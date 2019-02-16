
!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2018  Georg Umgiesser
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

c mercury routines
c
c contents :
c
c revision log :
c
c 15.05.2016    ggu     started mercury from bio3d
c
c notes :
c
c********************************************************************
c********************************************************************
c********************************************************************
c********************************************************************
c********************************************************************
c
c notes :
c
c State variables used: (mercury) -> Donata, please adjourn
c
c State variables used: (mer
c Hg0           81      1
c Hg2           82      2
c Hg3           83      3
c
c
c
c********************************************************************

!====================================================================
	module mercury
!====================================================================

	implicit none

	integer, parameter :: npstate = 3	!pelagic state variables
	integer, parameter :: nsstate = 3	!sediment state variables

	real, save, allocatable :: emp(:,:,:)	!pelagic state vector
	real, save, allocatable :: ems(:,:)	!sediment state vector

!       real, save, allocatable :: eload(:,:,:)   !atmospheric loading

	double precision, save :: da_out(4)

	integer, save :: iubp,iubs

!====================================================================
	end module mercury
!====================================================================

        subroutine mercury_module

c general interface to mercury module

        implicit none

	include 'femtime.h'

        real dt

	dt = dt_act
        call mercury3d(it,dt)

        end

c********************************************************************

	subroutine mercury3d(it,dt)

c eco-model cosimo

	use mod_diff_visc_fric
	use levels
	use basin
	use mercury
 
	implicit none

	integer it	!time in seconds
	real dt		!time step in seconds


        character*10 what

	integer k,i,l,lmax
	integer imerc
        integer id,idc
	integer nintp,ivar
	integer nbc
	real t,s
	real u,v

	real epela(npstate)
	real esedi(nsstate)

        real epload(npstate) !loading for each node
        real loadsup(npstate) 

! next array specifies boundary conditions if non are given

	real, save :: epbound(npstate) = (/0.5,0.5,0.05/)	!default bound cond.
	real, save :: epinit(npstate) = (/0.05,0.05,0.05/)	!default initial cond.
	real, save :: esinit(nsstate) = (/0.0,0.05,0.05/)	!default initial cond.
	real, save :: eploadin(nsstate) = (/0.10,0.01,0.0/)	!atm load kg/day

        real tpstot(npstate)              !for mass test
        real tsstot(nsstate)

	integer, save, allocatable :: idmerc(:)

	logical bsurf,bbottom
	integer iunit
	integer j
        real qrad
        real dtday
        real area,vol
        real getpar
        logical has_output,next_output
        logical has_output_d,next_output_d

        integer mode
        real ai,lsurf

        logical boxtype
        logical bcheck
        logical bresi,breact,bdecay
        integer ie,ii
        integer kspec
        integer nvar
        double precision dtime0,dtime
        real d          !element tickness 
        real depth      !element depth
        !real cbod,nh3,krear,sod
        real vel
        real windspeed,tempair
        real tday,t0,tsec
        real stp
        real mass
        real wsink
        real flux,conz(nlvdi)

        integer nbnds

        real, save :: rkpar,difmol
        integer, save :: icall = 0

c------------------------------------------------------------------
c	initial and boundary conditions  [mg/l]			??
c	initial loadings                 [g/m**2/sec]		??
c------------------------------------------------------------------

	breact = .true.		!use reactor

        what = 'mercury'

c-------------------------------------------------------------------
c initialization
c-------------------------------------------------------------------

	if( icall .le. -1 ) return

	if( icall .eq. 0 ) then
	  imerc = nint(getpar('imerc'))
	  if( imerc .le. 0 ) icall = -1
	  if( icall .le. -1 ) return
	  icall = 1

          call get_first_dtime(dtime0)

c         --------------------------------------------------
c	  initialize state variables
c         --------------------------------------------------

	  allocate(emp(nlvdi,nkndi,npstate))
	  allocate(ems(nkndi,nsstate))

           do i=1,npstate
            emp(:,:,i) = epinit(i)
          end do

           do i=1,nsstate
            ems(:,i) = esinit(i)
          end do

c         --------------------------------------------------
c	  initial conditions (only for pelagic part)
c         --------------------------------------------------

	  nvar = npstate
	  call mercury_init_file(dtime0,nvar,nlvdi,nlv,nkn,epinit,emp)

c         --------------------------------------------------
c	  set boundary conditions for all state variables
c         --------------------------------------------------

          nbc = nbnds()
          allocate(idmerc(nbc))
          idmerc = 0

          nintp = 2
	  nvar = npstate
          call bnds_init_new(what,dtime0,nintp,nvar,nkn,nlv
     +                          ,epbound,idmerc)

c         --------------------------------------------------
c	  initialize eco model
c         --------------------------------------------------

	  call mercury_init

c         --------------------------------------------------
c	  parameters for transport/diffusion resolution
c         --------------------------------------------------

          rkpar=getpar('chpar')
          difmol=getpar('difmol')

c         --------------------------------------------------
c	  initialize output 
c         --------------------------------------------------

	  call mercury_init_file_output
	  call mercury_write_file_output(dtime0)

	  write(6,*) 'mercury model initialized...'

	end if

c-------------------------------------------------------------------
c normal call
c-------------------------------------------------------------------

	wsink = 0.

c------------------------------------------------------------------
c       compute loading/unit surface	!FIXME -> this could be done only once
c------------------------------------------------------------------

	call mercury_surface_loading(npstate,eploadin,loadsup)

c-------------------------------------------------------------------
c       time management
c-------------------------------------------------------------------

        t0 = 0.
        dtday = dt / 86400.

        tsec = it
        tday = it / 86400. + t0         !time in days, FEM 0 is day t0

c	-------------------------------------------------------------------
c	loop on elements for reactor
c	-------------------------------------------------------------------

        mode = +1               !new time level for volume and depth

	if( breact ) then	!use reactor ?

	do k=1,nkn		!loop on nodes

          lmax = ilhkv(k)
	  call get_light(k,qrad)
          call get_sediment_values(k,flux,conz) !returns flux at bottom and conz in wc

          depth=0               !reinitialize the depth at each node
          do l=1,lmax
        
            call dvanode(l,k,mode,d,vol,area)   !gets depth, volume and area
            call getts(l,k,t,s)                 !gets temp and salt
            call getuv(l,k,u,v)                 !gets velocities u/v
            vel = sqrt(u*u+v*v)

            id = 1000*k+l
	    bsurf = (l==1)
	    bbottom = (l==lmax)
                
            depth=depth+d   !compute the depth at the centre of the element
        
            !write(6,*) bsurf,bbottom, l,d,k,depth,area,vol, 'mercury'
            !write(6,*) epload(:),'epload' !

            epload(:) = loadsup(:)*area
            epela(:) = emp(l,k,:)
            boxtype=.true.                
	    call mercury_react(id,bsurf,bbottom,boxtype,dtday,vol,d,k,t
     +				,qrad,epela,epload)
	    emp(l,k,:) = epela
          end do

          epload(:)=0
          esedi(:)=ems(k,:)
          boxtype=.false.
	  call mercury_react(id,bsurf,bbottom,boxtype,dtday,vol,d,k,t
     +				,qrad,esedi,epload)
	  ems(k,:) = esedi

	end do

	end if	!breact

c	-------------------------------------------------------------------
c	advection and diffusion
c	-------------------------------------------------------------------

	dtime = it
	call bnds_read_new(what,idmerc,dtime)

!$OMP PARALLEL PRIVATE(i)
!$OMP DO SCHEDULE(DYNAMIC)	

	do i=1,npstate

          call scal_adv(what,i
     +                          ,emp(1,1,i),idmerc
     +                          ,rkpar,wsink
     +                          ,difhv,difv,difmol)

          call tsmass (emp(1,1,i),1,nlvdi,tpstot(i)) !mass control

	end do

!$OMP END DO NOWAIT	
!$OMP END PARALLEL	

        do i=1,nsstate
          call scalmass(ems(1,i),0.1,tsstot(i))   !mass ctrl sed
	end do

c	-------------------------------------------------------------------
c	write of results
c	-------------------------------------------------------------------

	call mercury_write_file_output(dtime)

c	-------------------------------------------------------------------
c	debug output
c	-------------------------------------------------------------------

c	-------------------------------------------------------------------
c	end of routine
c	-------------------------------------------------------------------

	end

c*************************************************************
c*************************************************************
c*************************************************************

	subroutine mercury_init

! initializes mercury routines

	implicit none

	end

c*************************************************************

	subroutine mercury_init_file(dtime,nvar,nlvddi,nlv,nkn,val0,val)

c initialization of mercury from file

        implicit none

        double precision dtime
        integer nvar
        integer nlvddi
        integer nlv
        integer nkn
        real val0(nvar)
        real val(nlvddi,nkn,nvar)

        call tracer_file_init('mercury init','mercin',dtime
     +                          ,nvar,nlvddi,nlv,nkn,val0,val)

        end

c*************************************************************

        subroutine mercury_init_file_output

        use basin
        use levels
        use mercury

        implicit none

        integer nvar,id
        logical has_output_d

        call init_output_d('itmcon','idtcon',da_out)
        if( has_output_d(da_out) ) then
          nvar = npstate + nsstate
          call shyfem_init_scalar_file('merc',nvar,.false.,id)
          da_out(4) = id
        end if

        write(6,*) 'merc init file output done'

        end 

c*************************************************************

	subroutine mercury_write_file_output(dtime)

	use basin
	use levels
	use mercury

	implicit none

	double precision dtime

	integer nvar,id,idc,i
	logical next_output_d

        if( next_output_d(da_out) ) then
          id = nint(da_out(4))
          do i=1,npstate
            idc = 250 + i
            call shy_write_scalar_record(id,dtime,idc,nlvdi
     +                                          ,emp(1,1,i))
          end do
          do i=1,nsstate
            idc = 270 + i
            call shy_write_scalar_record(id,dtime,idc,1
     +                                          ,ems(1,i))
          end do
        end if

	end

c*************************************************************
c*************************************************************
c*************************************************************

      subroutine merc_euler(nstate,dt,vol,c,cold,cds)

! new c is computed
! cold is returned which is just c before call

      implicit none

      integer nstate            !number of state variables
      real dt                  !time step [day]
      real vol            !volume [m**3]
      real c(nstate)            !state variable [mg/L] == [g/m**3]
      real cold(nstate)      !old state variable (return)
      real cds(nstate)      !source term [g/day]

      integer i
      real volold,volnew      !volume [m**3]
      real mass            !mass [g]
      real mder            !derivative [g/day]

      volold = vol
      volnew = vol

      do i=1,nstate
        cold(i) = c(i)
        mass = c(i) * volold
        mder = cds(i)
        c(i) = ( mass + dt * mder ) / volnew
        !write(88,*)dt,'dt-day'
        !write(88,*) mder,'mder'
        !write(88,*) cds(i),'cds'
        !write(88,*) i, c(i),nstate,'c(i)'
        !if(c(i).lt.0.00001) c(i)=0.00001
      end do

      end
      
c*************************************************************
c*************************************************************
c*************************************************************

	subroutine mercury_surface_loading(npstate,eploadin,loadsup)

	use basin

	implicit none

	integer npstate
	real eploadin(npstate)
	real loadsup(npstate)

	integer mode
	integer k,l,i
	real d,vol,area,areatot

        loadsup=0
        mode=+1
        l=1

        areatot=0
        do k=1,nkn
          call dvanode(l,k,mode,d,vol,area)       !gets depth, volume, area
          areatot=areatot+area
        end do

        do i=1,npstate
          loadsup(i)=eploadin(i)/areatot          !loading per surf unit
        end do

        !write(6,*) loadsup(1),areatot,'loadsup1,areatot'
        !write(6,*) loadsup(2),areatot,'loadsup2,areatot'
        !write(6,*) loadsup(3),areatot,'loadsup3,areatot'

	end

c*************************************************************

