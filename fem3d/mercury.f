c
c $Id: bio3d.f,v 1.33 2008-10-10 09:29:54 georg Exp $
c
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
c nh3		71	1
c no3		72	2
c opo4		73	3
c phyto		74	4
c cbod		75	5
c do		76	6
c on		77	7
c op		78	8
c zoo		79	9
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

	integer, save :: ia_out(4)
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

! next array specifies boundary conditions if non are given

	real, save :: epbound(npstate) = (/0.,0.,0./)	!default bound cond.
	real, save :: epinit(npstate) = (/0.,0.,0./)	!default initial cond.

        real tpstot(npstate)              !for mass test
        real tsstot(nsstate)

	integer, save, allocatable :: idmerc(:)

	logical bsurf,bbottom
	integer iunit
	integer j
	real qrad
	real area,vol
	real getpar
	logical has_output,next_output
	logical has_output_d,next_output_d

        integer mode
        real ai,lsurf

	logical bcheck
	logical bresi,breact,bdecay
	integer ie,ii
	integer kspec
	integer itanf,nvar
	double precision dtime0,dtime
	real d
	real cbod,nh3,krear,sod
	real vel
	real windspeed,tempair
	real tday,t0,tsec
	real stp
        real mass
	real wsink

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

          dtime0 = itanf

c         --------------------------------------------------
c	  initialize state variables
c         --------------------------------------------------

	  allocate(emp(nlvdi,nkndi,npstate))
	  allocate(ems(nkndi,nsstate))

	  emp = 0.
	  ems = 0.

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

	  call get_first_time(itanf)
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

c	-------------------------------------------------------------------
c	loop on elements for biological reactor
c	-------------------------------------------------------------------

        mode = +1               !new time level for volume and depth

	if( breact ) then	!use reactor ?

	do k=1,nkn		!loop on nodes

          lmax = ilhkv(k)
	  call get_light(k,qrad)
	  esedi = ems(k,:)

          do l=1,lmax
            call dvanode(l,k,mode,d,vol,area)   !gets depth, volume and area
            call getts(l,k,t,s)                 !gets temp and salt
            call getuv(l,k,u,v)                 !gets velocities u/v
            vel = sqrt(u*u+v*v)

            id = 1000*k+l
	    bsurf = (l==1)
	    bbottom = (l==lmax)

	    epela = emp(l,k,:)

	    call mercury_react(id,bsurf,bbottom,dt,vol,d,vel,t,s,qrad
     +				,epela,esedi)

	    emp(l,k,:) = epela
          end do

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

        subroutine mercury_react(id,bsurf,bbottom,dt,vol,d,vel,t,s,qrad
     +                          ,epela,esedi)

! reactor for mercury

	implicit none

	integer, parameter :: npstate = 3	!pelagic state variables
	integer, parameter :: nsstate = 3	!sediment state variables

	integer id		!id of node and level
	logical bsurf		!is on surface?
	logical bbottom		!is on bottom?
	real dt			!time step [s]
	real vol		!volume of box [m**3]
	real d			!depth of box [m]
	real vel		!current velocity [m/s]
	real t			!water temperature [C]
	real s			!salinity [psu]
	real qrad		!radiation [W/m**2] (in/out)
	real epela(npstate)	!pelagic state variables (in/out)
	real esedi(nsstate)	!sediment state variables (in/out)

! qrad is input on top
! on return qrad should be the radiation on bottom (decay)

! this is an empty stub that does nothing ... ok for testing

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

	integer ishyff,nvar,id
	logical has_output,has_output_d
	real getpar

	ishyff = nint(getpar('ishyff'))

	  call init_output('itmcon','idtcon',ia_out)
          if( ishyff == 1 ) ia_out = 0
	  if( has_output(ia_out) ) then
	    call open_scalar_file(ia_out,nlv,npstate,'mer')
	    iubp = ia_out(4)
	    call open_scalar_file(ia_out,1,nsstate,'mes')
	    iubs = ia_out(4)
	  end if

          call init_output_d('itmcon','idtcon',da_out)
          if( ishyff == 0 ) da_out = 0
          if( has_output_d(da_out) ) then
	    nvar = npstate + nsstate
            call shyfem_init_scalar_file('merc',nvar,.false.,id)
            da_out(4) = id
          end if

	end 

c*************************************************************

	subroutine mercury_write_file_output(dtime)

	use basin
	use levels
	use mercury

	implicit none

	double precision dtime

	integer nvar,id,idc,i
	logical next_output,next_output_d

	if( next_output(ia_out) ) then
	  ia_out(4) = iubp
	  do i=1,npstate
	    idc = 250 + i
	    call write_scalar_file(ia_out,idc,nlvdi,emp(1,1,i))
	  end do

	  ia_out(4) = iubs
	  do i=1,nsstate
	    idc = 270 + i
	    call write_scalar_file(ia_out,idc,1,ems(1,i))
	  end do
        end if

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

