
!--------------------------------------------------------------------------
!
!    Copyright (C) 2006,2008,2010,2012,2014-2015,2014-2015  Georg Umgiesser
!    Copyright (C) 2018-2019  Georg Umgiesser
!    Copyright (C) 2006  Francesca De Pascalis
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

c atoxi3d - Toxical routines from ARPAV - shell
c
c contents :
c
c revision log :
c
c 15.02.2006	ggu&fdp	new routine atoxi3d for ARPAV (from bio3d)
c 23.03.2006	ggu	ntot eliminated
c 23.03.2006	ggu	changed time step to real
c 17.04.2008	ggu	new open boundary conditions 
c 22.04.2008	ggu	advection parallelized
c 23.04.2008	ggu	call to bnds_set_def() changed
c 09.10.2008	ggu	new call to confop
c 23.03.2010	ggu	changed v6.1.1
c 30.03.2012	ggu	changed VERS_6_1_51
c 21.10.2014	ggu	converted to new boundary treatment
c 23.12.2014	ggu	changed VERS_7_0_11
c 19.01.2015	ggu	changed VERS_7_1_3
c 26.02.2015	ggu	changed VERS_7_1_5
c 17.07.2015	ggu	changed VERS_7_1_53
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 30.07.2015	ggu	changed VERS_7_1_83
c 03.04.2018	ggu	changed VERS_7_5_43
c 16.02.2019	ggu	changed VERS_7_5_60
c 21.03.2022	ggu	upgraded to da_out
c 22.03.2022	ggu	upgraded to new cmed
c
c notes :
c
c State variables used: (Wasp)
c
c todo :
c
c - * FIXME -> number of levels nlvdim, and not 1 (done)
c - wind speed and air temp is not available -> introduce (wmeteo)
c
c********************************************************************
c********************************************************************
c********************************************************************
c********************************************************************
c********************************************************************
c
c***********************************************

	subroutine atoxi3d

c toxi module ARPAV

	use mod_diff_visc_fric
	use levels
	use basin

	implicit none

	include 'param.h'
	include 'mkonst.h'

	integer nstate
	parameter( nstate = 1 )

	real e(nlvdim,nkndim,nstate)	!state vector
	real eb(nlvdim,nkndim,nstate)	!boundary vector of state vectors
	save e,eb

        real tstot(nstate)              !for mass test

        character*10 what

	integer k,i,l,lmax
	integer itoxi
        integer id
	integer nintp,ivar
	integer nbc
	real t,s
	real u,v
	real dt		!time step in seconds

	integer, save, allocatable :: idtoxi(:)

	real eaux(nstate)
	real einit(nstate)
	real ebound(nstate)
	save einit,ebound

	integer iunit
	integer j
	real rlux,rluxaux,itot,fday
	real dtt,dttday
	real area,vol
	real oxysat
	real getpar
	integer iround
	integer ieint,ipint
	integer nvar
	double precision dtime0,dtime

        integer mode
        real ai,lsurf

	logical bcheck
	logical bsedim
	integer ie,ii
	integer kspec
	integer idc,is
	real d
	real cbod,nh3,krear,sod
	real vel
	real windspeed,tempair
	real tday,t0,tsec
	real stp
        real mass
	real wsink

	integer nbnds

	double precision, save :: da_out(4)

	logical has_output_d,next_output_d

	integer, save :: icall = 0
	integer iespecial,inspecial
	save iespecial,inspecial
	real rkpar,difmol
	save rkpar,difmol

c------------------------------------------------------------------
c	initial and boundary conditions  [mg/l]			??
c
c 	 data einit /0.0, 0., 0.0, 0.0, 0.,   0.,0.,0.0,0.0/
 	 data einit /0.0/
 	 data ebound /0.0/
c
c------------------------------------------------------------------

        what = 'toxi'

c-------------------------------------------------------------------
c initialization
c-------------------------------------------------------------------

	if( icall .le. -1 ) return

	if( icall .eq. 0 ) then
	  itoxi = iround(getpar('itoxi'))
	  if( itoxi .le. 0 ) icall = -1
	  if( icall .le. -1 ) return
	  icall = 1

	  stop 'error stop atoxi3d: not adapted yet for new framework'

c         --------------------------------------------------
c	  initialize state variables with einit
c         --------------------------------------------------

	  do k=1,nkn		!loop on nodes
            lmax = ilhkv(k)
            do l=1,lmax
	      do is=1,nstate
	        e(l,k,is) = einit(is)
              end do
	    end do
          end do

c         --------------------------------------------------
c	  initialize state variables from external file
c         --------------------------------------------------

          call inicfil('toxi',e,nstate)

c         --------------------------------------------------
c	  set boundary conditions for all state variables
c         --------------------------------------------------

          nbc = nbnds()
          allocate(idtoxi(nbc))
          idtoxi = 0

	  call get_first_dtime(dtime0)
	  nintp = 2
	  nvar = nstate
          call bnds_init_new(what,dtime0,nintp,nvar,nkn,nlv
     +                          ,ebound,idtoxi)

	  !call bnds_init(what,tox3dn,nintp,nstate,nb3dim,toxiarr,ebound)
	  !call bnds_set_def(what,nb3dim,toxiarr) !nvar != 1
	  !call bnds_print('init of '//what,nb3dim,toxiarr)

c         --------------------------------------------------
c	  initialize eco model
c         --------------------------------------------------

	  call atoxi_ini

c         --------------------------------------------------
c	  parameters for transport/diffusion resolution
c         --------------------------------------------------

          rkpar=getpar('chpar')
          difmol=getpar('difmol')

c         --------------------------------------------------
c	  initialize output 
c         --------------------------------------------------

          call init_output_d('itmcon','idtcon',da_out)
	  !call increase_output_d(da_out)
          if( has_output_d(da_out) ) then
            nvar = nstate
            call shyfem_init_scalar_file('tox',nvar,.false.,id)
            da_out(4) = id
          end if

	  write(6,*) 'toxi model initialized...'

	end if

c-------------------------------------------------------------------
c normal call
c-------------------------------------------------------------------

	kspec = -100
	!kspec = 930
	bcheck = .true.
	bcheck = .false.
	wsink = 0.

c	-------------------------------------------------------------------
c	time management
c	-------------------------------------------------------------------

	t0 = 0.
	call get_act_dtime(dtime)
	call get_timestep(dt)
	tsec = dtime

c	-------------------------------------------------------------------
c	loop on nodes for biological reactor
c	-------------------------------------------------------------------

        mode = +1               !new time level for volume and depth

	do k=1,nkn		!loop on nodes

          lmax = ilhkv(k)
          !call getmeteo(k,tempair,windspeed)    !meteo FIXME
          !call wmeteo(tempair,windspeed)      !meteo FIXME
          rlux = 1.

          do l=1,lmax
            call dvanode(l,k,mode,d,vol,area)   !gets depth, volume and area
            call getts(l,k,t,s)                 !gets temp and salt
            call getuv(l,k,u,v)                 !gets velocities u/v
            vel = sqrt(u*u+v*v)

            id = 1000*k+l

	    do is=1,nstate
	      eaux(is) = e(l,k,is)
	    end do

	    !call atoxi(id,tsec,dt,d,t,eaux)

	    do is=1,nstate
	      e(l,k,is) = eaux(is)
	    end do
          end do

	end do

c	-------------------------------------------------------------------
c	advection and diffusion
c	-------------------------------------------------------------------

	if( bcheck ) call check_toxi('before advection',e)

	call bnds_read_new(what,idtoxi,dtime)

!$OMP PARALLEL PRIVATE(i)
!$OMP DO SCHEDULE(DYNAMIC)

	do is=1,nstate

          call scal_adv(what,is
     +                          ,e(1,1,is),idtoxi
     +                          ,rkpar,wsink
     +                          ,difhv,difv,difmol)

          call tsmass (e(1,1,is),1,nlvdim,tstot(is)) !mass control

	end do

!$OMP END DO NOWAIT
!$OMP END PARALLEL

	if( bcheck ) call check_toxi('after advection',e)
	!write(86,*) dtime,tstot(1)

c	-------------------------------------------------------------------
c	write of results (file BIO)
c	-------------------------------------------------------------------

        if( next_output_d(da_out) ) then
          id = nint(da_out(4))
          idc = 120
	  do is=1,nstate
	    idc = idc + 1
            call shy_write_scalar_record(id,dtime,idc,nlvdi,e(:,:,is))
	  end do
        end if

c	call toxi_av_shell(e)		!aver/min/max of state vars

	if( bcheck ) call check_toxi('at end',e)

c	-------------------------------------------------------------------
c	debug output
c	-------------------------------------------------------------------

c        k = 100
c        l = 1
c        call getts(l,k,t,s)
c        call writeet(95,it,k,l,e,t,s,nlvdim,nkndim,nstate)

c	-------------------------------------------------------------------
c	end of routine
c	-------------------------------------------------------------------

	end

c*************************************************************

	subroutine writeet(iunit,it,k,l,e,t,s,nlvdim,nknddi,nstate)

c formatted write for debug

	implicit none

	integer iunit
	integer it,k,l
	integer nlvdim,nknddi,nstate
	real e(nlvdim,nknddi,nstate)
	real t
	real s

	integer i

	write(iunit,'(i10,11f12.4)') it,
     +			(e(l,k,i),i=1,nstate),
     +			t,s

	end

c*************************************************************

	subroutine toxi_av_shell(e)

c computes and writes average/min/max of bio variables
c
c id = 260
c
c e(1) average	== 261
c e(1) min	== 262
c e(1) max	== 263
c e(2) average	== 264
c ...

	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

c parameter

	integer, parameter :: nstate = 1

	real e(nlvdi,nkn,nstate)	!state vector

c local
	integer itsmed
	integer id,nvar,idc,is
	double precision dtime
	double precision rr
c function
	real getpar
	logical has_output_d,is_over_output_d,next_output_d
c save
        double precision, save :: da_out(4)
        double precision, save, allocatable :: bioacu(:,:,:)
        real, save, allocatable :: biomin(:,:,:)
        real, save, allocatable :: biomax(:,:,:)
        real, save, allocatable :: raux(:,:)

	integer, save :: icall = 0
	integer, save :: nr = 0

	if( icall .lt. 0 ) return

	if( icall .eq. 0 ) then

          itsmed=nint(getpar('itsmed'))
          if( itsmed .le. 0 ) then
            icall = -1
            return
          end if

          call init_output_d('itmcon','idtcon',da_out)
          call increase_output_d(da_out)
          if( has_output_d(da_out) ) then
            nvar = nstate*3
            call shyfem_init_scalar_file('bav',nvar,.false.,id)
            da_out(4) = id
          end if

          allocate(biomin(nlvdi,nkn,nstate),biomax(nlvdi,nkn,nstate))
	  allocate(bioacu(nlvdi,nkn,nstate))
          allocate(raux(nlvdi,nkn))

	  do is=1,nstate
            call cmed_reset(nr,bioacu(:,:,is)
     +			,biomin(:,:,is),biomax(:,:,is))
	  end do

	  icall = 1
	end if

        if( .not. is_over_output_d(da_out) ) return

        nr = nr + 1
	do is=1,nstate
          call cmed_accum(e(:,:,is),bioacu(:,:,is)
     +			,biomin(:,:,is),biomax(:,:,is))
	end do

        if( .not. next_output_d(da_out) ) return

        id = nint(da_out(4))
        call get_act_dtime(dtime)
        rr=1./nr

        idc = 260
	do is=1,nstate
	  idc = idc + 1
          raux = bioacu(:,:,is) * rr
          call shy_write_scalar_record(id,dtime,idc,nlvdi,raux)
          raux = biomin(:,:,is)
          call shy_write_scalar_record(id,dtime,idc,nlvdi,raux)
          raux = biomax(:,:,is)
          call shy_write_scalar_record(id,dtime,idc,nlvdi,raux)
	end do

	do is=1,nstate
          call cmed_reset(nr,bioacu(:,:,is)
     +			,biomin(:,:,is),biomax(:,:,is))
	end do

	end

c*************************************************************

	subroutine check_toxi(title,e)

c checks bio vars

	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'

	integer nstate
	parameter( nstate = 1 )

	character*(*) title
	real e(nlvdim,nkndim,nstate)	!state vector


        character*20 text
	integer i

	write(6,*) 'check_toxi: ',title

        text = '*** bio check e     '
	do i=1,nstate
          write(text(18:19),'(i2)') i
          call check2Dr(nlvdim,nlv,nkn,e(1,1,i),0.,1.e+20,text,title)
	end do

	end

c*************************************************************

