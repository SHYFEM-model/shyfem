
!--------------------------------------------------------------------------
!
!    Copyright (C) 2020  Georg Umgiesser
!    Copyright (C) 2020  Rasa Idzelyte
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
!
! thermodynamic ice model
!
! revision log :
!
! 08.11.2020    ggu&riz	started new module shyice_model
! 13.11.2020    ggu	bug fix: cc is a fraction, not a percentage
! 13.11.2020    ggu	new number of variables is 46
! 14.11.2020    ggu	output routines finished
! 22.11.2020    ggu	some code for debugging added
!
!******************************************************************

!==================================================================
        module shyice_model
!==================================================================

	implicit none

	integer, parameter :: nvars = 46
	integer, save :: nkn_ice = 0
	logical, save :: bice = .false.
	integer, save :: icemod = 0
	double precision, save :: idtice = 0.
	double precision, save :: da_out(4) = 0.
	double precision, parameter :: iceth0 = 1.D-5
	double precision, parameter :: Tkelvin = 273.15

	integer, save :: kdebug = 1	!if > 0 writes icethicknes of node
	integer, save :: iunit1 = 0

	double precision, save, allocatable :: icevars(:,:)
	real, save, allocatable :: icethick(:)

!==================================================================
        contains
!==================================================================

	subroutine shyice_model_init(nkn)

! allocates ice model arrays

	integer nkn

	if( nkn == nkn_ice ) return

	if( nkn_ice > 0 ) then
	  deallocate(icevars)
	  deallocate(icethick)
	end if

	nkn_ice = nkn

	if( nkn == 0 ) return

	allocate(icevars(nvars,nkn))
	allocate(icethick(nkn))

	icevars = 0.
	icethick = 0.

	end subroutine shyice_model_init

!==================================================================
        end module shyice_model
!==================================================================

	subroutine shyice_init

! initializes ice model

	use shyice_model
	use basin

	implicit none

	logical bicefile
	integer k
	double precision dgetpar

	icemod = nint(dgetpar('icemod'))
	idtice = dgetpar('idtice')

	bice = icemod > 0		!we use the ice model
	if( .not. bice ) return

	if( idtice > 0 ) then
	  write(6,*) 'Cannot yet handle idtice /= 0'
	  stop 'error stop shyice_init: idtice /= 0'
	end if

	call meteo_has_ice_file(bicefile)
	if( bicefile ) then
	  write(6,*) 'cannot handle both ice model and ice file'
	  write(6,*) 'you can either use the ice forcing file'
	  write(6,*) 'or you can run the ice model'
	  stop 'error stop shyice_init: not both file and model'
	end if

	call shyice_model_init(nkn)
	call ice_check_nvars(nvars)

	do k=1,nkn
	  call ice_init_vars(icevars(:,k))
	end do

	write(6,*) 'shyice is active: ',icemod,idtice

	end

!*****************************************************************

	subroutine shyice_run_all_nodes

! run ice model for all nodes

	use shyice_model
	use basin

	implicit none

	integer k
	integer mode,layer
	real qs,ta,rh,twb,uw,cc,p
	real r,e,eeff
	real tw,salt
	real dt,hm
	double precision vars(nvars)

	real depnode

	if( .not. bice ) return

	mode = +1			!new time step
	layer = 1			!only first layer
	call get_timestep(dt)

	do k=1,nkn

	  call getts(layer,k,tw,salt)
	  hm = depnode(layer,k,mode)

	  call meteo_get_heat_values(k,qs,ta,rh,twb,uw,cc,p)
	  call get_pe_values(k,r,e,eeff)

	  call shyice_run(k,qs,ta,rh,uw,cc,p,r,hm,tw,salt,dt)

	end do

	end

!*****************************************************************

	subroutine shyice_run(k,qsrad,tair,rh,uw,cc,p,r,hm,tm,sm,dt)

! run ice model for one node

	use shyice_model
	use basin

	implicit none

	integer k
	real qsrad			!solar radiation [W/m**2]
	real tair			!air temperature [C]
	real rh				!relative humidity [%]
	real uw				!wind speed [m/s]
	real cc				!cloud cover [0-1]
	real p				!atmospheric pressure [mbar]
	real r				!precipitation [mm/day]
	real hm				!depth of mixed layer [m]
	real tm				!temperature of mixed layer [C]
	real sm				!salinity of mixed layer [psu]
	real dt				!time step [s]

	logical bdebug
	integer, parameter :: i0 = 0
	integer, save :: icall = 0
	integer i
	real sh
	double precision Cl,Fsd_cloud,P_rate,qa,qs
	double precision Ta,ua,Sw,deltat
	double precision hdm,tdm,sdm
	double precision iceth
	double precision vars(nvars)
	double precision atime,dtime

	if( .not. bice ) return

	!Cl = cc * 100.
	Cl = cc 
	Fsd_cloud = qsrad
	P_rate = r / 86400.
	call rh2sh(tair,p,rh,sh)
	qa = sh
	qs = sh
	Ta = tair + Tkelvin
	Ua = uw

	deltat = dt

	vars(:) = icevars(:,k)
	iceth = icethick(k)
	vars(1) = tm + Tkelvin		!temp in mixed layer is vars(1)
	Sw = sm
        hdm = hm
	!call ice_set_hmix(hdm)

!---------------------------------------------------------------
! next is for debugging ice model
! chose node you want to debug and time after which to debug
! inside the ice model, debug sections can be inserted like this:
!
!	if( bice_debug ) then
!	  !debug statements
!	end if

	call get_act_dtime(dtime)		!relative time
	call get_absolute_act_time(atime)	!absolute time
	call set_ice_debug(.false.)		!set debug to false

	bdebug = ( k .eq. 1108 ) !.and. dtime > 128514000 )	!do we have to debug? k .eq. 100000 .and. dtime > 279558900
	if ( bdebug ) then
	  call set_ice_debug(bdebug)
	end if 

	if( bdebug ) then
          call ice_debug2(i,k,Cl,Fsd_cloud,P_rate,qa,qs
     &			,Ta,Ua,Sw,deltat,hdm,iceth,0,nvars,vars)
	end if
	if ( bdebug ) then
	  write(501,*) tm
	end if
	
	!if ( tm < -1.0 ) then
	!  tm = -5.0
	!  vars(1) = tm + Tkelvin
	!end if
!---------------------------------------------------------------
	call ice_run(k,Cl,Fsd_cloud,P_rate,qa,qs
     +				,Ta,Ua,Sw,deltat,hdm
     +                          ,i0,vars,iceth)

	icevars(:,k) = vars(:)
	icethick(k) = iceth
	
	tm = vars(1) - Tkelvin
	if ( tm < -1.0 ) then
	  tm = -1.0
	end if
    
	sm = Sw
	if ( bdebug ) then
	  write(502,*) tm
	end if
	!if ( bdebug ) then
	!  if ( iceth > 1.D-5 .and. iceth < 0.05 ) then
	!    write(501,*) iceth,tm,Ta
	!  end if
	!end if
	if( bdebug ) then
	  icall = icall + 1
	  i = icall
	  hdm = hm
	  tdm = vars(1)
          call ice_debug2(i,k,Cl,Fsd_cloud,P_rate,qa,qs
     &			,Ta,Ua,Sw,deltat,hdm,iceth,1,nvars,vars)
	end if

	if( iceth > iceth0 ) then	!ice cover is either 0 or 1
	  call set_ice_cover(k,1.)
	else
	  call set_ice_cover(k,0.)
	end if

	end

!*****************************************************************

	subroutine shyice_get_tsurf(k,t0)

	use shyice_model

	implicit none

	integer k
	real t0

	t0 = icevars(2,k)

	end

!*****************************************************************

	subroutine shyice_is_active(bicemodel)

! is ice model running?

	use shyice_model

	implicit none

	logical bicemodel

	bicemodel = bice

	end

!*****************************************************************
! output routines
!*****************************************************************

	subroutine shyice_init_output

	use shyice_model

	implicit none

	integer id
	integer, parameter :: nvar = 1
	logical, parameter :: b2d = .true.

        logical has_output_d

	if( .not. bice ) return

        da_out = 0

        call init_output_d('itmcon','idtcon',da_out)
        if( has_output_d(da_out) ) then
          call shyfem_init_scalar_file('ice',nvar,b2d,id)
          da_out(4) = id
        end if

	if( kdebug > 0 ) then
	  iunit1 = 555
	  call find_unit(iunit1)
	  open(iunit1,file='icethickness.txt'
     +			,status='unknown',form='formatted')
	end if

	end

!*****************************************************************

	subroutine shyice_write_output

	use shyice_model

	implicit none

	integer id,idvar,nlvdi
	double precision dtime
	character*20 aline

	logical next_output_d

	if( .not. bice ) return

        idvar = 86	!ice thickness
	nlvdi = 1

        if( next_output_d(da_out) ) then
          call get_act_dtime(dtime)
          id = nint(da_out(4))
          call shy_write_scalar_record(id,dtime,idvar,nlvdi,icethick)
        end if

	if( kdebug > 0 ) then
	  call get_act_timeline(aline)
	  write(iunit1,*) aline,icethick(kdebug)
	end if

	end

!*****************************************************************
! debug routines
!*****************************************************************

        subroutine ice_debug(i,k,Cl,Fsd_cloud,P_rate,qa,Ta,Ua           &
     &                          ,iceth,hm,tm,sm)
        implicit none
        integer i,k
        double precision Cl,Fsd_cloud,P_rate,qa,Ta,Ua
        double precision iceth,hm,tm,sm

        write(765,'(i8,10e14.6)') i                                     &
     &                  ,Cl,Fsd_cloud,P_rate,qa,Ta,Ua                   &
     &                  ,iceth,hm,tm,sm

        end

        subroutine ice_debug2(i,k,Cl,Fsd_cloud,P_rate,qa,qs
     &				,Ta,Ua,Sw,deltat,hdm,iceth
     &				,j,nvars,vars)

	implicit none
        integer i,k,j
        double precision Cl,Fsd_cloud,P_rate,qa,qs
        double precision Ta,Ua,Sw,deltat,hdm,iceth
	integer nvars
	double precision vars(nvars)

	write(777,*) Cl,Fsd_cloud,P_rate,qa,qs
     +			,Ta,Ua,Sw,deltat,hdm,iceth
     +			,j,nvars,vars

	end

!*****************************************************************
! aux routines for compilation - to be deleted
!*****************************************************************
!
!	subroutine ice_check_nvars(n)
!	use shyice_model
!	integer n
!	end
!
!	subroutine ice_init_vars(vars)
!	use shyice_model
!	double precision vars(nvars)
!	vars = 0.
!	end
!
!	subroutine ice_run(CL,Fsd_cloud,P_rate,qa,qs
!     +				,Ta,Ua,Sw,deltat
!     +                          ,i,vars,iceth)
!	use shyice_model
!	double precision CL,Fsd_cloud,P_rate,qa,qs
!	double precision Ta,Ua,Sw,deltat
!	double precision vars(nvars),iceth
!	integer i
!	end
!
!	subroutine ice_set_hmix(hdm)
!	use shyice_model
!	double precision hdm
!	end
!
!*****************************************************************

