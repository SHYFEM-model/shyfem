
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
! 08.11.2020    ggu&riz	started new module
! 13.11.2020    ggu	bug fix: cc is a fraction, not a percentage
!
!******************************************************************

!==================================================================
        module shyice_model
!==================================================================

	implicit none

	integer, parameter :: nvars = 44
	integer, save :: nkn_ice = 0
	logical, save :: bice = .false.
	double precision, save :: idtice = -1.
	double precision, parameter :: iceth0 = 1.D-5
	double precision, parameter :: Tkelvin = 273.15

	integer, save :: iunit = 0

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

	idtice = dgetpar('idtice')

	bice = idtice >= 0		!we use the ice model
	if( .not. bice ) return

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

	write(6,*) 'shyice is active: idtice'

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
	integer, parameter :: i = 0
	real sh
	double precision CL,Fsd_cloud,P_rate,qa,qs
	double precision Ta,ua,Sw,deltat
	double precision iceth
	double precision vars(nvars)

	if( .not. bice ) return

	bdebug = ( mod(k,50) == 0 )
	bdebug = .true.
	bdebug = .false.

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
	call ice_set_hmix(hm)

	if( bdebug ) then
	  write(6,*) k
	  write(6,*) Cl,Fsd_cloud,P_rate,qa,Ta,Ua
	  write(6,*) deltat
	  write(6,*) iceth
	  write(6,*) hm,tm,sm
	  write(6,*) vars(:)
	end if

	call ice_run(CL,Fsd_cloud,P_rate,qa,qs
     +				,Ta,Ua,Sw,deltat
     +                          ,i,vars,iceth)

	icevars(:,k) = vars(:)
	icethick(k) = iceth
	tm = vars(1) - Tkelvin
	sm = Sw

	if( bdebug ) then
	  write(6,*) k
	  write(6,*) iceth
	  write(6,*) hm,tm,sm
	  write(6,*) vars(:)
	  stop
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

	subroutine shyice_init_output
	use shyice_model
	if( .not. bice ) return
	iunit = 456
	open(iunit,file='icethickness.txt'
     +			,status='unknown',form='formatted')
	end

	subroutine shyice_write_output
	use shyice_model
	character*20 aline
	if( .not. bice ) return
	call get_act_timeline(aline)
	write(iunit,*) aline,icethick(1)
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
!	subroutine ice_set_hmix(hm)
!	use shyice_model
!	real hm
!	end
!
!*****************************************************************

