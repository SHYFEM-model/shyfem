
!--------------------------------------------------------------------------
!
!    Copyright (C) 2015-2020  Georg Umgiesser
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

! revision log :
!
! 16.11.2015	ggu	changed VERS_7_3_14
! 19.02.2016	ggu	changed VERS_7_5_2
! 31.03.2017	ggu	changed VERS_7_5_24
! 09.05.2017	ggu	changed VERS_7_5_26
! 11.07.2017	ggu	changed VERS_7_5_30
! 05.12.2017	ggu	changed VERS_7_5_39
! 13.07.2018	ggu	changed VERS_7_4_1
! 25.10.2018	ggu	changed VERS_7_5_51
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62
! 31.01.2020	ggu	integrated parts from connectivity

!**************************************************************************

!==================================================================
        module mod_lagrange
!==================================================================

        implicit none

	!---------------------------------------------
	! internal parameters
	!---------------------------------------------

        integer, private, save :: nel_lagr = 0
        integer, private, save :: nlv_lagr = 0
        integer, private, save :: nbdy_lagr = 0

        integer*8 :: bm_kind			!indicates size of bitmap

        integer, save :: nbdymax = 0		!max number of particles
        integer, save :: ncust = 1		!number of custom properties

	!---------------------------------------------
	! particle info
	!---------------------------------------------

	type :: lagr_body
          double precision :: xi(3)             !internal coordinate
          double precision :: time		!time in seconds
	  real :: z				!z-coordinate
	  integer :: ie				!element number
	  integer :: l				!layer number
	end type lagr_body

	type :: lagr_entry
          integer(kind(bm_kind)) :: bitmap_in,bitmap_out  !for connectivity
	  double precision :: sinking		!sinking velocity
	  real, allocatable :: custom(:)	!custom properties
	  integer :: id				!id of particle
	  integer :: type			!type of particle
	  type(lagr_body) :: init		!initial properties
	  type(lagr_body) :: actual		!run time properties
	end type lagr_entry
	  
        type(lagr_entry), save, allocatable :: lgr_ar(:)

	!---------------------------------------------
	! parameters
	!---------------------------------------------

        logical, parameter :: blgrxi = .true.   !new version with xi coords

        logical, save :: blgrdebug = .false.
        logical, save :: blgrsurf = .false.
        logical, save :: blgr2d = .false.
        logical, save :: bconnect = .false.
        logical, save :: bcount = .false.	!counts particles in elements
        logical, save :: bsedim = .false.	!sediment 
        logical, save :: blarvae = .false.	!larvae 
        logical, save :: boilsim = .false.	!oil simulation
        logical, save :: bcompress = .false.	!compress particle numbers 
        logical, save :: bvdiff = .true.	!compute vertical diffusion
        logical, save :: bhdiff = .true.	!compute horizontal diffusion
        logical, save :: bbeach = .false.	!allow particle to beach on the shore

        integer, save :: ilagr                  !type of lagrangian simulation
        integer, save :: nbdy                   !total number of particles
        integer, save :: idbdy                  !max id used
	integer, save :: lunit                  !unit for messages
        integer, save :: ipvert                 !vertical release
        integer, save :: linbot                 !bottom layer for vert release
        integer, save :: lintop                 !surface layer for vert release
        real, save :: dripar                    !drifter parameter for inertia
        real, save :: stkpar                    !stokes drift parameter

        integer, save :: artype                 !special element type

        real, save :: azlgr                     !az parameter
        real, save :: tdecay                    !decay time - do not use
        real, save :: fall                      !vertical sinking velocity
        real, save :: rwhpar                    !horizontal diffusivity
        real, save :: lbeach                    !factor for particle beaching

	!---------------------------------------------
	! horizontal diffusivity
	!---------------------------------------------

        real, save, allocatable :: rwhvar(:)    !horizontal diffusivity (vary)

        !---------------------------------------------
        ! vertical diffusivity
        !---------------------------------------------

        real, save, allocatable :: wde(:,:)	!vertical diffusivity on elems
        real, save, allocatable :: dtvd(:)	!vertical random walk time step

	!---------------------------------------------
	! backtracking
	!---------------------------------------------

        integer, save :: nback
        logical, save :: bback

        integer, save, allocatable :: ie_back(:) !backtracking element number
        real, save, allocatable :: u_lag(:)      !backtracking x-velocity
        real, save, allocatable :: v_lag(:)      !backtracking y-velocity
        real, save, allocatable :: x_back(:)     !backtracking x-coordinate
        real, save, allocatable :: y_back(:)     !backtracking y-coordinate
        real, save, allocatable :: z_back(:)     !backtracking z-coordinate

	!--------------------------------------------------
	! fluxes and velocities
	!--------------------------------------------------

        real, save, allocatable :: flux2d(:,:)      !fluxes of sides
        real, save, allocatable :: flux3d(:,:,:)    !fluxes of sides (3d)
        real, save, allocatable :: vel_ie(:,:)      !velocities of sides
        real, save, allocatable :: vel3d_ie(:,:,:)  !velocities of sides (3d)
        real, save, allocatable :: dvert(:,:)

	!--------------------------------------------------
	! connectivity
	!--------------------------------------------------

        integer, save, allocatable :: i_count(:)
        double precision, save, allocatable :: t_count(:)
	real, save :: pld_day = 0.

!==================================================================
        contains
!==================================================================

        subroutine mod_lagrange_init(nel,nlv)

	integer nel,nlv

        if( nel == nel_lagr .and. nlv == nlv_lagr ) return

        if( nel > 0 .or. nlv > 0 ) then
          if( nel == 0 .or. nlv == 0 ) then
            write(6,*) 'nel,nlv: ',nel,nlv
	    stop 'error stop mod_lagrange_init: incompatible parameters'
          end if
        end if

        if( nel_lagr > 0 ) then
          deallocate(rwhvar)
          deallocate(ie_back)
          deallocate(u_lag)
          deallocate(v_lag)
          deallocate(x_back)
          deallocate(y_back)
          deallocate(z_back)
          deallocate(flux2d)
          deallocate(flux3d)
          deallocate(vel_ie)
          deallocate(vel3d_ie)
          deallocate(dvert)
          deallocate(i_count)
          deallocate(t_count)
          deallocate(wde)
          deallocate(dtvd)
	end if

        nel_lagr = nel
        nlv_lagr = nlv

        if( nel == 0 ) return

        allocate(rwhvar(nel))
        allocate(ie_back(nel))
        allocate(u_lag(nel))
        allocate(v_lag(nel))
        allocate(x_back(nel))
        allocate(y_back(nel))
        allocate(z_back(nel))

        allocate(flux2d(3,nel))
        allocate(flux3d(nlv,3,nel))
        allocate(vel_ie(3,nel))
        allocate(vel3d_ie(nlv,3,nel))
        allocate(dvert(3,nel))

        allocate(i_count(nel))
        allocate(t_count(nel))

        allocate(wde(0:nlv,nel))
        allocate(dtvd(nel))

        end subroutine mod_lagrange_init

!******************************************************************

        subroutine mod_lagrange_handle_alloc(nbdy)

	integer nbdy

        if( nbdy_lagr == 0 ) then
	  call mod_lagrange_init_body(max(nbdy,1024))
        else if( nbdy < nbdy_lagr/3 ) then
	  call mod_lagrange_init_body(nbdy_lagr/2)
	else if( nbdy > nbdy_lagr ) then
	  call mod_lagrange_init_body(max(nbdy,nbdy_lagr*2))
	end if

	end subroutine mod_lagrange_handle_alloc

!******************************************************************

        subroutine mod_lagrange_init_body(nbdy)

	integer nbdy

	integer ndim
	real time1,time2
	type(lagr_entry), allocatable :: paux(:)

        if( nbdy == nbdy_lagr ) return

	call cpu_time(time1)

	if( nbdy_lagr == 0 ) then		!first time
	  nbdy_lagr = nbdy
	  allocate(lgr_ar(nbdy))
        else if( nbdy == 0 ) then		!last time
	  deallocate(lgr_ar)
	else
          ndim = min(nbdy,nbdy_lagr)
          allocate(paux(nbdy))
          paux(1:ndim) = lgr_ar(1:ndim)
          call move_alloc(paux,lgr_ar)
	end if

	call cpu_time(time2)
	!write(lunit,*) 'alloc_lagr: ',nbdy,nbdy_lagr,time2-time1
	write(6,*) 'alloc_lagr: ',nbdy,nbdy_lagr,time2-time1

        nbdy_lagr = nbdy

        end subroutine mod_lagrange_init_body

!==================================================================
        end module mod_lagrange
!==================================================================

