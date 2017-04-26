
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

	!---------------------------------------------
	! particle info
	!---------------------------------------------

        type  :: lagr_entry
          sequence
          double precision :: xi(3)             !internal coordinate
          double precision :: sv                !sinking velocity
          real    :: xst                        !initial x-coordinate
          real    :: yst                        !initial y-coordinate
          real    :: zst                        !initial z-coordinate
          real    :: x                          !x-coordinate
          real    :: y                          !y-coordinate
          real    :: z                          !z-coordinate
          real    :: c                          !custom property for plot
          real    :: tin                        !time of release
          integer :: est                        !initial element number
          integer :: id                         !id of particle
          integer :: ty                         !type of particle
          integer :: ie                         !element number
          integer :: l                          !layer number
          integer :: dummy                      !dummy argument for sequence
          !integer(kind(bm_kind)) :: bitmap_in,bitmap_out  !uncomment for connectivity
        end type lagr_entry

        type(lagr_entry), save, allocatable :: lgr_ar(:)

	!---------------------------------------------
	! parameters
	!---------------------------------------------

        logical, parameter :: blgrxi = .true.    !new version with xi coords

        logical, save :: blgrdebug = .false.
        logical, save :: blgrsurf = .false.
        logical, save :: bconnect = .false.
        logical, save :: bcount = .false.	!counts particles in elements
        logical, save :: bsedim = .false.	!sediment 
        logical, save :: blarvae = .false.	!larvae 
        logical, save :: boilsim = .false.	!oil simulation
        logical, save :: bcompress = .false.	!compress particle numbers 

        integer, save :: ilagr                  !type of lagrangian simulation
        integer, save :: nbdy                   !total number of particles
        integer, save :: idbdy                  !max id used
        integer, save :: lunit                  !unit for messages
        integer, save :: ipvert                 !vertical release
        integer, save :: linbot                 !bottom layer for vert release
        integer, save :: lintop                 !surface layer for vert release

        integer, save :: artype                 !special element type

        real, save :: azlgr                     !az parameter
        real, save :: tdecay                    !decay time - do not use
        real, save :: fall                      !vertical sinking velocity
        real, save :: rwhpar                    !horizontal diffusivity

	!---------------------------------------------
	! horizontal diffusivity
	!---------------------------------------------

        real, save, allocatable :: rwhvar(:)    !horizontal diffusivity (vary)

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
	! conncectivity
	!--------------------------------------------------

        integer, save, allocatable :: i_count(:)
        double precision, save, allocatable :: t_count(:)

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

        end subroutine mod_lagrange_init

!******************************************************************

        subroutine mod_lagrange_handle_alloc(nbdy)

	integer nbdy

        if( nbdy_lagr == 0 ) then
	  call mod_lagrange_init_body(1024)
        else if( nbdy < nbdy_lagr/3 ) then
	  call mod_lagrange_init_body(nbdy_lagr/2)
	else if( nbdy > nbdy_lagr ) then
	  call mod_lagrange_init_body(nbdy_lagr*2)
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

