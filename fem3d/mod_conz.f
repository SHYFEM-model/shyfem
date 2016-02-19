!
! routines for generic concentration
!
! contents :
!
! subroutine mod_conz_init(ncs,nkn,nlv)
!
! revision log :
!
! 09.11.2015    ggu     restructured with global values
! 19.02.2016    ggu     set default decay to zero
!
!******************************************************************

!==================================================================
        module mod_conz
!==================================================================

        implicit none

        integer, private, save :: ncs_conz = 0
        integer, private, save :: nkn_conz = 0
        integer, private, save :: nlv_conz = 0

        real, allocatable, save :: conzv(:,:,:)
        real, allocatable, save :: cnv(:,:)

	integer, save :: iconz = 0
        integer, save :: icall_conz = 0
        integer, save :: ninfo = 0
        integer, save :: level = 0
        integer, save :: iprogr = 0

        logical, save :: binfo = .true.

        real, save :: cref,rkpar,difmol,contau

        integer, save, allocatable :: idconz(:)
        integer, save :: ia_out(4)

        real, save, allocatable :: cdefs(:)
        real, save, allocatable :: tauv(:)
        real, save, allocatable :: massv(:)

	character*4, save :: what = 'conz'

        integer, parameter :: ndim_tau = 0
        real, save :: taupar(ndim_tau)
!        data taupar /0.,0.,0.,0.,0.,0.,0./

!==================================================================
	contains
!==================================================================

        subroutine mod_conz_init(ncs,nkn,nlv)

        integer ncs
        integer nkn
        integer nlv

        if( ncs == ncs_conz .and. nkn == nkn_conz 
     +		.and. nlv == nlv_conz ) return

        if( ncs > 0 .or. nkn > 0 .or. nlv > 0 ) then
          if( ncs == 0 .or. nkn == 0 .or. nlv == 0 ) then
            write(6,*) 'ncs,nkn,nlv: ',ncs,nkn,nlv
            stop 'error stop mod_conz_init: incompatible parameters'
          end if
        end if

        if( nkn_conz > 0 ) then
          deallocate(conzv)
          deallocate(cnv)
        end if

        ncs_conz = ncs
        nkn_conz = nkn
        nlv_conz = nlv

        if( nkn == 0 ) return

        allocate(conzv(nlv,nkn,ncs))
        allocate(cnv(nlv,nkn))

        end subroutine mod_conz_init

!==================================================================
        end module mod_conz
!==================================================================

