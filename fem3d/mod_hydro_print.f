!
! module for hydrodynamic print values
!
! revision log :
!
! 17.06.2016    ggu     wprv is now from 1:nlv
!
!============================================================
	module mod_hydro_print
!============================================================

	implicit none

        integer, private, save :: nkn_hydro_print = 0
        integer, private, save :: nlv_hydro_print = 0

        real, allocatable, save :: uprv(:,:)
        real, allocatable, save :: vprv(:,:)
        real, allocatable, save :: upro(:,:)
        real, allocatable, save :: vpro(:,:)
        real, allocatable, save :: wprv(:,:)
        real, allocatable, save :: up0v(:), vp0v(:)
        real, allocatable, save :: xv(:,:)

!============================================================
	contains
!============================================================

        subroutine mod_hydro_print_init(nkn,nlv)

        integer nkn, nlv

	if( nkn == nkn_hydro_print .and. nlv == nlv_hydro_print ) return

        if( nlv > 0 .or. nkn > 0 ) then
          if( nlv == 0 .or. nkn == 0 ) then
            write(6,*) 'nlv,nkn: ',nlv,nkn
            stop 'error stop mod_hydro_print_init: incompatible params'
          end if
        end if

	if( nkn_hydro_print > 0 ) then
          deallocate(uprv)
          deallocate(vprv)
          deallocate(upro)
          deallocate(vpro)
          deallocate(wprv)

          deallocate(up0v)
          deallocate(vp0v)

          deallocate(xv)
        end if

        nkn_hydro_print = nkn
        nlv_hydro_print = nlv

        if( nkn == 0 ) return

        allocate(uprv(nlv,nkn))
        allocate(vprv(nlv,nkn))
        allocate(upro(nlv,nkn))
        allocate(vpro(nlv,nkn))
        allocate(wprv(nlv,nkn))

	allocate(up0v(nkn))
	allocate(vp0v(nkn))

	allocate(xv(3,nkn))

        end subroutine mod_hydro_print_init

!============================================================
        end module mod_hydro_print
!============================================================

