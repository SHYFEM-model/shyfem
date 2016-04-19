
!==================================================================
	module mod_plot3d
!==================================================================

	implicit none

	integer, private, save  :: nkn_plot3d = 0
	integer, private, save  :: nel_plot3d = 0
	integer, private, save  :: nlv_plot3d = 0
	integer, private, save  :: np_plot3d = 0

	real, allocatable, save :: fvlv(:,:)
	real, allocatable, save :: wauxv(:,:)
	real, allocatable, save :: het3v(:,:)
	real, allocatable, save :: p3(:,:)

!==================================================================
	contains
!==================================================================

	subroutine mod_plot3d_init(nkn,nel,nlv,np)

	integer nkn,nel,nlv,np

        if( nkn == nkn_plot3d .and. nel == nel_plot3d .and.
     +		nlv == nlv_plot3d .and. np == np_plot3d ) return

        if( nel > 0 .or. nkn > 0 .or. nlv > 0 .or. np > 0 ) then
          if( nel == 0 .or. nkn == 0 .or. nlv == 0 .or. np == 0 ) then
            write(6,*) 'nel,nkn,nlv,np: ',nel,nkn,nlv,np
            stop 'error stop mod_plot3d_init: incompatible parameters'
          end if
        end if

        if( nkn_plot3d > 0 ) then
          deallocate(fvlv)
          deallocate(wauxv)
          deallocate(het3v)
          deallocate(p3)
        end if

        nkn_plot3d = nkn
        nel_plot3d = nel
        nlv_plot3d = nlv
        np_plot3d = np

	if( nkn == 0 ) return

        allocate(fvlv(nlv,nkn))
        allocate(wauxv(0:nlv,nkn))
        allocate(het3v(nlv,nel))
        allocate(p3(nlv,np))

	end subroutine mod_plot3d_init

!==================================================================
	end module mod_plot3d
!==================================================================

