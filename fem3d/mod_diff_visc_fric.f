	module mod_diff_visc_fric

	implicit none

	integer, private, save :: nkn_diff_visc_fric = 0
	integer, private, save :: nel_diff_visc_fric = 0
	integer, private, save :: nlv_diff_visc_fric = 0

	real, allocatable, save :: rfricv(:)
	real, allocatable, save :: czv(:)
	real, allocatable, save :: difhv(:,:)
	real, allocatable, save :: visv(:,:)
	real, allocatable, save :: difv(:,:)

	contains

!************************************************************

        subroutine mod_diff_visc_fric_init(nkn,nel,nlv)

        integer nkn, nel, nlv

        if( nkn == nkn_diff_visc_fric .and. nel == nel_diff_visc_fric
     +      .and. nlv == nlv_diff_visc_fric ) return

        if( nel > 0 .or. nkn > 0 .or. nlv > 0 ) then
          if( nel == 0 .or. nkn == 0 .or. nlv == 0 ) then
            write(6,*) 'nel,nkn,nlv: ',nel,nkn,nlv
	    stop 'error stop mod_diff_visc_fric_init: incompatible params'
          end if
        end if

        if( nkn_diff_visc_fric > 0 ) then
          deallocate(rfricv)
          deallocate(czv)
          deallocate(difhv)
          deallocate(visv)
          deallocate(difv)
        end if

        nkn_diff_visc_fric = nkn
        nel_diff_visc_fric = nel
        nlv_diff_visc_fric = nlv

        if( nkn == 0 ) return

        allocate(rfricv(nel))
        allocate(czv(nel))
        allocate(difhv(nlv,nel))
        allocate(visv(0:nlv,nkn))
        allocate(difv(0:nlv,nkn))

        end subroutine mod_diff_visc_fric_init

!************************************************************

        end module mod_diff_visc_fric
