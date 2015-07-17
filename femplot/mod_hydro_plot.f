
!==================================================================
	module mod_hydro_plot
!==================================================================

	implicit none

	integer, private, save  :: nkn_hydro_plot = 0
	integer, private, save  :: nel_hydro_plot = 0

	real, allocatable, save :: uv(:)
	real, allocatable, save :: vv(:)
	real, allocatable, save :: uvnv(:)
	real, allocatable, save :: vvnv(:)
	real, allocatable, save :: usnv(:)
	real, allocatable, save :: vsnv(:)
	real, allocatable, save :: wsnv(:)

!==================================================================
	contains
!==================================================================

	subroutine mod_hydro_plot_init(nkn,nel)

	integer nkn,nel

        if( nkn == nkn_hydro_plot .and. nel == nel_hydro_plot ) return

        if( nel > 0 .or. nkn > 0 ) then
          if( nel == 0 .or. nkn == 0 ) then
            write(6,*) 'nel,nkn: ',nel,nkn
            stop 'error stop mod_hydro_plot_init: incompatible params'
          end if
        end if

        if( nkn_hydro_plot > 0 ) then
          deallocate(uv)
          deallocate(vv)
          deallocate(uvnv)
          deallocate(vvnv)
          deallocate(usnv)
          deallocate(vsnv)
          deallocate(wsnv)
        end if

        nel_hydro_plot = nel
        nkn_hydro_plot = nkn

	if( nkn == 0 ) return

        allocate(uv(nkn))
        allocate(vv(nkn))
        allocate(uvnv(nel))
        allocate(vvnv(nel))
        allocate(usnv(nel))
        allocate(vsnv(nel))
        allocate(wsnv(nkn))

	end subroutine mod_hydro_plot_init

!==================================================================
	end module mod_hydro_plot
!==================================================================

