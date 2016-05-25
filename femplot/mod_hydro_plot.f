
!==================================================================
	module mod_hydro_plot
!==================================================================

	implicit none

	integer, private, save  :: nkn_hydro_plot = 0
	integer, private, save  :: nel_hydro_plot = 0
	integer, private, save  :: nx_hydro_plot = 0
	integer, private, save  :: ny_hydro_plot = 0
	integer, private, save  :: nxy_limit = 1000	!0 for no limit

	real, allocatable, save :: uvnode(:)	!variable in x on node
	real, allocatable, save :: vvnode(:)	!variable in y on node
	real, allocatable, save :: uvelem(:)	!variable in x on element
	real, allocatable, save :: vvelem(:)	!variable in y on elemen
	real, allocatable, save :: utrans(:)	!transport in x on element
	real, allocatable, save :: vtrans(:)	!transport in y on element
	real, allocatable, save :: uvover(:)	!scalar for overlay
	real, allocatable, save :: wsnv(:)

	real, allocatable, save :: ureg(:,:)
	real, allocatable, save :: vreg(:,:)

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
          deallocate(uvnode)
          deallocate(vvnode)
          deallocate(uvelem)
          deallocate(vvelem)
          deallocate(utrans)
          deallocate(vtrans)
          deallocate(uvover)
          deallocate(wsnv)
        end if

        nel_hydro_plot = nel
        nkn_hydro_plot = nkn

	if( nkn == 0 ) return

        allocate(uvnode(nkn))
        allocate(vvnode(nkn))
        allocate(uvelem(nel))
        allocate(vvelem(nel))
        allocate(utrans(nel))
        allocate(vtrans(nel))
        allocate(uvover(nkn))
        allocate(wsnv(nkn))

	end subroutine mod_hydro_plot_init

!******************************************************************

	subroutine hydro_plot_regular(nx,ny)

	integer nx,ny

        if( nx == nx_hydro_plot .and. ny == ny_hydro_plot ) return

        if( nx > 0 .or. ny > 0 ) then
          if( nx == 0 .or. ny == 0 ) then
            write(6,*) 'nx,ny: ',nx,ny
	    stop 'error stop mod_hydro_plot_regular: incompatible params'
          end if
        end if

	if( nxy_limit > 0 .and. max(nx,ny) > nxy_limit ) then
	  write(6,*) 'nx,ny,nxy_limit: ',nx,ny,nxy_limit
	  write(6,*) 'limiting size of regular grid'
	  nx = min(nx,nxy_limit)
	  ny = min(ny,nxy_limit)
	end if

        if( nx_hydro_plot > 0 ) then
          deallocate(ureg)
          deallocate(vreg)
	end if

        nx_hydro_plot = nx
        ny_hydro_plot = ny

	if( nx == 0 ) return

        allocate(ureg(nx,ny))
        allocate(vreg(nx,ny))
	
	end subroutine hydro_plot_regular

!==================================================================
	end module mod_hydro_plot
!==================================================================

