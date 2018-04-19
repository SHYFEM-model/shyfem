
!==================================================================
	module mod_hydro_plot
!==================================================================

	implicit none

	integer, private, save  :: nkn_hydro_plot = 0
	integer, private, save  :: nel_hydro_plot = 0
	integer, private, save  :: nlv_hydro_plot = 0
	integer, private, save  :: np_hydro_plot = 0
	integer, private, save  :: nx_hydro_plot = 0
	integer, private, save  :: ny_hydro_plot = 0
	integer, private, save  :: nxy_limit = 1000	!0 for no limit

	logical, save :: bminmax = .false.
	logical, save :: berrintp = .false.
	logical, save :: bonelem
	logical, save :: bisreg
	logical, save :: bistrans

	real, save :: hydro_regpar(7) = 0.

	real, allocatable, save :: uvnode(:)	!variable in x on node
	real, allocatable, save :: vvnode(:)	!variable in y on node
	real, allocatable, save :: uvelem(:)	!variable in x on element
	real, allocatable, save :: vvelem(:)	!variable in y on elemen
	real, allocatable, save :: utrans(:)	!transport in x on element
	real, allocatable, save :: vtrans(:)	!transport in y on element
	real, allocatable, save :: uvspeed(:)	!speed on node (aux)
	real, allocatable, save :: uvdir(:)	!direction on node (aux)
	real, allocatable, save :: uvover(:)	!scalar for overlay
	real, allocatable, save :: wsnv(:)

	real, allocatable, save :: ureg(:,:)
	real, allocatable, save :: vreg(:,:)

        real, allocatable, save :: arfvlv(:)
        real, allocatable, save :: hetv(:)
        real, allocatable, save :: parray(:)

        logical, allocatable, save :: bwater(:)
        logical, allocatable, save :: bkwater(:)
        logical, allocatable, save :: bplot(:)
        logical, allocatable, save :: bkplot(:)

        real, allocatable, save :: fvlv(:,:)
        real, allocatable, save :: wauxv(:,:)
        real, allocatable, save :: het3v(:,:)
        real, allocatable, save :: p3(:,:)

!==================================================================
	contains
!==================================================================

	subroutine mod_hydro_plot_init(nkn,nel,nlv,np)

	integer nkn,nel,nlv,np

        if( nkn == nkn_hydro_plot .and. nel == nel_hydro_plot .and.
     +          nlv == nlv_hydro_plot .and. np == np_hydro_plot ) return

        if( nel > 0 .or. nkn > 0 .or. nlv > 0 .or. np > 0 ) then
          if( nel == 0 .or. nkn == 0 . or. nlv == 0 .or. np == 0 ) then
            write(6,*) 'nel,nkn,nlv,np: ',nel,nkn,nlv,np
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
          deallocate(uvspeed)
          deallocate(uvdir)
          deallocate(uvover)
          deallocate(wsnv)

          deallocate(arfvlv)
          deallocate(hetv)
          deallocate(parray)
          deallocate(bwater)
          deallocate(bkwater)
          deallocate(bplot)
          deallocate(bkplot)

          deallocate(fvlv)
          deallocate(wauxv)
          deallocate(het3v)
          deallocate(p3)
        end if

        nel_hydro_plot = nel
        nkn_hydro_plot = nkn
        nlv_hydro_plot = nlv
        np_hydro_plot = np

	if( nkn == 0 ) return

        allocate(uvnode(nkn))
        allocate(vvnode(nkn))
        allocate(uvelem(nel))
        allocate(vvelem(nel))
        allocate(utrans(nel))
        allocate(vtrans(nel))
        allocate(uvspeed(nkn))
        allocate(uvdir(nkn))
        allocate(uvover(nkn))
        allocate(wsnv(nkn))

        allocate(arfvlv(nkn))
        allocate(hetv(nel))
        allocate(parray(np))
        allocate(bwater(nel))
        allocate(bkwater(nkn))
        allocate(bplot(nel))
        allocate(bkplot(nkn))

        allocate(fvlv(nlv,nkn))
        allocate(wauxv(0:nlv,nkn))
        allocate(het3v(nlv,nel))
        allocate(p3(nlv,np))

	end subroutine mod_hydro_plot_init

!******************************************************************

	subroutine mod_hydro_plot_regular_init(nx,ny)

	integer nx,ny

        if( nx == nx_hydro_plot .and. ny == ny_hydro_plot ) return

        if( nx > 0 .or. ny > 0 ) then
          if( nx == 0 .or. ny == 0 ) then
            write(6,*) 'nx,ny: ',nx,ny
	    stop 'error stop mod_hydro_plot_regular_init: ' //
     +				'incompatible params'
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
	
	end subroutine mod_hydro_plot_regular_init

!******************************************************************

	subroutine mod_hydro_set_regpar(regpar)

	real regpar(7)

	hydro_regpar = regpar

	end subroutine mod_hydro_set_regpar

!******************************************************************

	subroutine mod_hydro_get_regpar(regpar)

	real regpar(7)

	regpar = hydro_regpar

	end subroutine mod_hydro_get_regpar

!==================================================================
	end module mod_hydro_plot
!==================================================================

