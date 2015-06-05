	module mod_levels

	implicit none

        !integer ilhkv(nkndim)
        !common /ilhkv/ilhkv
        !integer ilhv(neldim)
        !common /ilhv/ilhv
        !real hlv(nlvdim), hldv(nlvdim)
        !common /hlv/hlv, /hldv/hldv
        !integer ilmv(neldim)
        !common /ilmv/ilmv
        !integer ilmkv(nkndim)
        !common /ilmkv/ilmkv
        !save /ilhkv/,/ilhv/,/hlv/,/hldv/,/ilmv/,/ilmkv/

	integer, private, save :: nkn_levels = 0
	integer, private, save :: nel_levels = 0
	integer, private, save :: nlv_levels = 0

	integer, allocatable, save :: ilhkv(:)
	integer, allocatable, save :: ilhv(:)
	real, allocatable, save :: hlv(:), hldv(:)
	integer, allocatable, save :: ilmv(:)
	integer, allocatable, save :: ilmkv(:)

	contains

!************************************************************

        subroutine mod_levels_init(nkn,nel,nlv)

        integer nkn, nel, nlv

        if( nkn == nkn_levels .and. nel == nel_levels .and.
     +      nlv == nlv_levels ) return

        if( nkn_levels > 0 .or. nel_levels > 0 .or.
     +      nlv_levels > 0 ) then
          deallocate(ilhkv)
          deallocate(ilmkv)

          deallocate(ilhv)
          deallocate(ilmv)

          deallocate(hlv)
          deallocate(hldv)
        end if

        if( nkn == 0 ) return
        if( nel == 0 ) return
        if( nlv == 0 ) return

        nkn_levels = nkn
        nel_levels = nel
        nlv_levels = nlv

        allocate(ilhkv(nkn))
        allocate(ilmkv(nkn))

        allocate(ilhv(nel))
        allocate(ilmv(nel))

        allocate(hlv(nlv))
        allocate(hldv(nlv))

        end subroutine mod_levels_init

!************************************************************

        end module mod_levels
