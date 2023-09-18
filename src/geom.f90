
!==================================================================
	module geom
!==================================================================

	implicit none

	integer, private, save  :: nkn_geom = 0
	integer, private, save  :: nel_geom = 0
	integer, private, save  :: ngr_geom = 0
	integer, private, save  :: nlk_geom = 0

	integer, save :: maxlnk = 0

	integer, allocatable, save :: ilinkv(:)
	integer, allocatable, save :: lenkv(:)
	integer, allocatable, save :: lenkiiv(:)
	integer, allocatable, save :: linkv(:)

	integer, allocatable, save :: ieltv(:,:)
	integer, allocatable, save :: kantv(:,:)
	double precision, allocatable, save :: dxv(:)
	double precision, allocatable, save :: dyv(:)
        integer, allocatable, save :: auxv_iei(:,:)

!==================================================================
	contains
!==================================================================

	subroutine mod_geom_init(nkn,nel,ngr,nknex)

	integer nkn,nel,ngr
        integer, optional :: nknex

	integer nlk

        if( ngr == ngr_geom .and. nel == nel_geom .and. &
     &      nkn == nkn_geom ) return

        if( ngr > 0 .or. nel > 0 .or. nkn > 0 ) then
          if( ngr == 0 .or. nel == 0 .or. nkn == 0 ) then
            write(6,*) 'ngr,nel,nkn: ',ngr,nel,nkn
            stop 'error stop mod_geom_init: incompatible parameters'
          end if
        end if

        if( nkn_geom > 0 ) then
          deallocate(ilinkv)
          deallocate(lenkv)
          deallocate(lenkiiv)
          deallocate(linkv)
          deallocate(ieltv)
          deallocate(kantv)
          deallocate(dxv)
          deallocate(dyv)
          deallocate(auxv_iei)
        end if

	nlk = 3*nel + 2*nkn
	maxlnk = ngr

        ngr_geom = ngr
        nel_geom = nel
        nkn_geom = nkn
        nlk_geom = nlk

	if( nkn == 0 ) return

        if(present(nknex)) then
          allocate(kantv(2,nknex))
        else
          allocate(kantv(2,nkn))
        end if

        allocate(ilinkv(nkn+1))
        allocate(lenkv(nlk))
        allocate(lenkiiv(nlk))
        allocate(linkv(nlk))
        allocate(ieltv(3,nel))
        allocate(dxv(nkn))
        allocate(dyv(nkn))
        allocate(auxv_iei(3,nel))

	end subroutine mod_geom_init

!==================================================================
	end module geom
!==================================================================

