	module mod_hydro_print

	implicit none

        !real uprv(nlvdim,nkndim)
        !common /uprv/uprv
        !real vprv(nlvdim,nkndim)
        !common /vprv/vprv
        !real upro(nlvdim,nkndim)
        !common /upro/upro
        !real vpro(nlvdim,nkndim)
        !common /vpro/vpro
        !real wprv(0:nlvdim,nkndim)
        !common /wprv/wprv
        !save /uprv/,/vprv/,/upro/,/vpro/,/wprv/
        !real up0v(nkndim)
        !common /up0v/up0v
        !real vp0v(nkndim)
        !common /vp0v/vp0v
        !save /up0v/,/vp0v/
        !real xv(3,nkndim)
        !common /xv/xv
        !save /xv/

        integer, private, save :: nkn_hydro_print = 0
        integer, private, save :: nlv_hydro_print = 0

        real, allocatable, save :: uprv(:,:)
        real, allocatable, save :: vprv(:,:)
        real, allocatable, save :: upro(:,:)
        real, allocatable, save :: vpro(:,:)
        real, allocatable, save :: wprv(:,:)
        real, allocatable, save :: up0v(:), vp0v(:)
        real, allocatable, save :: xv(3,:)

!************************************************************

        subroutine mod_hydro_print_init(nkn,nlv)

        integer nkn, nlv

        if( nkn == nkn_hydro_print .and. 
     +      nlv == nlv_hydro_print ) return

        if( nkn_hydro_print > 0 .or.
     +      nlv_hydro_print > 0 ) then
          deallocate(uprv)
          deallocate(vprv)
          deallocate(upro)
          deallocate(vpro)
          deallocate(wprv)

          deallocate(up0v)
          deallocate(vp0v)

          deallocate(xv)
        end if

        if( nkn == 0 ) return
        if( nlv == 0 ) return

        nkn_hydro_print = nkn
        nlv_hydro_print = nlv

        allocate(uprv(nlv,nkn))
        allocate(vprv(nlv,nkn))
        allocate(upro(nlv,nkn))
        allocate(vpro(nlv,nkn))
        allocate(wprv(0:nlv,nkn))

	allocate(up0v(nkn))
	allocate(vp0v(nkn))

	allocate(xv(3,nkn))

        end subroutine mod_hydro_print_init

!************************************************************

        end module mod_hydro_print
