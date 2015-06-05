
	module hydro_vel

	implicit none

        !real ulov(nlvdim,neldim)
        !common /ulov/ulov
        !real ulnv(nlvdim,neldim)
        !common /ulnv/ulnv
        !real vlov(nlvdim,neldim)
        !common /vlov/vlov
        !real vlnv(nlvdim,neldim)
        !common /vlnv/vlnv
        !real wlov(0:nlvdim,nkndim)
        !common /wlov/wlov
        !real wlnv(0:nlvdim,nkndim)
        !common /wlnv/wlnv
        !save /ulov/,/ulnv/,/vlov/,/vlnv/,/wlov/,/wlnv/

	integer, private, save :: nkn_hydro_vel = 0
	integer, private, save :: nel_hydro_vel = 0
	integer, private, save :: nlv_hydro_vel = 0

	real, allocatable, save :: ulov(:,:)
	real, allocatable, save :: ulnv(:,:)
	real, allocatable, save :: vlov(:,:)
	real, allocatable, save :: vlnv(:,:)
	real, allocatable, save :: wlov(:,:)
	real, allocatable, save :: wlnv(:,:)

	contains

!************************************************************

        subroutine mod_hydro_vel_init(nkn,nel,nlv)

        integer nkn, nel, nlv

        if( nkn == nkn_hydro_vel .and. nel == nel_hydro_vel .and.
     +      nlv == nlv_hydro_vel ) return

        if( nkn_hydro_vel > 0 .or. nel_hydro_vel > 0 .or.
     +      nlv_hydro_vel > 0 ) then
          deallocate(ulov)
          deallocate(ulnv)
          deallocate(vlov)
          deallocate(vlnv)

          deallocate(wlov)
          deallocate(wlnv)
        end if

        if( nkn == 0 ) return
        if( nel == 0 ) return
        if( nlv == 0 ) return

        nkn_hydro_vel = nkn
        nel_hydro_vel = nel
        nlv_hydro_vel = nlv

        allocate(ulov(nlv,nel))
        allocate(ulnv(nlv,nel))

        allocate(vlov(nlv,nel))
        allocate(vlnv(nlv,nel))

        allocate(wlov(0:nlv,nkn))
        allocate(wlnv(0:nlv,nkn))

        end subroutine mod_hydro_vel_init

!************************************************************

        end module hydro_vel
