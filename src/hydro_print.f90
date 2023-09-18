
	module hydro_print

	implicit none

        integer, private, save :: nkn_hydro_print = 0
        integer, private, save :: nlv_hydro_print = 0

        double precision, allocatable, save :: uprv(:,:)
        double precision, allocatable, save :: vprv(:,:)

        double precision, allocatable, save :: q2Dprv(:,:)
        double precision, allocatable, save :: c2Dprv(:,:)
        double precision, allocatable, save :: upro(:,:)
        double precision, allocatable, save :: vpro(:,:)
        double precision, allocatable, save :: wprv(:,:)
        double precision, allocatable, save :: up0v(:), vp0v(:)
        double precision, allocatable, save :: xv(:,:)

!************************************************************

	contains

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

          deallocate(q2Dprv)
          deallocate(c2Dprv)
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

        allocate(q2Dprv(nlv,nkn))
        allocate(c2Dprv(nlv,nkn))
        allocate(upro(nlv,nkn))
        allocate(vpro(nlv,nkn))
        allocate(wprv(0:nlv,nkn))

	allocate(up0v(nkn))
	allocate(vp0v(nkn))

	allocate(xv(3,nkn))

        end subroutine mod_hydro_print_init

!******************************************************************

	subroutine getuv(l,k,u,v)

! accessor routine to get velocities u/v

	implicit none

        integer l,k
        double precision u,v

        u = uprv(l,k)
        v = vprv(l,k)

        end

!************************************************************

        end module hydro_print
