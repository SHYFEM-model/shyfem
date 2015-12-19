!
! module for solving system matrix
!
! revision log :
!
! 15.12.2015    ggu&deb     adapted to 3d case (poisson)
!
! notes :
!
! use 61 91 121 etc..  61 indicates 10E-6 of precision etc..
! best choice is 121 for iterative
! use 0 for direct solver - this should be a save choice if in doubt
!
!==================================================================
        module mod_system
!==================================================================

        implicit none

        integer, private, save :: nkn_system = 0
        integer, private, save :: nel_system = 0
        integer, private, save :: ngr_system = 0
        integer, private, save :: mbw_system = 0
        integer, private, save :: nlv_system = 0
        integer, private, save :: nkn_amat = 0
        integer, private, save :: mbw_amat = 0

	logical, save :: bsys3d = .false.

	integer, save :: iprec = 0
	integer, save :: nthpard = 8	! Number of threads for Pardiso solver

	integer, save :: n2zero = 0
	integer, save :: n3zero = 0
	integer, save :: n2max = 0
	integer, save :: n3max = 0

        double precision, allocatable, save :: vs1v(:)
        double precision, allocatable, save :: vs3v(:)
        integer, allocatable, save :: is2v(:)

        double precision, allocatable, save :: rvec2d(:)
        double precision, allocatable, save :: raux2d(:)
        double precision, allocatable, save :: rvec3d(:)
        double precision, allocatable, save :: raux3d(:)

        !double precision, allocatable, save :: coo(:)
        !integer, allocatable, save :: icoo(:)
        !integer, allocatable, save :: jcoo(:)
        !integer, allocatable, save :: ijp(:)

        integer, allocatable, save :: ng(:)		!grade of node k
        integer, allocatable, save :: ntg(:)		!cumulative grade inc k
        integer, allocatable, save :: n3g(:)		!grade of node k
        integer, allocatable, save :: nt3g(:)		!cumulative grade inc k
        integer, allocatable, save :: diag(:)		!relativ pos of diag el
        integer, allocatable, save :: iorder(:,:)	!all nodes in k
        integer, allocatable, save :: ijp_ie(:,:,:)	!pointer to entry

        integer, allocatable, save :: i2coo(:)
        integer, allocatable, save :: j2coo(:)
        double precision, allocatable, save :: c2coo(:)

        integer, allocatable, save :: i3coo(:)
        integer, allocatable, save :: j3coo(:)
        double precision, allocatable, save :: c3coo(:)

        integer, allocatable, save :: back3coo(:,:)

        double precision, allocatable, save :: amat(:)

        double precision, parameter :: d_tiny = tiny(1.d+0)
        double precision, parameter :: r_tiny = tiny(1.)

!==================================================================
	contains
!==================================================================

        subroutine mod_system_init(nkn,nel,ngr,mbw,nlv)

        integer  :: nkn
        integer  :: nel
        integer  :: ngr
        integer  :: mbw
        integer  :: nlv

	integer ngl2d
	integer ngl3d

        if( mbw == mbw_system .and. nel == nel_system .and.
     +      nkn == nkn_system .and. ngr == ngr_system ) return

        if( mbw > 0 .or. nel > 0 .or. nkn > 0 .or. ngr > 0) then
          if( mbw == 0 .or. nel == 0 .or. nkn == 0 .or. ngr == 0) then
            write(6,*) 'mbw,ngr,nel,nkn: ',mbw,ngr,nel,nkn
            stop 'error stop mod_system_init: incompatible parameters'
          end if
        end if

        if( nkn_system > 0 ) then
          deallocate(vs1v)			!next three only used in lp
          deallocate(vs3v)
          deallocate(is2v)

          deallocate(rvec2d)
          deallocate(raux2d)
          deallocate(rvec3d)
          deallocate(raux3d)

	  deallocate(ng,ntg,n3g,nt3g,diag)
	  deallocate(iorder)			!prob not needed
	  deallocate(ijp_ie)

          deallocate(i2coo)
          deallocate(j2coo)
          deallocate(c2coo)

          deallocate(i3coo)
          deallocate(j3coo)
          deallocate(c3coo)
          deallocate(back3coo)
        end if

        nkn_system = nkn
        nel_system = nel
        ngr_system = ngr
        mbw_system = mbw
        nlv_system = nlv

	n2max = 7*nkn
	n3max = 6*nkn*nlv + nkn*(2+3*nlv)
	if( .not. bsys3d ) n3max = 1

	ngl2d = nkn
	ngl3d = nkn*(nlv+2)
	if( .not. bsys3d ) ngl3d = 1

        if( nkn == 0 ) return

        allocate(vs1v(nkn))
        allocate(vs3v(nkn))
        allocate(is2v(nkn))

        allocate(rvec2d(ngl2d))
        allocate(raux2d(ngl2d))
        allocate(rvec3d(ngl3d))
        allocate(raux3d(ngl3d))

        allocate(ng(nkn))
        allocate(ntg(0:nkn))
        allocate(n3g(nkn))
        allocate(nt3g(0:nkn))
        allocate(diag(nkn))
        allocate(iorder(ngr+1,nkn))
        allocate(ijp_ie(3,3,nel))

        allocate(i2coo(n2max))
        allocate(j2coo(n2max))
        allocate(c2coo(n2max))

        allocate(i3coo(n3max))
        allocate(j3coo(n3max))
        allocate(c3coo(n3max))
        allocate(back3coo(4,n3max))

        end subroutine mod_system_init

c****************************************************************

        subroutine mod_system_amat_init(nkn,mbw)

        integer  :: nkn
        integer  :: mbw

        integer  :: mat

        if( mbw == mbw_amat .and. nkn == nkn_amat ) return

        if( mbw > 0 .or. nkn > 0 ) then
          if( mbw == 0 .or. nkn == 0 ) then
            write(6,*) 'mbw,nkn: ',mbw,nkn
            stop 'error stop mod_system_amat_init: incompatible params'
          end if
        end if

        if( nkn_amat > 0 ) then
          deallocate(amat)
	end if

        nkn_amat = nkn
        mbw_amat = mbw

	mat = nkn*(1+3*mbw)

	if( nkn == 0 ) return

        allocate(amat(mat))

	end subroutine mod_system_amat_init

!==================================================================
        end module mod_system
!==================================================================

