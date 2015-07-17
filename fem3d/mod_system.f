
!==================================================================
        module mod_system
!==================================================================

! use 61 91 121 etc..  61 indicates 10E-6 of precision etc..
! best choice is 121 for iterative
! use 0 for direct solver - this should be a save choice if in doubt

        implicit none

        integer, private, save :: nkn_system = 0
        integer, private, save :: nel_system = 0
        integer, private, save :: mbw_system = 0
        integer, private, save :: nkn_amat = 0
        integer, private, save :: mbw_amat = 0

	integer, save :: iprec = 0
	integer, save :: nnzero = 0
	integer, save :: csrdim = 0

        double precision, allocatable, save :: vs1v(:)
        double precision, allocatable, save :: vs2v(:)
        double precision, allocatable, save :: vs3v(:)
        integer, allocatable, save :: is2v(:)

        double precision, allocatable, save :: rvec(:)
        double precision, allocatable, save :: raux(:)

        double precision, allocatable, save :: coo(:)
        integer, allocatable, save :: icoo(:)
        integer, allocatable, save :: jcoo(:)
        integer, allocatable, save :: ijp(:)

        double precision, allocatable, save :: amat(:)

!==================================================================
	contains
!==================================================================

        subroutine mod_system_init(nkn,nel,mbw)

        integer  :: nkn
        integer  :: nel
        integer  :: mbw

        integer  :: csr,mat

        if( mbw == mbw_system .and. nel == nel_system .and.
     +      nkn == nkn_system ) return

        if( mbw > 0 .or. nel > 0 .or. nkn > 0 ) then
          if( mbw == 0 .or. nel == 0 .or. nkn == 0 ) then
            write(6,*) 'mbw,nel,nkn: ',mbw,nel,nkn
            stop 'error stop mod_system_init: incompatible parameters'
          end if
        end if

        if( nkn_system > 0 ) then
          deallocate(vs1v)
          deallocate(vs2v)
          deallocate(vs3v)
          deallocate(is2v)
          deallocate(rvec)
          deallocate(raux)
          deallocate(coo)
          deallocate(icoo)
          deallocate(jcoo)
          deallocate(ijp)
        end if

        nkn_system = nkn
        nel_system = nel
        mbw_system = mbw

	csr = 9 * nel
	csrdim = csr			!may be deleted later
	mat = nkn*(1+3*mbw)

        if( nkn == 0 ) return

        allocate(vs1v(nkn))
        allocate(vs2v(nkn))
        allocate(vs3v(nkn))
        allocate(is2v(nkn))
        allocate(rvec(nkn))
        allocate(raux(nkn))
        allocate(coo(csr))
        allocate(icoo(csr))
        allocate(jcoo(csr))
        allocate(ijp(mat))

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

