
!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2018  Georg Umgiesser
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

! module for solving system matrix
!
! revision log :
!
! 15.12.2015    ggu&dbf adapted to 3d case (poisson)
! 23.04.2018    ggu	new framework - introduced matrix type
!
! notes :
!
! use 61 91 121 etc..  61 indicates 10E-6 of precision etc..
! best choice is 121 for iterative
! use 0 for direct solver - this should be a save choice if in doubt
!
! structure of calls:
!
! system_initialize						simsys_spk.f
!	mod_system_init						mod_system.f
!	mod_system_insert_elem_index				mod_system.f
!	spk_initialize_system					subspk.f
! system_init							simsys_spk.f
!	spk_init_system						subspk.f
!		coo_init_new					subcoo.f
! system_solve							simsys_spk.f
!	spk_solve_system					subspk.f
!		coocsr						subcoo.f
!		csort						subcoo.f
!		ilu0,ilut,ilutp,ilud,iludp			spk_ilut.f
!		runrc						spk_itaux.f
! system_assemble						simsys_spk.f
!	-> c2coo, rvec2d					(mod_system.f)
! system_add_rhs						simsys_spk.f
!	-> rvec2d
! system_get						simsys_spk.f
!	-> rvec2d
!
! system_adjust_matrix_3d					simsys_spk.f
!	coo_adjust_3d						subcoo.f
!
!==================================================================
        module mod_system
!==================================================================

        implicit none

        type :: smatrix

	logical :: bglobal = .false.		!global matrix mpi

        integer :: nkn_system = 0
        integer :: nel_system = 0
        integer :: ngr_system = 0
        integer :: mbw_system = 0
        integer :: nlv_system = 0
        integer :: nkn_amat = 0
        integer :: mbw_amat = 0

	integer :: iprec = 0
	integer :: nthpard = 8	! Number of threads for Pardiso solver

	integer :: n2zero = 0
	integer :: n3zero = 0
	integer :: n2max = 0
	integer :: n3max = 0

	integer, allocatable :: nen3v_system(:,:)

        double precision, allocatable :: vs1v(:)
        double precision, allocatable :: vs3v(:)
        integer, allocatable :: is2v(:)

        double precision, allocatable :: rvec2d(:)
        double precision, allocatable :: raux2d(:)
        double precision, allocatable :: rvec3d(:)
        double precision, allocatable :: raux3d(:)

        !double precision, allocatable :: coo(:)
        !integer, allocatable :: icoo(:)
        !integer, allocatable :: jcoo(:)
        !integer, allocatable :: ijp(:)

        integer, allocatable :: ng(:)		!grade of node k
        integer, allocatable :: ntg(:)		!cumulative grade inc k
        integer, allocatable :: n3g(:)		!grade of node k
        integer, allocatable :: nt3g(:)		!cumulative grade inc k
        integer, allocatable :: diag(:)		!relativ pos of diag el
        integer, allocatable :: iorder(:,:)	!all nodes in k
        integer, allocatable :: ijp_ie(:,:,:)	!pointer to entry

        integer, allocatable :: i2coo(:)
        integer, allocatable :: j2coo(:)
        double precision, allocatable :: c2coo(:)

        integer, allocatable :: i3coo(:)
        integer, allocatable :: j3coo(:)
        double precision, allocatable :: c3coo(:)

        integer, allocatable :: back3coo(:,:)

        double precision, allocatable :: amat(:)

	end type smatrix

	logical, save :: bsys3d = .false.
	logical, save :: bsysexpl = .false.

        double precision, parameter :: d_tiny = tiny(1.d+0)
        double precision, parameter :: r_tiny = tiny(1.)

	type(smatrix), save, target  :: l_matrix	!local matrix
	type(smatrix), save, target  :: g_matrix	!global matrix

!==================================================================
	contains
!==================================================================

        subroutine mod_system_init(nkn,nel,ngr,mbw,nlv,matrix)

        integer  :: nkn
        integer  :: nel
        integer  :: ngr
        integer  :: mbw
        integer  :: nlv
	type(smatrix) :: matrix

	integer ngl2d
	integer ngl3d
	integer n2max
	integer n3max

        if( mbw == matrix%mbw_system .and. 
     +		nel == matrix%nel_system .and.
     +		nkn == matrix%nkn_system .and. 
     +		ngr == matrix%ngr_system ) return

        if( mbw > 0 .or. nel > 0 .or. nkn > 0 .or. ngr > 0) then
          if( mbw == 0 .or. nel == 0 .or. nkn == 0 .or. ngr == 0) then
            write(6,*) 'mbw,ngr,nel,nkn: ',mbw,ngr,nel,nkn
            stop 'error stop mod_system_init: incompatible parameters'
          end if
        end if

        if( matrix%nkn_system > 0 ) then
          deallocate(matrix%nen3v_system)

          deallocate(matrix%vs1v)		!next three only used in lp
          deallocate(matrix%vs3v)
          deallocate(matrix%is2v)

          deallocate(matrix%rvec2d)
          deallocate(matrix%raux2d)
          deallocate(matrix%rvec3d)
          deallocate(matrix%raux3d)

	  deallocate(matrix%ng)
	  deallocate(matrix%ntg)
	  deallocate(matrix%n3g)
	  deallocate(matrix%nt3g)
	  deallocate(matrix%diag)
	  deallocate(matrix%iorder)		!prob not needed
	  deallocate(matrix%ijp_ie)

          deallocate(matrix%i2coo)
          deallocate(matrix%j2coo)
          deallocate(matrix%c2coo)

          deallocate(matrix%i3coo)
          deallocate(matrix%j3coo)
          deallocate(matrix%c3coo)
          deallocate(matrix%back3coo)
        end if

        matrix%nkn_system = nkn
        matrix%nel_system = nel
        matrix%ngr_system = ngr
        matrix%mbw_system = mbw
        matrix%nlv_system = nlv

	n2max = 7*nkn
	n3max = 6*nkn*nlv + nkn*(2+3*nlv)
	if( .not. bsys3d ) n3max = 1

	matrix%n2max = n2max
	matrix%n3max = n3max

	ngl2d = nkn
	ngl3d = nkn*(nlv+2)
	if( .not. bsys3d ) ngl3d = 1

        if( nkn == 0 ) return

        allocate(matrix%nen3v_system(3,nel))

        allocate(matrix%vs1v(nkn))
        allocate(matrix%vs3v(nkn))
        allocate(matrix%is2v(nkn))

        allocate(matrix%rvec2d(ngl2d))
        allocate(matrix%raux2d(ngl2d))
        allocate(matrix%rvec3d(ngl3d))
        allocate(matrix%raux3d(ngl3d))

        allocate(matrix%ng(nkn))
        allocate(matrix%ntg(0:nkn))
        allocate(matrix%n3g(nkn))
        allocate(matrix%nt3g(0:nkn))
        allocate(matrix%diag(nkn))
        allocate(matrix%iorder(ngr+1,nkn))
        allocate(matrix%ijp_ie(3,3,nel))

        allocate(matrix%i2coo(n2max))
        allocate(matrix%j2coo(n2max))
        allocate(matrix%c2coo(n2max))

        allocate(matrix%i3coo(n3max))
        allocate(matrix%j3coo(n3max))
        allocate(matrix%c3coo(n3max))
        allocate(matrix%back3coo(4,n3max))

        end subroutine mod_system_init

!****************************************************************

        subroutine mod_system_amat_init(nkn,mbw,matrix)

        integer  :: nkn
        integer  :: mbw
	type(smatrix) :: matrix

        integer  :: mat

        if( mbw == matrix%mbw_amat .and. nkn == matrix%nkn_amat ) return

        if( mbw > 0 .or. nkn > 0 ) then
          if( mbw == 0 .or. nkn == 0 ) then
            write(6,*) 'mbw,nkn: ',mbw,nkn
            stop 'error stop mod_system_amat_init: incompatible params'
          end if
        end if

        if( matrix%nkn_amat > 0 ) then
          deallocate(matrix%amat)
	end if

        matrix%nkn_amat = nkn
        matrix%mbw_amat = mbw

	mat = nkn*(1+3*mbw)

	if( nkn == 0 ) return

        allocate(matrix%amat(mat))

	end subroutine mod_system_amat_init

!****************************************************************

        subroutine mod_system_insert_elem_index(nel,nen3v,matrix)

        integer  :: nel
        integer  :: nen3v(3,nel)
	type(smatrix) :: matrix

	if( nel /= matrix%nel_system ) then
	  write(6,*) 'nel,nel_system: ',nel,matrix%nel_system
	  stop 'error stop mod_system_insert_elem_index: dimensions'
	end if

	matrix%nen3v_system = nen3v

        end subroutine mod_system_insert_elem_index

!****************************************************************

        subroutine mod_system_set_local(matrix)

	type(smatrix) :: matrix
	matrix%bglobal = .false.

        end subroutine mod_system_set_local

!****************************************************************

        subroutine mod_system_set_global(matrix)

	type(smatrix) :: matrix
	matrix%bglobal = .true.

        end subroutine mod_system_set_global

!==================================================================
        end module mod_system
!==================================================================

!==================================================================
        module mod_system_interface
!==================================================================

! this is only needed for older compilers - but it does not hurt

        INTERFACE

        subroutine spk_solve_system(matrix,buse3d,nndim,n,z)
        use mod_system
        type(smatrix), target :: matrix
        logical buse3d
        integer nndim           !dimension of non zeros in system
        integer n               !dimension of system (unknowns)
        real z(n)               !first guess
        end

        subroutine coo_adjust_3d(matrix)
        use mod_system
        type(smatrix), target :: matrix
        end

        subroutine coo_init_new(matrix)
        use mod_system
        type(smatrix), target :: matrix
        end

        END INTERFACE

!==================================================================
        end module mod_system_interface
!==================================================================

