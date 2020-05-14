
!--------------------------------------------------------------------------
!
!    Copyright (C) 2019  Marco Bajo
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

! Paralution solver routines
!
! revision log :
!
! 18.10.2019	mbj	first version
!
!*************************************************************************

	subroutine para_initialize_system(matrix)

! Initialize system - to be called only once

	use mod_system
	use mod_system_interface

	implicit none

	type(smatrix) :: matrix

	integer n2max,n3max
	integer n2zero,n3zero

	  call coo_init_new(matrix) !construct pointers for coo matrix format

	  n2max = matrix%n2max
	  n3max = matrix%n3max
	  n2zero = matrix%n2zero
	  n3zero = matrix%n3zero

	  if( n2zero > n2max ) then
	    stop 'error stop para_initialize_system: non zero 2d max'
	  end if
	  if( n3zero > n3max ) then
	    stop 'error stop para_initialize_system: non zero 3d max'
	  end if

	  write(6,*) 'SOLVER: Paralution'
          print*, 'coo-matrix initialisation...'
          print*, 'Number of non-zeros 2d: ',n2zero,n2max
          print*, 'Number of non-zeros 3d: ',n3zero,n3max

	end

!*************************************************************************

	subroutine para_init_system(matrix)

! Initialize vector and matrix      

	use mod_system
	use, intrinsic :: ISO_C_BINDING, only : C_INT, C_PTR, C_DOUBLE, C_CHAR, C_NULL_CHAR, C_LOC

	implicit none

	type(smatrix) :: matrix

	integer, save :: icall = 0

	interface
         subroutine paralution_init() BIND(C)
           use, intrinsic :: ISO_C_BINDING, only : C_INT, C_PTR, C_DOUBLE, C_CHAR
         end subroutine
        end interface

	if (icall == 0) then
           call paralution_init()
	   icall = 1
        end if

	matrix%rvec2d = 0.
	matrix%raux2d = 0.
	matrix%rvec3d = 0.
	matrix%raux3d = 0.
	matrix%c2coo = 0.
	matrix%c3coo = 0.

      end

!*************************************************************************

      subroutine para_solve_system_csr(matrix,buse3d,nndim,nin,z)

! Solver routine with Paralution CUDA methods.

      use mod_system
      use, intrinsic :: ISO_C_BINDING, only : C_INT, C_PTR, C_DOUBLE, C_CHAR, C_NULL_CHAR, C_LOC

      implicit none

	type(smatrix), target :: matrix
	logical buse3d
	integer nndim		!dimension of non zeros in system
	integer nin		!dimension of system (unknowns)
	real z(nin)		!first guess

      real*8   guess(nin)

      integer ngl,nnzero

      real*8, allocatable :: raux(:)

      integer, allocatable :: iwork(:)

      type(smatrix), pointer :: mm

      interface
      subroutine paralution_fortran_solve_csr( n, m, nnzero, solver, mformat, preconditioner, pformat,&
       &                                     icsr, jcsr, csr, rvec, atol, rtol, div, maxiter, basis,&
       &                                     p, q, x, iter, resnorm, ierr ) BIND(C)

      use, intrinsic :: ISO_C_BINDING, only : C_INT, C_PTR, C_DOUBLE, C_CHAR

      integer(kind=C_INT), value, intent(in)  :: n, m, nnzero, maxiter, basis, p, q
      real(kind=C_DOUBLE), value, intent(in)  :: atol, rtol, div
      integer(kind=C_INT),        intent(out) :: iter, ierr
      real(kind=C_DOUBLE),        intent(out) :: resnorm
      type(C_PTR),         value, intent(in)  :: icsr, jcsr, csr, rvec
      type(C_PTR),         value              :: x
      character(kind=C_CHAR)                  :: solver, mformat, preconditioner, pformat

      end subroutine paralution_fortran_solve_csr
      end interface

      integer(kind=C_INT)   :: n, m, iter, ierr
      real(kind=C_DOUBLE)   :: resnorm

      integer(kind=C_INT), allocatable, target :: icsr(:), jcsr(:)
      real(kind=C_DOUBLE), allocatable, target :: csr(:)
      real(kind=C_DOUBLE), allocatable, target :: rvec(:), x(:)

      !Paralution parameters
      character(len=80) :: Solver,Op_mat_form,Prec,Prec_mat_form

      mm => matrix

      n = nin
      m = n  !n row and col are the same

      allocate(rvec(n),raux(n))
      allocate( icsr(n+1), jcsr(nndim) )
      allocate( csr(nndim) )
      allocate(iwork(max(n+1,2*nndim)))


      rvec = 0.
      raux = 0.

      ngl = n

!--------------------------------------------------
! CONVERSION AND SORTING
!--------------------------------------------------

      if( buse3d ) then
        nnzero = mm%n3zero
        call coocsr(ngl,nnzero,mm%c3coo,mm%i3coo,mm%j3coo,csr,jcsr,icsr)
        !write(6,*)'3D nnzero',nnzero,ngl
	!call coo_show(ngl,nnzero,i3coo,j3coo,c3coo)
	!call csr_show(ngl,nnzero,icsr,jcsr,csr)
	!call coo_print(ngl,nnzero,i3coo,j3coo,c3coo,rvec3d)
      else
        nnzero = mm%n2zero
        call coocsr(ngl,nnzero,mm%c2coo,mm%i2coo,mm%j2coo,csr,jcsr,icsr)
        !write(6,*)'2D nnzero',nnzero,ngl
	!call coo_show(ngl,nnzero,i2coo,j2coo,c2coo)
	!call csr_show(ngl,nnzero,icsr,jcsr,csr)
	!call coo_print(ngl,nnzero,i2coo,j2coo,c2coo,rvec)
      end if
      
      if( nnzero .gt. nndim .or. ngl+1 .gt. 2*nndim ) goto 99

      !call csort (ngl,csr,jcsr,icsr,iwork,.true.)
      call csort (ngl,csr,jcsr,icsr,.true.)	!new calling modus spk2

!--------------------------------------------------

      if( buse3d ) then
	rvec = mm%rvec3d
	raux = mm%raux3d
      else
	rvec = mm%rvec2d
	raux = mm%raux2d
      end if

      guess = z
      !guess = 0.
      x = guess

      ! Preconditioner matrix format (DENSE,CSR,MCSR,COO,DIA,ELL,HYB)
      Prec_mat_form = 'MCSR'
      ! Operation matrix format(DENSE,CSR,MCSR,COO,DIA,ELL,HYB)
      Op_mat_form = 'MCSR'   !fastest

      ! Preconditioner (None,Jacobi,MultiColoredGS,MultiColoredSGS,ILU,MultiColoredILU)
      Prec = 'Jacobi'
      !Prec = 'MultiColoredILU'      !much slower with big matrices
      ! Solver (CG,BiCGStab,GMRES,Fixed-Point)
      Solver = 'BiCGStab'

      ! Run paralution C function for CSR matrices
      call paralution_fortran_solve_csr( n, m, nnzero,                                       &
      &                                  trim(Solver) // C_NULL_CHAR,                        &
      &                                  trim(Op_mat_form) // C_NULL_CHAR,                   &
      &                                  trim(Prec) // C_NULL_CHAR,                          &
      &                                  trim(Prec_mat_form) // C_NULL_CHAR,                 &
      &                                  C_LOC(icsr),C_LOC(jcsr),C_LOC(csr), C_LOC(rvec),    &
      &                                  1e-8_C_DOUBLE, 1e-8_C_DOUBLE, 1e+8_C_DOUBLE, 2000,  &
      &                                  30, 0, 1, C_LOC(x), iter, resnorm, ierr )

      raux = x

!-----------------------------------------------------------------

      if( buse3d ) then
         mm%rvec3d(1:ngl) = raux(1:ngl)
      else
         mm%rvec2d(1:ngl) = raux(1:ngl)
      endif

      deallocate(csr,icsr,jcsr,iwork)
      deallocate(raux,rvec)

      return
 99   continue
      write(6,*) nnzero,ngl+1,nndim
      stop 'error stop para_solve_system: dimension iwork'
      end

!*************************************************************************

