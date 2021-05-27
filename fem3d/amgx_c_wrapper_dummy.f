!--------------------------------------------------------------------------
!
!    Copyright (C) 2021 Celia Laurent
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
!
! dummy functions for enabling PETSc without AmgX
!
! revision log :
!
! 20.04.2021	clr	adding this dummy file to remove use_AmgX pragma directive
!
!*******************************************************************
      subroutine CAmgX_GetInitSolver(AmgX_Solver,PETSC_COMM,
     +                               modestr,cfgfile)
       use iso_c_binding
        implicit none
           type(c_ptr) :: AmgX_Solver
           integer :: PETSC_COMM  ! MPI Communicator
           character(len=*,kind=c_char) :: modestr
           character(len=*,kind=c_char) :: cfgfile
           stop "ERROR CAmgX_dummy, AmgX wasn't enabled at compilation"
      end subroutine CAmgX_GetInitSolver
!*******************************************************************
      subroutine CAmgX_GetSolver(AmgX_Solver)
       use iso_c_binding
        implicit none
           type(c_ptr) :: AmgX_Solver
           stop "ERROR CAmgX_dummy, AmgX wasn't enabled at compilation"
      end subroutine CAmgX_GetSolver
!*******************************************************************
      subroutine CAmgX_Initialize(AmgX_Solver,PETSC_COMM,
     +                            modestr,cfgfile,perr)
#include "petsc/finclude/petsc.h"
       use iso_c_binding
        use petscdm
        use petscdmlabel
        implicit none
           type(c_ptr) :: AmgX_Solver
           integer :: PETSC_COMM  ! MPI Communicator
           character(len=*,kind=c_char) :: modestr
           character(len=*,kind=c_char) :: cfgfile
           PetscErrorCode :: perr     ! error flag

           stop "ERROR CAmgX_dummy, AmgX wasn't enabled at compilation"
      end subroutine CAmgX_Initialize
!*******************************************************************
      subroutine CAmgX_SetA(AmgX_Solver,A,perr)
#include "petsc/finclude/petsc.h"
       use iso_c_binding
        use petscdm
        use petscdmlabel
        use petscmat
        implicit none
           type(c_ptr) :: AmgX_Solver
           Mat  :: A       ! Matrix A
           PetscErrorCode :: perr     ! error flag
           stop "ERROR CAmgX_dummy, AmgX wasn't enabled at compilation"
      end subroutine CAmgX_SetA
!*******************************************************************
      subroutine CAmgX_Solve(AmgX_Solver,X,B,perr)
#include "petsc/finclude/petsc.h"
       use iso_c_binding
        use petscdm
        use petscdmlabel
        use petscvec
        implicit none
           type(c_ptr) :: AmgX_Solver
           Vec  :: X       ! solution
           Vec  :: B       ! rhs
           PetscErrorCode :: perr     ! error flag
           stop "ERROR CAmgX_dummy, AmgX wasn't enabled at compilation"
      end subroutine CAmgX_Solve
!*******************************************************************
      subroutine CAmgX_Finalize(AmgX_Solver,perr)
#include "petsc/finclude/petsc.h"
       use iso_c_binding
        use petscdm
        use petscdmlabel
        implicit none
           type(c_ptr) :: AmgX_Solver
           PetscErrorCode :: perr     ! error flag
           stop "ERROR CAmgX_dummy, AmgX wasn't enabled at compilation"
      end subroutine CAmgX_Finalize
!*******************************************************************
