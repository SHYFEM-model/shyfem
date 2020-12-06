!--------------------------------------------------------------------------
!
!    Copyright (C) 2020 Celia Laurent
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
! module for solving system matrix using petsc
!
! revision log :
!   
! 21.12.2020	cla	original implementation
!
! notes :
!
! structure of calls:
!
!   shyfem
!   ├──system_initialize
!   |   ├──petsc_global_initialize
!   |   |  └──petsc_initialize
!   |   ├──petsc_global_create_setup
!   |   |  ├──petsc_create_indexes
!   |   |  └──pestc_identify_non_zeros_and_ghosts
!   |   ├──allocate(zeta_system)
!   |   ├──zeta_system=petsc_system(PETSc_zeta_config,AmgX_zeta_config)
!   |   ├──system_zeta%create_objects
!   |   |  ├──MatCreate, MatSet, ...
!   |   |  └──VecCreate,VecSet, VecMPISetGhost ...
!   |   └──petsc_global_close_setup
!   ├──loop on time
!   |   └──hydro
!   |      ├──system_init
!   |      |     └──reset_zero_entries
!   |      ├──hydro_zeta
!   |      |  ├──loop on elements
!   |      |  |  ├──zeta_system%mat3x3(:,:) = …
!   |      |  |  ├──zeta_system%vecx3(:) = …
!   |      |  |  └──zeta_system%add_matvec_values(element_number)
!   |      |  |      ├──MatSetValues 
!   |      |  |      └──VecSetValues 
!   |      |  └──zeta_system%add_full_rhs(length,values(:))
!   |      |     └──VecSetValue  
!   |      ├──system_solve
!   |      |  ├──zeta_system%matvec_assemble
!   |      |  |     └──Mat/VecAssemblyBegin/End   
!   |      |  ├──IF (ITER ==1) 
!   |      |  |     └──zeta_system%init_solver 
!   |      |  |         └──zeta_system%init_solver_PETSc
!   |      |  |             ├──KSPCreate, KSPSet, KSPGetPC, ...
!   |      |  |             └──PCSetType, PCFactorSetMatSolverType, ...
!   |      |  └──zeta_system%solve
!   |      |     ├──KSPSetOperators
!   |      |     └──KSPSolve
!   |      └──system_get
!   |         └──zeta_system%get_solution
!   |            └──VecGhostUpdateBegin/End
!   |                ├──VecGetArrayF90
!   |                └──VecRestoreArrayReadF90
!   └──system_finalize
!       ├──zeta_system%destroy_solver
!       |   └──zeta_system%destroy_solver_PETSc
!       |      └──KSPDestroy
!       ├──zeta_system%destroy_matvecobjects
!       |   ├──VecDestroy
!       |   └──MatDestroy
!       ├──deallocate(zeta_system)
!       ├──petsc_global_finalize
!       └──PetscFinalize
!
!==================================================================
        module mod_petsc
!==================================================================
        use mod_petsc_system
        use mod_petsc_global
        private
        
        type(petsc_system),public,pointer :: zeta_system
        character(len=80),public :: PETSc_zeta_configfile
        character(len=80),public :: AmgX_zeta_configfile

        integer :: petsc_iter=1

        public ::
     +           system_initialize,
     +           system_init,
     +           system_set_explicit,
     +           system_solve,
     +           system_solve_3d,
     +           system_get,
     +           system_get_3d,
     +           system_finalize      
        contains


!******************************************************************

      subroutine system_initialize

! allocates data structure

        implicit none
        !-------------------------------------------------------------
        ! Initialize Petsc 
        !-------------------------------------------------------------         
        call petsc_global_initialize
        call petsc_global_create_setup
         !----------------------------------------------------
         ! init the petsc system of matrix, vectors and solver
         !----------------------------------------------------
        allocate(zeta_system)
        zeta_system=petsc_system(PETSc_zeta_configfile,
     +                           AmgX_zeta_configfile
     +                            ) !constructor
        call zeta_system%create_objects  

        call petsc_global_close_setup ! free sparsity pattern and ghost arrays

      end subroutine system_initialize

!******************************************************************

	subroutine system_init

        call zeta_system%reset_zero_entries

	end subroutine system_init

!******************************************************************

	subroutine system_set_explicit

	use mod_system

	implicit none

	bsysexpl = .true.

	end subroutine system_set_explicit

!******************************************************************

	subroutine system_solve(n,z)

! solves system - z is used for initial guess, not the final solution
#include "pragma_directives.h"

	use mod_system
	use mod_petsc_system
	use mod_system_interface
	use shympi
	use shympi_time

	implicit none

	integer, intent(in) :: n
	real, intent(in)    :: z(n)
	double precision t_start,t_end,t_passed

        !---------------------------------
        ! Petsc Begin/End Assembling :
        !---------------------------------
          call zeta_system%matvec_assemble

          if(petsc_iter>1)then
            continue
          else
            write(*,*)'init_solver'
            t_start = shympi_wtime()
            call zeta_system%init_solver
            t_end = shympi_wtime()
	    t_passed = t_end - t_start
	    call shympi_time_accum(shympi_t_init_solver,t_passed)
            write(*,*)'MPI_SOLVER_INIT_TIME=',t_passed
          endif
          write(6,*)'iter is ',petsc_iter
          petsc_iter=petsc_iter+1

	t_start = shympi_wtime()

        call zeta_system%solve

	t_end = shympi_wtime()
	t_passed = t_end - t_start
	call shympi_time_accum(shympi_t_solve,t_passed)

	end subroutine system_solve

!******************************************************************

        subroutine system_solve_3d(n,nlvdi,nlv,z)
	integer, intent(in) :: n,nlvdi,nlv
	real, intent(in)    :: z(nlvdi,n)
           stop 'CALL TO NON-IMPLEMENTED system_solve_3d'
        end

!******************************************************************

        subroutine system_get(n,z)
        ! copies solution back to z

        implicit none
	integer n
	real z(n)
          call zeta_system%get_solution(n,z)
        end subroutine system_get

!******************************************************************

        subroutine system_get_3d(n,nlvdi,nlv,z)
	integer, intent(in) :: n,nlvdi,nlv
	real, intent(in)    :: z(nlvdi,n)
           stop 'CALL TO NON-IMPLEMENTED system_get_3d'
        end

!******************************************************************

       subroutine system_finalize
        implicit none

          call zeta_system%destroy_solver
          call zeta_system%destroy_matvecobjects
          deallocate(zeta_system)
          call petsc_global_finalize

       end subroutine system_finalize
!==================================================================
      end module mod_petsc
!==================================================================
