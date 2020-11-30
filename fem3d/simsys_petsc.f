
!--------------------------------------------------------------------------
!
!    Copyright (C) 2009-2011,2014-2015,2017-2019  Georg Umgiesser
!    Copyright (C) 2015  Debora Bellafiore
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

! revision log :
!
! 04.06.2020	ggu	creating routine calling PETSc libray to solve the system of equations
!
!******************************************************************

!==================================================================
	module mod_system_global
!==================================================================

        implicit none 

!==================================================================
	end module mod_system_global
!==================================================================

        subroutine system_initialize

! allocates data structure

	use mod_system
	use mod_system_petsc
	use levels
	use basin
	use shympi

        implicit none

	call shympi_barrier

          write(6,*) '----------------------------------------'
          write(6,*) 'initializing petsc library ',my_id,nkn
          write(6,*) '----------------------------------------'
          call mod_system_petsc_init(petsc_zeta_solver)

        end

!******************************************************************

	subroutine system_init
	use mod_system_petsc

        call mod_system_petsc_zeroentries(petsc_zeta_solver)

	end

!******************************************************************

	subroutine system_set_explicit

	use mod_system

	implicit none

	bsysexpl = .true.

	end

!******************************************************************

	subroutine system_solve(n,z)

! solves system - z is used for initial guess, not the final solution
#include "pragma_directives.h"

	use mod_system
	use mod_system_petsc
	use mod_system_interface
	use shympi
	use shympi_time

	implicit none

	integer, intent(in) :: n
	real, intent(in)    :: z(n)
	double precision t_start,t_end,t_passed

	t_start = shympi_wtime()

        call mod_system_petsc_solve(petsc_zeta_solver)

	t_end = shympi_wtime()
	t_passed = t_end - t_start
	call shympi_time_accum(shympi_t_solve,t_passed)

	end

!******************************************************************

	subroutine system_solve_3d(n,nlvdi,nlv,z)


	end

!******************************************************************

        subroutine system_get(n,z)

! copies solution back to z

	use mod_system
	use mod_system_petsc, only: mod_system_petsc_get_solution,
     +                              petsc_zeta_solver

        implicit none

	integer n
	real z(n)

          call mod_system_petsc_get_solution(n,z, 
     +                                       petsc_zeta_solver)

        end

!******************************************************************

        subroutine system_get_3d(n,nlvdi,nlv,z)

        end

!******************************************************************

        subroutine system_finalize

	use shympi
	use mod_system_petsc

        implicit none

        call mod_system_petsc_finalize(petsc_zeta_solver)

        end
