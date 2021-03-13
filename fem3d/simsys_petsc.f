
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
! 13.03.2021	ggu&clr	newly written for petsc solver
!
!******************************************************************

! routines we need:
!
! system_initialize
! system_init
! system_set_explicit
! system_solve
! system_solve_3d		empty
! system_assemble
! system_assemble_3d		empty
! system_get
! system_get_3d			empty
! system_add_rhs
! system_adjust_matrix_3d	empty
! system_finalize
!
! not needed:
!
! system_setup_global
! system_solve_global

!==================================================================
	module mod_system_global
!==================================================================

	implicit none

!==================================================================
	end module mod_system_global
!==================================================================

        subroutine system_initialize

! allocates data structure

	use shympi
	use mod_petsc

        implicit none

        write(6,*) '----------------------------------------'
        write(6,*) 'initializing matrix inversion routines'
        write(6,*) 'using PETSC routines ',my_id,nkn
        write(6,*) '----------------------------------------'

	call system_initialize_petsc

        end

!******************************************************************

	subroutine system_init

! initializes arrays (set to zero)
!
! must be called before every assembly of matrix

	use shympi
	use mod_petsc

	implicit none

	call system_init_petsc

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

	use shympi
	use mod_petsc

	implicit none

	integer, intent(in) :: n
	real, intent(in)    :: z(n)

	call system_solve_petsc(n,z)

	end

!******************************************************************

	subroutine system_solve_3d(n,nlvdi,nlv,z)

! solves system - z is used for initial guess, not the final solution

	implicit none

	integer, intent(in) :: n,nlvdi,nlv
	real, intent(in)    :: z(nlvdi,n)

	stop 'error stop system_solve_3d: not yet implemented'

	end

!******************************************************************

	subroutine system_assemble(ie,kn,mass,rhs)

! assembles element matrix into system matrix

	use mod_system
	use shympi
	use mod_petsc

	implicit none

	integer ie
	integer kn(3)
	real mass(3,3)
	real rhs(3)

	call system_assemble_petsc(ie,kn,mass,rhs)

	end

!******************************************************************

	subroutine system_assemble_3d(ie,l,nlv,kn,mass,rhs)

! assembles element matrix into system matrix (3d version)

	implicit none

	integer ie,l,nlv
	integer kn(3)
	real mass(-1:1,3,3)
	real rhs(3)

	stop 'error stop system_assemble_3d: not yet implemented'
	
	end

!******************************************************************

        subroutine system_get(n,z)

! copies solution back to z

	use mod_system
	use mod_petsc

        implicit none

	integer n
	real z(n)

	call system_get_petsc(n,z)

        end

!******************************************************************

        subroutine system_get_3d(n,nlvdi,nlv,z)

! copies solution back to z (was system_adjust_3d)

        implicit none

	integer n,nlvdi,nlv
	real z(nlvdi,n)

	stop 'error stop system_get_3d: not yet implemented'

        end

!******************************************************************

        subroutine system_add_rhs(dt,n,array)

! adds right hand side to system array

	use mod_system
	use mod_petsc

        implicit none

        real dt
	integer n
        real array(n)

	call system_add_rhs_petsc(dt,n,array)

        end

!******************************************************************

        subroutine system_adjust_matrix_3d

        implicit none

	stop 'error stop system_adjust_3d: not yet implemented'

        end

!******************************************************************

        subroutine system_finalize

	use mod_system
	use mod_petsc

        implicit none

	call system_finalize_petsc

	end

!******************************************************************

