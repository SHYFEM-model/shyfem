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
! 21.12.2020	clr	original implementation
! 13.03.2021	ggu	bug fix in add_full_rhs()
! 20.04.2021	clr	alternative implementation to replace pragma directives use_PETSc/SPK/AmgX
!
! notes :
!
! structure of calls:
!
!==================================================================
      module mod_petsc_system
!==================================================================

#include "petsc/finclude/petsc.h"

       use petscvec
       use petscmat
       use petscksp
       use petscdm
       use petscdmlabel
       use mod_petsc_global
       use shympi, only : shympi_barrier,my_id,bmpi
   
       
       implicit none 

        !------------------------------------------
        type :: petsc_system ! components and solver of a system AX=B
           private          
           PetscScalar,public :: mat3x3(3,3)
           PetscScalar,public :: vecx3(3)
           integer :: PETSC_COMM  ! MPI Communicator
           Vec  :: X       ! solution vector
           Vec  :: X_loc   ! local solution vector including ghost nodes
           Vec  :: B       ! rhs
           Mat  :: A       ! Matrix A
           KSP  :: ksp     ! Krylov solver
           PC   :: pc      ! Preconditioner
           PetscScalar, pointer :: p_X_loc(:) ! pointer to local solution
           logical :: solver_is_initialized
           type(c_ptr) :: AmgX_Solver

           character(len=80),public :: AmgX_configfile
           logical :: use_AmgX   
        contains
              procedure :: create_objects           
              procedure :: init_solver
              procedure :: init_solver_PETSc
              procedure :: init_solver_AmgX
              procedure :: add_matvec_values
              procedure :: add_full_rhs
              procedure :: matvec_assemble
              procedure :: solve
              procedure :: solve_PETSc
              procedure :: solve_AmgX
              procedure :: get_solution
              procedure :: reset_zero_entries
              procedure :: destroy_matvecobjects
              procedure :: destroy_solver
              procedure :: destroy_solver_PETSc
              procedure :: destroy_solver_AmgX

        end type petsc_system

        interface petsc_system
           module procedure new_system
        end interface petsc_system
        !------------------------------------------
        abstract interface
             subroutine subroutine_of_sysobj (sysobj)
                Import :: petsc_system
                class(petsc_system),target :: sysobj
             end subroutine subroutine_of_sysobj
        end interface

         
        !------------------------------------------

        PetscInt, parameter ::  three=3
        integer, parameter :: petsc_id=0
        integer, parameter :: amgx_id=1
        character(len=4) :: AmgX_mode='dDDI'

       public :: create_objects,
     +           init_solver,
     +           add_matvec_values,
     +           add_full_rhs,
     +           matvec_assemble,
     +           solve,
     +           get_solution,
     +           reset_zero_entries,
     +           destroy_matvecobjects,
     +           destroy_solver
!==================================================================
       contains
!==================================================================

! *****************************************************
!       class constructor 
! *****************************************************

       type(petsc_system) function new_system(petscconfig,amgxconfig)
        character(len=80) :: petscconfig
        character(len=80) :: amgxconfig
        character(len=10) :: shyfem_solver='none'
        PetscBool opt_found

        new_system%AmgX_configfile=trim(amgxconfig) 
         if( bmpi ) then
            new_system%PETSC_COMM=PETSC_COMM_WORLD
         else
            new_system%PETSC_COMM=PETSC_COMM_SELF
         endif

         write(*,*)'new_system constructor'

         if (trim(petscconfig).ne.'NO_FILE_GIVEN') then
         write(6,*) 'reading petscrc file : ',petscconfig
           call PetscOptionsInsertFILE(
     +                     new_system%PETSC_COMM,
     +                     PETSC_NULL_INTEGER,
     +                     petscconfig, 
     +                     PETSC_TRUE,
     +                     perr )
           if (perr .ne. 0) then
              write(6,*)'Unable to read from config file ',
     +                    petscconfig
              stop
           endif
         else
           write(6,*)'NO petsc_zconfig file was given. '
         endif
#ifdef Verbose
        write(*,*)'Options Get String'
#endif

         call PetscOptionsGetString(
     +                 PETSC_NULL_OPTIONS,
     +                 PETSC_NULL_CHARACTER,
     +                 "-shyfem_solver",
     +                 shyfem_solver,
     +                 opt_found,
     +                 perr)
         call assert(perr.eq.0,'PETScOptionsGetString',perr)
         if(opt_found.neqv. .true.)then
            shyfem_solver='petsc'
         endif
       write(*,*)'read shyfem_solver :',trim(shyfem_solver),'.'
#ifdef Verbose
        write(*,*)'set solver id'
#endif
         if (trim(shyfem_solver)=='amgx') then
           write(*,*)'using shyfem_solver ',shyfem_solver,
     +      ' => AmgX routines '     
            if (trim(amgxconfig).eq.'NO_FILE_GIVEN') then
              write(*,*)'using shyfem_solver ',shyfem_solver,
     +          'requires an AmgX configuration file name in the .str'
              stop "ERROR, AmgX configuration file name is missing"
            endif
           new_system%use_AmgX=.True.
         elseif(trim(shyfem_solver)=='petsc') then                    
           write(*,*)'using shyfem_solver ',shyfem_solver,
     +      ' => PETSc routines '     
           new_system%use_AmgX=.False.
         else
           stop "shyfem_solver must be 'petsc' or 'amgx'"
         endif
#ifdef Verbose
        write(*,*)'new_system done'
#endif

       end function  new_system

! *****************************************************

       subroutine create_objects(self)

#include "petsc/finclude/petsc.h"
        implicit none
        class(petsc_system),target :: self

        PetscInt  :: rowStart_read,rowEnd_read

    
         write(6,*)'PETSc Create Objects'
         !-------------------------------------------------------------        
         ! Initialize PETSc Mass Matrix
         !-------------------------------------------------------------        
#ifdef Verbose
         write(6,*)'PETSc Create Matrix',nodes_loc,nodes_glob
#endif
         call MatCreate(self%PETSC_COMM,
     +              self%A,perr)
         call assert(perr.eq.0,'PETScMatCreate',perr)
         call MatSetSizes(self%A,
     +                      nodes_loc,nodes_loc,
     +                      nodes_glob,nodes_glob,perr)
         call assert(perr.eq.0,'PETScMatSetSize',perr)
         call MatSetType(self%A,MATAIJ,perr) ! matrix type MATAIJ is identical
                                               ! to MATSEQAIJ when constructed with 
                                               ! a single process communicator, and
                                               ! MATMPIAIJ otherwise
         ! to run on the GPU request MATAIJCUSPARSE in the options database
         ! to replace default MATAIJ Type of Matrix A 
         call assert(perr.eq.0,'PETScMatSetType',perr)
         call MatSetFromOptions(self%A,perr)  
         call assert(perr.eq.0,'PETScMatSetFromOptions',perr)
         if( bmpi ) then
             call MatMPIAIJSetPreallocation(self%A,
     +                      PETSC_DECIDE,d_nnz,
     +                      PETSC_DECIDE,o_nnz,
     +                      perr)
         call assert(perr.eq.0,'MatMPIAIJSetPreallocation',perr)
         else
             call MatSEQAIJSetPreallocation(self%A,
     +                      PETSC_DECIDE,d_nnz,
     +                      perr)
         call assert(perr.eq.0,'MatSEQAIJSetPreallocation',perr)
         endif
         call PetscObjectSetName(self%A,'A (Mat)',perr)
         ! --------------------------------------------------------
         call MatSetOption(self%A,
     +                     MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE,
     +                     perr)
         call assert(perr.eq.0,'MatSetFromOpt perr',perr)
         call MatSetUp(self%A,perr)
         call assert(perr.eq.0,'MatSetUp perr',perr)

         ! --------------------------------------------------------
         call MatGetOwnershipRange(self%A,
     +                            rowStart_read,rowEnd_read,perr)
         call assert(rowStart_read==rowStart,
     +                     'rowStart changed',perr)
         call assert(rowEnd_read==rowEnd,'rowStart changed',perr)
         call MatGetSize(self%A,Rsize,Csize,perr)
         call assert(perr.eq.0,'MatGetSize',perr)

#ifdef Verbose
         write(6,'(a,i3,4(a,2i8))')
     +           'PETSc : rank',my_id,' with num loc rows and cols ',
     +       nodes_loc,nodes_loc,' and num glob rows and cols',
     +       nodes_glob,nodes_glob,' owns rows',rowStart,rowEnd,
     +          ' in matrix of size',Rsize,Csize
#endif
         !-------------------------------------------------------------        
         ! Initialize PETSc Vectors
         !-------------------------------------------------------------        
#ifdef Verbose
         write(6,*)'PETSc Create rhs Vector B'
#endif
         call VecCreate(self%PETSC_COMM,self%B,perr)
         call assert(perr.eq.0,'VecCreate B',perr)
         call VecSetSizes(self%B,nodes_loc,nodes_glob,perr)
         call VecSetType(self%B,VECSTANDARD,perr) ! seq on one process and mpi on several
         ! to run on the GPU request VECCUDA in the options database
         ! to replace default VECSTANDARD Type of vector B
         call VecSetFromOptions(self%B,perr) 
         call VecSetOption(self%B,
     +           VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE,perr)
         call PetscObjectSetName(self%B,'B (rhs)',perr)

         !-------------------------------------------------------------        
#ifdef Verbose
         write(6,*)'PETSc Create global Vec X'
#endif
         call VecCreate(self%PETSC_COMM,self%X,perr)
         call assert(perr.eq.0,'VecCreate X',perr)
         call VecSetSizes(self%X,nodes_loc,nodes_glob,perr)
         call VecSetType(self%X,VECSTANDARD,perr)
         ! to run on the GPU request VECCUDA in the options database
         ! to replace default VECSTANDARD Type of vector X  
         call VecSetFromOptions(self%X,perr)  
         call PetscObjectSetName(self%X,'X (distributed)',perr)
         if(bmpi)then
#ifdef Verbose
           write(6,*)'PETSc set X ghosts and create Vec X_loc'
#endif
           call VecMPISetGhost(self%X,nghosts,ghosts,perr)
           call VecGhostGetLocalForm(self%X,self%X_loc,perr)
           call PetscObjectSetName(self%X_loc,'X_loc',perr)
         endif

         !-------------------------------------------------------------        
         call VecGetOwnershipRange(self%B, rowStart,rowEnd,perr)
#ifdef Verbose
         write(6,*)'PETSc : rank',my_id,
     +                       ' owns rows ',rowStart,rowEnd  
#endif
         call VecGetOwnershipRange(self%X, rowStart,rowEnd,perr)

          self%solver_is_initialized=.false.

#ifdef Verbose
         write(6,*)'PETSc done initializing'
         call shympi_barrier
#endif
       end subroutine create_objects

! ************************************************************************
! choose solver to init  
! ************************************************************************
       subroutine init_solver(self)
         class(petsc_system),target :: self
          if (self%use_AmgX) then
             call self%init_solver_AmgX
          else
             call self%init_solver_PETSc
          endif
       end subroutine init_solver

! ************************************************************************
! init the PETSc solver 
! ************************************************************************

       subroutine init_solver_PETSc(self)

         implicit none
         class(petsc_system),target :: self
       
          PetscReal rtol
          PetscBool opt_found
          character(len=10) :: opt_val
         !-------------------------------------------------------------        
         ! setup KSP environment and Linear Solver including conditioner
         !-------------------------------------------------------------        
          if(self%solver_is_initialized)then
             stop 'ERROR ksp solver already initialized'
          else
             self%solver_is_initialized=.true.
          endif
         write(6,*)'PETSc Create KSP Solver'
       
         call KSPCreate(self%PETSC_COMM,self%ksp,perr)
         call assert(perr.eq.0,'KSPCreate perr',perr)

         call KSPSetOperators(self%ksp,
     +                        self%A,self%A,
     +                        perr)
         call assert(perr.eq.0,'KSPSetOperators perr',perr)
       
         call KSPSetUp(self%ksp,perr)
         call assert(perr.eq.0,'KSPSetup perr ',perr)
         ! other ksp solvers can ben entered in the options database at run time
         call KSPSetType(self%ksp,KSPGMRES,perr)
         !call KSPSetType(self%ksp,KSPPREONLY,perr)

         rtol  = 1e-8           
         call KSPSetTolerances(self%ksp,
     +            rtol,
     +            PETSC_DEFAULT_REAL,
     +            PETSC_DEFAULT_REAL,
     +            PETSC_DEFAULT_INTEGER,perr)
         call KSPSetFromOptions(self%ksp,perr)
         call PetscOptionsGetString(
     +                 PETSC_NULL_OPTIONS,
     +                 PETSC_NULL_CHARACTER,
     +                 "-ksp_type",
     +                 opt_val,
     +                 opt_found,
     +                 perr)
         write(*,*)'kps_type=',opt_val
         if(opt_found.neqv. .true. .or. (
     +     trim(opt_val).ne.'preonly'.and.trim(opt_val).ne.'choleski')) 
     +       call KSPSetInitialGuessNonzero(self%ksp,PETSC_TRUE,perr)

         call KSPGetPC(self%ksp,self%pc,perr)
         call assert(perr.eq.0,'KSPGetPC perr ', perr)
         call PCSetType(self%pc,PCBJACOBI,perr)
         !call PCSetType(self%pc,PCLU,perr)
         call assert(perr.eq.0,'PCSetType perr ',perr)
         ! to run on the GPU request MATSOLVERCUSPARSE in the options database
         ! to set PCFactorSetMatSolverType
!        call PCFactorSetMatSolverType(self%pc,
!    +                  MATSOLVERMUMPS,perr);  

          call PCSetFromOptions(self%pc,perr)
#ifdef Verbose
         call shympi_barrier
#endif

       end subroutine init_solver_PETSc

! ************************************************************************
! init the AmgX solver 
! ************************************************************************
       subroutine init_solver_AmgX(self)

       use iso_c_binding

       implicit none
       class(petsc_system),target :: self
       character(len=len_trim(AmgX_mode)+1,kind=c_char) :: modestr
       character(len=len_trim(self%AmgX_configfile)+1,
     +           kind=c_char) :: cfgfile
       external :: CAmgX_GetInitSolver
       external :: CAmgX_GetSolver
       external :: CAmgX_Initialize
       modestr=trim(AmgX_mode) // c_null_char
       cfgfile=trim(self%AmgX_configfile) // c_null_char
        if(self%solver_is_initialized)then
           stop 'ERROR ksp solver already initialized'
        else
           self%solver_is_initialized=.true.
        endif

        write(6,*)'Initialize AmgX solver'
!        call CAmgX_GetInitSolver(self%AmgX_Solver,self%PETSC_COMM,
!     +                        modestr,cfgfile)

         call CAmgX_GetSolver(self%AmgX_Solver)
         call CAmgX_Initialize(self%AmgX_Solver,self%PETSC_COMM,
     +                       modestr,cfgfile,perr)

        write(6,*)'AmgX solver creation is done'
       end subroutine init_solver_AmgX

!****************************************************************       
! insert values of element matrix and rhs vector into system matrix and vector
!****************************************************************

       subroutine add_matvec_values(self,ie)
          implicit none
          class(petsc_system),target :: self
          integer, intent(in):: ie

#ifdef Verbose
          write(6,'(2(a,i3),3(a,f15.7),2f15.7,a,3i3)')
     +          'rank',my_id,' ele ',ie,
     +          ' adds matrix values min: ',minval(self%mat3x3),
     +           ' max:',maxval(self%mat3x3),' vec:',
     +            self%vecx3,' in row,col=',
     +           nodes_eleshy2block(:,ie)
#endif
          call MatSetValues(self%A, ! the matrix that is set
     +                 three,nodes_eleshy2block(:,ie),     ! the number of rows and their global indices 
     +                 three,nodes_eleshy2block(:,ie),     ! the number of columns and their global indices
     +                 self%mat3x3,         ! the block of values to be inserted
     +                 ADD_VALUES,   ! sum with matrix values
     +                 perr)

          call VecSetValues(self%B, ! the matrix that is set
     +                 three,nodes_eleshy2block(:,ie),     ! the number of elements and their global indices 
     +                 self%vecx3,         ! the block of values to be inserted
     +                 ADD_VALUES,   ! sum with matrix values
     +                 perr)

       end subroutine add_matvec_values


! ************************************************************************
! reset to zero the entries of the system matrix A and rhs B vector
! ************************************************************************

       subroutine reset_zero_entries(self)
          implicit none
          class(petsc_system),target :: self
#ifdef Verbose
         call shympi_barrier
         if(my_id==0) write(6,*)'PETSc Preassemble'
#endif
         call MatZeroEntries(self%A,perr)
         call VecZeroEntries(self%B,perr)
#ifdef Verbose
         call shympi_barrier
#endif
       end subroutine reset_zero_entries

! ************************************************************************
! assemble system matrix and rhs B vector after the insertion of the entry values
! ************************************************************************

       subroutine matvec_assemble(self)
          implicit none
          class(petsc_system),target :: self
#ifdef Verbose
         call shympi_barrier
          if(my_id==0) write(*,*)'PETSc assemble'
#endif
          call MatAssemblyBegin(self%A,MAT_FINAL_ASSEMBLY,perr)
          call MatAssemblyEnd(self%A,MAT_FINAL_ASSEMBLY,perr)
          call VecAssemblyBegin(self%B,perr)
          call VecAssemblyEnd(self%B,perr)
#ifdef Verbose
         call MatView(self%A,PETSC_VIEWER_STDOUT_WORLD,perr)
         call VecView(self%B,PETSC_VIEWER_STDOUT_WORLD,perr)
         call shympi_barrier
#endif
       end subroutine matvec_assemble

! ************************************************************************
! add an array of values to the rhs B vector (boundary conditions)
! ************************************************************************
       subroutine add_full_rhs(self,dt,n,array)

           use mod_system
           implicit none
             class(petsc_system),target :: self
	     real dt
             integer n
             real array(n)
             integer k
             PetscInt row
             PetscScalar val
#ifdef Verbose
         call shympi_barrier
#endif
             do k=1,n
                val=dt*array(k)
                row=nodes_shy2block(k)  
                if( row>=rowStart .and. row<rowEnd )then !only assemble inner nodes
#ifdef Verbose
                write(6,'(a,i3,a,f10.5,2(a,i3))')'rank',my_id,
     +                     ' sets vector  value ',array(k),' of node',
     +                         k,' in row ',row
#endif
                call VecSetValue(self%B,
     +                           row,val,
     +                           ADD_VALUES,
     +                           perr) 
                call assert(perr.eq.0,'VecSetVale _add_full_rhs perr',
     +                            perr)
                endif
             end do
#ifdef Verbose
         call shympi_barrier
#endif

       end subroutine add_full_rhs

! ************************************************************************
! solve the linear system of equations
! ************************************************************************
       subroutine solve(self)
         class(petsc_system),target :: self
          if (self%use_AmgX) then
             call self%solve_AmgX
          else
             call self%solve_PETSc
          endif
       end subroutine solve
!*******************************************************************
       subroutine solve_PETSc(self)
        use mod_system
        implicit none

          class(petsc_system),target :: self
          if(self%solver_is_initialized)then
             continue
          else
             stop 'ERROR ksp solver was not initialized'
          endif
#ifdef Verbose
         call shympi_barrier
              if(my_id==0) write(*,*)'PETSc solve system'
#endif
              ! set KSP solver
              call KSPSetOperators(self%ksp,self%A,
     +                  self%A,perr)
              call assert(perr.eq.0,
     +                          'KSPSetOperators perr',perr)

              ! solve
              call KSPSolve(self%ksp,self%B,self%X,perr)
              call assert(perr.eq.0,'KSPSolve perr ',perr)

#ifdef Verbose
              call VecView(self%X,PETSC_VIEWER_STDOUT_WORLD,perr)
              call shympi_barrier
              if(my_id==0) write(*,*)'PETSc system solved'
#endif

       end subroutine solve_PETSc
!*******************************************************************
       subroutine solve_AmgX(self)
        use mod_system

        use iso_c_binding

        implicit none

        class(petsc_system),target :: self
        integer(kind=c_int) :: iters
        integer(kind=c_int) ::  iter
        real(kind=c_double) :: residual
        external :: CAmgX_SetA
        external :: CAmgX_Solve
#ifdef Verbose
        if(my_id==0) write(*,*)'solve system'
#endif
        call CAmgX_SetA(self%AmgX_Solver,self%A,perr) ! AmgX Wrapper
        call CAmgX_Solve(self%AmgX_Solver,self%X,
     +                   self%B,perr)  ! AmgX Wrapper
!       call CAmgX_getiters(self%AmgX_Solver,iters,perr)
!       do iter=0,iters-1
!         call CAmgX_getresidual(self%AmgX_Solver,
!    +                           iter,residual,perr)
!         if(my_id==0)write(6,*)'iter',iter,' residual :',
!    +                            residual
!       enddo


       end subroutine solve_AmgX

! ************************************************************************
! copy petsc solution vector into shyfem solution vector
! ************************************************************************

       subroutine get_solution(self,n,z)
        use mod_system

        implicit none
        class(petsc_system),target :: self

             integer :: n
             real :: z(n)
             integer :: k,row,tmp

             PetscOffset, save :: offset
#ifdef Verbose
         call shympi_barrier
#endif
             offset=0
             if(bmpi)then
               call VecGhostUpdateBegin(self%X,
     +                           INSERT_VALUES,SCATTER_FORWARD,perr)
               call VecGhostUpdateEnd(self%X,
     +                           INSERT_VALUES,SCATTER_FORWARD,perr)   
               call VecGetArrayReadF90(self%X_loc,self%p_X_loc,perr)
             else
               call VecGetArrayReadF90(self%X,self%p_X_loc,perr)
             endif

             tmp=int(offset)

#ifdef Verbose
               if(bmpi)then
                call VecView(self%X_loc,PETSC_VIEWER_STDOUT_WORLD,
     +                           perr)
               else
                call VecView(self%X,PETSC_VIEWER_STDOUT_WORLD,
     +                           perr)
               endif
               write(6,*)'rank ',my_id,' has X n=',n,' offset ',tmp
#endif

           do k=1,n
              row=k
              z(k)=real(self%p_X_loc(row))
           enddo

           if(bmpi)then
           call VecRestoreArrayReadF90(self%X_loc,self%p_X_loc,perr)
           else
           call VecRestoreArrayReadF90(self%X,self%p_X_loc,perr)
           endif
#ifdef Verbose
         call shympi_barrier
#endif

       end subroutine get_solution

! ************************************************************************

       subroutine destroy_matvecobjects(self)
        use shympi
        implicit none
        class(petsc_system),target :: self
          call VecDestroy(self%B,perr) 
          call VecDestroy(self%X,perr)   
          !if(bmpi .or. GhostVec )then
          if(bmpi)  call VecDestroy(self%X_loc,perr)   
          !endif
          call MatDestroy(self%A,perr)

       end subroutine destroy_matvecobjects

! ************************************************************************

       subroutine destroy_solver(self)
        use shympi
        implicit none
        class(petsc_system),target :: self
          if (self%use_AmgX) then
             call self%destroy_solver_AmgX
          else
             call self%destroy_solver_PETSc
          endif
       end subroutine destroy_solver
! ************************************************************************
       subroutine destroy_solver_PETSc(self)
        use shympi
        implicit none
        class(petsc_system),target :: self
        !call PCDestroy(self%pc,perr)
        write(*,*)'Finalize KSP Solver '
        !call PCDestroy(self%pc,perr)
        call KSPDestroy(self%ksp,perr)
       end subroutine destroy_solver_PETSc
! ************************************************************************
       subroutine destroy_solver_AmgX(self)
        use shympi
        implicit none
        class(petsc_system),target :: self
        external :: CAmgX_Finalize
        write(*,*)'Finalize AmgX Solver '
        call CAmgX_Finalize(self%AmgX_Solver,perr) ! AmgX Wrapper
       end subroutine destroy_solver_AmgX

! ************************************************************************
! assert that the logical lcond is verified, otherwise stop the program
! ************************************************************************

       subroutine assert(lcond,msg,icase)
         logical,intent(in) :: lcond
         character(len=*), intent(in) :: msg
         integer, intent(in) :: icase

         if (.not.lcond) then
            write(*,*) 'ERROR MSG: ',msg, icase
            stop 'error in PETSc routines called by mod_petsc_system'
         endif
         return
       end subroutine assert
                 
!==================================================================
      end module mod_petsc_system
!==================================================================

