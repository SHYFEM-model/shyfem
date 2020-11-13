
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

! module for solving system matrix using petsc
!
! revision log :
!
! notes :
!
! structure of calls:
!
!==================================================================
        module mod_system_petsc
!==================================================================

#include "pragma_directives.h"
#include "petsc/finclude/petsc.h"

        use mpi
        use petscvec
        use petscmat
        use petscksp
        use petscdm
        use petscdmlabel
        
        implicit none 

        type :: syspetsc_matrix ! components and solver of a system AX=B
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

        end type syspetsc_matrix

        type(syspetsc_matrix), save, target  :: petsc_zeta_solver

        integer, save, allocatable :: nodes_shy2block(:) ! index returning node id in block enumeration when given an internal global node id
        PetscInt, save :: nodes_loc ! number of nodes of this MPI process, including ghost nodes
        PetscInt, save :: nodes_glob ! number of nodes the whole global domain (summ of the inner nodes of all processes)
        ! for every row, number of non zeros in columns corresponding to the inner nodes:
        PetscInt, allocatable :: d_nnz(:)
       ! integer, save, allocatable :: d_nnz(:)
        ! for every row, number of non zeros in columns corresponding to the ghost nodes:
        PetscInt, allocatable :: o_nnz(:)
        !  list of ghost nodes (in global-block numerotation) : 
        PetscInt, allocatable :: ghosts(:) 
        PetscInt, save :: nghosts ! number of ghost nodes
        PetscErrorCode :: perr     ! error flag
        PetscInt    :: rowStart ! first row owned by this process in the global matrix
        PetscInt    :: rowEnd   ! rowEnd-1 is the last row owned by this process in the global matrix
        PetscInt    :: Rsize    ! Number of rows in the global matrix
        PetscInt    :: Csize    ! Number of columns in the global matrix
        PetscInt,parameter ::  three=3
        integer, allocatable :: nodes_eleshy2block(:,:)
        PetscInt :: block3indexes(3)
        !Logical,parameter :: GhostVec=.true.
       !interface resize_array
       !  module procedure resize_int_1darray
       !  module procedure resize_int_2darray
       !end interface resize_array
!==================================================================
        contains
!==================================================================


! ************************************************************************
! init the petsc system of matrix, vectors and solver
! ************************************************************************

        subroutine mod_system_petsc_init(sysobj)

        use shympi

        implicit none

        type(syspetsc_matrix) :: sysobj
        PetscInt  :: rowStart_read,rowEnd_read
 
#ifdef _Debug
         call shympi_barrier
#endif


        !-------------------------------------------------------------        
        ! Initialize Petsc 
        !-------------------------------------------------------------        

         write(6,*) 'PETSc Initialization'
         call PetscInitialize(
     +       PETSC_NULL_CHARACTER,
     +      perr)
         if (perr .ne. 0) then
            write(6,*)'Unable to initialize PETSc'
            stop
         endif

         !-------------------------------------------------------------        
         ! compute number of local rows, create indexes of nodes to 
         ! reorder matrix/vector rows and columns by blocks, 
         ! identify the non-zeros of the matrix as well as the ghost nodes
         !-------------------------------------------------------------        
         call petsc_create_indexes
         write(6,*)'PETSc indexes created'


         call petsc_identify_non_zeros_and_ghosts
         write(6,*)'PETSc non zeros and ghosts identified'


         if( bmpi ) then
            sysobj%PETSC_COMM=PETSC_COMM_WORLD
         else
            sysobj%PETSC_COMM=PETSC_COMM_SELF
         endif
    


         !-------------------------------------------------------------        
         ! Initialize PETSc Mass Matrix
         !-------------------------------------------------------------        
#ifdef _Debug
         write(6,*)'PETSc Create Matrix',nodes_loc,nodes_glob
#endif
         call MatCreate(sysobj%PETSC_COMM,
     +              sysobj%A,perr)
         call MatSetType(sysobj%A,MATAIJ,perr) ! matrix type MATAIJ is identical
                                               ! to MATSEQAIJ when constructed with 
                                               ! a single process communicator, and
                                               ! MATMPIAIJ otherwise
         call MatSetSizes(sysobj%A,
     +                      nodes_loc,nodes_loc,
     +                      nodes_glob,nodes_glob,perr)
         if( bmpi ) then
             call MatMPIAIJSetPreallocation(sysobj%A,
     +                      PETSC_DECIDE,d_nnz,
     +                      PETSC_DECIDE,o_nnz,
     +                      perr)
         else
             call MatSEQAIJSetPreallocation(sysobj%A,
     +                      PETSC_DECIDE,d_nnz,
     +                      perr)
         endif
         call PetscObjectSetName(sysobj%A,'A (Mat)',perr)
         ! --------------------------------------------------------
         call MatSetOption(sysobj%A,
     +                     MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE,
     +                     perr)
         call MatSetFromOptions(sysobj%A,perr)  
         call petsc_assert(perr.eq.0,'MatSetFromOpt perr',perr)
         call MatSetUp(sysobj%A,perr)
         call petsc_assert(perr.eq.0,'MatSetUp perr',perr)

         ! --------------------------------------------------------
         call MatGetOwnershipRange(sysobj%A,
     +                            rowStart_read,rowEnd_read,perr)
         call petsc_assert(rowStart_read==rowStart,
     +                     'rowStart changed',perr)
         call petsc_assert(rowEnd_read==rowEnd,'rowStart changed',perr)
         call MatGetSize(sysobj%A,Rsize,Csize,perr)

#ifdef _Debug
         write(6,'(a,i3,4(a,2i8))')
     +           'PETSc : rank',my_id,' with num loc rows and cols ',
     +       nodes_loc,nodes_loc,' and num glob rows and cols',
     +       nodes_glob,nodes_glob,' owns rows',rowStart,rowEnd,
     +          ' in matrix of size',Rsize,Csize
#endif
         !-------------------------------------------------------------        
         ! Initialize PETSc Vectors
         !-------------------------------------------------------------        
#ifdef _Debug
         write(6,*)'PETSc Create rhs Vector B'
#endif
         call VecCreate(sysobj%PETSC_COMM,sysobj%B,perr)
         call VecSetType(sysobj%B,VECSTANDARD,perr) ! seq on one process and mpi on several
         call VecSetSizes(sysobj%B,nodes_loc,nodes_glob,perr)
         call VecSetOption(sysobj%B,
     +           VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE,perr)
         call PetscObjectSetName(sysobj%B,'B (rhs)',perr)

         !-------------------------------------------------------------        
#ifdef _Debug
         write(6,*)'PETSc Create global Vec X'
#endif
         call VecCreate(sysobj%PETSC_COMM,sysobj%X,perr)
         call VecSetType(sysobj%X,VECSTANDARD,perr)
         call VecSetSizes(sysobj%X,nodes_loc,nodes_glob,perr)
         call PetscObjectSetName(sysobj%X,'X (distributed)',perr)
         if(bmpi)then
#ifdef _Debug
           write(6,*)'PETSc set X ghosts and create Vec X_loc'
#endif
           call VecMPISetGhost(sysobj%X,nghosts,ghosts,perr)
           call VecGhostGetLocalForm(sysobj%X,sysobj%X_loc,perr)
           call PetscObjectSetName(sysobj%X_loc,'X_loc',perr)
         endif

         !-------------------------------------------------------------        
         call VecGetOwnershipRange(sysobj%B, rowStart,rowEnd,perr)
#ifdef _Debug
         write(6,*)'PETSc : rank',my_id,
     +                       ' owns rows ',rowStart,rowEnd  
#endif
         call VecGetOwnershipRange(sysobj%X, rowStart,rowEnd,perr)

          write(*,*)' Mat and Vec Created'

            deallocate(d_nnz)
            deallocate(o_nnz)
          if(bmpi) deallocate(ghosts)


#ifdef _Debug
         write(6,*)'PETSc done initializing'
         call shympi_barrier
#endif
        end subroutine mod_system_petsc_init

        
! ************************************************************************
! init the PETSc solver 
! ************************************************************************

        subroutine mod_system_petsc_init_PETSc_solver(sysobj)
         use shympi, only : shympi_barrier
         implicit none
          type(syspetsc_matrix) :: sysobj
             PetscReal rtol

         !-------------------------------------------------------------        
         ! setup KSP environment and Linear Solver including conditioner
         !-------------------------------------------------------------        
#ifdef _Debug
         call shympi_barrier
         write(6,*)'PETSc Create KSP Solver'
#endif
       
         call KSPCreate(sysobj%PETSC_COMM,sysobj%ksp,perr)
         call petsc_assert(perr.eq.0,'KSPCreate perr',perr)

         call KSPSetOperators(sysobj%ksp,
     +                        sysobj%A,sysobj%A,
     +                        perr)
         call petsc_assert(perr.eq.0,'KSPSetOperators perr',perr)
       
         call KSPSetUp(sysobj%ksp,perr)
         call petsc_assert(perr.eq.0,'KSPSetup perr ',perr)
         call KSPSetType(sysobj%ksp,KSPGMRES,perr)
         !call KSPSetType(sysobj%ksp,KSPPREONLY,perr)

         rtol  = 1e-8           
         call KSPSetTolerances(sysobj%ksp,
     +            rtol,
     +            PETSC_DEFAULT_REAL,
     +            PETSC_DEFAULT_REAL,
     +            PETSC_DEFAULT_INTEGER,perr)
         call KSPSetInitialGuessNonzero(sysobj%ksp,PETSC_TRUE,perr)
         call KSPSetFromOptions(sysobj%ksp,perr)

         call KSPGetPC(sysobj%ksp,sysobj%pc,perr)
         call petsc_assert(perr.eq.0,'KSPGetPC perr ', perr)
         call PCSetType(sysobj%pc,PCBJACOBI,perr)
         !call PCSetType(sysobj%pc,PCLU,perr)
         call petsc_assert(perr.eq.0,'PCSetType perr ',perr)
!        call PCFactorSetMatSolverType(sysobj%pc,
!    +                  MATSOLVERMUMPS,perr);  
          call PCSetFromOptions(sysobj%pc,perr)
#ifdef _Debug
         call shympi_barrier
#endif

        end subroutine mod_system_petsc_init_PETSc_solver

!****************************************************************       
! insert values of element matrix and rhs vector into system matrix and vector
!****************************************************************

        subroutine mod_system_petsc_setvalues(ie,sysobj)
          implicit none
          integer, intent(in):: ie
          type(syspetsc_matrix), intent(inout) :: sysobj

#ifdef _Debug
          write(6,'(2(a,i3),3(a,f15.7),2f15.7,a,3i3)')
     +          'rank',my_id,' ele ',ie,
     +          ' adds matrix values min: ',minval(sysobj%mat3x3),
     +           ' max:',maxval(sysobj%mat3x3),' vec:',
     +            sysobj%vecx3,' in row,col=',
     +           nodes_eleshy2block(:,ie)
#endif
          call MatSetValues(sysobj%A, ! the matrix that is set
     +                 three,nodes_eleshy2block(:,ie),     ! the number of rows and their global indices 
     +                 three,nodes_eleshy2block(:,ie),     ! the number of columns and their global indices
     +                 sysobj%mat3x3,         ! the block of values to be inserted
     +                 ADD_VALUES,   ! sum with matrix values
     +                 perr)

          call VecSetValues(sysobj%B, ! the matrix that is set
     +                 three,nodes_eleshy2block(:,ie),     ! the number of elements and their global indices 
     +                 sysobj%vecx3,         ! the block of values to be inserted
     +                 ADD_VALUES,   ! sum with matrix values
     +                 perr)

        end subroutine mod_system_petsc_setvalues


! ************************************************************************
! reset to zero the entries of the system matrix A and rhs B vector
! ************************************************************************

        subroutine mod_system_petsc_zeroentries(sysobj)
          use shympi
          implicit none
          type(syspetsc_matrix) :: sysobj
#ifdef _Debug
         call shympi_barrier
         if(my_id==0) write(6,*)'PETSc Preassemble'
#endif
         !/* MatResetPreallocation restores the memory required by users */
!        call MatResetPreallocation(sysobj%A,perr);CHKERRA(perr)
!        call MatSetOption(sysobj%A,
!    +                     MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE,
!    +                     perr);CHKERRA(perr)
         call MatZeroEntries(sysobj%A,perr)
         call VecZeroEntries(sysobj%B,perr)
#ifdef _Debug
         call shympi_barrier
#endif
        end subroutine mod_system_petsc_zeroentries

! ************************************************************************
! assemble system matrix and rhs B vector after the insertion of the entry values
! ************************************************************************

        subroutine mod_system_petsc_assemble(sysobj)
          use shympi
          implicit none
          type(syspetsc_matrix) :: sysobj
#ifdef _Debug
         call shympi_barrier
          if(my_id==0) write(*,*)'PETSc assemble'
#endif
          call MatAssemblyBegin(sysobj%A,MAT_FINAL_ASSEMBLY,perr)
          call MatAssemblyEnd(sysobj%A,MAT_FINAL_ASSEMBLY,perr)
          call VecAssemblyBegin(sysobj%B,perr)
          call VecAssemblyEnd(sysobj%B,perr)
#ifdef _Debug
         call MatView(sysobj%A,PETSC_VIEWER_STDOUT_WORLD,perr)
         call VecView(sysobj%B,PETSC_VIEWER_STDOUT_WORLD,perr)
         call shympi_barrier
#endif
        end subroutine mod_system_petsc_assemble

! ************************************************************************
! add an array of values to the rhs B vector (boundary conditions)
! ************************************************************************

!            allocate(setvalue_indexes(0:n-1)
!            setvalue_indexes(:)=-1 ! negative indexes will be ignored
!            do k=1,n
!               row=nodes_shy2block(k)  
!               if( row>=rowStart .and. row<rowEnd )then !only assemble inner nodes
!               setvalue_indexes(k-1)=row
!            enddo
        subroutine mod_system_petsc_setvec(n,array,sysobj)

           use mod_system
           use shympi, only : shympi_barrier
           implicit none

             integer n
             real array(n)
             type(syspetsc_matrix) :: sysobj
             integer k
             PetscInt row
             PetscScalar val
#ifdef _Debug
         call shympi_barrier
#endif
             do k=1,n
                val=array(k)
                row=nodes_shy2block(k)  
                if( row>=rowStart .and. row<rowEnd )then !only assemble inner nodes
#ifdef _Debug
                write(6,'(a,i3,a,f10.5,2(a,i3))')'rank',my_id,
     +                     ' sets vector  value ',array(k),' of node',
     +                         k,' in row ',row
#endif
                call VecSetValue(sysobj%B,
     +                           row,val,
     +                           ADD_VALUES,
     +                           perr) 
                call petsc_assert(perr.eq.0,'VecSetVale _setvec perr',
     +                            perr)
                endif
             end do
#ifdef _Debug
         call shympi_barrier
#endif

        end subroutine mod_system_petsc_setvec

! ************************************************************************
! solve the linear system of equations
! ************************************************************************

        subroutine mod_system_petsc_solve(sysobj)
        use mod_system
        use mod_system_interface
        use shympi
        implicit none

             type(syspetsc_matrix) :: sysobj
#ifdef _Debug
         call shympi_barrier
              if(my_id==0) write(*,*)'PETSc solve system'
#endif
              ! set KSP solver
              call KSPSetOperators(sysobj%ksp,sysobj%A,
     +                  sysobj%A,perr)
              call petsc_assert(perr.eq.0,
     +                          'KSPSetOperators perr',perr)

              ! solve
              call KSPSolve(sysobj%ksp,sysobj%B,sysobj%X,perr)
              call petsc_assert(perr.eq.0,'KSPSolve perr ',perr)

#ifdef _Debug
              call VecView(sysobj%X,PETSC_VIEWER_STDOUT_WORLD,perr)
              call shympi_barrier
              if(my_id==0) write(*,*)'PETSc system solved'
#endif

        end subroutine mod_system_petsc_solve

! ************************************************************************
! copy petsc solution vector into shyfem solution vector
! ************************************************************************

        subroutine mod_system_petsc_get_solution(n,z,sysobj)
        use mod_system
        use mod_system_interface
        use shympi

        implicit none

             integer :: n
             real :: z(n)
             type(syspetsc_matrix) :: sysobj
             integer :: k,row,tmp

             PetscOffset, save :: offset
#ifdef _Debug
         call shympi_barrier
#endif
             offset=0
             if(bmpi)then
               call VecGhostUpdateBegin(sysobj%X,
     +                           INSERT_VALUES,SCATTER_FORWARD,perr)
               call VecGhostUpdateEnd(sysobj%X,
     +                           INSERT_VALUES,SCATTER_FORWARD,perr)   
               call VecGetArrayReadF90(sysobj%X_loc,sysobj%p_X_loc,perr)
             else
               call VecGetArrayReadF90(sysobj%X,sysobj%p_X_loc,perr)
             endif

             tmp=int(offset)

#ifdef _Debug
               if(bmpi)then
                call VecView(sysobj%X_loc,PETSC_VIEWER_STDOUT_WORLD,
     +                           perr)
               else
                call VecView(sysobj%X,PETSC_VIEWER_STDOUT_WORLD,
     +                           perr)
               endif
               write(6,*)'rank ',my_id,' has X n=',n,' offset ',tmp
#endif

           do k=1,n
              row=k
              z(k)=real(sysobj%p_X_loc(row))
           enddo

           if(bmpi)then
           call VecRestoreArrayReadF90(sysobj%X_loc,sysobj%p_X_loc,perr)
           else
           call VecRestoreArrayReadF90(sysobj%X,sysobj%p_X_loc,perr)
           endif
#ifdef _Debug
         call shympi_barrier
#endif

        end subroutine mod_system_petsc_get_solution

! ************************************************************************
! Destroy the sysobj deallocating memory, destroying sub-objects and
! calling PetscFinalize
! ************************************************************************

        subroutine mod_system_petsc_finalize(sysobj)
        use shympi
        implicit none
        type(syspetsc_matrix) :: sysobj
        PetscBool :: Petsc_is_initialized
          call VecDestroy(sysobj%B,perr) 
          call VecDestroy(sysobj%X,perr)   
          !if(bmpi .or. GhostVec )then
          if(bmpi)  call VecDestroy(sysobj%X_loc,perr)   
          !endif
          call MatDestroy(sysobj%A,perr)
          call PetscInitialized(Petsc_is_initialized,perr)
          if (Petsc_is_initialized)then
             call PetscFinalize(perr)
          endif
          deallocate(nodes_shy2block)
          deallocate(nodes_eleshy2block)
 
        end subroutine mod_system_petsc_finalize

! ************************************************************************
! assert that the logical lcond is verified, otherwise stop the program
! ************************************************************************

      subroutine petsc_assert(lcond,msg,icase)
         logical,intent(in) :: lcond
         character(len=*), intent(in) :: msg
         integer, intent(in) :: icase

         if (.not.lcond) then
            write(*,*) msg, icase
            stop 'error in PETSc routines called by mod_system_petsc'
         endif
         return
      end subroutine petsc_assert

! ************************************************************************
! compute number or local rows and
! create indexes of nodes to reorder matrix/vector rows so that
! every rank owns a single block of rows containing only its own
! nodes.
! ************************************************************************

      subroutine petsc_create_indexes()
        use basin, only: ipv ! returns internal global node number
        use shympi

        implicit none

        integer  k,kk,id
        integer nkn_max
        integer, allocatable :: inner_nodes_list(:)
        integer, allocatable :: nodes_by_ranks(:,:)
        integer, allocatable :: nodes_glob2block(:) ! index returning node id in block enumeration when given an internal global node id
#ifdef _Debug
        call shympi_barrier
#endif
        nodes_glob=nkn_global 

        allocate(nodes_shy2block(nkn_local))
        nodes_shy2block(:)=-999      

        nkn_max = shympi_max(nkn_local)
        allocate(inner_nodes_list(nkn_max))
        allocate(nodes_by_ranks(nkn_max,0:n_threads-1))
        allocate(nodes_glob2block(nkn_global))

        inner_nodes_list(:)=-1 
        nodes_loc=0
        do k=1,nkn_local
            if( id_node(k) /= my_id ) cycle        !only enter inner nodes
            nodes_loc=nodes_loc+1
            inner_nodes_list(k)= ip_int_nodes(k,my_id+1)
        enddo
        call shympi_gather(inner_nodes_list,nodes_by_ranks)
        k=1
        do id=0,n_threads-1
          do kk=1,nkn_max
            if(nodes_by_ranks(kk,id)>0) then
              nodes_glob2block(nodes_by_ranks(kk,id))=k
#ifdef _Debug
              write(*,'(4(a,i3))')'PETSc rank',my_id,
     +              ' take node from id',id,' node glob_int index=',
     +              nodes_by_ranks(kk,id),' -> block index -1 =',k-1
#endif
              k=k+1
            else
              exit
            endif
          enddo
        enddo 
        if(bmpi)then
           do k=1,nkn_local
#ifdef _Debug
              write(*,'(4(a,i3))')'rank,',my_id,' shynode ',k,
     +                  ' (ext ',ipv(k),
     +                  ') -> block id:',nodes_glob2block(ipv(k))-1
#endif
              nodes_shy2block(k)=nodes_glob2block(ipv(k))-1
           enddo
        else
           do k=1,nkn_local
              nodes_shy2block(k)=k-1
           enddo
        endif

        ! compute expected rowStart,rowEnd:
        rowStart=nkn_global+1
        rowEnd=0
        do kk=1,nodes_loc
         if(inner_nodes_list(kk)>0)then
         rowStart=min(rowStart,nodes_glob2block(inner_nodes_list(kk))-1)! rowStart numeration starts from 0
         rowEnd=max(rowEnd,nodes_glob2block(inner_nodes_list(kk)))      ! rowEnd is excluded from the owned rows
         else
             exit
         endif
        enddo
#ifdef _Debug
        write(*,*)'PETSc computed nodes rowStart,rowEnd=',
     +             rowStart,rowEnd,' rank',my_id
#endif
        deallocate(nodes_glob2block)
        deallocate(inner_nodes_list)
        deallocate(nodes_by_ranks)
#ifdef _Debug
        call shympi_barrier
#endif
      end subroutine petsc_create_indexes


! ************************************************************************
! identify the non-zeros of the matrix and the ghost nodes
! of every process 
! ************************************************************************

      subroutine petsc_identify_non_zeros_and_ghosts

        use basin, only: nel,nen3v
        use shympi

        implicit none

        integer  k,ie,ie_mpi,numele,row,col,i,j,max_nlocnod
        logical,allocatable :: node_is_ghost(:)
        integer,allocatable :: local_nodes(:,:) 
        integer,allocatable :: numlocnod_per_row(:) 
#ifdef _Debug
        call shympi_barrier
#endif
        max_nlocnod=12
        allocate(local_nodes(max_nlocnod,rowStart:rowEnd-1))
        local_nodes(:,:)=-1
        allocate(numlocnod_per_row(rowStart:rowEnd-1))
        numlocnod_per_row(:)=0
        allocate(node_is_ghost(0:nkn_global-1))
        node_is_ghost(:)=.False.
        allocate(d_nnz(rowStart:rowEnd-1))
        d_nnz(:)=0
        allocate(o_nnz(rowStart:rowEnd-1))
        o_nnz(:)=0
       
        allocate(nodes_eleshy2block(3,nel))
!        do ie_mpi=1,nel
!          ie = ip_sort_elem(ie_mpi)
        do ie=1,nel
          do i=1,3
            nodes_eleshy2block(i,ie)=nodes_shy2block( nen3v(i,ie) )
          end do
          do i=1,3
            row=nodes_eleshy2block(i,ie)
            if( row>=rowStart.and.row<rowEnd )then        !only consider rows of inner nodes
              do j=1,3
                  col=nodes_eleshy2block(j,ie)
                  if(numlocnod_per_row(row)+1>=max_nlocnod)then
                      max_nlocnod=max_nlocnod+2
                      call resize_2darray(local_nodes,
     +                                  max_nlocnod,rowEnd-rowStart,   ! new sizes of dim 1 and 2
     +                                  -1) ! default value
                  endif
                  if( any(local_nodes(:,row)==col) )then
                        continue
                  else
                         numlocnod_per_row(row)=numlocnod_per_row(row)+1
                         local_nodes(numlocnod_per_row(row),row)=col
                         if(col>=rowStart .and. col<rowEnd )then ! diagonal block
                          d_nnz(row)=d_nnz(row)+1
                         else ! non diagonal block (ghost columns)
                          o_nnz(row)=o_nnz(row)+1
                          node_is_ghost(col)=.True.
                         endif
                  endif
              end do
            endif
          end do
        end do
        write(*,*)'PETSc done identifying ele per row, rank',my_id,nel,
     +    ' nrows,nlocele=',rowEnd-rowStart-1

#ifdef _Debug
        write(*,'(4(a,i2),2(a,4i2))')'rank',my_id,
     +             'maxval(d_nnz+o_nnz)=',maxval(d_nnz+o_nnz),
     +            ' maxval(d_nnz)=',maxval(d_nnz),
     +            ' maxval(o_nnz)=',maxval(o_nnz),
     +            ' d_nnz=',d_nnz(rowStart+15:rowStart+50:10),
     +            ' o_nnz=',o_nnz(rowStart+15:rowStart+50:10)
#endif

        !-------------------------------------------------------------
        ! save the number of ghost nodes: 'nghosts' of every process   
        ! and the list of those ghost nodes (in global-block
        ! numerotation) in array 'ghosts'                                                         
        !-------------------------------------------------------------
        nghosts=0
        if(bmpi)then
           do col=0,nkn_global-1
             if(node_is_ghost(col)) nghosts=nghosts+1
           enddo 
           allocate(ghosts(nghosts))
           nghosts=0
           do k=1,nkn_local
             col=nodes_shy2block(k)
             if(node_is_ghost(col)) then
               nghosts=nghosts+1
               ghosts(nghosts)=col
             endif
           enddo
           write(6,'(5(a,i8))')'PETSc : rank',my_id,
     +         ' has ',rowEnd-rowStart-1,
     +         ' inner nodes and ',nghosts,' ghost nodes, sum(d_nnz)=',
     +         sum(d_nnz),' while sum(o_nnz)=',sum(o_nnz)
        endif
        deallocate(node_is_ghost)
        deallocate(local_nodes)
        deallocate(numlocnod_per_row)
        write(*,*)'PETSc done identifying non-zero and ghosts, rank',
     +            my_id
#ifdef _Debug
        call shympi_barrier
#endif
      end subroutine petsc_identify_non_zeros_and_ghosts

        subroutine resize_1darray(array,newsize,default_val)
         implicit none
          integer, intent(in):: newsize,default_val
          integer, dimension(:), allocatable, intent(inout):: array
          integer, dimension(:), allocatable :: tmparray
          integer :: lowerbound,upperbound,minsize
          lowerbound=lbound(array,1)
          upperbound=ubound(array,1)
          allocate(tmparray( lowerbound : upperbound ))
          tmparray(:)=array(:)
          deallocate(array)
          allocate(array(lowerbound:lowerbound+newsize))
          array(:)=default_val
          minsize=min(upperbound-lowerbound,newsize)
          array(lowerbound:lowerbound+minsize)=
     +         tmparray(lowerbound:lowerbound+minsize)
          deallocate(tmparray)
        end subroutine resize_1darray
          

        subroutine resize_2darray(array,
     +                                news1,news2,default_val)
         implicit none
          integer, intent(in):: news1,news2
          integer, intent(in):: default_val
          integer, dimension(:,:), allocatable, intent(inout):: array
          integer, dimension(:,:), allocatable :: tmparray
          integer :: mins1,mins2,lbound1,lbound2,olds1,olds2
          lbound1=lbound(array,1)
          lbound2=lbound(array,2)
          olds1=ubound(array,1)-lbound1
          olds2=ubound(array,2)-lbound2
          allocate(tmparray( lbound1 : lbound1+olds1 ,
     +                       lbound2 : lbound2+olds2 ))
          tmparray(:,:)=array(:,:)
          deallocate(array)
          allocate(array( lbound1 : lbound1+news1 ,
     +                       lbound2 : lbound2+news2 ))
          array(:,:)=default_val
          mins1=min(olds1,news1)
          mins2=min(olds2,news2)
          array(lbound1:lbound1+mins1,lbound2:lbound2+mins2)=
     +         tmparray(lbound1:lbound1+mins1,lbound2:lbound2+mins2)
          deallocate(tmparray)
        end subroutine resize_2darray
                  
!==================================================================
        end module mod_system_petsc
!==================================================================

