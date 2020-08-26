
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

#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscksp.h"

        use mpi
        use petscvec
        use petscmat
        use petscksp
        use petscdm
        use petscdmlabel
        
        implicit none 

#define get_zeta(pos)  p_zeta_loc(offset + (pos))

        type :: syspetsc_matrix
          
           Vec  :: zeta   ! solution vector
           Vec  :: zeta_loc   ! local solution vector including ghost nodes
           Vec  :: rhs ! rhs
           Mat  :: MassMat
           KSP  :: ksp              ! Krylov solver
           PC   :: pc               ! Preconditioner
           !MatPartitioning :: part ! pratition object
           IS              :: is
           integer :: PETSC_COMM
          
        end type syspetsc_matrix

        type(syspetsc_matrix), save, target  :: petsc_zeta_solver

        integer, save, allocatable :: nodes_ext2block(:) ! index returning node id in block enumeration when given an external node id
        integer, save, allocatable :: nodes_block2ext(:) ! index returning node id in ext enumeration when given a block node id
        integer, save, allocatable :: nodes_int2block(:) ! index returning node id in block enumeration when given an internal node id
        integer, save, allocatable :: nodes_block2int(:) ! index returning node id in internal enumeration when given a block node id
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
        Logical,parameter :: Debug=.false.
        PetscScalar, pointer :: p_zeta_loc(:)
!==================================================================
        contains
!==================================================================


! ************************************************************************
! init the petsc system of matrix, vectors and solver
! ************************************************************************

        subroutine mod_system_petsc_init(sysobj)

        use basin, only: ipv ! returns external global node number
        use shympi
        implicit none

        type(syspetsc_matrix) :: sysobj
        PetscInt  :: rowStart_read,rowEnd_read
 
        PetscBool :: Petsc_is_initialized
        PetscReal rtol
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


         if( bmpi ) then
            sysobj%PETSC_COMM=PETSC_COMM_WORLD
         else
            sysobj%PETSC_COMM=PETSC_COMM_SELF
         endif

         !-------------------------------------------------------------        
         ! compute number of local rows, create indexes of nodes to 
         ! reorder matrix/vector rows and columns by blocks, 
         ! identify the non-zeros of the matrix as well as the ghost nodes
         !-------------------------------------------------------------        
         call petsc_create_indexes

         call petsc_identify_non_zeros_and_ghosts

         !-------------------------------------------------------------        
         ! Initialize PETSc Mass Matrix
         !-------------------------------------------------------------        
         if(Debug) write(6,*)'PETSc Create Matrix',nodes_loc,nodes_glob
         call MatCreate(sysobj%PETSC_COMM,
     +              sysobj%MassMat,perr)
         call MatSetType(sysobj%MassMat,MATAIJ,perr)
         call MatSetSizes(sysobj%MassMat,
     +                      nodes_loc,nodes_loc,
     +                      nodes_glob,nodes_glob,perr)
         call MatMPIAIJSetPreallocation(sysobj%MassMat,
     +                      PETSC_DECIDE,d_nnz,
     +                      PETSC_DECIDE,o_nnz,
     +                      perr)

         call PetscObjectSetName(sysobj%MassMat,'MassMat',perr)
         ! --------------------------------------------------------
         call MatSetOption(sysobj%MassMat,
     +                     MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE,
     +                     perr)
         call MatSetFromOptions(sysobj%MassMat,perr)  
         call petsc_assert(perr.eq.0,'MatSetFromOpt perr',perr)
         call MatSetUp(sysobj%MassMat,perr)
         call petsc_assert(perr.eq.0,'MatSetUp perr',perr)

         ! --------------------------------------------------------
         call MatGetOwnershipRange(sysobj%MassMat,
     +                            rowStart_read,rowEnd_read,perr)
         call petsc_assert(rowStart_read==rowStart,
     +                     'rowStart changed',perr)
         call petsc_assert(rowEnd_read==rowEnd,'rowStart changed',perr)
         call MatGetSize(sysobj%MassMat,Rsize,Csize,perr)

         if(Debug) write(6,'(a,i3,4(a,2i8))')
     +           'PETSc : rank',my_id,' with num loc rows and cols ',
     +       nodes_loc,nodes_loc,' and num glob rows and cols',
     +       nodes_glob,nodes_glob,' owns rows',rowStart,rowEnd,
     +          ' in matrix of size',Rsize,Csize

         !-------------------------------------------------------------        
         ! Initialize PETSc Vectors
         !-------------------------------------------------------------        
         if(Debug) write(6,*)'PETSc Create Vector rhs'

         call VecCreateMPI(sysobj%PETSC_COMM,nodes_loc,
     +                            nodes_glob,sysobj%rhs,perr)
         call petsc_assert(perr.eq.0,'veccreatesMPI perr',perr)

         call PetscObjectSetName(sysobj%rhs,'rhs',perr)

         !-------------------------------------------------------------        
         if(Debug) write(6,*)'PETSc Create global Vec zeta and zeta_loc'

         call VecCreateGhost(sysobj%PETSC_COMM,nodes_loc,nodes_glob,
     +       nghosts,ghosts,sysobj%zeta,perr)
         call VecGhostGetLocalForm(sysobj%zeta,sysobj%zeta_loc,perr)
         call PetscObjectSetName(sysobj%zeta,'zeta',perr)
         call PetscObjectSetName(sysobj%zeta_loc,'zeta_loc',perr)

         !-------------------------------------------------------------        
         call VecGetOwnershipRange(sysobj%rhs, rowStart,rowEnd,perr)
         if(Debug) write(6,*)'PETSc : rank',my_id,
     +                       ' owns rows ',rowStart,rowEnd  
         call VecGetOwnershipRange(sysobj%zeta, rowStart,rowEnd,perr)

         !-------------------------------------------------------------        
         ! setup KSP environment and Linear Solver including conditioner
         !-------------------------------------------------------------        
         if(Debug) write(6,*)'PETSc Create KSP Solver'
       
         call KSPCreate(sysobj%PETSC_COMM,sysobj%ksp,perr)
         call petsc_assert(perr.eq.0,'KSPCreate perr',perr)

         call KSPSetOperators(sysobj%ksp,
     +                        sysobj%MassMat,sysobj%MassMat,
     +                        perr)
         call petsc_assert(perr.eq.0,'KSPSetOperators perr',perr)
       
         if(Debug) write(6,*)'PETSc Setup KSP solver'
         call KSPSetUp(sysobj%ksp,perr)
         call petsc_assert(perr.eq.0,'KSPSetup perr ',perr)
         call KSPSetType(sysobj%ksp,KSPPREONLY,perr)
         call KSPSetFromOptions(sysobj%ksp,perr)

         rtol  = 1e-5           
         call KSPSetTolerances(sysobj%ksp,
     +            rtol,
     +            PETSC_DEFAULT_REAL,
     +            PETSC_DEFAULT_REAL,
     +            PETSC_DEFAULT_INTEGER,perr)


         if(Debug) write(6,*)'PETSc setup preconditionner'
         call KSPGetPC(sysobj%ksp,sysobj%pc,perr)
         call petsc_assert(perr.eq.0,'KSPGetPC perr ', perr)
         call PCSetType(sysobj%pc,PCLU,perr)
         call petsc_assert(perr.eq.0,'PCSetType perr ',perr)
         call PCFactorSetMatSolverType(sysobj%pc,
!     +                  MATSOLVERSUPERLU_DIST,perr);  
     +                  MATSOLVERMUMPS,perr);  

            deallocate(d_nnz)
            deallocate(o_nnz)
            deallocate(ghosts)

         if(Debug) write(6,*)'PETSc done initializing'
        end subroutine mod_system_petsc_init

!****************************************************************       
! insert values of element matrix and vector into system matrix and vector
!****************************************************************

        subroutine mod_system_petsc_setvalue(ie,kn,mass,rhs,sysobj)
          use shympi, only : my_id
          use basin, only: ipv ! returns external global node number
          implicit none
          type(syspetsc_matrix) :: sysobj
         
          integer ie
          integer kn(3)
          real mass(3,3)
          real rhs(3)
          integer kext(3)
          integer i,j,kr,kc
          PetscInt  row,col
          PetscScalar val
          
          do i =1,3
            kext(i)=ipv(kn(i))
          enddo
!         if(Debug)
!    +    write(6,'(10(a,i3))')'rank',my_id,' ie=',ie,
!    +                         ' kext=(',kext(1),
!    +                              ' ,',kext(2),
!    +                              ' ,',kext(3),
!    +' ) -> matrix rows :(',nodes_ext2block(kext(1))-1,
!    +                  ' ,',nodes_ext2block(kext(2))-1,
!    +                  ' ,',nodes_ext2block(kext(3))-1,
!    +                ' ), row range is ',rowStart,' ,',rowEnd
          do i =1,3
            kr = kext(i)
            row=nodes_ext2block(kr) - 1  ! -1 because petsc indexes start at zero  
            if( row>=rowStart .and. row<rowEnd )then        !only assemble rows of inner nodes
              do j=1,3
                kc = kext(j)
                col=nodes_ext2block(kc) -1
                val=mass(i,j)
                call MatSetValue(sysobj%MassMat,
     +                       row,col,val,
     +                       ADD_VALUES,perr)
               call petsc_assert(perr.eq.0,'MatSetValue perr',perr)
!             if(Debug) write(6,'(a,i3,a,f,2(a,2i3))')'rank',my_id,
!    +                     ' adds  value ',mass(i,j),' of ',kr,kc,
!    +                     ' in row,col=',row,col
               end do
               val=rhs(i)
               call VecSetValue(sysobj%rhs,
     +                          row,val,
     +                          ADD_VALUES,
     +                          perr)
               call petsc_assert(perr.eq.0,'VecSetValue ',perr)
            endif
          end do

        end subroutine mod_system_petsc_setvalue

! ************************************************************************
! reset to zero the entries of the system matrix and rhs vector
! ************************************************************************

        subroutine mod_system_petsc_zeroentries(sysobj)
          use shympi
          implicit none
          type(syspetsc_matrix) :: sysobj
         if(Debug .or. my_id==0) write(6,*)'PETSc Preassemble'
         !/* MatResetPreallocation restores the memory required by users */
!        call MatResetPreallocation(sysobj%MassMat,perr);CHKERRA(perr)
!        call MatSetOption(sysobj%MassMat,
!    +                     MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE,
!    +                     perr);CHKERRA(perr)
         if(Debug)write(6,*)'rank ',my_id,'  MatZeroEntries'
         call MatZeroEntries(sysobj%MassMat,perr)
         if(Debug) write(6,*)'rank ',my_id,'  VecZeroEntries'
         call VecZeroEntries(sysobj%rhs,perr)
         call shympi_barrier

        end subroutine mod_system_petsc_zeroentries

! ************************************************************************
! assemble system matrix and rhs vector after the insertion of the entry values
! ************************************************************************

        subroutine mod_system_petsc_assemble(sysobj)
          use shympi
          implicit none
          type(syspetsc_matrix) :: sysobj
          if(Debug .or. my_id==0) write(*,*)'PETSc assemble'
          call MatAssemblyBegin(sysobj%MassMat,MAT_FINAL_ASSEMBLY,perr)
          call MatAssemblyEnd(sysobj%MassMat,MAT_FINAL_ASSEMBLY,perr)
          if(Debug)
     +      call MatView(sysobj%MassMat,PETSC_VIEWER_STDOUT_WORLD,perr)
          call VecAssemblyBegin(sysobj%rhs,perr)
          call VecAssemblyEnd(sysobj%rhs,perr)
          if(Debug)
     +       call VecView(sysobj%rhs,PETSC_VIEWER_STDOUT_WORLD,perr)

        end subroutine mod_system_petsc_assemble

! ************************************************************************
! add an array of values to the rhs vector (boundary conditions)
! ************************************************************************

        subroutine mod_system_petsc_setvec(dt,n,array,sysobj)

           use mod_system
           use shympi, only : my_id
           use basin, only: ipv
           implicit none

             real dt
             integer n
             real array(n)
             type(syspetsc_matrix) :: sysobj
             integer k
             integer kext
             PetscInt row
             PetscScalar val
             if(Debug) 
     +       write(6,*)'rank',my_id,' set Vector Value for n=',n,array
             do k=1,n
                val=array(k)
                kext=ipv(k)
                row=nodes_ext2block(kext) - 1  ! -1 because petsc indexes start at zero  
                if( row>=rowStart .and. row<rowEnd )then !only assemble inner nodes
                if(Debug) write(6,'(a,i3,a,f10.5,2(a,i3))')'rank',my_id,
     +                     ' sets vector  value ',array(k),' of node',
     +                    kext,' in row ',row
                call VecSetValue(sysobj%rhs,
     +                           row,val,
     +                           ADD_VALUES,
     +                           perr) 
                call petsc_assert(perr.eq.0,'VecSetVale _setvec perr',
     +                            perr)
                endif
             end do

        end subroutine mod_system_petsc_setvec

! ************************************************************************
! solve the linear system of equations
! ************************************************************************

        subroutine mod_system_petsc_solve(n,z,sysobj)
        use mod_system
        use mod_system_interface
        use shympi
        implicit none

             integer n
             real z(n)
             type(syspetsc_matrix) :: sysobj

              if(Debug .or. my_id==0) write(*,*)'PETSc solve system'

              ! set KSP solver
              call KSPSetOperators(sysobj%ksp,sysobj%MassMat,
     +                  sysobj%MassMat,perr)
              call petsc_assert(perr.eq.0,
     +                          'KSPSetOperators perr',perr)

              ! solve
              call KSPSolve(sysobj%ksp,sysobj%rhs,sysobj%zeta,perr)
!             call petsc_assert(perr.eq.0,'KSPSolve perr ',perr)
              if(Debug)
     +          call VecView(sysobj%zeta,PETSC_VIEWER_STDOUT_WORLD,perr)


        end subroutine mod_system_petsc_solve

! ************************************************************************
! copy petsc solution vector into shyfem solution vector
! ************************************************************************

        subroutine mod_system_petsc_zeta2z(n,z,sysobj)
        use mod_system
        use mod_system_interface
        use shympi


        implicit none

             integer :: n
             real :: z(n)
             type(syspetsc_matrix) :: sysobj
             integer :: k,row,tmp

             PetscOffset, save :: offset
             offset=0

              call VecGhostUpdateBegin(sysobj%zeta,
     +                           INSERT_VALUES,SCATTER_FORWARD,perr)
              call VecGhostUpdateEnd(sysobj%zeta,
     +                           INSERT_VALUES,SCATTER_FORWARD,perr)   
              if(Debug) write(6,'(a,i3)',advance='no')'rank',my_id
              if(Debug)
     +          call VecView(sysobj%zeta_loc,PETSC_VIEWER_STDOUT_WORLD,
     +                           perr)

              if(Debug) write(6,*)' ' 

             call VecGetArrayF90(sysobj%zeta_loc,p_zeta_loc,perr)
             tmp=offset
             if(Debug)
     +          write(6,*)'rank ',my_id,' has zeta n=',n,' offset ',tmp

             do k=1,n
                row=nodes_int2block(k)-1
                !z(k)=get_zeta(row+1)
                z(k)=p_zeta_loc(k)
                if(Debug) write(6,'(5(a,i3),a,f13.9)')'rank ',my_id,
     +                ' returns z[',k,' ] = zeta[',row,' ] ; ( kext=',
     +                nodes_block2ext(row+1) ,' p_zeta_loc[',
     +               k,' ] =',z(k)
             enddo

            call VecRestoreArray(sysobj%zeta_loc,p_zeta_loc,offset,perr)

        end subroutine mod_system_petsc_zeta2z

! ************************************************************************
! Destroy the sysobj deallocating memory, destroying sub-objects and
! calling PetscFinalize
! ************************************************************************

        subroutine mod_system_petsc_finalize(sysobj)
        implicit none
        type(syspetsc_matrix) :: sysobj
        PetscBool :: Petsc_is_initialized
          call VecDestroy(sysobj%rhs,perr) 
          call VecDestroy(sysobj%zeta,perr)   
          call MatDestroy(sysobj%MassMat,perr)
          call PetscInitialized(Petsc_is_initialized,perr)
          if (Petsc_is_initialized)then
             call PetscFinalize(perr)
          endif
 
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
        use basin, only: ipv ! returns external global node number
        use shympi

        implicit none

        integer  k,kk,id
        integer nkn_max
        integer, allocatable :: inner_nodes_list(:)
        integer, allocatable :: nodes_by_ranks(:,:)
        nodes_glob=nkn_global 


        nkn_max = shympi_max(nkn_local)
        allocate(inner_nodes_list(nkn_max))
        allocate(nodes_by_ranks(nkn_max,0:n_threads-1))
        allocate(nodes_ext2block(nkn_global))
        allocate(nodes_block2ext(nkn_global))
        allocate(nodes_int2block(nkn_global))
        allocate(nodes_block2int(nkn_global))

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
              nodes_ext2block(nodes_by_ranks(kk,id))=k
              nodes_block2ext(k)=nodes_by_ranks(kk,id)
              if(Debug) write(*,'(4(a,i3))')'PETSc rank',my_id,
     +              ' take node from id',id,' node external index=',
     +              nodes_by_ranks(kk,id),' -> block index=',k
              k=k+1
            else
              exit
            endif
          enddo
        enddo 
        deallocate(inner_nodes_list)
        deallocate(nodes_by_ranks)

        ! compute expected rowStart,rowEnd:
        rowStart=nkn_global+1
        rowEnd=0
        do k=1,nkn_local
           if( id_node(k) /= my_id ) cycle        !only enter inner nodes
           do kk=1,nkn_global
             if(nodes_block2ext(kk)==ip_int_nodes(k,my_id+1))then
                rowStart=min(rowStart,kk-1) ! rowStart numeration starts from 0
                rowEnd=max(rowEnd,kk)       ! rowEnd is excluded from the owned rows
                exit
             endif
           enddo
        enddo
        if(Debug) write(6,*)'rank',my_id,
     +             ' computed nodes rowStart,rowEnd=',rowStart,rowEnd

        nodes_block2int(:)=0
        do k=1,nkn_local
            do kk=1,nkn_global
                if(ipv(k)==nodes_block2ext(kk))then
                   nodes_block2int(kk)=k
                   nodes_int2block(k)=kk
              if(Debug) write(*,'(a,i2,3(a,i4))')'PETSc rank',my_id,
     +              ' with node internal index=',k,' (kext ',
     +              nodes_block2ext(kk),') -> block index',
     +              kk-1
                  exit
                endif
            enddo 
        enddo 

      end subroutine petsc_create_indexes


! ************************************************************************
! identify the non-zeros of the matrix and the ghost nodes
! of every process 
! ************************************************************************

      subroutine petsc_identify_non_zeros_and_ghosts

        use basin, only: ipv, ! returns external global node number
     +                   nel,nen3v
        use shympi

        implicit none

        integer, allocatable :: non_zero_vector(:)
        integer kext(3)
        integer kn(3)
        integer  k,kk,id,ie,ie_mpi,row,col,kr,kc,i,j,checkrow
        logical,allocatable :: node_is_ghost(:)

        allocate(non_zero_vector(0:nkn_global-1))
        allocate(d_nnz(rowStart:rowEnd-1))
        d_nnz(:)=0
        allocate(o_nnz(rowStart:rowEnd-1))
        o_nnz(:)=0
        allocate(node_is_ghost(0:nkn_global-1))
        node_is_ghost(:)=.False.
        DO row=rowStart,rowEnd-1
          !-------------------------------------------------------------
          ! loop over elements to get rows and cols of non-zero matrix entries
          !-------------------------------------------------------------
          non_zero_vector(:)=0
          do ie_mpi=1,nel
            ie = ip_sort_elem(ie_mpi)
            do i=1,3
              kk=nen3v(i,ie)
              kn(i)=kk
              kext(i)=ipv(kn(i))
            end do
            do i=1,3
              kr = kext(i)
              checkrow=nodes_ext2block(kr) - 1  ! -1 because petsc indexes start at zero  
              if( row==checkrow )then        !only consider rows of inner nodes
                do j=1,3
                  kc = kext(j)
      
                  col=nodes_ext2block(kc) -1
                  non_zero_vector(col)=1
                 end do
              endif
            end do
          end do
         
          !-------------------------------------------------------------
          ! store in d_nnz, for every row, the number of non zeros in
          ! columns corresponding to the inner nodes 
          ! store in o_nnz, for every row, the number of non zeros in
          ! columns corresponding to the ghost nodes 
          !-------------------------------------------------------------
          do col=0,nkn_global-1
              if(non_zero_vector(col)>0)then
                 if(col>=rowStart .and. col<rowEnd )then ! diagonal block
                  d_nnz(row)=d_nnz(row)+1
                 else ! non diagonal block (ghost columns)
                  o_nnz(row)=o_nnz(row)+1
                  node_is_ghost(col)=.True.
                 endif
              endif
          enddo
        END DO

        !-------------------------------------------------------------
        ! save the number of ghost nodes: 'nghosts' of every process   
        ! and the list of those ghost nodes (in global-block
        ! numerotation) in array 'ghosts'                                                         
        !-------------------------------------------------------------
        nghosts=0
        do col=0,nkn_global-1
          if(node_is_ghost(col)) nghosts=nghosts+1
        enddo 
        allocate(ghosts(nghosts))
        if(Debug) write(6,'(2(a,i3),a)',advance='no')'rank',my_id,
     +      ' has nghost=',nghosts,', ghost list:'
        nghosts=0
        do k=1,nkn_local
        !do col=0,nkn_global-1
          col=nodes_int2block(k)-1
          if(node_is_ghost(col)) then
            nghosts=nghosts+1
            ghosts(nghosts)=col
            if(Debug) write(6,'(i3)',advance='no')ghosts(nghosts)
          endif
        enddo
        if(Debug)write(6,*)' '

      end subroutine petsc_identify_non_zeros_and_ghosts



!==================================================================
        end module mod_system_petsc
!==================================================================

