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
! 20.04.2021	clr	alternative implementation to replace pragma directives use_PETSc/SPK/AmgX
!
! notes :
!
! structure of calls:
!
!==================================================================
	module mod_petsc_global
!==================================================================
#include "petsc/finclude/petsc.h"
        use shympi, only: bmpi,my_id,nkn_global,shympi_barrier
        use petscdm
        use petscdmlabel
        implicit none

        PetscInt :: block3indexes(3)
        PetscInt,public, save :: nodes_loc ! number of nodes of this MPI process, including ghost nodes
        PetscInt,public, save :: nodes_glob ! number of nodes the whole global domain (summ of the inner nodes of all processes)
        PetscErrorCode :: perr     ! error flag
        PetscInt,public    :: rowStart ! first row owned by this process in the global matrix
        PetscInt,public    :: rowEnd   ! rowEnd-1 is the last row owned by this process in the global matrix
        PetscInt,public    :: Rsize    ! Number of rows in the global matrix
        PetscInt,public    :: Csize    ! Number of columns in the global matrix
        PetscInt,protected,public :: nghosts ! number of ghost nodes
        integer,protected,public, allocatable :: nodes_shy2block(:) ! index returning node id in block enumeration when given an internal global node id
        integer,protected,public, allocatable :: nodes_eleshy2block(:,:)
        PetscInt, allocatable :: d_nnz(:)
        ! for every row, number of non zeros in columns corresponding to the ghost nodes:
        PetscInt, allocatable :: o_nnz(:)
        !  list of ghost nodes (in global-block numerotation) :
        PetscInt, allocatable :: ghosts(:) 
        logical,allocatable :: node_is_ghost(:)

        public :: 
     +            petsc_global_initialize,
     +            petsc_global_create_setup,        
     +            petsc_global_close_setup        
      contains
     
      subroutine petsc_global_initialize
#include "petsc/finclude/petsc.h"

         call PetscInitialize(
     +                  PETSC_NULL_CHARACTER,
     +                  perr)
         if (perr .ne. 0) then
            write(6,*)'Unable to initialize PETSc'
            stop
         endif


      end subroutine petsc_global_initialize
!****************************************************
      subroutine petsc_global_create_setup

	use mod_system
	use levels
	use basin
	use shympi
        ! for every row, number of non zeros in columns corresponding to the inner nodes:
        PetscInt  :: rowStart_read,rowEnd_read
        integer :: k,col,ng
        PetscBool opt_found

         write(*,*)'petsc_create_indexes'
         call petsc_create_indexes

         allocate(node_is_ghost(0:nkn_global-1))
         node_is_ghost(:)=.False.
         allocate(d_nnz(rowStart:rowEnd-1))
         d_nnz(:)=0
         allocate(o_nnz(rowStart:rowEnd-1))
         o_nnz(:)=0

         write(*,*)'non_zero_and_ghosts'
         call petsc_identify_non_zeros_and_ghosts
         ! create list of ghost nodes (in global-block numerotation) in array 'ghosts'
         if(bmpi)then
            allocate(ghosts(nghosts))
            ng=0
            do k=1,nkn_local
              col=nodes_shy2block(k)
              if(node_is_ghost(col)) then
                ng=ng+1
                ghosts(ng)=col
              endif
            enddo
         endif
      end subroutine petsc_global_create_setup

! ***********************************************************

      subroutine petsc_global_close_setup
#ifdef Verbose
        write(*,*)'petsc_global_close_setup...'
#endif
         deallocate(d_nnz)
         deallocate(o_nnz)
         deallocate(node_is_ghost)
         if(bmpi) deallocate(ghosts)
#ifdef Verbose
        write(*,*)'petsc_global_close_setup done'
#endif

      end subroutine petsc_global_close_setup

! ***********************************************************

      subroutine petsc_global_finalize
        PetscBool :: Petsc_is_initialized
          call PetscInitialized(Petsc_is_initialized,perr)
          if (Petsc_is_initialized)then
             call PetscFinalize(perr)
          endif
          deallocate(nodes_shy2block)
          deallocate(nodes_eleshy2block)
          write(*,*)"PETSc Finalized" 
      end subroutine petsc_global_finalize

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
#ifdef Verbose
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
#ifdef Verbose
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
#ifdef Verbose
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
#ifdef Verbose
        write(*,*)'PETSc computed nodes rowStart,rowEnd=',
     +             rowStart,rowEnd,' rank',my_id
#endif
        deallocate(nodes_glob2block)
        deallocate(inner_nodes_list)
        deallocate(nodes_by_ranks)
#ifdef Verbose
        call shympi_barrier
#endif
      end subroutine petsc_create_indexes


! ************************************************************************
! identify the non-zeros of the matrix and the ghost nodes
! of every process 
! ************************************************************************

      subroutine petsc_identify_non_zeros_and_ghosts

        use basin, only: nel,nen3v

        implicit none

        integer  k,ie,ie_mpi,numele,row,col,i,j,max_nlocnod
        integer,allocatable :: local_nodes(:,:) 
        integer,allocatable :: numlocnod_per_row(:) 
#ifdef Verbose
        call shympi_barrier
        write(*,'(3(2(a,i10),a,i6))')
     +    ' before working  d_nnz+o_nnz has min=',
     +        minval(d_nnz+o_nnz),' max=',maxval(d_nnz+o_nnz),
     +    ', sum=',sum(d_nnz+o_nnz),
     +    ' ; d_nnz has min=',minval(d_nnz),' max=',maxval(d_nnz),
     +    ' sum=',sum(d_nnz),
     +    ' ; o_nnz has min=',minval(o_nnz),' max=',maxval(o_nnz),
     +    ' sum=',sum(o_nnz)
#endif
        max_nlocnod=12
        allocate(local_nodes(max_nlocnod,rowStart:rowEnd-1))
        local_nodes(:,:)=-1
        allocate(numlocnod_per_row(rowStart:rowEnd-1))
        numlocnod_per_row(:)=0
       
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

        !-------------------------------------------------------------
        ! save the number of ghost nodes: 'nghosts' of every process   
        !-------------------------------------------------------------
         nghosts=0
         if(bmpi)then
            do col=0,nkn_global-1
              if(node_is_ghost(col)) nghosts=nghosts+1
            enddo 
         endif
        deallocate(local_nodes)
        deallocate(numlocnod_per_row)

        write(*,'(a,i3,4(a,i6),3(2(a,i10),a,i6))')'PETSc rank=',my_id,
     +    ' has a number of inner nodes (nrows)=',rowEnd-rowStart-1,
     +    ' (rows ',rowStart,
     +    ' to',rowEnd-1,') and ',nghosts,
     +    ' ghost nodes ; d_nnz+o_nnz has min=',
     +        minval(d_nnz+o_nnz),' max=',maxval(d_nnz+o_nnz),
     +    ', sum=',sum(d_nnz+o_nnz),
     +    ' ; d_nnz has min=',minval(d_nnz),' max=',maxval(d_nnz),
     +    ' sum=',sum(d_nnz),
     +    ' ; o_nnz has min=',minval(o_nnz),' max=',maxval(o_nnz),
     +    ' sum=',sum(o_nnz)
#ifdef Verbose
        write(*,'(a,i4)')
     +  'PETSc done identifying non-zero and ghosts, rank',my_id
        call shympi_barrier
#endif
      end subroutine petsc_identify_non_zeros_and_ghosts

! ***********************************************************

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
!******************************************************************


!==================================================================
      end module mod_petsc_global
!==================================================================
