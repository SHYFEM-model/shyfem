  subroutine partPHG()  ! MAIN subroutine for partitioning

    use mpi_global_graph_elem
    use mod_geom
    use basin
    use levels
    !use shympi

    implicit none

    real(Zoltan_FLOAT) :: version

    integer i, j, k, id, numGlobalNeighbors, numGlobalVertices, nnbors, procID, num, nbr

    integer partNel, minNel, maxNel
    
    type(Zoltan_Struct), pointer :: zz_obj, zz_obj2
    integer :: st(MPI_STATUS_SIZE), ierr

    type(GRAPH_DATA) :: mygraph
    type(GRAPH_DATA), dimension(:), allocatable :: send_graph
    type (IELTV_GRID) :: ieltGrid
    integer, dimension(:), allocatable :: parts
    integer, dimension(:), allocatable :: partAssign
    integer, dimension(:), allocatable :: idx
    integer, dimension(:), allocatable :: subIarv
    integer, dimension(:), allocatable :: testIarv
    integer, dimension(:,:), allocatable :: struct, mystruct
    integer, dimension(:), allocatable :: fullNodesAssign
    integer, dimension(n_threads) :: recvbuffer, displs
    integer sendbuffer

    integer, dimension(2) :: send_count
    integer :: count_tag=300, id_tag=600

    integer(Zoltan_INT) :: global_part, length 
    integer(Zoltan_INT), dimension(:), allocatable :: partids,wgtidx
    real(Zoltan_FLOAT), dimension(:), allocatable :: partsizes
    real offset,weight0,weight8,weight15

    character*(20) filename
    character*(20) format_string
    character*(10) what
    integer ounit

!******************************************************************
!** control on number of process and number nodes
!******************************************************************

    call mpi_get_param

    if((neldi/n_threads) .lt. ngr) then
      if(my_id .eq. 0) then
        write(6,*)
        write(6,*)'error stop:'
        write(6,*)'n_threads, neldi, ngr',n_threads,neldi,ngr
        write(6,*)'the number of processes is too high compared with the number of nodes'
        write(6,*)'you can try a smaller number of processes'
      end if
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      stop
    end if

!******************************************************************
!** Initialize Zoltan
!******************************************************************

    ierr=Zoltan_Initialize(version)

!******************************************************************
!** Read graph from input file and distribute it 
!******************************************************************

    if(my_id .eq. 0) then

      ! what(partitioning type): on elements or on nodes
      what='elems'

      !Get the number of vertices
      allocate(ieltGrid%nNeighbor(neldi))

      call findNeighbor(neldi, ieltv, ieltGrid)

      !Get the number of vertices
      numGlobalVertices = neldi

      !Get the number of vertex neighbors
      numGlobalNeighbors = ieltGrid%totalNeighbor

      allocate(graph%vertexGID(numGlobalVertices))
      allocate(graph%nborIndex(numGlobalVertices + 1))
      allocate(graph%nborGID(numGlobalNeighbors))
      allocate(graph%nborProc(numGlobalNeighbors))

      graph%nborIndex(1) = 1
      do i=1,numGlobalVertices

        nnbors = ieltGrid%nNeighbor(i)       !number of neighbors to the vertex/node 
        graph%vertexGID(i) = i

        k=0
        do j=0,2
          if(ieltv(j+1,i) .ne. 0 .and. ieltv(j+1,i) .ne. -1) then
            graph%nborGID(graph%nborIndex(i) + j - k) = ieltv(j+1,i)
          else
            k = k+1
          end if
        end do

        graph%nborIndex(i+1) = graph%nborIndex(i) + nnbors

      end do

      !Assign each vertex to a process using a hash function
      do i=1,numGlobalNeighbors
        id = graph%nborGID(i)
        graph%nborProc(i) = mod(id, n_threads)
      end do

      !Create a sub graph for each process
      allocate(send_graph(n_threads), STAT = ierr)

      do i=1, numGlobalVertices
        id = graph%vertexGID(i)
        procID = mod(id, n_threads)+1
        send_graph(procID)%numMyVertices = send_graph(procID)%numMyVertices + 1
      end do

      allocate(idx(n_threads))

      do i=1, n_threads
        idx(i) = 0
        num = send_graph(i)%numMyVertices
        allocate(send_graph(i)%vertexGID(num))
        allocate(send_graph(i)%nborIndex(num + 1))
        do j=1, num + 1
          send_graph(i)%nborIndex(j) = 0
        end do
      end do

      do i=1, numGlobalVertices
        id = graph%vertexGID(i)
        nnbors = graph%nborIndex(i+1) - graph%nborIndex(i)
        procID = mod(id, n_threads)

        j = idx(procID+1)
                
        send_graph(procID+1)%vertexGID(j+1) = id
        send_graph(procID+1)%nborIndex(j+2) = send_graph(procID+1)%nborIndex(j+1) + nnbors

        idx(procID+1) = j+1
      end do

      do i=1,n_threads
        num = send_graph(i)%nborIndex(send_graph(i)%numMyVertices+1)
        allocate (send_graph(i)%nborGID(num))
        allocate (send_graph(i)%nborProc(num))

        send_graph(i)%numAllNbors = num
      end do

      do i=1,n_threads
        idx(i) = 0
      end do

      k=0
      do i=1,numGlobalVertices
        id = graph%vertexGID(i)
        nnbors = graph%nborIndex(i+1) - graph%nborIndex(i)
        procID = mod(id, n_threads)
        j = idx(procID+1)

        if (nnbors > 0) then
          do nbr=1,nnbors
            send_graph(procID+1)%nborGID(j+nbr) = graph%nborGID(i+nbr-1+k)
            send_graph(procID+1)%nborProc(j+nbr) = graph%nborProc(i+nbr-1+k)
          end do
          k=k+nnbors-1

          idx(procID+1) = j + nnbors
        endif
      end do

      deallocate (idx)


      ! Process zero sub-graph

      deallocate (graph%vertexGID)
      deallocate (graph%nborProc)

      graph = send_graph(1) !!!!!!!!!!!!!!!!

      ! Send other processes their subgraph
      do i=1, n_threads-1  
        send_count(1) = send_graph(i+1)%numMyVertices
        send_count(2) = send_graph(i+1)%numAllNbors

        call MPI_Send(send_count, 2, MPI_INTEGER, i, count_tag, MPI_COMM_WORLD, st, ierr)

        if (send_count(1) > 0) then

          call MPI_Send(send_graph(i+1)%vertexGID, send_count(1), MPI_INTEGER, i, id_tag, MPI_COMM_WORLD, st, ierr)

          call MPI_Send(send_graph(i+1)%nborIndex, send_count(1) + 1, MPI_INTEGER, i, id_tag + 1, MPI_COMM_WORLD, st, ierr)

          if (send_count(2) > 0) then
            call MPI_Send(send_graph(i+1)%nborGID, send_count(2), MPI_INTEGER, i, id_tag + 2, MPI_COMM_WORLD, st, ierr)

            call MPI_Send(send_graph(i+1)%nborProc, send_count(2), MPI_INTEGER, i, id_tag + 3, MPI_COMM_WORLD, st, ierr)
          end if
        end if
      end do

      do i=1,n_threads
        deallocate(send_graph(i)%nborIndex)
        deallocate(send_graph(i)%vertexGID)
        deallocate (send_graph(i)%nborGID)
        deallocate (send_graph(i)%nborProc)
      end do

    else ! rank /= 0
      call MPI_Recv(send_count, 2, MPI_INTEGER, 0, count_tag, MPI_COMM_WORLD, st, ierr)

      if (send_count(1) < 0) then
        call MPI_Finalize(ierr)
        stop
      endif

      graph%numMyVertices = send_count(1)
      graph%numAllNbors  = send_count(2)

      if (send_count(1) > 0) then
        allocate (graph%vertexGID(send_count(1)))
        allocate (graph%nborIndex((send_count(1) + 1)))

        if (send_count(2) > 0) then
          allocate (graph%nborGID(send_count(2)))
          allocate (graph%nborProc(send_count(2)))
        endif
      endif

      if (send_count(1) > 0) then

        call MPI_Recv(graph%vertexGID,send_count(1), MPI_INTEGER, 0, id_tag, MPI_COMM_WORLD, st, ierr)

        call MPI_Recv(graph%nborIndex,send_count(1) + 1, MPI_INTEGER, 0, id_tag + 1, MPI_COMM_WORLD, st, ierr)

        if (send_count(2) > 0) then
          call MPI_Recv(graph%nborGID,send_count(2), MPI_INTEGER, 0, id_tag + 2, MPI_COMM_WORLD, st, ierr)
          call MPI_Recv(graph%nborProc,send_count(2), MPI_INTEGER, 0, id_tag + 3, MPI_COMM_WORLD, st, ierr)
        endif
      endif

    endif

    call MPI_BCAST(numGlobalVertices, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)

    call MPI_BCAST(ilhv, numGlobalVertices, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)

    nullify(zz_obj)
    zz_obj => Zoltan_Create(MPI_COMM_WORLD)


    ! General parameters
    ierr=Zoltan_Set_Param(zz_obj, "DEBUG_LEVEL", "0")
    ierr=Zoltan_Set_Param(zz_obj, "LB_METHOD", "GRAPH")       !partitioning method
    ierr=Zoltan_Set_Param(zz_obj, "LB_APPROACH", "PARTITION") !partitioning approach
    ierr=Zoltan_Set_Param(zz_obj, "NUM_GID_ENTRIES", "1")     !global IDs are integers
    ierr=Zoltan_Set_Param(zz_obj, "NUM_LID_ENTRIES", "1")     !local IDs are integers
    ierr=Zoltan_Set_Param(zz_obj, "RETURN_LISTS", "EXPORT")      !export AND import lists
    ierr=Zoltan_Set_Param(zz_obj, "PHG_REFINEMENT_QUALITY", "10")      

    ! Graph parameters
    ierr=Zoltan_Set_Param(zz_obj, "CHECK_GRAPH", "1")
    ierr=Zoltan_Set_Param(zz_obj, "PHG_EDGE_SIZE_THRESHOLD", ".35")  ! 0-remove all, 1-remove none
    ierr=Zoltan_Set_Param(zz_obj, "OBJ_WEIGHT_DIM", "1")  ! 0-no weight, 1-weight on vertices
    

    ! Query functions - defined in simpleQueries.h (nel file di esempio)
    ierr=Zoltan_Set_Fn(zz_obj, ZOLTAN_NUM_OBJ_FN_TYPE, get_number_of_vertices)
    ierr=Zoltan_Set_Fn(zz_obj, ZOLTAN_OBJ_LIST_FN_TYPE, get_vertex_list)
    ierr=Zoltan_Set_Fn(zz_obj, ZOLTAN_NUM_EDGES_MULTI_FN_TYPE, get_num_edges_list)
    ierr=Zoltan_Set_Fn(zz_obj, ZOLTAN_EDGE_LIST_MULTI_FN_TYPE, get_edge_list)


    ierr = Zoltan_LB_Partition(zz_obj, changes, numGidEntries, numLidEntries, &
                               numImport, importGlobalGids, importLocalGids, importProcs, importToPart, &
                               numExport, exportGlobalGids, exportLocalGids, exportProcs, exportToPart)

    if(ierr .ne. ZOLTAN_OK) then
       write(6,*)'sorry...'
       call MPI_FINALIZE()
       stop
    end if

    allocate (parts(graph%numMyVertices+1))

    allocate(partAssign(numGlobalVertices))
    allocate(allPartAssign(numGlobalVertices))

    do i=1,graph%numMyVertices
      parts(i)=my_id
    end do

    do i=1,numExport
      parts(exportLocalGids(i)) = exportToPart(i)
    end do

    do i=1,numGlobalVertices 
       partAssign(i)=0
    end do

    j=1
    do i=1, graph%numMyVertices
      if(parts(i) .ne. my_id) then
        partAssign(exportGlobalGids(j)) = parts(i)
        j=j+1
      else
        partAssign(graph%vertexGID(i)) = parts(i)
      end if
    end do
  
    call MPI_ALLREDUCE(partAssign, allPartAssign, numGlobalVertices, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)

    ierr = Zoltan_LB_Free_Part(importGlobalGids, importLocalGids, importProcs, importToPart)
    ierr = Zoltan_LB_Free_Part(exportGlobalGids, exportLocalGids, exportProcs, exportToPart)

    call Zoltan_Destroy(zz_obj)



    partNel = 0
    do i=1,neldi
      if(allPartAssign(i) .eq. my_id) then
        partNel = partNel + 1
      end if
    end do

    call MPI_ALLREDUCE(partNel, minNel, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(partNel, maxNel, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)

    if(my_id .eq. 0) then
      write(6,*)'partitioning in',n_threads,'processes'
      write(6,*)'least number of elements for domain',minNel
      write(6,*)'greatest number of elements for domain',maxNel
    end if

    if(minNel .le. 0) then
      if(my_id .eq. 0) then
        write(6,*)'tried partitioning with',n_threads,' processes'
        write(6,*)'process with zero elements found:'
        write(6,*) 'trying with PHG COARSENING METHOD=IPM'
      end if
      repart = .true.
    else
      call shypart_mk_elem
    end if

    if(repart) then

      nullify(zz_obj2)
      zz_obj2 => Zoltan_Create(MPI_COMM_WORLD)
      !zz_obj2 => Zoltan_Copy(zz_obj)

      phg_method = 'IPM'

      ! General parameters
      ierr=Zoltan_Set_Param(zz_obj2, "DEBUG_LEVEL", "0")
      ierr=Zoltan_Set_Param(zz_obj2, "LB_METHOD", "GRAPH")       !partitioning method
      ierr=Zoltan_Set_Param(zz_obj2, "LB_APPROACH", "PARTITION") !partitioning approach
      ierr=Zoltan_Set_Param(zz_obj2, "NUM_GID_ENTRIES", "1")     !global IDs are integers
      ierr=Zoltan_Set_Param(zz_obj2, "NUM_LID_ENTRIES", "1")     !local IDs are integers
      ierr=Zoltan_Set_Param(zz_obj2, "RETURN_LISTS", "EXPORT")      !export AND import lists
      ierr=Zoltan_Set_Param(zz_obj2, "PHG_REFINEMENT_QUALITY", "10")      
      ierr=Zoltan_Set_Param(zz_obj2, "PHG_COARSENING_METHOD", "IPM")

      ! Graph parameters
      ierr=Zoltan_Set_Param(zz_obj2, "CHECK_GRAPH", "1")
      ierr=Zoltan_Set_Param(zz_obj2, "PHG_EDGE_SIZE_THRESHOLD", ".35")  ! 0-remove all, 1-remove none
      ierr=Zoltan_Set_Param(zz_obj2, "OBJ_WEIGHT_DIM", "1")  ! 0-no weight, 1-weight on vertices
    

      ! Query functions - defined in simpleQueries.h (nel file di esempio)
      ierr=Zoltan_Set_Fn(zz_obj2, ZOLTAN_NUM_OBJ_FN_TYPE, get_number_of_vertices)
      ierr=Zoltan_Set_Fn(zz_obj2, ZOLTAN_OBJ_LIST_FN_TYPE, get_vertex_list)
      ierr=Zoltan_Set_Fn(zz_obj2, ZOLTAN_NUM_EDGES_MULTI_FN_TYPE, get_num_edges_list)
      ierr=Zoltan_Set_Fn(zz_obj2, ZOLTAN_EDGE_LIST_MULTI_FN_TYPE, get_edge_list)

      deallocate(parts)
      deallocate(partAssign)
      deallocate(allPartAssign)

      ierr = Zoltan_LB_Partition(zz_obj2, changes2, numGidEntries2, numLidEntries2, &
                               numImport2, importGlobalGids2, importLocalGids2, importProcs2, importToPart2, &
                               numExport2, exportGlobalGids2, exportLocalGids2, exportProcs2, exportToPart2)

      if(ierr .ne. ZOLTAN_OK) then
        write(6,*)'sorry...'
        call MPI_FINALIZE()
        stop
      end if

      allocate (parts(graph%numMyVertices+1))

      allocate(partAssign(numGlobalVertices))
      allocate(allPartAssign(numGlobalVertices))

      do i=1,graph%numMyVertices
        parts(i)=my_id
      end do

      do i=1,numExport2
        parts(exportLocalGids2(i)) = exportToPart2(i)
      end do

      do i=1,numGlobalVertices 
         partAssign(i)=0
      end do

      j=1
      do i=1, graph%numMyVertices
        if(parts(i) .ne. my_id) then
          partAssign(exportGlobalGids2(j)) = parts(i)
          j=j+1
        else
          partAssign(graph%vertexGID(i)) = parts(i)
        end if
      end do
  
      call MPI_ALLREDUCE(partAssign, allPartAssign, numGlobalVertices, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)

      ierr = Zoltan_LB_Free_Part(importGlobalGids2, importLocalGids2, importProcs2, importToPart2)
      ierr = Zoltan_LB_Free_Part(exportGlobalGids2, exportLocalGids2, exportProcs2, exportToPart2)

      call Zoltan_Destroy(zz_obj2)

      partNel = 0
      do i=1,neldi
        if(allPartAssign(i) .eq. my_id) then
          partNel = partNel + 1
        end if
      end do

      call MPI_ALLREDUCE(partNel, minNel, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(partNel, maxNel, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)

      if(my_id .eq. 0) then
        write(6,*)'partitioning in',n_threads,'processes'
        write(6,*)'least number of elements for domain',minNel
        write(6,*)'greatest number of elements for domain',maxNel
      end if

      if(minNel .le. 0) then
        if(my_id .eq. 0) then
          write(6,*)'tried partitioning with',n_threads,' processes'
          write(6,*)'too many processes for the basin, try a smaller number of processes'
        end if
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        stop
      end if

      call shypart_mk_elem

    end if


    if(my_id .eq. 0) then
      deallocate(ieltGrid%nNeighbor)
      deallocate(send_graph)
    end if
    deallocate (graph%vertexGID)
    deallocate (graph%nborIndex)
    deallocate (graph%nborGID)
    deallocate (graph%nborProc)
    deallocate(parts)
    deallocate(partAssign)

    if(my_id .eq. 0) then
      if(n_threads .gt. 999) then
        format_string = "(A10,I4)"
        write(filename,format_string)'part_elems',n_threads
      else if(n_threads .gt. 99) then
        format_string = "(A10,I3)"
        write(filename,format_string)'part_elems',n_threads
      else if(n_threads .gt. 9) then
        format_string = "(A10,I2)"
        write(filename,format_string)'part_elems',n_threads
      else
        format_string = "(A10,I1)"
        write(filename,format_string)'part_elems',n_threads
      end if

      open(unit=1500, file=filename, action='write')
        write(1500,*),nkndi,neldi,n_threads,what
        write(1500,fmt="(i12,i12,i12,i12,i12,i12)"),allPartAssign
      close(1500)
    end if

    return

end subroutine

