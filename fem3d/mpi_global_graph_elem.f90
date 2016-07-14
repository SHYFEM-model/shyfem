
!#######################################################################!
!*****************  start module Global_Graph_Data_Ele  ****************!
!#######################################################################!

module Global_Graph_Data_Ele

  use mpi
  use zoltan
  use basin
  use shypart
  use communicationStruct
   ! Structure to hold graph

   ! MPI variables rankID and processes number
   !integer, save, public :: my_id, n_threads
   !INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: allPartAssign

   type GRAPH_DATA
     integer :: numMyVertices  ! total vertices in in my partition
     integer :: numAllNbors   ! total number of neighbors of my vertices
     integer(ZOLTAN_INT), dimension(:), allocatable :: vertexGID    ! global ID of each of my vertices
     integer, dimension(:),allocatable :: nborIndex    ! nborIndex[i] is location of start of neighbors for vertex i
     integer(ZOLTAN_INT), dimension(:), allocatable :: nborGID      ! nborGIDs[nborIndex[i]] is first neighbor of vertex i
     integer, dimension(:),allocatable :: nborProc     ! process owning each nbor in nborGID
   end type GRAPH_DATA

   type IELTV_GRID
     integer :: totalNeighbor   ! total number of neighbors
     integer, dimension(:), allocatable :: nNeighbor  ! array containing the number of neighbor for each element
   end type IELTV_GRID     

   !! Zoltan data to store in module
   LOGICAL :: changes
   INTEGER(Zoltan_INT) :: numGidEntries, numLidEntries
   INTEGER(Zoltan_INT) :: numImport, numExport
   INTEGER(Zoltan_INT), POINTER, DIMENSION(:) :: importGlobalGids, exportGlobalGids
   INTEGER(Zoltan_INT), POINTER, DIMENSION(:) :: importLocalGids, exportLocalGids
   INTEGER(Zoltan_INT), POINTER, DIMENSION(:) :: importProcs, exportProcs
   INTEGER(Zoltan_INT), POINTER, DIMENSION(:) :: importToPart, exportToPart

   LOGICAL :: changes2
   INTEGER(Zoltan_INT) :: numGidEntries2, numLidEntries2
   INTEGER(Zoltan_INT) :: numImport2, numExport2
   INTEGER(Zoltan_INT), POINTER, DIMENSION(:) :: importGlobalGids2, exportGlobalGids2
   INTEGER(Zoltan_INT), POINTER, DIMENSION(:) :: importLocalGids2, exportLocalGids2
   INTEGER(Zoltan_INT), POINTER, DIMENSION(:) :: importProcs2, exportProcs2
   INTEGER(Zoltan_INT), POINTER, DIMENSION(:) :: importToPart2, exportToPart2

   ! Declare a global Mesh data structure. 
   type(GRAPH_DATA) :: graph,graph2

contains

!######################################################################!
!***************  start get_number_of_vertices function  **************!
!######################################################################!

  ! Application defined query functions
  function get_number_of_vertices(mygraph, ierr) result (nVert)

    implicit none

    ! in
    integer(Zoltan_INT) :: mygraph(1)

    ! out
    integer, intent(out):: ierr
    integer :: nVert

    ierr = ZOLTAN_OK
    nVert = graph%numMyVertices
  end function get_number_of_vertices

!######################################################################!
!******************  start get_vertex_list subroutine  ****************!
!######################################################################!

  subroutine get_vertex_list(mygraph, sizeGID, sizeLID, globalID, localID, wgt_dim, obj_wgts, ierr)

    use levels

    implicit none

    ! in
    integer(Zoltan_INT), intent(in) :: mygraph(1)
    integer(Zoltan_INT), intent(in) :: sizeGID, sizeLID, wgt_dim
    
    ! out
    integer(Zoltan_INT), intent(out), dimension(*) :: globalID, localID  
    real(Zoltan_FLOAT), dimension(*) :: obj_wgts
    !real(Zoltan_FLOAT), intent(out), dimension(*) :: obj_wgts
    integer(Zoltan_INT), intent(out) :: ierr
    
    ! local
    integer :: i
   
    ierr = ZOLTAN_OK

    ! In this example, return the IDs of our vertices, but no weights.
    ! Zoltan will assume equally weighted vertices.

    do i=1,graph%numMyVertices
      globalID(i) = graph%vertexGID(i)
      localID(i) = i
      !obj_wgts(i) = 0.650528 + 0.349472 * ilhv(globalID(i))
      !obj_wgts(i) = 0.92 + 0.08 * ilhv(globalID(i))
      !obj_wgts(i) = 0.88 + 0.12 * ilhv(globalID(i))
      !obj_wgts(i) = 0.85 + 0.15 * ilhv(globalID(i))
      !obj_wgts(i) = 0.75 + 0.25 * ilhv(globalID(i))
      !obj_wgts(i) = 0.81 + 0.19 * ilhv(globalID(i))
      !obj_wgts(i) = 1.35 + 0.15 * ilhv(globalID(i))
      obj_wgts(i) = 1.28 + 0.15 * ilhv(globalID(i))

    end do
  end subroutine get_vertex_list

!######################################################################!
!****************  start get_num_edges_list subroutine  ***************!
!######################################################################!

  subroutine get_num_edges_list(mygraph, sizeGID, sizeLID, num_obj, globalID, localID, numEdges, ierr)

    implicit none
    
    ! in
    integer(Zoltan_INT) :: mygraph(1)
    integer(Zoltan_INT), intent(in) :: sizeGID, sizeLID, num_obj
    integer(Zoltan_INT), dimension(*), intent(in) :: globalID,localID
    
    ! out
    integer(Zoltan_INT), dimension(*), intent(out) :: numEdges
    integer(Zoltan_INT), intent(out) :: ierr

    ! local
    integer :: i, idx

    if ((sizeGID /= 1) .or. (sizeLID /= 1) .or. (num_obj /= graph%numMyVertices)) then
      ierr = ZOLTAN_FATAL
      return
    end if
    
    do i=1,num_obj
      idx = localID(i)
      numEdges(i) = graph%nborIndex(idx+1) - graph%nborIndex(idx)
    end do

    ierr = ZOLTAN_OK
  end subroutine get_num_edges_list

!######################################################################!
!*******************  start get_edge_list subroutine  *****************!
!######################################################################!

  subroutine get_edge_list(mygraph,sizeGID,sizeLID,num_obj,globalID,localID,num_edges,nborGID,nborProc,wgt_dim,ewgts,ierr)

    implicit none

    !in
    integer(Zoltan_INT), intent(in) :: mygraph(1)
    integer(Zoltan_INT), intent (in) :: sizeGID, sizeLID, num_obj, wgt_dim
    integer(Zoltan_INT), dimension(*), intent(in) :: globalID,localID, num_edges
    
    !out
    integer(Zoltan_INT), dimension(*), intent(out) :: nborGID, nborProc
    real(Zoltan_FLOAT), dimension(*) :: ewgts
    integer(Zoltan_INT), intent (out) :: ierr
    
    !local
    integer i, j, from, to

    ierr = ZOLTAN_OK

	if ((sizeGID /= 1).or.(sizeLID /= 1).or.(num_obj /= graph%numMyVertices).or.(wgt_dim /= 0)) then
	    ierr = ZOLTAN_FATAL
	    return
	end if

    do i=1, num_obj

      ! In this example, we are not setting edge weights.
      ! Zoltan will set each edge to weight 1.0.

      to = graph%nborIndex(localID(i)+1)
      from = graph%nborIndex(localID(i))
      if ((to - from) /= num_edges(i)) then
        ierr = ZOLTAN_FATAL
        exit
      end if
      do j=from+1, to
        nborGID(j) = graph%nborGID(j)
        nborProc(j) = graph%nborProc(j)
      end do
    end do
    
  end subroutine get_edge_list

  subroutine mpi_get_param

    use basin
    use mod_geom

    implicit none

    integer param(10)
    integer ierr

    if(my_id .eq. 0) then
      param(1) = nkndi
      param(2) = neldi
      param(3) = ngr
      param(4) = mbw
      param(5) = nkn
    end if

    call MPI_BCAST(param, 5, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)  

    if(my_id .ne. 0) then
      nkndi = param(1)
      neldi = param(2)
      ngr = param(3)
      mbw = param(4)
      nkn = param(5)
    end if 

 

  end subroutine mpi_get_param


!######################################################################!
!*******************  start findNeighbor subroutines  *****************!
!######################################################################!
!######################################################################!
!** This subroutine find the number of neighbors for each element and *!
!** inserts it into the array nNeighbor of the ieltvGrid structure  ***!
!######################################################################!

  subroutine findNeighbor(nel, ieltv, ieltGrid)

    implicit none

    ! arguments
    integer nel
    integer ieltv(3,1)
    type (IELTV_GRID) :: ieltGrid

    ! local
    integer ie, ii, ien
    logical bverbose
    integer i,j       
 
!-------------------------------------------------------------------
! initialize
!-------------------------------------------------------------------

    bverbose = .false.
    ieltGrid%totalNeighbor = 0

!-------------------------------------------------------------------
! loop over elements
!-------------------------------------------------------------------

    do ie=1,nel
      ieltGrid%nNeighbor(ie) = 0  ! initializes ieltGrid%nNeighbor(ie)

      do ii=1,3

        ien = ieltv(ii,ie)

        if( ien .gt. 0 ) then

          if( ien .gt. nel ) then
            write(6,*) 'ie,ii,ien,nel: ',ie,ii,ien,nel
            stop 'error stop findNeighbor: corrupt data structure of ieltv'
          end if

          ieltGrid%totalNeighbor = ieltGrid%totalNeighbor + 1
          ieltGrid%nNeighbor(ie) = ieltGrid%nNeighbor(ie) + 1

        end if

      end do

    end do

    if( bverbose ) then
      write(6,*) 'findNeighbor is ok'
      write(6,*) '  total neighbors =      ',ieltGrid%totalNeighbor
    end if

!-------------------------------------------------------------------
! end of routine
!-------------------------------------------------------------------

    return

  end subroutine


!*****************************************************************

	subroutine shypart_mk_elem

        use basin
        !use shympi
        !use Graph_Ele
        use mod_geom
        use commonStructures
        !use communicationStruct
        use levels

	implicit none

        integer i,nlk,ierr
        integer, dimension(:,:), allocatable :: mystruct
        integer h
        integer temp_nkn, temp_nel

        nel = 0
        do i=1,neldi
          if(allPartAssign(i) .eq. my_id) then
            nel = nel + 1
          end if
        end do

        write(6,*)'myElements =',nel,neldi,my_id

        call makeMyEle2(nel,neldi,nen3v)

        allocate(mystruct(3,myele%numberID))

        call makeMyStruct(myele, nen3v, mystruct, 3)

        call makeMyNodes2(myele, mystruct, mynodes, nen3v, nkndi)

        deallocate(mystruct)
        
        nkn =mynodes%numberID
        nel =myele%numberID
 
        if(nel .eq. 0) then
          ngr = 0
        end if


        temp_nkn=nkndi
        temp_nel=neldi
        call transfer_domain
        neldi=temp_nel
        nkndi=temp_nkn

        nlk = 3*nel+2*nkn

        call mod_geom_init(nkn,nel,ngr)

        if( ngr .gt. maxlnk )then
          write(6,*) 'ngr,maxlnk: ',ngr,maxlnk
          write(6,*) 'Please adjust maxlnk in links.h'
          stop 'error stop shypart_mk_elem: maxlnk'
        end if

!-------------------------------------------------------------
! make static arrays
!-------------------------------------------------------------

        repart=.false.
        call mklenk(nlk,nkn,nel,nen3v,ilinkv,lenkv,repart)
        call mklenkii(nlk,nkn,nel,nen3v,ilinkv,lenkv,lenkiiv,repart)
        call mklink(nkn,ilinkv,lenkv,linkv,repart)

        call mkkant(nkn,ilinkv,lenkv,linkv,kantv,repart)
        call mkielt(nkn,nel,ilinkv,lenkv,linkv,ieltv,repart)

        call MPI_ALLREDUCE(repart, repart, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)

        if((repart) .and. (phg_method .eq. 'AGG')) then
          nel = neldi
          nkn = nkndi
          deallocate(nen3v)
          !deallocate(hm3v)
          !deallocate(ipv)
          !deallocate(ipev)
          !deallocate(iarv)
          !deallocate(iarnv)
          !deallocate(xgv)
          !deallocate(ygv)
          if(my_id .eq. 0) deallocate(allNodesAssign)
          deallocate(mypart%mysend%node_send)
          deallocate(univocal_nodes%localID)
          deallocate(univocal_nodes%globalID)
          deallocate(numberElements)
          deallocate(numberNodes)
          deallocate(myele%globalID)
          deallocate(mynodes%globalID)
          allocate(nen3v(3,neldi))
          !allocate(hm3v(3,neldi))
          !allocate(ipv(nkndi))
          !allocate(ipev(neldi))
          !allocate(iarv(neldi))
          !allocate(iarnv(nkndi))
          !allocate(xgv(nkndi))
          !allocate(ygv(nkndi))
          nen3v=tempNen3v
          !hm3v=tempHm3v
          !ipv=tempIpv
          !ipev=tempIpev
          !iarv=tempIarv
          !iarnv=tempIarnv
          !xgv=tempXgv
          !ygv=tempYgv
        
          call mod_geom_init(nkndi,neldi,ngr)
          ieltv=temp_ieltv
          ilhv=temp_ilhv
          linkv=temp_linkv
          lenkv=temp_lenkv
          ilinkv=temp_ilinkv
          !lenkiiv=temp_lenkiiv
          !kantv=temp_kantv
          
          !call shypart_setup
          write(6,*)'recall partPHG',my_id
          return
          !call partPHG
        else if((repart) .and. (phg_method .eq. 'IPM')) then
          if(my_id .eq. 0) then
            write(6,*)'tried partitioning with',n_threads,' processes'
            write(6,*)'static arrays errors occurred for the basin, try a smaller number of processes'
          end if
          call MPI_Barrier(MPI_COMM_WORLD,ierr)
          stop
        end if

	end subroutine shypart_mk_elem


!#######################################################################!
!*********************  start makemyele subroutine  ********************!
!#######################################################################!
!#######################################################################!
!* This subroutine makes the myele structure. It contains informations *!
!* on the elements of a subprocess(global ID elements, rank ID process *!
!* neighbor, globalID neighbor elements, etc.)                         *!
!#######################################################################!
!######################################################################################################!
! n = local ID of the Elements number(numMyVertices)                                                   !
!                        <---------- n --------->                                                      !
! myele%globalID         |||||||||||||||||||||||| globalID of the element                              !
! myele%rankID(1,*)      |||||||||||||||||||||||| contains the rankID of the process to export data    !   
! myele%rankID(2,*)      |||||||||||||||||||||||| contains the rankID of the process to export data    !  
! myele%rankID(3,*)      |||||||||||||||||||||||| contains the rankID of the process to export data    !  
! myele%neighborID(1,*)  |||||||||||||||||||||||| contains the global neighbor ID to which import data !
! myele%neighborID(2,*)  |||||||||||||||||||||||| contains the global neighbor ID to which import data !
! myele%neighborID(3,*)  |||||||||||||||||||||||| contains the global neighbor ID to which import data !
!######################################################################################################!

  subroutine makeMyEle2(numMyVertices, numGlobalVertices, struct)

    use commonStructures
    use mod_geom

    implicit none

    ! input
    ! rank process, number of my elements, number global elements, number global of neighbor
    integer numMyVertices, numGlobalVertices

    ! allPartAssign contains the ID of the process to which is assigned the element
    ! the number of the neighbors an element "nneighbor" is equal to
    integer, dimension(3,numGlobalVertices) :: struct 
!    integer, dimension(27*numMyVertices) :: neighbor, sort_neighbor   
    integer, dimension(:), allocatable :: neighbor, sort_neighbor   
    integer, dimension(:), allocatable :: sortRank, rank   
    integer :: st(MPI_STATUS_SIZE), ierr
    integer error, STAT
 
    ! local 
    integer eleneighbor, myneighbor, nneighbor, globalID, locID, neighbID
    integer x,p,s,i,j,k,h,n,ele,t,length,nprocesses,node,maxSends,maxReceives,sizev


    sizev = (ngr-2)*3*numMyVertices
    allocate(neighbor(ngr*3*numMyVertices),STAT=error)
        if (error .ne. 0) then
          write(6,*)'error: could not allocate memory for array',error
          stop
        endif
    allocate(sort_neighbor(ngr*3*numMyVertices),STAT=error)
        if (error .ne. 0) then
          write(6,*)'error: could not allocate memory for array',error
          stop
        endif
    allocate(rank(sizev),STAT=error)
        if (error .ne. 0) then
          write(6,*)'error: could not allocate memory for array',error
          stop
        endif
    allocate(sortRank(sizev),STAT=error)
        if (error .ne. 0) then
          write(6,*)'error: could not allocate memory for array',error
          stop
        endif

    do i=1,sizev
       rank(i) = -1
       sortRank(i) = 0
    end do

    do i=1, size(neighbor)
      neighbor(i) = 0
    end do

    myneighbor = 0
    k = 0
    myele%numberID=0
    do i=1,numGlobalVertices
       if(allPartAssign(i) .eq. my_id) then

          do h=1,3
            n = struct(h,i)
            eleneighbor=ilinkv(n+1)-ilinkv(n)
            do j=1,eleneighbor
              ele = lenkv(ilinkv(n) + j)
              if (ele .ne. 0) then
                if(allPartAssign(ele) .ne. my_id) then
                  k = k + 1
                  rank(k) = allPartAssign(ele)
                  myneighbor = myneighbor + 1
                  neighbor(myneighbor) = ele
                end if
              end if
            end do
          end do

          myele%numberID=myele%numberID+1
       end if
    end do

    if(myele%numberID .ne. numMyVertices)then
      write(*,*)'Error stop makemyele'
      stop
    end if

    write(6,*)'numberElements',my_id

    allocate(numberElements(n_threads),STAT=error)
        if (error .ne. 0) then
          write(6,*)'error: could not allocate memory for array',error
          stop
        endif

    call MPI_ALLGATHER(myele%numberID, 1, MPI_INTEGER, numberElements, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    if(myneighbor .gt. 0)then
      call remove_dups(neighbor,length,sort_neighbor)   
    else
      length = 0
    end if

    deallocate(neighbor)

    if(k .gt. 0)then
      call remove_dups0(rank,nprocesses,sortRank)   
    else
      nprocesses = 0
    end if

    deallocate(rank)

    write(6,*)'2numberElements',my_id
    if( .not. allocated(mypart%mysend%process)) then
      allocate(mypart%mysend%process(nprocesses),STAT=error)
    else
      deallocate(mypart%mysend%process)
      allocate(mypart%mysend%process(nprocesses),STAT=error)
    end if
    if( .not. allocated(mypart%myreceive%process)) then
      allocate(mypart%myreceive%process(nprocesses),STAT=error)
    else
      deallocate(mypart%myreceive%process)
      allocate(mypart%myreceive%process(nprocesses),STAT=error)
    end if

    mypart%mysend%sends = nprocesses
    mypart%myreceive%receives = nprocesses

    do i=1,nprocesses
       mypart%mysend%process(i) = sortRank(i)
       mypart%myreceive%process(i) = sortRank(i)
    end do

    deallocate(sortRank)

    if(myele%numberID .ne. 0) then
      allocate(myele%globalID(myele%numberID),STAT=error)
        if (error .ne. 0) then
          write(6,*)'error: could not allocate memory for array',error
          stop
        endif
    end if

    locID=1
    do i=1,numGlobalVertices
       if(allPartAssign(i) .eq. my_id) then
          myele%globalID(locID) = i
          locID = locID + 1
       end if
    end do

    deallocate(sort_neighbor)

    return

  end subroutine


!#######################################################################!
!********************  start makemynodes subroutine  *******************!
!#######################################################################!
!#######################################################################!
! This subroutine makes the mynodes structure. It contains informations !
!***  on the nodes of a subprocess(global ID nodes, rank ID process  ***!
!*****          neighbor, globalID neighbor nodes, etc.)           *****!
!#######################################################################!
!########################################################################################################!
! n = local ID of the node number                                                                        !
!                          <---------- n --------->                                                      !
! mynodes%globalID         |||||||||||||||||||||||| globalID of the node                                 !
! mynodes%rankID(1,*)      |||||||||||||||||||||||| contains the rankID of the process to export data    !   
! mynodes%rankID(2,*)      |||||||||||||||||||||||| contains the rankID of the process to export data    !  
! mynodes%rankID(3,*)      |||||||||||||||||||||||| contains the rankID of the process to export data    !  
! mynodes%neighborID(1,*)  |||||||||||||||||||||||| contains the global neighbor ID to which import data !
! mynodes%neighborID(2,*)  |||||||||||||||||||||||| contains the global neighbor ID to which import data !
! mynodes%neighborID(3,*)  |||||||||||||||||||||||| contains the global neighbor ID to which import data !
!########################################################################################################!

  subroutine makeMyNodes2(myele, mystruct, mynodes, struct, nkn)

    use commonStructures
    use mod_geom
    use levels
    implicit none

    ! input
    ! rank process, number of my elements

    integer nkn
    integer mystruct(3,1)
    integer, dimension(:,:) :: struct
    integer, dimension(:), allocatable :: neighbor, sort_neighbor, sort_nodes
    integer, dimension(:), allocatable :: fullNodesAssign, discard_nodes
    integer, dimension(n_threads) :: recvbuffer, displs
    type (COMMUNICATION_INFO) :: myele,mynodes,temp
    integer :: st(MPI_STATUS_SIZE), ierr 
    ! local 
    integer i, j, k, n, x, h, s, length, sendbuffer, node,nneighbor
    integer counter,nodiTotali
    integer, dimension(:), allocatable :: temp_discard_node
    logical debug, take_nodes
    integer, dimension(:), allocatable :: control

    !
    integer ounit,error
    character*(20) filename
    character*(20) format_string
    integer, dimension(:), allocatable :: nnodes
    integer, dimension(:), allocatable :: nnodesAssign
    integer sendbuffer2
    integer, dimension(n_threads) :: recvbuffer2, displs2
    integer maxNodes

    integer, dimension(:,:),allocatable :: data_recv
    integer, dimension(:),allocatable :: sreqArray
    integer, dimension(:),allocatable :: rreqArray

    if(.not. allocated(sreqArray)) then
      allocate(sreqArray(mypart%mysend%sends))
    else
      deallocate(sreqArray)
      allocate(sreqArray(mypart%mysend%sends))
    end if
    if(.not. allocated(rreqArray)) then
      allocate(rreqArray(mypart%myreceive%receives))
    else
      deallocate(rreqArray)
      allocate(rreqArray(mypart%myreceive%receives))
    end if
    

    debug = .false.

    if(.not. allocated(temp%globalID)) then
      allocate(temp%globalID(3*myele%numberID))
    else
      deallocate(temp%globalID)
      allocate(temp%globalID(3*myele%numberID))
    end if
    if(.not. allocated(sort_nodes)) then
      allocate(sort_nodes(3*myele%numberID))
    else
      deallocate(sort_nodes)
      allocate(sort_nodes(3*myele%numberID))
    end if

    n=0
    do i=1, myele%numberID
      do j=1,3
        n = n + 1
        temp%globalID(n) = mystruct(j,i)
      end do
    end do

    if(n .gt. 0)then
      call remove_dups(temp%globalID,length,sort_nodes)
    else
      length = 0
    end if

    mynodes%numberID = length

    write(6,*)'myNodes =',mynodes%numberID,nkndi,my_id

    call MPI_ALLREDUCE(mynodes%numberID, maxNodes, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)

    allocate(mynodes%globalID(mynodes%numberID))
 
    do i=1,mynodes%numberID
      mynodes%globalID(i) = sort_nodes(i)
    end do

    allocate(numberNodes(n_threads))

    call MPI_ALLGATHER(mynodes%numberID, 1, MPI_INTEGER, numberNodes, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    call spread_nodes_to_neighbor(sreqArray,rreqArray,data_recv,maxNodes)

    mypart%mysend%maxItems=maxNodes
    mypart%myreceive%maxItems=maxNodes

    totalNodes = numberNodes(1)
    do i=2,n_threads
      totalNodes = totalNodes + numberNodes(i)
    end do

    if(my_id .eq. 0)then
      write(6,*)'totalNodes=',totalNodes
    end if

    if(.not. allocated(allNodesAssign)) then
      allocate(allNodesAssign(n_threads,nkndi))
    else
      deallocate(allNodesAssign)
      allocate(allNodesAssign(n_threads,nkndi))
    end if

    if(.not. allocated(fullNodesAssign)) then
      allocate(fullNodesAssign(totalNodes))
    else
      deallocate(fullNodesAssign)
      allocate(fullNodesAssign(totalNodes))
    end if

    do i=1,n_threads
      recvbuffer(i) = numberNodes(i)
    end do

    sendbuffer = numberNodes(my_id+1)

    displs(1) = 0
    do i=2,n_threads
      displs(i) = displs(i-1) + numberNodes(i-1)
    end do

    call MPI_ALLGATHERV(sort_nodes, sendbuffer, MPI_INTEGER, fullNodesAssign, recvbuffer, displs, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    call makeNodesAssign(fullNodesAssign, nkn, n_threads)
    deallocate(fullNodesAssign)

    allocate(temp_discard_node(3*myele%numberID))

    do i=1,3*myele%numberID
       temp_discard_node(i) = 0
    end do

    allocate(mypart%mysend%node_temp(mypart%mysend%sends))
    do h=1,mypart%mysend%sends
       allocate(mypart%mysend%node_temp(h)%items(mynodes%numberID))
       do i=1,mynodes%numberID
          mypart%mysend%node_temp(h)%items(i)=-1
       end do
       mypart%mysend%node_temp(h)%numItems=0
    end do

    call waitAny(sreqArray)
    call waitAny(rreqArray)
    
    deallocate(sreqArray)
    deallocate(rreqArray)

    if(myele%numberID .gt. 0) then
       !allocate(control(mypart%myreceive%receives))
       !control =0
       counter =0
       do h=1,mypart%myreceive%receives
          !take_nodes = .false.
          do i=1,mynodes%numberID
             do n=1,numberNodes(mypart%myreceive%process(h)+1)
                if(mynodes%globalID(i) .eq. data_recv(n,h) ) then
                   mypart%mysend%node_temp(h)%numItems=mypart%mysend%node_temp(h)%numItems+1
                   mypart%mysend%node_temp(h)%items(mypart%mysend%node_temp(h)%numItems) = i

                   !if( (mod(control(h), 2) .eq. 0) .and. (my_id .gt. mypart%myreceive%process(h)) ) then
                   !   take_nodes = .true.
                   !else if( mod(control(h), 2) .eq. 1 .and. (my_id .lt. mypart%myreceive%process(h)) ) then
                   !   take_nodes = .true.
                   !end if
                   !control(h) = control(h) + 1

                   !if( take_nodes ) then
                   !   counter = counter + 1
                   !   temp_discard_node(counter) = mynodes%globalID(i)
                   !   take_nodes = .false.
                   !end if

                   if(my_id .gt. mypart%myreceive%process(h)) then
                      counter = counter + 1
                      temp_discard_node(counter) = mynodes%globalID(i)
 !                     discard_nodes(counter) = mynodes%globalID(i)
                   end if

                end if
             end do
          end do
       end do
    else
      counter = 0
    end if
    
    allocate(mypart%mysend%node_send(mypart%mysend%sends))
    do h=1,mypart%myreceive%receives
      mypart%mysend%node_send(h)%numItems=mypart%mysend%node_temp(h)%numItems
      if(.not. allocated(mypart%mysend%node_send(h)%items))then
        allocate(mypart%mysend%node_send(h)%items(mypart%mysend%node_send(h)%numItems))
      else
        deallocate(mypart%mysend%node_send(h)%items)
        allocate(mypart%mysend%node_send(h)%items(mypart%mysend%node_send(h)%numItems))
      end if
      do j=1,mypart%mysend%node_send(h)%numItems
        mypart%mysend%node_send(h)%items(j) = mypart%mysend%node_temp(h)%items(j)
      end do
    end do

    deallocate(data_recv)

    allocate(univocal_nodes%localID(mynodes%numberID))
    allocate(univocal_nodes%globalID(mynodes%numberID))

    allocate(discard_nodes(3*myele%numberID))


    if(counter .gt. 0)then
      call remove_dups(temp_discard_node,length,discard_nodes)
    else
      length = 0
    end if

    !length = counter

    counter = 0
    if(length .gt. 0) then
       n=1
       do i=1,mynodes%numberID
          do j=1,length
             if(discard_nodes(n) .eq. mynodes%globalID(i)) then
                n=n+1
                exit
             end if
          end do 
          if(j .eq. (length+1))then
             counter = counter + 1
             univocal_nodes%localID(counter) = i
             univocal_nodes%globalID(counter) = mynodes%globalID(i)
          end if
       end do
    else
       do i=1,mynodes%numberID
          counter = counter + 1
          univocal_nodes%localID(counter) = i
          univocal_nodes%globalID(counter) = mynodes%globalID(i)
       end do 
    end if

    univocal_nodes%numberID = counter

    !write(6,*),mynodes%totalID-mynodes%numberID
    !write(6,*),mypart%mysend%sends


!    do h=1,mypart%myreceive%receives
!       do i=1,mypart%mysend%node_send(h)%numItems
!          write(ounit,*),mypart%mysend%node_send(h)%items(i),mypart%mysend%process(h)
!          write(ounit,*),mypart%mysend%node_send(h)%numItems
!       end do
!    end do
!    close(ounit)

!     call sleep(8)

    debug = .false.

    if(debug) then

       if(my_id .gt. 99) then
         format_string = "(A9,I3)"
         write(filename,format_string)'globalID_',my_id
       else if(my_id .gt. 9) then
         format_string = "(A9,I2)"
         write(filename,format_string)'globalID_',my_id
       else
         format_string = "(A9,I1)"
         write(filename,format_string)'globalID_',my_id
       end if
   
       ounit=10+my_id 
       open(unit=ounit, file=filename, action='write')

       !write(ounit,*),myele%globalID
       !write(ounit,*),mynodes%globalID
       !write(ounit,*),univocal_nodes%globalID
       write(ounit,*)'send_receive',mypart%mysend%maxItems,mypart%myreceive%maxItems
       write(ounit,*)'nodes',numberNodes,totalNodes
       write(ounit,*)'ele_GID',myele%globalID
       write(ounit,*)'node_GID',mynodes%globalID
       write(ounit,*)'univocal_LID',univocal_nodes%numberID
       write(ounit,*)'unvocal_GID',univocal_nodes%globalID
       write(ounit,*)'send',mypart%mysend%process
       write(ounit,*)'receive',mypart%myreceive%process
       do h=1,mypart%mysend%sends
          write(ounit,*)'id,h,sends_num',my_id,h,mypart%mysend%sends
          write(ounit,*)'id,neighbor,numItems',my_id,mypart%mysend%process(h),mypart%mysend%node_send(h)%numItems
          write(ounit,*)'items',mypart%mysend%node_send(h)%items
          !write(ounit,*),mypart%mysend%node_receive(h)%numItems,mypart%mysend%node_receive(h)%items
       end do
       if(my_id .eq. 0)then
         write(ounit,*)'allNodesAssign',allNodesAssign
       end if

       write(ounit,*)length,counter,length+counter,mynodes%numberID,my_id

       call MPI_ALLREDUCE(univocal_nodes%numberID, nodiTotali, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

       write(ounit,*),counter,mynodes%numberID, nodiTotali, my_id

       allocate(nnodes(n_threads))
       allocate(nnodesAssign(nodiTotali))

       call MPI_AllGATHER(univocal_nodes%numberID, 1, MPI_INTEGER, nnodes, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

       do i=1,n_threads
          recvbuffer2(i) = nnodes(i)
       end do

       sendbuffer2 = nnodes(my_id+1)

       displs2(1) = 0
       do i=2,n_threads
         displs2(i) = displs2(i-1) + nnodes(i-1)
       end do

       call MPI_ALLGATHERV(univocal_nodes%globalID,sendbuffer2,MPI_INTEGER, &
                nnodesAssign,recvbuffer2,displs2,MPI_INTEGER,MPI_COMM_WORLD,ierr)

       do i=1,n_threads
          if( i .ne. (my_id+1)) then
             do j=1,univocal_nodes%numberID
                do k=1, nnodes(i)
                   if(univocal_nodes%globalID(j) .eq. nnodesAssign(displs2(i)+k)) then
                      write(ounit,*)'error: univocal_nodes',my_id
                      write(ounit,*),univocal_nodes%globalID(j),i-1,my_id
                   end if
                end do
             end do
          end if
       end do

       close(ounit)

    end if

    if(my_id .ne. 0) then
      deallocate(allNodesAssign)
    end if

    do h=1,mypart%mysend%sends
       deallocate(mypart%mysend%node_temp(h)%items)
    end do

    !deallocate(control)
    deallocate(mypart%mysend%node_temp)
    deallocate(discard_nodes)
    deallocate(temp_discard_node)
    deallocate (sort_nodes)
    deallocate(temp%globalID)

    return

  end subroutine

  
  subroutine makeNodesAssign(fullNodesAssign, nkn, n_threads)

    use commonStructures

    implicit none

    integer, dimension(totalNodes) :: fullNodesAssign
!    integer, dimension(n_threads,nkn) :: allNodesAssign
!    integer totalNodes, i, j, node, nkn, n_threads
    integer i, j, node, nkn, n_threads
    integer, dimension(n_threads) :: offset


    offset(1) = 0
    do i=1, n_threads
      do j=1,nkn
        allNodesAssign(i,j) = -1
      end do
      if(i .gt. 1) then
        offset(i) = offset(i-1) + numberNodes(i-1)
      end if
      do j=1, numberNodes(i)
        node = fullNodesAssign(offset(i)+j)
        allNodesAssign(i,node) = 1 
      end do
    end do
   
    return

  end subroutine


!######################################################################!
!*******************  start makeMyStruct subroutine  ******************!
!######################################################################!
!######################################################################!
!***  This subroutine build a integer substructure for the mapping  ***!
!***  between main structure and local structure to a process  ********!
!######################################################################!

  subroutine makeMyStruct(infoStruct, struct, mystruct, mysize)

!   infoStruct is a data structure with information on the comunication
!   that allows to make a mapping between a element/node of the struct
!   with a element of mystruct
!   struct is the structure needs to convert in substructure
!   mysize is the dimension of the structure for each element/node
!   my struct is the created substructure of the process

    implicit none
    ! input arguments
    integer mysize, mode
!    integer struct(mysize,1)
    type (COMMUNICATION_INFO) :: infoStruct
    integer, dimension(mysize,1) :: struct
    integer, dimension(mysize,1) :: mystruct
    integer :: st(MPI_STATUS_SIZE), ierr

    ! output arguments
!    type (MAKE_STRUCT) :: mystruct,struct

    ! local 
    integer n, j

    do n=1, infoStruct%numberID
      do j=1, mysize
        mystruct(j,n) = struct(j,infoStruct%globalID(n))
      end do
    end do


    return

  end subroutine


  subroutine transfer_domain

    !use basin 

    implicit none

    !write(6,*)'transfer_domain',myele%numberID,mynodes%numberID,neldi,nkndi,my_id

    integer, allocatable :: nen3v_new(:,:)
    integer, allocatable :: ipev_new(:)
    integer, allocatable :: ipv_new(:)
    integer, allocatable :: iarv_new(:)
    integer, allocatable :: iarnv_new(:)

    real, allocatable :: xgv_new(:)
    real, allocatable :: ygv_new(:)
    real, allocatable :: hm3v_new(:,:)

    integer own_nodes(nkndi)

    integer i,ii,ie,k,ierr

!   ----------------------------------
!   allocate aux arrays
!   ----------------------------------


    allocate(nen3v_new(3,nel))
    allocate(ipev_new(nel))
    allocate(ipv_new(nkn))
    allocate(iarv_new(nel))
    allocate(iarnv_new(nkn))

    allocate(xgv_new(nkn))
    allocate(ygv_new(nkn))
    allocate(hm3v_new(3,nel))


    own_nodes = 0
    do i=1,nkn
      k = mynodes%globalID(i)
      own_nodes(k) = i
      ipv_new(i) = ipv(k)
      iarnv_new(i) = iarnv(k)
      xgv_new(i) = xgv(k)
      ygv_new(i) = ygv(k)
    end do

    do i=1,nel
      ie = myele%globalID(i)
      do ii=1,3
        k = nen3v(ii,ie)
        if( own_nodes(k) <= 0 ) then
          write(6,*)'error stop transfer_domain: internal error'
          write(6,*) ie,k,own_nodes(k),my_id
          stop
        end if
        nen3v_new(ii,i) = own_nodes(k)
      end do
      hm3v_new(:,i) = hm3v(:,ie)
      ipev_new(i) = ipev(ie)
      iarv_new(i) = iarv(ie)
    end do

    call basin_init(nkn,nel)

    if(nel .ne. 0) then
    nen3v = nen3v_new
    ipev = ipev_new
    ipv = ipv_new
    iarv = iarv_new
    iarnv = iarnv_new
    xgv = xgv_new
    ygv = ygv_new
    hm3v = hm3v_new
    end if

!   ----------------------------------
!   deallocate temp arrays
!   ----------------------------------

    deallocate(nen3v_new)
    deallocate(ipev_new)
    deallocate(ipv_new)
    deallocate(iarv_new)
    deallocate(iarnv_new)

    deallocate(xgv_new)
    deallocate(ygv_new)
    deallocate(hm3v_new)

!   ----------------------------------
!   end routine
!   ----------------------------------

    return    

  end subroutine transfer_domain


!##########################################################################!
!**************  start spread_nodes_to_neighbor subroutine  ***************!
!##########################################################################!

   subroutine spread_nodes_to_neighbor(srequests, rrequests, data_receive, maxNodes)

        use MPI_Utility
        use commonStructures

        implicit none

        ! arguments
        integer, dimension(:,:),allocatable :: data_receive
        integer, dimension(:) :: srequests,rrequests
        integer maxNodes

        ! MPI
        integer ierr

        ! local
        integer i, comm_loop, tag, newtag
        integer sendNodes, receiveNodes

        comm_loop = mypart%mysend%sends
        allocate(data_receive(maxNodes,comm_loop))

        newtag = next_mpi_tag()

        do i = 1, comm_loop

           sendNodes = mynodes%numberID
           !write(6,*),numberNodes,maxNodes,my_id
           !write(6,*),mypart%mysend%process(i),my_id
           receiveNodes = numberNodes(mypart%mysend%process(i)+1)

           srequests(i) = MPI_REQUEST_NULL
           rrequests(i) = MPI_REQUEST_NULL

           tag = newtag + mypart%mysend%process(i) + my_id

           ! Non-blocking sends
           call MPI_Isend(mynodes%globalID(1),sendNodes,MPI_INTEGER, &
                mypart%mysend%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
           if(ierr .ne. MPI_SUCCESS) then
              write(6,*)'ierror:',ierr
           end if

           tag = newtag + my_id + mypart%mysend%process(i)
           ! Non-blocking receives
           call MPI_Irecv(data_receive(1,i),receiveNodes,MPI_INTEGER, &
                mypart%myreceive%process(i),tag, MPI_COMM_WORLD,rrequests(i),ierr)
           if(ierr .ne. MPI_SUCCESS) then
              write(6,*)'ierror:',ierr
           end if

        end do

        return

     end subroutine spread_nodes_to_neighbor

!##########################################################################!
!***************************  start waitAny  ******************************!
!##########################################################################!
!##########################################################################!
!* This subroutine is necessary when there are some pending non blocking  *!
!* communications, the subroutine take in input the requests vector       *!
!##########################################################################!

     subroutine waitAny(requests)

        implicit none

        ! argument
        integer, dimension(:) :: requests

        !local
        integer counter,ierr,i
        integer status(MPI_STATUS_SIZE),rindex

!        write(6,*)'size_requests',size(requests),my_id

        do i=1,size(requests)

           call MPI_WaitAny(size(requests), requests(1), rindex, status,ierr)
           if(ierr .ne. MPI_SUCCESS) then
              call abort
           end if

!           write(6,*)'count myrank',status(1)/4,my_id
!           write(6,*)'cancelled myrank',status(2),my_id
!           write(6,*)'sender myrank',status(3),my_id
!           write(6,*)'tag myrank',status(4),my_id
!           write(6,*)'error myrank',status(5),my_id
!           write(6,*)'index',rindex,my_id

        end do

        return
        
        end subroutine


!##########################################################################!
!*******************  start remove_dups0 subroutine  **********************!
!##########################################################################!
!##########################################################################!
!* This subroutine sort the elements of a vector and delete the duplicate *!
!* the elements must be integer positive included the 0                   *!
!##########################################################################!

  subroutine remove_dups0(array,length,res)

    implicit none
    
    integer, dimension(:) :: array         ! The input
    integer :: res(size(array))  ! The output
    integer :: length            ! The number of unique elements in output
    integer :: i, j, k, h, temp

    if( size(array) .le. 0) then
      length = 0
      return
    end if

    do i=1,size(array)
      res(i) = -1
    end do

    k = 1
    res(1) = array(1)
    outer: do i=2,size(array)
       if(array(i) .lt. 0) cycle outer
       do j=1,k
          if ((res(j) == array(i)) .or. (array(i) .lt. 0)) then
             ! Found a match so start looking again
             cycle outer
          end if
       end do
       ! No match found so add it to the output
       k = k + 1
       res(k) = array(i)

       h=k
       do while (res(h) .lt. res(h-1))
         temp = res(h-1)
         res(h-1) = res(h)
         res(h) = temp
         if(h .eq. 1) exit
         h = h-1
       end do

    end do outer

    length = k

    return

  end subroutine remove_dups0


!##########################################################################!
!*******************  start remove_dups subroutine  ***********************!
!##########################################################################!
!##########################################################################!
!* This subroutine sort the elements of a vector and delete the duplicate *!
!* the elements must be integer positive except the 0                     *!
!##########################################################################!

  subroutine remove_dups(array,length,res)

    implicit none
    
    integer, dimension(:) :: array         ! The input
    integer :: res(size(array))  ! The output
    integer :: length            ! The number of unique elements in output
    integer :: i, j, k, h, temp

    if( size(array) .le. 0) then
      length = 0
      return
    end if

    do i=1,size(res)
    !do i=1,size(array)
      res(i) = -1
    end do

    k = 1
    res(1) = array(1)
    outer: do i=2,size(array)
       if(array(i) .eq. 0) cycle outer
       do j=1,k
          if ((res(j) == array(i)) .or. (array(i) .eq. 0)) then
             ! Found a match so start looking again
             cycle outer
          end if
       end do
       ! No match found so add it to the output
       k = k + 1
       res(k) = array(i)

       h=k
       do while (res(h) .lt. res(h-1))
         temp = res(h-1)
         res(h-1) = res(h)
         res(h) = temp
         if(h .eq. 1) exit
         h = h-1
       end do

    end do outer

    length = k

    return

  end subroutine remove_dups


!######################################################################!
!******************  end Module Global_Graph_Data_Ele  ****************!
!######################################################################!

  end module Global_Graph_Data_Ele
