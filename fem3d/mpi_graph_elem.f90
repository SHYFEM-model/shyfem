
!#######################################################################!
!*******************  start module Global_Graph_Data  ******************!
!#######################################################################!

module Graph_Ele

  use mpi
  use communicationStruct
  use basin
  use shympi

contains


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

  subroutine makeMyEle(numMyVertices, numGlobalVertices, struct)

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
    integer, dimension(:), allocatable :: sortRank, rank,receives,sends   
    integer, dimension(:), allocatable :: part_receives, part_sends   
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
    allocate(receives(sizev))
    allocate(sends(sizev))
    allocate(part_receives(sizev))
    allocate(part_sends(sizev))
    allocate(rank(sizev))
    allocate(sortRank(sizev))

    do i=1,sizev
       receives(i) = 0
       sends(i) = 0
       part_receives(i) = 0
       part_sends(i) = 0
       rank(i) = -1
       sortRank(i) = 0
    end do

    do i=1, size(neighbor)
      neighbor(i) = 0
    end do

    myneighbor = 0
    k = 0
!    locID=1
    myele%numberID=0
    do i=1,numGlobalVertices
       if(allPartAssign(i) .eq. my_id) then


!          k=1
!          do j=0,nneighbor-1
!             if(allPartAssign(neighbID) .ne. my_id) then
!                myele%rankID(k,locID) = allPartAssign(neighbID)
!                myele%neighborID(k,locID) = neighbID
!                k= k+1
!             end if
!          end do


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
!          locID = locID + 1
       end if
    end do

    if(myele%numberID .ne. numMyVertices)then
      write(*,*)'Error stop makemyele'
      stop
    end if

    allocate(numberElements(n_threads))

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

    allocate(mypart%mysend%process(nprocesses),STAT=error)
    allocate(mypart%myreceive%process(nprocesses),STAT=error)
    allocate(mypart%mysend%part_send(nprocesses),STAT=error)
    allocate(mypart%myreceive%part_receive(nprocesses),STAT=error)
    allocate(mypart%mysend%all_send(nprocesses),STAT=error)
    allocate(mypart%myreceive%all_receive(nprocesses),STAT=error)

    mypart%mysend%sends = nprocesses
    mypart%myreceive%receives = nprocesses

    do i=1,nprocesses
       mypart%mysend%process(i) = sortRank(i)
       mypart%myreceive%process(i) = sortRank(i)
    end do

    deallocate(sortRank)

    myele%totalID = numMyVertices + length

    if(myele%totalID .ne. 0) then
      allocate(myele%globalID(myele%totalID),STAT=error)
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

    do i=numMyVertices+1,myele%totalID
      myele%globalID(i) = sort_neighbor(i-numMyVertices)
    end do

    maxSends=0
    maxReceives=0
    do h=1,nprocesses

       t = 0
       k = 0
       s = 0     
       p = 0
       receives = 0
       sends = 0
       part_receives = 0
       part_sends = 0

 
       do i=numMyVertices+1,myele%totalID

          if(allPartAssign(myele%globalID(i)) .eq. mypart%mysend%process(h)) then
             t = t +1
             receives(t) = i
          end if

       end do

       do ele=1,myele%numberID
          do n=1,3

             globalID=ieltv(n,myele%globalID(ele))
             if(globalID .gt. 0) then
             if(allPartAssign(globalID) .eq. mypart%mysend%process(h)) then

                do j=1,t
                   if(myele%globalID(receives(j)) .eq. globalID) then
                      k = k + 1
                      part_receives(k) = receives(j)
                   end if
                end do

                s = s + 1
                part_sends(s) = ele

             end if
             end if

          end do
       end do

       do ele=1,myele%numberID

          do j=1,3
                   
             node=struct(j,myele%globalID(ele))
             eleneighbor=ilinkv(node+1)-ilinkv(node)
                 
             do x=1,eleneighbor

                 if(lenkv(ilinkv(node)+x) .ne. 0) then
                    if(allPartAssign(lenkv(ilinkv(node)+x)) .eq. mypart%mysend%process(h)) then
                       p = p + 1
                       sends(p) = ele
                    end if
                 end if
             end do

          end do

       end do


       allocate(mypart%myreceive%all_receive(h)%items(t))

       mypart%myreceive%all_receive(h)%numItems = t

       do n=1,t
          mypart%myreceive%all_receive(h)%items(n) = receives(n)
       end do

       if(k .gt. 0) then
          call remove_dups(part_receives,length,sort_neighbor)
       else
         length = 0
       end if

       allocate(mypart%myreceive%part_receive(h)%items(length))

       mypart%myreceive%part_receive(h)%numItems = length

       do n=1,length
          mypart%myreceive%part_receive(h)%items(n) = sort_neighbor(n)
       end do

       if(s .gt. 0) then
          call remove_dups(part_sends,length,sort_neighbor)
       else
         length = 0
       end if

       allocate(mypart%mysend%part_send(h)%items(length))

       mypart%mysend%part_send(h)%numItems = length

       do n=1,length
          mypart%mysend%part_send(h)%items(n) = sort_neighbor(n)
       end do

       if(p .gt. 0) then
          call remove_dups(sends,length,sort_neighbor)
       else
         length = 0
       end if

       allocate(mypart%mysend%all_send(h)%items(length))

       mypart%mysend%all_send(h)%numItems = length

       do n=1,length
          mypart%mysend%all_send(h)%items(n) = sort_neighbor(n)
       end do

       if(length .gt. maxSends) then
          maxSends = length
       end if

       if(t .gt. maxReceives) then
          maxReceives = t
       end if

    end do

    mypart%mysend%maxItems = maxSends
    mypart%myreceive%maxItems = maxReceives

    deallocate(sort_neighbor)
    deallocate(receives)
    deallocate(sends)
    deallocate(part_receives)
    deallocate(part_sends)

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

  subroutine makeMyNodes(myele, mystruct, mynodes, struct, nkn)

    use commonStructures

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
    logical debug

    !
    integer ounit,error
    character*(20) filename
    character*(20) format_string
    integer, dimension(:), allocatable :: nnodes
    integer, dimension(:), allocatable :: nnodesAssign
    integer sendbuffer2
    integer, dimension(n_threads) :: recvbuffer2, displs2

    debug = .false.

    allocate(temp%globalID(3*myele%numberID))
    allocate(sort_nodes(3*myele%numberID))

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

    allocate(numberNodes(n_threads))

    call MPI_ALLGATHER(mynodes%numberID, 1, MPI_INTEGER, numberNodes, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    totalNodes = numberNodes(1)
    do i=2,n_threads
      totalNodes = totalNodes + numberNodes(i)
    end do

    if(my_id .eq. 0)then
      write(6,*)'totalNodes=',totalNodes
    end if

    allocate(allNodesAssign(n_threads,nkn))
    allocate(fullNodesAssign(totalNodes))


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

!    allocate(mynodes%rankID(mysize,mynodes%numberID))
!    allocate(mynodes%neighborID(mysize,mynodes%numberID))

!    mynodes%neighborID(:,:)=-1
!    mynodes%rankID(:,:)=-1

!    do i=1, length
!      n = ilinkv(sort_nodes(i)+1)-ilinkv(sort_nodes(i))
!      s=1
!      do k=1,n
!        node = linkv(ilinkv(sort_nodes(i))+ k)
!        if (AllNodesAssign(my_id+1,node) .ne. 1) then
!           do h=1,n_threads
!             if(AllNodesAssign(h,node) .eq. 1) then
!               mynodes%rankID(s,i) = h
!               mynodes%neighborID(s,i) = node
!               s = s +1
!               exit
!             end if
!           end do
!         end if
!       end do
!    end do

    allocate(neighbor(3*(myele%totalID-myele%numberID)))
    allocate(sort_neighbor(3*(myele%totalID-myele%numberID)))

    do i=1,3*(myele%totalID-myele%numberID)
      neighbor(i) = 0
    end do

    nneighbor = 0
    do i=myele%numberID+1,myele%totalID
      do j=1,3
        node = struct(j,myele%globalID(i))
        if(allNodesAssign(my_id+1,node) .ne. 1) then
         nneighbor = nneighbor + 1 
         neighbor(nneighbor) = node 
        end if
      end do
    end do

    if(nneighbor .gt. 0)then
      call remove_dups(neighbor,length,sort_neighbor)
    else
      length = 0
    end if

    mynodes%totalID = mynodes%numberID + length

    allocate(mynodes%globalID(mynodes%totalID))
 
    do i=1,mynodes%numberID
      mynodes%globalID(i) = sort_nodes(i)
    end do

    do i=mynodes%numberID+1,mynodes%totalID
      mynodes%globalID(i) = sort_neighbor(i-mynodes%numberID)
    end do


    allocate(temp_discard_node(3*myele%numberID))


    do i=1,3*myele%numberID
       temp_discard_node(i) = 0
    end do

       allocate(mypart%mysend%node_send(mypart%mysend%sends))
       allocate(mypart%mysend%node_temp(mypart%mysend%sends))
       do h=1,mypart%mysend%sends
          allocate(mypart%mysend%node_temp(h)%items(mypart%myreceive%all_receive(h)%numItems*3))
          do i=1,mypart%myreceive%all_receive(h)%numItems*3
             mypart%mysend%node_temp(h)%items(i)=0
          end do
          mypart%mysend%node_temp(h)%numItems=0
       end do

       if(myele%numberID .gt. 0) then

          counter = 0
          do i=myele%numberID, myele%totalID
              do j=1,3
                 k=struct(j,myele%globalID(i))
!loop1:           do n=1,mynodes%numberID                
                  do n=1,mynodes%numberID                
                    if(k .eq. mynodes%globalID(n)) then 
                       do h=1,mypart%myreceive%receives
                          do s=1,mypart%myreceive%all_receive(h)%numItems
                             if(mypart%myreceive%all_receive(h)%items(s) .eq. i) then

                                mypart%mysend%node_temp(h)%numItems=mypart%mysend%node_temp(h)%numItems+1

                                mypart%mysend%node_temp(h)%items(mypart%mysend%node_temp(h)%numItems) = n

                                if(my_id .gt. mypart%myreceive%process(h)) then
                                   counter = counter + 1
                                   temp_discard_node(counter) = k
!                                   exit
!                                   exit loop1
                                end if
                             end if
                          end do
                       end do
                    end if
                 end do
!                 end do loop1
              end do                             
          end do

      end if

    if(myele%numberID .gt. 0) then
    do h=1,mypart%mysend%sends
       if(mypart%mysend%node_temp(h)%numItems .gt. 0)then
          allocate(mypart%mysend%node_send(h)%items(mypart%myreceive%all_receive(h)%numItems*3),STAT=error)
          call remove_dups(mypart%mysend%node_temp(h)%items,length,mypart%mysend%node_send(h)%items)
       else
          length = 0
       end if
          mypart%mysend%node_send(h)%numItems=length
          deallocate(mypart%mysend%node_temp(h)%items)
    end do
    end if

    allocate(discard_nodes(size(temp_discard_node)))

    if(counter .gt. 0)then
      call remove_dups(temp_discard_node,length,discard_nodes)
    else
      length = 0
    end if

    allocate(univocal_nodes%localID(mynodes%numberID))
    allocate(univocal_nodes%globalID(mynodes%numberID))

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

!    if(my_id .gt. 9) then
!       format_string = "(A9,I2)"
!       write(filename,format_string)'globalID_',my_id
!    else
!       format_string = "(A9,I1)"
!       write(filename,format_string)'globalID_',my_id
!    end if
   
!    ounit=10+my_id 
!    open(unit=ounit, file=filename, action='write')
!    do i=1,univocal_nodes%numberID
!       write(ounit,*),univocal_nodes%globalID(i)
!    end do
!    close(ounit)

!    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    debug = .false.

!    call MPI_Barrier(MPI_COMM_WORLD,ierr)

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
       !write(ounit,*),discard_nodes
       !write(ounit,*),size(discard_nodes)
       !write(ounit,*),my_id

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
                nnodesAssign, recvbuffer2, displs2, MPI_INTEGER, MPI_COMM_WORLD, ierr)

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

    write(6,*)'myNodes is= ',mynodes%numberID,my_id

    if(my_id .ne. 0) then
      deallocate(allNodesAssign)
    end if

    deallocate(discard_nodes)
    deallocate (sort_nodes)
    deallocate (neighbor)
    deallocate (sort_neighbor)
    deallocate(temp%globalID)
    deallocate(temp_discard_node)

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


  subroutine makeMyNen3v(myele,mynodes, mystruct, struct)

    implicit none
  
    integer i,j,h,k
    integer mystruct(3,1), struct(3,1)
    type (COMMUNICATION_INFO) :: myele,mynodes,temp
    
    do i=1,myele%totalID
      do j=1,3
        k = mystruct(j,i)
        do h=1,mynodes%totalID
          if(mynodes%globalID(h) .eq. k ) then
            struct(j,i) = h
            exit
          end if
        end do
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


!######################################################################!
!*****************  start makeMyReal3DStruct subroutine  **************!
!######################################################################!
!######################################################################!
!***  This subroutine build a real substructure for the mapping  ***!
!***  between main structure and local structure to a process  ********!
!######################################################################!

  subroutine makeMyReal3DStruct(infoStruct, struct, mystruct, mysize, mysize2)

!   infoStruct is a data structure with information on the comunication
!   that allows to make a mapping between a element/node of the struct
!   with a element of mystruct
!   struct is the structure needs to convert in substructure
!   mysize is the dimension of the structure for each element/node
!   my struct is the created substructure of the process

    implicit none
    ! input arguments
    integer mysize, mysize2, mode
    type (COMMUNICATION_INFO) :: infoStruct
    real, dimension(mysize, mysize2, 1) :: struct
    real, dimension(mysize, mysize2, 1) :: mystruct

    ! local 
    integer n, j, i

    do n=1, infoStruct%totalID
      do j=1, mysize
        do i=1, mysize2
        mystruct(i,j,n) = struct(i,j,infoStruct%globalID(n))
        end do
      end do
    end do

    return

  end subroutine


!######################################################################!
!*****************  start makeMyRealStruct subroutine  ****************!
!######################################################################!
!######################################################################!
!***  This subroutine build a real substructure for the mapping  ***!
!***  between main structure and local structure to a process  ********!
!######################################################################!

  subroutine makeMyRealStruct(infoStruct, struct, mystruct, mysize)

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
    real, dimension(mysize,1) :: struct
    real, dimension(mysize,1) :: mystruct

    ! output arguments
!    type (MAKE_STRUCT) :: mystruct,struct

    ! local 
    integer n, j

    do n=1, infoStruct%totalID
      do j=1, mysize
        mystruct(j,n) = struct(j,infoStruct%globalID(n))
      end do
    end do

    return

  end subroutine


!######################################################################!
!*****************  start makeMyReal0Struct subroutine  ***************!
!######################################################################!
!######################################################################!
!***  This subroutine build a real substructure for the mapping  ***!
!***  between main structure and local structure to a process  ********!
!######################################################################!

  subroutine makeMyReal0Struct(infoStruct, struct, mystruct, mysize)

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
!    real, allocatable, dimension(:,:) :: struct
    real, dimension(0:mysize,1) :: struct
    real, dimension(0:mysize,1) :: mystruct

    ! output arguments
!    type (MAKE_STRUCT) :: mystruct,struct

    ! local 
    integer n, j

    do n=1, infoStruct%totalID
      do j=0, mysize
        mystruct(j,n) = struct(j,infoStruct%globalID(n))
      end do
    end do


    return

  end subroutine


!######################################################################!
!*****************  start makeMyDPStruct subroutine  ****************!
!######################################################################!
!######################################################################!
! This subroutine build a double precision substructure for the mapping!
!***  between main structure and local structure to a process  ********!
!######################################################################!

  subroutine makeMyDPStruct(infoStruct, struct, mystruct, mysize)

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
    double precision, dimension(mysize,1) :: struct
    double precision, dimension(mysize,1) :: mystruct

    ! output arguments
!    type (MAKE_STRUCT) :: mystruct,struct

    ! local 
    integer n, j

    do n=1, infoStruct%totalID
      do j=1, mysize
        mystruct(j,n) = struct(j,infoStruct%globalID(n))
      end do
    end do

    return

  end subroutine


!######################################################################!
!******************  start makeMySubStruct subroutine  ****************!
!######################################################################!
!######################################################################!
!***  This subroutine build a integer substructure for the mapping  ***!
!***  between main structure and local structure to a process  ********!
!######################################################################!

  subroutine makeMySubStruct(infoStruct, struct, mystruct)

!   infoStruct is a data structure with information on the comunication
!   that allows to make a mapping between a element/node of the struct
!   with a element of mystruct
!   struct is the structure needs to convert in substructure
!   mysize is the dimension of the structure for each element/node
!   my struct is the created substructure of the process

    implicit none
    ! input arguments
    integer mode
    integer, dimension(1) :: struct
    type (COMMUNICATION_INFO) :: infoStruct

    ! output arguments
    integer, dimension(1) :: mystruct
    ! local 
    integer n, j


    do n=1, infoStruct%totalID
        mystruct(n) = struct(infoStruct%globalID(n))
    end do

    return

  end subroutine


!######################################################################!
!****************  start makeMyRealSubStruct subroutine  **************!
!######################################################################!
!######################################################################!
!***  This subroutine build a integer substructure for the mapping  ***!
!***  between main structure and local structure to a process  ********!
!######################################################################!

  subroutine makeMyRealSubStruct(infoStruct, struct, mystruct)

!   infoStruct is a data structure with information on the comunication
!   that allows to make a mapping between a element/node of the struct
!   with a element of mystruct
!   struct is the structure needs to convert in substructure
!   mysize is the dimension of the structure for each element/node
!   my struct is the created substructure of the process
!   is node-based

    implicit none
    ! input arguments
    integer mode
    real, dimension(1) :: struct
    type (COMMUNICATION_INFO) :: infoStruct

    ! output arguments
    real, dimension(1) :: mystruct

    ! local 
    integer n, j

    do n=1, infoStruct%totalID
        mystruct(n) = struct(infoStruct%globalID(n))
    end do

    return

  end subroutine


!######################################################################!
!****************  start makeMyDPSubStruct subroutine  **************!
!######################################################################!
!######################################################################!
!** This subroutine build a D precision substructure for the mapping **!
!***  between main structure and local structure to a process  ********!
!######################################################################!

  subroutine makeMyDPSubStruct(infoStruct, struct, mystruct)

!   infoStruct is a data structure with information on the comunication
!   that allows to make a mapping between a element/node of the struct
!   with a element of mystruct
!   struct is the structure needs to convert in substructure
!   mysize is the dimension of the structure for each element/node
!   my struct is the created substructure of the process
!   is node-based

    implicit none
    ! input arguments
    integer mode
    double precision, dimension(1) :: struct
    type (COMMUNICATION_INFO) :: infoStruct

    ! output arguments
    double precision, dimension(1) :: mystruct

    ! local 
    integer n, j

    do n=1, infoStruct%totalID
        mystruct(n) = struct(infoStruct%globalID(n))
    end do

    return

  end subroutine



!##########################################################################!
!**********************  start setStruct subroutine  **********************!
!##########################################################################!
!##########################################################################!
!* This subroutine sets the two-dimensional substructures of int type for *!
!* the subprocesses each subprocess receve the full structures and through ! 
!* the use of the makeMyStruct subroutine creates its own substructure ****!
!##########################################################################!

  subroutine setStruct(fullStruct, infoStruct, mystruct, mysize, mode)

    !arguments
    !struct is the struct needs to be converted
    !mystruct is the substucture converted
    ! mode indicates the dimensions of the structure:
    ! 1 => struct(mysize,neldim)
    ! 2 => struct(mysize,nkndim)

    use basin

    implicit none

    integer mysize, mode, tag, length,i
    integer, allocatable, dimension(:,:) :: struct
    integer, dimension(:,:) :: mystruct
    integer fullStruct(mysize,1)
    integer :: st(MPI_STATUS_SIZE), ierr
    type (COMMUNICATION_INFO) :: infoStruct

    if(mode .eq. 1) then
       length = neldi * mysize
       allocate(struct(mysize, neldi))
    else if (mode .eq. 2) then
       allocate(struct(mysize, nkndi))
       length = nkndi * mysize
    else
      write(6,*)'Error in SetSubStruct, mode =', mode
    end if

    if (my_id .eq. 0)then
      struct=fullstruct
    end if

    call MPI_BCAST(struct, length, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)

    call makeMyStruct(infoStruct, struct, mystruct, mysize)

    deallocate(struct)

  end subroutine


!##########################################################################!
!*********************  start setSubStruct subroutine  ********************!
!##########################################################################!
!##########################################################################!
!* This subroutine sets the unidimensional substructures of int type for  *!
!* the subprocesses each subprocess receve the full structures and through ! 
!* the use of the makeMySubStruct subroutine creates its own substructure !
!##########################################################################!

  subroutine setSubStruct(fullStruct, infoStruct, mystruct, mode)

    !arguments
    !struct is the struct needs to be converted
    !mystruct is the substucture converted
    ! mode indicates the dimensions of the structure:
    ! 1 => struct(neldim)
    ! 2 => struct(nkndim)

    use basin

    implicit none

    integer, dimension(:) :: mystruct
    integer mode, tag, i, length,j
    integer, allocatable, dimension(:) :: struct
    integer fullStruct(1)
    integer :: st(MPI_STATUS_SIZE), ierr
    type (COMMUNICATION_INFO) :: infoStruct

    if(mode .eq. 1) then
       length = neldi
       allocate(struct(neldi))
    else if (mode .eq. 2) then
       allocate(struct(nkndi))
       length = nkndi
    else
      write(6,*)'Error in SetSubStruct, mode =', mode
    end if

    if (my_id .eq. 0)then
      struct=fullstruct
    end if

    call MPI_BCAST(struct, length, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)

    call makeMySubStruct(infoStruct, struct, mystruct)

    deallocate(struct)


  end subroutine


!##########################################################################!
!********************  start setReal3DStruct subroutine  ******************!
!##########################################################################!
!##########################################################################!
!This subroutine sets the three-dimensional substructures of real type for !
!* the subprocesses each subprocess receve the full structures and through ! 
!* the use of the makeMyRealStruct subroutine creates its own substructure !
!##########################################################################!

  subroutine setReal3DStruct(fullStruct, infoStruct, mystruct, mysize, mysize2, mode)

    !arguments
    !struct is the struct needs to be converted
    !mystruct is the substucture converted
    ! mode indicates the dimensions of the structure:
    ! 1 => struct(mysize,neldim)
    ! 2 => struct(mysize,nkndim)

    use basin

    implicit none

    integer mysize, mysize2, mode, tag, length, i, j, k
    real, allocatable, dimension(:,:,:) :: struct
    real, dimension(:,:,:) :: mystruct
    real, dimension(:,:,:) :: fullStruct
    integer :: st(MPI_STATUS_SIZE), ierr
    type (COMMUNICATION_INFO) :: infoStruct

    if(mode .eq. 1) then
       length = neldi * mysize * mysize2
       allocate(struct(mysize, mysize2, neldi))
    else if (mode .eq. 2) then
       allocate(struct(mysize, mysize2, nkndi))
       length = nkndi * mysize * mysize2
    else
      write(6,*)'Error in SetSubStruct, mode =', mode
    end if


    if (my_id .eq. 0)then
      struct=fullstruct
    end if

    call MPI_BCAST(struct, length, MPI_REAL, 0, MPI_COMM_WORLD,ierr)

    call makeMyReal3DStruct(infoStruct, struct, mystruct, mysize, mysize2)

    deallocate(struct)

  end subroutine



!##########################################################################!
!********************  start setRealStruct subroutine  ********************!
!##########################################################################!
!##########################################################################!
!* This subroutine sets the two-dimensional substructures of real type for !
!* the subprocesses each subprocess receve the full structures and through ! 
!* the use of the makeMyRealStruct subroutine creates its own substructure !
!##########################################################################!

  subroutine setRealStruct(fullStruct, infoStruct, mystruct, mysize, mode)

    !arguments
    !struct is the struct needs to be converted
    !mystruct is the substucture converted
    ! mode indicates the dimensions of the structure:
    ! 1 => struct(mysize,neldim)
    ! 2 => struct(mysize,nkndim)

    use basin

    implicit none

    integer mysize, mode, tag, length, i, j, k
    real, allocatable, dimension(:,:) :: struct
    real, dimension(:,:) :: mystruct
    real, dimension(:,:) :: fullStruct
    integer :: st(MPI_STATUS_SIZE), ierr
    type (COMMUNICATION_INFO) :: infoStruct

    if(mode .eq. 1) then
       length = neldi * mysize
       allocate(struct(mysize, neldi))
    else if (mode .eq. 2) then
       allocate(struct(mysize, nkndi))
       length = nkndi * mysize
    else
      write(6,*)'Error in SetSubStruct, mode =', mode
    end if


    if (my_id .eq. 0)then
      struct=fullstruct
    end if

    call MPI_BCAST(struct, length, MPI_REAL, 0, MPI_COMM_WORLD,ierr)

    call makeMyRealStruct(infoStruct, struct, mystruct, mysize)

    deallocate(struct)

  end subroutine


!#############################################################################!
!********************  start setRealSubStruct subroutine  ********************!
!#############################################################################!
!#############################################################################!
!* This subroutine sets the unidimensional substructures of real type for ****!
!* the subprocesses each subprocess receve the full structures and through ***! 
!* the use of the makeMyRealSubStruct subroutine creates its own substructure !
!#############################################################################!

  subroutine setRealSubStruct(fullStruct, infoStruct, mystruct, mode)

    !arguments
    !struct is the struct needs to be converted
    !mystruct is the substucture converted
    ! mode indicates the dimensions of the structure:
    ! 1 => struct(neldim)
    ! 2 => struct(nkndim)

    use basin

    implicit none


    real, dimension(:) :: mystruct
    integer mode,tag, i, length, j, k
    real, allocatable, dimension(:) :: struct
    real, dimension(:) :: fullStruct
    integer :: st(MPI_STATUS_SIZE), ierr
    type (COMMUNICATION_INFO) :: infoStruct

    if(mode .eq. 1) then
       length = neldi
       allocate(struct(neldi))
    else if (mode .eq. 2) then
       allocate(struct(nkndi))
       length = nkndi
    else
      write(6,*)'Error in SetSubStruct, mode =', mode
    end if

    if (my_id .eq. 0)then
      struct=fullstruct
    end if

    call MPI_BCAST(struct, length, MPI_REAL, 0, MPI_COMM_WORLD,ierr)

    call makeMyRealSubStruct(infoStruct, struct, mystruct)

    deallocate(struct)

    return

  end subroutine


!##########################################################################!
!********************  start setReal0Struct subroutine  *******************!
!##########################################################################!
!##########################################################################!
!* This subroutine sets the two-dimensional substructures of real type for !
!* the subprocesses each subprocess receve the full structures and through ! 
!* the use of the makeMyRealStruct subroutine creates its own substructure !
!##########################################################################!

  subroutine setReal0Struct(fullStruct, infoStruct, mystruct, mysize, mode)

    !arguments
    !struct is the struct needs to be converted
    !mystruct is the substucture converted
    ! mode indicates the dimensions of the structure:
    ! 1 => struct(mysize,neldim)
    ! 2 => struct(mysize,nkndim)

    use basin

    implicit none

    integer mysize, mode, tag, length
    real, allocatable, dimension(:,:) :: struct
    real, dimension(:,:) :: mystruct
    real, dimension(:,:) :: fullStruct
    integer :: st(MPI_STATUS_SIZE), ierr
    type (COMMUNICATION_INFO) :: infoStruct
    integer i, j, k


    if(mode .eq. 1) then
       length = neldi * (mysize +1)
       allocate(struct(0:mysize, neldi))
    else if (mode .eq. 2) then
       allocate(struct(0:mysize, nkndi))
       length = nkndi * (mysize +1)
    else
      write(6,*)'Error in SetSubStruct, mode =', mode
    end if


    if (my_id .eq. 0)then
      struct=fullstruct
    end if

    call MPI_BCAST(struct, length, MPI_REAL, 0, MPI_COMM_WORLD,ierr)

    call makeMyReal0Struct(infoStruct, struct, mystruct, mysize)

    deallocate(struct)

    return

  end subroutine


!##########################################################################!
!********************  start setDPStruct subroutine  **********************!
!##########################################################################!
!##########################################################################!
!* This subroutine sets the two-dimensional substructures of double       *!
!* precision type for the subprocesses. Each subprocess receve the full   *!
!* structures and through the use of the makeMyDPStruct subroutine creates !
!* its own substructure                                                   *!
!##########################################################################!

  subroutine setDPStruct(fullStruct, infoStruct, mystruct, mysize, mode)

    !arguments
    !struct is the struct needs to be converted
    !mystruct is the substucture converted
    ! mode indicates the dimensions of the structure:
    ! 1 => struct(mysize,neldim)
    ! 2 => struct(mysize,nkndim)

    use basin

    implicit none

    integer mysize, mode, tag, length, i, j, k
    double precision, allocatable, dimension(:,:) :: struct
    double precision, dimension(:,:) :: mystruct
    double precision, dimension(:,:) :: fullStruct
    integer :: st(MPI_STATUS_SIZE), ierr
    type (COMMUNICATION_INFO) :: infoStruct

    if(mode .eq. 1) then
       length = neldi * mysize
       allocate(struct(mysize, neldi))
    else if (mode .eq. 2) then
       allocate(struct(mysize, nkndi))
       length = nkndi * mysize
    else
      write(6,*)'Error in SetSubStruct, mode =', mode
    end if

    if (my_id .eq. 0)then
      struct=fullstruct
    end if

    call MPI_BCAST(struct, length, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,ierr)

    call makeMyDPStruct(infoStruct, struct, mystruct, mysize)

    deallocate(struct)


  end subroutine


!#############################################################################!
!********************  start setDPSubStruct subroutine  **********************!
!#############################################################################!
!#############################################################################!
!* This subroutine sets the two-dimensional substructures of double precision !
!* type for the subprocesses. Each subprocess receve the full structures and *!
!* through the use of the makeMyDPStruct subroutine creates its own          *!
!* substructure.                                                             *!
!#############################################################################!

  subroutine setDPSubStruct(fullStruct, infoStruct, mystruct, mode)

    !arguments
    !struct is the struct needs to be converted
    !mystruct is the substucture converted
    ! mode indicates the dimensions of the structure:
    ! 1 => struct(neldim)
    ! 2 => struct(nkndim)

    use basin

    implicit none


    double precision, dimension(:) :: mystruct
    integer mode,tag, i, length, j, k
    double precision, allocatable, dimension(:) :: struct
    double precision, dimension(:) :: fullStruct
    integer :: st(MPI_STATUS_SIZE), ierr
    type (COMMUNICATION_INFO) :: infoStruct

    if(mode .eq. 1) then
       length = neldi
       allocate(struct(neldi))
    else if (mode .eq. 2) then
       allocate(struct(nkndi))
       length = nkndi
    else
      write(6,*)'Error in SetSubStruct, mode =', mode
    end if

    if (my_id .eq. 0)then
      struct=fullstruct
    end if

    call MPI_BCAST(struct, length, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,ierr)

    call makeMyDPSubStruct(infoStruct, struct, mystruct)

    deallocate(struct)


  end subroutine


!#############################################################################!
!******************  start distributeVariables subroutine  *******************!
!#############################################################################!
!#############################################################################!
!** This subroutine update the required variables for each subprocess  *******!
!** necessary to run the conz3d routine ********************!
!#############################################################################!

     subroutine distributeVariables(    &
                       & dt,rkpar,difmol,itvd    &
                       & ,itvdv,robs,wsink,azpar    &
                       & ,adpar,aapar    &
                       & )



        implicit none

        include 'nlevel.h'

        !common variables
        !integer :: nlvdi,nlv
        !common /level/ nlvdi,nlv

        !arguments
        real difmol
        integer itvd
        integer itvdv
        real robs
        real wsink
        real dt,rkpar,azpar,adpar,aapar                        !$$azpar

        !local
        integer, dimension(6) :: tempintv
        real, dimension(8) :: temprealv

        !MPI variables
        integer :: st(MPI_STATUS_SIZE), ierr

        if (my_id .eq. 0) then
           tempintv(1) = nlvdi
           tempintv(2) = nlv
           tempintv(3) = itvd
           tempintv(4) = itvdv
           temprealv(1) = difmol
           temprealv(2) = robs
           temprealv(3) = wsink
           temprealv(4) = dt
           temprealv(5) = rkpar
           temprealv(6) = azpar
           temprealv(7) = adpar
           temprealv(8) = aapar
        end if



        call MPI_BCAST(temprealv, 8, MPI_REAL, 0, MPI_COMM_WORLD,ierr)
        call MPI_BCAST(tempintv, 4, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)

        if(my_id .ne. 0) then
           nlvdi = tempintv(1)        
           nlv = tempintv(2)        
           itvd = tempintv(3)        
           itvdv = tempintv(4)        
           difmol = temprealv(1)
           robs = temprealv(2)
           wsink = temprealv(3)
           dt = temprealv(4)
           rkpar = temprealv(5)
           azpar = temprealv(6)
           adpar = temprealv(7)
           aapar = temprealv(8)
        end if


        return

!       end subroutine


     end subroutine


!#############################################################################!
!***********************  start countLevel subroutine  ***********************!
!#############################################################################!
!#############################################################################!
!* This subroutine count the number of levels of each process and the total  *!
!* level for each process                                                    *!
!#############################################################################!

  subroutine countLevel(myele, ilhv, numlevels, totalLevels)

    use commonStructures

    implicit none

    include 'param.h'

    integer i
    integer, dimension(:) :: ilhv
    type (COMMUNICATION_INFO) :: myele
    integer :: st(MPI_STATUS_SIZE), ierr

    integer numlevels, totalLevels    

    numlevels = 0
    do i=1,myele%numberID
      numlevels = numlevels + ilhv(i)
    end do

    allocate(numberLevels(n_threads))
    call MPI_ALLGATHER(numlevels,1,MPI_INTEGER,numberLevels,1,MPI_INTEGER, MPI_COMM_WORLD,ierr)

    totalLevels = numberLevels(1)
    maxlevelsproc = numberLevels(1)
    do i=2,n_threads
      if (numberLevels(i) > maxlevelsproc) maxlevelsproc = numberLevels(i)
      totalLevels = totalLevels + numberLevels(i)
    end do

    return

  end subroutine


!#############################################################################!
!***********************  start recountLevel subroutine  ***********************!
!#############################################################################!
!#############################################################################!
!* This subroutine count the number of levels of each process                *!
!#############################################################################!

  subroutine recountLevel(myele, ilhv)

    use commonStructures

    implicit none

    include 'param.h'

    integer i
    integer, dimension(:) :: ilhv
    type (COMMUNICATION_INFO) :: myele
    integer :: st(MPI_STATUS_SIZE), ierr

    integer numlevels    

    numlevels = 0
    do i=1,myele%numberID
      numlevels = numlevels + ilhv(i)
    end do

    if ( numlevels .ge. (1.20 * numberLevels(my_id+1)) ) then
      write(6,*)'my old levels =', numberLevels(my_id+1), my_id
      write(6,*)'my new levels =', numlevels, my_id
    end if

    return

  end subroutine

!#############################################################################!
!**************************  start rebuildStructures  ************************!
!#############################################################################!
!#############################################################################!
!******** This subroutine update the variables itend and it to allow *********!
!***************** the execution of the loop on the time steps ***************!
!#############################################################################!

  subroutine rebuildIntSubStructures(struct, size2, sizeArray, newStruct, mode)

    use commonStructures

    implicit none

    include 'femtime.h'
    include 'nlevel.h'

    integer size2, mode
    integer :: count_tag, i
!    integer :: nlvdi,nlv
    integer :: st(MPI_STATUS_SIZE), ierr
!    integer itanf,itend,idt,nits,niter,it
    integer :: sendbuffer
    integer, dimension(size2) :: struct, newStruct
    integer, dimension(n_threads) :: displsProc, sizeArray, recvbuffer, displs
    integer, allocatable, dimension(:) :: fullStruct
!    common /femtim/ itanf,itend,idt,nits,niter,it
!    common /level/ nlvdi,nlv


    if(my_id .eq. 0) then
      if(mode .eq. 1) then  !structure based on elements
        allocate(fullStruct(size2)) 
      else if(mode .eq. 2) then  !structure based on nodes
        allocate(fullStruct(totalNodes))
      else
        write(6,*)'error in rebuildRealStructures', mode, my_id
        stop
      end if
    end if

    sendbuffer = sizeArray(my_id+1)

    displs(1) = 0
    do i=2,n_threads
      displs(i) = displs(i-1) + sizeArray(i-1)
    end do

    call MPI_GATHERV(struct, sendbuffer, MPI_INTEGER, fullStruct, sizeArray, displs, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    if(my_id .eq. 0) then
      displsProc(1) = 0
      do i=2,n_threads
        displsProc(i) = displsProc(i-1) + sizeArray(i-1)
      end do
      if(mode .eq. 1) then
        call rebuildIntSubStruct(fullStruct, displsProc, newStruct, size2)
      else
        call rebuildIntSubStruct2(fullStruct, displsProc, newStruct, size2)
      end if
      deallocate(fullStruct)
    end if

  end subroutine


  subroutine rebuildIntSubStruct(fullStruct, displsProc, newStruct, size2)

    use commonStructures

    implicit none

    integer :: place, i, process, size2
    integer :: st(MPI_STATUS_SIZE), ierr
    integer, dimension(:) :: countProc(n_threads), displsProc(1)
    integer, dimension(:) :: fullStruct, newStruct

    ! countProc(1) = counter for items of the master process
    do i=1,n_threads
      countProc(i)=1
    end do

    do i=1, size2
      process = allPartAssign(i)
      place = displsProc(process+1) + countProc(process+1)
      newStruct(i) = fullStruct(place)
      countProc(process+1) = countProc(process+1) + 1
    end do
   
    return 

  end subroutine


  subroutine rebuildIntSubStruct2(fullStruct, displsProc, newStruct, size2)

    use commonStructures

    implicit none

    integer :: place, i, j, k, size2
    integer :: st(MPI_STATUS_SIZE), ierr
    integer, dimension(:) :: countProc(n_threads), displsProc(1)
    integer, dimension(:) :: fullStruct, newStruct
    logical bexit 

    ! countProc(1) = counter for items of the master process
    do i=1,n_threads
      countProc(i)=1
    end do

    do i=1, size2
      bexit = .false.
      do j=1,n_threads
        if(allNodesAssign(j,i) .eq. 1) then
          place = displsProc(j) + countProc(j)
          countProc(j) = countProc(j) + 1
          if(bexit .eqv. .true.) then
            cycle
          end if
          newStruct(i) = fullStruct(place)
          bexit = .true.
        end if
      end do
    end do
   
    return 

  end subroutine


!#############################################################################!
!**************************  start rebuildStructures  ************************!
!#############################################################################!
!#############################################################################!
!******** This subroutine update the variables itend and it to allow *********!
!***************** the execution of the loop on the time steps ***************!
!#############################################################################!

  subroutine rebuildRealSubStructures(struct, size2, sizeArray, newStruct, mode)

    use commonStructures

    implicit none

    include 'femtime.h'
    include 'nlevel.h'

    integer size2, mode
    integer :: count_tag, i
!    integer :: nlvdi,nlv
    integer :: st(MPI_STATUS_SIZE), ierr
!    integer itanf,itend,idt,nits,niter,it
    integer :: sendbuffer
    real, dimension(size2) :: struct, newStruct
    integer, dimension(n_threads) :: displsProc, sizeArray, recvbuffer, displs
    real, allocatable, dimension(:) :: fullStruct
!    common /femtim/ itanf,itend,idt,nits,niter,it
!    common /level/ nlvdi,nlv


    if(my_id .eq. 0) then
      if(mode .eq. 1) then  !structure based on elements
        allocate(fullStruct(size2)) 
      else if(mode .eq. 2) then  !structure based on nodes
        allocate(fullStruct(totalNodes))
      else
        write(6,*)'error in rebuildRealStructures', mode, my_id
        stop
      end if
    end if

    sendbuffer = sizeArray(my_id+1)

    displs(1) = 0
    do i=2,n_threads
      displs(i) = displs(i-1) + sizeArray(i-1)
    end do

    call MPI_GATHERV(struct, sendbuffer, MPI_REAL, fullStruct, sizeArray, displs, MPI_REAL, 0, MPI_COMM_WORLD, ierr)

    if(my_id .eq. 0) then
      displsProc(1) = 0
      do i=2,n_threads
        displsProc(i) = displsProc(i-1) + sizeArray(i-1)
      end do
      if(mode .eq. 1) then
        call rebuildRealSubStruct(fullStruct, displsProc, newStruct, size2)
      else
        call rebuildRealSubStruct2(fullStruct, displsProc, newStruct, size2)
      end if
      deallocate(fullStruct)
    end if

  end subroutine


!#############################################################################!
!**************************  start rebuildStructures  ************************!
!#############################################################################!
!#############################################################################!
!******** This subroutine update the variables itend and it to allow *********!
!***************** the execution of the loop on the time steps ***************!
!#############################################################################!

  subroutine rebuildRealStructures(struct, size1, size2, sizeArray, newStruct, mode)

    use commonStructures

    implicit none

    include 'femtime.h'
    include 'nlevel.h'

    integer size1, size2, mode
    integer :: count_tag, i
!    integer :: nlvdi,nlv
    integer :: st(MPI_STATUS_SIZE), ierr
!    integer itanf,itend,idt,nits,niter,it
    integer :: sendbuffer
    real, dimension(size1,1) :: struct, newStruct
    integer, dimension(n_threads) :: displsProc, sizeArray, recvbuffer, displs
    real, allocatable, dimension(:,:) :: fullStruct
!    common /femtim/ itanf,itend,idt,nits,niter,it
!    common /level/ nlvdi,nlv


    if(my_id .eq. 0) then
      if(mode .eq. 1) then  !structure based on elements
        allocate(fullStruct(size1,size2)) 
      else if(mode .eq. 2) then  !structure based on nodes
        allocate(fullStruct(size1,totalNodes))
      else
        write(6,*)'error in rebuildRealStructures', mode, my_id
        stop
      end if
      do i=1,n_threads
        recvbuffer(i) = sizeArray(i) * size1
      end do
    end if

    sendbuffer = size1 * sizeArray(my_id+1)

    displs(1) = 0
    do i=2,n_threads
      displs(i) = displs(i-1) + size1 * sizeArray(i-1)
    end do

    call MPI_GATHERV(struct, sendbuffer, MPI_REAL, fullStruct, recvbuffer, displs, MPI_REAL, 0, MPI_COMM_WORLD, ierr)

    if(my_id .eq. 0) then
      displsProc(1) = 0
      do i=2,n_threads
        displsProc(i) = displsProc(i-1) + sizeArray(i-1)
      end do
      if(mode .eq. 1) then
        call rebuildRealStruct(fullStruct, displsProc, newStruct, size1, size2)
      else
        call rebuildRealStruct2(fullStruct, displsProc, newStruct, size1, size2)
      end if
      deallocate(fullStruct)
    end if

  end subroutine



  subroutine rebuildDP_ccle(struct, size1, size2, sizeArray, newStruct, mode)

    use commonStructures

    implicit none

    include 'femtime.h'
    include 'nlevel.h'

    integer size1, size2, mode
    integer :: count_tag, i
!    integer :: nlvdi,nlv
    integer :: st(MPI_STATUS_SIZE), ierr
!    integer itanf,itend,idt,nits,niter,it
    integer :: sendbuffer
    double precision, dimension(size1,3,nel) :: struct
    double precision, dimension(size1,3,neldi) :: newStruct
    integer, dimension(n_threads) :: displsProc, sizeArray, recvbuffer, displs
    double precision, allocatable, dimension(:,:,:) :: fullStruct
!    common /femtim/ itanf,itend,idt,nits,niter,it
!    common /level/ nlvdi,nlv


    if(my_id .eq. 0) then
      if(mode .eq. 1) then  !structure based on elements
        allocate(fullStruct(size1,3,size2)) 
      else 
        write(6,*)'error in rebuildDP_ccle', mode, my_id
        stop
      end if
      do i=1,n_threads
        recvbuffer(i) = sizeArray(i) * size1 * 3
      end do
    end if

    sendbuffer = size1 * 3 * sizeArray(my_id+1)

    displs(1) = 0
    do i=2,n_threads
      displs(i) = displs(i-1) + size1 * 3 * sizeArray(i-1)
    end do

    call MPI_GATHERV(struct, sendbuffer,MPI_DOUBLE_PRECISION,fullStruct, &
                recvbuffer, displs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    if(my_id .eq. 0) then
      displsProc(1) = 0
      do i=2,n_threads
        displsProc(i) = displsProc(i-1) + sizeArray(i-1)
      end do
      if(mode .eq. 1) then
        call rebuild_ccle(fullStruct, displsProc, newStruct, size1, size2)
      end if
      deallocate(fullStruct)
    end if

  end subroutine


  subroutine rebuild_ccle(fullStruct, displsProc, newStruct, size1, size2)

    use commonStructures

    implicit none

    integer :: place, i, process, j, size1, size2,k
    integer :: st(MPI_STATUS_SIZE), ierr
    integer, dimension(:) :: countProc(n_threads), displsProc(1)
    double precision, dimension(:,:,:) :: fullStruct, newStruct

    ! countProc(1) = counter for items of the master process
    do i=1,n_threads
      countProc(i)=1
    end do

!    open(unit=1500, file="testddxv2.txt", action='write')
    do i=1, size2
      process = allPartAssign(i)
      place = displsProc(process+1) + countProc(process+1)
!      if(i .gt. 100 .and. i .lt. 200) then
!        write(1500,*)'ie =',i
!      end if
      do j=1,3
        do k=1,size1
          newStruct(k,j,i) = fullStruct(k,j,place)
        end do
      end do
      countProc(process+1) = countProc(process+1) + 1
    end do
!    close(1500)
    

  end subroutine



!#############################################################################!
!**************************  start rebuildStructures  ************************!
!#############################################################################!
!#############################################################################!
!******** This subroutine update the variables itend and it to allow *********!
!***************** the execution of the loop on the time steps ***************!
!#############################################################################!

  subroutine rebuildReal0Structures(struct, size1, size2, sizeArray, newStruct, mode)

    use commonStructures

    implicit none

    include 'femtime.h'
    include 'nlevel.h'

    integer size1, size2, mode
    integer :: count_tag, i
!    integer :: nlvdi,nlv
    integer :: st(MPI_STATUS_SIZE), ierr
!    integer itanf,itend,idt,nits,niter,it
    integer :: sendbuffer
    real, dimension(0:size1,1) :: struct, newStruct
    integer, dimension(n_threads) :: displsProc, sizeArray, recvbuffer, displs
    real, allocatable, dimension(:,:) :: fullStruct
!    common /femtim/ itanf,itend,idt,nits,niter,it
!    common /level/ nlvdi,nlv


    if(my_id .eq. 0) then
      if(mode .eq. 1) then  !structure based on elements
        allocate(fullStruct(0:size1,size2)) 
      else if(mode .eq. 2) then  !structure based on nodes
        allocate(fullStruct(0:size1,totalNodes))
      else
        write(6,*)'error in rebuildRealStructures', mode, my_id
        stop
      end if
      do i=1,n_threads
        recvbuffer(i) = sizeArray(i) * (size1+1)
      end do
    end if

    sendbuffer = (size1+1) * sizeArray(my_id+1)

    displs(1) = 0
    do i=2,n_threads
      displs(i) = displs(i-1) + (size1+1) * sizeArray(i-1)
    end do

    call MPI_GATHERV(struct, sendbuffer, MPI_REAL, fullStruct, recvbuffer, displs, MPI_REAL, 0, MPI_COMM_WORLD, ierr)

    if(my_id .eq. 0) then
      displsProc(1) = 0
      do i=2,n_threads
        displsProc(i) = displsProc(i-1) + sizeArray(i-1)
      end do
      if(mode .eq. 1) then
        call rebuildReal0Struct(fullStruct, displsProc, newStruct, size1, size2)
      else
        call rebuildReal0Struct2(fullStruct, displsProc, newStruct, size1, size2)
      end if
      deallocate(fullStruct)
    end if

    return


  end subroutine


  subroutine rebuildDP0Structures(struct, size1, size2, sizeArray, newStruct, mode)

    use commonStructures

    implicit none

    include 'femtime.h'
    include 'nlevel.h'

    integer size1, size2, mode
    integer :: count_tag, i
!    integer :: nlvdi,nlv
    integer :: st(MPI_STATUS_SIZE), ierr
!    integer itanf,itend,idt,nits,niter,it
    integer :: sendbuffer
    double precision, dimension(0:size1,1) :: struct, newStruct
    integer, dimension(n_threads) :: displsProc, sizeArray, recvbuffer, displs
    double precision, allocatable, dimension(:,:) :: fullStruct
!    common /femtim/ itanf,itend,idt,nits,niter,it
!    common /level/ nlvdi,nlv


    if(my_id .eq. 0) then
      if(mode .eq. 1) then  !structure based on elements
        allocate(fullStruct(0:size1,size2)) 
      else if(mode .eq. 2) then  !structure based on nodes
        allocate(fullStruct(0:size1,totalNodes))
      else
        write(6,*)'error in rebuildRealStructures', mode, my_id
        stop
      end if
      do i=1,n_threads
        recvbuffer(i) = sizeArray(i) * (size1+1)
      end do
    end if

    sendbuffer = (size1+1) * sizeArray(my_id+1)

    displs(1) = 0
    do i=2,n_threads
      displs(i) = displs(i-1) + (size1+1) * sizeArray(i-1)
    end do

    call MPI_GATHERV(struct,sendbuffer,MPI_DOUBLE_PRECISION,fullStruct, &
                recvbuffer,displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

    if(my_id .eq. 0) then
      displsProc(1) = 0
      do i=2,n_threads
        displsProc(i) = displsProc(i-1) + sizeArray(i-1)
      end do
      if(mode .eq. 1) then
        call rebuildDP0Struct(fullStruct, displsProc, newStruct, size1, size2)
      else
        call rebuildDP0Struct2(fullStruct, displsProc, newStruct, size1, size2)
      end if
      deallocate(fullStruct)
    end if

    return


  end subroutine


!#############################################################################!
!*************************  start rebuildRealStruct2 *************************!
!#############################################################################!
!#############################################################################!
!******** This subroutine rebuild a complete structure after a parallel ******!
!**** subroutine, significant for the root process only. double precision ****!
!#############################################################################!

  subroutine rebuildReal0Struct2(fullStruct, displsProc, newStruct, size1, size2)

    use commonStructures

    implicit none

    integer :: place, i, j, k, size1, size2
    integer :: st(MPI_STATUS_SIZE), ierr
    integer, dimension(:) :: countProc(n_threads), displsProc(1)
    real, dimension(:,:) :: fullStruct, newStruct
    logical bexit 

    ! countProc(1) = counter for items of the master process
    do i=1,n_threads
      countProc(i)=1
    end do

    do i=1, size2
      bexit = .false.
      do j=1,n_threads
        if(allNodesAssign(j,i) .eq. 1) then
          place = displsProc(j) + countProc(j)
          countProc(j) = countProc(j) + 1
          if(bexit .eqv. .true.) then
            cycle
          end if
          !do k=0, size1
          !  newStruct(k,i) = fullStruct(k,place)
            newStruct(:,i) = fullStruct(:,place)
          !end do
          bexit = .true.
        end if
      end do
    end do
   
    return 

  end subroutine



!#############################################################################!
!*************************  start rebuildRealStruct2 *************************!
!#############################################################################!
!#############################################################################!
!******** This subroutine rebuild a complete structure after a parallel ******!
!**** subroutine, significant for the root process only. double precision ****!
!#############################################################################!

  subroutine rebuildDP0Struct2(fullStruct, displsProc, newStruct, size1, size2)

    use commonStructures

    implicit none

    integer :: place, i, j, k, size1, size2
    integer :: st(MPI_STATUS_SIZE), ierr
    integer, dimension(:) :: countProc(n_threads), displsProc(1)
    double precision, dimension(:,:) :: fullStruct, newStruct
    logical bexit 

    ! countProc(1) = counter for items of the master process
    do i=1,n_threads
      countProc(i)=1
    end do

    do i=1, size2
      bexit = .false.
      do j=1,n_threads
        if(allNodesAssign(j,i) .eq. 1) then
          place = displsProc(j) + countProc(j)
          countProc(j) = countProc(j) + 1
          if(bexit .eqv. .true.) then
            cycle
          end if
          !do k=0, size1
          !  newStruct(k,i) = fullStruct(k,place)
            newStruct(:,i) = fullStruct(:,place)
          !end do
          bexit = .true.
        end if
      end do
    end do
   
    return 

  end subroutine



!#############################################################################!
!*************************  start rebuildRealStruct2 *************************!
!#############################################################################!
!#############################################################################!
!******** This subroutine rebuild a complete structure after a parallel ******!
!**** subroutine, significant for the root process only. double precision ****!
!#############################################################################!

  subroutine rebuildRealStruct2(fullStruct, displsProc, newStruct, size1, size2)

    use commonStructures

    implicit none

    integer :: place, i, j, k, size1, size2
    integer :: st(MPI_STATUS_SIZE), ierr
    integer, dimension(:) :: countProc(n_threads), displsProc(1)
    real, dimension(:,:) :: fullStruct, newStruct
    logical bexit 

    ! countProc(1) = counter for items of the master process
    do i=1,n_threads
      countProc(i)=1
    end do

    do i=1, size2
      bexit = .false.
      do j=1,n_threads
        if(allNodesAssign(j,i) .eq. 1) then
          place = displsProc(j) + countProc(j)
          countProc(j) = countProc(j) + 1
          if(bexit .eqv. .true.) then
            cycle
          end if
          do k=1, size1
            newStruct(k,i) = fullStruct(k,place)
          end do
          bexit = .true.
        end if
      end do
    end do
   
    return 

  end subroutine



!#############################################################################!
!*************************  start rebuildRealStruct  *************************!
!#############################################################################!
!#############################################################################!
!******** This subroutine rebuild a complete structure after a parallel ******!
!**** subroutine, significant for the root process only. double precision ****!
!#############################################################################!

  subroutine rebuildReal0Struct(fullStruct, displsProc, newStruct, size1, size2)

    use commonStructures

    implicit none

    integer :: place, i, process, j, size1, size2
    integer :: st(MPI_STATUS_SIZE), ierr
    integer, dimension(:) :: countProc(n_threads), displsProc(1)
    real, dimension(:,:) :: fullStruct, newStruct

    ! countProc(1) = counter for items of the master process
    do i=1,n_threads
      countProc(i)=1
    end do

!    open(unit=1500, file="testddxv2.txt", action='write')
    do i=1, size2
      process = allPartAssign(i)
      place = displsProc(process+1) + countProc(process+1)
!      if(i .gt. 100 .and. i .lt. 200) then
!        write(1500,*)'ie =',i
!      end if
      do j=0,size1
        newStruct(j,i) = fullStruct(j,place)
!        if(i .gt. 100 .and. i .lt. 200) then
!          write(1500,*)'newstruct',newStruct(j,i)
!        end if
      end do
      countProc(process+1) = countProc(process+1) + 1
    end do
!    close(1500)
    

  end subroutine


!#############################################################################!
!*************************  start rebuildRealStruct  *************************!
!#############################################################################!
!#############################################################################!
!******** This subroutine rebuild a complete structure after a parallel ******!
!**** subroutine, significant for the root process only. double precision ****!
!#############################################################################!

  subroutine rebuildDP0Struct(fullStruct, displsProc, newStruct, size1, size2)

    use commonStructures

    implicit none

    integer :: place, i, process, j, size1, size2
    integer :: st(MPI_STATUS_SIZE), ierr
    integer, dimension(:) :: countProc(n_threads), displsProc(1)
    double precision, dimension(:,:) :: fullStruct, newStruct

    ! countProc(1) = counter for items of the master process
    do i=1,n_threads
      countProc(i)=1
    end do

!    open(unit=1500, file="testddxv2.txt", action='write')
    do i=1, size2
      process = allPartAssign(i)
      place = displsProc(process+1) + countProc(process+1)
!      if(i .gt. 100 .and. i .lt. 200) then
!        write(1500,*)'ie =',i
!      end if
      do j=0,size1
        newStruct(j,i) = fullStruct(j,place)
!        if(i .gt. 100 .and. i .lt. 200) then
!          write(1500,*)'newstruct',newStruct(j,i)
!        end if
      end do
      countProc(process+1) = countProc(process+1) + 1
    end do
!    close(1500)
    

  end subroutine



!#############################################################################!
!*************************  start rebuildRealStruct  *************************!
!#############################################################################!
!#############################################################################!
!******** This subroutine rebuild a complete structure after a parallel ******!
!**** subroutine, significant for the root process only. double precision ****!
!#############################################################################!

  subroutine rebuildRealStruct(fullStruct, displsProc, newStruct, size1, size2)

    use commonStructures

    implicit none

    integer :: place, i, process, j, size1, size2
    integer :: st(MPI_STATUS_SIZE), ierr
    integer, dimension(:) :: countProc(n_threads), displsProc(1)
    real, dimension(:,:) :: fullStruct, newStruct

    ! countProc(1) = counter for items of the master process
    do i=1,n_threads
      countProc(i)=1
    end do

!    open(unit=1500, file="testddxv2.txt", action='write')
    do i=1, size2
      process = allPartAssign(i)
      place = displsProc(process+1) + countProc(process+1)
!      if(i .gt. 100 .and. i .lt. 200) then
!        write(1500,*)'ie =',i
!      end if
      do j=1,size1
        newStruct(j,i) = fullStruct(j,place)
!        if(i .gt. 100 .and. i .lt. 200) then
!          write(1500,*)'newstruct',newStruct(j,i)
!        end if
      end do
      countProc(process+1) = countProc(process+1) + 1
    end do
!    close(1500)
    

  end subroutine


!#############################################################################!
!**************************  start rebuildStructures  ************************!
!#############################################################################!
!#############################################################################!
!******** This subroutine update the variables itend and it to allow *********!
!***************** the execution of the loop on the time steps ***************!
!#############################################################################!

  subroutine rebuildDPSubStructures(struct, size2, sizeArray, newStruct, mode)

    use commonStructures

    implicit none

    include 'femtime.h'
    include 'nlevel.h'

    integer size2, mode
    integer :: count_tag, i
    integer :: st(MPI_STATUS_SIZE), ierr
    integer :: sendbuffer
    double precision, dimension(1) :: struct, newStruct
    integer, dimension(n_threads) :: displsProc, sizeArray, recvbuffer, displs
    double precision, allocatable, dimension(:) :: fullStruct

    if(my_id .eq. 0) then
      if(mode .eq. 1) then  !structure based on elements
        allocate(fullStruct(size2)) 
      else if(mode .eq. 2) then  !structure based on nodes
        allocate(fullStruct(totalNodes))
      else
        write(6,*)'error in rebuildDPStructures', mode, my_id
        stop
      end if
      do i=1,n_threads
        recvbuffer(i) = sizeArray(i)
      end do
    end if

    sendbuffer = sizeArray(my_id+1)

    displs(1) = 0
    do i=2,n_threads
      displs(i) = displs(i-1) + sizeArray(i-1)
    end do

    call MPI_GATHERV(struct,sendbuffer,MPI_DOUBLE_PRECISION,fullStruct, &
                recvbuffer,displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

    if(my_id .eq. 0) then
      displsProc(1) = 0
      do i=2,n_threads
        displsProc(i) = displsProc(i-1) + sizeArray(i-1)
      end do
      if(mode .eq. 1) then
        call rebuildDPSubStruct(fullStruct, displsProc, newStruct, size2)
      else
        call rebuildDPSubStruct2(fullStruct, displsProc, newStruct, size2)
      end if
    end if

    return

  end subroutine

!#############################################################################!
!*************************  start rebuildRealStruct  *************************!
!#############################################################################!
!#############################################################################!
!******** This subroutine rebuild a complete structure after a parallel ******!
!**** subroutine, significant for the root process only. double precision ****!
!#############################################################################!

  subroutine rebuildDPSubStruct(fullStruct, displsProc, newStruct, size2)

    use commonStructures

    implicit none

    integer :: place, i, process, size2
    integer :: st(MPI_STATUS_SIZE), ierr
    integer, dimension(:) :: countProc(n_threads), displsProc(1)
    double precision, dimension(:) :: fullStruct, newStruct

    ! countProc(1) = counter for items of the master process
    do i=1,n_threads
      countProc(i)=1
    end do

    do i=1, size2
      process = allPartAssign(i)
      place = displsProc(process+1) + countProc(process+1)
      newStruct(i) = fullStruct(place)
      countProc(process+1) = countProc(process+1) + 1
    end do
   
    return 

  end subroutine

!#############################################################################!
!*************************  start rebuildRealSubStruct2 **********************!
!#############################################################################!
!#############################################################################!
!******** This subroutine rebuild a complete structure after a parallel ******!
!**** subroutine, significant for the root process only. double precision ****!
!#############################################################################!

  subroutine rebuildDPSubStruct2(fullStruct, displsProc, newStruct, size2)

    use commonStructures

    implicit none

    integer :: place, i, j, k, size2
    integer :: st(MPI_STATUS_SIZE), ierr
    integer, dimension(:) :: countProc(n_threads), displsProc(1)
    double precision, dimension(:) :: fullStruct, newStruct
    logical bexit 

    ! countProc(1) = counter for items of the master process
    do i=1,n_threads
      countProc(i)=1
    end do

    do i=1, size2
      bexit = .false.
      do j=1,n_threads
        if(allNodesAssign(j,i) .eq. 1) then
          place = displsProc(j) + countProc(j)
          countProc(j) = countProc(j) + 1
          if(bexit .eqv. .true.) then
            cycle
          end if
          newStruct(i) = fullStruct(place)
          bexit = .true.
        end if
      end do
    end do
   
    return 

  end subroutine

!#############################################################################!
!**************************  start rebuildStructures  ************************!
!#############################################################################!
!#############################################################################!
!******** This subroutine update the variables itend and it to allow *********!
!***************** the execution of the loop on the time steps ***************!
!#############################################################################!

  subroutine rebuildDPStructures(struct, size1, size2, sizeArray, newStruct, mode)

    use commonStructures

    implicit none

    include 'femtime.h'
    include 'nlevel.h'

    integer size1, size2, mode
    integer :: count_tag, i
!    integer :: nlvdi,nlv
    integer :: st(MPI_STATUS_SIZE), ierr
!    integer itanf,itend,idt,nits,niter,it
    integer :: sendbuffer
    double precision, dimension(size1,1) :: struct, newStruct
    integer, dimension(n_threads) :: displsProc, sizeArray, recvbuffer, displs
    double precision, allocatable, dimension(:,:) :: fullStruct
!    common /femtim/ itanf,itend,idt,nits,niter,it
!    common /level/ nlvdi,nlv

    if(my_id .eq. 0) then
      if(mode .eq. 1) then  !structure based on elements
        allocate(fullStruct(size1,size2)) 
      else if(mode .eq. 2) then  !structure based on nodes
        allocate(fullStruct(size1,totalNodes))
      else
        write(6,*)'error in rebuildDPStructures', mode, my_id
        stop
      end if
      do i=1,n_threads
        recvbuffer(i) = sizeArray(i) * size1
      end do
    end if

    sendbuffer = size1 * sizeArray(my_id+1)

    displs(1) = 0
    do i=2,n_threads
      displs(i) = displs(i-1) + size1 * sizeArray(i-1)
    end do

    call MPI_GATHERV(struct,sendbuffer,MPI_DOUBLE_PRECISION,fullStruct, &
                recvbuffer,displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

    if(my_id .eq. 0) then
      displsProc(1) = 0
      do i=2,n_threads
        displsProc(i) = displsProc(i-1) + sizeArray(i-1)
      end do
      if(mode .eq. 1) then
        call rebuildDPStruct(fullStruct, displsProc, newStruct, size1, size2)
      else
        call rebuildDPStruct2(fullStruct, displsProc, newStruct, size1, size2)
      end if
    end if

    return

  end subroutine


!#############################################################################!
!*************************  start rebuildDPStruct  ***************************!
!#############################################################################!
!#############################################################################!
!******** This subroutine rebuild a complete structure after a parallel ******!
!**** subroutine, significant for the root process only. double precision ****!
!#############################################################################!

  subroutine rebuildDPStruct(fullStruct, displsProc, newStruct, size1, size2)

    use commonStructures

    implicit none

    integer :: place, i, process, j, size1, size2
    integer :: st(MPI_STATUS_SIZE), ierr
    integer, dimension(:) :: countProc(n_threads), displsProc(1)
    double precision, dimension(:,:) :: fullStruct, newStruct

    ! countProc(1) = counter for items of the master process
    do i=1,n_threads
      countProc(i)=1
    end do

    do i=1, size2
      process = allPartAssign(i)
      place = displsProc(process+1) + countProc(process+1)
      do j=1,size1
        newStruct(j,i) = fullStruct(j,place)
      end do
      countProc(process+1) = countProc(process+1) + 1
    end do
   
    return 

  end subroutine


!#############################################################################!
!*************************  start rebuildDPStruct2  ***************************!
!#############################################################################!
!#############################################################################!
!******** This subroutine rebuild a complete structure after a parallel ******!
!**** subroutine, significant for the root process only. double precision ****!
!#############################################################################!

  subroutine rebuildDPStruct2(fullStruct, displsProc, newStruct, size1, size2)

    use commonStructures

    implicit none

    integer :: place, i, j, k, size1, size2
    integer :: st(MPI_STATUS_SIZE), ierr
    integer, dimension(:) :: countProc(n_threads), displsProc(1)
    double precision, dimension(:,:) :: fullStruct, newStruct
    logical bexit

    ! countProc(1) = counter for items of the master process
    do i=1,n_threads
      countProc(i)=1
    end do

    do i=1, size2
      bexit = .false.
      do j=1,n_threads
        if(allNodesAssign(j,i) .eq. 1) then
          place = displsProc(j) + countProc(j)
          countProc(j) = countProc(j) + 1
          if(bexit .eqv. .true.) then
            cycle
          end if
          do k=1, size1
            newStruct(k,i) = fullStruct(k,place)
          end do
          bexit = .true.
        end if
      end do
    end do
   
    return 

  end subroutine


!#############################################################################!
!*************************  start rebuildRealStruct  *************************!
!#############################################################################!
!#############################################################################!
!******** This subroutine rebuild a complete structure after a parallel ******!
!**** subroutine, significant for the root process only. double precision ****!
!#############################################################################!

  subroutine rebuildRealSubStruct(fullStruct, displsProc, newStruct, size2)

    use commonStructures

    implicit none

    integer :: place, i, process, size2
    integer :: st(MPI_STATUS_SIZE), ierr
    integer, dimension(:) :: countProc(n_threads), displsProc(1)
    real, dimension(:) :: fullStruct, newStruct

    ! countProc(1) = counter for items of the master process
    do i=1,n_threads
      countProc(i)=1
    end do

    do i=1, size2
      process = allPartAssign(i)
      place = displsProc(process+1) + countProc(process+1)
      newStruct(i) = fullStruct(place)
      countProc(process+1) = countProc(process+1) + 1
    end do
   
    return 

  end subroutine


!#############################################################################!
!*************************  start rebuildRealSubStruct2 **********************!
!#############################################################################!
!#############################################################################!
!******** This subroutine rebuild a complete structure after a parallel ******!
!**** subroutine, significant for the root process only. double precision ****!
!#############################################################################!

  subroutine rebuildRealSubStruct2(fullStruct, displsProc, newStruct, size2)

    use commonStructures

    implicit none

    integer :: place, i, j, k, size2
    integer :: st(MPI_STATUS_SIZE), ierr
    integer, dimension(:) :: countProc(n_threads), displsProc(1)
    real, dimension(:) :: fullStruct, newStruct
    logical bexit 

    ! countProc(1) = counter for items of the master process
    do i=1,n_threads
      countProc(i)=1
    end do

    do i=1, size2
      bexit = .false.
      do j=1,n_threads
        if(allNodesAssign(j,i) .eq. 1) then
          place = displsProc(j) + countProc(j)
          countProc(j) = countProc(j) + 1
          if(bexit .eqv. .true.) then
            cycle
          end if
          newStruct(i) = fullStruct(place)
          bexit = .true.
        end if
      end do
    end do
   
    return 

  end subroutine


   subroutine gatherStruct_r(array)

      use basin
      use levels
      use commonStructures

      implicit none

      real array(nlvdi,nkn)
      real, allocatable, dimension(:,:) :: temp,tempArray

      if(my_id .eq. 0) then
        allocate(tempArray(nlvdi,nkndi))
      end if
      allocate(temp(nlvdi,mynodes%totalID))



      call rebuildRealStructures(array,nlvdi,nkndi,numberNodes,tempArray, 2)

      if(my_id .ne. 0) then
        call setRealStruct(temp, mynodes, temp, nlvdi, 2)
      else
        call setRealStruct(tempArray, mynodes, temp, nlvdi, 2)
      end if

      array = temp

      return

   end subroutine





!######################################################################!
!*******************  start MergeSort subroutine  *********************!
!######################################################################!
!######################################################################!
!* This subroutine implement a mergesort algorithm for integer vector *!
!######################################################################!

  recursive subroutine MergeSort(array,length,T)

   implicit none

   !argument
   integer, intent(in) :: length
   integer, dimension(length), intent(in out) :: array

   integer, dimension((length+1)/2) :: T

   integer :: NA,NB,V

   if (length < 2) return
   if (length == 2) then
      if (array(1) > array(2)) then
         V = array(1)
         array(1) = array(2)
         array(2) = V
      endif
      return
   endif
   NA=(length+1)/2
   NB=length-NA

   call MergeSort(array,NA,T)
   call MergeSort(array(NA+1),NB,T)

   if (array(NA) > array(NA+1)) then
      T(1:NA)=array(1:NA)
      call Merge(T,NA,array(NA+1),NB,array,length)
   endif
   return

  end subroutine MergeSort

  subroutine Merge(A,NA,B,NB,C,NC)

   integer, intent(in) :: NA,NB,NC         ! Normal usage: NA+NB = NC
   integer, intent(in out) :: A(NA)        ! B overlays C(NA+1:NC)
   integer, intent(in)     :: B(NB)
   integer, intent(in out) :: C(NC)

   integer :: I,J,K

   I = 1; J = 1; K = 1;
   do while(I <= NA .and. J <= NB)
      if (A(I) <= B(J)) then
         C(K) = A(I)
         I = I+1
      else
         C(K) = B(J)
         J = J+1
      endif
      K = K + 1
   enddo
   do while (I <= NA)
      C(K) = A(I)
      I = I + 1
      K = K + 1
   enddo
   return

  end subroutine merge


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
!********************  end Module Global_Graph_Data  ******************!
!######################################################################!

  end module Graph_Ele
