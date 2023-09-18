
!#######################################################################!
!*******************  start module Global_Graph_Data  ******************!
!#######################################################################!

module mpi_graph_elem

  !use mpi
  use mpi_communication_struct
  use basin
  use timing
  use shympi


    logical, dimension(:,:),allocatable :: rank_neighbor

contains


!#######################################################################!
!*********************  start makeelems subroutine  ********************!
!#######################################################################!
!#######################################################################!
!* This subroutine makes the domain%elems structure. It contains informations *!
!* on the elements of a subprocess(global ID elements, rank ID process *!
!* neighbor, globalID neighbor elements, etc.)                         *!
!#######################################################################!

  subroutine makeMyEle(numMyVertices, numGlobalVertices, struct)

    use mpi_common_struct
    use geom
    use levels
    use shympi

    implicit none

    ! input
    ! rank process, number of my elements, number global elements, number global of neighbor
    integer numMyVertices, numGlobalVertices

    ! allPartAssign contains the ID of the process to which is assigned the element
    ! the number of the neighbors an element "nneighbor" is equal to
    integer, dimension(3,numGlobalVertices) :: struct 
    integer, dimension(:), allocatable :: neighbor, sort_neighbor   
    integer, dimension(:), allocatable :: sortRank, rank
    integer, dimension(:), allocatable :: part_receives, part_sends   
    integer :: st(MPI_STATUS_SIZE), ierr
    integer error, STAT
 
    ! local 
    integer eleneighbor, myneighbor, nneighbor, globalID, locID, neighbID
    integer x,p,s,i,j,k,h,n,ele,t,length,nprocesses,node,maxSends,maxReceives,sizev
    integer ele_n,jj
    integer numlevels, totalLevels
    double precision time1
#ifdef DEBUGON
    integer,dimension(:),allocatable :: GID_order
#endif
    allocate(rank_neighbor(0:n_threads-1,numGlobalVertices))


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
    allocate(part_receives(sizev))
    allocate(part_sends(sizev))
    allocate(rank(sizev))
    allocate(sortRank(sizev))

    do i=1,sizev
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
    domain%elems%domainID=0
    do i=1,numGlobalVertices
       if(allPartAssign(i) .eq. my_id) then

          do h=1,3

            ele_n=ieltv(h,i)
            if(ele_n .ne. 0) then
              if (my_id .ne. allPartAssign(ele_n)) then
                rank_neighbor(allPartAssign(ele_n),ele_n)=.true.
              end if
            end if


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

          domain%elems%domainID=domain%elems%domainID+1
       end if
    end do

    if(domain%elems%domainID .ne. numMyVertices)then
      write(*,*)'Error stop makeelems'
      stop
    end if

    allocate(numberElements(n_threads))

    if(ln_timing) time1 = shympi_wtime()

    call MPI_ALLGATHER(domain%elems%domainID, 1, MPI_INTEGER, numberElements, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    if(ln_timing) comm_time = comm_time + shympi_wtime() - time1

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

    allocate(domain%process(nprocesses),STAT=error)
    allocate(domain%ele_send(nprocesses),STAT=error)
    allocate(domain%ele_receive(nprocesses),STAT=error)
    allocate(domain%ele_temp(nprocesses),STAT=error)

    domain%exchanges = nprocesses

    do i=1,nprocesses
       domain%process(i) = sortRank(i)
    end do

    deallocate(sortRank)

    domain%elems%totalID = numMyVertices + length

    if(domain%elems%totalID .ne. 0) then
#ifdef DEBUGON
      allocate(domain%elems%mapID(domain%elems%totalID))
#endif
      allocate(domain%elems%globalID(domain%elems%totalID),STAT=error)
        if (error .ne. 0) then
          write(6,*)'error: could not allocate memory for array',error
          stop
        endif
    end if

    locID=1
    do i=1,numGlobalVertices
       if(allPartAssign(i) .eq. my_id) then
          domain%elems%globalID(locID) = i
          locID = locID + 1
       end if
    end do

    do i=numMyVertices+1,domain%elems%totalID
      domain%elems%globalID(i) = sort_neighbor(i-numMyVertices)
    end do

#ifdef DEBUGON
    allocate(GID_order(domain%elems%totalID))
    call remove_dups(domain%elems%globalID,length,GID_order)

    do i=1,domain%elems%totalID
myloop:   do j=1,domain%elems%totalID
            if(domain%elems%globalID(j).eq.GID_order(i)) then
              domain%elems%mapID(i)=j
              exit myloop
            end if
          end do myloop
    end do
#endif
    do h=1,nprocesses

       k = 0
       s = 0     
       part_receives = 0
       part_sends = 0

       do i=numMyVertices+1,domain%elems%totalID

          globalID=domain%elems%globalID(i)
          if(globalID .gt. 0) then
            if((allPartAssign(globalID) .eq. domain%process(h)) .and. &
           &  (rank_neighbor(domain%process(h),domain%elems%globalID(i)))) then
              k = k +1
              part_receives(k) = globalID
            end if
          end if

       end do

       do ele=1,domain%elems%domainID
          do n=1,3

            globalID=ieltv(n,domain%elems%globalID(ele))
            if(globalID .gt. 0) then
              if(allPartAssign(globalID) .eq. domain%process(h)) then
                s = s + 1
                part_sends(s) = domain%elems%globalID(ele)
              end if
            end if

          end do
       end do

       allocate(domain%ele_temp(h)%items(k))
       domain%ele_temp(h)%numItems = k

       if(k.gt.0) then
         call remove_dups(part_receives,length,domain%ele_temp(h)%items)
       else
        length=0
       end if

       allocate(domain%ele_receive(h)%items(length))
       domain%ele_receive(h)%numItems = length
       do n=1,length
loop1:   do jj=domain%elems%domainID,domain%elems%totalID
           if(domain%elems%globalID(jj).eq.domain%ele_temp(h)%items(n)) then
             domain%ele_receive(h)%items(n)=jj
             exit loop1
           end if
         end do loop1
       end do

       deallocate(domain%ele_temp(h)%items)
       allocate(domain%ele_temp(h)%items(s))
       domain%ele_temp(h)%numItems = s
       if(s .gt. 0) then
         call remove_dups(part_sends,length,domain%ele_temp(h)%items)
       else
         length = 0
       end if

       allocate(domain%ele_send(h)%items(length))

       domain%ele_send(h)%numItems = length

       do n=1,length
loop2:   do jj=1,domain%elems%domainID
           if(domain%elems%globalID(jj).eq.domain%ele_temp(h)%items(n)) then
             domain%ele_send(h)%items(n) = jj
             exit loop2
           end if
         end do loop2
       end do

    end do

    deallocate(domain%ele_temp)
    deallocate(sort_neighbor)
    !deallocate(receives)
    !deallocate(sends)
    deallocate(part_receives)
    deallocate(part_sends)

    return

  end subroutine

!#######################################################################!
!********************  start makenodes subroutine  *******************!
!#######################################################################!
!#######################################################################!
! This subroutine makes the domain%nodes structure. It contains informations !
!***  on the nodes of a subprocess(global ID nodes, rank ID process  ***!
!*****          neighbor, globalID neighbor nodes, etc.)           *****!
!#######################################################################!

  subroutine makeMyNodes(mystruct, struct, nkn)

    use mpi_common_struct
    use geom
    use levels

    implicit none


    logical debug
    ! input

    integer nkn
    integer mystruct(3,domain%elems%domainID)
    integer, dimension(:,:) :: struct
    integer, dimension(:), allocatable :: neighbor, sort_neighbor, sort_nodes
    integer, dimension(:), allocatable :: fullNodesAssign
    integer, dimension(n_threads) :: recvbuffer, displs
    type (MAPPING_INFO) :: temp
    integer :: st(MPI_STATUS_SIZE), ierr 
    ! local 
    integer i, j, k, n, x, h, s, length, sendbuffer, node,nneighbor
    integer counter

    !
    integer ounit,error
    character*(20) fmt1
    character*(20) fmt2
    integer sendbuffer2
    integer, dimension(n_threads) :: recvbuffer2, displs2

    integer nnodbound,nbound_num
    integer, dimension(:), allocatable :: nodbound, sort_nodbound

    integer tempvar,ele,eleneighbor,rank,numnodes,n_node

    integer jj   
    character*(20) filename
    character*(20) format_string 
    double precision time1

    counter = 0

    allocate(temp%globalID(3*domain%elems%domainID))
    allocate(sort_nodes(3*domain%elems%domainID))

    n=0
    do i=1, domain%elems%domainID
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

    domain%nodes%domainID = length

    allocate(numberNodes(n_threads))
    allocate(procNodes(n_threads))

    if(ln_timing) time1 = shympi_wtime()

    call MPI_ALLGATHER(domain%nodes%domainID, 1, MPI_INTEGER, numberNodes, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    if(ln_timing) comm_time = comm_time + shympi_wtime() - time1

    procNodes=numberNodes


    totalNodes = numberNodes(1)
    do i=2,n_threads
      totalNodes = totalNodes + numberNodes(i)
    end do

    if(my_id .eq. 0)then
      write(6,*)'totalNodes=',totalNodes
    end if

    allocate(allNodesAssign(nkn,n_threads))
    allocate(fullNodesAssign(totalNodes))


    do i=1,n_threads
      recvbuffer(i) = numberNodes(i)
    end do

    sendbuffer = numberNodes(my_id+1)

    displs(1) = 0
    do i=2,n_threads
      displs(i) = displs(i-1) + numberNodes(i-1)
    end do

    if(ln_timing) time1 = shympi_wtime()

    call MPI_ALLGATHERV(sort_nodes, sendbuffer, MPI_INTEGER, fullNodesAssign, recvbuffer, displs, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(size(sort_nodes),maxNodes,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD, ierr)

    if(ln_timing) comm_time = comm_time + shympi_wtime() - time1

    call makeNodesAssign(fullNodesAssign, nkn)
    deallocate(fullNodesAssign)

    allocate(neighbor(3*(domain%elems%totalID-domain%elems%domainID)))
    allocate(sort_neighbor(3*(domain%elems%totalID-domain%elems%domainID)))

    do i=1,3*(domain%elems%totalID-domain%elems%domainID)
      neighbor(i) = 0
    end do

    nneighbor = 0
    do i=domain%elems%domainID+1,domain%elems%totalID
      do j=1,3
        node = struct(j,domain%elems%globalID(i))
        if(allNodesAssign(node,my_id+1) .ne. 1) then
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

    domain%nodes%totalID = domain%nodes%domainID + length

    allocate(nodbound(3*domain%elems%domainID))
    allocate(sort_nodbound(3*domain%elems%domainID))

    nodbound = 0

    nnodbound=0
    do i=1,domain%elems%domainID
      do j=1,3
        node = struct(j,domain%elems%globalID(i))
        do k=my_id+2,n_threads
          if(allNodesAssign(node,k) .eq. 1) then
            nnodbound = nnodbound + 1
            nodbound(nnodbound) = node
            exit
          end if
        end do

      end do
    end do

    if(nnodbound .gt. 0) then
      call remove_dups(nodbound,nbound_num,sort_nodbound)
    else
      nbound_num = 0
    end if

    allocate(domain%nodes%globalID(domain%nodes%totalID))
 

    counter=0 
    do i=1,domain%nodes%domainID
      do j=1,nbound_num
        if(sort_nodbound(j) .eq. sort_nodes(i)) exit
        if(j .eq. nbound_num ) then
          counter = counter + 1
          domain%nodes%globalID(counter) = sort_nodes(i)
        end if
      end do
    end do

    if(nbound_num .eq. 0) then 
      do counter=1,domain%nodes%domainID
        domain%nodes%globalID(counter) = sort_nodes(counter)
      end do
      counter = counter - 1
    else
      do i=counter+1,domain%nodes%domainID
        domain%nodes%globalID(i) = sort_nodbound(i-counter)
      end do
    end if 

    if(counter .ne. (domain%nodes%domainID-nbound_num))then
       write(6,*)'error in mod_bound:',counter,nbound_num,domain%nodes%domainID
       stop
    end if

    domain%nodes%univocalID = counter
    nkn_inner = counter

    do i=domain%nodes%domainID+1,domain%nodes%totalID
      domain%nodes%globalID(i) = sort_neighbor(i-domain%nodes%domainID)
    end do

    call setUnivocal

    allocate(domain%node_send(domain%exchanges))
    allocate(domain%node_temp(domain%exchanges))
    do h=1,domain%exchanges
      allocate(domain%node_temp(h)%items((domain%nodes%totalID-domain%nodes%domainID)*8))
      do i=1,(domain%nodes%totalID-domain%nodes%domainID)*8
        domain%node_temp(h)%items(i)=0
      end do
      domain%node_temp(h)%numItems=0
    end do

      do i=domain%elems%domainID+1, domain%elems%totalID
        do j=1,3
          k=struct(j,domain%elems%globalID(i))
          do n=1,domain%nodes%domainID                
            if(k .eq. domain%nodes%globalID(n)) then
              do h=1,domain%exchanges
                if(allNodesAssign(k,domain%process(h)+1) .eq. 1) then
                    domain%node_temp(h)%numItems=domain%node_temp(h)%numItems+1
                    domain%node_temp(h)%items(domain%node_temp(h)%numItems) = n
                end if
              end do
            end if
          end do
        end do                             
      end do

    domain%maxItems=0

    do h=1,domain%exchanges
       if(domain%node_temp(h)%numItems .gt. 0)then
          call remove_dups(domain%node_temp(h)%items,length,domain%node_send(h)%items &
       ,domain%node_temp(h)%numItems)
       else
          length = 0
       end if
       domain%node_send(h)%numItems=length
       domain%maxItems=max(domain%maxItems,length)
       deallocate(domain%node_temp(h)%items)

       do i=1,domain%node_send(h)%numItems-1
         j=i
         do while(domain%nodes%globalID(domain%node_send(h)%items(j+1)) .lt. & 
         domain%nodes%globalID(domain%node_send(h)%items(j)) .and. (j .ge. 1))
           tempvar = domain%node_send(h)%items(j+1)
           domain%node_send(h)%items(j+1) = domain%node_send(h)%items(j)
           domain%node_send(h)%items(j) = tempvar
           if(j .eq. 1) exit
           j = j-1
         end do
       end do
    end do

    allocate(domain%halo_temp(domain%exchanges))
    allocate(domain%halo_temp2(domain%exchanges))
    allocate(domain%halo_temp3(domain%exchanges))
    allocate(domain%halo_temp4(domain%exchanges))
    allocate(domain%halo_recv_e(domain%exchanges))
    allocate(domain%halo_send_e(domain%exchanges))
    do h=1,domain%exchanges
      domain%halo_recv_e(h)%numItems=0
      domain%halo_send_e(h)%numItems=0
      k=0
      s=0
      allocate(domain%halo_temp(h)%items((domain%nodes%totalID-domain%nodes%domainID)*8))
      allocate(domain%halo_temp2(h)%items((domain%nodes%totalID-domain%nodes%domainID)*8))
      allocate(domain%halo_temp3(h)%items((domain%nodes%totalID-domain%nodes%domainID)*8))
      allocate(domain%halo_temp4(h)%items((domain%nodes%totalID-domain%nodes%domainID)*8))
      domain%halo_temp(h)%items=0
      domain%halo_temp2(h)%items=0
      domain%halo_temp3(h)%items=0
      domain%halo_temp4(h)%items=0
      rank=domain%process(h)
      do i=1,domain%node_send(h)%numItems
        node=domain%nodes%globalID(domain%node_send(h)%items(i))
      
        eleneighbor=ilinkv(node+1)-ilinkv(node)
        do j=1,eleneighbor
        ele = lenkv(ilinkv(node) + j)
        if (ele .gt. 0) then
          if(allPartAssign(ele) .ne. my_id) then
            !if (rank_neighbor(rank,ele)) then
            if (allPartAssign(ele).eq.rank) then
              s=s+1
!              do n=domain%elems%domainID,domain%elems%totalID                
!                if(ele .eq. domain%elems%globalID(n)) then
                  domain%halo_temp(h)%items(s)=ele
!                end if
!              end do
            end if
          else
            if(allNodesAssign(node,(rank+1)) .eq. 1) then
              k = k + 1
!              do n=1,domain%elems%domainID              
!                if(ele .eq. domain%elems%globalID(n)) then
                  domain%halo_temp2(h)%items(k)=ele
!                end if
!              end do
            end if
          end if
        end if
        end do
      end do
      domain%halo_temp(h)%numItems=s
      domain%halo_temp2(h)%numItems=k
      if(s .gt. 0)then
        call remove_dups(domain%halo_temp(h)%items,length,domain%halo_temp3(h)%items &
                        ,domain%halo_temp(h)%numItems)
        domain%halo_temp3(h)%numItems=length
        domain%maxItems=max(domain%maxItems,length)
      else
        domain%halo_temp3(h)%numItems=0
        domain%halo_temp3(h)%items=-1
      end if
      if(k .gt. 0)then

        call remove_dups(domain%halo_temp2(h)%items,length,domain%halo_temp4(h)%items &
                        ,domain%halo_temp2(h)%numItems)
        domain%halo_temp4(h)%numItems=length
        domain%maxItems=max(domain%maxItems,length)
      else
        domain%halo_temp4(h)%numItems=0
        domain%halo_temp4(h)%items=-1
      end if
    end do

    do h=1,domain%exchanges
      allocate(domain%halo_recv_e(h)%items(domain%halo_temp3(h)%numItems))
      allocate(domain%halo_send_e(h)%items(domain%halo_temp4(h)%numItems))
      domain%halo_recv_e(h)%numItems=domain%halo_temp3(h)%numItems
      domain%halo_send_e(h)%numItems=domain%halo_temp4(h)%numItems

      do i=1,domain%halo_temp3(h)%numItems
loop1:  do k=domain%elems%domainID,domain%elems%totalID
          if(domain%halo_temp3(h)%items(i).eq.domain%elems%globalID(k)) then
            domain%halo_recv_e(h)%items(i)=k
            exit loop1
          end if
        end do loop1
      end do

      do i=1,domain%halo_temp4(h)%numItems
loop2:  do k=1,domain%elems%domainID
          if(domain%halo_temp4(h)%items(i).eq.domain%elems%globalID(k)) then
            domain%halo_send_e(h)%items(i)=k
            exit loop2
          end if
        end do loop2
      end do
    end do

    allocate(domain%halo_recv_k(domain%exchanges))
    allocate(domain%halo_send_k(domain%exchanges))
    do h=1,domain%exchanges
      domain%halo_recv_k(h)%numItems=0
      domain%halo_send_k(h)%numItems=0
      k=0
      s=0
      domain%halo_temp(h)%items=0
      domain%halo_temp2(h)%items=0
      rank=domain%process(h)
      do i=1,domain%node_send(h)%numItems
        node=domain%nodes%globalID(domain%node_send(h)%items(i))
      
        numnodes=ilinkv(node+1)-ilinkv(node)
        do j=1,numnodes
        n_node = linkv(ilinkv(node) + j)
        if (n_node .gt. 0) then
          if(allNodesAssign(n_node,(my_id+1)) .ne. 1) then
            if(allNodesAssign(n_node,(rank+1)) .eq. 1) then
              s=s+1
!              do n=domain%nodes%domainID,domain%nodes%totalID                
!                if(n_node .eq. domain%nodes%globalID(n)) then
                  domain%halo_temp(h)%items(s)=n_node
!                end if
!              end do
            end if
          else
            if(allNodesAssign(n_node,(rank+1)) .ne. 1) then
              k = k + 1
!              do n=1,domain%nodes%domainID              
!                if(n_node .eq. domain%nodes%globalID(n)) then
                  domain%halo_temp2(h)%items(k)=n_node
!                end if
!              end do
            end if
          end if
        end if
        end do
      end do
      domain%halo_temp(h)%numItems=s
      domain%halo_temp2(h)%numItems=k
      if(s .gt. 0)then
        call remove_dups(domain%halo_temp(h)%items,length,domain%halo_temp3(h)%items &
                        ,domain%halo_temp(h)%numItems)
        domain%maxItems=max(domain%maxItems,length)
        domain%halo_temp3(h)%numItems=length
        domain%maxItems=max(domain%maxItems,length)
      else
        domain%halo_temp3(h)%numItems=0
        domain%halo_temp3(h)%items=-1
      end if
      if(k .gt. 0)then
        call remove_dups(domain%halo_temp2(h)%items,length,domain%halo_temp4(h)%items &
                        ,domain%halo_temp2(h)%numItems)
        domain%halo_temp4(h)%numItems=length
        domain%maxItems=max(domain%maxItems,length)
      else
        domain%halo_temp4(h)%numItems=0
        domain%halo_temp4(h)%items=-1
      end if
    end do

    do h=1,domain%exchanges
      allocate(domain%halo_recv_k(h)%items(domain%halo_temp3(h)%numItems))
      allocate(domain%halo_send_k(h)%items(domain%halo_temp4(h)%numItems))
      domain%halo_recv_k(h)%numItems=domain%halo_temp3(h)%numItems
      domain%halo_send_k(h)%numItems=domain%halo_temp4(h)%numItems

      do i=1,domain%halo_temp3(h)%numItems
loop3:  do k=domain%nodes%domainID,domain%nodes%totalID
            if(domain%halo_temp3(h)%items(i).eq.domain%nodes%globalID(k)) then
            domain%halo_recv_k(h)%items(i)=k
            exit loop3
          end if
        end do loop3
      end do

      do i=1,domain%halo_temp4(h)%numItems
loop4:  do k=1,domain%nodes%domainID
          if(domain%halo_temp4(h)%items(i).eq.domain%nodes%globalID(k)) then
            domain%halo_send_k(h)%items(i)=k
            exit loop4
          end if
        end do loop4
      end do
    end do

!                       gmica 216 mpi process debug printing
!    do h=1,domain%exchanges   
!    if (my_id.eq.170.and.domain%process(h).eq.143) then
!      do jj=1,domain%halo_recv_k(h)%numitems  
!        if (domain%halo_recv_k(h)%items(jj).gt.0) then    
!        !write(6,*) 'recv',ipv(domain%nodes%globalID(domain%halo_recv_k(h)%items(jj)))
!        write(6,*) 'recv',domain%nodes%globalID(domain%halo_recv_k(h)%items(jj))
!        !write(6,*) 'recv',domain%halo_recv_k(h)%items(jj)
!        end if  
!      end do          
!    end if          
!    end do    
!
!    do h=1,domain%exchanges   
!    if (my_id.eq.143.and.domain%process(h).eq.170) then
!      do jj=1,domain%halo_send_k(h)%numitems
!        if (domain%halo_send_k(h)%items(jj).gt.0) then    
!        !write(6,*) 'send',ipv(domain%nodes%globalID(domain%halo_send_k(h)%items(jj)))
!        write(6,*) 'send',domain%nodes%globalID(domain%halo_send_k(h)%items(jj))
!        !write(6,*) 'send',domain%halo_send_k(h)%items(jj)
!        end if  
!      enddo  
!    end if          
!    end do    

         domain%pos=0
myloop:  do i=1,domain%exchanges
            if(domain%process(i).lt.my_id) then
              domain%pos=i
            else
              exit myloop
            end if
         end do myloop

!    if(my_id.eq.0) then
!      open(unit=150, file="nen3v.txt", action='write')
!      open(unit=151, file="ilinkv.txt", action='write')
!      open(unit=152, file="lenkv.txt", action='write')
!      open(unit=153, file="linkv.txt", action='write')
!      open(unit=154, file="extel.txt", action='write')
!      open(unit=155, file="extkn.txt", action='write')
!      do i=1,30
!        write(150,*) nen3v(:,i)
!      end do
!      do i=1,31
!        write(151,*)ilinkv(i)
!      end do
!      write(152,*)'    k    ln    1°e   2°e   3°e   4°e   5°e   6°e   7°e   8°e'
!      write(153,*)'    k    ln    1°k   2°k   3°k   4°k   5°k   6°k   7°k   8°k'
!      do i=1,30
!        n=ilinkv(i+1)-ilinkv(i)
!        fmt2='(A,I1,A3)'
!        write(fmt1,fmt2)'(',n+2,'I6)'
!        write(152,fmt1)i,n,lenkv(ilinkv(i)+1:ilinkv(i+1))
!        write(153,fmt1)i,n,linkv(ilinkv(i)+1:ilinkv(i+1))
!      end do
!      do i=1,30
!        write(154,*) i,ipev(i)
!        write(155,*) i,ipv(i)
!      end do
!      close(150)
!      close(151)
!      close(152)
!      close(153)
!      close(154)
!      close(155)
!    end if

    deallocate(domain%node_temp)
    deallocate(domain%halo_temp)
    deallocate(domain%halo_temp2)
    deallocate(domain%halo_temp3)
    deallocate(domain%halo_temp4)
    deallocate(rank_neighbor)
    deallocate(allNodesAssign)

    deallocate (sort_nodes)
    deallocate (neighbor)
    deallocate (sort_neighbor)
    deallocate(temp%globalID)

    return

  end subroutine

!**************************************************************

     subroutine setup_bnd_neighbors(nbndo,nopnod,nopnodes)


        use basin
        use mpi_common_struct
        use mpi_communication_struct

	implicit none

        integer nbndo
        integer i,ierr,nsend,nrecv        
	integer ib		!nodal index
	integer kn		!node number to insert
        integer nopnod(:)
        integer nopnodes(:,:)
        type (commItem), dimension(:), allocatable :: temp,temp2
        integer length,counter
        integer :: send(n_threads),recv(n_threads)
	integer gID,nb,j,k,rank
        integer :: req(n_threads)
        integer st(MPI_STATUS_SIZE)
        double precision time1

        allocate(temp(n_threads))
        allocate(temp2(n_threads))
        do j=1,n_threads
          allocate(temp(j)%items(bounds%nneighbors))
          temp(j)%items=0
          temp(j)%numItems=0
          temp2(j)%numItems=0
        end do
        do ib=1,nbndo
	  nb = nopnod(ib)
	  do j=1,nb
	    kn = nopnodes(j,ib)
          !if(my_id.eq.170)write(6,*)'kn:',kn,ipv(kn),ib,j,nkn
            if(kn.gt.nkn) then
              gID=domain%nodes%globalID(kn)
              rank=univocalNodesAssign(gID)
              temp(rank+1)%prank=rank
              temp(rank+1)%numItems=temp(rank+1)%numItems+1
              temp(rank+1)%items(temp(rank+1)%numItems)=gID
          if(my_id.eq.143)write(6,*)'counter' &
             ,temp(rank+1)%numItems,gID,ipv(kn)
            end if
          end do
        end do

        do j=1,n_threads
          if(temp(j)%numItems.gt.0) then
            call remove_dups(temp(j)%items,length,temp2(j)%items &
              ,temp(j)%numItems)
          else
            length=0
          end if
          temp2(j)%numItems=length
        end do

        counter=0
        send=0
        do j=1,n_threads
          if(temp2(j)%numItems.gt.0) then
            counter=counter+1
            send(j)=send(j)+temp2(j)%numItems
          end if
        end do
        bounds%recv=counter

        if(ln_timing) time1 = shympi_wtime()

        call mpi_Alltoall(send,1,MPI_INTEGER,recv, &
               & 1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

        if(ln_timing) comm_time = comm_time + shympi_wtime() - time1

        nrecv=0
        do j=1,n_threads
          if(recv(j).ne.0) nrecv= nrecv+1
        end do
        bounds%send=nrecv

        allocate(bounds%snode(nrecv))

        nrecv=0
        do j=1,n_threads
          if(recv(j).gt.0) then
            nrecv=nrecv+1
            bounds%snode(nrecv)%numItems=recv(j)
            bounds%snode(nrecv)%prank=j-1
            allocate(bounds%snode(nrecv)%items(recv(j)))
          end if
        end do

        allocate(bounds%rnode(counter))
        counter=0
        do j=1,n_threads
          if(temp2(j)%numItems.gt.0) then
            counter=counter+1
            allocate(bounds%rnode(counter)%items(temp2(j)%numItems))
            do ib=1,temp2(j)%numItems
loop1:        do i=domain%nodes%domainID,domain%nodes%totalID
                if(temp2(j)%items(ib).eq.domain%nodes%globalID(i)) then
                  bounds%rnode(counter)%items(ib)=i
                  exit loop1
                end if
              end do loop1
            end do
            bounds%rnode(counter)%numItems=temp2(j)%numItems
            bounds%rnode(counter)%prank=temp(j)%prank
          end if
        end do

        if(ln_timing) time1 = shympi_wtime()

        do j=1,bounds%send
          rank=bounds%snode(j)%prank
          nrecv=bounds%snode(j)%numItems
          call MPI_IRecv(bounds%snode(j)%items,nrecv,MPI_INTEGER,rank &
     &        ,MPI_ANY_TAG,MPI_COMM_WORLD,req(rank+1),ierr)
        end do

        counter=0
        do j=1,n_threads
          if(temp2(j)%numItems.gt.0) then
            counter=counter+1
            rank=temp(j)%prank
            nsend=temp2(j)%numItems
            call MPI_Send(temp2(j)%items,nsend,MPI_INTEGER, &
     &                 rank,rank,MPI_COMM_WORLD,ierr)
          end if
        end do

        do j=1,bounds%send
          rank=bounds%snode(j)%prank
          CALL MPI_Wait(req(rank+1), st, ierr)
        end do

        if(ln_timing) comm_time = comm_time + shympi_wtime() - time1

        do j=1,bounds%send
          counter=1
loop2:    do i=1,domain%nodes%domainID
            if(domain%nodes%globalID(i).eq. &
     &          bounds%snode(j)%items(counter)) then
              bounds%snode(j)%items(counter)=i
              counter=counter+1
              exit loop2
            end if
          end do loop2
        end do

        do j=1,bounds%recv
          write(6,*)'recv_bn',ipv(bounds%rnode(j)%items(:)), &
                               bounds%rnode(j)%prank
        end do

        do j=1,bounds%send
          write(6,*)'send_bn',ipv(bounds%snode(j)%items(:)), &
                               bounds%snode(j)%prank
        end do

      return

    end subroutine

!########################################################################################################!

   subroutine setUnivocal

     implicit none
     integer i,j,k,test,n
     integer total_nod,ierr
     integer, allocatable,dimension(:) :: nnodes,nnodesAssign
     integer, dimension(n_threads) :: sendbuffer,recvbuffer, displs
     double precision temp(maxNodes)
     integer, allocatable,dimension(:) :: reorderNodes
     integer r
     double precision time1


     if(ln_timing) time1 = shympi_wtime()

     call MPI_ALLREDUCE(domain%nodes%univocalID, total_nod, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)


     allocate(nnodes(n_threads))
     allocate(nnodesAssign(total_nod))

      if(size(nnodesAssign) .ne. nkndi) then
       write(6,*)'error: univocal_nodes dimension',nkndi,size(nnodesAssign)
     end if

    !write(6,*)'checkMyNodes:',nkndi,total_nod

     call MPI_AllGATHER(domain%nodes%univocalID, 1, MPI_INTEGER, nnodes, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

     if(ln_timing) comm_time = comm_time + shympi_wtime() - time1

     if(bmpi_debug) then
     do i=1,n_threads
        recvbuffer(i) = nnodes(i)
     end do

     sendbuffer = nnodes(my_id+1)

     displs(1) = 0
     do i=2,n_threads
       displs(i) = displs(i-1) + nnodes(i-1)
     end do

     if(ln_timing) time1 = shympi_wtime()

     call MPI_ALLGATHERV(domain%nodes%globalID, sendbuffer, MPI_INTEGER, &
          nnodesAssign, recvbuffer, displs, MPI_INTEGER, MPI_COMM_WORLD, ierr)

     if(ln_timing) comm_time = comm_time + shympi_wtime() - time1

       do i=1,n_threads
          if( i .ne. (my_id+1)) then
             do j=1,domain%nodes%univocalID
                do k=1, nnodes(i)
                   if(domain%nodes%globalID(j) .eq. nnodesAssign(displs(i)+k)) then
                      write(6,*)'error: univocal_nodes',my_id
                      write(6,*),domain%nodes%globalID(j),i-1,my_id
                   end if
                end do
             end do
          end if
       end do
     end if

     allocate(univocalNodesAssign(nkndi))
     !allocate(reorderNodes(nkndi))

     nnodesAssign=-1

     do i=1,domain%nodes%domainID
        nnodesAssign(domain%nodes%globalID(i))=my_id
     end do
       
     if(ln_timing) time1 = shympi_wtime()

      call MPI_ALLREDUCE(nnodesAssign, univocalNodesAssign, nkndi, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)

     if(ln_timing) comm_time = comm_time + shympi_wtime() - time1




     !do i=1,nkndi
     !  do j=1,nkndi
     !    if(nnodesAssign(j) .eq. i) then
     !      test = 0
     !      k = 1
     !      do while(test .lt. j)
     !        test=nnodes(k)+test
     !        k = k+1
     !      end do
     !    end if
     !  end do
     !  univocalNodesAssign(i) = k-2
    !end do

     numberNodes = nnodes
     !totalNodes = nkndi

     deallocate(nnodes)
     !deallocate(nnodesAssign)

     if(my_id .eq. 0) then
       allocate(scatterNodes(totalNodes))
       k=0
       !r=0
       do i=1,n_threads
         n=0
         do j=1,nkndi
           if(allNodesAssign(j,i) .eq. 1) then
             if(univocalNodesAssign(j) .eq. (i-1)) then
               k=k+1
               !r=r+1
               scatterNodes(k) = j
               !reorderNodes(r) = j
             else
               n=n+1
               temp(n) = j
             end if
           end if
         end do
         do j=1,n
           k=k+1
           scatterNodes(k)=temp(j)
         end do
       end do


       !open(unit=1500, file="reorder.txt", action='write')
       !   write(6,*)'UnivocalNodes and ReorderNodes:'
       !   do i=1,nkndi
       !     write(1500,*),univocalNodesAssign(i),i,reorderNodes(i) 
       !   end do
       !close(1500)

     end if


   end subroutine

!########################################################################################################!

  
  subroutine makeNodesAssign(fullNodesAssign, nkn)

    use mpi_common_struct

    implicit none

    integer, dimension(totalNodes) :: fullNodesAssign
    integer i, j, node, nkn
    integer, dimension(n_threads) :: offset

    allNodesAssign = -1

    offset(1) = 0
    do i=1, n_threads
      if(i .gt. 1) then
        offset(i) = offset(i-1) + numberNodes(i-1)
      end if
      do j=1, numberNodes(i)
        node = fullNodesAssign(offset(i)+j)
        allNodesAssign(node,i) = 1 
      end do
    end do
   
    return

  end subroutine


!-----------------------------------------------------------------!

   subroutine tvd_mpi_init(itvd)

     use tvd

     implicit none

     integer itvd

     include 'param.h'

     itvd_type = itvd

     if( itvd_type .eq. 2 ) then
       !call tvd_upwind_init
       call tvd_mpi_upwind_init
       call setup_upwind_mpi
     end if


     if( itvd .eq. 0 ) then
       write(6,*) 'no horizontal TVD scheme used'
     else
       write(6,*) 'horizontal TVD scheme initialized: ',itvd
     end if

   end subroutine

!-----------------------------------------------------------------!

   subroutine setup_upwind_mpi

     use tvd
     !use mpi_communication

     implicit none

     integer ie,k,ii,j,ienew,n,rank,nrank,ierr,h,length
     !type (tvdInfo) :: tvdmpi
     type (tvdInfo) :: temp,temp2
     !logical :: rrecv(n_threads),rsend(n_threads)
     integer :: rrecv(n_threads)
     integer :: rsend(n_threads)
     integer tag,newtag,nsend,nrecv
     integer recv_data,send_data
     integer :: req(n_threads)
     integer st(MPI_STATUS_SIZE)
     double precision time1

     allocate(tvdmpi%eleID(3,3,domain%elems%domainID))
     allocate(tvdmpi%myele(3,3,domain%elems%domainID))

     allocate(temp%relem(n_threads))
     allocate(temp2%relem(n_threads))
     do rank=0,n_threads-1
       allocate(temp%relem(rank+1)%items(domain%elems%domainID))
       allocate(temp2%relem(rank+1)%items(domain%elems%domainID))
       temp%relem(rank+1)%numItems=0
       temp%relem(rank+1)%items=0
       temp2%relem(rank+1)%numItems=0
       temp2%relem(rank+1)%items=0
     end do


     n=0
     temp%recv=0
     do ie=1,nel
       !if(allPartAssign(ie).eq.my_id) then      
       n=n+1
       do ii=1,3
         do j=1,3
           ienew=ietvdup(j,ii,ie)
           if(ienew.ne.0) then
               rank=allPartAssign(ienew)
               if(rank.eq.my_id) then
                 do k=1,domain%elems%domainID
                   if(ienew.eq.domain%elems%globalID(k)) then
                     tvdmpi%eleID(j,ii,n)=k
                     !ietvdup(j,ii,ie)=k       !gmica valutare se decommentare qui
                   end if
                 end do
               else
                 temp%relem(rank+1)%numItems= &
                 temp%relem(rank+1)%numItems +1
                 temp%relem(rank+1)%items(temp%relem &
                 (rank+1)%numItems)=ienew
                 tvdmpi%eleID(j,ii,n)=ienew
               end if
               tvdmpi%myele(j,ii,n)=rank
             end if
         end do
       end do
       !end if
     end do

     k=0
     do rank=0,n_threads-1
       if(temp%relem(rank+1)%numItems.gt.0) k=k+1
     end do

     allocate(tvdmpi%relem(k))
     k=0
     do rank=0,n_threads-1
       if(temp%relem(rank+1)%numItems.gt.0) then
         k=k+1
         call remove_dups(temp%relem(rank+1)%items, &
         length,temp2%relem(rank+1)%items)
         if(.not. allocated(tvdmpi%relem(k)%items)) then
           allocate(tvdmpi%relem(k)%items(length))
         end if
         do j=1,length
           tvdmpi%relem(k)%items(j)=temp2%relem(rank+1)%items(j)
         end do
         tvdmpi%relem(k)%numItems=length
         tvdmpi%relem(k)%prank=rank
       end if
     end do


!     do k=1,tvdmpi%recv
!       allocate(tvdmpi%relem(k)%x(tvdmpi%relem(k)%numItems))
!       allocate(tvdmpi%relem(k)%y(tvdmpi%relem(k)%numItems))
!     end do
!
!     do ie=1,nel
!       do ii=1,3
!         do j=1,3
!           ienew=ietvdup(j,ii,ie)
!           if(ienew.ne.0) then
!             rank=allPartAssign(ienew)
!             if(rank.ne.my_id) then
!               do k=1,tvdmpi%recv
!loop1:           do h=1,tvdmpi%relem(k)%numItems
!                   if(tvdmpi%relem(k)%items(h).eq.ienew) then
!                     tvdmpi%relem(k)%x(h)=tvdupx(j,ii,ie)
!                     tvdmpi%relem(k)%y(h)=tvdupy(j,ii,ie)
!                     exit loop1
!                   end if
!                 end do loop1
!               end do
!             end if
!           end if
!         end do
!       end do
!     end do
!
!    if(my_id.eq.0) write(6,*)'myx:',tvdmpi%relem(1)%prank,tvdmpi%relem(1)%x
!


     !call remove_dups_s(temp%relem,length,tvdmpi%relem,k)

     tvdmpi%recv=k

     nrank=1
     rrecv=0
     rsend=0
     do k=0,n_threads-1
       do j=1,tvdmpi%recv
         if(tvdmpi%relem(j)%prank.eq.k) then
           rsend(k+1) = tvdmpi%relem(j)%numItems
           exit
         else
           rsend(k+1) = 0
         end if
       end do
     end do

     if(ln_timing) time1 = shympi_wtime()

     call mpi_Alltoall(rsend,1,MPI_INTEGER,rrecv, &
     & 1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

     if(ln_timing) comm_time = comm_time + shympi_wtime() - time1
  
    nrecv=0
    do j=1,n_threads
      if(rrecv(j).ne.0) nrecv= nrecv+1
    end do
    tvdmpi%send=nrecv
    allocate(tvdmpi%selem(nrecv))
    nrecv=0
    do j=1,n_threads
      if(rrecv(j).ne.0) then
        nrecv= nrecv+1
        tvdmpi%selem(nrecv)%numItems=rrecv(j)
        tvdmpi%selem(nrecv)%prank=j-1
        allocate(tvdmpi%selem(nrecv)%items(rrecv(j)))
      end if
    end do

    if(ln_timing) time1 = shympi_wtime()
 
    do j=1,tvdmpi%send
      rank=tvdmpi%selem(j)%prank
      nrecv=tvdmpi%selem(j)%numItems
        call MPI_IRecv(tvdmpi%selem(j)%items,nrecv,MPI_INTEGER, rank, MPI_ANY_TAG, &
                      & MPI_COMM_WORLD,req(rank+1),ierr)
    end do

    
    do j=1,tvdmpi%recv
      rank=tvdmpi%relem(j)%prank
      nsend=tvdmpi%relem(j)%numItems
      call MPI_Send(tvdmpi%relem(j)%items,nsend,MPI_INTEGER,rank,rank, &
        &                  MPI_COMM_WORLD,ierr)
    end do

     do j=1,tvdmpi%send
       rank=tvdmpi%selem(j)%prank
       CALL MPI_Wait(req(rank+1), st, ierr)
     end do

     if(ln_timing) comm_time = comm_time + shympi_wtime() - time1

     do j=1,tvdmpi%recv
       deallocate(temp2%relem(j)%items)
     end do

     do j=1,n_threads
       deallocate(temp%relem(j)%items)
     end do
     deallocate(temp%relem)
     deallocate(temp2%relem)
     
     return

   end subroutine

!-----------------------------------------------------------------!

	subroutine find_elem_from_old_mpi(ieold,xp,yp,ielem)

! finds element for point (xp,yp) starting from ieold

        use mpi_communication_struct
	use basin, only : neldi
        use regular

	implicit none

	integer ieold
	double precision xp,yp
	integer ielem	!element number on return

	integer iem,iep

!-------------------------------------------------------------
! check if old element is given -> if not test all elements
!-------------------------------------------------------------

	if( ieold .le. 0 .or. ieold .gt. neldi ) then
	  call find_element(xp,yp,ielem)
	  return
	end if

!-------------------------------------------------------------
! check if in old element
!-------------------------------------------------------------

	if( in_element(ieold,xp,yp) ) then
	  ielem = ieold
	  return
	end if

!-------------------------------------------------------------
! start from old element going upwards and downwards
!-------------------------------------------------------------

	iem = ieold-1
	if( iem .lt. 1 ) iem = neldi		!BUG_27.01.2011
	iep = ieold+1
	if( iep .gt. neldi ) iep = 1		!BUG_27.01.2011

	do while( iem .ne. ieold .and. iep .ne. ieold )
	  if( in_element(iem,xp,yp) ) then
	    ielem = iem
	    return
	  end if
	  iem = iem - 1
	  if( iem .lt. 1 ) iem = neldi

	  if( in_element(iep,xp,yp) ) then
	    ielem = iep
	    return
	  end if
	  iep = iep + 1
	  if( iep .gt. neldi ) iep = 1
	end do

!-------------------------------------------------------------
! no element found
!-------------------------------------------------------------

	ielem = 0

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------

    end subroutine

   subroutine tvd_mpi_upwind_init

! initializes position of upwind node
!
! sets position and element of upwind node

     use tvd
     use basin
     use evgeom

     implicit none

     include 'param.h'


     logical bsphe
     integer inode
     integer isphe
     integer ie,ii,j,k
     integer ienew,ienew2
     double precision x,y
     double precision r

     double precision xc,yc,xd,yd,xu,yu
     double precision dlat0,dlon0                    !center of projection

     integer ieext,gie

     write(6,*) 'setting up tvd mpi upwind information...'

     call get_coords_ev(isphe)
     bsphe = isphe .eq. 1
     inode = 0

     do ie=1,nel

       gie=domain%elems%globalID(ie)

       if ( bsphe ) call ev_make_center(gie,dlon0,dlat0)

       do ii=1,3

         k = nen3v(ii,gie)
         xc = xgv(k)
         yc = ygv(k)
         if ( bsphe ) call ev_g2c(xc,yc,xc,yc,dlon0,dlat0)

         j = mod(ii,3) + 1
         k = nen3v(j,gie)
         xd = xgv(k)
         yd = ygv(k)
         if ( bsphe ) call ev_g2c(xd,yd,xd,yd,dlon0,dlat0)

         xu = 2*xc - xd
         yu = 2*yc - yd
         if ( bsphe ) call ev_c2g(xu,yu,xu,yu,dlon0,dlat0)
         x = xu
         y = yu

         call find_elem_from_old_mpi(gie,x,y,ienew)

         tvdupx(j,ii,ie) = x
         tvdupy(j,ii,ie) = y
         ietvdup(j,ii,ie) = ienew

         j = mod(ii+1,3) + 1
         k = nen3v(j,gie)
         xd = xgv(k)
         yd = ygv(k)
	 if ( bsphe ) call ev_g2c(xd,yd,xd,yd,dlon0,dlat0)

         xu = 2*xc - xd
	 yu = 2*yc - yd
	 if ( bsphe ) call ev_c2g(xu,yu,xu,yu,dlon0,dlat0)
	 x = xu
	 y = yu

	 call find_elem_from_old_mpi(gie,x,y,ienew)

         tvdupx(j,ii,ie) = x
         tvdupy(j,ii,ie) = y
         ietvdup(j,ii,ie) = ienew

         tvdupx(ii,ii,ie) = 0.
         tvdupy(ii,ii,ie) = 0.
         ietvdup(ii,ii,ie) = 0

       end do
     end do

     write(6,*) '...tvd mpi upwind setup done (itvd=2)'

   end subroutine

!-----------------------------------------------------------------!
!
!   subroutine find_elem_mpi(ieold,xp,yp,ielem)
!
!! finds element for point (xp,yp) starting from ieold
!
!     use basin, only : nkn,nel,ngr,mbw
!
!     implicit none
!
!     integer ieold
!     double precision xp,yp
!     integer ielem	!element number on return
!     
!     logical in_element
!     integer iem,iep
!     
!!-------------------------------------------------------------
!! check if old element is given -> if not test all elements
!!-------------------------------------------------------------
!
!     if( ieold .le. 0 .or. ieold .gt. nel ) then
!       call find_element(xp,yp,ielem)
!       return
!     end if
!     
!!-------------------------------------------------------------
!! check if in old element
!!-------------------------------------------------------------
!
!     if( in_element(ieold,xp,yp) ) then
!       ielem = ieold
!       return
!     end if
!     
!!-------------------------------------------------------------
!! start from old element going upwards and downwards
!!-------------------------------------------------------------
!
!     iem = ieold-1
!     if( iem .lt. 1 ) iem = nel		!BUG_27.01.2011
!     iep = ieold+1
!     if( iep .gt. nel ) iep = 1		!BUG_27.01.2011
!     
!     do while( iem .ne. ieold .and. iep .ne. ieold )
!       if( in_element(iem,xp,yp) ) then
!         ielem = iem
!         return
!       end if
!       iem = iem - 1
!       if( iem .lt. 1 ) iem = nel
!     
!       if( in_element(iep,xp,yp) ) then
!         ielem = iep
!         return
!       end if
!       iep = iep + 1
!       if( iep .gt. nel ) iep = 1
!     end do
!
!!-------------------------------------------------------------
!! no element found
!!-------------------------------------------------------------
!
!     ielem = 0
!
!!-------------------------------------------------------------
!! end of routine
!!-------------------------------------------------------------
!
!     end subroutine
!
!!-----------------------------------------------------------------!

  subroutine makeMyNen3v(elems,nodes, mystruct, struct)

    implicit none
  
    integer i,j,h,k
    integer mystruct(3,1), struct(3,1)
    type (MAPPING_INFO) :: elems,nodes,temp
    
    do i=1,elems%totalID
      do j=1,3
        k = mystruct(j,i)
        do h=1,nodes%totalID
          if(nodes%globalID(h) .eq. k ) then
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

    use basin

    implicit none
    ! input arguments
    integer mysize, mode
!    integer struct(mysize,1)
    type (MAPPING_INFO) :: infoStruct
    integer, dimension(mysize,neldi) :: struct
    integer, dimension(mysize,infoStruct%totalID) :: mystruct
    integer :: st(MPI_STATUS_SIZE), ierr

    ! output arguments
!    type (MAKE_STRUCT) :: mystruct,struct

    ! local 
    integer n, j

    do n=1, infoStruct%domainID
      do j=1, mysize
        mystruct(j,n) = struct(j,infoStruct%globalID(n))
      end do
    end do


    return

  end subroutine


  subroutine transfer_domain

    !use basin 

    implicit none

    !write(6,*)'transfer_domain',elems%domainID,domain%nodes%domainID,neldi,nkndi,my_id

    integer, allocatable :: nen3v_new(:,:)
    integer, allocatable :: ipev_new(:)
    integer, allocatable :: ipv_new(:)
    integer, allocatable :: iarv_new(:)
    integer, allocatable :: iarnv_new(:)

    double precision, allocatable :: xgv_new(:)
    double precision, allocatable :: ygv_new(:)
    double precision, allocatable :: hm3v_new(:,:)

    integer own_nodes(nkndi)

    integer i,ii,ie,k,ierr

!   ----------------------------------
!   allocate aux arrays
!   ----------------------------------


    allocate(nen3v_new(3,domain%elems%totalID))
    allocate(ipev_new(nel))
    allocate(ipv_new(domain%nodes%totalID))
    allocate(iarv_new(nel))
    allocate(iarnv_new(nkn))

    allocate(xgv_new(domain%nodes%totalID))
    allocate(ygv_new(domain%nodes%totalID))
#ifdef DEBUGON
    allocate(hm3v_new(3,nel_local))
#else
    allocate(hm3v_new(3,nel))
#endif


    own_nodes = 0
    do i=1,nkn
      k = domain%nodes%globalID(i)
      own_nodes(k) = i
      ipv_new(i) = ipv(k)
      iarnv_new(i) = iarnv(k)
      xgv_new(i) = xgv(k)
      ygv_new(i) = ygv(k)
    end do
    do i=nkn+1,domain%nodes%totalID
      k = domain%nodes%globalID(i)
      own_nodes(k) = i
      ipv_new(i) = ipv(k)
      xgv_new(i) = xgv(k)
      ygv_new(i) = ygv(k)
    end do

    do i=1,domain%elems%totalID
      ie = domain%elems%globalID(i)
#ifdef DEBUGON
      hm3v_new(:,i) = hm3v(:,ie)
#endif
      do ii=1,3
        k = nen3v(ii,ie)
        if( own_nodes(k) <= 0 ) then
          write(6,*)'error stop transfer_domain: internal error'
          write(6,*) ie,k,own_nodes(k),my_id
          stop
        end if
        nen3v_new(ii,i) = own_nodes(k)
      end do
    end do

    do i=1,nel
      ie = domain%elems%globalID(i)
#ifndef DEBUGON
      hm3v_new(:,i) = hm3v(:,ie)
#endif
      ipev_new(i) = ipev(ie)
      iarv_new(i) = iarv(ie)
    end do

    call basin_init(nkn,nel,nkn_local,nel_local)

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

!******************************************************************

	subroutine shympi_setup

        use basin
        use geom
        use mpi_communication_struct
        use tvd
        use links
        use para

	implicit none

        integer i,nlkdi,ierr
        integer, dimension(:,:), allocatable :: mystruct
        integer temp_nkn,temp_nel,itvd
	double precision xmin,xmax,ymin,ymax

        temp_nel = 0
        do i=1,neldi
          if(allPartAssign(i) .eq. my_id) then
            temp_nel = temp_nel + 1
          end if
        end do
        nel=neldi

        nlkdi = 3*neldi+2*nkndi

        allocate(ilinkv(nkndi+1))
        allocate(lenkv(nlkdi))
        allocate(linkv(nlkdi))
        allocate(ieltv(3,neldi))
        allocate(total_ieltv(3,neldi))
        allocate(tilinkv(nkndi+1))
        allocate(tlinkv(nlkdi))

        call mklenk(nlkdi,nkndi,neldi,nen3v,ilinkv,lenkv)

        call mklink(nkndi,ilinkv,lenkv,linkv)

        call mkielt(nkndi,neldi,ilinkv,lenkv,linkv,ieltv)

        tilinkv = ilinkv
        tlinkv = linkv
        total_ieltv = ieltv

        nel =temp_nel

        call makeMyEle(nel,neldi,nen3v)

        allocate(mystruct(3,domain%elems%totalID))

        call makeMyStruct(domain%elems, nen3v, mystruct, 3)

        call makeMyNodes(mystruct, nen3v, nkndi)
     
        !nel=neldi 
        call mod_tvd_init(nel) 
        itvd = nint(getpar("itvd"))
        call tvd_mpi_init(itvd)

        nkn =domain%nodes%domainID
        nel =domain%elems%domainID
 
        if(nel .eq. 0) then
          ngr = 0
        end if

        nel_local=domain%elems%totalID
        nkn_local=domain%nodes%totalID

        temp_nkn=nkndi
        temp_nel=neldi
        call transfer_domain
        neldi=temp_nel
        nkndi=temp_nkn


        deallocate(ilinkv)
        deallocate(linkv)
        deallocate(lenkv)
        deallocate(ieltv)
       
	call bas_get_minmax(xmin,ymin,xmax,ymax)

	domain%nodes%xmin=xmin
	domain%nodes%xmax=xmax
	domain%nodes%ymin=ymin
	domain%nodes%ymax=ymax

	write(6,*)'domain minmax:',xmin,xmax,ymin,ymax
 
	return

	end subroutine shympi_setup

!#############################################################################!
!***********************  start countLevel subroutine  ***********************!
!#############################################################################!
!#############################################################################!
!* This subroutine count the number of levels of each process and the total  *!
!* level for each process                                                    *!
!#############################################################################!

  subroutine countLevel(numlevels, totalLevels)

    use mpi_common_struct
    use levels

    implicit none

    include 'param.h'

    integer i
    integer :: st(MPI_STATUS_SIZE), ierr

    integer numlevels, totalLevels
    double precision time1 

    numlevels = 0
    do i=1,domain%elems%domainID
      numlevels = numlevels + ilhv(i)
    end do

    allocate(numberLevels(n_threads))

    if(ln_timing) time1 = shympi_wtime()

    call MPI_ALLGATHER(numlevels,1,MPI_INTEGER,numberLevels,1,MPI_INTEGER, MPI_COMM_WORLD,ierr)

    if(ln_timing) comm_time = comm_time + shympi_wtime() - time1

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

  subroutine recountLevel(elems, ilhv)

    use mpi_common_struct

    implicit none

    include 'param.h'

    integer i
    integer, dimension(:) :: ilhv
    type (MAPPING_INFO) :: elems
    integer :: st(MPI_STATUS_SIZE), ierr

    integer numlevels    

    numlevels = 0
    do i=1,elems%domainID
      numlevels = numlevels + ilhv(i)
    end do

    if ( numlevels .ge. (1.20 * numberLevels(my_id+1)) ) then
      write(6,*)'my old levels =', numberLevels(my_id+1), my_id
      write(6,*)'my new levels =', numlevels, my_id
    end if

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

  subroutine remove_dups0(array,length,res,dimout)

    implicit none
    
    integer, dimension(:) :: array         ! The input
    !integer :: res(size(array))  ! The output
    integer :: length            ! The number of unique elements in output
    integer :: i, j, k, h, temp
    integer, optional :: dimout
    integer,allocatable,dimension(:) :: res

    if (present(dimout) .and. (.not. allocated(res))) then
      allocate(res(dimout))
    else if (.not. allocated(res)) then
       write(6,*)'error in remove_dups0'
    end if

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
         if(h .eq. 2) exit
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

  subroutine remove_dups(array,length,res,dimout)

    implicit none
    
    integer, dimension(:) :: array         ! The input
    !integer :: res(size(array))  ! The output
    integer :: length            ! The number of unique elements in output
    integer :: i, j, k, h, temp
    integer, optional :: dimout
    integer,allocatable,dimension(:) :: res

    if (present(dimout) .and. (.not. allocated(res))) then
      allocate(res(dimout))
    else if (.not. allocated(res)) then
       write(6,*)'error in remove_dups'
    end if

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
         if(h .eq. 2) exit
         h = h-1
       end do

    end do outer

    length = k

    return

  end subroutine remove_dups


  subroutine remove_dups_s(array,length,res,dimout)

    implicit none
    
    type (commItem), dimension(:) :: array         ! The input
    !integer :: res(size(array))  ! The output
    integer :: length            ! The number of unique elements in output
    integer :: i, j, k, h, temp
    integer, optional :: dimout
    type (commItem),allocatable,dimension(:) :: res

    if (present(dimout) .and. (.not. allocated(res))) then
      allocate(res(dimout))
    else if (.not. allocated(res)) then
       write(6,*)'error in remove_dups'
    end if

    if( size(array) .le. 0) then
      length = 0
      return
    end if

    do i=1,size(res)
    !do i=1,size(array)
      res(i)%numItems = 0
    end do

    k = 0
    !res(1) = array(1)
    outer: do i=1,size(array)
       if(array(i)%numItems .gt. 0) then
         k = k + 1
         res(k) = array(i)
         res(k)%prank=i-1
       end if
    end do outer

    length = k

    return

  end subroutine remove_dups_s


!######################################################################!
!********************  end Module Global_Graph_Data  ******************!
!######################################################################!

  end module mpi_graph_elem
