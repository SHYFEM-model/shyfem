module mpi_communication

   use mpi_communication_struct
   use mpi_utility
   use shympi
   use mpi_common_struct

   interface exchange_struct_2D
      module procedure exchange_struct_2D_i  &
                       ,exchange_struct_2D_r &
                       ,exchange_struct_2D_d
   end interface exchange_struct_2D

   interface rebuild_struct_2D_sum
      module procedure rebuild_struct_2D_sum_i  &
                       ,rebuild_struct_2D_sum_r &
                       ,rebuild_struct_2D_sum_d
   end interface rebuild_struct_2D_sum

   interface exchange_struct_3D
      module procedure exchange_struct_3D_r &
                       ,exchange_struct_3D_d
   end interface exchange_struct_3D

   interface rebuild_struct_3D_sum
      module procedure rebuild_struct_3D_sum_r &
                       ,rebuild_struct_3D_sum_d
   end interface rebuild_struct_3D_sum



   contains

     subroutine spread_nodes_to_neighbor(srequests, rrequests, data_receive, maxNodes)

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
        if(.not. allocated(data_receive)) then
          allocate(data_receive(maxNodes,comm_loop))
        else
          deallocate(data_receive)
          allocate(data_receive(maxNodes,comm_loop))
        end if

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
                mypart%mysend%process(i),tag, MPI_COMM_WORLD,srequests(i),ierr)
           if(ierr .ne. MPI_SUCCESS) then
              write(6,*)'ierror:',ierr
           end if

           tag = newtag + my_id + mypart%mysend%process(i) 
           ! Non-blocking receives
           call MPI_Irecv(data_receive(1,i),receiveNodes,MPI_INTEGER, &
                mypart%myreceive%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
           if(ierr .ne. MPI_SUCCESS) then
              write(6,*)'ierror:',ierr
           end if

        end do

        return

     end subroutine spread_nodes_to_neighbor

!##########################################################################!
!**********************  start send_recv_real_data  ***********************!
!##########################################################################!
!##########################################################################!
!* This subroutine is usefull for send or/and receve real data to/from  ***!
!* others processes                                                       *!
!##########################################################################!

     subroutine exchange_struct_2D_r(struct, mode, srequests, rrequests, data_send, data_receive)
        
        use basin

        implicit none

        ! arguments 
        integer mode
        real struct(nkn)

        !common
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

        ! local
        integer i,j,sendElements,receiveElements, comm_loop,k, tag, newtag
        real, dimension(:,:),allocatable :: data_send, data_receive 

        ! MPI
        integer ierr, data_type
        integer, dimension(:) :: srequests,rrequests

        comm_loop = mypart%mysend%sends
        data_type = MPI_REAL

        newtag = next_mpi_tag()

        allocate(data_send(mypart%mysend%maxItems*100,comm_loop))
        allocate(data_receive(mypart%myreceive%maxItems*100,comm_loop))

        if(mode .eq. 1) then

           do i = 1, comm_loop

              sendElements = mypart%mysend%all_send(i)%numItems
              receiveElements = mypart%myreceive%all_receive(i)%numItems

              do k=1,sendElements
                 data_send(k,i) = struct(mypart%mysend%all_send(i)%items(k))
              end do

!              if(my_id .eq. 2 .and. mypart%mysend%process(i) .eq. 1 ) then
!                 open(unit=1900, file="2data_send1.txt", action='write')
!                 do k=1,sendElements
!                    write(1900,*),k,mypart%mysend%all_send(i)%items(k),myele%globalID(mypart%mysend%all_send(i)%items(k))
!                    do j=1,size1                        
!                       write(1900,*), data_send(i,j,k)
!                    end do
!                 end do
!                 close(1900)
!                 stop
!              end if


              srequests(i) = MPI_REQUEST_NULL
              rrequests(i) = MPI_REQUEST_NULL

              tag = newtag + mypart%mysend%process(i) + my_id 

              ! Non-blocking sends
              call MPI_Isend(data_send(1,i),sendElements,data_type, &
                mypart%mysend%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              tag = newtag + my_id + mypart%mysend%process(i) 
              ! Non-blocking receives
              call MPI_Irecv(data_receive(1,i),receiveElements,data_type, &
                mypart%myreceive%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

           end do

        else if(mode .eq. 2) then

           do i = 1, comm_loop

              sendElements = mypart%mysend%part_send(i)%numItems
              receiveElements = mypart%myreceive%part_receive(i)%numItems

              do k=1,sendElements
                 data_send(k,i) = struct(mypart%mysend%part_send(i)%items(k))
              end do

              srequests(i) = MPI_REQUEST_NULL
              rrequests(i) = MPI_REQUEST_NULL

              tag = newtag + mypart%mysend%process(i) + my_id 

              ! Non-blocking sends
              call MPI_Isend(data_send(1,i),sendElements,data_type, &
                mypart%mysend%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              tag = newtag + my_id + mypart%mysend%process(i) 
              ! Non-blocking receives
              call MPI_Irecv(data_receive(1,i),receiveElements,data_type, &
                mypart%myreceive%process(i),tag, MPI_COMM_WORLD,rrequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

           end do

        else if(mode .eq. 3) then

           do i = 1, comm_loop

              sendElements = mypart%mysend%node_send(i)%numItems
              receiveElements = mypart%mysend%node_send(i)%numItems 

              do k=1,sendElements
                 data_send(k,i) = struct(mypart%mysend%node_send(i)%items(k))
              end do


              srequests(i) = MPI_REQUEST_NULL
              rrequests(i) = MPI_REQUEST_NULL

              tag = newtag + mypart%mysend%process(i) + my_id 

              ! Non-blocking sends
              call MPI_Isend(data_send(1,i),sendElements,data_type, &
                mypart%mysend%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              tag = newtag + my_id + mypart%mysend%process(i) 
              ! Non-blocking receives
              call MPI_Irecv(data_receive(1,i),receiveElements,data_type, &
                mypart%myreceive%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

           end do

        else
           write(6,*)'error in send_recv_int_data', mode
           stop
        end if

        return


    end subroutine

!##########################################################################!
!**********************  start send_recv_dp_data  ***********************!
!##########################################################################!
!##########################################################################!
!* This subroutine is usefull for send or/and receve double precision data !
!* to/from others processes                                               *!
!##########################################################################!

     subroutine exchange_struct_2D_d(struct, mode, srequests, rrequests, data_send, data_receive)
        
        use basin

        implicit none

        ! arguments 
        integer mode
        double precision struct(nkn)

        !common
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

        ! local
        integer i,j,sendElements,receiveElements, comm_loop,k, tag, newtag
        double precision, dimension(:,:),allocatable :: data_send, data_receive 

        ! MPI
        integer ierr, data_type
        integer, dimension(:) :: srequests,rrequests

        comm_loop = mypart%mysend%sends
        data_type = MPI_DOUBLE_PRECISION

        newtag = next_mpi_tag()

        allocate(data_send(mypart%mysend%maxItems*100,comm_loop))
        allocate(data_receive(mypart%myreceive%maxItems*100,comm_loop))

        if(mode .eq. 1) then

           do i = 1, comm_loop

              sendElements = mypart%mysend%all_send(i)%numItems
              receiveElements = mypart%myreceive%all_receive(i)%numItems

              do k=1,sendElements
                 data_send(k,i) = struct(mypart%mysend%all_send(i)%items(k))
              end do

!              if(my_id .eq. 2 .and. mypart%mysend%process(i) .eq. 1 ) then
!                 open(unit=1900, file="2data_send1.txt", action='write')
!                 do k=1,sendElements
!                    write(1900,*),k,mypart%mysend%all_send(i)%items(k),myele%globalID(mypart%mysend%all_send(i)%items(k))
!                    do j=1,size1                        
!                       write(1900,*), data_send(i,j,k)
!                    end do
!                 end do
!                 close(1900)
!                 stop
!              end if


              srequests(i) = MPI_REQUEST_NULL
              rrequests(i) = MPI_REQUEST_NULL

              tag = newtag + mypart%mysend%process(i) + my_id 

              ! Non-blocking sends
              call MPI_Isend(data_send(1,i),sendElements,data_type, &
                mypart%mysend%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              tag = newtag + my_id + mypart%mysend%process(i) 
              ! Non-blocking receives
              call MPI_Irecv(data_receive(1,i),receiveElements,data_type, &
                mypart%myreceive%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

           end do

        else if(mode .eq. 2) then

           do i = 1, comm_loop

              sendElements = mypart%mysend%part_send(i)%numItems
              receiveElements = mypart%myreceive%part_receive(i)%numItems

              do k=1,sendElements
                 data_send(k,i) = struct(mypart%mysend%part_send(i)%items(k))
              end do

              srequests(i) = MPI_REQUEST_NULL
              rrequests(i) = MPI_REQUEST_NULL

              tag = newtag + mypart%mysend%process(i) + my_id 

              ! Non-blocking sends
              call MPI_Isend(data_send(1,i),sendElements,data_type, &
                mypart%mysend%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              tag = newtag + my_id + mypart%mysend%process(i) 
              ! Non-blocking receives
              call MPI_Irecv(data_receive(1,i),receiveElements,data_type, &
                mypart%myreceive%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

           end do

        else if(mode .eq. 3) then

           do i = 1, comm_loop

              sendElements = mypart%mysend%node_send(i)%numItems
              receiveElements = mypart%mysend%node_send(i)%numItems 

              do k=1,sendElements
                 data_send(k,i) = struct(mypart%mysend%node_send(i)%items(k))
              end do


              srequests(i) = MPI_REQUEST_NULL
              rrequests(i) = MPI_REQUEST_NULL

              tag = newtag + mypart%mysend%process(i) + my_id 

              ! Non-blocking sends
              call MPI_Isend(data_send(1,i),sendElements,data_type, &
                mypart%mysend%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              tag = newtag + my_id + mypart%mysend%process(i) 
              ! Non-blocking receives
              call MPI_Irecv(data_receive(1,i),receiveElements,data_type, &
                mypart%myreceive%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

           end do

        else
           write(6,*)'error in send_recv_int_data', mode
           stop
        end if

        return


    end subroutine



!##########################################################################!
!**********************  start send_recv_int_data  ***********************!
!##########################################################################!
!##########################################################################!
!* This subroutine is usefull for send or/and receve integer data to/from *!
!* others processes                                                       *!
!##########################################################################!

     subroutine exchange_struct_2D_i(struct, mode, srequests, rrequests, data_send, data_receive)
        
        use basin

        implicit none

        ! arguments 
        integer mode
        integer struct(nkn)

        !common
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

        ! local
        integer i,j,sendElements,receiveElements, comm_loop,k, tag, newtag
        integer, dimension(:,:),allocatable :: data_send, data_receive 

        ! MPI
        integer ierr, data_type
        integer, dimension(:) :: srequests,rrequests

        comm_loop = mypart%mysend%sends
        data_type = MPI_INTEGER

        newtag = next_mpi_tag()

        allocate(data_send(mypart%mysend%maxItems*100,comm_loop))
        allocate(data_receive(mypart%myreceive%maxItems*100,comm_loop))

        if(mode .eq. 1) then

           do i = 1, comm_loop

              sendElements = mypart%mysend%all_send(i)%numItems
              receiveElements = mypart%myreceive%all_receive(i)%numItems

              do k=1,sendElements
                 data_send(k,i) = struct(mypart%mysend%all_send(i)%items(k))
              end do

!              if(my_id .eq. 2 .and. mypart%mysend%process(i) .eq. 1 ) then
!                 open(unit=1900, file="2data_send1.txt", action='write')
!                 do k=1,sendElements
!                    write(1900,*),k,mypart%mysend%all_send(i)%items(k),myele%globalID(mypart%mysend%all_send(i)%items(k))
!                    do j=1,size1                        
!                       write(1900,*), data_send(i,j,k)
!                    end do
!                 end do
!                 close(1900)
!                 stop
!              end if


              srequests(i) = MPI_REQUEST_NULL
              rrequests(i) = MPI_REQUEST_NULL

              tag = newtag + mypart%mysend%process(i) + my_id 

              ! Non-blocking sends
              call MPI_Isend(data_send(1,i),sendElements,data_type, &
                mypart%mysend%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              tag = newtag + my_id + mypart%mysend%process(i) 
              ! Non-blocking receives
              call MPI_Irecv(data_receive(1,i),receiveElements,data_type, &
                mypart%myreceive%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

           end do

        else if(mode .eq. 2) then

           do i = 1, comm_loop

              sendElements = mypart%mysend%part_send(i)%numItems
              receiveElements = mypart%myreceive%part_receive(i)%numItems

              do k=1,sendElements
                 data_send(k,i) = struct(mypart%mysend%part_send(i)%items(k))
              end do

              srequests(i) = MPI_REQUEST_NULL
              rrequests(i) = MPI_REQUEST_NULL

              tag = newtag + mypart%mysend%process(i) + my_id 

              ! Non-blocking sends
              call MPI_Isend(data_send(1,i),sendElements,data_type, &
                mypart%mysend%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              tag = newtag + my_id + mypart%mysend%process(i) 
              ! Non-blocking receives
              call MPI_Irecv(data_receive(1,i),receiveElements,data_type, &
                mypart%myreceive%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

           end do

        else if(mode .eq. 3) then

           do i = 1, comm_loop

              sendElements = mypart%mysend%node_send(i)%numItems
              receiveElements = mypart%mysend%node_send(i)%numItems 

              do k=1,sendElements
                 data_send(k,i) = struct(mypart%mysend%node_send(i)%items(k))
              end do


              srequests(i) = MPI_REQUEST_NULL
              rrequests(i) = MPI_REQUEST_NULL

              tag = newtag + mypart%mysend%process(i) + my_id 

              ! Non-blocking sends
              call MPI_Isend(data_send(1,i),sendElements,data_type, &
                mypart%mysend%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              tag = newtag + my_id + mypart%mysend%process(i) 
              ! Non-blocking receives
              call MPI_Irecv(data_receive(1,i),receiveElements,data_type, &
                mypart%myreceive%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

           end do

        else
           write(6,*)'error in send_recv_int_data', mode
           stop
        end if

        return


    end subroutine

!##########################################################################!
!**********************  start send_recv_real_data  ***********************!
!##########################################################################!
!##########################################################################!
!* This subroutine is usefull for send or/and receve real data to/from    *!
!* others processes                                                       *!
!##########################################################################!

     subroutine exchange_struct_3D_r(struct, size1, mode, srequests, rrequests, data_send, data_receive)
        
        implicit none

        ! arguments 
        integer size1, mode
        real struct(size1,1)

        !common
        !integer itanf,itend,idt,nits,niter,it
        !common /femtim/ itanf,itend,idt,nits,niter,it

        ! local
        integer i,j,sendElements,receiveElements, comm_loop,k, tag, newtag
        real, dimension(:,:,:),allocatable :: data_send, data_receive 

        ! MPI
        integer ierr, data_type
        integer, dimension(:) :: srequests,rrequests

        comm_loop = mypart%mysend%sends
        data_type = MPI_REAL

        newtag = next_mpi_tag()

        allocate(data_send(size1,mypart%mysend%maxItems,comm_loop))
        allocate(data_receive(size1,mypart%myreceive%maxItems,comm_loop))

        if(mode .eq. 1) then

           do i = 1, comm_loop

              sendElements = mypart%mysend%all_send(i)%numItems
              receiveElements = mypart%myreceive%all_receive(i)%numItems

              do k=1,sendElements
                 do j=1,size1
                    data_send(j,k,i) = struct(j,mypart%mysend%all_send(i)%items(k))
                 end do
              end do

              srequests(i) = MPI_REQUEST_NULL
              rrequests(i) = MPI_REQUEST_NULL

              tag = newtag + mypart%mysend%process(i) + my_id 

              ! Non-blocking sends
              call MPI_Isend(data_send(1,1,i),sendElements*size1,data_type, &
                mypart%mysend%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              tag = newtag + my_id + mypart%mysend%process(i) 
              ! Non-blocking receives
              call MPI_Irecv(data_receive(1,1,i),receiveElements*size1,data_type, &
                mypart%myreceive%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

           end do

        else if(mode .eq. 2) then

           do i = 1, comm_loop

              sendElements = mypart%mysend%part_send(i)%numItems
              receiveElements = mypart%myreceive%part_receive(i)%numItems

              do k=1,sendElements
                 do j=1,size1
                    data_send(j,k,i) = struct(j,mypart%mysend%part_send(i)%items(k))
                 end do
              end do

              srequests(i) = MPI_REQUEST_NULL
              rrequests(i) = MPI_REQUEST_NULL

              tag = newtag + mypart%mysend%process(i) + my_id 

              ! Non-blocking sends
              call MPI_Isend(data_send(1,1,i),sendElements*size1,data_type, &
                mypart%mysend%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              tag = newtag + my_id + mypart%mysend%process(i) 
              ! Non-blocking receives
              call MPI_Irecv(data_receive(1,1,i),receiveElements*size1,data_type, &
                mypart%myreceive%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

           end do

        else if(mode .eq. 3) then

           do i = 1, comm_loop

              sendElements = mypart%mysend%node_send(i)%numItems
              receiveElements = mypart%mysend%node_send(i)%numItems 

              do k=1,sendElements
                 do j=1,size1
                    data_send(j,k,i) = struct(j,mypart%mysend%node_send(i)%items(k))
                 end do
              end do


              srequests(i) = MPI_REQUEST_NULL
              rrequests(i) = MPI_REQUEST_NULL

              tag = newtag + mypart%mysend%process(i) + my_id 

              ! Non-blocking sends
              call MPI_Isend(data_send(1,1,i),sendElements*size1,data_type, &
                mypart%mysend%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              tag = newtag + my_id + mypart%mysend%process(i) 
              ! Non-blocking receives
              call MPI_Irecv(data_receive(1,1,i),receiveElements*size1,data_type, &
                mypart%myreceive%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

           end do

        else
           write(6,*)'error in send_recv_real_data', mode
           stop
        end if

        return


    end subroutine


!##########################################################################!
!**********************  start send_recv_real_data  ***********************!
!##########################################################################!
!##########################################################################!
!* This subroutine is usefull for send or/and receve real data to/from    *!
!* others processes                                                       *!
!##########################################################################!

     subroutine exchange_struct_3D_d(struct, size1, mode, srequests, rrequests, data_send, data_receive)
        
        implicit none

        ! arguments 
        integer size1, mode
        double precision struct(size1,1)

        !common
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

        ! local
        integer i,j,sendElements,receiveElements, comm_loop,k, tag, newtag
        double precision, dimension(:,:,:),allocatable :: data_send, data_receive 

        ! MPI
        integer ierr, data_type
        integer, dimension(:) :: srequests,rrequests

        comm_loop = mypart%mysend%sends
        data_type = MPI_DOUBLE_PRECISION

        newtag = next_mpi_tag()

        allocate(data_send(size1,mypart%mysend%maxItems,comm_loop))
        allocate(data_receive(size1,mypart%myreceive%maxItems,comm_loop))

        if(mode .eq. 1) then

           do i = 1, comm_loop

              sendElements = mypart%mysend%all_send(i)%numItems
              receiveElements = mypart%myreceive%all_receive(i)%numItems

              do k=1,sendElements
                 do j=1,size1
                    data_send(j,k,i) = struct(j,mypart%mysend%all_send(i)%items(k))
                 end do
              end do

              srequests(i) = MPI_REQUEST_NULL
              rrequests(i) = MPI_REQUEST_NULL

              tag = newtag + mypart%mysend%process(i) + my_id 

              ! Non-blocking sends
              call MPI_Isend(data_send(1,1,i),sendElements*size1,data_type, &
                mypart%mysend%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              tag = newtag + my_id + mypart%mysend%process(i) 
              ! Non-blocking receives
              call MPI_Irecv(data_receive(1,1,i),receiveElements*size1,data_type, &
                mypart%myreceive%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

           end do

        else if(mode .eq. 2) then

           do i = 1, comm_loop

              sendElements = mypart%mysend%part_send(i)%numItems
              receiveElements = mypart%myreceive%part_receive(i)%numItems

              do k=1,sendElements
                 do j=1,size1
                    data_send(j,k,i) = struct(j,mypart%mysend%part_send(i)%items(k))
                 end do
              end do

              srequests(i) = MPI_REQUEST_NULL
              rrequests(i) = MPI_REQUEST_NULL

              tag = newtag + mypart%mysend%process(i) + my_id 

              ! Non-blocking sends
              call MPI_Isend(data_send(1,1,i),sendElements*size1,data_type, &
                mypart%mysend%process(i),tag, MPI_COMM_WORLD,srequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              tag = newtag + my_id + mypart%mysend%process(i) 
              ! Non-blocking receives
              call MPI_Irecv(data_receive(1,1,i),receiveElements*size1,data_type, &
                mypart%myreceive%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

           end do

        else if(mode .eq. 3) then

           do i = 1, comm_loop

              sendElements = mypart%mysend%node_send(i)%numItems
              receiveElements = mypart%mysend%node_send(i)%numItems 

              do k=1,sendElements
                 do j=1,size1
                    data_send(j,k,i) = struct(j,mypart%mysend%node_send(i)%items(k))
                 end do
              end do


              srequests(i) = MPI_REQUEST_NULL
              rrequests(i) = MPI_REQUEST_NULL

              tag = newtag + mypart%mysend%process(i) + my_id 

              ! Non-blocking sends
              call MPI_Isend(data_send(1,1,i),sendElements*size1,data_type, &
                mypart%mysend%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              tag = newtag + my_id + mypart%mysend%process(i) 
              ! Non-blocking receives
              call MPI_Irecv(data_receive(1,1,i),receiveElements*size1,data_type, &
                mypart%myreceive%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

           end do

        else
           write(6,*)'error in send_recv_real_data', mode
           stop
        end if

        return


    end subroutine


!##########################################################################!
!**********************  start send_recv_dp_data  ***********************!
!##########################################################################!
!##########################################################################!
!* This subroutine is usefull for send or/and receve double precision data !
!* to/from others processes                                                !
!##########################################################################!

     subroutine send_recv_dp_data(struct, size1, mode, srequests, rrequests, data_send, data_receive)

        implicit none

        ! arguments 
        integer size1, mode
        double precision struct(size1,1)

        !common
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

        ! local
        integer i,j,sendElements,receiveElements, comm_loop,k, tag, newtag
        double precision, dimension(:,:,:),allocatable :: data_send, data_receive 

        ! MPI
        integer ierr, data_type
        integer, dimension(:) :: srequests,rrequests

        comm_loop = mypart%mysend%sends
        data_type = MPI_DOUBLE_PRECISION

        newtag = next_mpi_tag()

        allocate(data_send(size1,mypart%mysend%maxItems,comm_loop))
        allocate(data_receive(size1,mypart%myreceive%maxItems,comm_loop))


        if(mode .eq. 1) then

           do i = 1, comm_loop

              sendElements = mypart%mysend%all_send(i)%numItems
              receiveElements = mypart%myreceive%all_receive(i)%numItems

              do k=1,sendElements
                 do j=1,size1
                    data_send(j,k,i) = struct(j,mypart%mysend%all_send(i)%items(k))
                 end do
              end do

              srequests(i) = MPI_REQUEST_NULL
              rrequests(i) = MPI_REQUEST_NULL

              tag = newtag + mypart%mysend%process(i) + my_id 

              ! Non-blocking sends
              call MPI_Isend(data_send(1,1,i),sendElements*size1,data_type, &
                mypart%mysend%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              tag = newtag + my_id + mypart%mysend%process(i) 
              ! Non-blocking receives
              call MPI_Irecv(data_receive(1,1,i),receiveElements*size1,data_type, &
                mypart%myreceive%process(i),tag, MPI_COMM_WORLD,rrequests(i), ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

           end do

        else if(mode .eq. 2) then

           do i = 1, comm_loop

              sendElements = mypart%mysend%part_send(i)%numItems
              receiveElements = mypart%myreceive%part_receive(i)%numItems

              do k=1,sendElements
                 do j=1,size1
                    data_send(j,k,i) = struct(j,mypart%mysend%part_send(i)%items(k))
                 end do
              end do

              srequests(i) = MPI_REQUEST_NULL
              rrequests(i) = MPI_REQUEST_NULL

              tag = newtag + mypart%mysend%process(i) + my_id 

              ! Non-blocking sends
              call MPI_Isend(data_send(1,1,i),sendElements*size1,data_type, &
                 mypart%mysend%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              tag = newtag + my_id + mypart%mysend%process(i) 
              ! Non-blocking receives
              call MPI_Irecv(data_receive(1,1,i),receiveElements*size1,data_type, &
                mypart%myreceive%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

           end do

        else if(mode .eq. 3) then

           do i = 1, comm_loop

              sendElements = mypart%mysend%node_send(i)%numItems
              receiveElements = mypart%mysend%node_send(i)%numItems 

              do k=1,sendElements
                 do j=1,size1
                    data_send(j,k,i) = struct(j,mypart%mysend%node_send(i)%items(k))
                 end do
              end do


              srequests(i) = MPI_REQUEST_NULL
              rrequests(i) = MPI_REQUEST_NULL

              tag = newtag + mypart%mysend%process(i) + my_id 

              ! Non-blocking sends
              call MPI_Isend(data_send(1,1,i), sendElements*size1,data_type, &
                mypart%mysend%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              tag = newtag + my_id + mypart%mysend%process(i) 
              ! Non-blocking receives
              call MPI_Irecv(data_receive(1,1,i),receiveElements*size1,data_type, &
                mypart%myreceive%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

           end do

        else
           write(6,*)'error in send_recv_dp_data', mode
           stop
        end if

        return


    end subroutine


!##########################################################################!
!**********************  start send_recv_real3D_data  *********************!
!##########################################################################!
!##########################################################################!
!* This subroutine is usefull for send or/and receve real data to/from    *!
!* others processes                                                       *!
!##########################################################################!

     subroutine send_recv_real3D_data(struct,size1,size2,mode,srequests,rrequests,data_send,data_receive)

        implicit none

        ! arguments 
        integer size1,size2,mode
        real struct(size2,size1,1)

        !common
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

        ! local
        integer i,j,sendElements,receiveElements, comm_loop,k,n, tag, newtag
        real, dimension(:,:,:,:),allocatable :: data_send, data_receive 

        ! MPI
        integer ierr, data_type
        integer, dimension(:) :: srequests,rrequests


        comm_loop = mypart%mysend%sends
        data_type = MPI_REAL

        newtag = next_mpi_tag()

        allocate(data_send(size2,size1,mypart%mysend%maxItems,comm_loop))
        allocate(data_receive(size2,size1,mypart%myreceive%maxItems,comm_loop))


        if(mode .eq. 1) then

           do i = 1, comm_loop

              sendElements = mypart%mysend%all_send(i)%numItems
              receiveElements = mypart%myreceive%all_receive(i)%numItems

              do k=1,sendElements
                 do j=1,size1
                    do n=1,size2
                       data_send(n,j,k,i) = struct(n,j,mypart%mysend%all_send(i)%items(k))
                    end do
                 end do
              end do

              srequests(i) = MPI_REQUEST_NULL
              rrequests(i) = MPI_REQUEST_NULL

              tag = newtag + mypart%mysend%process(i) + my_id 

              ! Non-blocking sends
              call MPI_Isend(data_send(1,1,1,i),sendElements*size1*size2,data_type, &
                mypart%mysend%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              tag = newtag + my_id + mypart%mysend%process(i) 
              ! Non-blocking receives
              call MPI_Irecv(data_receive(1,1,1,i),receiveElements*size1*size2,data_type, &
                mypart%myreceive%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

           end do

        else if(mode .eq. 2) then

           do i = 1, comm_loop

              sendElements = mypart%mysend%part_send(i)%numItems
              receiveElements = mypart%myreceive%part_receive(i)%numItems

              do k=1,sendElements
                 do j=1,size1
                    do n=1,size2
                       data_send(n,j,k,i) = struct(n,j,mypart%mysend%part_send(i)%items(k))
                    end do
                 end do
              end do

              srequests(i) = MPI_REQUEST_NULL
              rrequests(i) = MPI_REQUEST_NULL

              tag = newtag + mypart%mysend%process(i) + my_id 

              ! Non-blocking sends
              call MPI_Isend(data_send(1,1,1,i),sendElements*size1*size2,data_type, &
                 mypart%mysend%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              tag = newtag + my_id + mypart%mysend%process(i) 
              ! Non-blocking receives
              call MPI_Irecv(data_receive(1,1,1,i),receiveElements*size1*size2,data_type, &
                 mypart%myreceive%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

           end do

        else if(mode .eq. 3) then

           do i = 1, comm_loop

              sendElements = mypart%mysend%node_send(i)%numItems
              receiveElements = mypart%mysend%node_send(i)%numItems 

              do k=1,sendElements
                 do j=1,size1
                    do n=1,size2
                       data_send(n,j,k,i) = struct(n,j,mypart%mysend%part_send(i)%items(k))
                    end do
                 end do
              end do

              srequests(i) = MPI_REQUEST_NULL
              rrequests(i) = MPI_REQUEST_NULL

              tag = newtag + mypart%mysend%process(i) + my_id 

              ! Non-blocking sends
              call MPI_Isend(data_send(1,1,1,i),sendElements*size1*size2,data_type, &
                mypart%mysend%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              tag = newtag + my_id + mypart%mysend%process(i) 
              ! Non-blocking receives
              call MPI_Irecv(data_receive(1,1,1,i),receiveElements*size1*size2, &
                 data_type,mypart%myreceive%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

           end do


        else
           write(6,*)'error in send_recv_real3D_data', mode
           stop
        end if

        return


    end subroutine


!##########################################################################!
!***************************  start waitAll  ******************************!
!##########################################################################!
!##########################################################################!
!* This subroutine is necessary when there are some pending non blocking  *!
!* communications                                                         *!
!##########################################################################!

     subroutine waitAll(struct, size1, mode, requests)

        implicit none

        ! arguments 
        integer size1, mode
        real struct(size1,1)

        ! local
!        integer i,j,sendElements,receiveElements, comm_loop,k, tag
!        real, dimension(:,:),allocatable :: data_send, data_receive 

        ! MPI
        integer, dimension(:) :: requests
        integer ierr
        integer counter,n

        !integer, dimension(:), allocatable :: requests, statuses

        ! communication info
        !type (myPart) mypart

        !   requests(i) = MPI_REQUEST_NULL
        !end do

        ! Wait for all non-blocking communications to complete
        allocate(statuses(MPI_STATUS_SIZE * size(requests)))

        call mpi_waitall(size(requests), requests, statuses, ierr)
        if(ierr .ne. MPI_SUCCESS) then
           call abort
        end if
        do n=0,mypart%mysend%sends-1
!              write(6,*)'sender myrank',statuses(MPI_STATUS_SIZE*n+1)%MPI_SOURCE,my_id
!              write(6,*)'tag myrank',statuses(MPI_STATUS_SIZE*n+1)%MPI_TAG,my_id
!              write(6,*)'error myrank',statuses(MPI_STATUS_SIZE*n+1)%MPI_ERROR,my_id
              write(6,*)'count myrank',n,statuses(MPI_STATUS_SIZE*n+1),my_id
              write(6,*)'cancelled myrank',n,statuses(MPI_STATUS_SIZE*n+2),my_id
              write(6,*)'sender myrank',n,statuses(MPI_STATUS_SIZE*n+3),my_id
              write(6,*)'tag myrank',n,statuses(MPI_STATUS_SIZE*n+4),my_id
              write(6,*)'error myrank',n,statuses(MPI_STATUS_SIZE*n+5),my_id
!              call MPI_GET_COUNT(statuses(MPI_STATUS_SIZE*n+1), MPI_REAL, counter,ierr)
!              write(6,*)'counter:',counter,my_id
        end do
        deallocate(statuses)

        return

     end subroutine


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
!**********************  start send_recv_test  ***********************!
!##########################################################################!
!##########################################################################!
!* This subroutine is usefull for send or/and receve real data to/from    *!
!* others processes                                                       *!
!##########################################################################!

     subroutine send_recv_test(srequests,rrequests,data_receive)

        
        implicit none

        ! arguments 
        integer, dimension(:) :: srequests,rrequests

        !common
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

        ! local
        integer i,j,sendElements,receiveElements, comm_loop,k, tag, newtag
        !real, dimension(:),allocatable :: data_send, data_receive 
        real, dimension(:) :: data_receive 
        real struct

        ! MPI
        integer ierr, data_type

        struct = 1

        comm_loop = mypart%mysend%sends
        data_type = MPI_REAL

        newtag = next_mpi_tag()

!        allocate(data_send(comm_loop))
!        allocate(data_receive(comm_loop))

!        write(6,*)'max_send',mypart%mysend%maxItems,my_id
!        write(6,*)'max_receive',mypart%myreceive%maxItems,my_id

           do i = 1, comm_loop

              !sendElements = mypart%mysend%all_send(i)%numItems
              !receiveElements = mypart%myreceive%all_receive(i)%numItems
              sendElements = 1 
              receiveElements = 1

!              data_send(i) = struct


              srequests(i) = MPI_REQUEST_NULL
              rrequests(i) = MPI_REQUEST_NULL

              tag = newtag + mypart%mysend%process(i)! + my_id 

              ! Non-blocking sends
              !call MPI_Isend(data_send(i), sendElements, data_type, mypart%mysend%process(i), tag, MPI_COMM_WORLD, srequests(i), ierr)
              call MPI_Isend(struct,sendElements,data_type, &
                mypart%mysend%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              write(6,*)'isend',sendElements,mypart%mysend%process(i),tag,my_id
 
              tag = newtag + my_id ! + mypart%mysend%process(i) 
              ! Non-blocking receives
              !call MPI_Irecv(data_receive(i), receiveElements, data_type, mypart%myreceive%process(i), tag, MPI_COMM_WORLD, rrequests(i), ierr)
              call MPI_Irecv(data_receive(i),receiveElements,data_type, &
                mypart%myreceive%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              write(6,*)'irecv',receiveElements,mypart%myreceive%process(i),tag,my_id

!              struct = data_receive(i)

           end do

!        call waitAny(srequests)
!        call MPI_Barrier(MPI_COMM_WORLD,ierr)
!        call waitAny(rrequests)
!        call MPI_Barrier(MPI_COMM_WORLD,ierr)


!        deallocate(data_send)
!        deallocate(data_receive)

        return

    end subroutine

!================================================================

    subroutine rebuild_struct_2D_max_i(struct,data_receive,mode)


        implicit none

        !arguments
        integer mode
        integer struct(1)
        integer, dimension(:,:) :: data_receive

        !local
        integer receiveElements,i,j,k

        if(mode .eq. 1) then


           do i=1,mypart%myreceive%receives

              receiveElements = mypart%myreceive%all_receive(i)%numItems

              do k=1,receiveElements
                 struct(mypart%myreceive%all_receive(i)%items(k)) = &
                  max(struct(mypart%myreceive%part_receive(i)%items(k)), data_receive(k,i))
              end do

           end do

        else if( mode .eq. 2) then

           do i=1,mypart%myreceive%receives

           receiveElements = mypart%myreceive%part_receive(i)%numItems

              do k=1,receiveElements
                 struct(mypart%myreceive%part_receive(i)%items(k)) = &
                 max(struct(mypart%myreceive%part_receive(i)%items(k)),data_receive(k,i))
              end do
           end do

        else if( mode .eq. 3) then

           do i=1,mypart%myreceive%receives

           receiveElements = mypart%mysend%node_send(i)%numItems

              do k=1,receiveElements
                 struct(mypart%mysend%node_send(i)%items(k)) = max(struct(mypart%mysend%node_send(i)%items(k)),data_receive(k,i))
              end do
           end do

        else
           write(6,*)'error in rebuild_int_max1, mode=',mode
        end if
        

        return


    end subroutine

!================================================================

    subroutine rebuild_struct_2D_max_r(struct,data_receive,mode)


        implicit none

        !arguments
        integer mode
        real struct(1)
        real, dimension(:,:) :: data_receive

        !local
        integer receiveElements,i,j,k

        if(mode .eq. 1) then


           do i=1,mypart%myreceive%receives

              receiveElements = mypart%myreceive%all_receive(i)%numItems

              do k=1,receiveElements
                 struct(mypart%myreceive%all_receive(i)%items(k)) = &
                max(struct(mypart%myreceive%part_receive(i)%items(k)), data_receive(k,i))
              end do

           end do

        else if( mode .eq. 2) then

           do i=1,mypart%myreceive%receives

           receiveElements = mypart%myreceive%part_receive(i)%numItems

              do k=1,receiveElements
                 struct(mypart%myreceive%part_receive(i)%items(k)) = &
                max(struct(mypart%myreceive%part_receive(i)%items(k)),data_receive(k,i))
              end do
           end do

        else if( mode .eq. 3) then

           do i=1,mypart%myreceive%receives

           receiveElements = mypart%mysend%node_send(i)%numItems

              do k=1,receiveElements
                 struct(mypart%mysend%node_send(i)%items(k)) = max(struct(mypart%mysend%node_send(i)%items(k)), data_receive(k,i))
              end do
           end do

        else
           write(6,*)'error in rebuild_int_max1, mode=',mode
        end if
        

        return


    end subroutine

!================================================================

    subroutine rebuild_struct_2D_min_i(struct,data_receive,mode)


        implicit none

        !arguments
        integer mode
        integer struct(1)
        integer, dimension(:,:) :: data_receive

        !local
        integer receiveElements,i,j,k

        if(mode .eq. 1) then


           do i=1,mypart%myreceive%receives

              receiveElements = mypart%myreceive%all_receive(i)%numItems

              do k=1,receiveElements
                 struct(mypart%myreceive%all_receive(i)%items(k)) = &
                min(struct(mypart%myreceive%part_receive(i)%items(k)), data_receive(k,i))
              end do

           end do

        else if( mode .eq. 2) then

           do i=1,mypart%myreceive%receives

           receiveElements = mypart%myreceive%part_receive(i)%numItems

              do k=1,receiveElements
                 struct(mypart%myreceive%part_receive(i)%items(k)) = &
                min(struct(mypart%myreceive%part_receive(i)%items(k)),data_receive(k,i))
              end do
           end do

        else if( mode .eq. 3) then

           do i=1,mypart%myreceive%receives

           receiveElements = mypart%mysend%node_send(i)%numItems

              do k=1,receiveElements
                 struct(mypart%mysend%node_send(i)%items(k)) = &
                min(struct(mypart%mysend%node_send(i)%items(k)), data_receive(k,i))
              end do
           end do

        else
           write(6,*)'error in rebuild_int_min1, mode=',mode
        end if
        

        return


    end subroutine

!================================================================

    subroutine rebuild_struct_2D_min_r(struct,data_receive,mode)


        implicit none

        !arguments
        integer mode
        real struct(1)
        real, dimension(:,:) :: data_receive

        !local
        integer receiveElements,i,j,k

        if(mode .eq. 1) then


           do i=1,mypart%myreceive%receives

              receiveElements = mypart%myreceive%all_receive(i)%numItems

              do k=1,receiveElements
                 struct(mypart%myreceive%all_receive(i)%items(k)) = &
                min(struct(mypart%myreceive%part_receive(i)%items(k)), data_receive(k,i))
              end do

           end do

        else if( mode .eq. 2) then

           do i=1,mypart%myreceive%receives

           receiveElements = mypart%myreceive%part_receive(i)%numItems

              do k=1,receiveElements
                 struct(mypart%myreceive%part_receive(i)%items(k)) = &
                min(struct(mypart%myreceive%part_receive(i)%items(k)),data_receive(k,i))
              end do
           end do

        else if( mode .eq. 3) then

           do i=1,mypart%myreceive%receives

           receiveElements = mypart%mysend%node_send(i)%numItems

              do k=1,receiveElements
                 struct(mypart%mysend%node_send(i)%items(k)) = &
                min(struct(mypart%mysend%node_send(i)%items(k)), data_receive(k,i))
              end do
           end do

        else
           write(6,*)'error in rebuild_int_min1, mode=',mode
        end if
        

        return


    end subroutine

!================================================================

    subroutine rebuild_struct_2D_sum_r(struct,data_receive,mode)


        implicit none

        !arguments
        integer mode
        real struct(1)
        real, dimension(:,:) :: data_receive

        !local
        integer receiveElements,i,j,k

        if(mode .eq. 1) then


           do i=1,mypart%myreceive%receives

                 receiveElements = mypart%myreceive%all_receive(i)%numItems

              do k=1,receiveElements
                 struct(mypart%myreceive%all_receive(i)%items(k)) = &
                struct(mypart%myreceive%all_receive(i)%items(k)) + data_receive(k,i)
              end do

           end do

        else if( mode .eq. 2) then

           do i=1,mypart%myreceive%receives

              receiveElements = mypart%myreceive%part_receive(i)%numItems

              do k=1,receiveElements
                 struct(mypart%myreceive%part_receive(i)%items(k)) = &
                struct(mypart%myreceive%part_receive(i)%items(k)) + data_receive(k,i)
              end do
           end do

        else if( mode .eq. 3) then

           do i=1,mypart%myreceive%receives

              receiveElements = mypart%mysend%node_send(i)%numItems

              do k=1,receiveElements
                 struct(mypart%mysend%node_send(i)%items(k)) = &
                struct(mypart%mysend%node_send(i)%items(k)) + data_receive(k,i)
              end do
           end do

        else
           write(6,*)'error in Rebuild_struct_2D_sum, mode=',mode
        end if
        

        return


    end subroutine

!================================================================

    subroutine rebuild_struct_2D_sum_d(struct,data_receive,mode)


        implicit none

        !arguments
        integer mode
        double precision struct(1)
        double precision, dimension(:,:) :: data_receive

        !local
        integer receiveElements,i,j,k

        if(mode .eq. 1) then


           do i=1,mypart%myreceive%receives

              receiveElements = mypart%myreceive%all_receive(i)%numItems

              do k=1,receiveElements
                 struct(mypart%myreceive%all_receive(i)%items(k)) = &
                struct(mypart%myreceive%all_receive(i)%items(k)) + data_receive(k,i)
              end do

           end do

        else if( mode .eq. 2) then

           do i=1,mypart%myreceive%receives

           receiveElements = mypart%myreceive%part_receive(i)%numItems

              do k=1,receiveElements
                 struct(mypart%myreceive%part_receive(i)%items(k)) = &
                struct(mypart%myreceive%part_receive(i)%items(k)) + data_receive(k,i)
              end do
           end do

        else if( mode .eq. 3) then

           do i=1,mypart%myreceive%receives

           receiveElements = mypart%mysend%node_send(i)%numItems

              do k=1,receiveElements
                 struct(mypart%mysend%node_send(i)%items(k)) = &
                struct(mypart%mysend%node_send(i)%items(k)) + data_receive(k,i)
              end do
           end do

        else
           write(6,*)'error in rebuild_int_min1, mode=',mode
        end if
        

        return


    end subroutine



!================================================================

    subroutine rebuild_struct_2D_sum_i(struct,data_receive,mode)


        implicit none

        !arguments
        integer mode
        integer struct(1)
        integer, dimension(:,:) :: data_receive

        !local
        integer receiveElements,i,j,k

        if(mode .eq. 1) then


           do i=1,mypart%myreceive%receives

              receiveElements = mypart%myreceive%all_receive(i)%numItems

              do k=1,receiveElements
                 struct(mypart%myreceive%all_receive(i)%items(k)) = &
                struct(mypart%myreceive%all_receive(i)%items(k)) + data_receive(k,i)
              end do

           end do

        else if( mode .eq. 2) then

           do i=1,mypart%myreceive%receives

           receiveElements = mypart%myreceive%part_receive(i)%numItems

              do k=1,receiveElements
                 struct(mypart%myreceive%part_receive(i)%items(k)) = &
                struct(mypart%myreceive%part_receive(i)%items(k)) + data_receive(k,i)
              end do
           end do

        else if( mode .eq. 3) then

           do i=1,mypart%myreceive%receives

           receiveElements = mypart%mysend%node_send(i)%numItems

              do k=1,receiveElements
                 struct(mypart%mysend%node_send(i)%items(k)) = &
                struct(mypart%mysend%node_send(i)%items(k)) + data_receive(k,i)
              end do
           end do

        else
           write(6,*)'error in rebuild_int_min1, mode=',mode
        end if
        

        return


    end subroutine

!================================================================

    subroutine rebuild_struct_3D_sum_r(struct,size1,data_receive,mode)


        implicit none

        !arguments
        integer size1,mode
        real struct(size1,1)
        real, dimension(:,:,:) :: data_receive

        !local
        integer receiveElements,i,j,k

        if(mode .eq. 1) then


           do i=1,mypart%myreceive%receives

              receiveElements = mypart%myreceive%all_receive(i)%numItems

              do k=1,receiveElements
                 do j=1,size1
                    struct(j,mypart%myreceive%all_receive(i)%items(k)) = &
                struct(j,mypart%myreceive%all_receive(i)%items(k)) + data_receive(j,k,i)
                 end do
              end do

           end do

        else if( mode .eq. 2) then

           do i=1,mypart%myreceive%receives

           receiveElements = mypart%myreceive%part_receive(i)%numItems

              do k=1,receiveElements
                 do j=1,size1
                    struct(j,mypart%myreceive%part_receive(i)%items(k)) = &
                struct(j,mypart%myreceive%part_receive(i)%items(k)) + data_receive(j,k,i)
                 end do
              end do
           end do

        else if( mode .eq. 3) then

           do i=1,mypart%myreceive%receives

           receiveElements = mypart%mysend%node_send(i)%numItems

              do k=1,receiveElements
                 do j=1,size1
                    struct(j,mypart%mysend%node_send(i)%items(k)) = &
                struct(j,mypart%mysend%node_send(i)%items(k)) + data_receive(j,k,i)
                 end do
              end do
           end do

        else
           write(6,*)'error in rebuild_RealHalo, mode=',mode
        end if
        

        return


    end subroutine

!================================================================

    subroutine rebuild_struct_3D_sum_d(struct,size1,data_receive,mode)


        implicit none

        !arguments
        integer size1,mode
        double precision struct(size1,1)
        double precision, dimension(:,:,:) :: data_receive

        !local
        integer receiveElements,i,j,k

        if(mode .eq. 1) then


           do i=1,mypart%myreceive%receives

              receiveElements = mypart%myreceive%all_receive(i)%numItems

              do k=1,receiveElements
                 do j=1,size1
                    struct(j,mypart%myreceive%all_receive(i)%items(k)) = &
                struct(j,mypart%myreceive%all_receive(i)%items(k)) + data_receive(j,k,i)
                 end do
              end do

           end do

        else if( mode .eq. 2) then

           do i=1,mypart%myreceive%receives

           receiveElements = mypart%myreceive%part_receive(i)%numItems

              do k=1,receiveElements
                 do j=1,size1
                    struct(j,mypart%myreceive%part_receive(i)%items(k)) = &
                struct(j,mypart%myreceive%part_receive(i)%items(k)) + data_receive(j,k,i)
                 end do
              end do
           end do

        else if( mode .eq. 3) then

           do i=1,mypart%myreceive%receives

           receiveElements = mypart%mysend%node_send(i)%numItems

              do k=1,receiveElements
                 do j=1,size1
                    struct(j,mypart%mysend%node_send(i)%items(k)) = &
                struct(j,mypart%mysend%node_send(i)%items(k)) + data_receive(j,k,i)
                 end do
              end do
           end do

        else
           write(6,*)'error in rebuild_RealHalo, mode=',mode
        end if
        

        return


    end subroutine

!================================================================

    subroutine rebuild_struct_3D_2xSum_r(struct,size1,data_receive,mode)


        implicit none

        !arguments
        integer size1,mode
        real struct(size1,1)
        real, dimension(:,:,:) :: data_receive

        !local
        integer receiveElements,i,j,k

        if(mode .eq. 1) then


           do i=1,mypart%myreceive%receives

              receiveElements = mypart%myreceive%all_receive(i)%numItems

              do k=1,receiveElements
                 do j=1,size1
                    struct(j,mypart%myreceive%all_receive(i)%items(k)) = &
                2 * (struct(j,mypart%myreceive%all_receive(i)%items(k)) + data_receive(j,k,i))
                 end do
              end do

           end do

        else if( mode .eq. 2) then

           do i=1,mypart%myreceive%receives

           receiveElements = mypart%myreceive%part_receive(i)%numItems

              do k=1,receiveElements
                 do j=1,size1
                    struct(j,mypart%myreceive%part_receive(i)%items(k)) = &
                2 * (struct(j,mypart%myreceive%part_receive(i)%items(k)) + data_receive(j,k,i))
                 end do
              end do
           end do

        else if( mode .eq. 3) then

           do i=1,mypart%myreceive%receives

           receiveElements = mypart%mysend%node_send(i)%numItems

              do k=1,receiveElements
                 do j=1,size1
                    struct(j,mypart%mysend%node_send(i)%items(k)) = &
                2 * (struct(j,mypart%mysend%node_send(i)%items(k)) + data_receive(j,k,i))
                 end do
              end do
           end do

        else
           write(6,*)'error in rebuild_RealHalo, mode=',mode
        end if
        

        return


    end subroutine

!================================================================

    subroutine rebuild_struct_3D_2xSum_d(struct,size1,data_receive,mode)


        implicit none

        !arguments
        integer size1,mode
        double precision struct(size1,1)
        double precision, dimension(:,:,:) :: data_receive

        !local
        integer receiveElements,i,j,k

        if(mode .eq. 1) then


           do i=1,mypart%myreceive%receives

              receiveElements = mypart%myreceive%all_receive(i)%numItems

              do k=1,receiveElements
                 do j=1,size1
                    struct(j,mypart%myreceive%all_receive(i)%items(k)) = &
                2 * (struct(j,mypart%myreceive%all_receive(i)%items(k)) + data_receive(j,k,i))
                 end do
              end do

           end do

        else if( mode .eq. 2) then

           do i=1,mypart%myreceive%receives

           receiveElements = mypart%myreceive%part_receive(i)%numItems

              do k=1,receiveElements
                 do j=1,size1
                    struct(j,mypart%myreceive%part_receive(i)%items(k)) = &
                2 * (struct(j,mypart%myreceive%part_receive(i)%items(k)) + data_receive(j,k,i))
                 end do
              end do
           end do

        else if( mode .eq. 3) then

           do i=1,mypart%myreceive%receives

           receiveElements = mypart%mysend%node_send(i)%numItems

              do k=1,receiveElements
                 do j=1,size1
                    struct(j,mypart%mysend%node_send(i)%items(k)) = &
                2 * (struct(j,mypart%mysend%node_send(i)%items(k)) + data_receive(j,k,i))
                 end do
              end do
           end do

        else
           write(6,*)'error in rebuild_RealHalo, mode=',mode
        end if
        

        return


    end subroutine

!******************************************************************

    subroutine rebuild_DPHalo(struct,size1,data_receive,mode)


        implicit none

        !arguments
        integer size1,mode
        double precision struct(size1,1)
        double precision, dimension(:,:,:) :: data_receive

        !local
        integer receiveElements,i,j,k

        if(mode .eq. 1) then


           do i=1,mypart%myreceive%receives

              receiveElements = mypart%myreceive%all_receive(i)%numItems

              do k=1,receiveElements
                 do j=1,size1
                    struct(j,mypart%myreceive%all_receive(i)%items(k)) =  data_receive(j,k,i)
                 end do
              end do

           end do

        else if( mode .eq. 2) then

           do i=1,mypart%myreceive%receives

           receiveElements = mypart%myreceive%part_receive(i)%numItems

              do k=1,receiveElements
                 do j=1,size1
                    struct(j,mypart%myreceive%part_receive(i)%items(k)) =  data_receive(j,k,i)
                 end do
              end do
           end do

        else if( mode .eq. 3) then

           do i=1,mypart%myreceive%receives

           receiveElements = mypart%mysend%node_send(i)%numItems

              do k=1,receiveElements
                 do j=1,size1
                    struct(j,mypart%mysend%node_send(i)%items(k)) = &
                struct(j,mypart%mysend%node_send(i)%items(k)) + data_receive(j,k,i)
                 end do
              end do
           end do

        else
           write(6,*)'error in rebuild_DPHalo, mode=',mode
        end if
        

        return


    end subroutine


!******************************************************************

    subroutine rebuild_RealHalo(struct,size1,data_receive,mode)


        implicit none

        !arguments
        integer size1,mode
        real struct(size1,1)
        real, dimension(:,:,:) :: data_receive

        !local
        integer receiveElements,i,j,k

        if(mode .eq. 1) then


           do i=1,mypart%myreceive%receives

              receiveElements = mypart%myreceive%all_receive(i)%numItems

              do k=1,receiveElements
                 do j=1,size1
                    struct(j,mypart%myreceive%all_receive(i)%items(k)) =  data_receive(j,k,i)
                 end do
              end do

           end do

        else if( mode .eq. 2) then

           do i=1,mypart%myreceive%receives

           receiveElements = mypart%myreceive%part_receive(i)%numItems

              do k=1,receiveElements
                 do j=1,size1
                    struct(j,mypart%myreceive%part_receive(i)%items(k)) =  data_receive(j,k,i)
                 end do
              end do
           end do

        else if( mode .eq. 3) then

           do i=1,mypart%myreceive%receives

           receiveElements = mypart%mysend%node_send(i)%numItems

              do k=1,receiveElements
                 do j=1,size1
                    struct(j,mypart%mysend%node_send(i)%items(k)) = &
                struct(j,mypart%mysend%node_send(i)%items(k)) + data_receive(j,k,i)
                 end do
              end do
           end do

        else
           write(6,*)'error in rebuild_RealHalo, mode=',mode
        end if
        

        return


    end subroutine


    subroutine rebuild_Real3D_Halo(struct,size1,size2,data_receive,mode)


        implicit none

        !arguments
        integer size1,size2,mode
        real struct(size1,size2,1)
        real, dimension(:,:,:,:) :: data_receive

        !local
        integer receiveElements,i,j,k,n

        if(mode .eq. 1) then

           do i=1,mypart%myreceive%receives

              receiveElements = mypart%myreceive%all_receive(i)%numItems



!              if(my_id .eq. 1 .and. mypart%myreceive%process(i) .eq. 2 ) then
!                 open(unit=1900, file="1data_receive2.txt", action='write')
!              end if


              do k=1,receiveElements


!              if(my_id .eq. 1 .and. mypart%myreceive%process(i) .eq. 2 ) then
!                    write(1900,*),k,mypart%myreceive%all_receive(i)%items(k),myele%globalID(mypart%myreceive%all_receive(i)%items(k))
!              end if


                 do j=1,size1
                    do n=1,size2
                       struct(n,j,mypart%myreceive%all_receive(i)%items(k)) =  data_receive(n,j,k,i)


!              if(my_id .eq. 1 .and. mypart%myreceive%process(i) .eq. 2 ) then
!                          write(1900,*), data_receive(n,j,k,i)
!              end if


                    end do
                 end do
              end do


!              if(my_id .eq. 1 .and. mypart%myreceive%process(i) .eq. 2 ) then
!                 close(1900)
!                 stop
!              end if


           end do

        else if( mode .eq. 2) then

           do i=1,mypart%myreceive%receives

           receiveElements = mypart%myreceive%part_receive(i)%numItems

              do k=1,receiveElements
                 do j=1,size1
                    do n=1,size2
                       struct(n,j,mypart%myreceive%part_receive(i)%items(k)) =  data_receive(n,j,k,i)
                    end do
                 end do
              end do
           end do

        else
           write(6,*)'error in rebuild_Real3D_Halo, mode=',mode
        end if
        

        return


    end subroutine


end module mpi_communication
