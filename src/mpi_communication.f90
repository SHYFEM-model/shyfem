!-------------------------------------------------------------------------------------
   module mpi_communication
!-------------------------------------------------------------------------------------

   use mpi_communication_struct
   use mpi_common_struct
   use timing

!-------------------------------------------------------------------------------------

   interface exchange_struct_2d
      module procedure exchange_struct_2d_i  &
                       ,exchange_struct_2d_r !&
!                       ,exchange_struct_2d_d
   end interface exchange_struct_2d

   interface rebuild_struct_2d_sum
      module procedure rebuild_struct_2d_sum_i  &
                       ,rebuild_struct_2d_sum_r !&
!                       ,rebuild_struct_2d_sum_d
   end interface rebuild_struct_2d_sum

   interface exchange_struct_3d
      module procedure exchange_struct_3d_i &
                       ,exchange_struct_3d_r !&
!                       ,exchange_struct_3d_d
   end interface exchange_struct_3d

   interface rebuild_struct_3d_sum
      module procedure rebuild_struct_3d_sum_r !&
!                       ,rebuild_struct_3d_sum_d
   end interface rebuild_struct_3d_sum

   interface rebuild_struct_3d
      module procedure rebuild_struct_3d_r !& 
!                       ,rebuild_struct_3d_d
   end interface rebuild_struct_3d

   interface rebuild_struct_2d
      module procedure rebuild_struct_2d_i &
                       ,rebuild_struct_2d_r !&
!                       ,rebuild_struct_2d_d
   end interface rebuild_struct_2d

!-------------------------------------------------------------------------------------
   contains
!-------------------------------------------------------------------------------------

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
        double precision time1

        if(ln_timing) time1 = MPI_Wtime()

        comm_loop = domain%exchanges
        if(.not. allocated(data_receive)) then
          allocate(data_receive(maxNodes,comm_loop))
        else
          deallocate(data_receive)
          allocate(data_receive(maxNodes,comm_loop))
        end if
        data_receive = 0.

        newtag = next_mpi_tag()

        do i = 1, comm_loop

           sendNodes = domain%nodes%domainID
           receiveNodes = numberNodes(domain%process(i)+1)

           srequests(i) = MPI_REQUEST_NULL
           rrequests(i) = MPI_REQUEST_NULL

           tag = newtag + domain%process(i) + my_id 

           ! Non-blocking sends
           call MPI_Isend(domain%nodes%globalID(1),sendNodes,MPI_INTEGER, &
                domain%process(i),tag, MPI_COMM_WORLD,srequests(i),ierr)
           if(ierr .ne. MPI_SUCCESS) then
              write(6,*)'ierror:',ierr
           end if

           tag = newtag + my_id + domain%process(i) 
           ! Non-blocking receives
           call MPI_Irecv(data_receive(1,i),receiveNodes,MPI_INTEGER, &
                domain%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
           if(ierr .ne. MPI_SUCCESS) then
              write(6,*)'ierror:',ierr
           end if

        end do

        if(ln_timing) comm_time = comm_time + MPI_Wtime() - time1

        return

     end subroutine spread_nodes_to_neighbor

!##########################################################################!
!**********************  start send_recv_real_data  ***********************!
!##########################################################################!
!##########################################################################!
!* This subroutine is usefull for send or/and receve double precision data to/from  ***!
!* others processes                                                       *!
!##########################################################################!

     subroutine exchange_struct_2d_r(struct, mode, srequests, rrequests, data_send, data_receive)
        
        use basin

        implicit none

        ! arguments 
        integer mode
        double precision struct(nkn_local)

        !common
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

        ! local
        integer i,j,sendElements,receiveElements, comm_loop,k, tag, newtag
        double precision, dimension(:,:),allocatable :: data_send, data_receive 

        ! MPI
        integer ierr, data_type
        integer, dimension(:) :: srequests,rrequests
        double precision time1

        if(ln_timing) time1 = MPI_Wtime()

        comm_loop = domain%exchanges
        data_type = MPI_DOUBLE_PRECISION

        newtag = next_mpi_tag()

        if(.not. allocated(data_send)) then
          allocate(data_send(domain%maxItems,comm_loop))
          allocate(data_receive(domain%maxItems,comm_loop))
        end if
        data_receive = 0.

        if(mode .eq. 2) then

           do i = 1, comm_loop

              sendElements = domain%ele_send(i)%numItems
              receiveElements = domain%ele_receive(i)%numItems

              do k=1,sendElements
                 data_send(k,i) = struct(domain%ele_send(i)%items(k))
              end do

              srequests(i) = MPI_REQUEST_NULL
              rrequests(i) = MPI_REQUEST_NULL

              tag = newtag + domain%process(i) + my_id 

              ! Non-blocking sends
              call MPI_Isend(data_send(1,i),sendElements,data_type, &
                domain%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              tag = newtag + my_id + domain%process(i) 
              ! Non-blocking receives
              call MPI_Irecv(data_receive(1,i),receiveElements,data_type, &
                domain%process(i),tag, MPI_COMM_WORLD,rrequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

           end do

        else if(mode .eq. 3) then

           do i = 1, comm_loop

              sendElements = domain%node_send(i)%numItems
              receiveElements = domain%node_send(i)%numItems 

              do k=1,sendElements
                 data_send(k,i) = struct(domain%node_send(i)%items(k))
              end do


              srequests(i) = MPI_REQUEST_NULL
              rrequests(i) = MPI_REQUEST_NULL

              tag = newtag + domain%process(i) + my_id 

              ! Non-blocking sends
              call MPI_Isend(data_send(1,i),sendElements,data_type, &
                domain%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              tag = newtag + my_id + domain%process(i) 
              ! Non-blocking receives
              call MPI_Irecv(data_receive(1,i),receiveElements,data_type, &
                domain%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

           end do

        else if(mode .eq. 4) then

           do i = 1, comm_loop

              sendElements = domain%halo_send_e(i)%numItems
              receiveElements = domain%halo_recv_e(i)%numItems

              do k=1,sendElements
                 data_send(k,i) = struct(domain%halo_send_e(i)%items(k))
              end do

              srequests(i) = MPI_REQUEST_NULL
              rrequests(i) = MPI_REQUEST_NULL

              tag = newtag + domain%process(i) + my_id 

              ! Non-blocking sends
              call MPI_Isend(data_send(1,i),sendElements,data_type, &
                domain%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              tag = newtag + my_id + domain%process(i) 
              ! Non-blocking receives
              call MPI_Irecv(data_receive(1,i),receiveElements,data_type, &
                domain%process(i),tag, MPI_COMM_WORLD,rrequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

           end do

        else if(mode .eq. 5) then

           do i = 1, comm_loop

              sendElements = domain%halo_send_k(i)%numItems
              receiveElements = domain%halo_recv_k(i)%numItems

              do k=1,sendElements
                 data_send(k,i) = struct(domain%halo_send_k(i)%items(k))
              end do

              srequests(i) = MPI_REQUEST_NULL
              rrequests(i) = MPI_REQUEST_NULL

              tag = newtag + domain%process(i) + my_id 

              ! Non-blocking sends
              call MPI_Isend(data_send(1,i),sendElements,data_type, &
                domain%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              tag = newtag + my_id + domain%process(i) 
              ! Non-blocking receives
              call MPI_Irecv(data_receive(1,i),receiveElements,data_type, &
                domain%process(i),tag, MPI_COMM_WORLD,rrequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

           end do

        else
           write(6,*)'error in send_recv_int_data', mode
           stop
        end if

        if(ln_timing) comm_time = comm_time + MPI_Wtime() - time1

        return


    end subroutine

!##########################################################################!
!**********************  start send_recv_dp_data  ***********************!
!##########################################################################!
!##########################################################################!
!* This subroutine is usefull for send or/and receve double precision data !
!* to/from others processes                                               *!
!##########################################################################!

!     subroutine exchange_struct_2d_d(struct, mode, srequests, rrequests, data_send, data_receive)
!        
!        use basin
!
!        implicit none
!
!        ! arguments 
!        integer mode
!        double precision struct(nkn_local)
!
!        !common
!        integer itanf,itend,idt,nits,niter,it
!        common /femtim/ itanf,itend,idt,nits,niter,it
!
!        ! local
!        integer i,j,sendElements,receiveElements, comm_loop,k, tag, newtag
!        double precision, dimension(:,:),allocatable :: data_send, data_receive 
!
!        ! MPI
!        integer ierr, data_type
!        integer, dimension(:) :: srequests,rrequests
!
!        comm_loop = domain%exchanges
!        data_type = MPI_DOUBLE_PRECISION
!
!        newtag = next_mpi_tag()
!
!        if(.not. allocated(data_send)) then
!          allocate(data_send(domain%maxItems,comm_loop))
!          allocate(data_receive(domain%maxItems,comm_loop))
!        end if
!        data_receive = 0.
!
!        if(mode .eq. 2) then
!
!           do i = 1, comm_loop
!
!              sendElements = domain%ele_send(i)%numItems
!              receiveElements = domain%ele_receive(i)%numItems
!
!              do k=1,sendElements
!                 data_send(k,i) = struct(domain%ele_send(i)%items(k))
!              end do
!
!              srequests(i) = MPI_REQUEST_NULL
!              rrequests(i) = MPI_REQUEST_NULL
!
!              tag = newtag + domain%process(i) + my_id 
!
!              ! Non-blocking sends
!              call MPI_Isend(data_send(1,i),sendElements,data_type, &
!                domain%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
!              if(ierr .ne. MPI_SUCCESS) then
!                 write(6,*)'ierror:',ierr
!              end if
!
!              tag = newtag + my_id + domain%process(i) 
!              ! Non-blocking receives
!              call MPI_Irecv(data_receive(1,i),receiveElements,data_type, &
!                domain%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
!              if(ierr .ne. MPI_SUCCESS) then
!                 write(6,*)'ierror:',ierr
!              end if
!
!           end do
!
!        else if(mode .eq. 3) then
!
!           do i = 1, comm_loop
!
!              sendElements = domain%node_send(i)%numItems
!              receiveElements = domain%node_send(i)%numItems 
!
!              do k=1,sendElements
!                 data_send(k,i) = struct(domain%node_send(i)%items(k))
!              end do
!
!
!              srequests(i) = MPI_REQUEST_NULL
!              rrequests(i) = MPI_REQUEST_NULL
!
!              tag = newtag + domain%process(i) + my_id 
!
!              ! Non-blocking sends
!              call MPI_Isend(data_send(1,i),sendElements,data_type, &
!                domain%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
!              if(ierr .ne. MPI_SUCCESS) then
!                 write(6,*)'ierror:',ierr
!              end if
!
!              tag = newtag + my_id + domain%process(i) 
!              ! Non-blocking receives
!              call MPI_Irecv(data_receive(1,i),receiveElements,data_type, &
!                domain%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
!              if(ierr .ne. MPI_SUCCESS) then
!                 write(6,*)'ierror:',ierr
!              end if
!
!           end do
!
!        else if(mode .eq. 4) then
!
!           do i = 1, comm_loop
!
!              sendElements = domain%halo_send_e(i)%numItems
!              receiveElements = domain%halo_recv_e(i)%numItems
!
!              do k=1,sendElements
!                 data_send(k,i) = struct(domain%halo_send_e(i)%items(k))
!              end do
!
!              srequests(i) = MPI_REQUEST_NULL
!              rrequests(i) = MPI_REQUEST_NULL
!
!              tag = newtag + domain%process(i) + my_id 
!
!              ! Non-blocking sends
!              call MPI_Isend(data_send(1,i),sendElements,data_type, &
!                domain%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
!              if(ierr .ne. MPI_SUCCESS) then
!                 write(6,*)'ierror:',ierr
!              end if
!
!              tag = newtag + my_id + domain%process(i) 
!              ! Non-blocking receives
!              call MPI_Irecv(data_receive(1,i),receiveElements,data_type, &
!                domain%process(i),tag, MPI_COMM_WORLD,rrequests(i),ierr)
!              if(ierr .ne. MPI_SUCCESS) then
!                 write(6,*)'ierror:',ierr
!              end if
!
!           end do
!
!        else if(mode .eq. 5) then
!
!           do i = 1, comm_loop
!
!              sendElements = domain%halo_send_k(i)%numItems
!              receiveElements = domain%halo_recv_k(i)%numItems
!
!              do k=1,sendElements
!                 data_send(k,i) = struct(domain%halo_send_k(i)%items(k))
!              end do
!
!              srequests(i) = MPI_REQUEST_NULL
!              rrequests(i) = MPI_REQUEST_NULL
!
!              tag = newtag + domain%process(i) + my_id 
!
!              ! Non-blocking sends
!              call MPI_Isend(data_send(1,i),sendElements,data_type, &
!                domain%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
!              if(ierr .ne. MPI_SUCCESS) then
!                 write(6,*)'ierror:',ierr
!              end if
!
!              tag = newtag + my_id + domain%process(i) 
!              ! Non-blocking receives
!              call MPI_Irecv(data_receive(1,i),receiveElements,data_type, &
!                domain%process(i),tag, MPI_COMM_WORLD,rrequests(i),ierr)
!              if(ierr .ne. MPI_SUCCESS) then
!                 write(6,*)'ierror:',ierr
!              end if
!
!           end do
!
!        else
!           write(6,*)'error in send_recv_int_data', mode
!           stop
!        end if
!
!        return
!
!
!    end subroutine



!##########################################################################!
!**********************  start send_recv_int_data  ***********************!
!##########################################################################!
!##########################################################################!
!* This subroutine is usefull for send or/and receve integer data to/from *!
!* others processes                                                       *!
!##########################################################################!

     subroutine exchange_struct_2d_i(struct, mode, srequests, rrequests, data_send, data_receive)
        
        use basin

        implicit none

        ! arguments 
        integer mode
        integer struct(:)

        !common
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

        ! local
        integer i,j,sendElements,receiveElements, comm_loop,k, tag, newtag
        integer, dimension(:,:),allocatable :: data_send, data_receive 

        ! MPI
        integer ierr, data_type
        integer, dimension(:) :: srequests,rrequests
        double precision time1

        if(ln_timing) time1 = MPI_Wtime()

        comm_loop = domain%exchanges
        data_type = MPI_INTEGER

        newtag = next_mpi_tag()

        if(.not. allocated(data_send)) then
          allocate(data_send(domain%maxItems,comm_loop))
          allocate(data_receive(domain%maxItems,comm_loop))
        end if
        data_receive = 0.

        if(mode .eq. 2) then

           do i = 1, comm_loop

              sendElements = domain%ele_send(i)%numItems
              receiveElements = domain%ele_receive(i)%numItems

              do k=1,sendElements
                 data_send(k,i) = struct(domain%ele_send(i)%items(k))
              end do

              srequests(i) = MPI_REQUEST_NULL
              rrequests(i) = MPI_REQUEST_NULL

              tag = newtag + domain%process(i) + my_id 

              ! Non-blocking sends
              call MPI_Isend(data_send(1,i),sendElements,data_type, &
                domain%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              tag = newtag + my_id + domain%process(i) 
              ! Non-blocking receives
              call MPI_Irecv(data_receive(1,i),receiveElements,data_type, &
                domain%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

           end do

        else if(mode .eq. 3) then

           do i = 1, comm_loop

              sendElements = domain%node_send(i)%numItems
              receiveElements = domain%node_send(i)%numItems 

              do k=1,sendElements
                 data_send(k,i) = struct(domain%node_send(i)%items(k))
              end do


              srequests(i) = MPI_REQUEST_NULL
              rrequests(i) = MPI_REQUEST_NULL

              tag = newtag + domain%process(i) + my_id 

              ! Non-blocking sends
              call MPI_Isend(data_send(1,i),sendElements,data_type, &
                domain%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              tag = newtag + my_id + domain%process(i) 
              ! Non-blocking receives
              call MPI_Irecv(data_receive(1,i),receiveElements,data_type, &
                domain%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

           end do

        else if(mode .eq. 4) then

           do i = 1, comm_loop

              sendElements = domain%halo_send_e(i)%numItems
              receiveElements = domain%halo_recv_e(i)%numItems

              do k=1,sendElements
                 data_send(k,i) = struct(domain%halo_send_e(i)%items(k))
              end do

              srequests(i) = MPI_REQUEST_NULL
              rrequests(i) = MPI_REQUEST_NULL

              tag = newtag + domain%process(i) + my_id 

              ! Non-blocking sends
              call MPI_Isend(data_send(1,i),sendElements,data_type, &
                domain%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              tag = newtag + my_id + domain%process(i) 
              ! Non-blocking receives
              call MPI_Irecv(data_receive(1,i),receiveElements,data_type, &
                domain%process(i),tag, MPI_COMM_WORLD,rrequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

           end do

        else if(mode .eq. 5) then

           do i = 1, comm_loop

              sendElements = domain%halo_send_k(i)%numItems
              receiveElements = domain%halo_recv_k(i)%numItems

              do k=1,sendElements
                 data_send(k,i) = struct(domain%halo_send_k(i)%items(k))
              end do

              srequests(i) = MPI_REQUEST_NULL
              rrequests(i) = MPI_REQUEST_NULL

              tag = newtag + domain%process(i) + my_id 

              ! Non-blocking sends
              call MPI_Isend(data_send(1,i),sendElements,data_type, &
                domain%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              tag = newtag + my_id + domain%process(i) 
              ! Non-blocking receives
              call MPI_Irecv(data_receive(1,i),receiveElements,data_type, &
                domain%process(i),tag, MPI_COMM_WORLD,rrequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

           end do

        else
           write(6,*)'error in send_recv_int_data', mode
           stop
        end if

        if(ln_timing) comm_time = comm_time + MPI_Wtime() - time1

        return


    end subroutine

!##########################################################################!
!**********************  start send_recv_real_data  ***********************!
!##########################################################################!
!##########################################################################!
!* This subroutine is usefull for send or/and receve double precision data to/from    *!
!* others processes                                                       *!
!##########################################################################!

     subroutine exchange_struct_3d_r(struct, size1,size2, mode, srequests, rrequests, data_send, data_receive)
        
        implicit none

        ! arguments 
        integer size1, size2, mode
        double precision struct(size1,size2)

        !common
        !integer itanf,itend,idt,nits,niter,it
        !common /femtim/ itanf,itend,idt,nits,niter,it

        ! local
        integer i,j,sendElements,receiveElements, comm_loop,k, tag, newtag
        double precision, dimension(:,:,:),allocatable :: data_send, data_receive 

        ! MPI
        integer ierr, data_type
        !integer, dimension(:) :: srequests,rrequests
        integer srequests(domain%exchanges)
        integer rrequests(domain%exchanges)
        double precision time1

        if(ln_timing) time1 = MPI_Wtime()

        comm_loop = domain%exchanges
        data_type = MPI_DOUBLE_PRECISION

        newtag = next_mpi_tag()
        
        if(.not. allocated(data_send)) then
          allocate(data_send(size1,domain%maxItems,comm_loop))
          allocate(data_receive(size1,domain%maxItems,comm_loop))
        end if

        data_receive = 0.
        if(mode .eq. 2) then

           do i = 1, comm_loop

              sendElements = domain%ele_send(i)%numItems
              receiveElements = domain%ele_receive(i)%numItems

              do k=1,sendElements
                 do j=1,size1
                    data_send(j,k,i) = struct(j,domain%ele_send(i)%items(k))
                 end do
              end do

              srequests(i) = MPI_REQUEST_NULL
              rrequests(i) = MPI_REQUEST_NULL

              tag = newtag + domain%process(i) + my_id 

              ! Non-blocking sends
              call MPI_Isend(data_send(1,1,i),sendElements*size1,data_type, &
                domain%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              tag = newtag + my_id + domain%process(i) 
              ! Non-blocking receives
              call MPI_Irecv(data_receive(1,1,i),receiveElements*size1,data_type, &
                domain%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

           end do

        else if(mode .eq. 3) then

           do i = 1, comm_loop

              sendElements = domain%node_send(i)%numItems
              receiveElements = domain%node_send(i)%numItems 

              do k=1,sendElements
                 do j=1,size1
                    data_send(j,k,i) = struct(j,domain%node_send(i)%items(k))
                 end do
              end do


              srequests(i) = MPI_REQUEST_NULL
              rrequests(i) = MPI_REQUEST_NULL

              tag = newtag + domain%process(i) + my_id 

              ! Non-blocking sends
              call MPI_Isend(data_send(1,1,i),sendElements*size1,data_type, &
                domain%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              tag = newtag + my_id + domain%process(i) 
              ! Non-blocking receives
              call MPI_Irecv(data_receive(1,1,i),receiveElements*size1,data_type, &
                domain%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

           end do

        else if(mode .eq. 4) then

           do i = 1, comm_loop

              sendElements = domain%halo_send_e(i)%numItems
              receiveElements = domain%halo_recv_e(i)%numItems

              do k=1,sendElements
                 do j=1,size1
                    data_send(j,k,i) = struct(j,domain%halo_send_e(i)%items(k))
                 end do
              end do

              srequests(i) = MPI_REQUEST_NULL
              rrequests(i) = MPI_REQUEST_NULL

              tag = newtag + domain%process(i) + my_id 

              ! Non-blocking sends
              call MPI_Isend(data_send(1,1,i),sendElements*size1,data_type, &
                domain%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              tag = newtag + my_id + domain%process(i) 
              ! Non-blocking receives
              call MPI_Irecv(data_receive(1,1,i),receiveElements*size1,data_type, &
                domain%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

           end do

        else if(mode .eq. 5) then

           do i = 1, comm_loop

              sendElements = domain%halo_send_k(i)%numItems
              receiveElements = domain%halo_recv_k(i)%numItems

              do k=1,sendElements
                 do j=1,size1
                    data_send(j,k,i) = struct(j,domain%halo_send_k(i)%items(k))
                 end do
              end do

              srequests(i) = MPI_REQUEST_NULL
              rrequests(i) = MPI_REQUEST_NULL

              tag = newtag + domain%process(i) + my_id 

              ! Non-blocking sends
              call MPI_Isend(data_send(1,1,i),sendElements*size1,data_type, &
                domain%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              tag = newtag + my_id + domain%process(i) 
              ! Non-blocking receives
              call MPI_Irecv(data_receive(1,1,i),receiveElements*size1,data_type, &
                domain%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

           end do

        else
           write(6,*)'error in send_recv_real_data', mode
           stop
        end if

        if(ln_timing) comm_time = comm_time + MPI_Wtime() - time1

        return


    end subroutine


     subroutine exchange_struct_3d_i(struct, size1,size2, mode, srequests, rrequests, data_send, data_receive)
        
        implicit none

        ! arguments 
        integer size1, size2, mode
        integer struct(size1,size2)

        !common
        !integer itanf,itend,idt,nits,niter,it
        !common /femtim/ itanf,itend,idt,nits,niter,it

        ! local
        integer i,j,sendElements,receiveElements, comm_loop,k, tag, newtag
        integer, dimension(:,:,:),allocatable :: data_send, data_receive 

        ! MPI
        integer ierr, data_type
        !integer, dimension(:) :: srequests,rrequests
        integer srequests(domain%exchanges)
        integer rrequests(domain%exchanges)
        double precision time1

        if(ln_timing) time1 = MPI_Wtime()

        comm_loop = domain%exchanges
        data_type = MPI_DOUBLE_PRECISION

        newtag = next_mpi_tag()
        
        if(.not. allocated(data_send)) then
          allocate(data_send(size1,domain%maxItems,comm_loop))
          allocate(data_receive(size1,domain%maxItems,comm_loop))
        end if
        data_receive = 0.

        if(mode .eq. 2) then

           do i = 1, comm_loop

              sendElements = domain%ele_send(i)%numItems
              receiveElements = domain%ele_receive(i)%numItems

              do k=1,sendElements
                 do j=1,size1
                    data_send(j,k,i) = struct(j,domain%ele_send(i)%items(k))
                 end do
              end do

              srequests(i) = MPI_REQUEST_NULL
              rrequests(i) = MPI_REQUEST_NULL

              tag = newtag + domain%process(i) + my_id 

              ! Non-blocking sends
              call MPI_Isend(data_send(1,1,i),sendElements*size1,data_type, &
                domain%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              tag = newtag + my_id + domain%process(i) 
              ! Non-blocking receives
              call MPI_Irecv(data_receive(1,1,i),receiveElements*size1,data_type, &
                domain%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

           end do

        else if(mode .eq. 3) then

           do i = 1, comm_loop

              sendElements = domain%node_send(i)%numItems
              receiveElements = domain%node_send(i)%numItems 

              do k=1,sendElements
                 do j=1,size1
                    data_send(j,k,i) = struct(j,domain%node_send(i)%items(k))
                 end do
              end do


              srequests(i) = MPI_REQUEST_NULL
              rrequests(i) = MPI_REQUEST_NULL

              tag = newtag + domain%process(i) + my_id 

              ! Non-blocking sends
              call MPI_Isend(data_send(1,1,i),sendElements*size1,data_type, &
                domain%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              tag = newtag + my_id + domain%process(i) 
              ! Non-blocking receives
              call MPI_Irecv(data_receive(1,1,i),receiveElements*size1,data_type, &
                domain%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

           end do

        else if(mode .eq. 4) then

           do i = 1, comm_loop

              sendElements = domain%halo_send_e(i)%numItems
              receiveElements = domain%halo_recv_e(i)%numItems

              do k=1,sendElements
                 do j=1,size1
                    data_send(j,k,i) = struct(j,domain%halo_send_e(i)%items(k))
                 end do
              end do

              srequests(i) = MPI_REQUEST_NULL
              rrequests(i) = MPI_REQUEST_NULL

              tag = newtag + domain%process(i) + my_id 

              ! Non-blocking sends
              call MPI_Isend(data_send(1,1,i),sendElements*size1,data_type, &
                domain%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              tag = newtag + my_id + domain%process(i) 
              ! Non-blocking receives
              call MPI_Irecv(data_receive(1,1,i),receiveElements*size1,data_type, &
                domain%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

           end do

        else if(mode .eq. 5) then

           do i = 1, comm_loop

              sendElements = domain%halo_send_k(i)%numItems
              receiveElements = domain%halo_recv_k(i)%numItems

              do k=1,sendElements
                 do j=1,size1
                    data_send(j,k,i) = struct(j,domain%halo_send_k(i)%items(k))
                 end do
              end do

              srequests(i) = MPI_REQUEST_NULL
              rrequests(i) = MPI_REQUEST_NULL

              tag = newtag + domain%process(i) + my_id 

              ! Non-blocking sends
              call MPI_Isend(data_send(1,1,i),sendElements*size1,data_type, &
                domain%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              tag = newtag + my_id + domain%process(i) 
              ! Non-blocking receives
              call MPI_Irecv(data_receive(1,1,i),receiveElements*size1,data_type, &
                domain%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

           end do

        else
           write(6,*)'error in send_recv_real_data', mode
           stop
        end if

        if(ln_timing) comm_time = comm_time + MPI_Wtime() - time1

        return


    end subroutine


!##########################################################################!
!**********************  start send_recv_real_data  ***********************!
!##########################################################################!
!##########################################################################!
!* This subroutine is usefull for send or/and receve double precision data to/from    *!
!* others processes                                                       *!
!##########################################################################!

!     subroutine exchange_struct_3d_d(struct, size1,size2, mode, srequests, rrequests, data_send, data_receive)
!        
!        implicit none
!
!        ! arguments 
!        integer size1,size2, mode
!        double precision struct(size1,size2)
!
!        !common
!        integer itanf,itend,idt,nits,niter,it
!        common /femtim/ itanf,itend,idt,nits,niter,it
!
!        ! local
!        integer i,j,sendElements,receiveElements, comm_loop,k, tag, newtag
!        double precision, dimension(:,:,:),allocatable :: data_send, data_receive 
!
!        ! MPI
!        integer ierr, data_type
!        integer, dimension(:) :: srequests,rrequests
!
!        comm_loop = domain%exchanges
!        data_type = MPI_DOUBLE_PRECISION
!
!        newtag = next_mpi_tag()
!
!        if(.not. allocated(data_send)) then
!          allocate(data_send(size1,domain%maxItems,comm_loop))
!          allocate(data_receive(size1,domain%maxItems,comm_loop))
!        end if
!        data_receive = 0.
!
!        if(mode .eq. 2) then
!
!           do i = 1, comm_loop
!
!              sendElements = domain%ele_send(i)%numItems
!              receiveElements = domain%ele_receive(i)%numItems
!
!              do k=1,sendElements
!                 do j=1,size1
!                    data_send(j,k,i) = struct(j,domain%ele_send(i)%items(k))
!                 end do
!              end do
!
!              srequests(i) = MPI_REQUEST_NULL
!              rrequests(i) = MPI_REQUEST_NULL
!
!              tag = newtag + domain%process(i) + my_id 
!
!              ! Non-blocking sends
!              call MPI_Isend(data_send(1,1,i),sendElements*size1,data_type, &
!                domain%process(i),tag, MPI_COMM_WORLD,srequests(i),ierr)
!              if(ierr .ne. MPI_SUCCESS) then
!                 write(6,*)'ierror:',ierr
!              end if
!
!              tag = newtag + my_id + domain%process(i) 
!              ! Non-blocking receives
!              call MPI_Irecv(data_receive(1,1,i),receiveElements*size1,data_type, &
!                domain%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
!              if(ierr .ne. MPI_SUCCESS) then
!                 write(6,*)'ierror:',ierr
!              end if
!
!           end do
!
!        else if(mode .eq. 3) then
!
!           do i = 1, comm_loop
!
!              sendElements = domain%node_send(i)%numItems
!              receiveElements = domain%node_send(i)%numItems 
!
!              do k=1,sendElements
!                 do j=1,size1
!                    data_send(j,k,i) = struct(j,domain%node_send(i)%items(k))
!                 end do
!              end do
!
!
!              srequests(i) = MPI_REQUEST_NULL
!              rrequests(i) = MPI_REQUEST_NULL
!
!              tag = newtag + domain%process(i) + my_id 
!
!              ! Non-blocking sends
!              call MPI_Isend(data_send(1,1,i),sendElements*size1,data_type, &
!                domain%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
!              if(ierr .ne. MPI_SUCCESS) then
!                 write(6,*)'ierror:',ierr
!              end if
!
!              tag = newtag + my_id + domain%process(i) 
!              ! Non-blocking receives
!              call MPI_Irecv(data_receive(1,1,i),receiveElements*size1,data_type, &
!                domain%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
!              if(ierr .ne. MPI_SUCCESS) then
!                 write(6,*)'ierror:',ierr
!              end if
!
!           end do
!
!        else if(mode .eq. 4) then
!
!           do i = 1, comm_loop
!
!              sendElements = domain%halo_send_e(i)%numItems
!              receiveElements = domain%halo_recv_e(i)%numItems
!
!              do k=1,sendElements
!                 do j=1,size1
!                    data_send(j,k,i) = struct(j,domain%halo_send_e(i)%items(k))
!                 end do
!              end do
!
!              srequests(i) = MPI_REQUEST_NULL
!              rrequests(i) = MPI_REQUEST_NULL
!
!              tag = newtag + domain%process(i) + my_id 
!
!              ! Non-blocking sends
!              call MPI_Isend(data_send(1,1,i),sendElements*size1,data_type, &
!                domain%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
!              if(ierr .ne. MPI_SUCCESS) then
!                 write(6,*)'ierror:',ierr
!              end if
!
!              tag = newtag + my_id + domain%process(i) 
!              ! Non-blocking receives
!              call MPI_Irecv(data_receive(1,1,i),receiveElements*size1,data_type, &
!                domain%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
!              if(ierr .ne. MPI_SUCCESS) then
!                 write(6,*)'ierror:',ierr
!              end if
!
!           end do
!
!        else if(mode .eq. 5) then
!
!           do i = 1, comm_loop
!
!              sendElements = domain%halo_send_k(i)%numItems
!              receiveElements = domain%halo_recv_k(i)%numItems
!
!              do k=1,sendElements
!                 do j=1,size1
!                    data_send(j,k,i) = struct(j,domain%halo_send_k(i)%items(k))
!                 end do
!              end do
!
!              srequests(i) = MPI_REQUEST_NULL
!              rrequests(i) = MPI_REQUEST_NULL
!
!              tag = newtag + domain%process(i) + my_id 
!
!              ! Non-blocking sends
!              call MPI_Isend(data_send(1,1,i),sendElements*size1,data_type, &
!                domain%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
!              if(ierr .ne. MPI_SUCCESS) then
!                 write(6,*)'ierror:',ierr
!              end if
!
!              tag = newtag + my_id + domain%process(i) 
!              ! Non-blocking receives
!              call MPI_Irecv(data_receive(1,1,i),receiveElements*size1,data_type, &
!                domain%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
!              if(ierr .ne. MPI_SUCCESS) then
!                 write(6,*)'ierror:',ierr
!              end if
!
!           end do
!
!        else
!           write(6,*)'error in send_recv_real_data', mode
!           stop
!        end if
!
!        return
!
!
!    end subroutine


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
        double precision time1

        if(ln_timing) time1 = MPI_Wtime()

        comm_loop = domain%exchanges
        data_type = MPI_DOUBLE_PRECISION

        newtag = next_mpi_tag()

        if(.not. allocated(data_send)) then
          allocate(data_send(size1,domain%maxItems,comm_loop))
          allocate(data_receive(size1,domain%maxItems,comm_loop))
        end if
        data_receive = 0.

        if(mode .eq. 2) then

           do i = 1, comm_loop

              sendElements = domain%ele_send(i)%numItems
              receiveElements = domain%ele_receive(i)%numItems

              do k=1,sendElements
                 do j=1,size1
                    data_send(j,k,i) = struct(j,domain%ele_send(i)%items(k))
                 end do
              end do

              srequests(i) = MPI_REQUEST_NULL
              rrequests(i) = MPI_REQUEST_NULL

              tag = newtag + domain%process(i) + my_id 

              ! Non-blocking sends
              call MPI_Isend(data_send(1,1,i),sendElements*size1,data_type, &
                 domain%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              tag = newtag + my_id + domain%process(i) 
              ! Non-blocking receives
              call MPI_Irecv(data_receive(1,1,i),receiveElements*size1,data_type, &
                domain%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

           end do

        else if(mode .eq. 3) then

           do i = 1, comm_loop

              sendElements = domain%node_send(i)%numItems
              receiveElements = domain%node_send(i)%numItems 

              do k=1,sendElements
                 do j=1,size1
                    data_send(j,k,i) = struct(j,domain%node_send(i)%items(k))
                 end do
              end do


              srequests(i) = MPI_REQUEST_NULL
              rrequests(i) = MPI_REQUEST_NULL

              tag = newtag + domain%process(i) + my_id 

              ! Non-blocking sends
              call MPI_Isend(data_send(1,1,i), sendElements*size1,data_type, &
                domain%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              tag = newtag + my_id + domain%process(i) 
              ! Non-blocking receives
              call MPI_Irecv(data_receive(1,1,i),receiveElements*size1,data_type, &
                domain%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

           end do

        else
           write(6,*)'error in send_recv_dp_data', mode
           stop
        end if

        if(ln_timing) comm_time = comm_time + MPI_Wtime() - time1

        return


    end subroutine


!##########################################################################!
!**********************  start send_recv_real3d_data  *********************!
!##########################################################################!
!##########################################################################!
!* This subroutine is usefull for send or/and receve double precision data to/from    *!
!* others processes                                                       *!
!##########################################################################!

     subroutine send_recv_real3d_data(struct,size1,size2,mode,srequests,rrequests,data_send,data_receive)

        implicit none

        ! arguments 
        integer size1,size2,mode
        double precision struct(size2,size1,1)

        !common
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

        ! local
        integer i,j,sendElements,receiveElements, comm_loop,k,n, tag, newtag
        double precision, dimension(:,:,:,:),allocatable :: data_send, data_receive 

        ! MPI
        integer ierr, data_type
        integer, dimension(:) :: srequests,rrequests
        double precision time1

        if(ln_timing) time1 = MPI_Wtime()


        comm_loop = domain%exchanges
        data_type = MPI_DOUBLE_PRECISION

        newtag = next_mpi_tag()

        if(.not. allocated(data_send)) then
          allocate(data_send(size2,size1,domain%maxItems,comm_loop))
          allocate(data_receive(size2,size1,domain%maxItems,comm_loop))
        end if
        data_receive = 0.

        if(mode .eq. 2) then

           do i = 1, comm_loop

              sendElements = domain%ele_send(i)%numItems
              receiveElements = domain%ele_receive(i)%numItems

              do k=1,sendElements
                 do j=1,size1
                    do n=1,size2
                       data_send(n,j,k,i) = struct(n,j,domain%ele_send(i)%items(k))
                    end do
                 end do
              end do

              srequests(i) = MPI_REQUEST_NULL
              rrequests(i) = MPI_REQUEST_NULL

              tag = newtag + domain%process(i) + my_id 

              ! Non-blocking sends
              call MPI_Isend(data_send(1,1,1,i),sendElements*size1*size2,data_type, &
                 domain%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              tag = newtag + my_id + domain%process(i) 
              ! Non-blocking receives
              call MPI_Irecv(data_receive(1,1,1,i),receiveElements*size1*size2,data_type, &
                 domain%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

           end do

        else if(mode .eq. 3) then

           do i = 1, comm_loop

              sendElements = domain%node_send(i)%numItems
              receiveElements = domain%node_send(i)%numItems 

              do k=1,sendElements
                 do j=1,size1
                    do n=1,size2
                       data_send(n,j,k,i) = struct(n,j,domain%node_send(i)%items(k))
                    end do
                 end do
              end do

              srequests(i) = MPI_REQUEST_NULL
              rrequests(i) = MPI_REQUEST_NULL

              tag = newtag + domain%process(i) + my_id 

              ! Non-blocking sends
              call MPI_Isend(data_send(1,1,1,i),sendElements*size1*size2,data_type, &
                domain%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              tag = newtag + my_id + domain%process(i) 
              ! Non-blocking receives
              call MPI_Irecv(data_receive(1,1,1,i),receiveElements*size1*size2, &
                 data_type,domain%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

           end do


        else
           write(6,*)'error in send_recv_real3d_data', mode
           stop
        end if

        if(ln_timing) comm_time = comm_time + MPI_Wtime() - time1

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
        double precision struct(size1,1)

        ! local
!        integer i,j,sendElements,receiveElements, comm_loop,k, tag
!        double precision, dimension(:,:),allocatable :: data_send, data_receive 

        ! MPI
        integer, dimension(:) :: requests
        integer ierr
        integer counter,n
        double precision time1

        if(ln_timing) time1 = MPI_Wtime()

        !integer, dimension(:), allocatable :: requests, statuses

        ! Wait for all non-blocking communications to complete
        allocate(statuses(MPI_STATUS_SIZE, size(requests)))

        call mpi_waitall(size(requests), requests, statuses, ierr)
        if(ierr .ne. MPI_SUCCESS) then
           call abort
        end if
        do n=0,domain%exchanges-1
           write(6,*)'count myrank',n,statuses(MPI_STATUS_SIZE,n+1),my_id
           write(6,*)'cancelled myrank',n,statuses(MPI_STATUS_SIZE,n+2),my_id
           write(6,*)'sender myrank',n,statuses(MPI_STATUS_SIZE,n+3),my_id
           write(6,*)'tag myrank',n,statuses(MPI_STATUS_SIZE,n+4),my_id
           write(6,*)'error myrank',n,statuses(MPI_STATUS_SIZE,n+5),my_id
        end do
        deallocate(statuses)

        if(ln_timing) comm_time = comm_time + MPI_Wtime() - time1

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

	include 'femtime.h'

        ! argument
        integer, dimension(:) :: requests

        !local
        integer counter,ierr,i
        integer status(MPI_STATUS_SIZE),rindex
	integer,save :: icall
	data icall /0/
        double precision time1

        if(ln_timing) time1 = MPI_Wtime()

	!write(6,*),'icall waitany:',icall
	icall=icall+1
        !write(6,*)'size_requests',size(requests),requests

        do i=1,size(requests)
           call MPI_WaitAny(size(requests), requests, rindex, status,ierr)
           if(ierr .ne. MPI_SUCCESS) then
              call abort
           end if
        end do

        if(ln_timing) comm_time = comm_time + MPI_Wtime() - time1

        return

        end subroutine


!##########################################################################!
!**********************  start send_recv_test  ***********************!
!##########################################################################!
!##########################################################################!
!* This subroutine is usefull for send or/and receve double precision data to/from    *!
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
        !double precision, dimension(:),allocatable :: data_send, data_receive 
        double precision, dimension(:) :: data_receive 
        double precision struct

        ! MPI
        integer ierr, data_type
        double precision time1

        if(ln_timing) time1 = MPI_Wtime()

        struct = 1

        comm_loop = domain%exchanges
        data_type = MPI_DOUBLE_PRECISION

        newtag = next_mpi_tag()

           do i = 1, comm_loop

              sendElements = 1 
              receiveElements = 1

              srequests(i) = MPI_REQUEST_NULL
              rrequests(i) = MPI_REQUEST_NULL

              tag = newtag + domain%process(i)! + my_id 

              ! Non-blocking sends
              call MPI_Isend(struct,sendElements,data_type, &
                domain%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              write(6,*)'isend',sendElements,domain%process(i),tag,my_id
 
              tag = newtag + my_id
              ! Non-blocking receives
              call MPI_Irecv(data_receive(i),receiveElements,data_type, &
                domain%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
              if(ierr .ne. MPI_SUCCESS) then
                 write(6,*)'ierror:',ierr
              end if

              write(6,*)'irecv',receiveElements,domain%process(i),tag,my_id

           end do

        if(ln_timing) comm_time = comm_time + MPI_Wtime() - time1

        return

    end subroutine

!================================================================

    subroutine rebuild_struct_2d_max_i(struct,size1,data_receive,mode)


        implicit none

        !arguments
        integer mode,size1
        integer struct(size1)
        integer, dimension(:,:) :: data_receive

        !local
        integer receiveElements,i,j,k
        double precision time1

        if(ln_timing) time1 = MPI_Wtime()

        if( mode .eq. 2) then

           do i=1,domain%exchanges

           receiveElements = domain%ele_receive(i)%numItems

              do k=1,receiveElements
                 struct(domain%ele_receive(i)%items(k)) = &
                 max(struct(domain%ele_receive(i)%items(k)),data_receive(k,i))
              end do
           end do

        else if( mode .eq. 3) then

           do i=1,domain%exchanges

           receiveElements = domain%node_send(i)%numItems

              do k=1,receiveElements
                 struct(domain%node_send(i)%items(k)) = max(struct(domain%node_send(i)%items(k)),data_receive(k,i))
              end do
           end do

        else
           write(6,*)'error in rebuild_int_max1, mode=',mode
        end if

        if(ln_timing) comm_time = comm_time + MPI_Wtime() - time1
        

        return


    end subroutine

!================================================================

    subroutine rebuild_struct_2d_max_r(struct,size1,data_receive,mode)


        implicit none

        !arguments
        integer mode,size1
        double precision struct(size1)
        double precision, dimension(:,:) :: data_receive

        !local
        integer receiveElements,i,j,k
        double precision time1

        if(ln_timing) time1 = MPI_Wtime()

        if( mode .eq. 2) then

           do i=1,domain%exchanges

           receiveElements = domain%ele_receive(i)%numItems

              do k=1,receiveElements
                 struct(domain%ele_receive(i)%items(k)) = &
                max(struct(domain%ele_receive(i)%items(k)),data_receive(k,i))
              end do
           end do

        else if( mode .eq. 3) then

           do i=1,domain%exchanges

           receiveElements = domain%node_send(i)%numItems

              do k=1,receiveElements
                 struct(domain%node_send(i)%items(k)) = max(struct(domain%node_send(i)%items(k)), data_receive(k,i))
              end do
           end do

        else
           write(6,*)'error in rebuild_int_max1, mode=',mode
        end if

        if(ln_timing) comm_time = comm_time + MPI_Wtime() - time1

        return


    end subroutine

!================================================================

    subroutine rebuild_struct_2d_min_i(struct,size1,data_receive,mode)


        implicit none

        !arguments
        integer mode,size1
        integer struct(size1)
        integer, dimension(:,:) :: data_receive

        !local
        integer receiveElements,i,j,k
        double precision time1

        if(ln_timing) time1 = MPI_Wtime()

        if( mode .eq. 2) then

           do i=1,domain%exchanges

           receiveElements = domain%ele_receive(i)%numItems

              do k=1,receiveElements
                 struct(domain%ele_receive(i)%items(k)) = &
                min(struct(domain%ele_receive(i)%items(k)),data_receive(k,i))
              end do
           end do

        else if( mode .eq. 3) then

           do i=1,domain%exchanges

           receiveElements = domain%node_send(i)%numItems

              do k=1,receiveElements
                 struct(domain%node_send(i)%items(k)) = &
                min(struct(domain%node_send(i)%items(k)), data_receive(k,i))
              end do
           end do

        else
           write(6,*)'error in rebuild_int_min1, mode=',mode
        end if
        
        if(ln_timing) comm_time = comm_time + MPI_Wtime() - time1

        return


    end subroutine

!================================================================

    subroutine rebuild_struct_2d_min_r(struct,size1,data_receive,mode)


        implicit none

        !arguments
        integer mode,size1
        double precision struct(size1)
        double precision, dimension(:,:) :: data_receive

        !local
        integer receiveElements,i,j,k
        double precision time1

        if(ln_timing) time1 = MPI_Wtime()

        if( mode .eq. 2) then

           do i=1,domain%exchanges

           receiveElements = domain%ele_receive(i)%numItems

              do k=1,receiveElements
                 struct(domain%ele_receive(i)%items(k)) = &
                min(struct(domain%ele_receive(i)%items(k)),data_receive(k,i))
              end do
           end do

        else if( mode .eq. 3) then

           do i=1,domain%exchanges

           receiveElements = domain%node_send(i)%numItems

              do k=1,receiveElements
                 struct(domain%node_send(i)%items(k)) = &
                min(struct(domain%node_send(i)%items(k)), data_receive(k,i))
              end do
           end do

        else
           write(6,*)'error in rebuild_int_min1, mode=',mode
        end if
        
        if(ln_timing) comm_time = comm_time + MPI_Wtime() - time1

        return


    end subroutine

!================================================================

    subroutine rebuild_struct_2d_sum_r(struct,size1,data_receive,mode)


        implicit none

        !arguments
        integer mode,size1
        double precision struct(size1)
        double precision, dimension(:,:) :: data_receive

        !local
        integer receiveElements,i,j,k
        double precision time1

        if(ln_timing) time1 = MPI_Wtime()

        if( mode .eq. 2) then

           do i=1,domain%exchanges

              receiveElements = domain%ele_receive(i)%numItems

              do k=1,receiveElements
                 struct(domain%ele_receive(i)%items(k)) = &
                struct(domain%ele_receive(i)%items(k)) + data_receive(k,i)
              end do
           end do

        else if( mode .eq. 3) then

           do i=1,domain%exchanges

              receiveElements = domain%node_send(i)%numItems

              do k=1,receiveElements
                 struct(domain%node_send(i)%items(k)) = &
                struct(domain%node_send(i)%items(k)) + data_receive(k,i)
              end do
           end do

        else
           write(6,*)'error in Rebuild_struct_2d_sum, mode=',mode
        end if
        
        if(ln_timing) comm_time = comm_time + MPI_Wtime() - time1

        return


    end subroutine

!================================================================

!    subroutine rebuild_struct_2d_sum_d(struct,size1,data_receive,mode)
!
!
!        implicit none
!
!        !arguments
!        integer mode,size1
!        double precision struct(size1)
!        double precision, dimension(:,:) :: data_receive
!
!        !local
!        integer receiveElements,i,j,k
!
!        if( mode .eq. 2) then
!
!           do i=1,domain%exchanges
!
!           receiveElements = domain%ele_receive(i)%numItems
!
!              do k=1,receiveElements
!                 struct(domain%ele_receive(i)%items(k)) = &
!                struct(domain%ele_receive(i)%items(k)) + data_receive(k,i)
!              end do
!           end do
!
!        else if( mode .eq. 3) then
!
!           do i=1,domain%exchanges
!
!           receiveElements = domain%node_send(i)%numItems
!
!              do k=1,receiveElements
!                 struct(domain%node_send(i)%items(k)) = &
!                struct(domain%node_send(i)%items(k)) + data_receive(k,i)
!              end do
!           end do
!
!        else
!           write(6,*)'error in rebuild_int_min1, mode=',mode
!        end if
!        
!
!        return
!
!
!    end subroutine



!================================================================

    subroutine rebuild_struct_2d_sum_i(struct,size1,data_receive,mode)


        implicit none

        !arguments
        integer mode,size1
        integer struct(size1)
        integer, dimension(:,:) :: data_receive

        !local
        integer receiveElements,i,j,k
        double precision time1

        if(ln_timing) time1 = MPI_Wtime()

        if( mode .eq. 2) then

           do i=1,domain%exchanges

           receiveElements = domain%ele_receive(i)%numItems

              do k=1,receiveElements
                 struct(domain%ele_receive(i)%items(k)) = &
                struct(domain%ele_receive(i)%items(k)) + data_receive(k,i)
              end do
           end do

        else if( mode .eq. 3) then

           do i=1,domain%exchanges

           receiveElements = domain%node_send(i)%numItems

              do k=1,receiveElements
                 struct(domain%node_send(i)%items(k)) = &
                struct(domain%node_send(i)%items(k)) + data_receive(k,i)
              end do
           end do

        else
           write(6,*)'error in rebuild_int_min1, mode=',mode
        end if
        
        if(ln_timing) comm_time = comm_time + MPI_Wtime() - time1

        return


    end subroutine

!================================================================

    subroutine rebuild_struct_3d_sum_r(struct,size1,size2,data_receive,mode)


        implicit none

        !arguments
        integer size1,size2,mode
        double precision struct(size1,size2)
        double precision, dimension(:,:,:) :: data_receive

        !local
        integer receiveElements,i,j,k
        double precision time1

        if(ln_timing) time1 = MPI_Wtime()

        if( mode .eq. 2) then

           do i=1,domain%exchanges

           receiveElements = domain%ele_receive(i)%numItems

              do k=1,receiveElements
                 do j=1,size1
                    struct(j,domain%ele_receive(i)%items(k)) = &
                struct(j,domain%ele_receive(i)%items(k)) + data_receive(j,k,i)
                 end do
              end do
           end do

        else if( mode .eq. 3) then

           do i=1,domain%exchanges

           receiveElements = domain%node_send(i)%numItems

              do k=1,receiveElements
                 do j=1,size1
                    struct(j,domain%node_send(i)%items(k)) = &
                struct(j,domain%node_send(i)%items(k)) + data_receive(j,k,i)
                 end do
              end do
           end do

        else
           write(6,*)'error in rebuild_RealHalo, mode=',mode
        end if
        
        if(ln_timing) comm_time = comm_time + MPI_Wtime() - time1

        return


    end subroutine

!================================================================

    subroutine rebuild_struct_3d_i(struct,size1,size2,data_receive,mode)


        implicit none

        !arguments
        integer size1,size2,mode
        integer struct(size1,size2)
        integer, dimension(:,:,:) :: data_receive

        !local
        integer receiveElements,i,j,k
        double precision time1

        if(ln_timing) time1 = MPI_Wtime()

        if( mode .eq. 2) then

           do i=1,domain%exchanges

           receiveElements = domain%halo_recv_e(i)%numItems

              do k=1,receiveElements
                 do j=1,size1
                    struct(j,domain%halo_recv_e(i)%items(k)) = &
                        data_receive(j,k,i)
                 end do
              end do
           end do

        else if( mode .eq. 3) then

           do i=1,domain%exchanges

           receiveElements = domain%halo_recv_k(i)%numItems

              do k=1,receiveElements
                 do j=1,size1
                    struct(j,domain%halo_recv_k(i)%items(k)) = &
                        data_receive(j,k,i)
                 end do
              end do
           end do

        else
           write(6,*)'error in rebuild_RealHalo, mode=',mode
        end if
        
        if(ln_timing) comm_time = comm_time + MPI_Wtime() - time1

        return


    end subroutine

!================================================================

    subroutine rebuild_struct_3d_r(struct,size1,size2,data_receive,mode)


        implicit none

        !arguments
        integer size1,size2,mode
        double precision struct(size1,size2)
        double precision, dimension(:,:,:) :: data_receive

        !local
        integer receiveElements,i,j,k
        double precision time1

        if(ln_timing) time1 = MPI_Wtime()

        if( mode .eq. 2) then

           do i=1,domain%exchanges

           receiveElements = domain%halo_recv_e(i)%numItems

              do k=1,receiveElements
                 do j=1,size1
                    struct(j,domain%halo_recv_e(i)%items(k)) = &
                        data_receive(j,k,i)
                 end do
              end do
           end do

        else if( mode .eq. 3) then

           do i=1,domain%exchanges

           receiveElements = domain%halo_recv_k(i)%numItems

              do k=1,receiveElements
                 do j=1,size1
                    struct(j,domain%halo_recv_k(i)%items(k)) = &
                        data_receive(j,k,i)
                 end do
              end do
           end do

        else
           write(6,*)'error in rebuild_RealHalo, mode=',mode
        end if
        
        if(ln_timing) comm_time = comm_time + MPI_Wtime() - time1

        return


    end subroutine

!================================================================

!    subroutine rebuild_struct_3d_d(struct,size1,size2,data_receive,mode)
!
!
!        implicit none
!
!        !arguments
!        integer size1,size2,mode
!        double precision struct(size1,size2)
!        double precision, dimension(:,:,:) :: data_receive
!
!        !local
!        integer receiveElements,i,j,k
!
!        if( mode .eq. 2) then
!
!           do i=1,domain%exchanges
!
!           receiveElements = domain%halo_recv_e(i)%numItems
!
!              do k=1,receiveElements
!                 do j=1,size1
!                    struct(j,domain%halo_recv_e(i)%items(k)) = &
!                        data_receive(j,k,i)
!                 end do
!              end do
!           end do
!
!        else if( mode .eq. 3) then
!
!           do i=1,domain%exchanges
!
!           receiveElements = domain%halo_recv_k(i)%numItems
!
!              do k=1,receiveElements
!                 do j=1,size1
!                    struct(j,domain%halo_recv_k(i)%items(k)) = &
!                        data_receive(j,k,i)
!                 end do
!              end do
!           end do
!
!        else
!           write(6,*)'error in rebuild_RealHalo, mode=',mode
!        end if
!        
!
!        return
!
!
!    end subroutine

!================================================================

    subroutine rebuild_struct_2d_i(struct,size1,data_receive,mode)


        implicit none

        !arguments
        integer size1,mode
        integer struct(size1)
        integer, dimension(:,:) :: data_receive

        !local
        integer receiveElements,i,j,k
        double precision time1

        if(ln_timing) time1 = MPI_Wtime()

        if( mode .eq. 2) then

           do i=1,domain%exchanges

           receiveElements = domain%halo_recv_e(i)%numItems

              do k=1,receiveElements
                struct(domain%halo_recv_e(i)%items(k)) = &
                        data_receive(k,i)
              end do
           end do

        else if( mode .eq. 3) then

           do i=1,domain%exchanges

           receiveElements = domain%halo_recv_k(i)%numItems

              do k=1,receiveElements
                struct(domain%halo_recv_k(i)%items(k)) = &
                        data_receive(k,i)
              end do
           end do

        else
           write(6,*)'error in rebuild_RealHalo, mode=',mode
        end if
        
        if(ln_timing) comm_time = comm_time + MPI_Wtime() - time1

        return


    end subroutine

!================================================================

    subroutine rebuild_struct_2d_r(struct,size1,data_receive,mode)


        implicit none

        !arguments
        integer size1,mode
        double precision struct(size1)
        double precision, dimension(:,:) :: data_receive

        !local
        integer receiveElements,i,j,k
        double precision time1

        if(ln_timing) time1 = MPI_Wtime()

        if( mode .eq. 2) then

           do i=1,domain%exchanges

           receiveElements = domain%halo_recv_e(i)%numItems

              do k=1,receiveElements
                struct(domain%halo_recv_e(i)%items(k)) = &
                        data_receive(k,i)
              end do
           end do

        else if( mode .eq. 3) then

           do i=1,domain%exchanges

           receiveElements = domain%halo_recv_k(i)%numItems

              do k=1,receiveElements
                struct(domain%halo_recv_k(i)%items(k)) = &
                        data_receive(k,i)
              end do
           end do

        else
           write(6,*)'error in rebuild_RealHalo, mode=',mode
        end if
        
        if(ln_timing) comm_time = comm_time + MPI_Wtime() - time1

        return


    end subroutine

!================================================================

!    subroutine rebuild_struct_2d_d(struct,size1,data_receive,mode)
!
!
!        implicit none
!
!        !arguments
!        integer size1,mode
!        double precision struct(size1)
!        double precision, dimension(:,:) :: data_receive
!
!        !local
!        integer receiveElements,i,j,k
!
!        if( mode .eq. 2) then
!
!           do i=1,domain%exchanges
!
!           receiveElements = domain%halo_recv_e(i)%numItems
!
!              do k=1,receiveElements
!                struct(domain%halo_recv_e(i)%items(k)) = &
!                        data_receive(k,i)
!              end do
!           end do
!
!        else if( mode .eq. 3) then
!
!           do i=1,domain%exchanges
!
!           receiveElements = domain%halo_recv_k(i)%numItems
!
!              do k=1,receiveElements
!                struct(domain%halo_recv_k(i)%items(k)) = &
!                        data_receive(k,i)
!              end do
!           end do
!
!        else
!           write(6,*)'error in rebuild_RealHalo, mode=',mode
!        end if
!        
!
!        return
!
!
!    end subroutine
!
!================================================================

!    subroutine rebuild_struct_3d_sum_d(struct,size1,size2,data_receive,mode)
!
!
!        implicit none
!
!        !arguments
!        integer size1,size2,mode
!        double precision struct(size1,size2)
!        double precision, dimension(:,:,:) :: data_receive
!
!        !local
!        integer receiveElements,i,j,k
!
!        if( mode .eq. 2) then
!
!           do i=1,domain%exchanges
!
!           receiveElements = domain%ele_receive(i)%numItems
!
!              do k=1,receiveElements
!                 do j=1,size1
!                    struct(j,domain%ele_receive(i)%items(k)) = &
!                struct(j,domain%ele_receive(i)%items(k)) + data_receive(j,k,i)
!                 end do
!              end do
!           end do
!
!        else if( mode .eq. 3) then
!
!           do i=1,domain%exchanges
!
!           receiveElements = domain%node_send(i)%numItems
!
!              do k=1,receiveElements
!                 do j=1,size1
!                    struct(j,domain%node_send(i)%items(k)) = &
!                struct(j,domain%node_send(i)%items(k)) + data_receive(j,k,i)
!                 end do
!              end do
!           end do
!
!        else
!           write(6,*)'error in rebuild_RealHalo, mode=',mode
!        end if
!        
!
!        return
!
!
!    end subroutine

!================================================================

    subroutine rebuild_struct_3d_2xSum_r(struct,size1,size2,data_receive,mode)


        implicit none

        !arguments
        integer size1,size2,mode
        double precision struct(size1,size2)
        double precision, dimension(:,:,:) :: data_receive

        !local
        integer receiveElements,i,j,k
        double precision time1

        if(ln_timing) time1 = MPI_Wtime()

        if( mode .eq. 2) then

           do i=1,domain%exchanges

           receiveElements = domain%ele_receive(i)%numItems

              do k=1,receiveElements
                 do j=1,size1
                    struct(j,domain%ele_receive(i)%items(k)) = &
                2 * (struct(j,domain%ele_receive(i)%items(k)) + data_receive(j,k,i))
                 end do
              end do
           end do

        else if( mode .eq. 3) then

           do i=1,domain%exchanges

           receiveElements = domain%node_send(i)%numItems

              do k=1,receiveElements
                 do j=1,size1
                    struct(j,domain%node_send(i)%items(k)) = &
                2 * (struct(j,domain%node_send(i)%items(k)) + data_receive(j,k,i))
                 end do
              end do
           end do

        else
           write(6,*)'error in rebuild_RealHalo, mode=',mode
        end if
        
        if(ln_timing) comm_time = comm_time + MPI_Wtime() - time1

        return


    end subroutine

!================================================================

    subroutine rebuild_struct_3d_2xSum_d(struct,size1,size2,data_receive,mode)


        implicit none

        !arguments
        integer size1,size2,mode
        double precision struct(size1,size2)
        double precision, dimension(:,:,:) :: data_receive

        !local
        integer receiveElements,i,j,k
        double precision time1

        if(ln_timing) time1 = MPI_Wtime()

        if( mode .eq. 2) then

           do i=1,domain%exchanges

           receiveElements = domain%ele_receive(i)%numItems

              do k=1,receiveElements
                 do j=1,size1
                    struct(j,domain%ele_receive(i)%items(k)) = &
                2 * (struct(j,domain%ele_receive(i)%items(k)) + data_receive(j,k,i))
                 end do
              end do
           end do

        else if( mode .eq. 3) then

           do i=1,domain%exchanges

           receiveElements = domain%node_send(i)%numItems

              do k=1,receiveElements
                 do j=1,size1
                    struct(j,domain%node_send(i)%items(k)) = &
                2 * (struct(j,domain%node_send(i)%items(k)) + data_receive(j,k,i))
                 end do
              end do
           end do

        else
           write(6,*)'error in rebuild_RealHalo, mode=',mode
        end if
        
        if(ln_timing) comm_time = comm_time + MPI_Wtime() - time1

        return


    end subroutine

!******************************************************************

    subroutine rebuild_DPHalo(struct,size1,size2,data_receive,mode)


        implicit none

        !arguments
        integer size1,size2,mode
        double precision struct(size1,size2)
        double precision, dimension(:,:,:) :: data_receive

        !local
        integer receiveElements,i,j,k
        double precision time1

        if(ln_timing) time1 = MPI_Wtime()

        if( mode .eq. 2) then

           do i=1,domain%exchanges

           receiveElements = domain%ele_receive(i)%numItems

              do k=1,receiveElements
                 do j=1,size1
                    struct(j,domain%ele_receive(i)%items(k)) =  data_receive(j,k,i)
                 end do
              end do
           end do

        else if( mode .eq. 3) then

           do i=1,domain%exchanges

           receiveElements = domain%node_send(i)%numItems

              do k=1,receiveElements
                 do j=1,size1
                    struct(j,domain%node_send(i)%items(k)) = &
                struct(j,domain%node_send(i)%items(k)) + data_receive(j,k,i)
                 end do
              end do
           end do

        else
           write(6,*)'error in rebuild_DPHalo, mode=',mode
        end if
        
        if(ln_timing) comm_time = comm_time + MPI_Wtime() - time1

        return


    end subroutine


!******************************************************************

    subroutine rebuild_RealHalo(struct,size1,data_receive,mode)


        implicit none

        !arguments
        integer size1,mode
        double precision struct(size1,1)
        double precision, dimension(:,:,:) :: data_receive

        !local
        integer receiveElements,i,j,k
        double precision time1

        if(ln_timing) time1 = MPI_Wtime()

        if( mode .eq. 2) then

           do i=1,domain%exchanges

           receiveElements = domain%ele_receive(i)%numItems

              do k=1,receiveElements
                 do j=1,size1
                    struct(j,domain%ele_receive(i)%items(k)) =  data_receive(j,k,i)
                 end do
              end do
           end do

        else if( mode .eq. 3) then

           do i=1,domain%exchanges

           receiveElements = domain%node_send(i)%numItems

              do k=1,receiveElements
                 do j=1,size1
                    struct(j,domain%node_send(i)%items(k)) = &
                struct(j,domain%node_send(i)%items(k)) + data_receive(j,k,i)
                 end do
              end do
           end do

        else
           write(6,*)'error in rebuild_RealHalo, mode=',mode
        end if
        
        if(ln_timing) comm_time = comm_time + MPI_Wtime() - time1

        return


    end subroutine


    subroutine rebuild_Real3d_Halo(struct,size1,size2,data_receive,mode)


        implicit none

        !arguments
        integer size1,size2,mode
        double precision struct(size1,size2,1)
        double precision, dimension(:,:,:,:) :: data_receive

        !local
        integer receiveElements,i,j,k,n
        double precision time1

        if(ln_timing) time1 = MPI_Wtime()

        if( mode .eq. 2) then

           do i=1,domain%exchanges

           receiveElements = domain%ele_receive(i)%numItems

              do k=1,receiveElements
                 do j=1,size1
                    do n=1,size2
                       struct(n,j,domain%ele_receive(i)%items(k)) =  data_receive(n,j,k,i)
                    end do
                 end do
              end do
           end do

        else
           write(6,*)'error in rebuild_Real3d_Halo, mode=',mode
        end if
        
        if(ln_timing) comm_time = comm_time + MPI_Wtime() - time1

        return


    end subroutine

!######################################################################!
!*******************  start next_mpi_tag function  ********************!
!######################################################################!
!######################################################################!
!* This function computes the next tag for mpi communications *********!
!######################################################################!

   integer function next_mpi_tag()

      implicit none
      include 'femtime.h'
      integer, save::last_tag=0, tag_it=0

      if(tag_it .ne. niter) then
        tag_it = niter
        last_tag = 0
      end if
      last_tag = last_tag + n_threads + 1

      next_mpi_tag = last_tag

   end function next_mpi_tag

   !##########################################################################!

     subroutine send_coo(struct, srequests, rrequests, data_send, &
     &                 data_receive,index_coo,size1)
        
        use basin

        implicit none

        integer node,size1
        ! arguments 
        double precision struct(size1)
        integer index_coo(size1)

        ! local
        integer i,j,n,sendElements,receiveElements, comm_loop,k, tag, newtag
        double precision, dimension(:,:),allocatable :: data_send, data_receive 

        ! MPI
        integer ierr, data_type
        integer, dimension(:) :: srequests,rrequests
        double precision time1

        if(ln_timing) time1 = MPI_Wtime()

        comm_loop = domain%exchanges
        data_type = MPI_DOUBLE_PRECISION

        newtag = next_mpi_tag()

        if(.not. allocated(data_send)) then
          allocate(data_send(domain%maxItems,comm_loop))
          allocate(data_receive(domain%maxItems,comm_loop))
        end if
        data_receive = 0.

        do i = 1, comm_loop

           sendElements = domain%node_send(i)%numItems
           receiveElements = domain%node_send(i)%numItems 

           do k=1,sendElements
             node=domain%node_send(i)%items(k)
             do n=1,size1
               if(node.eq.index_coo(n).and.(struct(n).gt.0)) exit
             end do
             data_send(k,i) = struct(n)
           end do


           srequests(i) = MPI_REQUEST_NULL
           rrequests(i) = MPI_REQUEST_NULL

           tag = newtag + domain%process(i) + my_id 

           ! Non-blocking sends
           call MPI_Isend(data_send(1,i),sendElements,data_type, &
             domain%process(i),tag,MPI_COMM_WORLD,srequests(i),ierr)
           if(ierr .ne. MPI_SUCCESS) then
              write(6,*)'ierror:',ierr
           end if

           tag = newtag + my_id + domain%process(i) 
           ! Non-blocking receives
           call MPI_Irecv(data_receive(1,i),receiveElements,data_type, &
             domain%process(i),tag,MPI_COMM_WORLD,rrequests(i),ierr)
           if(ierr .ne. MPI_SUCCESS) then
              write(6,*)'ierror:',ierr
           end if

         end do

        if(ln_timing) comm_time = comm_time + MPI_Wtime() - time1

         return

      end subroutine

!================================================================

    subroutine recv_coo(struct,index_coo,size1,data_receive)


        implicit none

        !arguments
        integer size1,node
        double precision struct(size1)
        double precision temp(size1)
        integer index_coo(size1)
        double precision, dimension(:,:) :: data_receive
        
        !local
        integer receiveElements,i,j,k,n
        double precision time1

        if(ln_timing) time1 = MPI_Wtime()
        !character*26 formato

        !formato='(A,E24.17,E24.17,I8,I8)'
        temp = 0.d0

         do i=1,domain%pos

            receiveElements = domain%node_send(i)%numItems
            
            do k=1,receiveElements
              node =domain%node_send(i)%items(k)
              do n=1,size1
                if(node.eq.index_coo(n).and.(struct(n).gt.0)) exit
              end do

                temp(n) = temp(n) + data_receive(k,i)

            end do
         end do

         temp(:) = temp(:) + struct(:)

         do i=domain%pos+1,domain%exchanges
            
            receiveElements = domain%node_send(i)%numItems
            
            do k=1,receiveElements
              node =domain%node_send(i)%items(k)
              do n=1,size1
                if(node.eq.index_coo(n).and.(struct(n).gt.0)) exit
              end do

                temp(n) = temp(n) + data_receive(k,i)

            end do
         end do

          struct = temp

!
!                if(((my_id.eq.12).and.(node.eq.5496)).or.((my_id.eq.13).and.(node.eq.4045)).or.((my_id.eq.14).and.(node.eq.5470))) then
!              write(6,formato)'part_ccoo:',struct(n),data_receive(k,i),my_id,domain%process(i)
!                end if
!

!                if(((my_id.eq.12).and.(node.eq.5496)).or.((my_id.eq.13).and.(node.eq.4045)).or.((my_id.eq.14).and.(node.eq.5470))) then
!              write(6,formato)'part_ccoo:',struct(n),data_receive(k,i),my_id,domain%process(i)
!                end if
!


        if(ln_timing) comm_time = comm_time + MPI_Wtime() - time1

        return


    end subroutine

!-------------------------------------------------------------------------------------
   end module mpi_communication
!-------------------------------------------------------------------------------------

