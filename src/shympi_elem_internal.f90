!
! mpi routines
!
! contents :
!
! revision log :
!
! 24.11.2015    ggu     project started
! 22.06.2016    ggu     added sum option to shympi_reduce
!
!******************************************************************

        module shympi_internal

        use mpi_common_struct
        use mpi_communication
        use timing 

        implicit none

        contains

!********************************

	subroutine shympi_init_internal(my_id,n_threads)

	implicit none

	integer my_id,n_threads

	integer ierr

        call MPI_INIT( ierr )
        call MPI_COMM_RANK( MPI_COMM_WORLD, my_id, ierr )
        call MPI_COMM_SIZE( MPI_COMM_WORLD, n_threads, ierr )

	end subroutine shympi_init_internal

!******************************************************************

	subroutine shympi_barrier_internal

	implicit none

	integer ierr

	call MPI_BARRIER( MPI_COMM_WORLD, ierr)

	end subroutine shympi_barrier_internal

!******************************************************************

        subroutine shympi_abort_internal

        implicit none

        integer ierr,ierr_code

        ierr_code = 33
	call MPI_ABORT(MPI_COMM_WORLD,ierr_code,ierr)

        end subroutine shympi_abort_internal

!******************************************************************

        subroutine shympi_finalize_internal

        implicit none

        integer ierr

	call MPI_FINALIZE(ierr)

        return

        end subroutine shympi_finalize_internal

!******************************************************************

        subroutine shympi_get_status_size_internal(size)

        implicit none

        integer size

	size = MPI_STATUS_SIZE

        end subroutine shympi_get_status_size_internal

!******************************************************************

	subroutine shympi_syncronize_internal

	implicit none

	integer ierr

	!flush(6)
	call MPI_BARRIER( MPI_COMM_WORLD, ierr)

	end subroutine shympi_syncronize_internal

!******************************************************************
!******************************************************************
!******************************************************************

        subroutine shympi_gather_i_internal(val)

	implicit none

        integer val

        integer ierr
        double precision time1

        if(ln_timing) time1 = MPI_WTIME()

        call MPI_GATHER (val,1,MPI_INT,ival,1,MPI_INT,0,MPI_COMM_WORLD,ierr) 

        if(ln_timing) comm_time = comm_time + MPI_WTIME() - time1

        end subroutine shympi_gather_i_internal

!*******************************

        subroutine shympi_bcast_i_internal(val,id)

	implicit none

        integer val
        integer, optional :: id

        integer ierr
        double precision time1

        if(ln_timing) time1 = MPI_WTIME()

        if(present(id)) then
          call MPI_BCAST(val,1,MPI_INT,id,MPI_COMM_WORLD,ierr)
        else
          call MPI_BCAST(val,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
        end if

        if(ln_timing) comm_time = comm_time + MPI_WTIME() - time1

        end subroutine shympi_bcast_i_internal

!*******************************

        subroutine shympi_bcast_r_internal(val,id)

	implicit none

        double precision val
        integer, optional :: id

        integer ierr
        double precision time1

        if(ln_timing) time1 = MPI_WTIME()

        if(present(id)) then
          call MPI_BCAST(val,1,MPI_DOUBLE_PRECISION,id,MPI_COMM_WORLD,ierr)
        else
          call MPI_BCAST(val,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        end if

        if(ln_timing) comm_time = comm_time + MPI_WTIME() - time1

        end subroutine shympi_bcast_r_internal

!*******************************

	subroutine shympi_reduce_r_internal(what,val)

	implicit none

	character*(*) what
	double precision val

        integer ierr
	double precision valout
        double precision time1

        if(ln_timing) time1 = MPI_WTIME()

        if( what == 'min' ) then
	  call MPI_ALLREDUCE(val,valout,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
	  val = valout
        else if( what == 'max' ) then
	  call MPI_ALLREDUCE(val,valout,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
	  val = valout
        else if( what == 'sum' ) then
	  call MPI_ALLREDUCE(val,valout,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
	  val = valout
        else
          write(6,*) 'what = ',what
          stop 'error stop shympi_reduce_r_internal: not ready'
        end if

        if(ln_timing) comm_time = comm_time + MPI_WTIME() - time1

	end subroutine shympi_reduce_r_internal

!*******************************

	subroutine shympi_reduce_d_internal(what,val)

	implicit none

	character*(*) what
	double precision val

        integer ierr
	double precision valout
        double precision time1

        if(ln_timing) time1 = MPI_WTIME()

        if( what == 'min' ) then
	  call MPI_ALLREDUCE(val,valout,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
	  val = valout
        else if( what == 'max' ) then
	  call MPI_ALLREDUCE(val,valout,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
	  val = valout
        else if( what == 'sum' ) then
	  call MPI_ALLREDUCE(val,valout,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
	  val = valout
        else
          write(6,*) 'what = ',what
          stop 'error stop shympi_reduce_r_internal: not ready'
        end if

        if(ln_timing) comm_time = comm_time + MPI_WTIME() - time1

	end subroutine shympi_reduce_d_internal

!******************************************************************

	subroutine shympi_reduce_i_internal(what,val)

	implicit none

	character*(*) what
	integer val

        integer ierr
	integer valout
        double precision time1

        if(ln_timing) time1 = MPI_WTIME()

        if( what == 'min' ) then
	  call MPI_ALLREDUCE(val,valout,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
	  val = valout
        else if( what == 'max' ) then
	  call MPI_ALLREDUCE(val,valout,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
	  val = valout
        else if( what == 'sum' ) then
	  call MPI_ALLREDUCE(val,valout,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
	  val = valout
        else
          write(6,*) 'what = ',what
          stop 'error stop shympi_reduce_i_internal: not ready'
        end if

        if(ln_timing) comm_time = comm_time + MPI_WTIME() - time1

	end subroutine shympi_reduce_i_internal

!*****************************************************************

        subroutine shympi_ex_3d_nodes_sum_r_internal(array)

        use basin
        use levels

        implicit none

        double precision array(nlvdi,nkn)

        double precision, dimension(:,:,:),allocatable :: data_send_array
        double precision, dimension(:,:,:),allocatable :: data_recv_array
        integer, dimension(:),allocatable :: sreqArray
        integer, dimension(:),allocatable :: rreqArray

        if(.not.allocated(sreqArray)) then
          allocate(sreqArray(domain%exchanges))
          allocate(rreqArray(domain%exchanges))
        end if

        call exchange_struct_3d_r(array,nlvdi,nkn,3,sreqArray,rreqArray,  &
     &                   data_send_array,data_recv_array)

        call waitAny(sreqArray)
        call waitAny(rreqArray)

        call rebuild_struct_3d_sum_r(array,nlvdi,nkn,data_recv_array,3)

        deallocate(data_send_array)
        deallocate(data_recv_array)

        deallocate(sreqArray)
        deallocate(rreqArray)

        end subroutine shympi_ex_3d_nodes_sum_r_internal

!*****************************************************************

!        subroutine shympi_ex_3d_nodes_sum_d_internal(array)
!
!        use basin
!        use levels
!
!        implicit none
!
!        double precision array(nlvdi,nkn)
!
!        double precision, dimension(:,:,:),allocatable::data_send_array
!        double precision, dimension(:,:,:),allocatable::data_recv_array
!        integer, dimension(:),allocatable :: sreqArray
!        integer, dimension(:),allocatable :: rreqArray
!
!        if(.not.allocated(sreqArray)) then
!          allocate(sreqArray(domain%exchanges))
!          allocate(rreqArray(domain%exchanges))
!        end if
!        
!        call exchange_struct_3d(array,nlvdi,nkn,3,sreqArray,rreqArray, &
!     &                   data_send_array,data_recv_array)
!
!        call waitAny(sreqArray)
!        call waitAny(rreqArray)
!
!        call rebuild_struct_3d_sum_d(array,nlvdi,nkn,data_recv_array,3)
!
!        deallocate(data_send_array)
!        deallocate(data_recv_array)
!
!        deallocate(sreqArray)
!        deallocate(rreqArray)
!
!        end subroutine shympi_ex_3d_nodes_sum_d_internal

!*****************************************************************

        subroutine shympi_ex_2d_nodes_sum_i_internal(array)

        use basin

        implicit none

        integer array(nkn)


        integer, dimension(:,:),allocatable :: data_send_array
        integer, dimension(:,:),allocatable :: data_recv_array
        integer, dimension(:),allocatable :: sreqArray
        integer, dimension(:),allocatable :: rreqArray

        if(.not.allocated(sreqArray)) then
          allocate(sreqArray(domain%exchanges))
          allocate(rreqArray(domain%exchanges))
        end if

        call exchange_struct_2d(array,3,sreqArray,rreqArray,  &
     &                   data_send_array,data_recv_array)


        call waitAny(sreqArray)
        call waitAny(rreqArray)

        deallocate(sreqArray)
        deallocate(rreqArray)

        call rebuild_struct_2d_sum(array,nkn,data_recv_array,3)


        deallocate(data_send_array)
        deallocate(data_recv_array)


        end subroutine shympi_ex_2d_nodes_sum_i_internal

!*****************************************************************

        subroutine shympi_ex_halo_2d_nodes_i_internal(array)

        use basin

        implicit none

        integer array(nkn_local)


        integer, dimension(:,:),allocatable :: data_send_array
        integer, dimension(:,:),allocatable :: data_recv_array
        integer, dimension(:),allocatable :: sreqArray
        integer, dimension(:),allocatable :: rreqArray

        if(.not.allocated(sreqArray)) then
          allocate(sreqArray(domain%exchanges))
          allocate(rreqArray(domain%exchanges))
        end if

        call exchange_struct_2d(array,5,sreqArray,rreqArray,  &
     &                   data_send_array,data_recv_array)


        call waitAny(sreqArray)
        call waitAny(rreqArray)

        deallocate(sreqArray)
        deallocate(rreqArray)

        call rebuild_struct_2d(array,nkn_local,data_recv_array,3)


        deallocate(data_send_array)
        deallocate(data_recv_array)


        end subroutine 

!*****************************************************************

        subroutine shympi_ex_halo_2d_elems_i_internal(array)

        use basin

        implicit none

        integer array(nel_local)


        integer, dimension(:,:),allocatable :: data_send_array
        integer, dimension(:,:),allocatable :: data_recv_array
        integer, dimension(:),allocatable :: sreqArray
        integer, dimension(:),allocatable :: rreqArray

        if(.not.allocated(sreqArray)) then
          allocate(sreqArray(domain%exchanges))
          allocate(rreqArray(domain%exchanges))
        end if

        call exchange_struct_2d(array,4,sreqArray,rreqArray, &
     &                   data_send_array,data_recv_array)


        call waitAny(sreqArray)
        call waitAny(rreqArray)

        deallocate(sreqArray)
        deallocate(rreqArray)

        call rebuild_struct_2d(array,nel_local,data_recv_array,2)


        deallocate(data_send_array)
        deallocate(data_recv_array)


        end subroutine 

!*****************************************************************

        subroutine shympi_ex_2d_nodes_sum_r_internal(array)

        use basin

        implicit none

        double precision array(nkn)


        double precision, dimension(:,:),allocatable :: data_send_array
        double precision, dimension(:,:),allocatable :: data_recv_array
        integer, dimension(:),allocatable :: sreqArray
        integer, dimension(:),allocatable :: rreqArray

        if(.not.allocated(sreqArray)) then
          allocate(sreqArray(domain%exchanges))
          allocate(rreqArray(domain%exchanges))
        end if

        call exchange_struct_2d(array,3,sreqArray,rreqArray, &
     &                   data_send_array,data_recv_array)


        call waitAny(sreqArray)
        call waitAny(rreqArray)

        deallocate(sreqArray)
        deallocate(rreqArray)

        call rebuild_struct_2d_sum(array,nkn,data_recv_array,3)


        deallocate(data_send_array)
        deallocate(data_recv_array)


        end subroutine shympi_ex_2d_nodes_sum_r_internal

!*****************************************************************

        subroutine shympi_exchange_coo_internal(array,index_coo,size1)

        use mpi_communication
        use basin

        implicit none
        integer size1
        double precision array(size1)
        integer index_coo(size1)


        double precision, dimension(:,:),allocatable :: data_send_array
        double precision, dimension(:,:),allocatable :: data_recv_array
        integer, dimension(:),allocatable :: sreqArray
        integer, dimension(:),allocatable :: rreqArray

        if(.not.allocated(sreqArray)) then
          allocate(sreqArray(domain%exchanges))
          allocate(rreqArray(domain%exchanges))
        end if

        call send_coo(array,sreqArray,rreqArray, &
     &              data_send_array,data_recv_array,index_coo,size1)


        call waitAny(sreqArray)
        call waitAny(rreqArray)

        deallocate(sreqArray)
        deallocate(rreqArray)

        call recv_coo(array,index_coo,size1,data_recv_array)


        deallocate(data_send_array)
        deallocate(data_recv_array)


        end subroutine 

!*****************************************************************

        subroutine shympi_ex_halo_2d_nodes_r_internal(array)

        use basin

        implicit none

        double precision array(nkn_local)


        double precision, dimension(:,:),allocatable :: data_send_array
        double precision, dimension(:,:),allocatable :: data_recv_array
        integer, dimension(:),allocatable :: sreqArray
        integer, dimension(:),allocatable :: rreqArray

        if(.not.allocated(sreqArray)) then
          allocate(sreqArray(domain%exchanges))
          allocate(rreqArray(domain%exchanges))
        end if

        call exchange_struct_2d(array,5,sreqArray,rreqArray, &
     &                   data_send_array,data_recv_array)


        call waitAny(sreqArray)
        call waitAny(rreqArray)

        deallocate(sreqArray)
        deallocate(rreqArray)

        call rebuild_struct_2d(array,nkn_local,data_recv_array,3)


        deallocate(data_send_array)
        deallocate(data_recv_array)


        end subroutine 

!*****************************************************************

        subroutine shympi_ex_halo_2d_elems_r_internal(array)

        use basin
        implicit none

        double precision array(nel_local)


        double precision, dimension(:,:),allocatable :: data_send_array
        double precision, dimension(:,:),allocatable :: data_recv_array
        integer, dimension(:),allocatable :: sreqArray
        integer, dimension(:),allocatable :: rreqArray

        if(.not.allocated(sreqArray)) then
          allocate(sreqArray(domain%exchanges))
          allocate(rreqArray(domain%exchanges))
        end if

        call exchange_struct_2d(array,4,sreqArray,rreqArray, &
     &                   data_send_array,data_recv_array)


        call waitAny(sreqArray)
        call waitAny(rreqArray)

        deallocate(sreqArray)
        deallocate(rreqArray)

        call rebuild_struct_2d(array,nel_local,data_recv_array,2)


        deallocate(data_send_array)
        deallocate(data_recv_array)


        end subroutine 

!*****************************************************************

        subroutine shympi_ex_2d_nodes_sum_d_internal(array)

        use basin
        implicit none

        double precision array(nkn)

        double precision, dimension(:,:),allocatable :: data_send_array
        double precision, dimension(:,:),allocatable :: data_recv_array
        integer, dimension(:),allocatable :: sreqArray
        integer, dimension(:),allocatable :: rreqArray

        if(.not.allocated(sreqArray)) then
          allocate(sreqArray(domain%exchanges))
          allocate(rreqArray(domain%exchanges))
        end if

        call exchange_struct_2d(array,3,sreqArray,rreqArray, &
     &                   data_send_array,data_recv_array)

        call waitAny(sreqArray)
        call waitAny(rreqArray)

        deallocate(sreqArray)
        deallocate(rreqArray)

        call rebuild_struct_2d_sum(array,nkn,data_recv_array,3)

        deallocate(data_send_array)
        deallocate(data_recv_array)

        end subroutine shympi_ex_2d_nodes_sum_d_internal

!*****************************************************************

        subroutine shympi_ex_halo_2d_nodes_d_internal(array)

        use basin
        implicit none

        double precision array(nkn_local)

        double precision, dimension(:,:),allocatable :: data_send_array
        double precision, dimension(:,:),allocatable :: data_recv_array
        integer, dimension(:),allocatable :: sreqArray
        integer, dimension(:),allocatable :: rreqArray

        if(.not.allocated(sreqArray)) then
          allocate(sreqArray(domain%exchanges))
          allocate(rreqArray(domain%exchanges))
        end if

        call exchange_struct_2d(array,5,sreqArray,rreqArray, &
     &                   data_send_array,data_recv_array)

        call waitAny(sreqArray)
        call waitAny(rreqArray)

        deallocate(sreqArray)
        deallocate(rreqArray)

        call rebuild_struct_2d(array,nkn_local,data_recv_array,3)

        deallocate(data_send_array)
        deallocate(data_recv_array)

        end subroutine 

!*****************************************************************

        subroutine shympi_ex_halo_2d_elems_d_internal(array)

        use basin
        implicit none

        double precision array(nel_local)

        double precision, dimension(:,:),allocatable :: data_send_array
        double precision, dimension(:,:),allocatable :: data_recv_array
        integer, dimension(:),allocatable :: sreqArray
        integer, dimension(:),allocatable :: rreqArray

        if(.not.allocated(sreqArray)) then
          allocate(sreqArray(domain%exchanges))
          allocate(rreqArray(domain%exchanges))
        end if

        call exchange_struct_2d(array,4,sreqArray,rreqArray, &
     &                   data_send_array,data_recv_array)

        call waitAny(sreqArray)
        call waitAny(rreqArray)

        deallocate(sreqArray)
        deallocate(rreqArray)

        call rebuild_struct_2d(array,nel_local,data_recv_array,2)

        deallocate(data_send_array)
        deallocate(data_recv_array)

        end subroutine 

!*****************************************************************

        subroutine shympi_ex_halo_3d2_nodes_i_internal(dim2,array)

        use basin
        use levels

        implicit none

        integer dim2
        integer array(dim2,nkn_local)

        integer, dimension(:,:,:),allocatable :: data_send_array
        integer, dimension(:,:,:),allocatable :: data_recv_array
        integer, dimension(:),allocatable :: sreqArray
        integer, dimension(:),allocatable :: rreqArray

        if(.not.allocated(sreqArray)) then
          allocate(sreqArray(domain%exchanges))
          allocate(rreqArray(domain%exchanges))
        end if

        call exchange_struct_3d_i(array,dim2,nkn_local,5,sreqArray, &
     &                  rreqArray,data_send_array,data_recv_array)

        call waitAny(sreqArray)
        call waitAny(rreqArray)

        call rebuild_struct_3d_i(array,dim2,nkn_local,  &
     &                  data_recv_array,3)

        deallocate(data_send_array)
        deallocate(data_recv_array)

        deallocate(sreqArray)
        deallocate(rreqArray)

        end subroutine 

!*****************************************************************

        subroutine shympi_ex_halo_3d_nodes_r_internal(array)

        use basin
        use levels

        implicit none

        double precision array(nlvdi,nkn_local)

        double precision, dimension(:,:,:),allocatable :: data_send_array
        double precision, dimension(:,:,:),allocatable :: data_recv_array
        integer, dimension(:),allocatable :: sreqArray
        integer, dimension(:),allocatable :: rreqArray

        if(.not.allocated(sreqArray)) then
          allocate(sreqArray(domain%exchanges))
          allocate(rreqArray(domain%exchanges))
        end if

        call exchange_struct_3d_r(array,nlvdi,nkn_local,5,sreqArray, &
     &                  rreqArray,data_send_array,data_recv_array)

        call waitAny(sreqArray)
        call waitAny(rreqArray)

        call rebuild_struct_3d_r(array,nlvdi,nkn_local, &
     &                  data_recv_array,3)

        deallocate(data_send_array)
        deallocate(data_recv_array)

        deallocate(sreqArray)
        deallocate(rreqArray)

        end subroutine 

!*****************************************************************

        subroutine shympi_ex_halo_3d0_nodes_r_internal(array)

        use basin
        use levels

        implicit none

        double precision array(0:nlvdi,nkn_local)

        double precision, dimension(:,:,:),allocatable :: data_send_array
        double precision, dimension(:,:,:),allocatable :: data_recv_array
        integer, dimension(:),allocatable :: sreqArray
        integer, dimension(:),allocatable :: rreqArray

        if(.not.allocated(sreqArray)) then
          allocate(sreqArray(domain%exchanges))
          allocate(rreqArray(domain%exchanges))
        end if

        call exchange_struct_3d_r(array(0:nlvdi,:),nlvdi+1,nkn_local,5,sreqArray, &
     &                  rreqArray,data_send_array,data_recv_array)

        call waitAny(sreqArray)
        call waitAny(rreqArray)

        call rebuild_struct_3d_r(array(0:nlvdi,:),nlvdi+1,nkn_local, &
     &                  data_recv_array,3)

        deallocate(data_send_array)
        deallocate(data_recv_array)

        deallocate(sreqArray)
        deallocate(rreqArray)

        end subroutine 

!*****************************************************************

!        subroutine shympi_ex_halo_3d_nodes_d_internal(array)
!
!        use basin
!        use levels
!
!        implicit none
!
!        double precision array(nlvdi,nkn_local)
!
!        double precision, dimension(:,:,:),allocatable :: &
!     &          data_send_array,data_recv_array
!        integer, dimension(:),allocatable :: sreqArray
!        integer, dimension(:),allocatable :: rreqArray
!
!        if(.not.allocated(sreqArray)) then
!          allocate(sreqArray(domain%exchanges))
!          allocate(rreqArray(domain%exchanges))
!        end if
!
!        call exchange_struct_3d_d(array,nlvdi,nkn_local,5,sreqArray, &
!     &          rreqArray,data_send_array,data_recv_array)
!
!        call waitAny(sreqArray)
!        call waitAny(rreqArray)
!
!        call rebuild_struct_3d_d(array,nlvdi,nkn_local,data_recv_array,3)
!
!        deallocate(data_send_array)
!        deallocate(data_recv_array)
!
!        deallocate(sreqArray)
!        deallocate(rreqArray)
!
!        end subroutine 

!*****************************************************************

        subroutine shympi_ex_halo_3d_elems_r_internal(array)

        use basin
        use levels

        implicit none

        double precision array(nlvdi,nel_local)

        double precision, dimension(:,:,:),allocatable :: data_send_array
        double precision, dimension(:,:,:),allocatable :: data_recv_array
        integer, dimension(:),allocatable :: sreqArray
        integer, dimension(:),allocatable :: rreqArray

        if(.not.allocated(sreqArray)) then
          allocate(sreqArray(domain%exchanges))
          allocate(rreqArray(domain%exchanges))
        end if

        call exchange_struct_3d_r(array,nlvdi,nel_local,4,sreqArray,  &
     &                    rreqArray,data_send_array,data_recv_array)

        call waitAny(sreqArray)
        call waitAny(rreqArray)

        call rebuild_struct_3d_r(array,nlvdi,nel_local,  &
     &                  data_recv_array,2)

        deallocate(data_send_array)
        deallocate(data_recv_array)

        deallocate(sreqArray)
        deallocate(rreqArray)

        end subroutine 

!*****************************************************************

        subroutine shympi_ex_halo_4d_elems_r_internal(newdim,array)

        use basin

        implicit none

        integer newdim
        double precision array(newdim,nel_local)

        double precision, dimension(:,:,:),allocatable :: data_send_array
        double precision, dimension(:,:,:),allocatable :: data_recv_array
        integer, dimension(:),allocatable :: sreqArray
        integer, dimension(:),allocatable :: rreqArray

        if(.not.allocated(sreqArray)) then
          allocate(sreqArray(domain%exchanges))
          allocate(rreqArray(domain%exchanges))
        end if

        call exchange_struct_3d_r(array,newdim,nel_local,4,sreqArray,  &
     &                    rreqArray,data_send_array,data_recv_array)

        call waitAny(sreqArray)
        call waitAny(rreqArray)

        call rebuild_struct_3d_r(array,newdim,nel_local,  &
     &                  data_recv_array,2)

        deallocate(data_send_array)
        deallocate(data_recv_array)

        deallocate(sreqArray)
        deallocate(rreqArray)

        end subroutine 

!*****************************************************************

        subroutine shympi_ex_halo_3d_elems_r3_internal(array)

        use basin

        implicit none

        double precision array(3,nel_local)

        double precision, dimension(:,:,:),allocatable :: data_send_array
        double precision, dimension(:,:,:),allocatable :: data_recv_array
        integer, dimension(:),allocatable :: sreqArray
        integer, dimension(:),allocatable :: rreqArray

        if(.not.allocated(sreqArray)) then
          allocate(sreqArray(domain%exchanges))
          allocate(rreqArray(domain%exchanges))
        end if

        call exchange_struct_3d_r(array,3,nel_local,4,sreqArray,  &
     &                    rreqArray,data_send_array,data_recv_array)

        call waitAny(sreqArray)
        call waitAny(rreqArray)

        call rebuild_struct_3d_r(array,3,nel_local,  &
     &                  data_recv_array,2)

        deallocate(data_send_array)
        deallocate(data_recv_array)

        deallocate(sreqArray)
        deallocate(rreqArray)

        end subroutine 

!*****************************************************************

!        subroutine shympi_ex_halo_3d_elems_d_internal(array)
!
!        use basin
!        use levels
!
!        implicit none
!
!        double precision array(nlvdi,nel_local)
!
!        double precision, dimension(:,:,:),allocatable :: data_send_array,data_recv_array
!        integer, dimension(:),allocatable :: sreqArray
!        integer, dimension(:),allocatable :: rreqArray
!
!        if(.not.allocated(sreqArray)) then
!          allocate(sreqArray(domain%exchanges))
!          allocate(rreqArray(domain%exchanges))
!        end if
!
!        call exchange_struct_3d_d(array,nlvdi,nel_local,4,sreqArray, &
!     &                  rreqArray,data_send_array,data_recv_array)
!
!        call waitAny(sreqArray)
!        call waitAny(rreqArray)
!
!        call rebuild_struct_3d_d(array,nlvdi,nel_local,data_recv_array,2)
!
!        deallocate(data_send_array)
!        deallocate(data_recv_array)
!
!        deallocate(sreqArray)
!        deallocate(rreqArray)
!
!        end subroutine 

!*****************************************************************

        subroutine shympi_ex_2d_nodes_min_i_internal(array)

        use basin
        implicit none

        integer array(nkn)


        integer, dimension(:,:),allocatable :: data_send_array
        integer, dimension(:,:),allocatable :: data_recv_array
        integer, dimension(:),allocatable :: sreqArray
        integer, dimension(:),allocatable :: rreqArray

        if(.not.allocated(sreqArray)) then
          allocate(sreqArray(domain%exchanges))
          allocate(rreqArray(domain%exchanges))
        end if

        call exchange_struct_2d_i(array,3,sreqArray,rreqArray,  &
     &                   data_send_array,data_recv_array)


        call waitAny(sreqArray)
        call waitAny(rreqArray)

        deallocate(sreqArray)
        deallocate(rreqArray)

        call rebuild_struct_2d_min_i(array,nkn,data_recv_array,3)


        deallocate(data_send_array)
        deallocate(data_recv_array)


        end subroutine shympi_ex_2d_nodes_min_i_internal

!*****************************************************************

        subroutine shympi_ex_2d_nodes_min_r_internal(array)

        use basin
        implicit none

        double precision array(nkn)


        double precision, dimension(:,:),allocatable :: data_send_array
        double precision, dimension(:,:),allocatable :: data_recv_array
        integer, dimension(:),allocatable :: sreqArray
        integer, dimension(:),allocatable :: rreqArray

        if(.not.allocated(sreqArray)) then
          allocate(sreqArray(domain%exchanges))
          allocate(rreqArray(domain%exchanges))
        end if

        call exchange_struct_2d_r(array,3,sreqArray,rreqArray,  &
     &                   data_send_array,data_recv_array)


        call waitAny(sreqArray)
        call waitAny(rreqArray)

        deallocate(sreqArray)
        deallocate(rreqArray)

        call rebuild_struct_2d_min_r(array,nkn,data_recv_array,3)


        deallocate(data_send_array)
        deallocate(data_recv_array)


        end subroutine shympi_ex_2d_nodes_min_r_internal

!*****************************************************************

        subroutine shympi_ex_2d_nodes_max_i_internal(array)

        use basin
        implicit none

        !integer, dimension(:) :: array
        integer array(nkn)


        integer, dimension(:,:),allocatable :: data_send_array
        integer, dimension(:,:),allocatable :: data_recv_array
        integer, dimension(:),allocatable :: sreqArray
        integer, dimension(:),allocatable :: rreqArray

        if(.not.allocated(sreqArray)) then
          allocate(sreqArray(domain%exchanges))
          allocate(rreqArray(domain%exchanges))
        end if

        call exchange_struct_2d_i(array,3,sreqArray,rreqArray,  &
     &                   data_send_array,data_recv_array)


        call waitAny(sreqArray)
        call waitAny(rreqArray)

        deallocate(sreqArray)
        deallocate(rreqArray)

        call rebuild_struct_2d_max_i(array,nkn,data_recv_array,3)

        deallocate(data_send_array)
        deallocate(data_recv_array)


        end subroutine shympi_ex_2d_nodes_max_i_internal

!*****************************************************************

        subroutine shympi_ex_2d_nodes_max_r_internal(array)

        use basin
        implicit none

        !integer, dimension(:) :: array
        double precision array(nkn)


        double precision, dimension(:,:),allocatable :: data_send_array
        double precision, dimension(:,:),allocatable :: data_recv_array
        integer, dimension(:),allocatable :: sreqArray
        integer, dimension(:),allocatable :: rreqArray

        if(.not.allocated(sreqArray)) then
          allocate(sreqArray(domain%exchanges))
          allocate(rreqArray(domain%exchanges))
        end if

        call exchange_struct_2d_r(array,3,sreqArray,rreqArray,  &
     &                   data_send_array,data_recv_array)


        call waitAny(sreqArray)
        call waitAny(rreqArray)

        deallocate(sreqArray)
        deallocate(rreqArray)

        call rebuild_struct_2d_max_r(array,nkn,data_recv_array,3)

        deallocate(data_send_array)
        deallocate(data_recv_array)


        end subroutine shympi_ex_2d_nodes_max_r_internal

!******************************************************************

       subroutine send_halo_int_r(array,dim2,dim1,sreqv,rreqv,what)


         implicit none

         integer dim2, dim1
         double precision array(dim2,dim1)
         character*(2) what
         integer sreqv(1),rreqv(1)
         

         if(what == 'ut') then
         call exchange_struct_3d(array,dim2,dim1,2,sreqv,rreqv,  &
     &                        data_send_ut ,data_recv_ut)

         else if(what == 'vt') then
         call exchange_struct_3d(array,dim2,dim1,2,sreqv,rreqv,  &
     &                         data_send_vt ,data_recv_vt)
         else
           write(6,*)'what different from ut or vt in send_halo_int'
           stop
         end if

        end

!******************************************************************


        subroutine recv_halo_int_r(array,dim2,dim1,what)


         implicit none

         integer dim2, dim1
         double precision array(dim2,dim1)
         character*(2) what


         if(what == 'ut') then
           call waitAny(sreq_ut)
           call waitAny(rreq_ut)
           call rebuild_RealHalo(array,dim2,data_recv_ut,2)
           deallocate(data_send_ut)
           deallocate(data_recv_ut)
           deallocate(sreq_ut)
           deallocate(rreq_ut)
         else if(what == 'vt') then
           call waitAny(sreq_vt)
           call waitAny(rreq_vt)
           call rebuild_RealHalo(array,dim2,data_recv_vt,2)
           deallocate(data_send_vt)
           deallocate(data_recv_vt)
           deallocate(sreq_vt)
           deallocate(rreq_vt)
         else
           write(6,*)'what is different from uv or vt in recv_halo'
           stop
         end if

       end subroutine


!******************************************************************

       subroutine send_halo_int_d(array,dim2,dim1)


         implicit none

         integer dim2, dim1
         double precision array(dim2,dim1)
       
         if(.not.allocated(sreq)) then
           allocate(sreq(domain%exchanges))
           allocate(rreq(domain%exchanges))
         end if

         call exchange_struct_3d(array,dim2,dim1,2,sreq,rreq,data_send_d  &
     &                          ,data_recv_d)

        end

!******************************************************************


        subroutine recv_halo_int_d(array,dim2,dim1)


         implicit none

         integer dim2, dim1
         double precision array(dim2,dim1)

         call waitAny(sreq)
         call waitAny(rreq)
         call rebuild_DPHalo(array,dim2,dim1,data_recv_d,2)

         deallocate(data_send_d)
         deallocate(data_recv_d)
         deallocate(sreq)
         deallocate(rreq)

       end subroutine

!*****************************************************************
!*****************************************************************
!*****************************************************************

        function shympi_wtime_internal()

        use mpi_common_struct
        implicit none

	double precision shympi_wtime_internal

	shympi_wtime_internal = MPI_WTIME()

        end function shympi_wtime_internal

!******************************************************************


        end module shympi_internal

