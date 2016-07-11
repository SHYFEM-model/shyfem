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

        module shympi_aux

        use mpi

        implicit none

        !include 'mpif.h'

        end module shympi_aux

!********************************

	subroutine shympi_init_internal(my_id,n_threads)

	use shympi_aux

	implicit none

	integer my_id,n_threads

	integer ierr

        call MPI_INIT( ierr )
        call MPI_COMM_RANK( MPI_COMM_WORLD, my_id, ierr )
        call MPI_COMM_SIZE( MPI_COMM_WORLD, n_threads, ierr )

	end subroutine shympi_init_internal

!******************************************************************

	subroutine shympi_barrier_internal

	use shympi_aux

	implicit none

	integer ierr

	call MPI_BARRIER( MPI_COMM_WORLD, ierr)

	end subroutine shympi_barrier_internal

!******************************************************************

        subroutine shympi_abort_internal

	use shympi_aux

        implicit none

        integer ierr,ierr_code

        ierr_code = 33
	call MPI_ABORT(MPI_COMM_WORLD,ierr_code,ierr)

        end subroutine shympi_abort_internal

!******************************************************************

        subroutine shympi_finalize_internal

	use shympi_aux

        implicit none

        integer ierr

	call MPI_BARRIER( MPI_COMM_WORLD, ierr)
	call MPI_FINALIZE(ierr)

        end subroutine shympi_finalize_internal

!******************************************************************

        function shympi_wtime_internal()

	use shympi_aux

        implicit none

	double precision shympi_wtime_internal

	shympi_wtime_internal = MPI_WTIME()

        end function shympi_wtime_internal

!******************************************************************

        subroutine shympi_get_status_size_internal(size)

	use shympi_aux

        implicit none

        integer size

	size = MPI_STATUS_SIZE

        end subroutine shympi_get_status_size_internal

!******************************************************************

	subroutine shympi_syncronize_internal

	use shympi_aux

	implicit none

	integer ierr

	flush(6)
	call MPI_BARRIER( MPI_COMM_WORLD, ierr)

	end subroutine shympi_syncronize_internal

!******************************************************************
!******************************************************************
!******************************************************************

        subroutine shympi_gather_i_internal(val)

	use shympi_aux

	use shympi

	implicit none

        integer val

        integer ierr

        call MPI_GATHER (val,1,MPI_INT
     +                  ,ival,1,MPI_INT
     +                  ,0,MPI_COMM_WORLD,ierr)

        end subroutine shympi_gather_i_internal

!*******************************

        subroutine shympi_bcast_i_internal(val)

	use shympi_aux

	implicit none

        integer val

        integer ierr

        call MPI_BCAST(val,1,MPI_INT,0,MPI_COMM_WORLD,ierr)

        end subroutine shympi_bcast_i_internal

!*******************************

	subroutine shympi_reduce_r_internal(what,val)

	use shympi_aux

	implicit none

	character*(*) what
	real val

        integer ierr
	real valout

        if( what == 'min' ) then
	  call MPI_ALLREDUCE(val,valout,1,MPI_REAL,MPI_MIN
     +				,MPI_COMM_WORLD,ierr)
	  val = valout
        else if( what == 'max' ) then
	  call MPI_ALLREDUCE(val,valout,1,MPI_REAL,MPI_MAX
     +				,MPI_COMM_WORLD,ierr)
	  val = valout
        else if( what == 'sum' ) then
	  call MPI_ALLREDUCE(val,valout,1,MPI_REAL,MPI_SUM
     +				,MPI_COMM_WORLD,ierr)
	  val = valout
        else
          write(6,*) 'what = ',what
          stop 'error stop shympi_reduce_r_internal: not ready'
        end if

	end subroutine shympi_reduce_r_internal

!******************************************************************

	subroutine shympi_reduce_i_internal(what,val)

	use shympi_aux

	implicit none

	character*(*) what
	integer val

        integer ierr
	integer valout

        if( what == 'min' ) then
	  call MPI_ALLREDUCE(val,valout,1,MPI_INTEGER,MPI_MIN
     +				,MPI_COMM_WORLD,ierr)
	  val = valout
        else if( what == 'max' ) then
	  call MPI_ALLREDUCE(val,valout,1,MPI_INTEGER,MPI_MAX
     +				,MPI_COMM_WORLD,ierr)
	  val = valout
        else if( what == 'sum' ) then
	  call MPI_ALLREDUCE(val,valout,1,MPI_INTEGER,MPI_SUM
     +				,MPI_COMM_WORLD,ierr)
	  val = valout
        else
          write(6,*) 'what = ',what
          stop 'error stop shympi_reduce_i_internal: not ready'
        end if

	end subroutine shympi_reduce_i_internal

!*****************************************************************

	subroutine shympi_setup

        use basin
        use shympi
        use Graph_Ele
        use mod_geom
        use communicationStruct

	implicit none

        integer i,nlkdi,ierr
        integer, dimension(:,:), allocatable :: mystruct
        integer temp_nkn,temp_nel

        nel = 0
        do i=1,neldi
          if(allPartAssign(i) .eq. my_id) then
            nel = nel + 1
          end if
        end do

        write(6,*)'myElements =',nel,neldi,my_id

        nlkdi = 3*neldi+2*nkndi

        allocate(ilinkv(nkndi+1))
        allocate(lenkv(nlkdi))
        allocate(linkv(nlkdi))
        allocate(ieltv(3,neldi))
        allocate(total_ieltv(3,neldi))

        call mklenk(nlkdi,nkndi,neldi,nen3v,ilinkv,lenkv)

        call mklink(nkndi,ilinkv,lenkv,linkv)

        call mkielt(nkndi,neldi,ilinkv,lenkv,linkv,ieltv)

        total_ieltv = ieltv

        call makeMyEle(nel,neldi,nen3v)

        allocate(mystruct(3,myele%totalID))

        call makeMyStruct(myele, nen3v, mystruct, 3)

        call makeMyNodes(myele, mystruct, mynodes, nen3v, nkndi)
        
        nkn =mynodes%numberID
        nel =myele%numberID
 
        if(nel .eq. 0) then
          ngr = 0
        end if
        !write(6,*)'myNodes =',nkn,nkndi,my_id

        temp_nkn=nkndi
        temp_nel=neldi
        call transfer_domain
        neldi=temp_nel
        nkndi=temp_nkn

        nel_tot=myele%totalID
 
        !deallocate(nen3v)
        !allocate(nen3v(3,myele%totalID))

        !call makeMyNen3v(myele, mynodes, mystruct, nen3v)

        !call sleep(4)
        !call shympi_barrier

        deallocate(ilinkv)
        deallocate(linkv)
        deallocate(lenkv)
        deallocate(ieltv)

        call shympi_alloc
        
        !nkn =mynodes%numberID
        !nel =myele%numberID

        !call shympi_alloc

        !call sleep(5)
        
        !stop


	end subroutine shympi_setup

!******************************************************************

        subroutine shympi_ex_3d_nodes_sum_r_internal(array)

        use MPI_Communications
        use basin
        use levels

        implicit none

        real array(nlvdi,nkn)

        real, dimension(:,:,:),allocatable :: data_send_array
        real, dimension(:,:,:),allocatable :: data_recv_array
        integer, dimension(:),allocatable :: sreqArray
        integer, dimension(:),allocatable :: rreqArray

        allocate(sreqArray(mypart%mysend%sends))
        allocate(rreqArray(mypart%myreceive%receives))

        call exchange_struct_3D_r(array,nlvdi,3,sreqArray,rreqArray,
     +                   data_send_array,data_recv_array)

        call waitAny(sreqArray)
        call waitAny(rreqArray)

        call rebuild_struct_3d_sum_r(array,nlvdi,data_recv_array,3)

        deallocate(data_send_array)
        deallocate(data_recv_array)

        deallocate(sreqArray)
        deallocate(rreqArray)

        end subroutine shympi_ex_3d_nodes_sum_r_internal

!******************************************************************

        subroutine shympi_ex_3d_nodes_sum_d_internal(array)

        use MPI_Communications
        use basin
        use levels

        implicit none

        double precision array(nlvdi,nkn)

        double precision, dimension(:,:,:),allocatable::data_send_array
        double precision, dimension(:,:,:),allocatable::data_recv_array
        integer, dimension(:),allocatable :: sreqArray
        integer, dimension(:),allocatable :: rreqArray

        allocate(sreqArray(mypart%mysend%sends))
        allocate(rreqArray(mypart%myreceive%receives))

        call exchange_struct_3D(array,nlvdi,3,sreqArray,rreqArray,
     +                   data_send_array,data_recv_array)

        call waitAny(sreqArray)
        call waitAny(rreqArray)

        call rebuild_struct_3D_sum_d(array,nlvdi,data_recv_array,3)

        deallocate(data_send_array)
        deallocate(data_recv_array)

        deallocate(sreqArray)
        deallocate(rreqArray)

        end subroutine shympi_ex_3d_nodes_sum_d_internal

!*****************************************************************

        subroutine shympi_ex_2d_nodes_sum_r_internal(array)

        use MPI_Communications
        use basin

        implicit none

        real array(nkn)


        real, dimension(:,:),allocatable :: data_send_array
        real, dimension(:,:),allocatable :: data_recv_array
        integer, dimension(:),allocatable :: sreqArray
        integer, dimension(:),allocatable :: rreqArray

        allocate(sreqArray(mypart%mysend%sends))
        allocate(rreqArray(mypart%myreceive%receives))

        call exchange_struct_2d(array,3,sreqArray,rreqArray,
     +                   data_send_array,data_recv_array)


        call waitAny(sreqArray)
        call waitAny(rreqArray)

        deallocate(sreqArray)
        deallocate(rreqArray)

        call rebuild_struct_2d_sum(array,data_recv_array,3)


        deallocate(data_send_array)
        deallocate(data_recv_array)


        end subroutine shympi_ex_2d_nodes_sum_r_internal

!*****************************************************************

        subroutine shympi_ex_2d_nodes_sum_d_internal(array)

        use MPI_Communications
        use basin

        implicit none

        double precision array(nkn)

        double precision, dimension(:,:),allocatable :: data_send_array
        double precision, dimension(:,:),allocatable :: data_recv_array
        integer, dimension(:),allocatable :: sreqArray
        integer, dimension(:),allocatable :: rreqArray

        allocate(sreqArray(mypart%mysend%sends))
        allocate(rreqArray(mypart%myreceive%receives))

        call exchange_struct_2d(array,3,sreqArray,rreqArray,
     +                   data_send_array,data_recv_array)

        call waitAny(sreqArray)
        call waitAny(rreqArray)

        deallocate(sreqArray)
        deallocate(rreqArray)

        call rebuild_struct_2d_sum(array,data_recv_array,3)

        deallocate(data_send_array)
        deallocate(data_recv_array)

        end subroutine shympi_ex_2d_nodes_sum_d_internal

!*****************************************************************

        subroutine shympi_ex_2d_nodes_min_i_internal(array)

        use MPI_Communications
        use basin

        implicit none

        integer array(nkn)


        integer, dimension(:,:),allocatable :: data_send_array
        integer, dimension(:,:),allocatable :: data_recv_array
        integer, dimension(:),allocatable :: sreqArray
        integer, dimension(:),allocatable :: rreqArray

        allocate(sreqArray(mypart%mysend%sends))
        allocate(rreqArray(mypart%myreceive%receives))

        call exchange_struct_2d_i(array,3,sreqArray,rreqArray,
     +                   data_send_array,data_recv_array)


        call waitAny(sreqArray)
        call waitAny(rreqArray)

        deallocate(sreqArray)
        deallocate(rreqArray)

        call rebuild_struct_2d_min_i(array,data_recv_array,3)


        deallocate(data_send_array)
        deallocate(data_recv_array)


        end subroutine shympi_ex_2d_nodes_min_i_internal

!*****************************************************************

        subroutine shympi_ex_2d_nodes_min_r_internal(array)

        use MPI_Communications
        use basin

        implicit none

        real array(nkn)


        real, dimension(:,:),allocatable :: data_send_array
        real, dimension(:,:),allocatable :: data_recv_array
        integer, dimension(:),allocatable :: sreqArray
        integer, dimension(:),allocatable :: rreqArray

        allocate(sreqArray(mypart%mysend%sends))
        allocate(rreqArray(mypart%myreceive%receives))

        call exchange_struct_2d_r(array,3,sreqArray,rreqArray,
     +                   data_send_array,data_recv_array)


        call waitAny(sreqArray)
        call waitAny(rreqArray)

        deallocate(sreqArray)
        deallocate(rreqArray)

        call rebuild_struct_2d_min_r(array,data_recv_array,3)


        deallocate(data_send_array)
        deallocate(data_recv_array)


        end subroutine shympi_ex_2d_nodes_min_r_internal

!*****************************************************************

        subroutine shympi_ex_2d_nodes_max_i_internal(array)

        use MPI_Communications
        use basin

        implicit none

        !integer, dimension(:) :: array
        integer array(nkn)


        integer, dimension(:,:),allocatable :: data_send_array
        integer, dimension(:,:),allocatable :: data_recv_array
        integer, dimension(:),allocatable :: sreqArray
        integer, dimension(:),allocatable :: rreqArray

        allocate(sreqArray(mypart%mysend%sends))
        allocate(rreqArray(mypart%myreceive%receives))

        call exchange_struct_2d_i(array,3,sreqArray,rreqArray,
     +                   data_send_array,data_recv_array)


        call waitAny(sreqArray)
        call waitAny(rreqArray)

        deallocate(sreqArray)
        deallocate(rreqArray)

        call rebuild_struct_2d_max_i(array,data_recv_array,3)

        deallocate(data_send_array)
        deallocate(data_recv_array)


        end subroutine shympi_ex_2d_nodes_max_i_internal

!*****************************************************************

        subroutine shympi_ex_2d_nodes_max_r_internal(array)

        use MPI_Communications
        use basin

        implicit none

        !integer, dimension(:) :: array
        real array(nkn)


        real, dimension(:,:),allocatable :: data_send_array
        real, dimension(:,:),allocatable :: data_recv_array
        integer, dimension(:),allocatable :: sreqArray
        integer, dimension(:),allocatable :: rreqArray

        allocate(sreqArray(mypart%mysend%sends))
        allocate(rreqArray(mypart%myreceive%receives))

        call exchange_struct_2d_r(array,3,sreqArray,rreqArray,
     +                   data_send_array,data_recv_array)


        call waitAny(sreqArray)
        call waitAny(rreqArray)

        deallocate(sreqArray)
        deallocate(rreqArray)

        call rebuild_struct_2d_max_r(array,data_recv_array,3)

        deallocate(data_send_array)
        deallocate(data_recv_array)


        end subroutine shympi_ex_2d_nodes_max_r_internal

!******************************************************************

!******************************************************************
!******************************************************************
!******************************************************************

