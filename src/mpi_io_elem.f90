
!==================================================================
        module mpi_io_admin
!==================================================================

        use shympi
        use timing

        double precision, allocatable, dimension(:,:) :: inTempv
        double precision, allocatable, dimension(:,:) :: inSaltv
        integer, allocatable, dimension(:) :: inIlhkv
        double precision, allocatable, dimension(:) :: inHkv
        double precision, allocatable, dimension(:) :: inHev

        double precision, allocatable, dimension(:) :: outZnv
        double precision, allocatable, dimension(:) :: outHev
        double precision, allocatable, dimension(:,:) :: outZenv,outSaux
        integer, allocatable, dimension(:) :: outIlhv,outIlhkv
        double precision, allocatable, dimension(:,:) :: outUtlnv,outVtlnv
        double precision, allocatable, dimension(:,:) :: outTempv

        !double precision, allocatable,save :: hlv_fem_temp(:)
        double precision, allocatable,save :: hk_fem_temp(:)
        double precision, allocatable,save :: he_fem_temp(:)
        integer, allocatable,save :: ilhkv_fem_temp(:)
        !integer,save :: hlv_dim_new
	integer,save :: iminreg,imaxreg,jminreg,jmaxreg
	integer,save :: nxreg,nyreg
	double precision,save :: x0reg,y0reg

        INTERFACE rebuild_2d_nodes
        	MODULE PROCEDURE  &
     			  rebuild_2d_nodes_i  &
     			  ,rebuild_2d_nodes_r !&
!     			  ,rebuild_2d_nodes_d
        END INTERFACE

        INTERFACE rebuild_3d_elems3
        	MODULE PROCEDURE  & 
     			  rebuild_3d_elems3_r &
     			  ,rebuild_3d_elems3_i
        END INTERFACE

        INTERFACE rebuild_2d_elems
        	MODULE PROCEDURE &
     			  rebuild_2d_elems_i &
     			  ,rebuild_2d_elems_r
        END INTERFACE


        INTERFACE rebuild_3d_nodes0
        	MODULE PROCEDURE & 
     			  rebuild_3d_nodes0_r
        END INTERFACE

        INTERFACE rebuild_3d_nodes3
        	MODULE PROCEDURE & 
     			  rebuild_3d_nodes3_i
        END INTERFACE

        INTERFACE rebuild_3d_nodes
        	MODULE PROCEDURE & 
     			  rebuild_3d_nodes_r &
     &			  ,rebuild_3d_nodes_r2 !&
!     			  ,rebuild_3d_nodes_d
        END INTERFACE

        INTERFACE rebuild_3d_elems
        	MODULE PROCEDURE & 
     			  rebuild_3d_elems_r
        END INTERFACE

        contains


!******************************************************************

        subroutine rebuild_2d_nodes_i(array, reb_array)

        use basin
        use mpi_common_struct

        implicit none

        integer array(nkn)
        integer, allocatable, optional, dimension(:) :: reb_array

        if(my_id .eq. 0) then
          if(.not. allocated(reb_array)) then
            allocate(reb_array(nkndi))
          end if
          call rebuild_integ_2d(array, nkndi, numberNodes, 2,reb_array)
        else
          call rebuild_integ_2d(array, nkndi, numberNodes, 2)
        end if

        end subroutine rebuild_2d_nodes_i

!******************************************************************

        subroutine rebuild_2d_nodes_r(array, reb_array)

        use basin
        use mpi_common_struct

        implicit none

        double precision array(nkn)
        double precision, allocatable, optional, dimension(:) :: reb_array

        if(my_id .eq. 0) then
          if(.not. allocated(reb_array)) then
            allocate(reb_array(nkndi))
          end if
          call rebuild_real_2d(array, nkndi, numberNodes, 2,reb_array)
        else
          call rebuild_real_2d(array, nkndi, numberNodes, 2)
        end if


        end subroutine rebuild_2d_nodes_r

!******************************************************************

!        subroutine rebuild_2d_nodes_d(array, reb_array)
!
!        use basin
!        use mpi_common_struct
!
!        implicit none
!
!        double precision array(nkn)
!        double precision,allocatable,optional,dimension(:) :: reb_array
!
!        if(my_id .eq. 0) then
!          if(.not. allocated(reb_array)) then
!            allocate(reb_array(nkndi))
!          end if
!          call rebuild_dp_2d(array, nkndi, numberNodes, 2,reb_array)
!        else
!          call rebuild_dp_2d(array, nkndi, numberNodes, 2)
!        end if
!
!
!        end subroutine rebuild_2d_nodes_d

!******************************************************************

        subroutine rebuild_3d_elems3_r(array, reb_array)

        use basin
        use mpi_common_struct
        implicit none

        double precision array(3,nel)
        double precision, allocatable, optional, dimension(:,:) :: reb_array

        if(my_id .eq. 0) then
          if(.not. allocated(reb_array)) then
            allocate(reb_array(3,neldi))
          end if
          call rebuild_real_3d(array,3, neldi, numberElements,1,reb_array)
        else
          call rebuild_real_3d(array,3, neldi, numberElements,1)
        end if

        end subroutine rebuild_3d_elems3_r

!******************************************************************
        subroutine rebuild_3d_elems3_i(array, reb_array)

        use basin
        use mpi_common_struct
        implicit none

        integer array(3,nel)
        integer, allocatable, optional, dimension(:,:) :: reb_array

        if(my_id .eq. 0) then
          if(.not. allocated(reb_array)) then
            allocate(reb_array(3,neldi))
          end if
          call rebuild_integ_3d(array,3, neldi, numberElements,1,reb_array)
        else
          call rebuild_integ_3d(array,3, neldi, numberElements,1)
        end if

        end subroutine rebuild_3d_elems3_i

!******************************************************************

        subroutine rebuild_2d_elems_i(array,reb_array)

        use basin
        use mpi_common_struct
        implicit none

        integer array(nel)
        integer, allocatable, optional, dimension(:) :: reb_array

        if(my_id .eq. 0) then
          if(.not. allocated(reb_array)) then
            allocate(reb_array(neldi))
          end if
          call rebuild_integ_2d(array, neldi, numberElements, 1,reb_array)
        else
          call rebuild_integ_2d(array, neldi, numberElements, 1)
        end if

        end subroutine rebuild_2d_elems_i

!******************************************************************

        subroutine rebuild_2d_elems_r(array,reb_array)

        use basin
        use mpi_common_struct
        implicit none

        double precision array(nel)
        double precision, allocatable, optional, dimension(:) :: reb_array

        if(my_id .eq. 0) then
          if(.not. allocated(reb_array)) then
            allocate(reb_array(neldi))
          end if
          call rebuild_real_2d(array, neldi, numberElements, 1,reb_array)
        else
          call rebuild_real_2d(array, neldi, numberElements, 1)
        end if

        end subroutine rebuild_2d_elems_r

!******************************************************************

        subroutine rebuild_3d_nodes0_r(array, reb_array)

        use basin
        use levels
        use mpi_common_struct

        implicit none

        double precision array(nlvdi+1,nkn)
        double precision, allocatable, optional, dimension(:,:) :: reb_array

        if(my_id .eq. 0) then
          if(.not. allocated(reb_array)) then
            allocate(reb_array(nlvdi+1,nkndi))
          end if
          call rebuild_real_3d(array,nlvdi+1,nkndi,numberNodes, 2,reb_array)
        else
          call rebuild_real_3d(array,nlvdi+1,nkndi,numberNodes, 2)
        end if

        end subroutine rebuild_3d_nodes0_r


!******************************************************************

        subroutine rebuild_3d_nodes_r2(dim2, array, reb_array)

        use basin
        use mpi_common_struct

        implicit none

        integer dim2
        double precision array(dim2,nkn)
        double precision, allocatable, optional, dimension(:,:) :: reb_array

        if(my_id .eq. 0) then
          if(.not. allocated(reb_array)) then
            allocate(reb_array(dim2,nkndi))
          else
            deallocate(reb_array)
            allocate(reb_array(dim2,nkndi))
          end if
          call rebuild_real_3d(array,dim2,nkndi,numberNodes, 2,reb_array)
        else
          call rebuild_real_3d(array,dim2,nkndi,numberNodes, 2)
        end if

        end subroutine rebuild_3d_nodes_r2

!******************************************************************

        subroutine rebuild_3d_nodes_r(array, reb_array)

        use basin
        use levels
        use mpi_common_struct

        implicit none

        double precision array(nlvdi,nkn)
        double precision, allocatable, optional, dimension(:,:) :: reb_array

        if(my_id .eq. 0) then
          if(.not. allocated(reb_array)) then
            allocate(reb_array(nlvdi,nkndi))
          end if
          call rebuild_real_3d(array,nlvdi,nkndi,numberNodes, 2,reb_array)
        else
          call rebuild_real_3d(array,nlvdi,nkndi,numberNodes, 2)
        end if

        end subroutine rebuild_3d_nodes_r


!******************************************************************

!        subroutine rebuild_3d_nodes_d(array, reb_array)
!
!        use basin
!        use levels
!        use mpi_common_struct
!
!        implicit none
!
!        double precision array(nlvdi,nkn)
!        double precision, allocatable, optional, dimension(:,:) :: reb_array
!
!        if(my_id .eq. 0) then
!          if(.not. allocated(reb_array)) then
!            allocate(reb_array(nlvdi,nkndi))
!          end if
!          call rebuild_dp_3d(array,nlvdi,nkndi,numberNodes, 2,reb_array)
!        else
!          call rebuild_dp_3d(array,nlvdi,nkndi,numberNodes, 2)
!        end if
!
!        end subroutine rebuild_3d_nodes_d


!******************************************************************

        subroutine rebuild_3d_elems_r(array, reb_array)

        use basin
        use levels
        use mpi_common_struct

        implicit none

        double precision array(nlvdi,nel)
        double precision, allocatable, optional, dimension(:,:) :: reb_array

        if(my_id .eq. 0) then
          if(.not. allocated(reb_array)) then
            allocate(reb_array(nlvdi,neldi))
          end if
          call rebuild_real_3d(array,nlvdi, neldi, numberElements, 1,reb_array)
        else
          call rebuild_real_3d(array,nlvdi, neldi, numberElements, 1)
        end if

        end subroutine rebuild_3d_elems_r

!******************************************************************

        subroutine rebuild_3d_nodes3_i(array, reb_array)

        use basin
        use levels
        use mpi_common_struct

        implicit none

        integer array(3,nkn)
        integer, allocatable, optional, dimension(:,:) :: reb_array

        if(my_id .eq. 0) then
          if(.not. allocated(reb_array)) then
            allocate(reb_array(nlvdi,nkndi))
          end if
          call rebuild_integ_3d(array,3,nkndi,numberNodes, 2,reb_array)
        else
          call rebuild_integ_3d(array,3,nkndi,numberNodes, 2)
        end if

        end subroutine rebuild_3d_nodes3_i

!******************************************************************

!******************************************************************

        subroutine rebuild_structures(ilhv,znv,zenv,utlnv,vtlnv)

        use basin
        !use hydro_admin
        use levels ,only: nlvdi

        implicit none
        integer ilhv(nel)
        double precision znv(nkn)
        double precision zenv(3,nel)
        double precision utlnv(nlvdi,nel)
        double precision vtlnv(nlvdi,nel)

        if(my_id .eq. 0) then
          call rebuild_2d_nodes(znv,outZnv)
          call rebuild_3d_elems3(zenv,outZenv)
          call rebuild_2d_elems(ilhv,outIlhv)
          call rebuild_3d_elems(utlnv,outUtlnv)
          call rebuild_3d_elems(vtlnv,outVtlnv)
        else
          call rebuild_2d_nodes(znv)
          call rebuild_3d_elems3(zenv)
          call rebuild_2d_elems(ilhv)
          call rebuild_3d_elems(utlnv)
          call rebuild_3d_elems(vtlnv)
        end if

        end subroutine rebuild_structures

!******************************************************************

        subroutine rebuild_ous_header

        use basin
        use levels
        use depth

        implicit none
        
        if(my_id .eq. 0) then
          call rebuild_2d_elems(ilhv,outIlhv)
          call rebuild_2d_elems(hev,outHev)
        else
          call rebuild_2d_elems(ilhv)
          call rebuild_2d_elems(hev)
        end if

        end subroutine rebuild_ous_header

!******************************************************************

        subroutine rebuild_nos_header

        use basin
        use levels
        use depth

        implicit none
        
        if(my_id .eq. 0) then
          call rebuild_2d_nodes(ilhkv, outIlhkv)
          call rebuild_2d_elems(hev,outHev)
        else
          call rebuild_2d_nodes(ilhkv)
          call rebuild_2d_elems(hev)
        end if

        

        end subroutine rebuild_nos_header

!******************************************************************

        subroutine rebuild_scalar(dim2,array)

        implicit none

        integer dim2
        double precision array(dim2,*)
        
        character*(5) what

        if(my_id .eq. 0) then
          call rebuild_3d_nodes_r2(dim2,array,outTempv)
        else
          call rebuild_3d_nodes_r2(dim2,array)
        end if
        

        end subroutine rebuild_scalar

!******************************************************************

        subroutine compute_my_tracer(fullStruct,myStruct)

        use basin
        use levels
        use mpi_communication_struct

        implicit none

        integer i
        double precision fullStruct(nlvdi,nkndi)
        double precision myStruct(nlvdi,nkndi)

        do i=1,nkn
          myStruct(:,i) = fullStruct(:,domain%nodes%globalID(i))
        end do

        end subroutine compute_my_tracer

!******************************************************************

        subroutine shympi_ts_init(ilhkv_fem,hk_fem,he_fem,hlv_fem, &
                                        nkn_fem,nel_fem,nlv_fem)

        use levels, only: ilhkv,hlv,ilhv,nlvdi
        use depth
        use basin, only: neldi,nkndi

        implicit none
        !argument
        integer nkn_fem,nel_fem,nlv_fem
        integer, allocatable, dimension(:) :: ilhkv_fem
        double precision, allocatable, dimension (:) :: hk_fem, he_fem, hlv_fem
        double precision, allocatable, dimension (:) :: new_hlv_fem

        integer ierr, hlv_dim
        double precision time1

          allocate(inTempv(nlvdi,nkndi))
          allocate(inSaltv(nlvdi,nkndi))

          if(my_id .eq. 0) then
            call rebuild_2d_nodes(ilhkv_fem, inIlhkv)
            call rebuild_2d_nodes(hk_fem, inHkv)
            call rebuild_2d_elems(he_fem,inHev)
          else
            call rebuild_2d_nodes(ilhkv_fem)
            call rebuild_2d_nodes(hk_fem)
            call rebuild_2d_elems(he_fem)
          end if

          if(my_id .ne. 0) then
            allocate(inIlhkv(nkndi))
            allocate(inHkv(nkndi))
            allocate(inHev(neldi))
          end if

          !hlv_dim_new=size(hlv_fem)
          hlv_dim = shympi_max(size(hlv_fem))

          if(size(hlv_fem) .lt. hlv_dim) then
            deallocate(hlv_fem)
            allocate(hlv_fem(hlv_dim))
            hlv_fem = 0.
            !allocate(new_hlv_fem(hlv_dim))
          end if
          !if(.not. allocated(new_hlv_fem)) then
          allocate(new_hlv_fem(size(hlv_fem)))
          !end if

          if(ln_timing) time1 = shympi_wtime()

          call MPI_Bcast(inIlhkv(1),nkndi,MPI_INTEGER,0,MPI_COMM_WORLD ,ierr)

          call MPI_Bcast(inHkv(1),nkndi,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD ,ierr)

          call MPI_Bcast(inHev(1),neldi,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD ,ierr)

          call MPI_Allreduce(hlv_fem(1),new_hlv_fem(1),hlv_dim,MPI_DOUBLE_PRECISION, &
                                MPI_MAX,MPI_COMM_WORLD,ierr)

          if(ln_timing) comm_time = comm_time + shympi_wtime() - time1

          nkn_fem = nkndi
          nel_fem = neldi
          nlv_fem = hlv_dim


          allocate(ilhkv_fem_temp(size(ilhkv_fem)))
          ilhkv_fem_temp = ilhkv_fem
          deallocate(ilhkv_fem)
          allocate(ilhkv_fem(nkndi))
          ilhkv_fem = inIlhkv 

          allocate(hk_fem_temp(size(hk_fem)))
          hk_fem_temp = hk_fem
          deallocate(hk_fem)
          allocate(hk_fem(nkndi))
          hk_fem = inHkv

          allocate(he_fem_temp(size(he_fem)))
          he_fem_temp = he_fem
          deallocate(he_fem)
          allocate(he_fem(neldi))
          he_fem=inHev

          hlv_fem = new_hlv_fem

        return

        end subroutine

!******************************************************************

        subroutine shympi_ts_finalize(ilhkv_fem,he_fem,hk_fem,hlv_fem, nkn)

        use levels, only: ilhkv,hlv,ilhv,nlvdi
        use depth

        implicit none

        integer nkn,i
        integer, allocatable,dimension(:) :: ilhkv_fem
        double precision, allocatable, dimension(:) :: hk_fem
        double precision, allocatable, dimension(:) :: he_fem
        double precision, allocatable, dimension(:) :: hlv_fem

        deallocate(he_fem)
        allocate(he_fem(size(he_fem_temp)))
        he_fem = he_fem_temp
        deallocate(he_fem_temp)

        deallocate(hk_fem)
        allocate(hk_fem(size(hk_fem_temp)))
        hk_fem = hk_fem_temp
        deallocate(hk_fem_temp)

        deallocate(ilhkv_fem)
        allocate(ilhkv_fem(size(ilhkv_fem_temp)))
        ilhkv_fem = ilhkv_fem_temp
        deallocate(ilhkv_fem_temp)
        
        return

        end subroutine

!******************************************************************

        subroutine rebuild_integ_2d(struct, size2, sizeArray, mode, newStruct)

        use mpi_common_struct

        implicit none

        integer size2, mode
        integer :: count_tag, i
        integer :: st(MPI_STATUS_SIZE), ierr
        integer :: sendbuffer
        integer,optional, dimension(size2) :: newStruct
        integer, dimension(size2) :: struct
        integer, dimension(n_threads) :: displsProc, sizeArray, recvbuffer, displs
        integer, allocatable, dimension(:) :: fullStruct
        double precision time1

        if(my_id .eq. 0) then
          allocate(fullStruct(size2)) 
        end if

        sendbuffer = sizeArray(my_id+1)

        displs(1) = 0
        do i=2,n_threads
          displs(i) = displs(i-1) + sizeArray(i-1)
        end do

        if(ln_timing) time1 = shympi_wtime()

        call MPI_GATHERV(struct, sendbuffer, MPI_INTEGER, fullStruct, &
                sizeArray, displs, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

        if(ln_timing) comm_time = comm_time + shympi_wtime() - time1

        if(my_id .eq. 0) then
          displsProc(1) = 0
          do i=2,n_threads
            displsProc(i) = displsProc(i-1) + sizeArray(i-1)
          end do
          if(mode .eq. 1) then
            call rebuild_integ_2d_int(fullStruct, displsProc, newStruct, allPartAssign, size2)
          else
            call rebuild_integ_2d_int(fullStruct, displsProc, newStruct, univocalNodesAssign, size2)
          end if
          deallocate(fullStruct)
        end if

      end subroutine

!******************************************************************

      subroutine rebuild_integ_2d_int(fullStruct, displsProc, newStruct, itemMap, size2)

        implicit none

        integer :: place, i, process, size2
        integer :: st(MPI_STATUS_SIZE), ierr
        integer, dimension(:) :: countProc(n_threads), displsProc(n_threads)
        integer, dimension(:) :: fullStruct, newStruct
        integer itemMap(size2)

        ! countProc(1) = counter for items of the master process
        do i=1,n_threads
          countProc(i)=1
        end do

        do i=1, size2
          process = itemMap(i)
          place = displsProc(process+1) + countProc(process+1)
          newStruct(i) = fullStruct(place)
          countProc(process+1) = countProc(process+1) + 1
        end do
   
        return 

      end subroutine

!******************************************************************

      subroutine rebuild_real_2d(struct, size2, sizeArray, mode, newStruct)

        use mpi_common_struct

        implicit none

        integer size2, mode
        integer :: count_tag, i
        integer :: st(MPI_STATUS_SIZE), ierr
        integer :: sendbuffer
        double precision,optional, dimension(size2) :: newStruct
        double precision, dimension(size2) :: struct
        integer, dimension(n_threads) :: displsProc, sizeArray, recvbuffer, displs
        double precision, allocatable, dimension(:) :: fullStruct
        double precision time1

        if(my_id .eq. 0) then
          allocate(fullStruct(size2)) 
        end if

        sendbuffer = sizeArray(my_id+1)

        displs(1) = 0
        do i=2,n_threads
          displs(i) = displs(i-1) + sizeArray(i-1)
        end do

        if(ln_timing) time1 = shympi_wtime()

        call MPI_GATHERV(struct, sendbuffer, MPI_DOUBLE_PRECISION, fullStruct, &
                sizeArray, displs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

        if(ln_timing) comm_time = comm_time + shympi_wtime() - time1

        if(my_id .eq. 0) then
          displsProc(1) = 0
          do i=2,n_threads
            displsProc(i) = displsProc(i-1) + sizeArray(i-1)
          end do
          if(mode .eq. 1) then
            call rebuild_real_2d_int(fullStruct, displsProc, newStruct, allPartAssign, size2)
          else
            call rebuild_real_2d_int(fullStruct, displsProc, newStruct, univocalNodesAssign, size2)
          end if
          deallocate(fullStruct)
        end if

      end subroutine

!******************************************************************

      subroutine rebuild_real_2d_int(fullStruct, displsProc, newStruct, itemMap, size2)

        implicit none

        integer :: place, i, process, size2
        integer :: st(MPI_STATUS_SIZE), ierr
        integer, dimension(:) :: countProc(n_threads), displsProc(n_threads)
        double precision, dimension(:) :: fullStruct, newStruct
        integer itemMap(size2)

        ! countProc(1) = counter for items of the master process
        do i=1,n_threads
          countProc(i)=1
        end do

        do i=1, size2
          process = itemMap(i)
          place = displsProc(process+1) + countProc(process+1)
          newStruct(i) = fullStruct(place)
          countProc(process+1) = countProc(process+1) + 1
        end do
   
        return 

      end subroutine

!******************************************************************

      subroutine rebuild_dp_2d(struct, size2, sizeArray, mode, newStruct)

        use mpi_common_struct

        implicit none

        integer size2, mode
        integer :: count_tag, i
        integer :: st(MPI_STATUS_SIZE), ierr
        integer :: sendbuffer
        double precision,optional, dimension(size2) :: newStruct
        double precision, dimension(size2) :: struct
        integer, dimension(n_threads) :: displsProc, sizeArray, recvbuffer, displs
        double precision, allocatable, dimension(:) :: fullStruct
        double precision time1

        if(my_id .eq. 0) then
          allocate(fullStruct(size2)) 
        end if

        sendbuffer = sizeArray(my_id+1)

        displs(1) = 0
        do i=2,n_threads
          displs(i) = displs(i-1) + sizeArray(i-1)
        end do

        if(ln_timing) time1 = shympi_wtime()

        call MPI_GATHERV(struct, sendbuffer, MPI_DOUBLE_PRECISION, &
            fullStruct,sizeArray, displs, MPI_DOUBLE_PRECISION, 0, &
            MPI_COMM_WORLD, ierr)

        if(ln_timing) comm_time = comm_time + shympi_wtime() - time1

        if(my_id .eq. 0) then
          displsProc(1) = 0
          do i=2,n_threads
            displsProc(i) = displsProc(i-1) + sizeArray(i-1)
          end do
          if(mode .eq. 1) then
            call rebuild_dp_2d_int(fullStruct, displsProc, newStruct, allPartAssign, size2)
          else
            call rebuild_dp_2d_int(fullStruct, displsProc, newStruct, univocalNodesAssign, size2)
          end if
          deallocate(fullStruct)
        end if

      end subroutine

!******************************************************************

      subroutine rebuild_dp_2d_int(fullStruct, displsProc, newStruct, itemMap, size2)

        implicit none

        integer :: place, i, process, size2
        integer :: st(MPI_STATUS_SIZE), ierr
        integer, dimension(:) :: countProc(n_threads), displsProc(n_threads)
        double precision, dimension(:) :: fullStruct, newStruct
        integer itemMap(size2)

        ! countProc(1) = counter for items of the master process
        do i=1,n_threads
          countProc(i)=1
        end do

        do i=1, size2
          process = itemMap(i)
          place = displsProc(process+1) + countProc(process+1)
          newStruct(i) = fullStruct(place)
          countProc(process+1) = countProc(process+1) + 1
        end do
   
        return 

      end subroutine

!******************************************************************

      subroutine rebuild_integ_3d(struct, size1, size2, sizeArray, mode, newStruct)

        use mpi_common_struct

        implicit none

        integer size1, size2, mode
        integer :: count_tag, i
        integer :: st(MPI_STATUS_SIZE), ierr
        integer :: sendbuffer
        integer, optional, dimension(size1,size2) :: newStruct
        integer, dimension(size1,size2) :: struct
        integer, dimension(n_threads) :: displsProc, sizeArray, recvbuffer, displs
        integer, allocatable, dimension(:,:) :: fullStruct
        double precision time1

        if(my_id .eq. 0) then
          allocate(fullStruct(size1,size2)) 
          do i=1,n_threads
            recvbuffer(i) = sizeArray(i) * size1
          end do
        end if

        sendbuffer = size1 * sizeArray(my_id+1)

        displs(1) = 0
        do i=2,n_threads
          displs(i) = displs(i-1) + size1 * sizeArray(i-1)
        end do

        if(ln_timing) time1 = shympi_wtime()

        call MPI_GATHERV(struct, sendbuffer, MPI_INTEGER, fullStruct, &
                recvbuffer, displs, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

        if(ln_timing) comm_time = comm_time + shympi_wtime() - time1

        if(my_id .eq. 0) then
          displsProc(1) = 0
          do i=2,n_threads
            displsProc(i) = displsProc(i-1) + sizeArray(i-1)
          end do
          if(mode .eq. 1) then
            call rebuild_integ_3d_int(fullStruct,displsProc,newStruct,allPartAssign,size1,size2)
          else
            call rebuild_integ_3d_int(fullStruct,displsProc,newStruct,univocalNodesAssign,size1,size2)
          end if
          deallocate(fullStruct)
        end if

      end subroutine

!******************************************************************

      subroutine rebuild_integ_3d_int(fullStruct, displsProc, newStruct, itemMap, size1, size2)

        implicit none

        integer :: place, i, process, j, size1, size2
        integer :: st(MPI_STATUS_SIZE), ierr
        integer, dimension(:) :: countProc(n_threads), displsProc(n_threads)
        integer, dimension(:,:) :: fullStruct, newStruct
        integer itemMap(size2)
    
        ! countProc(1) = counter for items of the master process
        do i=1,n_threads
          countProc(i)=1
        end do

        do i=1, size2
          process = itemMap(i)
          place = displsProc(process+1) + countProc(process+1)
          do j=1,size1
            newStruct(j,i) = fullStruct(j,place)
          end do
          countProc(process+1) = countProc(process+1) + 1
        end do

      end subroutine

!******************************************************************

      subroutine rebuild_real_3d(struct, size1, size2, sizeArray, mode, newStruct)

        use mpi_common_struct

        implicit none

        integer size1, size2, mode
        integer :: count_tag, i
        integer :: st(MPI_STATUS_SIZE), ierr
        integer :: sendbuffer
        double precision, optional, dimension(size1,size2) :: newStruct
        double precision, dimension(size1,size2) :: struct
        integer, dimension(n_threads) :: displsProc, sizeArray, recvbuffer, displs
        double precision, allocatable, dimension(:,:) :: fullStruct
        double precision time1

        if(my_id .eq. 0) then
          allocate(fullStruct(size1,size2)) 
          do i=1,n_threads
            recvbuffer(i) = sizeArray(i) * size1
          end do
        end if

        sendbuffer = size1 * sizeArray(my_id+1)

        displs(1) = 0
        do i=2,n_threads
          displs(i) = displs(i-1) + size1 * sizeArray(i-1)
        end do

        if(ln_timing) time1 = shympi_wtime()

        call MPI_GATHERV(struct, sendbuffer, MPI_DOUBLE_PRECISION, fullStruct, &
                recvbuffer, displs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

        if(ln_timing) comm_time = comm_time + shympi_wtime() - time1

        if(my_id .eq. 0) then
          displsProc(1) = 0
          do i=2,n_threads
            displsProc(i) = displsProc(i-1) + sizeArray(i-1)
          end do
          if(mode .eq. 1) then
            call rebuild_real_3d_int(fullStruct,displsProc,newStruct,allPartAssign,size1,size2)
          else
            call rebuild_real_3d_int(fullStruct,displsProc,newStruct,univocalNodesAssign,size1,size2)
          end if
          deallocate(fullStruct)
        end if

      end subroutine

!******************************************************************

      subroutine rebuild_real_3d_int(fullStruct, displsProc, newStruct, itemMap, size1, size2)

        implicit none

        integer :: place, i, process, j, size1, size2
        integer :: st(MPI_STATUS_SIZE), ierr
        integer, dimension(:) :: countProc(n_threads), displsProc(n_threads)
        double precision, dimension(:,:) :: fullStruct, newStruct
        integer itemMap(size2)
    
        ! countProc(1) = counter for items of the master process
        do i=1,n_threads
          countProc(i)=1
        end do

        do i=1, size2
          process = itemMap(i)
          place = displsProc(process+1) + countProc(process+1)
          do j=1,size1
            newStruct(j,i) = fullStruct(j,place)
          end do
          countProc(process+1) = countProc(process+1) + 1
        end do

      end subroutine

!******************************************************************

      subroutine rebuild_dp_3d(struct, size1, size2, sizeArray, mode, newStruct)

        use mpi_common_struct

        implicit none

        integer size1, size2, mode
        integer :: count_tag, i
        integer :: st(MPI_STATUS_SIZE), ierr
        integer :: sendbuffer
        double precision, optional, dimension(size1,size2) :: newStruct
        double precision, dimension(size1,size2) :: struct
        integer, dimension(n_threads) :: displsProc, sizeArray, recvbuffer, displs
        double precision, allocatable, dimension(:,:) :: fullStruct
        double precision time1

        if(my_id .eq. 0) then
          allocate(fullStruct(size1,size2)) 
          do i=1,n_threads
            recvbuffer(i) = sizeArray(i) * size1
          end do
        end if

        sendbuffer = size1 * sizeArray(my_id+1)

        displs(1) = 0
        do i=2,n_threads
          displs(i) = displs(i-1) + size1 * sizeArray(i-1)
        end do

        if(ln_timing) time1 = shympi_wtime()

        call MPI_GATHERV(struct, sendbuffer, MPI_DOUBLE_PRECISION, &
          fullStruct, recvbuffer, displs, MPI_DOUBLE_PRECISION, 0, &
          MPI_COMM_WORLD, ierr)

        if(ln_timing) comm_time = comm_time + shympi_wtime() - time1

        if(my_id .eq. 0) then
          displsProc(1) = 0
          do i=2,n_threads
            displsProc(i) = displsProc(i-1) + sizeArray(i-1)
          end do
          if(mode .eq. 1) then
            call rebuild_dp_3d_int(fullStruct,displsProc,newStruct,allPartAssign,size1,size2)
          else
            call rebuild_dp_3d_int(fullStruct,displsProc,newStruct,univocalNodesAssign,size1,size2)
          end if
          deallocate(fullStruct)
        end if

      end subroutine

!******************************************************************

      subroutine rebuild_dp_3d_int(fullStruct, displsProc, newStruct, itemMap, size1, size2)

        implicit none

        integer :: place, i, process, j, size1, size2
        integer :: st(MPI_STATUS_SIZE), ierr
        integer, dimension(:) :: countProc(n_threads), displsProc(n_threads)
        double precision, dimension(:,:) :: fullStruct, newStruct
        integer itemMap(size2)
    
        ! countProc(1) = counter for items of the master process
        do i=1,n_threads
          countProc(i)=1
        end do

        do i=1, size2
          process = itemMap(i)
          place = displsProc(process+1) + countProc(process+1)
          do j=1,size1
            newStruct(j,i) = fullStruct(j,place)
          end do
          countProc(process+1) = countProc(process+1) + 1
        end do

      end subroutine

!******************************************************************

        subroutine rebuild_ev(array, reb_array)

        use basin
        use mpi_common_struct
        use evgeom

        implicit none

        double precision array(evdim,nel)
        double precision, allocatable, optional, dimension(:,:) :: reb_array

        if(my_id .eq. 0) then
          if(.not. allocated(reb_array)) then
            allocate(reb_array(evdim,neldi))
          end if
          call rebuild_ev_3d(array,evdim, neldi, numberElements,1,reb_array)
        else
          call rebuild_ev_3d(array,evdim, neldi, numberElements,1)
        end if

        end subroutine rebuild_ev

!******************************************************************

      subroutine rebuild_ev_3d(struct, size1, size2, sizeArray, mode, newStruct)

        use mpi_common_struct

        implicit none

        integer size1, size2, mode
        integer :: count_tag, i
        integer :: st(MPI_STATUS_SIZE), ierr
        integer :: sendbuffer
        double precision, optional, dimension(size1,size2) :: newStruct
        double precision, dimension(size1,size2) :: struct
        integer, dimension(n_threads) :: displsProc, sizeArray, recvbuffer, displs
        double precision, allocatable, dimension(:,:) :: fullStruct
        double precision time1

        if(my_id .eq. 0) then
          allocate(fullStruct(size1,size2)) 
          do i=1,n_threads
            recvbuffer(i) = sizeArray(i) * size1
          end do
        end if

        sendbuffer = size1 * sizeArray(my_id+1)

        displs(1) = 0
        do i=2,n_threads
          displs(i) = displs(i-1) + size1 * sizeArray(i-1)
        end do

        if(ln_timing) time1 = shympi_wtime()

        call MPI_GATHERV(struct, sendbuffer, MPI_DOUBLE_PRECISION      &
                , fullStruct, recvbuffer, displs, MPI_DOUBLE_PRECISION &
                , 0, MPI_COMM_WORLD, ierr)

        if(ln_timing) comm_time = comm_time + shympi_wtime() - time1

        if(my_id .eq. 0) then
          displsProc(1) = 0
          do i=2,n_threads
            displsProc(i) = displsProc(i-1) + sizeArray(i-1)
          end do
          if(mode .eq. 1) then
            call rebuild_ev_3d_int(fullStruct,displsProc,newStruct,allPartAssign,size1,size2)
          else
            call rebuild_ev_3d_int(fullStruct,displsProc,newStruct,univocalNodesAssign,size1,size2)
          end if
          deallocate(fullStruct)
        end if

      end subroutine

!******************************************************************

      subroutine rebuild_ev_3d_int(fullStruct, displsProc, newStruct, itemMap, size1, size2)

        implicit none

        integer :: place, i, process, j, size1, size2
        integer :: st(MPI_STATUS_SIZE), ierr
        integer, dimension(:) :: countProc(n_threads), displsProc(n_threads)
        double precision, dimension(:,:) :: fullStruct, newStruct
        integer itemMap(size2)
    
        ! countProc(1) = counter for items of the master process
        do i=1,n_threads
          countProc(i)=1
        end do

        do i=1, size2
          process = itemMap(i)
          place = displsProc(process+1) + countProc(process+1)
          do j=1,size1
            newStruct(j,i) = fullStruct(j,place)
          end do
          countProc(process+1) = countProc(process+1) + 1
        end do

      end subroutine

!**************************************************************
!**************************************************************
! opening of default simulation
!**************************************************************
!**************************************************************
!**************************************************************

        function ifemopp(ext,form,status)

! opens file with default name (run) and extension given for fem model
! returns with error code

! ext   extension (with dot)
! form  formatted ?
! status open status

        implicit none

	integer ifemopp
        character*(*) ext,status,form

	character*80 file,defdir,defnam

        call def_makep(ext,file)

	ifemopp=ifileop(0,file,form,status)

	end function

!**************************************************************
!**************************************************************
!**************************************************************

        subroutine def_makep(ext,file)

! makes file with defaults supplied
!
! ext   extension (with dot)
! file  created file name (return)

        use para
        use defnames
        implicit none

        character*(*) ext,file
        character*80 dir,name

	name = def_nam
        call getfnm('datdir',dir)	! this has to be deleted
        call getfnm('runnam',name)	! this has to be deleted

	call mknamep(dir,name,ext,file)

	end subroutine

!**************************************************************

        subroutine mknamep(dir,name,ext,file)

! makes file name given its constituents
!
! dir   directory
! name  name
! ext   extension (with or without dot)
! file  created file name (return)

        use utility
        use defnames
        implicit none

! arguments
        character*(*) dir,name,ext,file
! local
        integer nall,nstart,nend,naux

        nall=1
        file=' '

        nstart=ichafs(dir)
        nend=ichanm(dir)
        if(nend.gt.0) then
		file(nall:)=dir(nstart:nend)
        	nall=nall+nend-nstart+1
	end if

        nstart=ichafs(name)
        nend=ichanm(name)
        if(nend.gt.0) then
		file(nall:)=name(nstart:nend)
       		nall=nall+nend-nstart+1
	end if

        call add_rank(file)

	call add_extension(file,ext,.false.)

	return

	end subroutine

!***************************************************************

	function ifileop(iunit,file,format,status)

! opens file
!
! iunit		unit number to be opened
! file		file name to be opened
! format	f*ormatted, u*nformatted or b*inary
! status	o*ld, n*ew, s*cratch  or u*nknown
! ifileo	effective unit number the file is opened
!
! routine tries to open file at iunit. If this is not possible
! ...higher unit numbers are tested and the file is eventually
! ...opened at one of these units. The effective unit number
! ...is passed back in ifileo. -1 is returned if no file is opened.
! use iunit <= 0 to use standard unit number
! iunit < 0 does not complain if not existing

        use fil

	implicit none

	integer ifileop
	integer iunit
	character*(*) file,format,status

	integer iu,nfile,ios
	logical found,error,opened,ex,od,bquiet
	character*1 cf,cs
	character*15 form,stat,access

        integer iustd
        save iustd
        data iustd / 20 /

!----------------------------------------------------------------
! initialize basic things
!----------------------------------------------------------------

	ifileop=-1
	iu=iunit
	bquiet = ( iu < 0 )		!if unit negative do not complain
	if( iu .le. 0 ) iu = iustd + my_id

!----------------------------------------------------------------
! check key words
!----------------------------------------------------------------

	cf=format(1:1)
	cs=status(1:1)

	access = 'sequential'

	if( cf .eq. 'f' .or. cf .eq. 'F' ) then
		form='formatted'
	else if( cf .eq. 'u' .or. cf .eq. 'U' ) then
		form='unformatted'
	else if( cf .eq. 'b' .or. cf .eq. 'B' ) then
		form='unformatted'
!lahey#		access = 'transparent'
	else
		write(6,*) 'format keyword not recognized :',format
		return
	end if

	if( cs .eq. 'o' .or. cs .eq. 'O' ) then
		stat='old'
	else if( cs .eq. 'n' .or. cs .eq. 'N' ) then
! for VAX      change to	stat='new'
! for DOS/UNIX change to	stat='unknown'
!		stat='new'
		stat='unknown'
	else if( cs .eq. 's' .or. cs .eq. 'S' ) then
		stat='scratch'
	else if( cs .eq. 'u' .or. cs .eq. 'U' ) then
		stat='unknown'
	else
		write(6,*) 'status keyword not recognized :',status
		return
	end if

!----------------------------------------------------------------
! check if file exists (in case of status=old)
!----------------------------------------------------------------

	inquire(file=file,exist=ex)
	if(.not.ex.and.stat.eq.'old') then
	  if( .not. bquiet ) then
	    write(6,*) 'file does not exist : ',trim(file)
	  end if
	  return
	end if

!----------------------------------------------------------------
! find unit where to open file
!----------------------------------------------------------------

	call find_unit(iu)
	found=iu.gt.0

!----------------------------------------------------------------
! open file and check error
!----------------------------------------------------------------

	if(found) then
          	open(	 unit=iu &
     			,file=file &
     			,form=form &
     			,status=stat &
     			,access=access &
     			,iostat=ios &
     	            )
	  if(ios.ne.0) then
	    nfile=ichanm0(file)
	    if(nfile.le.0) nfile=1
	        write(6,*) 'error opening file : ',file(1:nfile)
		write(6,*) 'unit : ',iu,'  iostat : ',ios
		write(6,*) 'error : ',mod(ios,256)
		inquire(file=file,opened=opened)
		if( opened ) then
		  write(6,*) '...the file is already open...'
		  write(6,*) 'If you are using gfortran or pgf90'
		  write(6,*) 'please remember that you can open'
		  write(6,*) 'a file only once. You will have to'
		  write(6,*) 'copy the file to files with different'
		  write(6,*) 'names and open these files instead'
		end if
	  else
		rewind(iu)
		ifileop=iu
                if( iunit .le. 0 ) iustd=iu ! remember for next time
                !write(10,*) 'ifileo: ',iu,iunit,'  ',file(1:40)
	  end if
	end if

!----------------------------------------------------------------
! end of routine
!----------------------------------------------------------------

	end function


!**************************************************************

        subroutine add_rank(file)

        implicit none

	character*80 file
        character*(80) filename
        character*(20) format_string
        character*10 rank_str
        integer length

        if(my_id .le. 9) then
           write(rank_str,'(A4,I1)'),'_000',my_id
        else if(my_id .le. 99) then
           write(rank_str,'(A3,I2)'),'_00',my_id
        else if(my_id .le. 999) then
           write(rank_str,'(A2,I3)'),'_0',my_id
        else
           write(rank_str,'(A1,I4)'),'_',my_id
        end if 

        length=len_trim(file)
        write(format_string,'(A2,I2,A5)'),'(A',length,',A5)'

        write(filename,format_string),file,rank_str

        file=filename
 

        end subroutine

!******************************************************************

        subroutine adjust_name(file)

        use defnames
        implicit none

        character*80 file 

        call cutext(file)

        call add_rank(file)

	call add_extension(file,'rst',.false.)

        return

        end subroutine

!******************************************************************

        subroutine cutext(file)

        implicit none

        integer offset
        character*80 file,filename,myformat


        offset=index(file,'.')
          offset = offset-1
        if(offset .lt. 10) then
          write(myformat,'(A2,I1,A1)'),'(A',offset,')'
        else
          write(myformat,'(A2,I2,A1)'),'(A',offset,')'
        end if
        write(filename,myformat),file

        file=filename

        return

        end subroutine

!******************************************************************

        subroutine up_ielt_mpi(inodv)

        use basin
        use links

        implicit none

        integer ii,ie,ierr
        integer, allocatable, dimension(:) :: ginodv
        integer, allocatable, dimension(:,:) :: gnen3v
        integer, allocatable, dimension(:,:) :: temp_nen3v
        integer inodv(nkn)
        double precision time1

!-------------------------------------------------------------
! prepare update ieltv mpi
!-------------------------------------------------------------

        if(bmpi) then
          allocate(gnen3v(3,nel))
          do ie=1,nel
            do ii=1,3
                gnen3v(ii,ie)=domain%nodes%globalID(nen3v(ii,ie))
            end do
          end do

          if(my_id.eq.0) then
            allocate(temp_nen3v(3,nel_local))
            temp_nen3v=nen3v
            deallocate(nen3v)
            allocate(nen3v(3,neldi))
            call rebuild_3d_elems3(gnen3v,nen3v)
            call rebuild_2d_nodes_i(inodv,ginodv)
          else
            call rebuild_3d_elems3(gnen3v)
            call rebuild_2d_nodes_i(inodv)
          end if
        end if

        if(bmpi.and.(my_id.eq.0)) then                
          call update_ielt(neldi,ginodv,total_ieltv)
          deallocate(nen3v)
          allocate(nen3v(3,nel_local))
          nen3v=temp_nen3v
        end if

        if(ln_timing) time1 = shympi_wtime()

        if (bmpi) call MPI_BCAST(total_ieltv, 3*neldi, MPI_INTEGER, 0, &
                        MPI_COMM_WORLD,ierr)

        if(ln_timing) comm_time = comm_time + shympi_wtime() - time1

        return

        end subroutine

!==================================================================
        end module mpi_io_admin
!==================================================================

