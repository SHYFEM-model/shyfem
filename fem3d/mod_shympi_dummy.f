
!--------------------------------------------------------------------------
!
!    Copyright (C) 2015,2017-2019  Georg Umgiesser
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

! mpi dummy routines
!
! contents :
!
! revision log :
!
! 24.11.2015	ggu	project started
! 05.12.2017	ggu	changed VERS_7_5_39
! 07.12.2017	ggu	changed VERS_7_5_40
! 24.01.2018	ggu	changed VERS_7_5_41
! 10.04.2018	ggu	adjourned with new function calls
! 19.04.2018	ggu	changed VERS_7_5_45
! 26.04.2018	ggu	changed VERS_7_5_46
! 11.05.2018	ggu	new function calls for globals and internal pointers
! 06.07.2018	ggu	changed VERS_7_5_48
! 16.02.2019	ggu	changed VERS_7_5_60
! 10.06.2020	ggu	adjournments from node integrated
! 09.04.2021    clr     bug fix in shympi_bcast_array_r() -> real arg
! 17.04.2021    clr     new shympi_exchange_array_3(), check_external_numbers()
! 23.04.2021    clr     formal change in MODULE PROCEDURE declarations for meson compatibility
! 08.06.2021    ggu     parameters in shympi_exchange_array_3() were integer
! 25.06.2021    ggu     in shympi_init() check if basin has been read
! 20.10.2021    ggu     do not stop when reading grd file
!
!******************************************************************

!==================================================================
        module shympi
!==================================================================

	implicit none

	public

	logical, save :: bmpi = .false.
	logical, save :: bmpi_debug = .false.

	integer,save :: n_threads = 1
	integer,save :: my_id = 0
	integer,save :: my_unit = 0

        integer,save :: status_size = 0          !ngr of total basin

        integer,save :: ngr_global = 0          !ngr of total basin
        integer,save :: nlv_global = 0          !nlv of total basin

	integer,save :: nkn_global = 0		!total basin
	integer,save :: nel_global = 0
	integer,save :: nkn_local = 0		!this domain
	integer,save :: nel_local = 0
	integer,save :: nkn_unique = 0		!this domain unique
	integer,save :: nel_unique = 0
	integer,save :: nkn_inner = 0		!only proper, no ghost
	integer,save :: nel_inner = 0

        integer,save :: nk_max = 0              !max of nkn of all domains
        integer,save :: ne_max = 0              !max of nel of all domains
        integer,save :: nn_max = 0              !max of nkn/nel of all domains

	integer,save :: n_ghost_areas = 0
	integer,save :: n_ghost_nodes_max = 0
	integer,save :: n_ghost_elems_max = 0
	integer,save :: n_ghost_max = 0
	integer,save :: n_buffer = 0

        integer,save,allocatable :: nkn_domains(:)
        integer,save,allocatable :: nel_domains(:)
        integer,save,allocatable :: nkn_cum_domains(:)
        integer,save,allocatable :: nel_cum_domains(:)

	integer,save,allocatable :: ghost_areas(:,:)
	integer,save,allocatable :: ghost_nodes_in(:,:)
	integer,save,allocatable :: ghost_nodes_out(:,:)
	integer,save,allocatable :: ghost_elems_in(:,:)
	integer,save,allocatable :: ghost_elems_out(:,:)

	!integer,save,allocatable :: i_buffer_in(:,:)
	!integer,save,allocatable :: i_buffer_out(:,:)
	!real,save,allocatable    :: r_buffer_in(:,:)
	!real,save,allocatable    :: r_buffer_out(:,:)
	
	integer,save,allocatable :: id_node(:)
	integer,save,allocatable :: id_elem(:,:)

        integer,save,allocatable :: ip_sort_node(:)     !pointer to sorted node
        integer,save,allocatable :: ip_sort_elem(:)     !pointer to sorted elem

        integer,save,allocatable :: ip_ext_node(:)      !global external nums
        integer,save,allocatable :: ip_ext_elem(:)
        integer,save,allocatable :: ip_int_node(:)      !global internal nums
        integer,save,allocatable :: ip_int_elem(:)
        integer,save,allocatable :: ip_int_nodes(:,:)   !all global int nums
        integer,save,allocatable :: ip_int_elems(:,:)
        integer,save,allocatable :: nen3v_global(:,:)
	real,save,allocatable :: hlv_global(:)

        type communication_info
          integer, public :: numberID
          integer, public :: totalID
          !contains the global ID of the elements belonging to the process
          integer, public, dimension(:), allocatable :: globalID,localID
          !contains the rank ID of the process to which belong the global ID
          !elements with which to communicate
          integer, public, dimension(:,:), allocatable :: rankID
          integer, public, dimension(:,:), allocatable :: neighborID
        end type communication_info

	type (communication_info), SAVE :: univocal_nodes
        integer, allocatable, save, dimension(:) :: allPartAssign

!---------------------

        INTERFACE shympi_exchange_3d_node
        MODULE PROCEDURE  shympi_exchange_3d_node_r
     +                   ,shympi_exchange_3d_node_d
     +                   ,shympi_exchange_3d_node_i
        END INTERFACE

        INTERFACE shympi_exchange_3d0_node
       	MODULE PROCEDURE   shympi_exchange_3d0_node_r
!     +                   ,shympi_exchange_3d0_node_d
!     +                   ,shympi_exchange_3d0_node_i
        END INTERFACE

        INTERFACE shympi_exchange_2d_node
        MODULE PROCEDURE  shympi_exchange_2d_node_r
     +                   ,shympi_exchange_2d_node_d
     +                   ,shympi_exchange_2d_node_i
        END INTERFACE

        INTERFACE shympi_exchange_3d_elem
        MODULE PROCEDURE   shympi_exchange_3d_elem_r
!     +                   ,shympi_exchange_3d_elem_d
!     +                   ,shympi_exchange_3d_elem_i
        END INTERFACE

        INTERFACE shympi_exchange_2d_elem
        MODULE PROCEDURE  shympi_exchange_2d_elem_r
     +                   ,shympi_exchange_2d_elem_d
     +                   ,shympi_exchange_2d_elem_i
        END INTERFACE

!---------------------

        INTERFACE shympi_check_elem
        MODULE PROCEDURE  shympi_check_2d_elem_r
     +                   ,shympi_check_2d_elem_d
     +                   ,shympi_check_2d_elem_i
     +			 ,shympi_check_3d_elem_r
!     +                   ,shympi_check_3d_elem_d
!     +                   ,shympi_check_3d_elem_i
        END INTERFACE

        INTERFACE shympi_check_node
        MODULE PROCEDURE  shympi_check_2d_node_r
     +                   ,shympi_check_2d_node_d
     +                   ,shympi_check_2d_node_i
     +			 ,shympi_check_3d_node_r
!     +                   ,shympi_check_3d_node_d
!     +                   ,shympi_check_3d_node_i
        END INTERFACE

        INTERFACE shympi_check_2d_node
        MODULE PROCEDURE  shympi_check_2d_node_r
     +                   ,shympi_check_2d_node_d
     +                   ,shympi_check_2d_node_i
        END INTERFACE

        INTERFACE shympi_check_2d_elem
        MODULE PROCEDURE  shympi_check_2d_elem_r
     +                   ,shympi_check_2d_elem_d
     +                   ,shympi_check_2d_elem_i
        END INTERFACE

        INTERFACE shympi_check_3d_node
        MODULE PROCEDURE  shympi_check_3d_node_r
!     +                   ,shympi_check_3d_node_d
!     +                   ,shympi_check_3d_node_i
        END INTERFACE

        INTERFACE shympi_check_3d0_node
        MODULE PROCEDURE  shympi_check_3d0_node_r
!     +                   ,shympi_check_3d0_node_d
!     +                   ,shympi_check_3d0_node_i
        END INTERFACE

        INTERFACE shympi_check_3d_elem
        MODULE PROCEDURE  shympi_check_3d_elem_r
!     +                   ,shympi_check_3d_elem_d
!     +                   ,shympi_check_3d_elem_i
        END INTERFACE

        INTERFACE shympi_check_array
        MODULE PROCEDURE  shympi_check_array_i
     +                   ,shympi_check_array_r
     +                   ,shympi_check_array_d
        END INTERFACE

!---------------------

        INTERFACE shympi_gather
        MODULE PROCEDURE  shympi_gather_scalar_i
     +                   ,shympi_gather_array_2d_i
     +                   ,shympi_gather_array_2d_r
     +                   ,shympi_gather_array_2d_d
        END INTERFACE

        INTERFACE shympi_gather_and_sum
        MODULE PROCEDURE  shympi_gather_and_sum_i
     +                   ,shympi_gather_and_sum_r
     +                   ,shympi_gather_and_sum_d
        END INTERFACE

        INTERFACE shympi_bcast
        MODULE PROCEDURE  shympi_bcast_scalar_i
     +                   ,shympi_bcast_array_r
        END INTERFACE

        INTERFACE shympi_collect_node_value
        MODULE PROCEDURE   shympi_collect_node_value_2d_i
     +                    ,shympi_collect_node_value_2d_r
     +                    ,shympi_collect_node_value_3d_r
!     +                    ,shympi_collect_node_value_2d_i
        END INTERFACE

        INTERFACE shympi_reduce
        MODULE PROCEDURE shympi_reduce_r
!     +                   ,shympi_reduce_i
        END INTERFACE

!---------------------

        INTERFACE shympi_min
        MODULE PROCEDURE   shympi_min_r
     +			  ,shympi_min_i
     +			  ,shympi_min_d
     +			  ,shympi_min_0_r
     +			  ,shympi_min_0_i
     +			  ,shympi_min_0_d
        END INTERFACE

        INTERFACE shympi_max
        MODULE PROCEDURE   shympi_max_r
     +			  ,shympi_max_i
!     +			  ,shympi_max_d
     +			  ,shympi_max_0_r
     +			  ,shympi_max_0_i
!     +			  ,shympi_max_0_d
        END INTERFACE

        INTERFACE shympi_sum
        MODULE PROCEDURE   shympi_sum_r
     +			  ,shympi_sum_i
     +			  ,shympi_sum_d
     +			  ,shympi_sum_0_r
     +			  ,shympi_sum_0_i
     +			  ,shympi_sum_0_d
        END INTERFACE

!---------------------

        INTERFACE shympi_exchange_array
        MODULE PROCEDURE   shympi_exchange_array_2d_r
     +			  ,shympi_exchange_array_2d_i
     +			  ,shympi_exchange_array_3d_r
     +			  ,shympi_exchange_array_3d_i
        END INTERFACE

        INTERFACE shympi_get_array
        MODULE PROCEDURE   shympi_get_array_2d_r
     +			  ,shympi_get_array_2d_i
!     +			  ,shympi_get_array_3d_r
!     +			  ,shympi_get_array_3d_i
        END INTERFACE

        INTERFACE shympi_getvals
        MODULE PROCEDURE   shympi_getvals_2d_node_r
     +			  ,shympi_getvals_2d_node_i
     +			  ,shympi_getvals_3d_node_r
     +			  ,shympi_getvals_3d_node_i
        END INTERFACE

!---------------------

        INTERFACE shympi_exchange_and_sum_3d_nodes
        MODULE PROCEDURE   shympi_exchange_and_sum_3d_nodes_r
     +			  ,shympi_exchange_and_sum_3d_nodes_d
        END INTERFACE

        INTERFACE shympi_exchange_and_sum_2d_nodes
        MODULE PROCEDURE   shympi_exchange_and_sum_2d_nodes_r
     +			  ,shympi_exchange_and_sum_2d_nodes_d
        END INTERFACE

        INTERFACE shympi_exchange_2d_nodes_min
        MODULE PROCEDURE  shympi_exchange_2d_nodes_min_i
     +			  ,shympi_exchange_2d_nodes_min_r
        END INTERFACE

        INTERFACE shympi_exchange_2d_nodes_max
        MODULE PROCEDURE   shympi_exchange_2d_nodes_max_i
     +			  ,shympi_exchange_2d_nodes_max_r
        END INTERFACE

!==================================================================
        contains
!==================================================================

	subroutine shympi_init(b_want_mpi)

	use basin

	logical b_want_mpi

	integer ierr,size
	character*10 cunit
	character*80 file

        !-----------------------------------------------------
        ! initializing (next call returns n_threads == 1)
        !-----------------------------------------------------

	call shympi_init_internal(my_id,n_threads)
	!call check_part_basin('nodes')

	if( .not. basin_has_read_basin() ) then
	  write(6,*) 'grd file has been read: ',nkn,nel,ngr
	  if( nkn == 0 ) then
	    stop 'error stop shympi_init: basin has not been initialized'
	  end if
	end if

	ngr_global = ngr

	nkn_global = nkn
	nel_global = nel
	nkn_local = nkn
	nel_local = nel
	nkn_inner = nkn
	nel_inner = nel
	nkn_unique = nkn
	nel_unique = nel

	bmpi = n_threads > 1		!this is always false

	if( b_want_mpi ) then
	 if( shympi_is_master() ) then
          !write(6,*) 'mpi thread: ',my_id
          write(6,*) 'program wants mpi but only one thread available'
	  write(6,*) 'the program has not been compiled with mpi support'
	 end if
         call shympi_stop('error stop shympi_init')
	end if

        !-----------------------------------------------------
        ! allocate important arrays
        !-----------------------------------------------------

	call shympi_get_status_size_internal(size)

        allocate(nkn_domains(n_threads))
        allocate(nel_domains(n_threads))
        allocate(nkn_cum_domains(0:n_threads))
        allocate(nel_cum_domains(0:n_threads))

        nkn_domains(1) = nkn
        nel_domains(1) = nel
	nk_max = nkn
	ne_max = nel
	nn_max = max(nkn,nel)

	!write(79,*) 'domains: '
     +	!	,nkn_domains,nel_domains,nk_max,ne_max,nn_max
	!stop

        nkn_cum_domains(0) = 0
        nkn_cum_domains(1) = nkn
        nel_cum_domains(0) = 0
        nel_cum_domains(1) = nel

	call shympi_alloc_global(nkn,nel,nen3v,ipv,ipev)

        !-----------------------------------------------------
        ! next is needed if program is not running in mpi mode
        !-----------------------------------------------------

	if( .not. bmpi ) then
	  call shympi_alloc_id(nkn,nel)
          call shympi_alloc_sort(nkn,nel)
	  call mpi_sort_index(nkn,nel)
	end if

        !-----------------------------------------------------
        ! output to terminal
        !-----------------------------------------------------

	if( bmpi ) then
	  write(cunit,'(i10)') my_id
	  cunit = adjustl(cunit)
	  file = 'mpi_debug_' // trim(cunit) // '.txt'
	  call shympi_get_new_unit(my_unit)
	  open(unit=my_unit,file=file,status='unknown')
	  !open(newunit=my_unit,file=file,status='unknown')
	  write(my_unit,*) 'shympi initialized: ',my_id,n_threads
	else if( bmpi_debug ) then
          write(6,*) 'shympi initialized: ',my_id,n_threads
          write(6,*) 'shympi is not running in mpi mode'
	end if

	flush(6)

        !-----------------------------------------------------
        ! finish up
        !-----------------------------------------------------

	call shympi_barrier_internal
	call shympi_syncronize_initial
	call shympi_syncronize_internal

        !-----------------------------------------------------
        ! end of routine
        !-----------------------------------------------------

	end subroutine shympi_init

!******************************************************************

	subroutine shympi_alloc_id(nk,ne)

	integer nk,ne

	!write(6,*) 'shympi_alloc_id: ',nk,ne

	allocate(id_node(nk))
	allocate(id_elem(0:2,ne))

	id_node = my_id
	id_elem = my_id

	end subroutine shympi_alloc_id

!******************************************************************

        subroutine shympi_alloc_sort(nk,ne)

        integer nk,ne

	if( allocated(ip_sort_node) ) deallocate(ip_sort_node)
	if( allocated(ip_sort_elem) ) deallocate(ip_sort_elem)

!	write(6,*) 'allocating sort arrays: ',nk,ne

        allocate(ip_sort_node(nk))
        allocate(ip_sort_elem(ne))

!	write(6,*) 'allocated sort arrays: '
!     +			,size(ip_sort_node),size(ip_sort_elem)

        ip_sort_node = 0
        ip_sort_elem = 0

        end subroutine shympi_alloc_sort

!******************************************************************

        subroutine shympi_alloc_global(nk,ne,nen3v,ipv,ipev)

        integer nk,ne
        integer nen3v(3,ne)
	integer ipv(nk)
	integer ipev(ne)

        integer i

        allocate(ip_ext_node(nk))
        allocate(ip_ext_elem(ne))
        allocate(ip_int_node(nk))
        allocate(ip_int_elem(ne))
        allocate(ip_int_nodes(nk,1))
        allocate(ip_int_elems(ne,1))
        allocate(nen3v_global(3,ne))

        do i=1,nk
          ip_int_node(i) = i
          ip_int_nodes(i,1) = i
        end do

        do i=1,ne
          ip_int_elem(i) = i
          ip_int_elems(i,1) = i
        end do

        nen3v_global = nen3v
        ip_ext_node = ipv
        ip_ext_elem = ipev

        end subroutine shympi_alloc_global

!******************************************************************

	subroutine shympi_alloc_ghost(n)

	use basin

	integer n

	allocate(ghost_areas(5,n_ghost_areas))
        allocate(ghost_nodes_out(n,n_ghost_areas))
        allocate(ghost_nodes_in(n,n_ghost_areas))
        allocate(ghost_elems_out(n,n_ghost_areas))
        allocate(ghost_elems_in(n,n_ghost_areas))

	ghost_areas = 0
        ghost_nodes_out = 0
        ghost_nodes_in = 0
        ghost_elems_out = 0
        ghost_elems_in = 0

	end subroutine shympi_alloc_ghost

!******************************************************************

        subroutine shympi_alloc_buffer(n)

        integer n
	integer,save,allocatable :: i_buffer_in(:,:)
	integer,save,allocatable :: i_buffer_out(:,:)
	real,save,allocatable    :: r_buffer_in(:,:)
	real,save,allocatable    :: r_buffer_out(:,:)

        if( n_buffer >= n ) return

        if( n_buffer > 0 ) then
          deallocate(i_buffer_in)
          deallocate(i_buffer_out)
          deallocate(r_buffer_in)
          deallocate(r_buffer_out)
        end if

        n_buffer = n

        allocate(i_buffer_in(n_buffer,n_ghost_areas))
        allocate(i_buffer_out(n_buffer,n_ghost_areas))

        allocate(r_buffer_in(n_buffer,n_ghost_areas))
        allocate(r_buffer_out(n_buffer,n_ghost_areas))

        end subroutine shympi_alloc_buffer

!******************************************************************
!******************************************************************
!******************************************************************

        function shympi_internal_node(k)

        integer shympi_internal_node
        integer k
        integer i

        do i=1,nkn_global
          if( ip_ext_node(i) == k ) exit
        end do
        if( i > nkn_global ) i = 0

        shympi_internal_node = i

        end function shympi_internal_node

!******************************************************************

        function shympi_internal_elem(ie)

        integer shympi_internal_elem
        integer ie
        integer i

        do i=1,nel_global
          if( ip_ext_elem(i) == ie ) exit
        end do
        if( i > nel_global ) i = 0

        shympi_internal_elem = i

        end function shympi_internal_elem

!******************************************************************

        function shympi_exist_node(k)

        logical shympi_exist_node
        integer k

        shympi_exist_node = ( shympi_internal_node(k) > 0 )

        end function shympi_exist_node

!******************************************************************

        function shympi_exist_elem(ie)

        logical shympi_exist_elem
        integer ie

        shympi_exist_elem = ( shympi_internal_elem(ie) > 0 )

        end function shympi_exist_elem

!******************************************************************
!******************************************************************
!******************************************************************

        subroutine shympi_set_hlv(nlv,hlv)

	integer nlv
	real hlv(nlv)

	if( .not. allocated(hlv_global) ) allocate(hlv_global(nlv))

	nlv_global = nlv
	hlv_global = hlv

        end subroutine shympi_set_hlv

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_get_new_unit(iunit)

	integer iunit

	integer iu,iumin,iumax,iostat
        logical opened

	return

	iumin = 20
	iumax = 1000

	do iu=iumin,iumax
          inquire (unit=iu, opened=opened, iostat=iostat)
          if (iostat.ne.0) cycle
          if (.not.opened) exit
	end do

	if( iu > iumax ) then
	  iu = 0
	  stop 'error stop shympi_get_new_unit: no new unit'
	end if

	iunit = iu

	end subroutine shympi_get_new_unit

!******************************************************************
!******************************************************************
!******************************************************************

	function shympi_partition_on_elements()

	logical shympi_partition_on_elements

	shympi_partition_on_elements = .false.

	end function shympi_partition_on_elements

!******************************************************************

        function shympi_partition_on_nodes()

        logical shympi_partition_on_nodes

        shympi_partition_on_nodes = .false.

        end function shympi_partition_on_nodes

!******************************************************************
!******************************************************************
!******************************************************************

	function shympi_is_inner_node(k)

	logical shympi_is_inner_node
	integer k

	shympi_is_inner_node = ( k <= nkn_inner )

	end function shympi_is_inner_node

!******************************************************************

	function shympi_is_inner_elem(ie)

	logical shympi_is_inner_elem
	integer ie

	shympi_is_inner_elem = ( ie <= nel_inner )

	end function shympi_is_inner_elem

!******************************************************************

	function shympi_is_unique_node(k)

	logical shympi_is_unique_node
	integer k

	shympi_is_unique_node = ( k <= nkn_unique )

	end function shympi_is_unique_node

!******************************************************************

	function shympi_is_unique_elem(ie)

	logical shympi_is_unique_elem
	integer ie

	shympi_is_unique_elem = ( ie <= nel_unique )

	end function shympi_is_unique_elem

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_barrier
	end subroutine shympi_barrier

!******************************************************************

	subroutine shympi_syncronize
	end subroutine shympi_syncronize

!******************************************************************

	subroutine shympi_finalize
	end subroutine shympi_finalize

!******************************************************************

	subroutine shympi_exit(ierr)

	integer ierr

	call exit(ierr)
	stop

	end subroutine shympi_exit

!******************************************************************

	subroutine shympi_stop(text)

	character*(*) text

	write(6,*) 'error stop shympi_stop: ',trim(text)
	stop

	end subroutine shympi_stop

!******************************************************************

	subroutine shympi_abort

	call exit(33)
	stop

	end subroutine shympi_abort

!******************************************************************

	function shympi_wtime()

	double precision shympi_wtime
	double precision shympi_wtime_internal

	shympi_wtime = 0.

	end function shympi_wtime

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_exchange_3d_node_i(val)

	use basin
	use levels

	integer val(nlvdi,nkn)
	logical, parameter :: belem = .false.

	end subroutine shympi_exchange_3d_node_i

!*******************************

	subroutine shympi_exchange_3d_node_r(val)

	use basin
	use levels

	real val(nlvdi,nkn)
	logical, parameter :: belem = .false.

	end subroutine shympi_exchange_3d_node_r

!*******************************

	subroutine shympi_exchange_3d_node_d(val)

	use basin
	use levels

	double precision val(nlvdi,nkn)
	logical, parameter :: belem = .false.

	end subroutine shympi_exchange_3d_node_d

!******************************************************************

	subroutine shympi_exchange_3d0_node_r(val)

	use basin
	use levels

	real val(0:nlvdi,nkn)
	logical, parameter :: belem = .false.

	end subroutine shympi_exchange_3d0_node_r

!******************************************************************

	subroutine shympi_exchange_3d_elem_r(val)

	use basin
	use levels

	real val(nlvdi,nel)
	logical, parameter :: belem = .true.

	end subroutine shympi_exchange_3d_elem_r

!******************************************************************

	subroutine shympi_exchange_2d_node_i(val)

	use basin
	use levels

	integer val(nkn)
	logical, parameter :: belem = .false.

	end subroutine shympi_exchange_2d_node_i

!*******************************

	subroutine shympi_exchange_2d_node_r(val)

	use basin
	use levels

	real val(nkn)
	logical, parameter :: belem = .false.

	end subroutine shympi_exchange_2d_node_r

!*******************************

	subroutine shympi_exchange_2d_node_d(val)

	use basin
	use levels

	double precision val(nkn)
	logical, parameter :: belem = .false.

	end subroutine shympi_exchange_2d_node_d

!******************************************************************

	subroutine shympi_exchange_2d_elem_i(val)

	use basin
	use levels

	integer val(nel)
	logical, parameter :: belem = .true.

	end subroutine shympi_exchange_2d_elem_i

!*******************************

	subroutine shympi_exchange_2d_elem_r(val)

	use basin
	use levels

	real val(nel)
	logical, parameter :: belem = .true.

	end subroutine shympi_exchange_2d_elem_r

!*******************************

	subroutine shympi_exchange_2d_elem_d(val)

	use basin
	use levels

	double precision val(nel)
	logical, parameter :: belem = .true.

	end subroutine shympi_exchange_2d_elem_d

!******************************************************************
!******************************************************************
!******************************************************************

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_check_2d_node_i(val,text)

	use basin

	integer val(nkn)
	character*(*) text

	end subroutine shympi_check_2d_node_i

!*******************************

	subroutine shympi_check_2d_node_r(val,text)

	use basin

	real val(nkn)
	character*(*) text

	end subroutine shympi_check_2d_node_r

!*******************************

	subroutine shympi_check_2d_node_d(val,text)

	use basin

	double precision val(nkn)
	character*(*) text

	end subroutine shympi_check_2d_node_d

!******************************************************************

	subroutine shympi_check_3d_node_r(val,text)

	use basin
	use levels

	real val(nlvdi,nkn)
	character*(*) text

	end subroutine shympi_check_3d_node_r

!******************************************************************

	subroutine shympi_check_3d0_node_r(val,text)

	use basin
	use levels

	real val(0:nlvdi,nkn)
	character*(*) text

	end subroutine shympi_check_3d0_node_r

!******************************************************************

	subroutine shympi_check_2d_elem_i(val,text)

	use basin

	integer val(nel)
	character*(*) text

	end subroutine shympi_check_2d_elem_i

!*******************************

	subroutine shympi_check_2d_elem_r(val,text)

	use basin

	real val(nel)
	character*(*) text

	end subroutine shympi_check_2d_elem_r

!*******************************

	subroutine shympi_check_2d_elem_d(val,text)

	use basin

	double precision val(nel)
	character*(*) text

	end subroutine shympi_check_2d_elem_d

!******************************************************************

	subroutine shympi_check_3d_elem_r(val,text)

	use basin
	use levels

	real val(nlvdi,nel)
	character*(*) text

	end subroutine shympi_check_3d_elem_r

!******************************************************************

        subroutine shympi_check_array_i(n,a1,a2,text)

        integer n
        integer a1(n),a2(n)
        character*(*) text

        end subroutine shympi_check_array_i

!******************************************************************

        subroutine shympi_check_array_r(n,a1,a2,text)

        integer n
        real a1(n),a2(n)
        character*(*) text

        end subroutine shympi_check_array_r

!******************************************************************

        subroutine shympi_check_array_d(n,a1,a2,text)

        integer n
        double precision a1(n),a2(n)
        character*(*) text

        end subroutine shympi_check_array_d

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_gather_scalar_i(val,vals)

	integer val
	integer vals(n_threads)

	vals(1) = val

	end subroutine shympi_gather_scalar_i

!*******************************

        subroutine shympi_gather_array_2d_i(val,vals)

        integer val(:)
        integer vals(size(val),n_threads)

	vals(:,1) = val(:)

        end subroutine shympi_gather_array_2d_i

!*******************************

        subroutine shympi_gather_array_2d_r(val,vals)

        real val(:)
        real vals(size(val),n_threads)

	vals(:,1) = val(:)

        end subroutine shympi_gather_array_2d_r

!*******************************

        subroutine shympi_gather_array_2d_d(val,vals)

        double precision val(:)
        double precision vals(size(val),n_threads)

	vals(:,1) = val(:)

        end subroutine shympi_gather_array_2d_d

!*******************************

        subroutine shympi_gather_and_sum_i(val)

        integer val(:)

        end subroutine shympi_gather_and_sum_i

!*******************************

        subroutine shympi_gather_and_sum_r(val)

        real val(:)

        end subroutine shympi_gather_and_sum_r

!*******************************

        subroutine shympi_gather_and_sum_d(val)

        double precision val(:)

        end subroutine shympi_gather_and_sum_d

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_bcast_scalar_i(val)

	integer val

	end subroutine shympi_bcast_scalar_i

!*******************************

        subroutine shympi_bcast_array_r(val)

        real val(:)

        end subroutine shympi_bcast_array_r

!*******************************

        subroutine shympi_collect_node_value_2d_r(k,vals,val)

        integer k
        real vals(:)
        real val

        val = 0
        if( k > 0 .and. k <= nkn_unique ) val = vals(k)

        end subroutine shympi_collect_node_value_2d_r

!*******************************

        subroutine shympi_collect_node_value_2d_i(k,vals,val)

        integer k
        integer vals(:)
        integer val

        val = 0
        if( k > 0 .and. k <= nkn_unique ) val = vals(k)

        end subroutine shympi_collect_node_value_2d_i

!*******************************

        subroutine shympi_collect_node_value_3d_r(k,vals,val)

        integer k
        real vals(:,:)
        real val(:)

        real vaux(size(val),n_threads)

        val = 0
        if( k > 0 .and. k <= nkn_unique ) val(:) = vals(:,k)

        end subroutine shympi_collect_node_value_3d_r

!*******************************

!       subroutine shympi_collect_node_value_3d_i(k,vals,val)
!
!       integer k
!       integer vals(:)
!       integer val
!
!       val = 0
!       if( k > 0 .and. k <= nkn_unique ) val = vals(k)
!       if( bmpi ) val = shympi_sum(val)
!
!       end subroutine shympi_collect_node_value_3d_i

!*******************************

        subroutine shympi_find_node(ke,ki,id)

        use basin

        integer ke,ki,id

        integer k,ic,i

        do k=1,nkn_unique
          if( ipv(k) == ke ) exit
        end do
	if( k > nkn_unique ) k = 0

        ki = k
        id = 0

        end subroutine shympi_find_node

!*******************************

	subroutine shympi_reduce_r(what,vals,val)

	character*(*) what
	real vals(:)
	real val

	if( what == 'min' ) then
	  val = MINVAL(vals)
	else if( what == 'max' ) then
	  val = MAXVAL(vals)
	else if( what == 'sum' ) then
	  val = SUM(vals)
	else
	  write(6,*) 'what = ',what
	  stop 'error stop shympi_reduce_r: not ready'
	end if

	end subroutine shympi_reduce_r

!******************************************************************
!******************************************************************
!******************************************************************

        subroutine shympi_get_array_2d_r(n,vals,val_out)

        use basin
        use levels

        integer n
        real vals(:)
        real val_out(n)

	val_out = vals

        end subroutine shympi_get_array_2d_r

!*******************************

        subroutine shympi_get_array_2d_i(n,vals,val_out)

        use basin
        use levels

        integer n
        integer vals(:)
        integer val_out(n)

	val_out = vals

        end subroutine shympi_get_array_2d_i

!******************************************************************
!******************************************************************
!******************************************************************

        subroutine shympi_exchange_array_3d_r(vals,val_out)

        real vals(:,:)
        real val_out(:,:)

	val_out = vals

        end subroutine shympi_exchange_array_3d_r

!*******************************

        subroutine shympi_exchange_array_3d_i(vals,val_out)

        integer vals(:,:)
        integer val_out(:,:)

	val_out = vals

        end subroutine shympi_exchange_array_3d_i

!*******************************

        subroutine shympi_exchange_array_3(vals,val_out)

        real vals(:,:)
        real val_out(:,:)

	val_out = vals

        end subroutine shympi_exchange_array_3

!*******************************

        subroutine shympi_exchange_array_2d_r(vals,val_out)

        real vals(:)
        real val_out(:)

	val_out = vals

        end subroutine shympi_exchange_array_2d_r

!*******************************

        subroutine shympi_exchange_array_2d_i(vals,val_out)

        integer vals(:)
        integer val_out(:)

	val_out = vals

        end subroutine shympi_exchange_array_2d_i

!******************************************************************
!******************************************************************
!******************************************************************

        subroutine shympi_getvals_2d_node_r(kind,vals,val)

        use basin
        use levels

        integer kind(2)
        real vals(nkn)
        real val

	val = vals(kind(1))

        end subroutine shympi_getvals_2d_node_r

!*******************************

        subroutine shympi_getvals_2d_node_i(kind,vals,val)

        use basin
        use levels

        integer kind(2)
        integer vals(nkn)
        integer val

	val = vals(kind(1))

        end subroutine shympi_getvals_2d_node_i

!*******************************

        subroutine shympi_getvals_3d_node_r(kind,vals,val)

        use basin
        use levels

        integer kind(2)
        real vals(nlvdi,nkn)
        real val(nlvdi)

	val(:) = vals(:,kind(1))

        end subroutine shympi_getvals_3d_node_r

!*******************************

        subroutine shympi_getvals_3d_node_i(kind,vals,val)

        use basin
        use levels

        integer kind(2)
        integer vals(nlvdi,nkn)
        integer val(nlvdi)

	val(:) = vals(:,kind(1))

        end subroutine shympi_getvals_3d_node_i

!******************************************************************
!******************************************************************
!******************************************************************

	function shympi_min_d(vals)

	double precision shympi_min_d
	double precision vals(:)
	double precision val

	val = MINVAL(vals)

	shympi_min_d = val

	end function shympi_min_d

!******************************************************************

	function shympi_min_r(vals)

	real shympi_min_r
	real vals(:)
	real val

	val = MINVAL(vals)

	shympi_min_r = val

	end function shympi_min_r

!******************************************************************

	function shympi_min_i(vals)

	integer shympi_min_i
	integer vals(:)
	integer val

	val = MINVAL(vals)

	shympi_min_i = val

	end function shympi_min_i

!******************************************************************

	function shympi_min_0_d(val)

	double precision shympi_min_0_d
	double precision val

	shympi_min_0_d = val

	end function shympi_min_0_d

!******************************************************************

	function shympi_min_0_r(val)

	real shympi_min_0_r
	real val

	shympi_min_0_r = val

	end function shympi_min_0_r

!******************************************************************

	function shympi_min_0_i(val)

	integer shympi_min_0_i
	integer val

	shympi_min_0_i = val

	end function shympi_min_0_i

!******************************************************************

	function shympi_max_r(vals)

	real shympi_max_r
	real vals(:)
	real val

	val = MAXVAL(vals)

	shympi_max_r = val

	end function shympi_max_r

!******************************************************************

	function shympi_max_i(vals)

	integer shympi_max_i
	integer vals(:)
	integer val

	val = MAXVAL(vals)

	shympi_max_i = val

	end function shympi_max_i

!******************************************************************

	function shympi_max_0_i(val)

! routine for val that is scalar

	integer shympi_max_0_i
	integer val

	shympi_max_0_i = val

	end function shympi_max_0_i

!******************************************************************

	function shympi_max_0_r(val)

! routine for val that is scalar

	real shympi_max_0_r
	real val

	shympi_max_0_r = val

	end function shympi_max_0_r

!******************************************************************

        function shympi_sum_d(vals)

        double precision shympi_sum_d
        double precision vals(:)
        double precision val

        val = SUM(vals)

        shympi_sum_d = val

        end function shympi_sum_d

!******************************************************************

	function shympi_sum_r(vals)

	real shympi_sum_r
	real vals(:)
	real val

	val = SUM(vals)

	shympi_sum_r = val

	end function shympi_sum_r

!******************************************************************

	function shympi_sum_i(vals)

	integer shympi_sum_i
	integer vals(:)
	integer val

	val = SUM(vals)

	shympi_sum_i = val

	end function shympi_sum_i

!******************************************************************

        function shympi_sum_0_d(val)

        double precision shympi_sum_0_d
        double precision val

        shympi_sum_0_d = val

        end function shympi_sum_0_d

!******************************************************************

	function shympi_sum_0_r(val)

	real shympi_sum_0_r
	real val

	shympi_sum_0_r = val

	end function shympi_sum_0_r

!******************************************************************

	function shympi_sum_0_i(val)

	integer shympi_sum_0_i
	integer val

	shympi_sum_0_i = val

	end function shympi_sum_0_i

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_exchange_and_sum_3d_nodes_r(a)
	real a(:,:)
	end subroutine shympi_exchange_and_sum_3d_nodes_r

	subroutine shympi_exchange_and_sum_3d_nodes_d(a)
	double precision a(:,:)
	end subroutine shympi_exchange_and_sum_3d_nodes_d

	subroutine shympi_exchange_and_sum_2d_nodes_r(a)
	real a(:)
	end subroutine shympi_exchange_and_sum_2d_nodes_r

	subroutine shympi_exchange_and_sum_2d_nodes_d(a)
	double precision a(:)
	end subroutine shympi_exchange_and_sum_2d_nodes_d

	subroutine shympi_exchange_2d_nodes_min_i(a)
	integer a(:)
	end subroutine shympi_exchange_2d_nodes_min_i

	subroutine shympi_exchange_2d_nodes_max_i(a)
	integer a(:)
	end subroutine shympi_exchange_2d_nodes_max_i

	subroutine shympi_exchange_2d_nodes_min_r(a)
	real a(:)
	end subroutine shympi_exchange_2d_nodes_min_r

	subroutine shympi_exchange_2d_nodes_max_r(a)
	real a(:)
	end subroutine shympi_exchange_2d_nodes_max_r

!******************************************************************
!******************************************************************
!******************************************************************

	function shympi_output()

	logical shympi_output

	shympi_output = my_id == 0

	end function shympi_output

!******************************************************************

	function shympi_is_master()

	logical shympi_is_master

	shympi_is_master = my_id == 0

	end function shympi_is_master

!******************************************************************

	subroutine shympi_comment(text)

	character*(*) text

	!if( bmpi .and. bmpi_debug .and. my_id == 0 ) then
	if( bmpi_debug .and. my_id == 0 ) then
	  write(6,*) 'shympi_comment: ' // trim(text)
	  write(299,*) 'shympi_comment: ' // trim(text)
	end if

	end subroutine shympi_comment

!******************************************************************

	subroutine shympi_parallel_code(text)

	character*(*) text

	text = 'serial'

	end subroutine shympi_parallel_code

!******************************************************************

	function shympi_is_parallel()

	logical shympi_is_parallel

	shympi_is_parallel = bmpi

	end function shympi_is_parallel

!******************************************************************

        function shympi_can_parallel()

        logical shympi_can_parallel

        shympi_can_parallel = .false.

        end function shympi_can_parallel

!******************************************************************

	subroutine shympi_univocal_nodes

	use basin, only : nkn

	implicit none

	integer ierr

	integer i

	univocal_nodes%numberID = nkn

	allocate(univocal_nodes%localID(nkn))

	do i=1,nkn
	  univocal_nodes%localID(i) = i
	end do

	end subroutine shympi_univocal_nodes

!******************************************************************

        subroutine check_part_basin(what)

        use basin

        implicit none

        character*(5) what
        integer pnkn,pnel,pn_threads,ierr
        integer control
        character*(14) filename
        character*(11) pwhat
        integer i

        call shympi_get_filename(filename,what)

        write(6,*) filename,what

        open(unit=108, file=filename, form="formatted"
     +   , iostat=control, status="old", action="read")
        
        if(control .ne. 0) then
          if(my_id .eq. 0) then
            write(6,*)'error stop: partitioning file not found'
          end if
          call shympi_barrier
        stop
        end if

        read(unit=108, fmt="(i12,i12,i12,A12)") pnkn,pnel,pn_threads
     +                  ,pwhat

        if(pnkn .ne. nkndi .or. pnel .ne. neldi .or. pn_threads
     &          .ne. n_threads .or. pwhat .ne. what) then
         if(my_id .eq. 0) then
          write(6,*)'basin file does not match'
          write(6,*)'partitioning file:nkn,nel,n_threads,partitioning'
     +          ,pnkn,pnel,pn_threads,pwhat
          write(6,*)'basin in str file:nkn,nel,n_threads,partitioning'
     +          ,nkndi,neldi,n_threads,what
         end if
         call shympi_barrier
         stop
        end if

        if(what .eq. 'nodes') then
          allocate(allPartAssign(nkndi))
          read(unit=108,fmt="(i12,i12,i12,i12,i12,i12)")
     +          (allPartAssign(i),i=1,nkndi)
        else 
          write(6,*)'error partitioning file on nodes'
          stop
        end if

        close(108)

        return

        end subroutine check_part_basin        

!******************************************************************

        subroutine shympi_get_filename(filename,what)

          implicit none

          character*(5) what
          character*(14) filename,format_string

          if(n_threads .gt. 999) then
            format_string = "(A5,A5,I4)"
            write(filename,format_string)'part_',what,n_threads
          else if(n_threads .gt. 99) then
            format_string = "(A5,A5,I3)"
            write(filename,format_string)'part_',what,n_threads
          else if(n_threads .gt. 9) then
            format_string = "(A5,A5,I2)"
            write(filename,format_string)'part_',what,n_threads
          else
            format_string = "(A5,A5,I1)"
            write(filename,format_string)'part_',what,n_threads
          end if

          return

        end subroutine shympi_get_filename

!******************************************************************
!******************************************************************
!******************************************************************

        subroutine check_external_numbers
        implicit none
	end subroutine

!==================================================================
        end module shympi
!==================================================================

!******************************************************************
!******************************************************************
!******************************************************************
! next are internal routines - needed to be able to compile
!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_init_internal(my_id,n_threads)
	implicit none
	integer ierr
	integer my_id,n_threads
	my_id = 0
	n_threads = 1
	end subroutine shympi_init_internal

!******************************************************************

	subroutine shympi_barrier_internal
	end subroutine shympi_barrier_internal

!******************************************************************

	subroutine shympi_syncronize_internal
	end subroutine shympi_syncronize_internal

!******************************************************************

	subroutine shympi_finalize_internal
	end subroutine shympi_finalize_internal

!******************************************************************

	subroutine shympi_abort_internal(ierr)
	implicit none
	integer ierr
	end subroutine shympi_abort_internal

!******************************************************************

        subroutine shympi_get_status_size_internal(size)
        implicit none
        integer size
        size = 1
        end subroutine shympi_get_status_size_internal

!******************************************************************

	subroutine shympi_syncronize_initial
	end subroutine shympi_syncronize_initial

!******************************************************************

        function shympi_wtime_internal()
        implicit none
        double precision shympi_wtime_internal
        shympi_wtime_internal = 0.
        end function shympi_wtime_internal

!******************************************************************
!******************************************************************
!******************************************************************

        subroutine shympi_allgather_i_internal(n,no,val,vals)
        use shympi
        implicit none
	integer n,no
        integer val(n)
        integer vals(no,1)
	if(n/=no) stop 'error stop shympi_allgather_internal: n/=no'
	vals(:,1) = val
        end subroutine shympi_allgather_i_internal

!******************************************************************

        subroutine shympi_allgather_r_internal(n,no,val,vals)
        use shympi
        implicit none
	integer n,no
        real val(n)
        real vals(no,1)
	if(n/=no) stop 'error stop shympi_allgather_internal: n/=no'
	vals(:,1) = val
        end subroutine shympi_allgather_r_internal

!******************************************************************

        subroutine shympi_allgather_d_internal(n,no,val,vals)
        use shympi
        implicit none
	integer n,no
        double precision val(n)
        double precision vals(no,1)
	if(n/=no) stop 'error stop shympi_allgather_internal: n/=no'
	vals(:,1) = val
        end subroutine shympi_allgather_d_internal

!******************************************************************

        subroutine shympi_bcast_i_internal(n,val)
        implicit none
	integer n
        integer val
        end subroutine shympi_bcast_i_internal

!******************************************************************
!******************************************************************
!******************************************************************

        subroutine shympi_reduce_r_internal(what,val)
        implicit none
        character*(*) what
        real val
        end subroutine shympi_reduce_r_internal

!******************************************************************

        subroutine shympi_reduce_i_internal(what,val)
        implicit none
        character*(*) what
        integer val
        end subroutine shympi_reduce_i_internal


!******************************************************************

        subroutine shympi_ex_3d_nodes_sum_r_internal(array)
        implicit none
        real array(:,:)
        end subroutine

!******************************************************************

        subroutine shympi_ex_3d_nodes_sum_d_internal(array)
        implicit none
        double precision array(:,:)
        end subroutine

!******************************************************************
        subroutine shympi_ex_2d_nodes_sum_r_internal(array)
        implicit none
        real array(:)
        end subroutine

!******************************************************************

        subroutine shympi_ex_2d_nodes_sum_d_internal(array)
        implicit none
        double precision array(:)
        end subroutine

!******************************************************************

        subroutine shympi_ex_2d_nodes_min_i_internal(array)
        implicit none
        integer array(:)
        end subroutine

!******************************************************************

        subroutine shympi_ex_2d_nodes_min_r_internal(array)
        implicit none
        real array(:)
        end subroutine

!******************************************************************

        subroutine shympi_ex_2d_nodes_max_i_internal(array)
        implicit none
        integer array(:)
        end subroutine

!******************************************************************

        subroutine shympi_ex_2d_nodes_max_r_internal(array)
        implicit none
        real array(:)
        end subroutine

!******************************************************************
!******************************************************************
!******************************************************************

!******************************************************************
!******************************************************************
!******************************************************************

