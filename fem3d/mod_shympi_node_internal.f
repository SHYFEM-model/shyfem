
!--------------------------------------------------------------------------
!
!    Copyright (C) 2015-2019  Georg Umgiesser
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

! mpi routines
!
! contents :
!
! revision log :
!
! 24.11.2015	ggu	project started
! 22.06.2016	ggu	added sum option to shympi_reduce
! 07.12.2017	ggu	changed VERS_7_5_40
! 24.01.2018	ggu	changed VERS_7_5_41
! 10.04.2018	ggu	code to exchange arrays
! 19.04.2018	ggu	changed VERS_7_5_45
! 26.04.2018	ggu	changed VERS_7_5_46
! 11.05.2018	ggu	bug fix in exchange arrays for zeta levels
! 16.02.2019	ggu	changed VERS_7_5_60
!
!******************************************************************

!==================================================================
        module shympi_aux
!==================================================================

! this module is only used inside this file

        use mpi

        implicit none

        !include 'mpif.h'

	integer, save :: n_p_threads = 1	!total number of threads
	integer, save :: my_p_id = 0		!id of this thread
	logical, save :: bpmaster = .true.	!is this the master?
	logical, save :: bpmpi = .false.	!use mpi? (threads > 1)

!==================================================================
        end module shympi_aux
!==================================================================

!******************************************************************

	subroutine shympi_error(routine,what,ierr)

	use shympi_aux

	implicit none

	character*(*) routine,what
	integer ierr

	integer eslen,iserr
	character*80 estring

	if( ierr /= 0 ) then
	  eslen = 0
	  iserr = 0
	  estring = ' '
	  !call MPI_ERROR_STRING(ierr,estring,eslen,iserr)
	  if( eslen > 80 .or. iserr /= 0 ) then
	    estring = 'no description for error'
	  end if
	  write(6,*) 'error in routine ',trim(routine)
	  write(6,*) 'ierr = ',ierr,' while doing ',trim(what)
	  write(6,*) 'error: ',trim(estring)
	  stop 'error stop shympi_error'
	end if

	end subroutine shympi_error

!******************************************************************

	subroutine shympi_init_internal(my_id,n_threads)

	use shympi_aux

	implicit none

	integer my_id,n_threads

	integer ierr,iberr
	integer required,provided

	required = MPI_THREAD_MULTIPLE
	required = MPI_THREAD_SERIALIZED

        !call MPI_INIT( ierr )
        call MPI_INIT_THREAD( required, provided, ierr )
	!write(6,*) 'thread safety: ',required, provided
	!write(6,*) 'initializing MPI: ',ierr

	call MPI_BARRIER( MPI_COMM_WORLD, iberr)
	call shympi_error('shympi_init_internal','init',ierr)
        call MPI_COMM_RANK( MPI_COMM_WORLD, my_id, ierr )
	call shympi_error('shympi_init_internal','rank',ierr)
        call MPI_COMM_SIZE( MPI_COMM_WORLD, n_threads, ierr )
	call shympi_error('shympi_init_internal','size',ierr)
	call MPI_BARRIER( MPI_COMM_WORLD, ierr)
	call shympi_error('shympi_init_internal','barrier',ierr)

	n_p_threads = n_threads
	my_p_id = my_id
	bpmaster = ( my_p_id == 0 )
	bpmpi = ( n_p_threads > 1 )

        !call shympi_finalize_internal
	!stop
	!write(6,*) 'MPI internally initialized: ',my_id,n_threads

	end subroutine shympi_init_internal

!******************************************************************

	subroutine shympi_barrier_internal

	use shympi_aux

	implicit none

	integer ierr

	call MPI_BARRIER( MPI_COMM_WORLD, ierr)

	end subroutine shympi_barrier_internal

!******************************************************************

        subroutine shympi_abort_internal(ierr_code)

	use shympi_aux

        implicit none

        integer ierr,ierr_code

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

	subroutine shympi_syncronize_initial

	use shympi_aux

        implicit none

	integer my_id,nt
	integer root,ierr,i
	integer count
	integer local
	integer, allocatable :: buf(:)

        call MPI_COMM_RANK( MPI_COMM_WORLD, my_id, ierr )
        call MPI_COMM_SIZE( MPI_COMM_WORLD, nt, ierr )

	allocate(buf(nt))

	root = my_id
	count = 1
	root = 0

	if( my_id == root ) then
	  do i=1,nt
	    buf(i) = i
	  end do
	end if

	call MPI_SCATTER (buf,count,MPI_INTEGER
     +			,local,count,MPI_INTEGER
     +			,root,MPI_COMM_WORLD,ierr)

	local = local * 2
	root = 0

	call MPI_GATHER (local,count,MPI_INTEGER
     +			,buf,count,MPI_INTEGER
     +			,root,MPI_COMM_WORLD,ierr)

	call MPI_BARRIER( MPI_COMM_WORLD, ierr)

	if( my_id == root ) then
	  !write(6,*) 'mpi sync: ',nt,root,(buf(i),i=1,nt)
	  !write(6,*) 'mpi sync: ',nt,root
	end if

	deallocate(buf)

	end subroutine shympi_syncronize_initial

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_exchange_internal_i(belem,n0,nlvddi,n,il
     +						,g_in,g_out,val)

	use shympi_aux

	use shympi

	implicit none

	logical belem
	integer n0,nlvddi,n
	integer il(n)
	integer g_in(n_ghost_max,n_ghost_areas)
	integer g_out(n_ghost_max,n_ghost_areas)
	integer val(n0:nlvddi,n)

	integer tag,ir,ia,id
	integer i,k,nc,ierr
	integer nb
	integer iout,iin
	integer nbs(2,n_ghost_areas)
	integer status(status_size,2*n_threads)
	integer request(2*n_threads)
	integer, allocatable :: buffer_in(:,:)
	integer, allocatable :: buffer_out(:,:)

        tag=121
	ir = 0

	iout = 2
	iin = 3
	if( belem ) then
	  iout = 4
	  iin = 5
	end if

	nb = (nlvddi-n0+1) * n_ghost_max
	call shympi_alloc_buffer(nb)
	allocate(buffer_in(nb,n_ghost_areas))
	allocate(buffer_out(nb,n_ghost_areas))

	do ia=1,n_ghost_areas
	  nc = ghost_areas(iout,ia)
	  call count_buffer(n0,nlvddi,n,nc,il,g_out(:,ia),nb)
	  nbs(1,ia) = nb
	  nc = ghost_areas(iin,ia)
	  call count_buffer(n0,nlvddi,n,nc,il,g_in(:,ia),nb)
	  nbs(2,ia) = nb
	end do

cccgguccc!$OMP CRITICAL

	do ia=1,n_ghost_areas
	  ir = ir + 1
	  id = ghost_areas(1,ia)
	  nc = ghost_areas(iout,ia)
	  nb = nbs(1,ia)
          call MPI_Irecv(buffer_out(:,ia),nb,MPI_INTEGER,id
     +	          ,tag,MPI_COMM_WORLD,request(ir),ierr)
	end do

	do ia=1,n_ghost_areas
	  ir = ir + 1
	  id = ghost_areas(1,ia)
	  nc = ghost_areas(iin,ia)
	  nb = nbs(2,ia)
	  call to_buffer_i(n0,nlvddi,n,nc,il
     +		,g_in(:,ia),val,nb,buffer_in(:,ia))
          call MPI_Isend(buffer_in(:,ia),nb,MPI_INTEGER,id
     +	          ,tag,MPI_COMM_WORLD,request(ir),ierr)
	end do

        call MPI_WaitAll(ir,request,status,ierr)

cccgguccc!$OMP END CRITICAL

	do ia=1,n_ghost_areas
	  id = ghost_areas(1,ia)
	  nc = ghost_areas(iout,ia)
	  nb = nbs(1,ia)
	  call from_buffer_i(n0,nlvddi,n,nc,il
     +		,g_out(:,ia),val,nb,buffer_out(:,ia))
	end do

	end subroutine shympi_exchange_internal_i

!******************************************************************

	subroutine shympi_exchange_internal_r(belem,n0,nlvddi,n,il
     +						,g_in,g_out,val)

	use shympi_aux

	use shympi

	implicit none

	logical belem
	integer n0,nlvddi,n
	integer il(n)
	integer g_in(n_ghost_max,n_ghost_areas)
	integer g_out(n_ghost_max,n_ghost_areas)
	real val(n0:nlvddi,n)

	integer tag,ir,ia,id
	integer i,k,nc,ierr
	integer nb
	integer iout,iin
	integer nbs(2,n_ghost_areas)
	integer status(status_size,2*n_threads)
	integer request(2*n_threads)
	real, allocatable :: buffer_in(:,:)
	real, allocatable :: buffer_out(:,:)

        tag=122
	ir = 0

	iout = 2
	iin = 3
	if( belem ) then
	  iout = 4
	  iin = 5
	end if

	nb = (nlvddi-n0+1) * n_ghost_max
	call shympi_alloc_buffer(nb)
	allocate(buffer_in(nb,n_ghost_areas))
	allocate(buffer_out(nb,n_ghost_areas))

	do ia=1,n_ghost_areas
	  nc = ghost_areas(iout,ia)
	  call count_buffer(n0,nlvddi,n,nc,il,g_out(:,ia),nb)
	  nbs(1,ia) = nb
	  nc = ghost_areas(iin,ia)
	  call count_buffer(n0,nlvddi,n,nc,il,g_in(:,ia),nb)
	  nbs(2,ia) = nb
	end do

cccgguccc!$OMP CRITICAL

	do ia=1,n_ghost_areas
	  ir = ir + 1
	  id = ghost_areas(1,ia)
	  nc = ghost_areas(iout,ia)
	  nb = nbs(1,ia)
          call MPI_Irecv(buffer_out(:,ia),nb,MPI_REAL,id
     +	          ,tag,MPI_COMM_WORLD,request(ir),ierr)
	end do

	do ia=1,n_ghost_areas
	  ir = ir + 1
	  id = ghost_areas(1,ia)
	  nc = ghost_areas(iin,ia)
	  nb = nbs(2,ia)
	  call to_buffer_r(n0,nlvddi,n,nc,il
     +		,g_in(:,ia),val,nb,buffer_in(:,ia))
          call MPI_Isend(buffer_in(:,ia),nb,MPI_REAL,id
     +	          ,tag,MPI_COMM_WORLD,request(ir),ierr)
	end do

        call MPI_WaitAll(ir,request,status,ierr)

cccgguccc!$OMP END CRITICAL

	do ia=1,n_ghost_areas
	  id = ghost_areas(1,ia)
	  nc = ghost_areas(iout,ia)
	  nb = nbs(1,ia)
	  call from_buffer_r(n0,nlvddi,n,nc,il
     +		,g_out(:,ia),val,nb,buffer_out(:,ia))
	end do

	end subroutine shympi_exchange_internal_r

!******************************************************************

	subroutine shympi_exchange_internal_d(belem,n0,nlvddi,n,il
     +						,g_in,g_out,val)

	use shympi_aux
	use shympi
	
	logical belem
	integer n0,nlvddi,n
	integer il(n)
	integer g_in(n_ghost_max,n_ghost_areas)
	integer g_out(n_ghost_max,n_ghost_areas)
	double precision val(n0:nlvddi,n)

        stop 'error stop shympi_exchange_internal_d: not ready'

	end subroutine shympi_exchange_internal_d

!******************************************************************
!******************************************************************
!******************************************************************

        subroutine shympi_allgather_i_internal(n,no,val,vals)

	use shympi_aux
	use shympi

	implicit none

	integer n,no
        integer val(n)
        integer vals(no,n_threads)

        integer ierr

	if( bpmpi ) then
          call MPI_ALLGATHER (val,n,MPI_INTEGER
     +                  ,vals,no,MPI_INTEGER
     +                  ,MPI_COMM_WORLD,ierr)
	  call shympi_error('shympi_allgather_i_internal'
     +			,'gather',ierr)
	else
	  vals(:,1) = val(:)
	end if

        end subroutine shympi_allgather_i_internal

!******************************************************************

        subroutine shympi_allgather_r_internal(n,no,val,vals)

	use shympi_aux
	use shympi

	implicit none

	integer n,no
        real val(n)
        real vals(no,n_threads)

        integer ierr

	if( bpmpi ) then
          call MPI_ALLGATHER (val,n,MPI_REAL
     +                  ,vals,no,MPI_REAL
     +                  ,MPI_COMM_WORLD,ierr)
	  call shympi_error('shympi_allgather_i_internal'
     +			,'gather',ierr)
	else
	  vals(:,1) = val(:)
	end if

        end subroutine shympi_allgather_r_internal

!******************************************************************

        subroutine shympi_allgather_d_internal(n,no,val,vals)

	use shympi_aux
	use shympi

	implicit none

	integer n,no
        double precision val(n)
        double precision vals(no,n_threads)

        integer ierr

	if( bpmpi ) then
          call MPI_ALLGATHER (val,n,MPI_DOUBLE
     +                  ,vals,no,MPI_DOUBLE
     +                  ,MPI_COMM_WORLD,ierr)
	  call shympi_error('shympi_allgather_i_internal'
     +			,'gather',ierr)
	else
	  vals(:,1) = val(:)
	end if

        end subroutine shympi_allgather_d_internal

!******************************************************************

        subroutine shympi_bcast_i_internal(n,val)

	use shympi_aux
	use shympi

	implicit none

	integer n
        integer val(n)

        integer ierr

	if( bpmpi ) then
          call MPI_BCAST(val,n,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	  call shympi_error('shympi_bcast_i_internal','bcast',ierr)
	end if

        end subroutine shympi_bcast_i_internal

!******************************************************************

        subroutine shympi_bcast_r_internal(n,val)

	use shympi_aux
	use shympi

	implicit none

	integer n
        real val(n)

        integer ierr

	if( bpmpi ) then
          call MPI_BCAST(val,n,MPI_REAL,0,MPI_COMM_WORLD,ierr)
	  call shympi_error('shympi_bcast_r_internal','bcast',ierr)
	end if

        end subroutine shympi_bcast_r_internal

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_reduce_d_internal(what,val)

	use shympi_aux

	implicit none

	character*(*) what
	double precision val

        integer ierr
	double precision valout

	if( bpmpi ) then
         if( what == 'min' ) then
	  call MPI_ALLREDUCE(val,valout,1,MPI_DOUBLE,MPI_MIN
     +				,MPI_COMM_WORLD,ierr)
	  val = valout
         else if( what == 'max' ) then
	  call MPI_ALLREDUCE(val,valout,1,MPI_DOUBLE,MPI_MAX
     +				,MPI_COMM_WORLD,ierr)
	  val = valout
         else if( what == 'sum' ) then
	  call MPI_ALLREDUCE(val,valout,1,MPI_DOUBLE,MPI_SUM
     +				,MPI_COMM_WORLD,ierr)
	  val = valout
         else
          write(6,*) 'what = ',what
          stop 'error stop shympi_reduce_d_internal: not ready'
         end if
	else
	 ierr = 0
	end if

	call shympi_error('shympi_reduce_d_internal','reduce',ierr)

	end subroutine shympi_reduce_d_internal

!******************************************************************

	subroutine shympi_reduce_r_internal(what,val)

	use shympi_aux

	implicit none

	character*(*) what
	real val

        integer ierr
	real valout

	if( bpmpi ) then
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
	else
	 ierr = 0
	end if

	call shympi_error('shympi_reduce_r_internal','reduce',ierr)

	end subroutine shympi_reduce_r_internal

!******************************************************************

	subroutine shympi_reduce_i_internal(what,val)

	use shympi_aux

	implicit none

	character*(*) what
	integer val

        integer ierr
	integer valout

	if( bpmpi ) then
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
	else
	 ierr = 0
	end if

	call shympi_error('shympi_reduce_i_internal','reduce',ierr)

	end subroutine shympi_reduce_i_internal

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_getvals_internal_r(kind,nlvddi,n
     +						,val_in,val_out)

	use shympi_aux

	use shympi

	implicit none

	integer kind(2)
	integer nlvddi,n
	real val_in(nlvddi,n)
	real val_out(nlvddi)

	integer id,k,nb,lmax
	integer ierr

	id = kind(2) - 1
	k = kind(1)
	lmax = nlvddi
	nb = lmax

	if( my_id == id ) then
	  val_out(1:lmax) = val_in(1:lmax,k)
	end if

	!write(6,*) '========',id,k,nb,val_out(1)

        call MPI_BCAST(val_out,nb,MPI_REAL,id
     +	          ,MPI_COMM_WORLD,ierr)

	end subroutine shympi_getvals_internal_r

!******************************************************************

	subroutine shympi_getvals_internal_i(kind,nlvddi,n
     +						,val_in,val_out)

	use shympi_aux

	use shympi

	implicit none

	integer kind(2)
	integer nlvddi,n
	integer val_in(nlvddi,n)
	integer val_out(nlvddi)

	integer id,k,nb,lmax
	integer ierr

	id = kind(2) - 1
	k = kind(1)
	lmax = nlvddi
	nb = lmax

	if( my_id == id ) then
	  val_out(1:lmax) = val_in(1:lmax,k)
	end if

	!write(6,*) '========',id,k,nb,val_out(1)

        call MPI_BCAST(val_out,nb,MPI_INTEGER,id
     +	          ,MPI_COMM_WORLD,ierr)

	end subroutine shympi_getvals_internal_i

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_get_array_internal_r(nlvddi,n
     +						,val_in,val_out)

	use shympi_aux

	use shympi

	implicit none

	integer kind(2)
	integer nlvddi,n
	real val_in(nlvddi,*)
	real val_out(nlvddi,n)

	integer i,ir,ns,nb,tag,id
	integer ierr
	integer ip(0:n_threads)
	integer req(2*n_threads)
	integer status(status_size,2*n_threads)

        tag=123
	ir = 0

	if( n == nkn_global ) then
	  ip = nkn_cum_domains
	else if( n == nel_global ) then
	  ip = nel_cum_domains
	else
	  write(6,*) 'n,nkn_global,nel_global: ',n,nkn_global,nel_global
	  call shympi_stop('error stop shympi_get_array_internal_i:'//
     +				' size of out array')
	end if

cccgguccc!$OMP CRITICAL

	if( my_id == 0 ) then
	  do i=2,n_threads
	    id = i - 1
	    ir = ir + 1
	    ns = nlvddi*ip(i-1) + 1
	    nb = nlvddi*(ip(i) - ip(i-1))
            call MPI_Irecv(val_out(1,ns),nb,MPI_REAL,id
     +	          ,tag,MPI_COMM_WORLD,req(ir),ierr)
	  end do
	  nb = ip(1)
	  val_out(:,1:nb) = val_in(:,1:nb)
	else
	    i = my_id + 1
	    ir = ir + 1
	    nb = nlvddi*(ip(i) - ip(i-1))
            call MPI_Isend(val_in,nb,MPI_REAL,0
     +	          ,tag,MPI_COMM_WORLD,req(ir),ierr)
	end if

        call MPI_WaitAll(ir,req,status,ierr)

cccgguccc!$OMP END CRITICAL

	end subroutine shympi_get_array_internal_r

!******************************************************************

	subroutine shympi_get_array_internal_i(nlvddi,n
     +						,val_in,val_out)

	use shympi_aux

	use shympi

	implicit none

	integer nlvddi,n
	integer val_in(nlvddi,*)
	integer val_out(nlvddi,n)

	integer i,ir,ns,nb,tag,id
	integer ierr
	integer ip(0:n_threads)
	integer req(2*n_threads)
	integer status(status_size,2*n_threads)

        tag=124
	ir = 0

	if( n == nkn_global ) then
	  ip = nkn_cum_domains
	else if( n == nel_global ) then
	  ip = nel_cum_domains
	else
	  write(6,*) 'n,nkn_global,nel_global: ',n,nkn_global,nel_global
	  call shympi_stop('error stop shympi_get_array_internal_i:'//
     +				' size of out array')
	end if

cccgguccc!$OMP CRITICAL

	if( my_id == 0 ) then
	  do i=2,n_threads
	    id = i - 1
	    ir = ir + 1
	    ns = nlvddi*ip(i-1) + 1
	    nb = nlvddi*(ip(i) - ip(i-1))
            call MPI_Irecv(val_out(1,ns),nb,MPI_INTEGER,id
     +	          ,tag,MPI_COMM_WORLD,req(ir),ierr)
	  end do
	  nb = ip(1)
	  val_out(:,1:nb) = val_in(:,1:nb)
	else
	    i = my_id + 1
	    ir = ir + 1
	    nb = nlvddi*(ip(i) - ip(i-1))
            call MPI_Isend(val_in,nb,MPI_INTEGER,0
     +	          ,tag,MPI_COMM_WORLD,req(ir),ierr)
	end if

        call MPI_WaitAll(ir,req,status,ierr)

cccgguccc!$OMP END CRITICAL

	end subroutine shympi_get_array_internal_i

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_exchange_array_internal_r(nlin,nlout,nkin,nkout
     +						,val_in,val_out)

	use shympi_aux

	use shympi

	implicit none

	integer nlin,nlout,nkin,nkout
	real val_in(nlin,nkin)
	real val_out(nlout,nkout)

	integer i,ir,ns,ne,nb,tag,id,nn,n
	integer ierr
	integer ip(0:n_threads)
	integer req(2*n_threads)
	integer status(status_size,2*n_threads)
	real, allocatable :: valaux(:,:)

        tag=125
	ir = 0

	n = nkout
	if( n == nkn_global ) then
	  ip = nkn_cum_domains
	else if( n == nel_global ) then
	  ip = nel_cum_domains
	else
	  write(6,*) 'n,nkn_global,nel_global: ',n,nkn_global,nel_global
	  call shympi_stop('error stop shympi_exchange_array_internal_i:'
     +				//' size of outer array')
	end if

	!write(6,*) 'exchanging: ',nlvddi,n

cccgguccc!$OMP CRITICAL

	do i=1,n_threads
	  id = i - 1
	  if( id == my_id ) cycle
	  ir = ir + 1
	  ns = ip(i-1) + 1
	  ne = ip(i)
	  nn = ip(i) - ip(i-1)
	  nb = nlout*nn
	  !write(6,1000) 'receiving: ',my_id,id,ir,ns,ne,nb,nb/nlvddi
          call MPI_Irecv(val_out(1,ns),nb,MPI_REAL,id
     +	          ,tag,MPI_COMM_WORLD,req(ir),ierr)
	end do

	i = my_id + 1			!we always send from this id
	nn = ip(i) - ip(i-1)
	nb = nlout*nn
	!write(6,'(a,6i7)') 'internal: ',my_id,nn,nkin,nkout,nlin,nlout
	allocate(valaux(nlout,nn))
	valaux = 0.
	valaux(1:nlin,1:nn) = val_in(1:nlin,1:nn)

	do i=1,n_threads
	  id = i - 1
	  if( id == my_id ) cycle
	  ir = ir + 1
	  !write(6,1000) 'sending: ',my_id,id,ir,nb,nb/nlvddi
          call MPI_Isend(valaux,nb,MPI_REAL,id
     +	          ,tag,MPI_COMM_WORLD,req(ir),ierr)
	end do

	i = my_id + 1
	ns = ip(i-1) + 1
	ne = ip(i)
	nn = ip(i) - ip(i-1)
	!write(6,1000) 'copying: ',my_id,ns,ne,nn
	val_out(1:nlin,ns:ne) = val_in(1:nlin,1:nn)

        call MPI_WaitAll(ir,req,status,ierr)

cccgguccc!$OMP END CRITICAL

 1000	format(a,8i5)
	end subroutine shympi_exchange_array_internal_r

!******************************************************************

	subroutine shympi_exchange_array_internal_i(nlin,nlout,nkin,nkout
     +						,val_in,val_out)

	use shympi_aux

	use shympi

	implicit none

	integer nlin,nlout,nkin,nkout
	integer val_in(nlin,nkin)
	integer val_out(nlout,nkout)

	integer i,ir,ns,ne,nb,tag,id,nn,n
	integer ierr
	integer ip(0:n_threads)
	integer req(2*n_threads)
	integer status(status_size,2*n_threads)
	integer, allocatable :: valaux(:,:)

        tag=126
	ir = 0

	n = nkout
	if( n == nkn_global ) then
	  ip = nkn_cum_domains
	else if( n == nel_global ) then
	  ip = nel_cum_domains
	else
	  write(6,*) 'n,nkn_global,nel_global: ',n,nkn_global,nel_global
	  call shympi_stop('error stop shympi_exchange_array_internal_i:'
     +				//' size of outer array')
	end if

	!write(6,*) 'exchanging: ',nlvddi,n

cccgguccc!$OMP CRITICAL

	do i=1,n_threads
	  id = i - 1
	  if( id == my_id ) cycle
	  ir = ir + 1
	  ns = ip(i-1) + 1
	  ne = ip(i)
	  nn = ip(i) - ip(i-1)
	  nb = nlout*nn
	  !write(6,1000) 'receiving: ',my_id,id,ir,ns,ne,nb,nb/nlvddi
          call MPI_Irecv(val_out(1,ns),nb,MPI_INTEGER,id
     +	          ,tag,MPI_COMM_WORLD,req(ir),ierr)
	end do

	i = my_id + 1			!we always send from this id
	nn = ip(i) - ip(i-1)
	nb = nlout*nn
	allocate(valaux(nlout,nn))
	valaux = 0
	valaux(1:nlin,1:nn) = val_in(1:nlin,1:nn)

	do i=1,n_threads
	  id = i - 1
	  if( id == my_id ) cycle
	  ir = ir + 1
	  !write(6,1000) 'sending: ',my_id,id,ir,nb,nb/nlvddi
          call MPI_Isend(valaux,nb,MPI_INTEGER,id
     +	          ,tag,MPI_COMM_WORLD,req(ir),ierr)
	end do

	i = my_id + 1
	ns = ip(i-1) + 1
	ne = ip(i)
	nn = ip(i) - ip(i-1)
	!write(6,1000) 'copying: ',my_id,ns,ne,nn
	val_out(1:nlin,ns:ne) = val_in(1:nlin,1:nn)

        call MPI_WaitAll(ir,req,status,ierr)

cccgguccc!$OMP END CRITICAL

 1000	format(a,8i5)
	end subroutine shympi_exchange_array_internal_i

!******************************************************************
!******************************************************************
!******************************************************************


