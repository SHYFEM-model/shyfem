!
! interpolation routines from file
!
! revision log :
!
! 16.06.2014	ggu	time is now double
! 25.06.2014	ggu	various bug fixes
! 07.07.2014	ggu	first version finished
! 20.10.2014	ggu	deal with datetime in fem/ts files
! 07.01.2015	ggu	bug fix in iff_populate_records() -> handle holes
! 08.01.2015	ggu	bug fix for parallel: make variables local
! 05.02.2015	ggu	iff_read_and_interpolate() introduced for parallel bug
! 25.09.2015	ggu	prepared to interpolate from reg onto elements
! 29.09.2015	ggu	in iff_interpolate() do not interpolate with flag
! 18.12.2015	ggu	in iff_peek_next_record() adjust date only if ierr==0
! 15.02.2016	ggu	if dtime==-1 do not check time, changes in exffil
! 07.03.2016	ggu	bug fix in iff_read_header(): no time adjust if ierr/=0
! 01.04.2016	ggu	bug fix in iff_init(): open file with np and not nexp
! 05.05.2016	ggu	new functionality for 3d matrices
!
!****************************************************************
!
! info:
!
! iformat
!		-3  file finished (no more read or intp allowed)
!		-2  file closed (read EOF)
!		-1  no file given
!		 0  fem unformatted
!		 1  fem formatted
!		 2  fem unformatted direct
!		 3  ts formatted
!
! nintp		if 0 -> constant (do not read any more file)
!
! nvar		on entry: expected
!		on return: what has been read (should be checked)
! nexp		0 < nexp <= nkn
! lexp		vertical dimension, 0 if 2D field is expected
! np		number of points in file
! lmax		number of vertical points in file (lmax > 0)
!
! how to check for things:
!
! iform_ts = 3
! iform_none = -1
! iform_closed = -2
! iform_forget = -3
! iform_error = -9
!
! bnofile = iformat < 0				no file open
! bfem = iformat >= 0 .and. iformat <= 2	fem file format
! bts = iformat == iform_ts			ts file format
! bconst = nintp == 0				constant field
! b2d = lexp == 0				2D field expected
! bonepoint					only one point stored
!
! calling sequence:
!
! iff_init_global		intializes module
! iff_init(...,id)		initializes and gets file id
!	iff_populate_records
!		iff_read_next_record
!		iff_peek_next_record
!		iff_allocate_fem_data_structure
!		iff_space_interpolate
! iff_read_and_interpolate
!	iff_read_next_record
!	iff_space_interpolate
! iff_time_interpolate		interpolates for new time
!	iff_interpolate		interpolates in time
!
!****************************************************************

!================================================================
	module intp_fem_file
!================================================================

	implicit none

	type, private :: info
	  integer :: iunit = 0
	  integer :: id0 = 0		!take fdata from this id - not working
	  integer :: iformat = 0
	  integer :: nvers = 0
	  integer :: ntype = 0
	  integer :: nvar = 0
	  integer :: nintp = 0		!if 0 -> constant
	  integer :: ibc = 0		!may indicate number of open boundary
	  integer :: irec = 0
	  integer :: np = 0		!horizontal points in file
	  integer :: lmax = 0		!vertical layers in file
	  integer :: nexp = 0		!expected horizontal points
	  integer :: lexp = 0		!expected vertical points (0 for 2D)
	  integer :: ilast = 0		!last entry in time
	  integer :: ireg = 0		!regular grid
	  logical :: eof = .false.	!EOF encountered?

	  real :: flag = -999.		!flag for value not available
	  logical :: bneedall = .true.	!need values for all points

	  integer :: datetime(2) = 0	!date and time parameters
	  real :: regpar(7) = 0.	!parameters for regular grid

	  character*80 :: file = ' '
	  character*80 :: descript = ' '	!description of entry
	  logical :: bonepoint = .false.	!only one point stored
	  logical :: bfemdata = .false.		!fem data structure allocated
	  logical :: bfiledata = .false.	!file data structure allocated

	  integer, allocatable :: nodes(:)
	  double precision, allocatable :: time(:)
	  real, allocatable :: data(:,:,:,:)

	  double precision time_file		!maybe not needed
	  character*80, allocatable :: strings_file(:)
	  real, allocatable :: data_file(:,:,:)
	  real, allocatable :: hlv_file(:)
	  integer, allocatable :: ilhkv_file(:)
	  real, allocatable :: hd_file(:)
	end type info

	logical, parameter :: bassert = .true.

	integer, parameter :: iform_none = -1
	integer, parameter :: iform_closed = -2
	integer, parameter :: iform_forget = -3
	integer, parameter :: iform_error = -9
	integer, parameter :: iform_error_opening = -15
	integer, parameter :: iform_no_such_file = -11
	integer, parameter :: iform_ts = 3

	integer, parameter, private :: ndim = 300
	type(info), save, dimension(ndim) :: pinfo

	integer, save :: idlast = 0

	integer, save :: date_fem = 0
	integer, save :: time_fem = 0
	double precision, save :: atime0_fem = 0

	integer, save :: nkn_fem = 0
	integer, save :: nel_fem = 0
	integer, save :: nlv_fem = 0
	integer, save, allocatable :: ilhkv_fem(:)
	integer, save, allocatable :: ilhv_fem(:)
	real, save, allocatable :: hk_fem(:)
	real, save, allocatable :: he_fem(:)
	real, save, allocatable :: hlv_fem(:)

!================================================================
	contains
!================================================================

	subroutine iff_print_info(idp,iunit,bdebug)

	integer idp	!print info on this id, if 0 all info
	integer, optional :: iunit
	logical, optional :: bdebug

	integer id,ids,ide,iu
	logical debug
	integer ilast,ifirst
	character*38 name
	character*80 descrp

	iu = 6
	if( present(iunit) ) iu = iunit
	if( iu <= 0 ) iu = 6
	debug = .false.
	if( present(bdebug) ) debug = bdebug

	if( idp <= 0 ) then
	  ids = 1
	  ide = idlast
	else
	  ids = idp
	  ide = idp
	end if

	write(iu,*) 'iff_print_info:'
	write(iu,1010)
	do id=ids,ide
	  ilast = len_trim(pinfo(id)%file)
	  ifirst = max(1,ilast-38+1)
	  name = pinfo(id)%file(ifirst:ilast)
	  descrp = pinfo(id)%descript
	  write(iu,1000) id,pinfo(id)%ibc
     +			,pinfo(id)%iunit,pinfo(id)%nvar
     +			,pinfo(id)%nintp,pinfo(id)%iformat
     +			,descrp(1:10),name
	end do

	if( .not. debug ) return

	write(iu,*) 'debug information:'
	do id=ids,ide
	  write(iu,*) id,pinfo(id)%nvers,pinfo(id)%ntype,pinfo(id)%irec
	  write(iu,*) id,pinfo(id)%np,pinfo(id)%lmax,pinfo(id)%nexp
	  write(iu,*) id,pinfo(id)%ilast,pinfo(id)%bonepoint
	  if( pinfo(id)%bfemdata ) then
	  write(iu,*) id,'fem variables: nodes,time,data'
	  write(iu,*) id,pinfo(id)%nodes
	  write(iu,*) id,pinfo(id)%time
	  write(iu,*) id,pinfo(id)%data
	  end if
	  if( pinfo(id)%bfiledata ) then
	  write(iu,*) id,'file variables: hlv,time,data'
	  write(iu,*) id,pinfo(id)%hlv_file
	  write(iu,*) id,pinfo(id)%time_file
	  write(iu,*) id,pinfo(id)%data_file
	  end if
	  !write(iu,*) id,pinfo(id)%ilhkv_file
	  !write(iu,*) id,pinfo(id)%hd_file
	end do

	return
 1010	format('   id  ibc unit nvar intp form descript   file')
 1000	format(6i5,1x,a10,1x,a38)
	end subroutine iff_print_info

!****************************************************************

	subroutine iff_print_file_info(id)

	integer id

	call iff_print_info(id)

	end subroutine iff_print_file_info

!****************************************************************

	subroutine iff_print_boundary_info(ibc,iunit,bdebug)

	integer ibc
	integer iunit
	logical bdebug

	integer id

	do id=1,idlast
	  if( ibc .eq. pinfo(id)%ibc ) then
	    call iff_print_info(id,iunit,bdebug)
	  end if
	end do

	end subroutine iff_print_boundary_info

!****************************************************************

	subroutine iff_close_file(id)

	integer id

	pinfo(id)%iformat = iform_closed
	call iff_allocate_file_arrays(id,0,0,0)
	close(pinfo(id)%iunit)
	pinfo(id)%iunit = -2

	end subroutine iff_close_file

!****************************************************************

	subroutine iff_forget_file(id)

	integer id

	pinfo(id)%iformat = iform_forget
	call iff_delete_entry(id)
	close(pinfo(id)%iunit)
	pinfo(id)%iunit = -3

	end subroutine iff_forget_file

!****************************************************************

	subroutine iff_delete_entry(id)

	integer id

	deallocate(pinfo(id)%strings_file)

	if( pinfo(id)%bfemdata ) then
	  deallocate(pinfo(id)%time)
	  deallocate(pinfo(id)%data)
	end if

	if( pinfo(id)%bfiledata ) then
	  deallocate(pinfo(id)%hlv_file)
	  deallocate(pinfo(id)%data_file)
	  deallocate(pinfo(id)%ilhkv_file)
	  deallocate(pinfo(id)%hd_file)
	end if

	pinfo(id)%bfemdata = .false.
	pinfo(id)%bfiledata = .false.

	end subroutine iff_delete_entry

!****************************************************************

	subroutine iff_init_entry(id)

	integer id

	if( nkn_fem == 0 ) then
	  write(6,*) 'iff routines have not been initialized'
	  write(6,*) 'iff_init_global() must be called first'
	  stop 'error stop iff_init_entry: no initialization'
	end if

	pinfo(id)%iunit = 0

	end subroutine iff_init_entry

!****************************************************************
!****************************************************************
!****************************************************************

	subroutine iff_get_var_description(id,ivar,string)

	integer id
	integer ivar
	character*(*) string

	integer nvar

	nvar = pinfo(id)%nvar
	if( ivar < 1 .or. ivar > nvar ) goto 99

	string = pinfo(id)%strings_file(ivar)

	return
   99	continue
	write(6,*) 'error in parameter'
	write(6,*) 'ivar,nvar: ',ivar,nvar
	call iff_print_file_info(id)
	stop 'error stop iff_get_var_description'
	end subroutine iff_get_var_description

!****************************************************************

	subroutine iff_set_var_description(id,ivar,string)

	integer id
	integer ivar
	character*(*) string

	integer nvar

	nvar = pinfo(id)%nvar
	if( ivar < 1 .or. ivar > nvar ) goto 99

	pinfo(id)%strings_file(ivar) = string

	return
   99	continue
	write(6,*) 'error in parameter'
	write(6,*) 'ivar,nvar: ',ivar,nvar
	stop 'error stop iff_set_var_description'
	end subroutine iff_set_var_description

!****************************************************************

	subroutine iff_set_description(id,ibc,string)

	integer id
	integer ibc
	character*(*) string

	pinfo(id)%ibc = ibc
	pinfo(id)%descript = string

	end subroutine iff_set_description

!****************************************************************

	function iff_get_nvar(id)

	integer iff_get_nvar
	integer id

	iff_get_nvar = pinfo(id)%nvar

	end function iff_get_nvar

!****************************************************************

	function iff_is_onepoint(id)

	logical iff_is_onepoint
	integer id

	iff_is_onepoint = pinfo(id)%bonepoint

	end function iff_is_onepoint

!****************************************************************

	function iff_is_constant(id)

	logical iff_is_constant
	integer id

	iff_is_constant = pinfo(id)%nintp == 0

	end function iff_is_constant

!****************************************************************

	function iff_has_file(id)

	logical iff_has_file
	integer id

	iff_has_file = pinfo(id)%iformat /= iform_none

	end function iff_has_file

!****************************************************************
!****************************************************************
!****************************************************************

	subroutine iff_init_global(nkn,nel,nlv,ilhkv,ilhv
     +					,hkv,hev,hlv,date,time)

! passes params and arrays from fem needed for interpolation
!
! should be only called once at the beginning of the simulation

	integer nkn
	integer nel
	integer nlv
	integer ilhkv(nkn)
	integer ilhv(nel)
	real hkv(nkn)
	real hev(nel)
	real hlv(nlv)
	integer date,time

	integer i

	if( nkn_fem == 0 ) then
	  do i=1,ndim
	    pinfo(i)%iunit = 0
	  end do
	end if

	if( nkn /= nkn_fem .or. nel /= nel_fem 
     +				.or. nlv /= nlv_fem ) then
	  if( nkn_fem > 0 ) then
	    deallocate(ilhkv_fem,hk_fem)
	    deallocate(ilhv_fem,he_fem)
	    deallocate(hlv_fem)
	  end if
	  allocate(ilhkv_fem(nkn),hk_fem(nkn))
	  allocate(ilhv_fem(nel),he_fem(nel))
	  allocate(hlv_fem(nlv))
	end if

	nkn_fem = nkn
	nel_fem = nel
	nlv_fem = nlv
	ilhkv_fem = ilhkv
	ilhv_fem = ilhv
	hk_fem = hkv
	he_fem = hev
	hlv_fem = hlv

	call iff_init_global_date_internal(date,time)

	end subroutine iff_init_global

!****************************************************************

	subroutine iff_init_global_date_internal(date,time)

	integer date,time

	date_fem = date
	time_fem = time
	call dts_to_abs_time(date,time,atime0_fem)

	end subroutine iff_init_global_date_internal

!****************************************************************

	subroutine iff_init(dtime,file,nvar,nexp,lexp,nintp
     +					,nodes,vconst,id)

	double precision dtime	!initial time
	character*(*) file	!file name
	integer nvar		!expected number of variables (might change)
	integer nexp		!expected number of points
	integer lexp		!expected max vertical data (0 for 2D)
	integer nintp		!requested interpolation (2 linear, 4 cubic)
	integer nodes(nexp)	!internal nodes numbers if nexp < nkn
	real vconst(nvar)	!constant values in case no file is given
	integer id		!identification of file info (return)

! nvar is return value (if file can be read)

	integer iformat,iunit
	integer ierr,np,i
	integer nvar_orig
	integer datetime(2)
	integer ntype,itype(2)
	integer id0,ibc
	logical breg
	logical bok
	logical bts,bfem,bnofile,bfile,berror,bnosuchfile,boperr

	!---------------------------------------------------------
	! get new id for file
	!---------------------------------------------------------

	idlast = idlast + 1
	if( idlast > ndim ) then
	  stop 'error stop iff_init: too many files opened'
	end if
	id = idlast

	call iff_init_entry(id)

	!---------------------------------------------------------
	! check input parameters
	!---------------------------------------------------------

	if( nvar < 1 ) goto 97
	if( nexp < 1 ) goto 97
	!if( nexp > nkn_fem ) goto 97
	if( lexp < 0 ) goto 97
	if( nintp < 1 ) goto 97

	!---------------------------------------------------------
	! get info on file
	!---------------------------------------------------------

	nvar_orig = nvar

	call iff_get_file_info(file,np,nvar,ntype,iformat)

	bnofile = iformat == iform_none			!no file given
	bfile = .not. bnofile				!file has been given
	bnosuchfile = iformat == iform_no_such_file	!file not existing
	!bnofile = iformat < 0
	!berror = iformat < 0
	boperr = iformat == iform_error_opening		!error opening
	berror = iformat == iform_error			!error file

	id0 = 0
	if( boperr ) then	!see if we can take data from other file
	  id0 = iff_find_id_to_file(file)
	  id0 = 0		!not working
	  if( id0 > 0 ) then
	    boperr = .false.
	    iformat = pinfo(id0)%iformat
	    ntype = pinfo(id0)%ntype
	  end if
	end if

	bts = iformat == iform_ts			!file is TS
	bfem = iformat >= 0 .and. iformat <= 2

        call fem_file_make_type(ntype,2,itype)
        breg = itype(2) > 0

	if( file /= ' ' .and. bnofile ) goto 99
	if( bnosuchfile ) goto 99
	if( bfile .and. berror ) goto 93
	if( bfile .and. boperr ) goto 91
	if( bfile .and. np < 1 ) goto 96
	if( .not. breg .and. np > 1 .and. np /= nexp ) goto 96

	if( nvar <= 0 ) nvar = nvar_orig

	!---------------------------------------------------------
	! store information
	!---------------------------------------------------------

	pinfo(id)%iunit = -1
	pinfo(id)%nvar = nvar
	pinfo(id)%nintp = nintp
	pinfo(id)%file = file
	pinfo(id)%nexp = nexp
	pinfo(id)%lexp = lexp

	pinfo(id)%id0 = id0
	pinfo(id)%iformat = iformat
	pinfo(id)%ntype = ntype
	pinfo(id)%ireg = itype(2)

	!---------------------------------------------------------
	! get data description and allocate data structure
	!---------------------------------------------------------

	if( .not. breg .and. nexp > 0 
     +		.and. nexp /= nkn_fem .and. nexp /= nel_fem) then
	  allocate(pinfo(id)%nodes(nexp))	!lateral BC
	  pinfo(id)%nodes = nodes
	end if

	allocate(pinfo(id)%strings_file(nvar))

	if( id0 > 0 ) then
	  pinfo(id)%strings_file = pinfo(id0)%strings_file
	  pinfo(id)%bonepoint = pinfo(id0)%bonepoint
	else if( bfem ) then
          call fem_file_get_data_description(file
     +				,pinfo(id)%strings_file,ierr)
	  if( ierr /= 0 ) goto 98
	else if( bts ) then
	  pinfo(id)%bonepoint = .true.
	  pinfo(id)%strings_file = ' '
	else if( bnofile ) then		!constant
	  pinfo(id)%nintp = 0
	  pinfo(id)%ilast = 1
	  pinfo(id)%bonepoint = .true.
	  pinfo(id)%strings_file = ' '
	  call iff_allocate_fem_data_structure(id)
          pinfo(id)%time(1) = 0.
	  do i=1,nvar
            pinfo(id)%data(1,1,i,1) = vconst(i)
	  end do
	else
	  stop 'error stop iff_init: internal error (3)'
	end if

	!---------------------------------------------------------
	! finally open files
	!---------------------------------------------------------

	if( bnofile ) then
	  return
	else if( bfem ) then
	  !call fem_file_read_open(file,nexp,iunit,iformat)
	  call fem_file_read_open(file,np,iunit,iformat)
	else if( bts ) then
	  call ts_open_file(file,nvar,datetime,iunit)
	  pinfo(id)%datetime = datetime
	else
	  stop 'error stop iff_init: internal error (3)'
	end if

	if( iunit < 0 ) goto 99
	if( iunit == 0 ) goto 90
	pinfo(id)%iunit = iunit

	write(6,*) 'file opened: ',id,trim(file)

	!---------------------------------------------------------
	! populate data base
	!---------------------------------------------------------

	call iff_populate_records(id,dtime)

	!---------------------------------------------------------
	! end of routine
	!---------------------------------------------------------

	return
   90	continue
	write(6,*) 'error in opening file: ',trim(file)
	write(6,*) 'iformat = ',iformat
	stop 'error stop iff_init'
   91	continue
	write(6,*) 'error opening file: ',trim(file)
	id0 = iff_find_id_to_file(file)
	if( id0 > 0 ) then
	  ibc = pinfo(id0)%ibc
	  write(6,*) 'the file is already open on boundary ',ibc
	else 
	  write(6,*) 'iformat = ',iformat
	end if
	stop 'error stop iff_init'
   93	continue
	write(6,*) 'error in file: ',trim(file)
	write(6,*) 'iformat = ',iformat
	stop 'error stop iff_init'
   96	continue
	write(6,*) 'file does not contain expected data size'
	write(6,*) 'nexp,np: ',nexp,np
	call iff_print_file_info(id)
	stop 'error stop iff_init'
   97	continue
	write(6,*) 'error in input parameters of file: ',trim(file)
	write(6,*) 'nvar: ',nvar
	write(6,*) 'nexp,lexp: ',nexp,lexp
	write(6,*) 'nintp: ',nintp
	write(6,*) 'nkn_fem: ',nkn_fem
	call iff_print_file_info(id)
	stop 'error stop iff_init'
   98	continue
	write(6,*) 'error reading data description of file: ',trim(file)
	call iff_print_file_info(id)
	stop 'error stop iff_init'
   99	continue
	write(6,*) 'no such file: ',trim(file)
	write(6,*) 'iformat = ',iformat
	stop 'error stop iff_init'
	end subroutine iff_init

!****************************************************************

	subroutine iff_set_constant(id,vconst)

c (re-) sets constant if no file has been opened 

	integer id
	real vconst(pinfo(id)%nvar)

	integer iformat
	logical bnofile

	iformat = pinfo(id)%iformat
	bnofile = iformat < 0

	if( bnofile ) then
	  pinfo(id)%data_file(1,1,:) = vconst
	end if

	end subroutine iff_set_constant

!****************************************************************

	subroutine iff_get_file_info(file,np,nvar,ntype,iformat)

c coputes info on type of file
c
c	-1	no file
c	 0	unformatted
c	 1	formatted
c	 2	time series

	character*(*) file
	integer np
	integer nvar
	integer ntype
	integer iformat		!info on file type (return)

	integer il
	integer itype(2)
	logical filex

	np = 0
	nvar = 0
	ntype = 0
	iformat = iform_none

	if( file == ' ' ) return

	il = len_trim(file)

	if( .not. filex(file) ) then
	  iformat = iform_no_such_file
	  return
	end if

	call fem_file_test_formatted(file,np,nvar,ntype,iformat)

	if( nvar > 0 ) then
	  write(6,*) 'file is fem file with format: ',iformat
	  write(6,*) file(1:il)
	else
	  call ts_get_file_info(file,nvar)
	  if( nvar > 0 ) then
	    np = 1
	    ntype = 0
	    iformat = iform_ts
	    write(6,*) 'file is time series with columns: ',nvar
	    write(6,*) file(1:il)
	  else if( iformat == -77 ) then
	    !write(6,*) 'error opening file: ',file(1:il)
	    !write(6,*) '(maybe the file is already open?)'
	    iformat = iform_error_opening
	  else
	    write(6,*) 'cannot determine file format: ',file(1:il)
	    write(6,*) 'file is neither FEM file nor time series'
	    iformat = iform_error
	  end if
	end if

	if( ntype .gt. 0 ) then
	  call fem_file_make_type(ntype,2,itype)
	  if( itype(2) .gt. 0 ) then
	    write(6,*) 'file is regular file: ',itype(2)
	  end if
	end if

	end subroutine iff_get_file_info

!****************************************************************

	function iff_find_id_to_file(file)

! finds id given file, returns 0 if file is not open

	integer iff_find_id_to_file
	character*(*) file

	integer id

	do id=1,idlast
	  if( file == pinfo(id)%file ) exit
	end do
	if( id > idlast ) id = 0

	iff_find_id_to_file = id

	end function iff_find_id_to_file

!****************************************************************

	subroutine iff_populate_records(id,dtime0)

	integer id
	double precision dtime0

	!integer it,it2,idt,its,itold
	double precision dtime,dtime2
	!double precision dtimes,ddt
	double precision dtimefirst,dtimelast
	integer nintp,i
	logical bok,bts

        if( .not. iff_read_next_record(id,dtime) ) goto 99
	dtimefirst = dtime

        bok = iff_peek_next_record(id,dtime2)

	if( bok ) then				!at least two records
		nintp = pinfo(id)%nintp
		call iff_assert(nintp > 0,'nintp<=0')

                !ddt = dtime2 - dtime
                !if( ddt <= 0 ) goto 98
                !dtimes = dtime0 - nintp*ddt	!first record needed
		!if( dtime0 == -1 ) dtimes = dtime	!just take first

		!dtimeold = dtime
                !do while( dtime < dtimes )
                !        bok = iff_read_next_record(id,dtime)
                !        if( .not. bok ) goto 97
		!	ddt = dtime - dtimeold	!if time step changes
                !	if( ddt <= 0 ) goto 98
		!	dtimes = dtime0 - nintp*ddt
		!	dtimeold = dtime
                !end do

	        do
		  if( dtime0 == -1. ) exit	! no real time given
		  bok = iff_peek_next_record(id,dtime2)
                  if( .not. bok ) goto 97
		  if( dtime2 >= dtime0 ) exit
                  bok = iff_read_next_record(id,dtime)
                  if( .not. bok ) goto 97
		  dtimelast = dtime
		end do

		if( dtime0 /= -1. ) then
		  if( dtime0 < dtimefirst ) goto 91
		  if( dtime0 > dtime2 ) goto 91
		end if
		!write(6,*) 'populate: ',dtimefirst,dtime,dtime2,dtime0

		call iff_allocate_fem_data_structure(id)

                call iff_space_interpolate(id,1,dtime)
                do i=2,nintp
                        bok = iff_read_next_record(id,dtime)
                        if( .not. bok ) goto 96
                        call iff_space_interpolate(id,i,dtime)
                end do

		pinfo(id)%ilast = nintp
	else					!constant field
		pinfo(id)%nintp = 0
		pinfo(id)%ilast = 1
		call iff_allocate_fem_data_structure(id)
                call iff_space_interpolate(id,1,dtime)
		call iff_close_file(id)
        end if
		
	return
   91	continue
	call iff_print_file_info(id)
	write(6,*) 'cannot find time records'
	write(6,*) 'looking for it = ',dtime0
	write(6,*) 'first time found = ',dtimefirst
	stop 'error stop iff_populate_records: no time record found'
   96	continue
	call iff_print_file_info(id)
	write(6,*) 'cannot find enough time records'
	write(6,*) 'would need at least ',nintp
	stop 'error stop iff_populate_records: not enough records'
   97	continue
	call iff_print_file_info(id)
	write(6,*) 'cannot find time record'
	!write(6,*) 'looking at least for it = ',dtimes
	write(6,*) 'looking for it = ',dtime0
	write(6,*) 'last time found it = ',dtime
	stop 'error stop iff_populate_records: not enough records'
   98	continue
	call iff_print_file_info(id)
	write(6,*) 'time step less than 0'
	write(6,*) 'this happens at it = ',dtime
	stop 'error stop iff_populate_records: time step <= 0'
   99	continue
	call iff_print_file_info(id)
	write(6,*) 'error reading first record of file'
	stop 'error stop iff_populate_records: read error'
	end  subroutine iff_populate_records

!****************************************************************

	subroutine iff_allocate_fem_data_structure(id)

	integer id

	integer nexp,lexp,nvar,nintp
	logical bonepoint

	bonepoint = pinfo(id)%bonepoint
	nvar = pinfo(id)%nvar
	nexp = pinfo(id)%nexp
	lexp = max(1,pinfo(id)%lexp)
	nintp = max(1,pinfo(id)%nintp)

	if( bonepoint ) then	! time series or constant - store only one point
	  allocate(pinfo(id)%time(nintp))
	  allocate(pinfo(id)%data(1,1,nvar,nintp))
	else
	  allocate(pinfo(id)%time(nintp))
	  allocate(pinfo(id)%data(lexp,nexp,nvar,nintp))
	end if

	pinfo(id)%bfemdata = .true.

	pinfo(id)%time = 0.
	pinfo(id)%data = 0.

	end subroutine iff_allocate_fem_data_structure

!****************************************************************

        function iff_read_next_record(id,dtime)

	logical iff_read_next_record
	integer id
	double precision dtime

        if( iff_read_header(id,dtime) ) then
          call iff_read_data(id,dtime)
	  iff_read_next_record = .true.
	else
	  call iff_close_file(id)
	  iff_read_next_record = .false.
	end if

        end function iff_read_next_record

!****************************************************************

        function iff_peek_next_record(id,dtime)

! just gets new time stamp of next record

	logical iff_peek_next_record
	integer id
	double precision dtime

	integer iunit
	integer nvers,np,lmax,nvar,ntype
	integer datetime(2)
	integer ierr,iformat
	real f(1)
	logical bts,bnofile

	iunit = pinfo(id)%iunit
	iformat = pinfo(id)%iformat
	bnofile = iformat < 0
	bts = iformat == iform_ts

	iff_peek_next_record = .false.

	if( bnofile ) then
	  return
	else if( bts ) then
	  nvar = 0
	  call ts_peek_next_record(iunit,nvar,dtime,f,datetime,ierr)
	  if( ierr == 0 ) then
	    if( datetime(1) > 0 ) pinfo(id)%datetime = datetime
	    call iff_adjust_datetime(id,pinfo(id)%datetime,dtime)
	  end if
	else
          call fem_file_peek_params(iformat,iunit,dtime
     +                          ,nvers,np,lmax,nvar,ntype,datetime,ierr)
	  if( ierr == 0 ) then
	    pinfo(id)%datetime = datetime
	    call iff_adjust_datetime(id,pinfo(id)%datetime,dtime)
	  end if
	end if

	iff_peek_next_record = ( ierr == 0 )
	if( ierr .gt. 0 ) goto 99

	return
   99	continue
	write(6,*) 'read error in reading file header: ',ierr
	call iff_print_file_info(id)
	stop 'error stop iff_peek_next_record'
        end function iff_peek_next_record

!****************************************************************

        function iff_read_header(id,dtime)

	logical iff_read_header
	integer id
	double precision dtime

	logical bts,bnofile
	integer nvers,np,lmax,nvar,ntype
	integer datetime(2)
	integer iunit,ierr,iformat
	integer ivar
	real f(pinfo(id)%nvar)

	iff_read_header = .false.

	iunit = pinfo(id)%iunit
	ntype = pinfo(id)%ntype
	iformat = pinfo(id)%iformat
	bts = iformat == iform_ts
	bnofile = iformat < 0
	dtime = 0.
	datetime = 0

	if( bnofile ) return		!no header to read

	pinfo(id)%irec = pinfo(id)%irec + 1

	if( bts ) then
	  nvar = pinfo(id)%nvar
	  call ts_read_next_record(iunit,nvar,dtime,f,datetime,ierr)
	  if( ierr == 0 ) then
	    if( datetime(1) > 0 ) pinfo(id)%datetime = datetime
	    call iff_adjust_datetime(id,pinfo(id)%datetime,dtime)
	  end if
	else
          call fem_file_read_params(iformat,iunit,dtime
     +                          ,nvers,np,lmax,nvar,ntype,datetime,ierr)
	  if( ierr == 0 ) then
	    pinfo(id)%datetime = datetime
	    call iff_adjust_datetime(id,pinfo(id)%datetime,dtime)
	  end if
	end if

	if( ierr < 0 ) return
	if( ierr > 0 ) goto 99
	if( nvar /= pinfo(id)%nvar ) goto 98

	pinfo(id)%time_file = dtime

	if( bts ) then
	  call iff_allocate_file_arrays(id,nvar,1,1)
	  pinfo(id)%data_file(1,1,:) = f	!we have already read the data
	else
	  pinfo(id)%nvers = nvers
	  pinfo(id)%ntype = ntype
	  call iff_allocate_file_arrays(id,nvar,np,lmax)
	  call fem_file_read_2header(iformat,iunit,ntype,lmax
     +					,pinfo(id)%hlv_file
     +					,pinfo(id)%regpar
     +					,ierr)
	  if( ierr /= 0 ) goto 97
	end if

	iff_read_header = .true.

	return
   97	continue
	write(6,*) 'read error in reading hlv header'
	call iff_print_file_info(id)
	stop 'error stop iff_read_header'
   98	continue
	write(6,*) 'cannot change number of variables'
	call iff_print_file_info(id)
	write(6,*) 'nvar_old: ',pinfo(id)%nvar,' nvar_new: ',nvar
	stop 'error stop iff_read_header'
   99	continue
	write(6,*) 'read error in reading file header'
	call iff_print_file_info(id)
	stop 'error stop iff_read_header'
	end function iff_read_header

!****************************************************************

	subroutine iff_allocate_file_arrays(id,nvar,np,lmax)

! allocates file data structure
!
! if np/lmax are different from stored ones: first deallocate, then alloacte
! if np/lmax == 0: only deallocate

	integer id
	integer nvar
	integer np
	integer lmax

	integer lm

	!---------------------------------------------------------
	! check input params
	!---------------------------------------------------------

	if( np < 0 .or. lmax < 0 ) goto 99

	!---------------------------------------------------------
	! see if everything is the same
	!---------------------------------------------------------

	if( np == pinfo(id)%np .and. lmax == pinfo(id)%lmax ) return

	!---------------------------------------------------------
	! check consistency
	!---------------------------------------------------------

	if( pinfo(id)%np > 0 .and. pinfo(id)%lmax == 0 ) goto 98
	if( pinfo(id)%np == 0 .and. pinfo(id)%lmax > 0 ) goto 98
	if( np > 0 .and. lmax == 0 ) goto 98
	if( np == 0 .and. lmax > 0 ) goto 98

	!---------------------------------------------------------
	! deallocate old arrays
	!---------------------------------------------------------

	if( pinfo(id)%np > 0 .or. pinfo(id)%lmax > 0 ) then
	  deallocate(pinfo(id)%hlv_file)
	  deallocate(pinfo(id)%data_file)
	  deallocate(pinfo(id)%ilhkv_file)
	  deallocate(pinfo(id)%hd_file)
	end if

	!---------------------------------------------------------
	! allocate new arrays
	!---------------------------------------------------------

	pinfo(id)%bfiledata = .false.

	pinfo(id)%np = np
	pinfo(id)%lmax = lmax

	if( np == 0 .or. lmax == 0 ) return

	if( np > 0 .or. lmax > 0 ) then
	  allocate(pinfo(id)%hlv_file(lmax))
	  allocate(pinfo(id)%data_file(lmax,np,nvar))
	  allocate(pinfo(id)%ilhkv_file(np))
	  allocate(pinfo(id)%hd_file(np))
	end if

	pinfo(id)%bfiledata = .true.

	pinfo(id)%hlv_file = 0.
	pinfo(id)%data_file = 0.
	pinfo(id)%ilhkv_file = 0
	pinfo(id)%hd_file = 0.

	!---------------------------------------------------------
	! end of routine
	!---------------------------------------------------------

	return
   98	continue
	write(6,*) 'error in parameters: '
	write(6,*) 'stored parameters: '
	write(6,*) 'np = ',pinfo(id)%np,' lmax = ',pinfo(id)%lmax
	write(6,*) 'requested parameters: '
	write(6,*) 'np = ',np,' lmax = ',lmax
	stop 'error stop iff_allocate_file_arrays: internal error (1)'
   99	continue
	write(6,*) 'error in parameters: '
	write(6,*) 'np = ',np,' lmax = ',lmax
	stop 'error stop iff_allocate_file_arrays'
	end subroutine iff_allocate_file_arrays

!****************************************************************

        subroutine iff_read_data(id,dtime)

	integer id
	double precision dtime

	integer iunit,nvers,np,lmax
	integer nlvddi,nvar
	integer ierr,i,iformat
	logical bnofile,bts
	character*60 string

	iformat = pinfo(id)%iformat
	bts = iformat == iform_ts
	bnofile = iformat < 0

	if( bnofile ) return

	iunit = pinfo(id)%iunit
	nvers = pinfo(id)%nvers
	np = pinfo(id)%np
	lmax = pinfo(id)%lmax
	nvar = pinfo(id)%nvar

	nlvddi = lmax

	if( bts ) then
	  ! ts data has already been read
	else
	  do i=1,nvar
            call fem_file_read_data(iformat,iunit
     +                          ,nvers,np,lmax
     +                          ,string
     +                          ,pinfo(id)%ilhkv_file
     +                          ,pinfo(id)%hd_file
     +                          ,nlvddi
     +				,pinfo(id)%data_file(1,1,i)
     +				,ierr)
	    if( ierr /= 0 ) goto 99
	    if( string /= pinfo(id)%strings_file(i) ) goto 98
	  end do
	end if

	return
   98	continue
	write(6,*) 'string description has changed for var ',i
	write(6,*) 'time: ',dtime
	write(6,*) 'old: ',pinfo(id)%strings_file(i)
	write(6,*) 'new: ',string
	call iff_print_file_info(id)
	stop 'error stop iff_read_data'
   99	continue
	write(6,*) 'error reading data: ',ierr
	write(6,*) 'time: ',dtime
	call iff_print_file_info(id)
	stop 'error stop iff_read_data'
	end subroutine iff_read_data

!****************************************************************
!****************************************************************
!****************************************************************

	subroutine iff_space_interpolate(id,iintp,dtime)

c interpolates in space all variables in data set id

	integer id
	integer iintp
	double precision dtime

	integer nintp,np,nexp,lexp,ip
	integer ivar,nvar,ireg
	integer l,j,lfem,ipl
	logical bts,bdebug

	bdebug = .false.
	!bdebug = id == 6

        pinfo(id)%time(iintp) = dtime

        nintp = pinfo(id)%nintp
        nvar = pinfo(id)%nvar
        np = pinfo(id)%np
        nexp = pinfo(id)%nexp
        lexp = pinfo(id)%lexp
        ireg = pinfo(id)%ireg
	bts = pinfo(id)%iformat == iform_ts

	if( nintp > 0 .and. iintp > nintp ) goto 99
	if( nintp == 0 .and. iintp > 1 ) goto 99

	if( ireg > 0 ) then
	  if( lexp > 1 ) then
	    call iff_handle_regular_grid_3d(id,iintp)
	  else
	    call iff_handle_regular_grid_2d(id,iintp)
	  end if
	else if( bts ) then
          nvar = pinfo(id)%nvar
	  do ivar=1,nvar
	    pinfo(id)%data(1,1,ivar,iintp) = pinfo(id)%data_file(1,1,ivar)
	  end do
	else if( np == 1 ) then
	  do ip=1,nexp
	    call iff_handle_vertical(id,iintp,1,ip)
	  end do
	else if( np /= nexp ) then
	  goto 98
	else
	  do ip=1,np
	    call iff_handle_vertical(id,iintp,ip,ip)
	  end do
	end if

	if( bdebug ) then
	  write(6,*) 'iff_space_interpolate: data ---------------'
	  write(6,*) np,dtime
	  do j=1,nintp
	    write(6,*) 'iintp = ',j
	    do ivar=1,nvar
	      write(6,*) 'ivar = ',ivar
	      do ip=1,nexp
		ipl = ip
		if( nexp /= nkn_fem ) ipl = pinfo(id)%nodes(ip)
		lfem = ilhkv_fem(ipl)
	        write(6,*) 'node = ',ip,lfem,lexp
	        write(6,*) (pinfo(id)%data(l,ip,ivar,j),l=1,lexp)
	      end do
	    end do
	  end do
	  write(6,*) 'end iff_space_interpolate: data -----------'
	end if

	return
   96	continue
	write(6,*) 'regular grid only for 2d field'
	write(6,*) 'ireg,lexp: ',ireg,lexp
	!call iff_print_file_info(id)
	stop 'error stop iff_space_interpolate'
   97	continue
	write(6,*) 'cannot yet handle ntype > 0'
	!call iff_print_file_info(id)
	stop 'error stop iff_space_interpolate'
   98	continue
	write(6,*) 'error in number of points: ',np,nexp
	!call iff_print_file_info(id)
	stop 'error stop iff_space_interpolate'
   99	continue
	write(6,*) 'error in parameters: ',iintp,nintp
	!call iff_print_file_info(id)
	stop 'error stop iff_space_interpolate'
	end subroutine iff_space_interpolate

!****************************************************************

	subroutine iff_handle_regular_grid_2d(id,iintp)

	integer id
	integer iintp

	logical bneedall
	integer ivar,nvar
	integer nx,ny
	integer nexp
	integer ierr
	real x0,y0,dx,dy,flag

        nvar = pinfo(id)%nvar

	nx = nint(pinfo(id)%regpar(1))
	ny = nint(pinfo(id)%regpar(2))
	x0 = pinfo(id)%regpar(3)
	y0 = pinfo(id)%regpar(4)
	dx = pinfo(id)%regpar(5)
	dy = pinfo(id)%regpar(6)
	flag = pinfo(id)%regpar(7)
        nexp = pinfo(id)%nexp
	bneedall = pinfo(id)%bneedall
	pinfo(id)%flag = flag		!use this in time_interpolate

	if( nexp == nkn_fem ) then
	  do ivar=1,nvar
	    call intp_reg_nodes(nx,ny,x0,y0,dx,dy,flag
     +			,pinfo(id)%data_file(1,1,ivar)
     +			,pinfo(id)%data(1,1,ivar,iintp)
     +			,ierr
     +		    )
	    !if( bneedall .and. ierr .ne. 0 ) goto 99	!is handled later
	  end do
	else if( nexp == nel_fem ) then
	  do ivar=1,nvar
	    call intp_reg_elems(nx,ny,x0,y0,dx,dy,flag
     +			,pinfo(id)%data_file(1,1,ivar)
     +			,pinfo(id)%data(1,1,ivar,iintp)
     +			,ierr
     +		    )
	    !if( bneedall .and. ierr .ne. 0 ) goto 99	!is handled later
	  end do
	else
	  write(6,*) 'nexp,nkn,nel: ',nexp,nkn_fem,nel_fem
	  write(6,*) 'Cannot yet handle...'
	  stop 'error stop iff_handle_regular_grid_2d: nexp'
	end if

	return
   99	continue
	write(6,*) 'error interpolating from regular grid: '
	write(6,*) 'ierr =  ',ierr
	write(6,*) 'bneedall =  ',bneedall
	stop 'error stop iff_handle_regular_grid_2d: reg interpolate'
	end subroutine iff_handle_regular_grid_2d

!****************************************************************

	subroutine iff_handle_regular_grid_3d(id,iintp)

	integer id
	integer iintp

	logical bneedall
	integer ivar,nvar
	integer nx,ny
	integer nexp,lexp
	integer np,ip,l,lmax
	integer ierr
	real x0,y0,dx,dy,flag
	real, allocatable :: fr(:,:)
	real, allocatable :: data(:,:,:)
	real, allocatable :: data2dreg(:)
	real, allocatable :: data2dfem(:)
	real, allocatable :: hfem(:)

        nvar = pinfo(id)%nvar

	nx = nint(pinfo(id)%regpar(1))
	ny = nint(pinfo(id)%regpar(2))
	x0 = pinfo(id)%regpar(3)
	y0 = pinfo(id)%regpar(4)
	dx = pinfo(id)%regpar(5)
	dy = pinfo(id)%regpar(6)
	flag = pinfo(id)%regpar(7)
        nexp = pinfo(id)%nexp
        np = pinfo(id)%np
        lexp = pinfo(id)%lexp
        lmax = pinfo(id)%lmax		!vertical size of data
	bneedall = pinfo(id)%bneedall
	pinfo(id)%flag = flag		!use this in time_interpolate

	if( np /= nx*ny ) goto 95
	if( nexp /= nkn_fem .and. nexp /= nel_fem ) goto 98

	allocate(fr(4,nexp))
	allocate(data(lmax,nexp,nvar))
	allocate(data2dreg(np),data2dfem(nexp))
	allocate(hfem(nexp))

	call intp_reg_setup_fr(nx,ny,x0,y0,dx,dy,nexp,fr)
	call intp_reg_intp_fr(nx,ny,flag,pinfo(id)%hd_file
     +                          ,nexp,fr,hfem,ierr)

	do ivar=1,nvar
	  do l=1,lmax
	    data2dreg = pinfo(id)%data_file(l,:,ivar)
	    call intp_reg_intp_fr(nx,ny,flag,data2dreg
     +                          ,nexp,fr,data2dfem,ierr)
	    data(l,:,ivar) = data2dfem
	  end do
	end do

	do ip=1,nexp
	  call iff_interpolate_vertical_int(id,iintp
     +				,lmax,hfem(ip),data(:,ip,:),ip)
	end do

	return
   95	continue
	write(6,*) 'np,nx*ny: ',np,nx*ny
	stop 'error stop iff_handle_regular_grid_3d: internal error (1)'
   96	continue
	write(6,*) 'regular grid only for 2d field'
	write(6,*) 'lexp: ',lexp
	!call iff_print_file_info(id)
	stop 'error stop iff_handle_regular_grid_3d: not ready'
   98	continue
	write(6,*) 'nexp,nkn,nel: ',nexp,nkn_fem,nel_fem
	write(6,*) 'Cannot yet handle...'
	stop 'error stop iff_handle_regular_grid_3d: nexp'
   99	continue
	write(6,*) 'error interpolating from regular grid: '
	write(6,*) 'ierr =  ',ierr
	write(6,*) 'bneedall =  ',bneedall
	stop 'error stop iff_handle_regular_grid_3d: reg interpolate'
	end subroutine iff_handle_regular_grid_3d

!****************************************************************

	subroutine iff_handle_vertical(id,iintp,ip_from,ip_to)

	integer id
	integer iintp
	integer ip_from,ip_to

	integer lmax,lexp

        lmax = pinfo(id)%lmax
        lexp = pinfo(id)%lexp

	if( lmax == 1 ) then		!2D field given
	  call iff_distribute_vertical(id,iintp,ip_from,ip_to)
	else
	  if( lexp < 1 ) goto 99	!3D field given but 2D needed
	  if( lexp == 1 ) then
	    call iff_integrate_vertical(id,iintp,ip_from,ip_to)
	  else
	    call iff_interpolate_vertical(id,iintp,ip_from,ip_to)
	  end if
	end if

	return
   99	continue
	write(6,*) 'applying 3D data to 2D field'
	call iff_print_file_info(id)
	stop 'error stop iff_handle_vertical'
	end subroutine iff_handle_vertical

!****************************************************************

	subroutine iff_distribute_vertical(id,iintp,ip_from,ip_to)

! lmax must be 1

	integer id
	integer iintp
	integer ip_from,ip_to

	real data(pinfo(id)%nvar)

	data = pinfo(id)%data_file(1,ip_from,:)
	call iff_distribute_vertical_int(id,iintp,data,ip_to)

	end subroutine iff_distribute_vertical

!****************************************************************

	subroutine iff_distribute_vertical_int(id,iintp,data,ip_to)

! lmax must be 1

	integer id
	integer iintp
	real data(pinfo(id)%nvar)
	integer ip_to

	integer lfem,l,ipl,lexp
	integer ivar,nvar
	real value

        nvar = pinfo(id)%nvar
	lexp = pinfo(id)%lexp

	if( lexp <= 1 ) then		!2D
	  lfem = 1
	else				!3D
	  ipl = ip_to
	  if( pinfo(id)%nexp /= nkn_fem ) ipl = pinfo(id)%nodes(ip_to)
	  lfem = ilhkv_fem(ipl)
	end if

	do ivar=1,nvar
	  value = data(ivar)
	  do l=1,lfem
	    pinfo(id)%data(l,ip_to,ivar,iintp) = value
	  end do
	end do

	end subroutine iff_distribute_vertical_int

!****************************************************************

	subroutine iff_integrate_vertical(id,iintp,ip_from,ip_to)

! lexp/lfem must be 1

	integer id
	integer iintp
	integer ip_from,ip_to

	integer lmax
	real h
	real data(pinfo(id)%lmax,pinfo(id)%nvar)

	lmax = pinfo(id)%ilhkv_file(ip_from)
	h = pinfo(id)%hd_file(ip_from)
	data = pinfo(id)%data_file(:,ip_from,:)

	call iff_integrate_vertical_int(id,iintp,lmax,h,data,ip_to)

	end subroutine iff_integrate_vertical

!****************************************************************

	subroutine iff_integrate_vertical_int(id,iintp,lmax,h,data,ip_to)

! lexp/lfem must be 1

	integer id
	integer iintp
	integer lmax
	real h
	real data(pinfo(id)%lmax,pinfo(id)%nvar)
	integer ip_to

	integer l
	integer ivar,nvar
	integer nsigma
	double precision acum,htot
	real z
	real hsigma
	real value,hlayer
	real hl(pinfo(id)%lmax)

        nvar = pinfo(id)%nvar
	if( h < -990 ) h = pinfo(id)%hlv_file(lmax)	!take from hlv array
	z = 0.

	call compute_sigma_info(lmax,pinfo(id)%hlv_file,nsigma,hsigma)
	call get_layer_thickness(lmax,nsigma,hsigma,z,h
     +					,pinfo(id)%hlv_file,hl)

	do ivar=1,nvar
	  acum = 0.
	  htot = 0.
	  do l=1,lmax
	    value = data(l,ivar)
	    hlayer = hl(l)
	    acum = acum + value*hlayer
	    htot = htot + hlayer
	  end do
	  pinfo(id)%data(1,ip_to,ivar,iintp) = acum / htot
	end do

	end subroutine iff_integrate_vertical_int

!****************************************************************

	subroutine iff_interpolate_vertical(id,iintp,ip_from,ip_to)

c global lmax and lexp are > 1

	integer id
	integer iintp
	integer ip_from,ip_to

	integer lmax,lfem,ipl
	real h
	real data(pinfo(id)%lmax,pinfo(id)%nvar)

	lmax = pinfo(id)%ilhkv_file(ip_from)
	h = pinfo(id)%hd_file(ip_from)
	data = pinfo(id)%data_file(:,ip_from,:)

	ipl = ip_to
	if( pinfo(id)%nexp /= nkn_fem ) ipl = pinfo(id)%nodes(ip_to)
	lfem = ilhkv_fem(ipl)

	if( lmax <= 1 ) then
	  call iff_distribute_vertical_int(id,iintp,data(1,:),ip_to)
	else if( lfem <= 1 ) then
	  call iff_integrate_vertical_int(id,iintp,lmax,h,data,ip_to)
	else
	  call iff_interpolate_vertical_int(id,iintp,lmax,h,data,ip_to)
	end if

	end subroutine iff_interpolate_vertical

!****************************************************************

	subroutine iff_interpolate_vertical_int(id,iintp,lmax,h,data
     +						,ip_to)

c global lmax and lexp are > 1

	integer id
	integer iintp
	integer lmax
	real h
	real data(pinfo(id)%lmax,pinfo(id)%nvar)
	integer ip_to

	logical bcenter,bcons
	integer l,ipl,lfem
	integer ivar,nvar
	integer nsigma
	double precision acum,htot
	real z,hfem,hfile
	real hsigma
	real value,hlayer
	real hl(pinfo(id)%lmax)
	real hz_file(0:pinfo(id)%lmax+1)
	real val_file(pinfo(id)%lmax+1)
	real hl_fem(nlv_fem)
	real hz_fem(0:nlv_fem)
	real val_fem(nlv_fem)
	logical bdebug

	bcenter = .false.
	bcons = .false.
	bdebug = .true.
	bdebug = .false.

        nvar = pinfo(id)%nvar

	ipl = ip_to
	if( pinfo(id)%nexp /= nkn_fem ) ipl = pinfo(id)%nodes(ip_to)
	lfem = ilhkv_fem(ipl)

	if( lmax <= 1 ) then
	  call iff_distribute_vertical_int(id,iintp,data(1,:),ip_to)
	  return
	else if( lfem <= 1 ) then
	  call iff_integrate_vertical_int(id,iintp,lmax,h,data,ip_to)
	  return
	end if

	z = 0.

	hfem = hk_fem(ipl)
	hfile = h
	if( hfile < -990 ) hfile = pinfo(id)%hlv_file(lmax) !take from hlv array

	call compute_sigma_info(lmax,pinfo(id)%hlv_file,nsigma,hsigma)
	call get_layer_thickness(lmax,nsigma,hsigma,z,hfile
     +					,pinfo(id)%hlv_file,hl)
	call get_bottom_of_layer(bcenter,lmax,z,hl,hz_file(1))
	hz_file(0) = z

	call compute_sigma_info(lfem,hlv_fem,nsigma,hsigma)
	call get_layer_thickness(lfem,nsigma,hsigma,z,hfem
     +					,hlv_fem,hl_fem)
	call get_bottom_of_layer(bcenter,lfem,z,hl_fem,hz_fem(1))
	hz_fem(0) = z

	if( bdebug ) then
	  write(6,*) 'iff_interpolate_vertical: -------------------'
	  write(6,*) id
	  write(6,*) ip_to
	  write(6,*) hfile,hfem
	  write(6,*) lmax,lfem
	  write(6,*) 'hlv_file: ',(pinfo(id)%hlv_file(l),l=1,lmax)
	  write(6,*) 'hlv_fem: ',(hlv_fem(l),l=1,lfem)
	  write(6,*) (hl(l),l=1,lmax)
	  write(6,*) (hl_fem(l),l=1,lfem)
	  write(6,*) (hz_file(l),l=0,lmax)
	  write(6,*) (hz_fem(l),l=0,lfem)
	  write(6,*) 'end iff_interpolate_vertical: ----------------'
	end if

	do ivar=1,nvar
	  val_file(1:lmax) = data(1:lmax,ivar)
	  call intp_vert(bcons,lmax,hz_file,val_file,lfem,hz_fem,val_fem)
	  pinfo(id)%data(1:lfem,ip_to,ivar,iintp) = val_fem(1:lfem)
	  if( bdebug ) then
	    write(6,*) 'iff_interpolate_vertical - ivar = : ',ivar
	    write(6,*) (val_file(l),l=1,lmax)
	    write(6,*) (val_fem(l),l=1,lfem)
	    write(6,*) 'end iff_interpolate_vertical - ivar = : ',ivar
	  end if
	end do

	end subroutine iff_interpolate_vertical_int

!****************************************************************
!****************************************************************
!****************************************************************

	subroutine iff_time_interpolate(id,itact,ivar,ndim,ldim,value)

	integer id
	double precision itact
	integer ivar
	integer ndim		!horizontal dimension of value
	integer ldim		!vertical dimension of value
	real value(ldim,ndim)

	integer iv,nvar,iformat
	integer nintp,lexp,nexp
	integer ilast,ifirst
	double precision it,itlast,itfirst
	logical bok
	double precision t,tc

	!---------------------------------------------------------
	! set up parameters
	!---------------------------------------------------------

	if( id < 1 .or. id > idlast ) goto 95
        iformat = pinfo(id)%iformat
	if( iformat == -3 ) goto 96

        nintp = pinfo(id)%nintp
        nvar = pinfo(id)%nvar
        lexp = max(1,pinfo(id)%lexp)
        nexp = pinfo(id)%nexp
	t = itact

	!---------------------------------------------------------
	! loop until time window is centered over desidered time
	!---------------------------------------------------------

	if( iff_must_read(id,t) ) then
	  write(6,*) 'warning: reading data in iff_time_interpolate'
	  stop 'error stop iff_time_interpolate: internal error (1)'
	  !call iff_read_and_interpolate(id,t)
	end if

	!---------------------------------------------------------
	! some sanity checks
	!---------------------------------------------------------

	if( nintp > 0 ) then
          ilast = pinfo(id)%ilast
	  itlast = nint(pinfo(id)%time(ilast))
	  ifirst = mod(ilast,nintp) + 1
	  itfirst = nint(pinfo(id)%time(ifirst))
          if( itlast < itact ) goto 98
          if( itfirst > itact ) goto 98
	end if
	if( lexp > ldim ) goto 97
	if( nexp > ndim ) goto 97

	!---------------------------------------------------------
	! do the interpolation
	!---------------------------------------------------------

	if( ivar .eq. 0 ) then
	  do iv=1,nvar
            call iff_interpolate(id,t,iv,ndim,ldim,value)
	  end do
	else
          call iff_interpolate(id,t,ivar,ndim,ldim,value)
	end if

	!---------------------------------------------------------
	! end of routine
	!---------------------------------------------------------

	return
   95	continue
	write(6,*) 'id out of range: ',id,idlast
	call iff_print_file_info(0)
	stop 'error stop iff_time_interpolate: internal error (1)'
   96	continue
	write(6,*) 'file is closed... cannot interpolate anymore'
	call iff_print_file_info(id)
	stop 'error stop iff_time_interpolate: internal error (2)'
   97	continue
	write(6,*) 'incompatible dimensions'
	write(6,*) 'ldim,lexp: ',ldim,lexp
	write(6,*) 'ndim,nexp: ',ndim,nexp
	call iff_print_file_info(id)
	stop 'error stop iff_time_interpolate'
   98	continue
	write(6,*) 'file does not contain needed time value'
	write(6,*) 'looking for time: ',itact
	write(6,*) 'first time available: ',itfirst
	write(6,*) 'last time available: ',itlast
	call iff_print_file_info(id)
	stop 'error stop iff_time_interpolatei: time out of range'
	end subroutine iff_time_interpolate

!****************************************************************

	function iff_must_read(id,t)

c this routine determines if new data has to be read from file

	logical iff_must_read
	integer id
	double precision t		!time for which to interpolate

	integer ilast,nintp
	double precision tc

        nintp = pinfo(id)%nintp
        ilast = pinfo(id)%ilast			!index of last record

        iff_must_read = .false.
	if( pinfo(id)%eof ) return		!already at EOF

	tc = tcomp(t,nintp,ilast,pinfo(id)%time)

        iff_must_read = ( tc < t )

	end function iff_must_read

!****************************************************************

	subroutine iff_read_and_interpolate(id,t)

c this routine reads and interpolates new data - no parallel execution

	integer id
	double precision t		!time for which to interpolate

	logical bok
	integer ilast,nintp
	double precision itlast,it
	double precision tc		!check time

	bok = .true.
        nintp = pinfo(id)%nintp
        ilast = pinfo(id)%ilast			!index of last record
	itlast = nint(pinfo(id)%time(ilast))	!time of last record

	tc = tcomp(t,nintp,ilast,pinfo(id)%time)

        do while( tc < t )
          bok = iff_read_next_record(id,it)
	  if( .not. bok ) exit
	  if( it <= itlast ) goto 99
	  ilast = mod(ilast,nintp) + 1
	  call iff_space_interpolate(id,ilast,it)
	  tc = tcomp(t,nintp,ilast,pinfo(id)%time)
	  itlast = it
        end do

        pinfo(id)%eof = .not. bok
        pinfo(id)%ilast = ilast

	return
   99	continue
	write(6,*) 'time record not in increasing sequence'
	write(6,*) 'it,itlast: ',it,itlast
	call iff_print_file_info(id)
	stop 'error stop iff_read_and_interpolate'
	end subroutine iff_read_and_interpolate

!****************************************************************

        subroutine iff_interpolate(id,t,ivar,ndim,ldim,value)

c does the final interpolation in time

	integer id
	double precision t
	integer ivar
	integer ndim		!horizontal dimension of value
	integer ldim		!vertical dimension of value
	real value(ldim,ndim)

	integer nintp,lexp,nexp,ilast
	logical bonepoint,bconst,bnodes,b2d,bmulti,bflag
	integer ipl,lfem,i,l,ip,j,iflag
	real val,tr,flag
	double precision time(pinfo(id)%nintp)
	!real time(pinfo(id)%nintp)
	real vals(pinfo(id)%nintp)
	double precision rd_intp_neville
	real intp_neville

        nintp = pinfo(id)%nintp
        lexp = max(1,pinfo(id)%lexp)
        nexp = pinfo(id)%nexp
        ilast = pinfo(id)%ilast
        flag = pinfo(id)%flag
        bonepoint = pinfo(id)%bonepoint
	bconst = nintp == 0
	bmulti = .not. bonepoint
	b2d = lexp <= 1
	bnodes = pinfo(id)%nexp /= nkn_fem	!use node pointer

	tr = t		!real version of t
	iflag = 0

	if( bconst .or. bonepoint ) then
	  if( bconst ) then
	    val = pinfo(id)%data(1,1,ivar,1)
	  else
	    time = pinfo(id)%time
	    bflag = .false.
	    val = flag
	    do j=1,nintp
	      vals(j) = pinfo(id)%data(1,1,ivar,j)
	      if( vals(j) == flag ) bflag = .true.
	    end do
	    if( .not. bflag ) val = rd_intp_neville(nintp,time,vals,t)
	  end if
	  do i=1,nexp
	    do l=1,lexp
	      if( bmulti ) val = pinfo(id)%data(l,i,ivar,1)
	      value(l,i) = val
	      if( val == flag ) iflag = iflag + 1
	    end do
	  end do
	else
	  time = pinfo(id)%time
	  do i=1,nexp
	    if( b2d ) then
	      lfem = 1
	    else
	      ipl = i
	      if( nexp /= nkn_fem .and. nexp /= nel_fem ) then
		ipl = pinfo(id)%nodes(i)
	      end if
	      lfem = ilhkv_fem(ipl)
	    end if
	    val = -888.		!just for check
	    do l=1,lfem
	      bflag = .false.
	      val = flag
	      do j=1,nintp
	        vals(j) = pinfo(id)%data(l,i,ivar,j)
	        if( vals(j) == flag ) bflag = .true.
	      end do
	      if( .not. bflag ) val = rd_intp_neville(nintp,time,vals,t)
	      value(l,i) = val
	      if( val == flag ) iflag = iflag + 1
	    end do
	    do l=lfem+1,ldim
	      value(l,i) = val
	      if( val == flag ) iflag = iflag + 1
	    end do
	  end do
	end if

	if( iflag > 0 .and. pinfo(id)%bneedall ) then
	  write(6,*) 'flag values found: ',iflag
	  write(6,*) id,bconst,bonepoint
	  write(6,*) nintp, nexp,lexp
	  write(6,*) val,flag
	  call iff_print_file_info(id)
	  write(6,*) 'we need all values for interpolation'
	   stop 'error stop iff_interpolate: iflag'
	end if

	end subroutine iff_interpolate

!****************************************************************

	subroutine iff_adjust_datetime(id,datetime,dtime)

! here we have to transform fem_file time to simulation time

	integer id
	integer datetime(2)		!reference date of fem_file
	double precision dtime		!relative time of fem_file

	double precision atime0

	if( datetime(1) <= 0 ) return		!nothing to adjust

	if( datetime(1) > 0 .and. date_fem <= 0 ) then
	  write(6,*) 'file has absolute date but simulation has not'
	  write(6,*) 'please set the date variable in the STR file'
	  call iff_print_file_info(id)
	  stop 'error stop iff_adjust_datetime: date'
	end if

	call dts_to_abs_time(datetime(1),datetime(2),atime0)

	dtime = dtime + atime0 - atime0_fem

	end subroutine iff_adjust_datetime

!****************************************************************

	function tcomp(t,nintp,ilast,time)

	double precision tcomp
	double precision t
	integer nintp
	integer ilast
	double precision time(nintp)

	integer n1,n2

	!n1=1+nintp/2
	!n2 = n1
	!if( mod(nintp,2) /= 0 .and. nintp > 1 ) n2=n2+1		!odd

	if( nintp < 1 ) then
	  tcomp = t
	  return
	end if

	n1=1+mod(ilast+nintp/2,nintp)
	n2 = n1
	if( mod(nintp,2) /= 0 .and. nintp > 1 ) n2=mod(n2,nintp)+1	!odd

	tcomp = 0.5 * (time(n1)+time(n2))

	end function tcomp

!****************************************************************
!****************************************************************
!****************************************************************

	subroutine iff_get_value(id,ivar,iintp,l,k,val)

	integer id
	integer ivar
	integer iintp
	integer l
	integer k
	real val

	val =  pinfo(id)%data(l,k,ivar,iintp)

	end subroutine iff_get_value

!****************************************************************

	subroutine iff_get_file_value(id,ivar,l,k,val)

	integer id
	integer ivar
	integer l
	integer k
	real val

	val =  pinfo(id)%data_file(l,k,ivar)

	end subroutine iff_get_file_value

!****************************************************************

	subroutine iff_assert(bval,text)

	logical bval
	character*(*) text

	if( bassert ) then
	  if( .not. bval ) then
	    write(6,*) 'assertion violated'
	    write(6,*) text
	    stop 'error stop iff_assert: assertion violated'
	  end if
	end if

	end subroutine iff_assert

!****************************************************************

	subroutine iff_need_all_values(id,bneedall)

	integer id
	logical bneedall

	pinfo(id)%bneedall = bneedall

	end subroutine iff_need_all_values

!****************************************************************

	subroutine iff_flag_ok(id)

	integer id

	call iff_need_all_values(id,.false.)

	end subroutine iff_flag_ok

!****************************************************************

	subroutine iff_get_flag(id,flag)

	integer id
	real flag

	flag = pinfo(id)%flag

	end subroutine iff_get_flag

!================================================================
	end module intp_fem_file
!================================================================

!****************************************************************
! next are dummy routines that can be deleted somewhen...
!****************************************************************

        subroutine exffil(file,nintp,nvar,nsize,ndim,array)

	use intp_fem_file

c opens file and inititializes array
c
c everything needed is in array (unit, vars etc...)

        implicit none

        character*(*) file      !file name
        integer nintp           !grade of interpolation (2=linear,4=cubic)
        integer nvar            !how many vars (columns) to read/interpolate
        integer nsize           !number of data per variable (may be 0 -> 1)
        integer ndim            !dimension of array
        real array(ndim)        !array with all information

	integer nv,id
	integer nexp,lexp
	integer nodes(1)
	double precision dtime
	real vconst(nvar)

	dtime = -1
	nexp = 1
	lexp = 0
	vconst = 0.
	nodes = 0
	nv = nvar
	if( ndim < 2 ) stop 'error stop exffil: ndim'
	
	call iff_init(dtime,file,nvar,nexp,lexp,nintp
     +					,nodes,vconst,id)

	if( nv .ne. nvar ) then
	  write(6,*) 'nvar,nv: ',nvar,nv
	  stop 'error stop exffil: parameter mismatch'
	end if

	array(1) = id
	array(2) = nvar

	end

!****************************************************************

        subroutine exffils(file,ndim,array)

c opens file and inititializes array - simplified version

        implicit none

        character*(*) file      !file name
        integer ndim            !dimension of array
        real array(ndim)        !array with all information

        integer nintp           !grade of interpolation (2=linear,4=cubic)
        integer nvar            !how many columns to read/interpolate
        integer nsize           !number of data per variable

        nintp=2
        nvar=1
        nsize=0

        call exffil(file,nintp,nvar,nsize,ndim,array)

        end

!****************************************************************

        subroutine exfintp(array,t,rint)

	use intp_fem_file

	implicit none

        real array(*)           !array with information from set-up
        real t                  !t value for which to interpolate
        real rint(*)            !interpolated values

	integer id,ldim,ndim,ivar,nvar
	double precision dtime

	id = nint(array(1))
	nvar = nint(array(2))
	dtime = t
	ldim = 1
	ndim = 1

	do ivar=1,nvar
	  call iff_time_interpolate(id,dtime,ivar,ndim,ldim,rint(ivar))
	end do

        end

!****************************************************************

	subroutine iff_init_global_2d(nkn,nel,hkv,hev,date,time)

	use intp_fem_file

	implicit none

	integer nkn,nel
	real hkv(nkn)
	real hev(nel)
	integer date,time

	integer nlv
	integer ilhkv(nkn)
	integer ilhv(nel)
	real hlv(1)

	nlv = 1
	ilhkv = 1
	ilhv = 1
	hlv(1) = 10000.

	call iff_init_global(nkn,nel,nlv,ilhkv,ilhv
     +				,hkv,hev,hlv,date,time)

	end

!****************************************************************

	subroutine iff_init_global_date(date,time)

	use intp_fem_file

	implicit none

	integer date,time

	call iff_init_global_date_internal(date,time)

	end

!****************************************************************

