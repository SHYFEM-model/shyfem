!
! SHY management routines
!
! contents :
!
! revision log :
!
! 10.10.2015    ggu     started routine
!
!**************************************************************
!**************************************************************
!**************************************************************

!==================================================================
	module shyfile
!==================================================================

	implicit none

	integer, parameter, private :: shytype = 1617

	integer, parameter, private :: minvers = 11
	integer, parameter, private :: maxvers = 11

	integer, parameter, private ::  no_type = 0
	integer, parameter, private :: ous_type = 1
	integer, parameter, private :: nos_type = 2
	integer, parameter, private :: ext_type = 3
	integer, parameter, private :: flx_type = 4

	type, private :: entry

	  integer :: iunit
	  integer :: nvers
	  integer :: ftype
	  integer :: nkn,nel,nlv,nvar
	  integer :: date,time
	  character*80 :: title
	  character*80 :: femver
          integer, allocatable :: nen3v(:,:)
          integer, allocatable :: ipev(:)
          integer, allocatable :: ipv(:)
          integer, allocatable :: iarv(:)
          integer, allocatable :: iarnv(:)
          real, allocatable :: xgv(:)
          real, allocatable :: ygv(:)
          real, allocatable :: hm3v(:,:)
          integer, allocatable :: hlv(:)
          integer, allocatable :: ilhv(:)
          integer, allocatable :: ilhkv(:)

	  logical :: is_allocated
	end type entry

	integer, save, private :: idlast = 0
	integer, save, private :: ndim = 0
	type(entry), save, allocatable :: pentry(:)

!==================================================================
	contains
!==================================================================

        subroutine shy_init_alloc

        type(entry), allocatable :: paux(:)

        if( ndim == 0 ) then
          ndim = 10
          allocate(pentry(ndim))
          return
        else
          ndim = ndim*2
          allocate(paux(ndim))
          paux(1:ndim/2) = pentry(1:ndim/2)
          call move_alloc(paux,pentry)
        end if

        end subroutine shy_init_alloc

!******************************************************************

	subroutine shy_init_new_id(id)

	integer id

	idlast = idlast + 1
	if( idlast > ndim ) then
          call shy_init_alloc
	end if
	id = idlast

	call shy_init_id(id)
	
	end subroutine shy_init_new_id

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shy_init_id(id)

	integer id

	if( id > ndim ) then
	  stop 'error stop shy_init_id: ndim'
	end if

	pentry(id)%iunit = 0
	pentry(id)%nvers = 0
	pentry(id)%ftype = 0
	pentry(id)%nkn = 0
	pentry(id)%nel = 0
	pentry(id)%nlv = 0
	pentry(id)%nvar = 0
	pentry(id)%date = 0
	pentry(id)%time = 0
	pentry(id)%title = ' '
	pentry(id)%femver = ' '

	pentry(id)%is_allocated = .false.
	
	end subroutine shy_init_id

!******************************************************************

	subroutine shy_alloc_arrays(id)

	integer id

	integer nkn,nel,nlv

	if( pentry(id)%is_allocated ) return

	nkn = pentry(id)%nkn
	nel = pentry(id)%nel
	nlv = pentry(id)%nlv

	allocate(pentry(id)%nen3v(3,nel))
	allocate(pentry(id)%ipev(nel))
	allocate(pentry(id)%ipv(nkn))
	allocate(pentry(id)%iarv(nel))
	allocate(pentry(id)%iarnv(nkn))
	allocate(pentry(id)%xgv(nkn))
	allocate(pentry(id)%ygv(nkn))
	allocate(pentry(id)%hm3v(3,nel))
	allocate(pentry(id)%hlv(nlv))
	allocate(pentry(id)%ilhv(nel))
	allocate(pentry(id)%ilhkv(nkn))

	pentry(id)%is_allocated = .true.

	end subroutine shy_alloc_arrays

!******************************************************************

	subroutine shy_dealloc_arrays(id)

	integer id

	if( .not. pentry(id)%is_allocated ) return

	deallocate(pentry(id)%nen3v)
	deallocate(pentry(id)%ipev)
	deallocate(pentry(id)%ipv)
	deallocate(pentry(id)%iarv)
	deallocate(pentry(id)%iarnv)
	deallocate(pentry(id)%xgv)
	deallocate(pentry(id)%ygv)
	deallocate(pentry(id)%hm3v)
	deallocate(pentry(id)%hlv)
	deallocate(pentry(id)%ilhv)
	deallocate(pentry(id)%ilhkv)

	pentry(id)%is_allocated = .false.

	end subroutine shy_dealloc_arrays

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shy_get_file_unit(iu)

	integer iu

	integer, parameter :: iu_min = 20
	integer, parameter :: iu_max = 1000

	logical bopen
	integer ios

	do iu=iu_min,iu_max
          inquire (unit=iu, opened=bopen, iostat=ios)
          if ( ios /= 0 ) cycle
          if ( .not. bopen ) return
	end do

	iu = 0

	end subroutine shy_get_file_unit

!******************************************************************

	subroutine shy_init(iunit,id)

	integer iunit
	integer nversion
	integer id

	integer nvers
	integer idempty

	if( iunit .le. 0 ) then
	  write(6,*) 'Impossible unit: ',iunit
	  stop 'error stop shy_init: iunit'
	end if

	nvers = maxvers

	idempty = 0
	do id=1,idlast
	  if( pentry(id)%iunit == 0 ) idempty = id
	  if( pentry(id)%iunit == iunit ) then
	    write(6,*) 'unit already initialized: ',iunit
	    stop 'error stop shy_init: iunit used'
	  end if
	end do

	id = idempty
	if( id == 0 ) then
	  call shy_init_new_id(id)
	else
	  call shy_init_id(id)
	end if

	pentry(id)%iunit = iunit
	pentry(id)%nvers = nvers

	rewind(iunit)

	end subroutine shy_init

c************************************************************

	subroutine shy_close(id)

	integer id

	pentry(id)%iunit = 0
	call shy_dealloc_arrays(id)
	if( id == idlast ) idlast = idlast - 1

	end subroutine shy_close

c************************************************************

	function shy_are_compatible(id1,id2)

	logical shy_are_compatible
	integer id1,id2

	shy_are_compatible = .false.

	if( pentry(id1)%nkn /= pentry(id2)%nkn ) return 
	if( pentry(id1)%nel /= pentry(id2)%nel ) return 
	if( pentry(id1)%nlv /= pentry(id2)%nlv ) return 
	if( pentry(id1)%nvar /= pentry(id2)%nvar ) return 
	if( pentry(id1)%date /= pentry(id2)%date ) return 
	if( pentry(id1)%time /= pentry(id2)%time ) return 

	if( .not. all(pentry(id1)%nen3v==pentry(id2)%nen3v) ) return
	if( .not. all(pentry(id1)%xgv==pentry(id2)%xgv) ) return
	if( .not. all(pentry(id1)%ygv==pentry(id2)%ygv) ) return
	if( .not. all(pentry(id1)%hm3v==pentry(id2)%hm3v) ) return

	!still to be finished
	!deallocate(pentry(id)%ipev)
	!deallocate(pentry(id)%ipv)
	!deallocate(pentry(id)%iarv)
	!deallocate(pentry(id)%iarnv)
	!deallocate(pentry(id)%hlv)
	!deallocate(pentry(id)%ilhv)
	!deallocate(pentry(id)%ilhkv)

	shy_are_compatible = .true.

	end function shy_are_compatible

c************************************************************

	subroutine shy_clone(id_from,id_to)

	integer id_from,id_to

	integer iunit,nvers

	iunit = pentry(id_to)%iunit
	nvers = pentry(id_to)%nvers

	pentry(id_to) = pentry(id_from)

	pentry(id_to)%iunit = iunit
	pentry(id_to)%nvers = nvers

	end subroutine shy_clone

c************************************************************

	subroutine shy_info(id)

	integer id

        write(6,*) 'iunit     : ',pentry(id)%iunit
        write(6,*) 'nvers     : ',pentry(id)%nvers
        write(6,*) 'nkn,nel   : ',pentry(id)%nkn,pentry(id)%nel
        write(6,*) 'nlv,nvar  : ',pentry(id)%nlv,pentry(id)%nvar
        write(6,*) 'date,time : ',pentry(id)%date,pentry(id)%time
        write(6,*) 'title     : ',pentry(id)%title
        write(6,*) 'femver    : ',pentry(id)%femver

	end subroutine shy_info

c************************************************************
c************************************************************
c************************************************************

	function shy_is_shy_file(iunit)

	logical shy_is_shy_file
	integer iunit

	integer ntype,nvers,ios

	shy_is_shy_file = .false.
	if( iunit .le. 0 ) return

	read(iunit,iostat=ios) ntype,nvers

	if( ios /= 0 ) return
	if( ntype .ne. shytype ) return
	if( nvers .lt. minvers .or. nvers .gt. maxvers ) return

	shy_is_shy_file = .true.
	rewind(iunit)

	end function shy_is_shy_file

c************************************************************
c************************************************************
c************************************************************

	subroutine shy_get_ftype(id,ftype)
	integer id
	integer ftype
	ftype = pentry(id)%ftype
	end subroutine shy_get_ftype

	subroutine shy_set_ftype(id,ftype)
	integer id
	integer ftype
	pentry(id)%ftype = ftype
	end subroutine shy_set_ftype

c************************************************************

	subroutine shy_get_params(id,nkn,nel,nlv,nvar)
	integer id
	integer nkn,nel,nlv,nvar
	nkn = pentry(id)%nkn
	nel = pentry(id)%nel
	nlv = pentry(id)%nlv
	nvar = pentry(id)%nvar
	end subroutine shy_get_params

	subroutine shy_set_params(id,nkn,nel,nlv,nvar)
	integer id
	integer nkn,nel,nlv,nvar
	pentry(id)%nkn = nkn
	pentry(id)%nel = nel
	pentry(id)%nlv = nlv
	pentry(id)%nvar = nvar
	end subroutine shy_set_params

c************************************************************

	subroutine shy_get_date(id,date,time)
	integer id
	integer date,time
	date = pentry(id)%date
	time = pentry(id)%time
	end subroutine shy_get_date

	subroutine shy_set_date(id,date,time)
	integer id
	integer date,time
	pentry(id)%date = date
	pentry(id)%time = time
	end subroutine shy_set_date

c************************************************************

	subroutine shy_get_title(id,title)
	integer id
	character*(*) title
	title = pentry(id)%title
	end subroutine shy_get_title

	subroutine shy_set_title(id,title)
	integer id
	character*(*) title
	pentry(id)%title = title
	end subroutine shy_set_title

c************************************************************

	subroutine shy_get_femver(id,femver)
	integer id
	character*(*) femver
	femver = pentry(id)%femver
	end subroutine shy_get_femver

	subroutine shy_set_femver(id,femver)
	integer id
	character*(*) femver
	pentry(id)%femver = femver
	end subroutine shy_set_femver

c************************************************************

	subroutine shy_get_elemindex(id,nen3v)
	integer id
	integer nen3v(3,pentry(id)%nel)
	nen3v = pentry(id)%nen3v
	end subroutine shy_get_elemindex

	subroutine shy_set_elemindex(id,nen3v)
	integer id
	integer nen3v(3,pentry(id)%nel)
	pentry(id)%nen3v = nen3v
	end subroutine shy_set_elemindex

c************************************************************

	subroutine shy_get_coords(id,xgv,ygv)
	integer id
	integer xgv(pentry(id)%nkn), ygv(pentry(id)%nkn)
	xgv = pentry(id)%xgv
	ygv = pentry(id)%ygv
	end subroutine shy_get_coords

	subroutine shy_set_coords(id,xgv,ygv)
	integer id
	integer xgv(pentry(id)%nkn), ygv(pentry(id)%nkn)
	pentry(id)%xgv = xgv
	pentry(id)%ygv = ygv
	end subroutine shy_set_coords

c************************************************************

	subroutine shy_get_depth(id,hm3v)
	integer id
	integer hm3v(3,pentry(id)%nel)
	hm3v = pentry(id)%hm3v
	end subroutine shy_get_depth

	subroutine shy_set_depth(id,hm3v)
	integer id
	integer hm3v(3,pentry(id)%nel)
	pentry(id)%hm3v = hm3v
	end subroutine shy_set_depth

c************************************************************

	subroutine shy_get_layers(id,hlv)
	integer id
	integer hlv(pentry(id)%nlv)
	hlv = pentry(id)%hlv
	end subroutine shy_get_layers

	subroutine shy_set_layers(id,hlv)
	integer id
	integer hlv(pentry(id)%nlv)
	pentry(id)%hlv = hlv
	end subroutine shy_set_layers

c************************************************************

	subroutine shy_get_layerindex(id,ilhv,ilhkv)
	integer id
	integer ilhv(pentry(id)%nel), ilhkv(pentry(id)%nkn)
	ilhv = pentry(id)%ilhv
	ilhkv = pentry(id)%ilhkv
	end subroutine shy_get_layerindex

	subroutine shy_set_layerindex(id,ilhv,ilhkv)
	integer id
	integer ilhv(pentry(id)%nel), ilhkv(pentry(id)%nkn)
	pentry(id)%ilhv = ilhv
	pentry(id)%ilhkv = ilhkv
	end subroutine shy_set_layerindex

c************************************************************
c************************************************************
c************************************************************

	subroutine shy_read_header(id,ierr)

	integer id,ierr

	call shy_read_header_1(id,ierr)
	if( ierr /= 0 ) return
	call shy_read_header_2(id,ierr)

	end subroutine shy_read_header

c************************************************************

	subroutine shy_peek_header(id,ierr)

	integer id,ierr

	call shy_read_header_1(id,ierr)
	if( ierr /= 0 ) return
	rewind(pentry(id)%iunit)

	end subroutine shy_peek_header

c************************************************************

	subroutine shy_skip_header(id,ierr)

	integer id,ierr

	call shy_read_header_1(id,ierr)
	if( ierr /= 0 ) return
	call shy_skip_header_2(id,ierr)

	end subroutine shy_skip_header

c************************************************************

	subroutine shy_read_header_1(id,ierr)

	integer id,ierr

	integer ios,iunit
	integer ntype,nvers
	integer ftype
	integer nkn,nel,nlv,nvar
	integer date,time
	character*80 title
	character*80 femver

	iunit = pentry(id)%iunit

	ierr = 1
        read(iunit,iostat=ios) ntype,nvers
	if( ios /= 0 ) return

	ierr = 91
	if( ntype /= shytype ) return
	ierr = 92
	if( nvers < minvers .or. nvers > maxvers ) return
	pentry(id)%nvers = nvers

	ierr = 2
        read(iunit,iostat=ios) ftype
	if( ios /= 0 ) return
	call shy_set_ftype(id,ftype)

	ierr = 3
        read(iunit,iostat=ios) nkn,nel,nlv,nvar
	if( ios /= 0 ) return
	call shy_set_params(id,nkn,nel,nlv,nvar)

	ierr = 4
        read(iunit,iostat=ios) date,time
	if( ios /= 0 ) return
	call shy_set_date(id,date,time)

	ierr = 5
        read(iunit,iostat=ios) title
	if( ios /= 0 ) return
	call shy_set_title(id,title)

	ierr = 6
        read(iunit,iostat=ios) femver
	if( ios /= 0 ) return
	call shy_set_femver(id,femver)

	ierr = 0

	end subroutine shy_read_header_1

!**************************************************************

	subroutine shy_read_header_2(id,ierr)

	integer id,ierr

	integer iunit

	iunit = pentry(id)%iunit

        read(iunit,err=99) pentry(id)%nen3v
        read(iunit,err=99) pentry(id)%ipev
        read(iunit,err=99) pentry(id)%ipv
        read(iunit,err=99) pentry(id)%iarv
        read(iunit,err=99) pentry(id)%iarnv
        read(iunit,err=99) pentry(id)%xgv
        read(iunit,err=99) pentry(id)%ygv
        read(iunit,err=99) pentry(id)%hm3v
        read(iunit,err=99) pentry(id)%hlv
        read(iunit,err=99) pentry(id)%ilhv
        read(iunit,err=99) pentry(id)%ilhkv

	ierr = 0

	return
   99	continue
	ierr = 21
	return
	end subroutine shy_read_header_2

!**************************************************************

	subroutine shy_skip_header_2(id,ierr)

	integer id,ierr

	integer iunit,i

	iunit = pentry(id)%iunit

	do i=1,11
          read(iunit,err=99) 
	end do

	ierr = 0

	return
   99	continue
	ierr = 31
	return
	end subroutine shy_skip_header_2

!**************************************************************

	subroutine shy_read_record(id,dtime,ivar,n,m,lmax,nlvddi,c,ierr)

	integer id,ierr
	double precision dtime
	integer ivar
	integer n,m
	integer lmax
	integer nlvddi
	real c(nlvddi,*)

	integer iunit
	integer i,k,ie,l

	iunit = pentry(id)%iunit

	read(iunit,iostat=ierr) dtime,ivar,n,m,lmax
	if( ierr /= 0 ) return

	if( lmax <= 1 ) then
	  read(iunit,iostat=ierr) (c(1,i),i=1,n*m)
	else if( n == pentry(id)%nkn ) then
	  read(iunit,iostat=ierr) ((c(l,k)
     +			,l=1,m*pentry(id)%ilhkv(k)),k=1,n)
	else if( n == pentry(id)%nel ) then
	  read(iunit,iostat=ierr) ((c(l,ie)
     +			,l=1,m*pentry(id)%ilhv(ie)),ie=1,n)
	end if

	end subroutine shy_read_record

!**************************************************************

	subroutine shy_peek_record(id,dtime,ivar,n,m,lmax,ierr)

	integer id,ierr
	double precision dtime
	integer ivar
	integer n,m
	integer lmax

	integer iunit

	iunit = pentry(id)%iunit

	read(iunit,iostat=ierr) dtime,ivar,n,m,lmax
	if( ierr /= 0 ) return
	rewind(iunit,iostat=ierr)

	end subroutine shy_peek_record

!**************************************************************

	subroutine shy_skip_record(id,dtime,ivar,n,m,lmax,ierr)

	integer id,ierr
	double precision dtime
	integer ivar
	integer n,m
	integer lmax

	integer iunit

	iunit = pentry(id)%iunit

	read(iunit,iostat=ierr) dtime,ivar,n,m,lmax
	if( ierr /= 0 ) return
	read(iunit,iostat=ierr)

	end subroutine shy_skip_record

!**************************************************************

	subroutine shy_back_record(id,ierr)

	integer id,ierr

	integer iunit

	iunit = pentry(id)%iunit

	backspace(iunit,iostat=ierr)
	if( ierr /= 0 ) return
	backspace(iunit,iostat=ierr)

	end subroutine shy_back_record

!**************************************************************
!**************************************************************
!**************************************************************

	subroutine shy_write_header(id,ierr)

	integer id,ierr

	integer ios,iunit
	integer ntype,nvers
	integer ftype
	integer nkn,nel,nlv,nvar
	integer date,time
	character*80 title
	character*80 femver

	iunit = pentry(id)%iunit
	nvers = maxvers

	call shy_get_ftype(id,ftype)
	call shy_get_params(id,nkn,nel,nlv,nvar)
	call shy_get_date(id,date,time)
	call shy_get_title(id,title)
	call shy_get_femver(id,femver)

        write(iunit,err=99) shytype,nvers
        write(iunit,err=99) ftype
        write(iunit,err=99) nkn,nel,nlv,nvar
        write(iunit,err=99) date,time
        write(iunit,err=99) title
        write(iunit,err=99) femver

        write(iunit,err=99) pentry(id)%nen3v
        write(iunit,err=99) pentry(id)%ipev
        write(iunit,err=99) pentry(id)%ipv
        write(iunit,err=99) pentry(id)%iarv
        write(iunit,err=99) pentry(id)%iarnv
        write(iunit,err=99) pentry(id)%xgv
        write(iunit,err=99) pentry(id)%ygv
        write(iunit,err=99) pentry(id)%hm3v
        write(iunit,err=99) pentry(id)%hlv
        write(iunit,err=99) pentry(id)%ilhv
        write(iunit,err=99) pentry(id)%ilhkv

	ierr = 0

	return
   99	continue
	ierr = 51
	return
	end subroutine shy_write_header

!**************************************************************

	subroutine shy_write_record(id,dtime,ivar,n,m,lmax,nlvddi,c,ierr)

	integer id,ierr
	double precision dtime
	integer ivar
	integer n,m
	integer lmax
	integer nlvddi
	real c(nlvddi,*)

	integer iunit
	integer i,k,ie,l

	iunit = pentry(id)%iunit

	write(iunit,iostat=ierr) dtime,ivar,n,m,lmax
	if( ierr /= 0 ) return

	if( lmax <= 1 ) then
	  write(iunit,iostat=ierr) (c(1,i),i=1,n*m)
	else if( n == pentry(id)%nkn ) then
	  write(iunit,iostat=ierr) ((c(l,k)
     +			,l=1,m*pentry(id)%ilhkv(k)),k=1,n)
	else if( n == pentry(id)%nel ) then
	  write(iunit,iostat=ierr) ((c(l,ie)
     +			,l=1,m*pentry(id)%ilhv(ie)),ie=1,n)
	end if

	end subroutine shy_write_record

!**************************************************************


!**************************************************************

!==================================================================
	end module shyfile
!==================================================================


!**************************************************************
!**************************************************************
!**************************************************************


c************************************************************
c************************************************************
c************************************************************

	subroutine test_units

	implicit none

	integer iu,iumax,ios

	iumax = 1000000
	iu = 20
	iu = 0

	do
	  iu = iu + 100
	  open(iu,iostat=ios)
	  if( ios /= 0 ) exit
	  close(iu)
	  if( mod(iu,1000) .eq. 0 )  write(6,*) iu
	  if( iu >= iumax ) exit
	end do

	write(6,*) 'last unit tested: ',iu

	end

c************************************************************
c************************************************************
c************************************************************

	program shy_main
	!call test_units
	end

c************************************************************
c************************************************************
c************************************************************

