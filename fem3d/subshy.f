!
! SHY management routines
!
! contents :
!
! revision log :
!
! 10.10.2015    ggu     started routine
! 15.10.2015    ggu     completed basic parts
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
	integer, parameter, private :: eos_type = 3
	integer, parameter, private :: ext_type = 4
	integer, parameter, private :: flx_type = 5

	type, private :: entry

	  integer :: iunit
	  integer :: nvers
	  integer :: ftype
	  integer :: nkn,nel,npr,nlv,nvar
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
          real, allocatable :: hlv(:)
          integer, allocatable :: ilhv(:)
          integer, allocatable :: ilhkv(:)

	  logical :: is_allocated
	end type entry

	integer, save, private :: idlast = 0
	integer, save, private :: ndim = 0
	type(entry), save, allocatable :: pentry(:)

        INTERFACE shy_is_shy_file
        MODULE PROCEDURE shy_is_shy_file_by_name,shy_is_shy_file_by_unit
        END INTERFACE

        INTERFACE shy_init
        MODULE PROCEDURE shy_init_by_unit,shy_init_by_file
        END INTERFACE

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
	pentry(id)%npr = 0
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

	function shy_open_file(file,status)

	integer shy_open_file
	character*(*) file
	character*(*), optional :: status

	integer iunit,ios
	character*20 stat

	stat = 'unknown'
	if( present(status) ) stat = status

	shy_open_file = 0

	call shy_get_file_unit(iunit)
	if( iunit == 0 ) return

	open(iunit,file=file,status=stat,form='unformatted'
     +				,iostat=ios)
	if( ios /= 0 ) return

	shy_open_file = iunit

	end function shy_open_file

!******************************************************************

	subroutine shy_get_file_unit(iunit)

	integer iunit

	integer, parameter :: iu_min = 20
	integer, parameter :: iu_max = 1000

	logical bopen
	integer ios

	do iunit=iu_min,iu_max
          inquire (unit=iunit, opened=bopen, iostat=ios)
          if ( ios /= 0 ) cycle
          if ( .not. bopen ) return
	end do

	iunit = 0

	end subroutine shy_get_file_unit

!******************************************************************
!******************************************************************
!******************************************************************

	function shy_init_by_file(file)

	integer shy_init_by_file
	character*(*) file

	integer iunit

	shy_init_by_file = 0

	iunit = shy_open_file(file)
	if( iunit == 0 ) return

	shy_init_by_file = shy_init_by_unit(iunit)

	end function shy_init_by_file

!******************************************************************

	function shy_init_by_unit(iunit)

	integer shy_init_by_unit
	integer iunit

	integer id
	integer nvers
	integer idempty

	if( iunit .le. 0 ) then
	  write(6,*) 'Impossible unit: ',iunit
	  stop 'error stop shy_init_by_unit: iunit'
	end if

	nvers = maxvers
	shy_init_by_unit = 0

	idempty = 0
	do id=1,idlast
	  if( pentry(id)%iunit == 0 ) idempty = id
	  if( pentry(id)%iunit == iunit ) then
	    write(6,*) 'unit already initialized: ',iunit
	    stop 'error stop shy_init_by_unit: iunit used'
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

	shy_init_by_unit = id

	end function shy_init_by_unit

!************************************************************

	subroutine shy_close(id)

	integer id

	if( id <= 0 ) return

	close(pentry(id)%iunit)
	pentry(id)%iunit = 0
	call shy_dealloc_arrays(id)
	if( id == idlast ) idlast = idlast - 1

	end subroutine shy_close

!************************************************************

	function shy_are_compatible(id1,id2)

	logical shy_are_compatible
	integer id1,id2

	integer nk,ne
	character*20 serr

	shy_are_compatible = .false.
	nk = pentry(id1)%nkn
	ne = pentry(id1)%nel

	serr = 'parameters'
	if( pentry(id1)%nkn /= pentry(id2)%nkn ) goto 99 
	if( pentry(id1)%nel /= pentry(id2)%nel ) goto 99 
	if( pentry(id1)%nlv /= pentry(id2)%nlv ) goto 99 
	if( pentry(id1)%npr /= pentry(id2)%npr ) goto 99 
	if( pentry(id1)%nvar /= pentry(id2)%nvar ) goto 99 
	if( pentry(id1)%date /= pentry(id2)%date ) goto 99 
	if( pentry(id1)%time /= pentry(id2)%time ) goto 99 

	serr = 'nen3v'
	if( .not. all(pentry(id1)%nen3v==pentry(id2)%nen3v) ) goto 99
	serr = 'xgv'
	if( .not. all(pentry(id1)%xgv==pentry(id2)%xgv) ) goto 99
	serr = 'ygv'
	if( .not. all(pentry(id1)%ygv==pentry(id2)%ygv) ) goto 99
	serr = 'hm3v'
	!call shy_diff_internal_r(3*ne,pentry(id1)%hm3v,pentry(id2)%hm3v)
	!if( .not. all(pentry(id1)%hm3v==pentry(id2)%hm3v) ) goto 99

	!still to be finished
	!deallocate(pentry(id)%ipev)
	!deallocate(pentry(id)%ipv)
	!deallocate(pentry(id)%iarv)
	!deallocate(pentry(id)%iarnv)
	serr = 'hlv'
	if( .not. all(pentry(id1)%hlv==pentry(id2)%hlv) ) goto 99
	!deallocate(pentry(id)%hlv)
	!deallocate(pentry(id)%ilhv)
	!deallocate(pentry(id)%ilhkv)

	shy_are_compatible = .true.

	return
  99	continue
	write(6,*) 'error checking ',serr
	return
	end function shy_are_compatible

!************************************************************

	subroutine shy_diff_internal_r(n,r1,r2)
	integer n
	real r1(n),r2(n)
	integer i
	do i=1,n
	  write(6,*) r1(i),r2(i),r2(i)-r1(i)
	end do
	end

!************************************************************

	subroutine shy_clone(id_from,id_to)

	integer id_from,id_to

	integer iunit,nvers

	iunit = pentry(id_to)%iunit
	nvers = pentry(id_to)%nvers

	pentry(id_to) = pentry(id_from)

	pentry(id_to)%iunit = iunit
	pentry(id_to)%nvers = nvers

	end subroutine shy_clone

!************************************************************

	subroutine shy_info(id)

	integer id

	character*80 file

	call shy_get_filename(id,file)

        write(6,*) 'filename  : ',trim(file)
        write(6,*) 'iunit     : ',pentry(id)%iunit
        write(6,*) 'nvers     : ',pentry(id)%nvers
        write(6,*) 'nkn,nel   : ',pentry(id)%nkn,pentry(id)%nel
        write(6,*) 'npr,nlv   : ',pentry(id)%npr,pentry(id)%nlv
        write(6,*) 'nvar      : ',pentry(id)%nvar
        write(6,*) 'date,time : ',pentry(id)%date,pentry(id)%time
        write(6,*) 'title     : ',trim(pentry(id)%title)
        write(6,*) 'femver    : ',trim(pentry(id)%femver)

        write(6,*) 'hlv       : ',pentry(id)%hlv

	end subroutine shy_info

!************************************************************

	subroutine shy_get_filename(id,file)

	integer id
	character*(*) file

	integer iunit,ios

	file = ' '
	iunit = pentry(id)%iunit
	if( iunit == 0 ) return

	inquire(iunit,name=file,iostat=ios)

	end subroutine shy_get_filename

!************************************************************
!************************************************************
!************************************************************

	function shy_exist_file(file)

	logical shy_exist_file
	character*(*) file

	integer iunit

	shy_exist_file = .false.

	iunit = shy_open_file(file,'old')
	if( iunit .le. 0 ) return

	shy_exist_file = .true.
	close(iunit)

	end function shy_exist_file

!************************************************************

	function shy_is_shy_file_by_name(file)

	logical shy_is_shy_file_by_name
	character*(*) file

	integer iunit

	shy_is_shy_file_by_name = .false.

	iunit = shy_open_file(file)
	if( iunit .le. 0 ) return

	shy_is_shy_file_by_name = shy_is_shy_file_by_unit(iunit)
	close(iunit)

	end function shy_is_shy_file_by_name

!************************************************************

	function shy_is_shy_file_by_unit(iunit)

	logical shy_is_shy_file_by_unit
	integer iunit

	integer ntype,nvers,ios

	shy_is_shy_file_by_unit = .false.
	if( iunit .le. 0 ) return

	read(iunit,iostat=ios) ntype,nvers

	!write(6,*) 'shy_is_shy_file: ',ios,ntype,nvers

	if( ios /= 0 ) return
	if( ntype .ne. shytype ) return
	if( nvers .lt. minvers .or. nvers .gt. maxvers ) return

	shy_is_shy_file_by_unit = .true.
	rewind(iunit)

	end function shy_is_shy_file_by_unit

!************************************************************
!************************************************************
!************************************************************

	subroutine shy_convert_2d(id)

	integer id

	pentry(id)%nlv = 1

	deallocate(pentry(id)%hlv)
	allocate(pentry(id)%hlv(1))

	pentry(id)%hlv(1) = 10000.
	pentry(id)%ilhv = 1
	pentry(id)%ilhkv = 1

	end subroutine shy_convert_2d

!************************************************************

	subroutine shy_convert_1var(id)

	integer id

	pentry(id)%nvar = 1

	end subroutine shy_convert_1var

!************************************************************
!************************************************************

!************************************************************
!************************************************************
!************************************************************

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

!************************************************************

	subroutine shy_get_params(id,nkn,nel,npr,nlv,nvar)
	integer id
	integer nkn,nel,npr,nlv,nvar
	nkn = pentry(id)%nkn
	nel = pentry(id)%nel
	npr = pentry(id)%npr
	nlv = pentry(id)%nlv
	nvar = pentry(id)%nvar
	end subroutine shy_get_params

	subroutine shy_set_params(id,nkn,nel,npr,nlv,nvar)
	integer id
	integer nkn,nel,npr,nlv,nvar
	pentry(id)%nkn = nkn
	pentry(id)%nel = nel
	pentry(id)%npr = npr
	pentry(id)%nlv = nlv
	pentry(id)%nvar = nvar
	end subroutine shy_set_params

!************************************************************

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

!************************************************************

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

!************************************************************

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

!************************************************************

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

!************************************************************

	subroutine shy_get_coords(id,xgv,ygv)
	integer id
	real xgv(pentry(id)%nkn), ygv(pentry(id)%nkn)
	xgv = pentry(id)%xgv
	ygv = pentry(id)%ygv
	end subroutine shy_get_coords

	subroutine shy_set_coords(id,xgv,ygv)
	integer id
	real xgv(pentry(id)%nkn), ygv(pentry(id)%nkn)
	pentry(id)%xgv = xgv
	pentry(id)%ygv = ygv
	end subroutine shy_set_coords

!************************************************************

	subroutine shy_get_depth(id,hm3v)
	integer id
	real hm3v(3,pentry(id)%nel)
	hm3v = pentry(id)%hm3v
	end subroutine shy_get_depth

	subroutine shy_set_depth(id,hm3v)
	integer id
	real hm3v(3,pentry(id)%nel)
	pentry(id)%hm3v = hm3v
	end subroutine shy_set_depth

!************************************************************

	subroutine shy_get_layers(id,hlv)
	integer id
	real hlv(pentry(id)%nlv)
	hlv = pentry(id)%hlv
	end subroutine shy_get_layers

	subroutine shy_set_layers(id,hlv)
	integer id
	real hlv(pentry(id)%nlv)
	pentry(id)%hlv = hlv
	end subroutine shy_set_layers

!************************************************************

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

!************************************************************

	subroutine shy_get_extnumbers(id,ipev,ipv)
	integer id
	integer ipev(pentry(id)%nel), ipv(pentry(id)%nkn)
	ipev = pentry(id)%ipev
	ipv = pentry(id)%ipv
	end subroutine shy_get_extnumbers

	subroutine shy_set_extnumbers(id,ipev,ipv)
	integer id
	integer ipev(pentry(id)%nel), ipv(pentry(id)%nkn)
	pentry(id)%ipev = ipev
	pentry(id)%ipv = ipv
	end subroutine shy_set_extnumbers

!************************************************************

	subroutine shy_get_areacode(id,iarv,iarnv)
	integer id
	integer iarv(pentry(id)%nel), iarnv(pentry(id)%nkn)
	iarv = pentry(id)%iarv
	iarnv = pentry(id)%iarnv
	end subroutine shy_get_areacode

	subroutine shy_set_areacode(id,iarv,iarnv)
	integer id
	integer iarv(pentry(id)%nel), iarnv(pentry(id)%nkn)
	pentry(id)%iarv = iarv
	pentry(id)%iarnv = iarnv
	end subroutine shy_set_areacode

!************************************************************
!************************************************************
!************************************************************

	subroutine shy_read_header(id,ierr)

	integer id,ierr

	call shy_read_header_1(id,ierr)
	if( ierr /= 0 ) return
	call shy_alloc_arrays(id)
	call shy_read_header_2(id,ierr)

	end subroutine shy_read_header

!************************************************************

	subroutine shy_peek_header(id,ierr)

	integer id,ierr

	call shy_read_header_1(id,ierr)
	if( ierr /= 0 ) return
	rewind(pentry(id)%iunit)

	end subroutine shy_peek_header

!************************************************************

	subroutine shy_skip_header(id,ierr)

	integer id,ierr

	call shy_read_header_1(id,ierr)
	if( ierr /= 0 ) return
	call shy_skip_header_2(id,ierr)

	end subroutine shy_skip_header

!************************************************************

	subroutine shy_read_header_1(id,ierr)

	integer id,ierr

	integer ios,iunit
	integer ntype,nvers
	integer ftype
	integer nkn,nel,npr,nlv,nvar
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
        read(iunit,iostat=ios) nkn,nel,npr,nlv,nvar
	if( ios /= 0 ) return
	call shy_set_params(id,nkn,nel,npr,nlv,nvar)

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
	integer i,k,ie,l,j
	integer, allocatable :: il(:)

	iunit = pentry(id)%iunit

	ierr = 55
	if( .not. pentry(id)%is_allocated ) return

	read(iunit,iostat=ierr) dtime,ivar,n,m,lmax
	if( ierr /= 0 ) return

	allocate(il(n))
	if( n == pentry(id)%nkn ) then
	  il = pentry(id)%ilhkv
	else if( n == pentry(id)%nel ) then
	  il = pentry(id)%ilhv
	else
	  write(6,*) n,pentry(id)%nkn,pentry(id)%nel
	  write(6,*) 'cannot determine layer pointer'
	  call shy_info(id)
	  stop 'error stop shy_read_record: layer pointer'
	end if

	!write(6,*) id,ivar,n,m,lmax

	if( lmax <= 1 ) then
	  read(iunit,iostat=ierr) ( c(1,i),i=1,n*m )
	else if( m == 1 ) then
	  read(iunit,iostat=ierr) (( c(l,i)
     +			,l=1,il(i) )
     +			,i=1,n )
	else
	  read(iunit,iostat=ierr) (( c(l,i)
     +			,l=1,il(1+(i-1)/m) )
     +			,i=1,n*m )
	end if
	deallocate(il)

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
	if( ierr > 0 ) return
	backspace(iunit,iostat=ierr)

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
	integer nkn,nel,npr,nlv,nvar
	integer date,time
	character*80 title
	character*80 femver

	iunit = pentry(id)%iunit
	nvers = maxvers

	call shy_get_ftype(id,ftype)
	call shy_get_params(id,nkn,nel,npr,nlv,nvar)
	call shy_get_date(id,date,time)
	call shy_get_title(id,title)
	call shy_get_femver(id,femver)

        write(iunit,err=99) shytype,nvers
        write(iunit,err=99) ftype
        write(iunit,err=99) nkn,nel,npr,nlv,nvar
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
	integer i,k,ie,l,j
	integer, allocatable :: il(:)

	iunit = pentry(id)%iunit

	write(iunit,iostat=ierr) dtime,ivar,n,m,lmax
	if( ierr /= 0 ) return

	allocate(il(n))
	if( n == pentry(id)%nkn ) then
	  il = pentry(id)%ilhkv
	else if( n == pentry(id)%nel ) then
	  il = pentry(id)%ilhv
	else
	  write(6,*) n,pentry(id)%nkn,pentry(id)%nel
	  write(6,*) 'cannot determine layer pointer'
	  stop 'error stop shy_read_record: layer pointer'
	end if

	if( lmax <= 1 ) then
	  write(iunit,iostat=ierr) ( c(1,i),i=1,n*m )
	else if( m == 1 ) then
	  write(iunit,iostat=ierr) (( c(l,i)
     +			,l=1,il(i) )
     +			,i=1,n )
	else
	  write(iunit,iostat=ierr) (( c(l,i)
     +			,l=1,il(1+(i-1)/m) )
     +			,i=1,n*m )
	end if
	deallocate(il)

	end subroutine shy_write_record

!==================================================================
	end module shyfile
!==================================================================

!************************************************************
!************************************************************
!************************************************************

	subroutine test_shy
	end

!************************************************************

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

!************************************************************
!************************************************************
!************************************************************

	!program shy_main
	!call test_units
	!call test_shy
	!end

!************************************************************
!************************************************************
!************************************************************

