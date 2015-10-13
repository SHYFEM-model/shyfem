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

	type, private :: entry

	  integer :: iunit
	  integer :: nvers
	  integer :: itype
	  double precision, allocatable :: array(:)

	end type entry

	integer, save, private :: idlast = 0
	integer, save, private :: ndim = 0
	type(entry), save, allocatable :: pentry(:)

	integer, parameter, private :: minvers = 11
	integer, parameter, private :: maxvers = 11

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

	subroutine shy_init_id(id)

	integer id

	if( id > ndim ) then
	  stop 'error stop shy_init_id: ndim'
	end if

	pentry(id)%iunit = 0
	pentry(id)%nvers = 0
	pentry(id)%itype = 0
	
	end subroutine shy_init_id

!******************************************************************

	subroutine shy_dealloc_arrays(id)

	integer id

	if( allocated(pentry(id)%array) ) deallocate(pentry(id)%array)

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

	subroutine shy_init(iunit,nversion,id)

	integer iunit
	integer nversion
	integer id

	integer nvers
	integer idempty

	if( iunit .le. 0 ) then
	  write(6,*) 'Impossible unit: ',iunit
	  stop 'error stop shy_init: iunit'
	end if

	nvers = nversion
	if( nvers .le. 0 ) nvers = maxvers

	if( nvers < minvers .or. nvers > maxvers ) then
	  write(6,*) 'Impossible version number: ',nvers
	  stop 'error stop shy_init: nvers'
	end if

	nvers = maxvers	!always write with highest version

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

	end subroutine shy_close

c************************************************************


!==================================================================
	end module shyfile
!==================================================================


!**************************************************************
!**************************************************************
!**************************************************************


c
c        subroutine nos_init(iunit,nversion)
c        subroutine nos_close(iunit)
c        subroutine nos_check_dimension(iunit,nknddi,nelddi,nlvddi)
c
c        subroutine nos_get_date(iunit,date,time)
c        subroutine nos_set_date(iunit,date,time)
c        subroutine nos_get_title(iunit,title)
c        subroutine nos_set_title(iunit,title)
c        subroutine nos_get_femver(iunit,femver)
c        subroutine nos_set_femver(iunit,femver)
c        subroutine nos_get_params(iunit,nkn,nel,nlv,nvar)
c        subroutine nos_set_params(iunit,nkn,nel,nlv,nvar)

c        subroutine nos_clone_params(iu_from,iu_to)
c        subroutine nos_check_compatibility(iu1,iu2)
c
c	 subroutine nos_is_nos_file(iunit,nvers)
c
c        subroutine nos_read_header(iunit,nkn,nel,nlv,nvar,ierr)
c        subroutine nos_write_header(iunit,nkn,nel,nlv,nvar,ierr)
c        subroutine nos_read_header2(iu,ilhkv,hlv,hev,ierr)
c        subroutine nos_write_header2(iunit,ilhkv,hlv,hev,ierr)
c        subroutine nos_read_record(iu,it,ivar,nlvddi,ilhkv,c,ierr)
c        subroutine nos_write_record(iunit,it,ivar,nlvddi,ilhkv,c,ierr)
c
c        subroutine nos_back_record(iunit)
c        subroutine nos_skip_header(iunit,nvar,ierr)
c        subroutine nos_skip_record(iunit,it,ivar,ierr)
c

c************************************************************
c************************************************************
c************************************************************
c public routines
c************************************************************
c************************************************************
c************************************************************

	subroutine nos_init(iunit,nversion)

	implicit none

	include 'nosinf.h'

	integer iunit
	integer nversion

	integer n,nvers
	integer findnos

	!!call ininos

	if( iunit .le. 0 ) then
	  write(6,*) 'nos_init: Cannot initialize for this unit'
	  write(6,*) 'iunit = ',iunit
	  !call errnos(iunit,'nos_init','Impossible unit number.')
	end if

	nvers = nversion
	if( nvers .le. 0 ) nvers = maxvers

	if( nvers .gt. maxvers ) then
	  write(6,*) 'nos_init: Impossible version number'
	  write(6,*) 'nvers = ',nvers,'   maxvers = ',maxvers
	  !call errnos(iunit,'nos_init','Impossible version number.')
	end if

	if( nvers .lt. maxcomp ) then
	  write(6,*) 'nos_init: Old function call'
	  write(6,*) 'nvers = ',nvers,'   maxcomp = ',maxcomp
	  !call errnos(iunit,'nos_init','Old function call.')
	end if

	nvers = maxvers	!always write with highest version

	!n = findnos(iunit)
	if( n .ne. 0 ) then
	  !call errnos(iunit,'nos_init','Unit already open.')
	end if

	!n = findnos(0)
	if( n .eq. 0 ) then
	  !call errnos(iunit,'nos_init','No space left (ndim).')
	end if

	nosvar(0,n) = iunit
	nosvar(1,n) = nvers

	rewind(iunit)

	end

c************************************************************

	subroutine nos_close(iunit)

	implicit none

	integer iunit

	!call delnos(iunit)

	end

c************************************************************

	subroutine nos_check_dimension(iunit,nknddi,nelddi,nlvddi)

	implicit none

	integer iunit,nknddi,nelddi,nlvddi

	!call dimnos(iunit,nknddi,nelddi,nlvddi)

	end

c************************************************************
c************************************************************
c************************************************************

	subroutine nos_get_date(iunit,date,time)

	implicit none

	include 'nosinf.h'

	integer iunit
	integer date,time

	integer n

	!call findnos_err(iunit,'nos_get_date','Cannot find entry.',n)

	date = nosvar(6,n)
	time = nosvar(7,n)

	end

c************************************************************

	subroutine nos_set_date(iunit,date,time)

	implicit none

	include 'nosinf.h'

	integer iunit
	integer date,time

	integer n

	!call findnos_err(iunit,'nos_set_date','Cannot find entry.',n)

	nosvar(6,n) = date
	nosvar(7,n) = time

	end

c************************************************************

	subroutine nos_get_title(iunit,title)

	implicit none

	include 'nosinf.h'

	integer iunit
	character*(*) title

	integer n

	!call findnos_err(iunit,'nos_get_title','Cannot find entry.',n)

	title = noschar(1,n)

	end

c************************************************************

	subroutine nos_set_title(iunit,title)

	implicit none

	include 'nosinf.h'

	integer iunit
	character*(*) title

	integer n

	!call findnos_err(iunit,'nos_set_title','Cannot find entry.',n)

	noschar(1,n) = title

	end

c************************************************************

	subroutine nos_get_femver(iunit,femver)

	implicit none

	include 'nosinf.h'

	integer iunit
	character*(*) femver

	integer n

	!call findnos_err(iunit,'nos_get_femver','Cannot find entry.',n)

	femver = noschar(2,n)

	end

c************************************************************

	subroutine nos_set_femver(iunit,femver)

	implicit none

	include 'nosinf.h'

	integer iunit
	character*(*) femver

	integer n

	!call findnos_err(iunit,'nos_set_femver','Cannot find entry.',n)

	noschar(2,n) = femver

	end

c************************************************************

	subroutine nos_get_params(iunit,nkn,nel,nlv,nvar)

	implicit none

	include 'nosinf.h'

	integer iunit
	integer nkn,nel,nlv,nvar

	integer nvers

	!call getnos(iunit,nvers,nkn,nel,nlv,nvar)

	end

c************************************************************

	subroutine nos_set_params(iunit,nkn,nel,nlv,nvar)

	implicit none

	include 'nosinf.h'

	integer iunit
	integer nkn,nel,nlv,nvar

	!call setnos(iunit,0,nkn,nel,nlv,nvar)

	end

c************************************************************

	subroutine nos_clone_params(iu_from,iu_to)

c clones data from one to other file 
c
c second file must have already been opened and initialized with nos_init
c should be only used to write file -> nvers should be max version

	implicit none

	include 'nosinf.h'

	integer iu_from
	integer iu_to

	integer i,nf,nt

	!call findnos_err(iu_from,'nos_clone_params'
!     +				,'Cannot find entry.',nf)
	!call findnos_err(iu_to,'nos_clone_params'
!     +				,'Cannot find entry.',nt)

	do i=2,nitdim		!unit and version are not cloned
	  nosvar(i,nt) = nosvar(i,nf)
	end do
	do i=1,nchdim
	  noschar(i,nt) = noschar(i,nf)
	end do

	end

c************************************************************

	subroutine nos_check_compatibility(iu1,iu2)

c checks compatibility between two nos files
c
c second file must have already been opened and initialized with nos_init
c should be only used to write file -> nvers should be max version

	implicit none

	include 'nosinf.h'

	integer iu1,iu2

	integer i,n1,n2

	!call findnos_err(iu1,'nos_check_compatibility'
   !  +				,'Cannot find entry.',n1)
	!call findnos_err(iu2,'nos_check_compatibility'
!     +				,'Cannot find entry.',n2)

	!unit and version are not checked (0,1)
	!neither are date and time (6,7)

	do i=2,5
	  if( nosvar(i,n1) /= nosvar(i,n1) ) exit
	end do

	if( i > 5 ) return	!all ok

	write(6,*) 'the two nos files are not compatible'
	write(6,*) 'units: ',n1,n2
	do i=0,nitdim
	  write(6,*) i,nosvar(i,n1),nosvar(i,n2)
	end do
	stop 'error stop nos_check_compatibility: not compatible'

	end

c************************************************************
c************************************************************
c************************************************************

	subroutine nos_is_nos_file(iunit,nvers)

c checks if iunit is open on nos file - returns nvers
c
c nvers == 0	no nos file (ntype is different) or read error
c nvers < 0	version number is wrong
c nvers > 0	good nos file

	implicit none

	include 'nosinf.h'

	integer iunit,nvers

	integer ntype

	nvers = 0
	if( iunit .le. 0 ) return

	read(iunit,end=1,err=1) ntype,nvers

	if( ntype .ne. ftype ) nvers = 0
	if( nvers .le. 0 .or. nvers .gt. maxvers ) nvers = -abs(nvers)

	rewind(iunit)

    1	continue

	end

c************************************************************
c************************************************************
c************************************************************

	subroutine nos_read_header(iunit,nkn,nel,nlv,nvar,ierr)

c before this nos_init has to be called

	implicit none

	include 'nosinf.h'

	integer iunit
	integer nkn,nel,nlv,nvar
	integer ierr

	integer n,nvers
	integer ntype,irec
	integer date,time
	character*80 line

	!call ininos

	!call findnos_err(iunit,'nos_read_header','Cannot find entry.',n)

c first record - find out what version

	irec = 1
	read(iunit,end=91,err=99) ntype,nvers

c control version number and type of file

	if( ntype .ne. ftype ) goto 97
	if( nvers .le. 0 .or. nvers .gt. maxvers ) goto 98

c next records

	date = 0
	time = 0

	irec = 2
	if( nvers .eq. 1 ) then
	  read(iunit,err=99) line
	  read(iunit,err=99) nkn,nlv	!nvar is actually nlv
	  nel = 0
	  nvar = 1			!always only one variable
	else if( nvers .eq. 2 ) then
	  read(iunit,err=99)	 nkn,nlv,nvar
	  read(iunit,err=99)	 line
	  nel = 0
	else if( nvers .ge. 3 ) then
	  read(iunit,err=99)	 nkn,nel,nlv,nvar
	  read(iunit,err=99)	 line
	else
	   stop 'error stop nos_read_header: internal error (1)'
	end if

	!call setnos(iunit,nvers,nkn,nel,nlv,nvar)
	call nos_set_title(iunit,line)

	irec = 3
	if( nvers .ge. 4 ) then
	  read(iunit,err=99)	 date,time
	  call nos_set_date(iunit,date,time)
	end if

	irec = 4
	if( nvers .ge. 5 ) then
	  read(iunit,err=99)	 line
	  call nos_set_femver(iunit,line)
	end if

	ierr=0

	return
   99	continue
	write(6,*) 'nos_read_header: Error encountered while'
	write(6,*) 'reading record number ',irec
	write(6,*) 'of NOS file header'
	write(6,*) 'nvers = ',nvers
	ierr=99
	return
   98	continue
	write(6,*) 'nos_read_header: Version not recognized : ',nvers
	ierr=98
	return
   97	continue
	write(6,*) 'nos_read_header: Wrong type of file : ',ntype
	write(6,*) 'Expected ',ftype
	ierr=97
	return
   91	continue
	write(6,*) 'nos_read_header: File is empty'
	ierr=91
	return
	end

c********************************************************************

	subroutine nos_write_header(iunit,nkn,nel,nlv,nvar,ierr)

c writes first header of NOS file

	implicit none

	include 'nosinf.h'

	integer iunit
	integer nkn,nel,nlv,nvar
	integer ierr

	integer n,nvers
	integer date,time
	character*80 title,femver

	!call ininos

	!call findnos_err(iunit,'nos_write_header','Cannot find entry.',n)

	nvers = maxvers
	!call setnos(iunit,nvers,nkn,nel,nlv,nvar)

	call nos_get_title(iunit,title)
	call nos_get_date(iunit,date,time)
	call nos_get_femver(iunit,femver)

	write(iunit)		ftype,maxvers
	write(iunit)		nkn,nel,nlv,nvar
	write(iunit)		title
	write(iunit)		date,time
	write(iunit)		femver

	ierr=0

	end

c************************************************************

	subroutine nos_read_header2(iu,ilhkv,hlv,hev,ierr)

c reads second record of NOS file

	implicit none

	integer iu
	integer ilhkv(1)
	real hlv(1)
	real hev(1)
	integer ierr

	logical bdata
	integer iunit
	integer k,l,ie
	integer nvers,nkn,nel,nlv,nvar

	bdata = iu .gt. 0	!with negative unit number skip arrays
	iunit = abs(iu)

	!!call getnos(iunit,nvers,nkn,nel,nlv,nvar)

	if( .not .bdata ) then	!do not read arrays
	  nkn = 0
	  nel = 0
	  nlv = 0
	else if( nlv .le. 1 ) then
	  do k=1,nkn
	    ilhkv(k) = 1
	  end do
	  hlv(1) = 10000.

	  nkn = 0
	  nlv = 0
	end if

c read records

	if( nvers .eq. 1 ) then
c	  no second header for version 1
	else if( nvers .eq. 2 ) then
	  if( nlv .ne. 0 ) then
	    read(iunit,err=99) (ilhkv(k),k=1,nkn)
	  end if
	else if( nvers .ge. 3 ) then
	  read(iunit,err=99) (ilhkv(k),k=1,nkn)
	  read(iunit,err=99) (hlv(l),l=1,nlv)
	  read(iunit,err=99) (hev(ie),ie=1,nel)
	else
	   stop 'error stop nos_read_header2: internal error (1)'
	end if

	ierr = 0

	return
   99	continue
	write(6,*) 'nos_read_header2: Error encountered while'
	write(6,*) 'reading second part of NOS file header'
	ierr=99
	return
	end

c************************************************************

	subroutine nos_write_header2(iunit,ilhkv,hlv,hev,ierr)

c writes second record of NOS file

	implicit none

	integer iunit
	integer ilhkv(1)
	real hlv(1)
	real hev(1)
	integer ierr

	integer k,l,ie
	integer nvers,nkn,nel,nlv,nvar

	!call getnos(iunit,nvers,nkn,nel,nlv,nvar)

c only one layer

	if( nlv .le. 1 ) then
	  nlv = 0
	  nkn = 0
	end if

c write records

	write(iunit) (ilhkv(k),k=1,nkn)
	write(iunit) (hlv(l),l=1,nlv)
	write(iunit) (hev(ie),ie=1,nel)

	ierr = 0

	end

c************************************************************

	subroutine nos_read_record(iu,it,ivar,nlvddi,ilhkv,c,ierr)

c reads data record of NOS file

	implicit none

c arguments
	integer iu,it,ivar
	integer nlvddi
	integer ilhkv(1)
	real c(nlvddi,1)
	integer ierr
c local
	integer l,k,lmax
	integer nvers,nkn,nel,nlv,nvar
	integer iunit
	logical bdata

	bdata = iu .gt. 0	!with negative unit number only time record
	iunit = abs(iu)

	!call getnos(iunit,nvers,nkn,nel,nlv,nvar)

	!it = -1
	!ivar = 0
	lmax = min(nlv,nlvddi)

	if( nvers .eq. 1 ) then
	   ivar = 1
	   read(iunit,end=88,err=98) it
	   if( bdata ) then
	     read(iunit,end=99,err=99) (c(1,k),k=1,nkn)
	   else
	     read(iunit,end=99,err=99)
	   end if
	else if( nvers .ge. 2 ) then
	   if( nvers .ge. 4 ) then
	     read(iunit,end=88,err=98) it,ivar,lmax
	   else
	     read(iunit,end=88,err=98) it,ivar
	   end if
	   if( bdata ) then
	     if( lmax .le. 1 ) then
               read(iunit,end=99,err=99) (c(1,k),k=1,nkn)
	     else
               read(iunit,end=99,err=99) ((c(l,k),l=1,ilhkv(k)),k=1,nkn)
	     end if
	   else
	     read(iunit,end=99,err=99)
	   end if
	else
	   write(6,*) 'nvers = ',nvers,'  iunit = ',iunit
	   stop 'error stop nos_read_record: internal error (1)'
	end if

	ierr=0

	return
   88	continue
	ierr=-1
	return
   98	continue
	write(6,*) 'nvers,it,ivar,lmax: ',nvers,it,ivar,lmax
	write(6,*) 'nos_read_record: Error while reading'
	write(6,*) 'time record of NOS file'
	ierr=98
	return
   99	continue
	write(6,*) 'nos_read_record: Error while reading'
	write(6,*) 'data record of NOS file'
	write(6,*) 'it = ',it,'  ivar = ',ivar
	ierr=99
	return
	end

c************************************************************

	subroutine nos_write_record(iunit,it,ivar,nlvddi,ilhkv,c,ierr)

c writes data record of NOS file

	implicit none

c arguments
	integer iunit,it,ivar
	integer nlvddi
	integer ilhkv(1)
	real c(nlvddi,1)
	integer ierr
c local
	integer l,k,lmax
	integer nvers,nkn,nel,nlv,nvar

	!call getnos(iunit,nvers,nkn,nel,nlv,nvar)

	lmax = min(nlv,nlvddi)

	write(iunit) it,ivar,lmax

	if( lmax .le. 1 ) then
	  write(iunit) (c(1,k),k=1,nkn)
	else
	  write(iunit) ((c(l,k),l=1,ilhkv(k)),k=1,nkn)
	end if

	ierr=0

	return
	end

c************************************************************

	subroutine nos_peek_record(iu,it,ivar,ierr)

c peeks into data record of NOS file

	implicit none

c arguments
	integer iu,it,ivar
	integer ierr
c local
	integer l,k,lmax
	integer nvers,nkn,nel,nlv,nvar
	integer iunit,ios

	iunit = abs(iu)

	!call getnos(iunit,nvers,nkn,nel,nlv,nvar)

	lmax = nlv

	if( nvers .eq. 1 ) then
	   ivar = 1
	   read(iunit,iostat=ios) it
	else if( nvers .ge. 2 ) then
	   if( nvers .ge. 4 ) then
	     read(iunit,iostat=ios) it,ivar,lmax
	   else
	     read(iunit,iostat=ios) it,ivar
	   end if
	else
	   write(6,*) 'nvers = ',nvers,'  iunit = ',iunit
	   stop 'error stop nos_peek_record: internal error (1)'
	end if

	if( ios > 0 ) then
	  write(6,*) 'nos_peek_record: Error while reading'
	  write(6,*) 'time record of NOS file'
	  ierr=98
	  return
	end if

	backspace(iu)

	if( ios < 0 ) then
	  ierr=-1
	else
	  ierr=0
	end if

	end

c************************************************************
c************************************************************
c************************************************************

	subroutine nos_back_record(iunit)

c skips back one data record (contains two reads)

	implicit none

	integer iunit

	backspace(iunit)
	backspace(iunit)

	end

c************************************************************

	subroutine nos_skip_header(iunit,nvar,ierr)

	implicit none

	integer iunit,nvar,ierr

	integer nkn,nel,nlv
	integer ilhkv(1)
	real hlv(1)
	real hev(1)

	call nos_read_header(iunit,nkn,nel,nlv,nvar,ierr)
	if( ierr .ne. 0 ) return
	call nos_read_header2(-iunit,ilhkv,hlv,hev,ierr)

	end

c************************************************************

	subroutine nos_skip_record(iunit,it,ivar,ierr)

	implicit none

	integer iunit,it,ivar,ierr

	integer nlvddi
	integer ilhkv(1)
	real c(1,1)

	nlvddi = 1
	call nos_read_record(-iunit,it,ivar,nlvddi,ilhkv,c,ierr)

	end

c************************************************************
c************************************************************
c************************************************************
c old routines
c************************************************************
c************************************************************
c************************************************************

	subroutine rhnos	(iunit,nvers
     +				,nknddi,nelddi,nlvddi
     +				,nkn,nel,nlv,nvar
     +				,ilhkv,hlv,hev
     +				,title
     +				)

c reads all headers of NOS file
c
c nvers		on entry maximal version that can be read
c		-> must be an input, used to check the corectness
c		.. of the call parameters
c		on return actual version read

	implicit none

c arguments
	integer iunit,nvers
	integer nknddi,nelddi,nlvddi
	integer nkn,nel,nlv,nvar
	integer date,time
	integer ilhkv(1)
	real hlv(1)
	real hev(1)
	character*(*) title
c local
	integer ierr,l

        call rfnos(iunit,nvers,nkn,nel,nlv,nvar,title,ierr)
        if(ierr.ne.0) goto 99

	call nos_get_date(iunit,date,time)

        write(6,*) 'nvers    : ',nvers
        write(6,*) 'nkn,nel  : ',nkn,nel
        write(6,*) 'nlv,nvar : ',nlv,nvar
        write(6,*) 'date,time: ',date,time
        write(6,*) 'title    : ',title

        !call dimnos(iunit,nknddi,nelddi,nlvddi)

        call rsnos(iunit,ilhkv,hlv,hev,ierr)
        if(ierr.ne.0) goto 98

        write(6,*) 'Available levels: ',nlv
        write(6,*) (hlv(l),l=1,nlv)

	return
   98	continue
	write(6,*) 'ierr: ',ierr
	stop 'error stop rhnos: error reading second header'
   99	continue
	write(6,*) 'ierr: ',ierr
	stop 'error stop rhnos: error reading first header'
	end

c************************************************************

        subroutine whnos        (iunit,nvers
     +                          ,nkn,nel,nlv,nvar
     +				,ilhkv,hlv,hev
     +                          ,title
     +                          )

c writes all headers of NOS file
c
c nvers         on entry maximal version
c               -> must be an input, used to check the corectness
c               .. of the call parameters

        implicit none

c arguments
        integer iunit,nvers
        integer nkn,nel,nlv,nvar
	integer ilhkv(1)
	real hlv(1)
	real hev(1)
        character*(*) title

        integer ierr

        call wfnos(iunit,nvers,nkn,nel,nlv,nvar,title,ierr)
        if(ierr.ne.0) goto 1

	call wsnos(iunit,ilhkv,hlv,hev,ierr)
        if(ierr.ne.0) goto 1

	return
    1	continue
	write(6,*) 'ierr: ',ierr
	stop 'error stop whnos'
	end

c************************************************************
c************************************************************
c************************************************************
c compatibility
c************************************************************
c************************************************************
c************************************************************

	subroutine rfnos	(iunit,nvers
     +				,nkn,nel,nlv,nvar
     +				,title
     +				,ierr
     +				)

	implicit none

	integer iunit,nvers
	integer nkn,nel,nlv,nvar
	character*(*) title
	integer ierr

	write(6,*) 'rfnos: ',iunit,nvers
	call nos_init(iunit,nvers)
	call nos_read_header(iunit,nkn,nel,nlv,nvar,ierr)
	call nos_get_title(iunit,title)

	end

c************************************************************

	subroutine wfnos	(iunit,nvers
     +				,nkn,nel,nlv,nvar
     +				,title
     +				,ierr
     +				)

	implicit none

	integer iunit,nvers
	integer nkn,nel,nlv,nvar
	character*(*) title
	integer ierr

	call nos_init(iunit,nvers)
	call nos_set_title(iunit,title)
	call nos_write_header(iunit,nkn,nel,nlv,nvar,ierr)

	end

c************************************************************

	subroutine rsnos(iunit,ilhkv,hlv,hev,ierr)

	implicit none

	integer iunit
	integer ilhkv(1)
	real hlv(1)
	real hev(1)
	integer ierr

	call nos_read_header2(iunit,ilhkv,hlv,hev,ierr)

	end

c************************************************************

	subroutine wsnos(iunit,ilhkv,hlv,hev,ierr)

	implicit none

	integer iunit
	integer ilhkv(1)
	real hlv(1)
	real hev(1)
	integer ierr

	call nos_write_header2(iunit,ilhkv,hlv,hev,ierr)

	end

c************************************************************

	subroutine rdnos(iunit,it,ivar,nlvddi,ilhkv,c,ierr)

	implicit none

	integer iunit,it,ivar
	integer nlvddi
	integer ilhkv(1)
	real c(nlvddi,1)
	integer ierr

	call nos_read_record(iunit,it,ivar,nlvddi,ilhkv,c,ierr)

	end

c************************************************************

	subroutine wrnos(iunit,it,ivar,nlvddi,ilhkv,c,ierr)

	implicit none

	integer iunit,it,ivar
	integer nlvddi
	integer ilhkv(1)
	real c(nlvddi,1)
	integer ierr

	call nos_write_record(iunit,it,ivar,nlvddi,ilhkv,c,ierr)

	end

c************************************************************
c************************************************************
c************************************************************

	subroutine infnos(ivar,name)

c returns description of variable id

	implicit none

	integer ivar
	character*(*) name

	if( ivar .eq. 10 ) then
	  name = 'Concentration []'
	else if( ivar .eq. 11 ) then
	  name = 'Salinity [psu]'
	else if( ivar .eq. 12 ) then
	  name = 'Temperature [C]'
	else if( ivar .eq. 99 ) then
	  name = 'Water residence time [days]'
	else
	  name = ''
	end if

	end

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

