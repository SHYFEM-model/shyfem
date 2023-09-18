!
! $Id: subnos.f,v 1.15 2007-03-20 13:14:53 georg Exp $
!
! utility routines to read/write NOS file - file type 161
!
! contents :
!
!	 subroutine ininos
!	 subroutine setnos(iunit,nvers,nkn,nel,nlv,nvar)
!	 subroutine getnos(iunit,nvers,nkn,nel,nlv,nvar)
!	 subroutine delnos(iunit)
!	 subroutine dimnos(iunit,nknddi,nelddi,nlvddi)
!
!        subroutine errnos(iunit,routine,text)
!        subroutine findnos_err(iunit,routine,text,n)
!        function findnos(iunit)
!        subroutine infonos(iunit,iout)
!
!        subroutine nos_init(iunit,nversion)
!        subroutine nos_close(iunit)
!        subroutine nos_check_dimension(iunit,nknddi,nelddi,nlvddi)
!
!        subroutine nos_get_date(iunit,date,time)
!        subroutine nos_set_date(iunit,date,time)
!        subroutine nos_get_title(iunit,title)
!        subroutine nos_set_title(iunit,title)
!        subroutine nos_get_femver(iunit,femver)
!        subroutine nos_set_femver(iunit,femver)
!        subroutine nos_get_params(iunit,nkn,nel,nlv,nvar)
!        subroutine nos_set_params(iunit,nkn,nel,nlv,nvar)
!        subroutine nos_clone_params(iu_from,iu_to)
!        subroutine nos_check_compatibility(iu1,iu2)
!
!	 subroutine nos_is_nos_file(iunit,nvers)
!
!        subroutine nos_read_header(iunit,nkn,nel,nlv,nvar,ierr)
!        subroutine nos_write_header(iunit,nkn,nel,nlv,nvar,ierr)
!        subroutine nos_read_header2(iu,ilhkv,hlv,hev,ierr)
!        subroutine nos_write_header2(iunit,ilhkv,hlv,hev,ierr)
!        subroutine nos_read_record(iu,it,ivar,nlvddi,ilhkv,c,ierr)
!        subroutine nos_write_record(iunit,it,ivar,nlvddi,ilhkv,c,ierr)
!
!        subroutine nos_back_record(iunit)
!        subroutine nos_skip_header(iunit,nvar,ierr)
!        subroutine nos_skip_record(iunit,it,ivar,ierr)
!
! revision log :
!
! 19.08.1998	ggu	version 3 implemented
! 18.11.1998	ggu	subroutines setnos, getnos, dimnos added
! 26.01.1999	ggu	some comments
! 10.03.2004	ggu	in wfnos make title 80 chars long
! 28.09.2006	ggu	new routine delnos to close nos file explicitly
! 28.09.2006	ggu	check if nos file is already open
! 07.03.2007	ggu	new routines rhnos, whnos
! 10.02.2012	ggu	new routines to skip records
! 30.10.2012	ggu	new format for date and time (new accessor routines)
! 16.11.2012	ggu	in wrnos bugfix - call setnos first
! 02.12.2012	ggu	restructured
! 21.01.2013	ggu	code for next and back record
! 18.01.2014	ggu	restructured, new date,time,femver (version 4+5)
! 29.10.2014	ggu	new routine nos_is_nos_file()
! 22.09.2015	ggu	new routine nos_check_compatibility()
!
! notes :
!
! Usage writing:
!
!	open file
!	call nos_init
!	call nos_set_title	(not obligatory)
!	call nos_set_date	(not obligatory)
!	call nos_set_femver	(not obligatory)
!	call nos_write_header
!	call nos_write_header2
!	call nos_write_record
!	...
!	call nos_close
!
! Usage reading:
!
!	open file
!	call nos_init
!	call nos_read_header
!	call dimnos
!	call nos_get_title	(not obligatory)
!	call nos_get_date	(not obligatory)
!	call nos_get_femver	(not obligatory)
!	call nos_read_header2
!	call nos_read_record
!	...
!	call nos_close
!
! format of file:
!
! version 3 and greater
!
!	ftype,nvers
!	nkn,nel,nlv,nvar
!	title
!	date,time				(version 4)
!	femver					(version 5)
!
!	(ilhkv(k),k=1,nkn)			empty if nlv <= 1
!	(hlv(k),k=1,nlv)			empty if nlv <= 1
!	(hev(k),k=1,nel)
!	
!	it,ivar					(version <= 3)
!	it,ivar,lmax				(version >= 4)
!	(c(1,k),k=1,nkn)			if nlv <= 1 or lmax <= 1
!	((c(l,k),l=1,ilhkv(k)),k=1,nkn)		if nlv > 1
!
! version 2
!
!	ftype,nvers
!	nkn,nlv,nvar
!	title
!
!	(ilhkv(k),k=1,nkn)			only if nlv > 1
!	
!	it,ivar
!	((c(l,k),l=1,ilhkv(k)),k=1,nkn)		if nlv > 1
!	(c(1,k),k=1,nkn)			if nlv <= 1
!
! version 1
!
!	ftype,nvers
!	title
!	nkn,nvar
!
!	it
!	(c(i),i=1,nkn*nvar)
!
! variable types:
!
!	10	conz
!	11	salt
!	12	temp
!	15	oxygen
!	18	rms
!
! todo :
!
! - routine to close file (read and write) -> clnos()
! - remove data structure from nosvar when closed -> iunit=0
! - in setnos look for first available iunit
! - in setnos check if iunit is already open
!
!************************************************************
!************************************************************
!************************************************************
! internal routines
!************************************************************
!************************************************************
!************************************************************
!-------------------------------------------------------------
        module ionos
!-------------------------------------------------------------

        include 'nosinf.h'

!-------------------------------------------------------------
        contains
!-------------------------------------------------------------

	subroutine ininos

! sets up initial common block - internal routine

	implicit none

	integer i,n

	logical binit
	save binit
	data binit /.false./

	if( binit ) return

	binit = .true.

	nositem = 0

	do n=1,ndim
	  do i=0,nitdim
	    nosvar(i,n) = 0
	  end do
	  do i=1,nchdim
	    noschar(i,n) = ' '
	  end do
	end do

	end

!************************************************************

	subroutine setnos(iunit,nvers,nkn,nel,nlv,nvar)

! sets up parameter common block - internal routine

	implicit none

	integer iunit,nvers,nkn,nel,nlv,nvar

	integer n

! we do not check if unit has already been opened -> open with ifileo
! changed -> before calling this nos_init has to be called

	n = findnos(iunit)
	!if( n .eq. 0 ) then
	!  n = findnos(0)
	!end if

	if( n .eq. 0 ) then
	  call errnos(iunit,'setnos','Cannot find entry.')
	end if

	!nosvar(0,n) = iunit
	if( nvers .gt. 0 ) nosvar(1,n) = nvers
	if(   nkn .gt. 0 ) nosvar(2,n) = nkn
	if(   nel .gt. 0 ) nosvar(3,n) = nel
	if(   nlv .gt. 0 ) nosvar(4,n) = nlv
	if(  nvar .gt. 0 ) nosvar(5,n) = nvar

	end

!************************************************************

	subroutine getnos(iunit,nvers,nkn,nel,nlv,nvar)

! gets parameter common block - internal routine

	implicit none

	integer iunit,nvers,nkn,nel,nlv,nvar

	integer n

	n = findnos(iunit)
	if( n .eq. 0 ) then
	  call errnos(iunit,'getnos','File is not initialized.')
	end if

	nvers = nosvar(1,n)
	nkn   = nosvar(2,n)
	nel   = nosvar(3,n)
	nlv   = nosvar(4,n)
	nvar  = nosvar(5,n)

	end

!************************************************************

	subroutine delnos(iunit)

! closes nos file internal structure - internal routine
!
! please note that the file has still to be closed manually

	implicit none

	integer iunit

	integer n,i

	call findnos_err(iunit,'delnos','File is not open, cannot close.',n)

	do i=0,nitdim
	  nosvar(i,n) = 0
	end do
	do i=1,nchdim
	  noschar(i,n) = ' '
	end do

	end

!************************************************************

	subroutine dimnos(iunit,nknddi,nelddi,nlvddi)

! checks dimension of arrays - internal routine

	implicit none

	integer iunit,nknddi,nelddi,nlvddi

	integer nvers,nkn,nel,nlv,nvar

	call getnos(iunit,nvers,nkn,nel,nlv,nvar)

        if( nkn .gt. nknddi ) goto 99
        if( nel .gt. nelddi ) goto 99
        if( nlv .gt. nlvddi ) goto 99

	return
   99   continue
        write(6,*) 'nkn,nknddi : ',nkn,nknddi
        write(6,*) 'nel,nelddi : ',nel,nelddi
        write(6,*) 'nlv,nlvddi : ',nlv,nlvddi
        stop 'error stop dimnos: dimension error'
	end

!************************************************************
!************************************************************
!************************************************************

	subroutine errnos(iunit,routine,text)

! error routine for nos - internal routine

	implicit none

	integer iunit
	character*(*) routine,text

	write(6,*) 'For unit ',iunit,' in routine ',routine
	write(6,*) text
	stop 'error stop errnos'

	end

!************************************************************

	subroutine findnos_err(iunit,routine,text,n)

! finds entry for iunit -> returns it in n or stops with error

	implicit none

	integer iunit
	character*(*) routine,text
	integer n

	n = findnos(iunit)

	if( n .eq. 0 ) then
	  call errnos(iunit,routine,text)
	end if

	end

!************************************************************

	function findnos(iunit)

! finds entry for iunit - internal routine

	implicit none

	integer iunit
	integer findnos

	integer n

	do n=1,min(nositem+1,ndim)		!look at one entry more
	  if( nosvar(0,n) .eq. iunit ) goto 1
	end do
	n = 0
    1	continue

	if( n .gt. nositem ) nositem = n
	findnos = n

	end

!************************************************************

	subroutine infonos(iunit,iout)

! writes info for unit - internal routine

	implicit none

	integer iunit,iout

	integer n,i

	call findnos_err(iunit,'nos_info','Cannot find entry.',n)

	write(iout,*) 'iunit = ',iunit,' position = ',n
	
	do i=0,nitdim
	  write(iout,*) i,nosvar(i,n)
	end do

	end

!************************************************************
!************************************************************
!************************************************************
! public routines
!************************************************************
!************************************************************
!************************************************************

	subroutine nos_init(iunit,nversion)

	implicit none

	integer iunit
	integer nversion

	integer n,nvers

	call ininos

	if( iunit .le. 0 ) then
	  write(6,*) 'nos_init: Cannot initialize for this unit'
	  write(6,*) 'iunit = ',iunit
	  call errnos(iunit,'nos_init','Impossible unit number.')
	end if

	nvers = nversion
	if( nvers .le. 0 ) nvers = maxvers

	if( nvers .gt. maxvers ) then
	  write(6,*) 'nos_init: Impossible version number'
	  write(6,*) 'nvers = ',nvers,'   maxvers = ',maxvers
	  call errnos(iunit,'nos_init','Impossible version number.')
	end if

	if( nvers .lt. maxcomp ) then
	  write(6,*) 'nos_init: Old function call'
	  write(6,*) 'nvers = ',nvers,'   maxcomp = ',maxcomp
	  call errnos(iunit,'nos_init','Old function call.')
	end if

	nvers = maxvers	!always write with highest version

	n = findnos(iunit)
	if( n .ne. 0 ) then
	  call errnos(iunit,'nos_init','Unit already open.')
	end if

	n = findnos(0)
	if( n .eq. 0 ) then
	  call errnos(iunit,'nos_init','No space left (ndim).')
	end if

	nosvar(0,n) = iunit
	nosvar(1,n) = nvers

	rewind(iunit)

	end

!************************************************************

	subroutine nos_close(iunit)

	implicit none

	integer iunit

	call delnos(iunit)

	end

!************************************************************

	subroutine nos_check_dimension(iunit,nknddi,nelddi,nlvddi)

	implicit none

	integer iunit,nknddi,nelddi,nlvddi

	call dimnos(iunit,nknddi,nelddi,nlvddi)

	end

!************************************************************
!************************************************************
!************************************************************

	subroutine nos_get_date(iunit,date,time)

	implicit none

	integer iunit
	integer date,time

	integer n

	call findnos_err(iunit,'nos_get_date','Cannot find entry.',n)

	date = nosvar(6,n)
	time = nosvar(7,n)

	end

!************************************************************

	subroutine nos_set_date(iunit,date,time)

	implicit none

	integer iunit
	integer date,time

	integer n

	call findnos_err(iunit,'nos_set_date','Cannot find entry.',n)

	nosvar(6,n) = date
	nosvar(7,n) = time

	end

!************************************************************

	subroutine nos_get_title(iunit,title)

	implicit none

	integer iunit
	character*(*) title

	integer n

	call findnos_err(iunit,'nos_get_title','Cannot find entry.',n)

	title = noschar(1,n)

	end

!************************************************************

	subroutine nos_set_title(iunit,title)

	implicit none

	integer iunit
	character*(*) title

	integer n

	call findnos_err(iunit,'nos_set_title','Cannot find entry.',n)

	noschar(1,n) = title

	end

!************************************************************

	subroutine nos_get_femver(iunit,femver)

	implicit none

	integer iunit
	character*(*) femver

	integer n

	call findnos_err(iunit,'nos_get_femver','Cannot find entry.',n)

	femver = noschar(2,n)

	end

!************************************************************

	subroutine nos_set_femver(iunit,femver)

	implicit none

	integer iunit
	character*(*) femver

	integer n

	call findnos_err(iunit,'nos_set_femver','Cannot find entry.',n)

	noschar(2,n) = femver

	end

!************************************************************

	subroutine nos_get_params(iunit,nkn,nel,nlv,nvar)

	implicit none

	integer iunit
	integer nkn,nel,nlv,nvar

	integer nvers

	call getnos(iunit,nvers,nkn,nel,nlv,nvar)

	end

!************************************************************

	subroutine nos_set_params(iunit,nkn,nel,nlv,nvar)

	implicit none

	integer iunit
	integer nkn,nel,nlv,nvar

	call setnos(iunit,0,nkn,nel,nlv,nvar)

	end

!************************************************************

	subroutine nos_clone_params(iu_from,iu_to)

! clones data from one to other file 
!
! second file must have already been opened and initialized with nos_init
! should be only used to write file -> nvers should be max version

	implicit none

	integer iu_from
	integer iu_to

	integer i,nf,nt

	call findnos_err(iu_from,'nos_clone_params','Cannot find entry.',nf)
	call findnos_err(iu_to,'nos_clone_params','Cannot find entry.',nt)

	do i=2,nitdim		!unit and version are not cloned
	  nosvar(i,nt) = nosvar(i,nf)
	end do
	do i=1,nchdim
	  noschar(i,nt) = noschar(i,nf)
	end do

	end

!************************************************************

	subroutine nos_check_compatibility(iu1,iu2)

! checks compatibility between two nos files
!
! second file must have already been opened and initialized with nos_init
! should be only used to write file -> nvers should be max version

	implicit none

	integer iu1,iu2

	integer i,n1,n2

	call findnos_err(iu1,'nos_check_compatibility','Cannot find entry.',n1)
	call findnos_err(iu2,'nos_check_compatibility','Cannot find entry.',n2)

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

!************************************************************
!************************************************************
!************************************************************

	subroutine nos_is_nos_file(iunit,nvers)

! checks if iunit is open on nos file - returns nvers
!
! nvers == 0	no nos file (ntype is different) or read error
! nvers < 0	version number is wrong
! nvers > 0	good nos file

	implicit none

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

!************************************************************
!************************************************************
!************************************************************

	subroutine nos_read_header(iunit,nkn,nel,nlv,nvar,ierr)

! before this nos_init has to be called

	use shympi
        use timing

	implicit none

	integer iunit
	integer nkn,nel,nlv,nvar
	integer ierr

	integer n,nvers
	integer ntype,irec
	integer date,time
	character*80 line
        double precision time1

        if(ln_timing) time1 = shympi_wtime()

	call ininos

	call findnos_err(iunit,'nos_read_header','Cannot find entry.',n)

! first record - find out what version

	irec = 1
	read(iunit,end=91,err=99) ntype,nvers

! control version number and type of file

	if( ntype .ne. ftype ) goto 97
	if( nvers .le. 0 .or. nvers .gt. maxvers ) goto 98

! next records

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

	call setnos(iunit,nvers,nkn,nel,nlv,nvar)
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

        if(ln_timing) then
           io_time = io_time + shympi_wtime() - time1
           io_scalar_time = io_scalar_time + shympi_wtime() - time1
        endif

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

!********************************************************************

	subroutine nos_write_header(iunit,nkn,nel,nlv,nvar,ierr)

! writes first header of NOS file

	use shympi
        use timing

	implicit none

	integer iunit
	integer nkn,nel,nlv,nvar
	integer ierr

	integer n,nvers
	integer date,time
	character*80 title,femver
        double precision time1

        if(ln_timing) time1 = shympi_wtime()

	call ininos

	call findnos_err(iunit,'nos_write_header','Cannot find entry.',n)

	nvers = maxvers
	call setnos(iunit,nvers,nkn,nel,nlv,nvar)

	ierr=0
	if(.not.shympi_is_master()) return

	call nos_get_title(iunit,title)
	call nos_get_date(iunit,date,time)
	call nos_get_femver(iunit,femver)

	write(iunit)		ftype,maxvers
	write(iunit)		nkn,nel,nlv,nvar
	write(iunit)		title
	write(iunit)		date,time
	write(iunit)		femver

        if(ln_timing) then
           io_time = io_time + shympi_wtime() - time1
           io_scalar_time = io_scalar_time + shympi_wtime() - time1
        endif

	ierr=0

	end

!************************************************************

	subroutine nos_read_header2(iu,ilhkv,hlv,hev,ierr)

! reads second record of NOS file

        use timing
        use shympi

	implicit none

	integer iu
	integer ilhkv(*)
	double precision hlv(*)
	double precision hev(*)
        real, dimension(:),allocatable :: hl,he
	integer ierr

	logical bdata
	integer iunit
	integer k,l,ie
	integer nvers,nkn,nel,nlv,nvar
        double precision time1

        if(ln_timing) time1 = shympi_wtime()

	bdata = iu .gt. 0	!with negative unit number skip arrays
	iunit = abs(iu)

	call getnos(iunit,nvers,nkn,nel,nlv,nvar)

#ifdef SINGLEP
        allocate(hl(nlv))
        allocate(he(nel))
#endif

	if( .not. bdata ) then	!do not read arrays
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

! read records

	if( nvers .eq. 1 ) then
!	  no second header for version 1
	else if( nvers .eq. 2 ) then
	  if( nlv .ne. 0 ) then
	    read(iunit,err=99) (ilhkv(k),k=1,nkn)
	  end if
	else if( nvers .ge. 3 ) then
	  read(iunit,err=99) (ilhkv(k),k=1,nkn)
#ifdef SINGLEP
	  read(iunit,err=99) (hl(l),l=1,nlv)
	  read(iunit,err=99) (he(ie),ie=1,nel)
          hlv(1:nlv)=dble(hl(1:nlv))
          hev(1:nkn)=dble(he(1:nkn))
#else
	  read(iunit,err=99) (hlv(l),l=1,nlv)
	  read(iunit,err=99) (hev(ie),ie=1,nel)
#endif
	else
	   stop 'error stop nos_read_header2: internal error (1)'
	end if

	ierr = 0

        if(ln_timing) then
           io_time = io_time + shympi_wtime() - time1
           io_scalar_time = io_scalar_time + shympi_wtime() - time1
        endif

	return
   99	continue
	write(6,*) 'nos_read_header2: Error encountered while'
	write(6,*) 'reading second part of NOS file header'
	ierr=99
	return
	end

!************************************************************

	subroutine nos_write_header2(iunit,ilhkv,hlv,hev,ierr)

! writes second record of NOS file

	use shympi
        use timing

	implicit none

	integer iunit
	integer ilhkv(*)
	double precision hlv(*)
	double precision hev(*)
        real,dimension(:),allocatable :: hl
        real,dimension(:),allocatable :: he
	integer ierr

	integer k,l,ie
	integer nvers,nkn,nel,nlv,nvar
        double precision time1

        if(ln_timing) time1 = shympi_wtime()

	ierr = 0
	if(.not.shympi_is_master()) return

	call getnos(iunit,nvers,nkn,nel,nlv,nvar)
#ifdef SINGLEP
        allocate(hl(nlv))
        allocate(he(nel))
        hl=real(hlv(1:nlv))
        he=real(hev(1:nkn))
#endif

! only one layer

	if( nlv .le. 1 ) then
	  nlv = 0
	  nkn = 0
	end if

! write records

	write(iunit) (ilhkv(k),k=1,nkn)
#ifdef SINGLEP
	write(iunit) (hl(l),l=1,nlv)
	write(iunit) (he(ie),ie=1,nel)
#else
	write(iunit) (hlv(l),l=1,nlv)
	write(iunit) (hev(ie),ie=1,nel)
#endif

	ierr = 0

        if(ln_timing) then
           io_time = io_time + shympi_wtime() - time1
           io_scalar_time = io_scalar_time + shympi_wtime() - time1
        endif

	end

!************************************************************

	subroutine nos_read_record(iu,it,ivar,nlvddi,ilhkv,c,ierr)

! reads data record of NOS file

	use shympi
        use timing

	implicit none

! arguments
	integer iu,it,ivar
	integer nlvddi
	integer ilhkv(*)
	double precision c(nlvddi,*)
	real,dimension(:,:),allocatable :: rc
	integer ierr
! local
	integer l,k,lmax
	integer nvers,nkn,nel,nlv,nvar
	integer iunit
	logical bdata
        double precision time1

        if(ln_timing) time1 = shympi_wtime()

	bdata = iu .gt. 0	!with negative unit number only time record
	iunit = abs(iu)

	call getnos(iunit,nvers,nkn,nel,nlv,nvar)

#ifdef SINGLEP
        allocate (rc(nlvddi,nkn))
#endif
	!it = -1
	!ivar = 0
	lmax = min(nlv,nlvddi)

	if( nvers .eq. 1 ) then
	   ivar = 1
	   read(iunit,end=88,err=98) it
	   if( bdata ) then
#ifdef SINGLEP
	     read(iunit,end=99,err=99) (rc(1,k),k=1,nkn)
             c(1,1:nkn)=dble(rc(1,1:nkn))
#else
	     read(iunit,end=99,err=99) (c(1,k),k=1,nkn)
#endif
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
#ifdef SINGLEP
               read(iunit,end=99,err=99) (rc(1,k),k=1,nkn)
               c(1,1:nkn)=dble(rc(1,1:nkn))
#else
               read(iunit,end=99,err=99) (c(1,k),k=1,nkn)
#endif
	     else
#ifdef SINGLEP
               read(iunit,end=99,err=99) ((rc(l,k),l=1,ilhkv(k)),k=1,nkn)
               c(1:nlvddi,1:nkn)=dble(rc(1:nlvddi,1:nkn))
#else
               read(iunit,end=99,err=99) ((c(l,k),l=1,ilhkv(k)),k=1,nkn)
#endif
	     end if
	   else
	     read(iunit,end=99,err=99)
	   end if
	else
	   write(6,*) 'nvers = ',nvers,'  iunit = ',iunit
	   stop 'error stop nos_read_record: internal error (1)'
	end if

	ierr=0

        if(ln_timing) then
           io_time = io_time + shympi_wtime() - time1
           io_scalar_time = io_scalar_time + shympi_wtime() - time1
        endif

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

!************************************************************

	subroutine nos_write_record(iunit,it,ivar,nlvddi,ilhkv,c,ierr)

! writes data record of NOS file

	use shympi
        use timing

	implicit none

! arguments
	integer iunit,it,ivar
	integer nlvddi
	integer ilhkv(*)
	double precision c(nlvddi,*)
        real, dimension(:,:),allocatable :: rc
	integer ierr
! local
	integer l,k,lmax
	integer nvers,nkn,nel,nlv,nvar
        double precision time1

        if(ln_timing) time1 = shympi_wtime()

	ierr=0
	if(.not.shympi_is_master()) return

	call getnos(iunit,nvers,nkn,nel,nlv,nvar)
#ifdef SINGLEP
        allocate(rc(nlvddi,nkn))
        rc=real(c(1:nlvddi,1:nkn))
#endif
	lmax = min(nlv,nlvddi)

	write(iunit) it,ivar,lmax

	if( lmax .le. 1 ) then
#ifdef SINGLEP
	  write(iunit) (rc(1,k),k=1,nkn)
#else
	  write(iunit) (c(1,k),k=1,nkn)
#endif
	else
#ifdef SINGLEP
	  write(iunit) ((rc(l,k),l=1,ilhkv(k)),k=1,nkn)
#else
	  write(iunit) ((c(l,k),l=1,ilhkv(k)),k=1,nkn)
#endif
	end if

	ierr=0

        if(ln_timing) then
           io_time = io_time + shympi_wtime() - time1
           io_scalar_time = io_scalar_time + shympi_wtime() - time1
        endif

	return
	end

!************************************************************

	subroutine nos_peek_record(iu,it,ivar,ierr)

! peeks into data record of NOS file

	use shympi
        use timing

	implicit none

! arguments
	integer iu,it,ivar
	integer ierr
! local
	integer l,k,lmax
	integer nvers,nkn,nel,nlv,nvar
	integer iunit,ios
        double precision time1

        if(ln_timing) time1 = shympi_wtime()

	iunit = abs(iu)

	call getnos(iunit,nvers,nkn,nel,nlv,nvar)

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

        if(ln_timing) then
           io_time = io_time + shympi_wtime() - time1
           io_scalar_time = io_scalar_time + shympi_wtime() - time1
        endif

	end

!************************************************************
!************************************************************
!************************************************************

	subroutine nos_back_record(iunit)

! skips back one data record (contains two reads)

	implicit none

	integer iunit

	backspace(iunit)
	backspace(iunit)

	end

!************************************************************

	subroutine nos_skip_header(iunit,nvar,ierr)

	implicit none

	integer iunit,nvar,ierr

	integer nkn,nel,nlv
	integer ilhkv(1)
	double precision hlv(1)
	double precision hev(1)

	call nos_read_header(iunit,nkn,nel,nlv,nvar,ierr)
	if( ierr .ne. 0 ) return
	call nos_read_header2(-iunit,ilhkv,hlv,hev,ierr)

	end

!************************************************************

	subroutine nos_skip_record(iunit,it,ivar,ierr)

	implicit none

	integer iunit,it,ivar,ierr

	integer nlvddi
	integer ilhkv(1)
	double precision c(1,1)

	nlvddi = 1
	call nos_read_record(-iunit,it,ivar,nlvddi,ilhkv,c,ierr)

	end

!************************************************************
!************************************************************
!************************************************************
! old routines
!************************************************************
!************************************************************
!************************************************************

	subroutine rhnos (iunit,nvers,nknddi,nelddi,nlvddi,nkn,nel,nlv,nvar,ilhkv,hlv,hev,title)

! reads all headers of NOS file
!
! nvers		on entry maximal version that can be read
!		-> must be an input, used to check the corectness
!		.. of the call parameters
!		on return actual version read

	implicit none

! arguments
	integer iunit,nvers
	integer nknddi,nelddi,nlvddi
	integer nkn,nel,nlv,nvar
	integer date,time
	integer ilhkv(*)
	double precision hlv(*)
	double precision hev(*)
	character*(*) title
! local
	integer ierr,l

        call rfnos(iunit,nvers,nkn,nel,nlv,nvar,title,ierr)
        if(ierr.ne.0) goto 99

	call nos_get_date(iunit,date,time)

        write(6,*) 'nvers    : ',nvers
        write(6,*) 'nkn,nel  : ',nkn,nel
        write(6,*) 'nlv,nvar : ',nlv,nvar
        write(6,*) 'date,time: ',date,time
        write(6,*) 'title    : ',title

        call dimnos(iunit,nknddi,nelddi,nlvddi)

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

!************************************************************

        subroutine whnos(iunit,nvers,nkn,nel,nlv,nvar,ilhkv,hlv,hev,title)

! writes all headers of NOS file
!
! nvers         on entry maximal version
!               -> must be an input, used to check the corectness
!               .. of the call parameters

        implicit none

! arguments
        integer iunit,nvers
        integer nkn,nel,nlv,nvar
	integer ilhkv(*)
	double precision hlv(*)
	double precision hev(*)
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

!************************************************************
!************************************************************
!************************************************************
! compatibility
!************************************************************
!************************************************************
!************************************************************

	subroutine rfnos(iunit,nvers,nkn,nel,nlv,nvar,title,ierr)

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

!************************************************************

	subroutine wfnos(iunit,nvers,nkn,nel,nlv,nvar,title,ierr)

	implicit none

	integer iunit,nvers
	integer nkn,nel,nlv,nvar
	character*(*) title
	integer ierr

	call nos_init(iunit,nvers)
	call nos_set_title(iunit,title)
	call nos_write_header(iunit,nkn,nel,nlv,nvar,ierr)

	end

!************************************************************

	subroutine rsnos(iunit,ilhkv,hlv,hev,ierr)

	implicit none

	integer iunit
	integer ilhkv(*)
	double precision hlv(*)
	double precision hev(*)
	integer ierr

	call nos_read_header2(iunit,ilhkv,hlv,hev,ierr)

	end

!************************************************************

	subroutine wsnos(iunit,ilhkv,hlv,hev,ierr)

	implicit none

	integer iunit
	integer ilhkv(*)
	double precision hlv(*)
	double precision hev(*)
	integer ierr

	call nos_write_header2(iunit,ilhkv,hlv,hev,ierr)

	end

!************************************************************

	subroutine rdnos(iunit,it,ivar,nlvddi,ilhkv,c,ierr)

	implicit none

	integer iunit,it,ivar
	integer nlvddi
	integer ilhkv(*)
	double precision c(nlvddi,*)
	integer ierr

	call nos_read_record(iunit,it,ivar,nlvddi,ilhkv,c,ierr)

	end

!************************************************************

	subroutine wrnos(iunit,it,ivar,nlvddi,ilhkv,c,ierr)

        use shympi
        use mpi_io_admin

	implicit none

	integer iunit,it,ivar
	integer nlvddi
	integer ilhkv(*)
	double precision c(nlvddi,*)
	integer ierr

        if(shympi_partition_on_elements()) then
          call rebuild_scalar(nlvddi,c)
          if(shympi_is_master()) then
	    call nos_write_record(iunit,it,ivar,nlvddi,ilhkv,outTempv,ierr)
          end if
        else
	  call nos_write_record(iunit,it,ivar,nlvddi,ilhkv,c,ierr)
        end if

	end

!************************************************************
!************************************************************
!************************************************************

	subroutine infnos(ivar,name)

! returns description of variable id

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

!************************************************************

!-------------------------------------------------------------
        end module ionos
!-------------------------------------------------------------

