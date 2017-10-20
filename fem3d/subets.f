c
c $Id: subets.f,v 1.15 2007-03-20 13:14:53 georg Exp $
c
c utility routines to read/write ETS file - file type 161
c
c contents :
c
c	 subroutine iniets
c	 subroutine setets(iunit,nvers,nkn,nlv,nvar)
c	 subroutine getets(iunit,nvers,nkn,nlv,nvar)
c	 subroutine delets(iunit)
c	 subroutine dimets(iunit,nknddi,nlvddi)
c
c        subroutine errets(iunit,routine,text)
c        subroutine findets_err(iunit,routine,text,n)
c        function findets(iunit)
c        subroutine infoets(iunit,iout)
c
c        subroutine ets_init(iunit,nversion)
c        subroutine ets_close(iunit)
c        subroutine ets_check_dimension(iunit,nknddi,nlvddi)
c
c        subroutine ets_get_date(iunit,date,time)
c        subroutine ets_set_date(iunit,date,time)
c        subroutine ets_get_title(iunit,title)
c        subroutine ets_set_title(iunit,title)
c        subroutine ets_get_femver(iunit,femver)
c        subroutine ets_set_femver(iunit,femver)
c        subroutine ets_get_params(iunit,nkn,nlv,nvar)
c        subroutine ets_set_params(iunit,nkn,nlv,nvar)
c        subroutine ets_clone_params(iu_from,iu_to)
c
c	 subroutine ets_is_ets_file(iunit,nvers)
c
c	 subroutine ets_peek_header(iunit,nkn,nlv,nvar,ierr)
c        subroutine ets_read_header(iunit,nkn,nlv,nvar,ierr)
c        subroutine ets_write_header(iunit,nkn,nlv,nvar,ierr)
c        subroutine ets_read_header2(iu,ilhkv,hlv,hkv
c     +					,nodes,xg,yg,desc,ierr)
c        subroutine ets_write_header2(iunit,ilhkv,hlv,hkv
c     +					,nodes,xg,yg,desc,ierr)
c        subroutine ets_read_record(iu,it,ivar,nlvddi,ilhkv,c,ierr)
c        subroutine ets_write_record(iunit,it,ivar,nlvddi,ilhkv,c,ierr)
c
c        subroutine ets_next_record(iunit,it,ivar,ierr)
c        subroutine ets_back_record(iunit)
c        subroutine ets_skip_header(iunit,nvar,ierr)
c        subroutine ets_skip_record(iunit,it,ivar,ierr)
c
c revision log :
c
c 24.01.2014	ggu	copied from subnos.f
c
c notes :
c
c Usage writing:
c
c	open file
c	call ets_init
c	call ets_set_title	(not obligatory)
c	call ets_set_date	(not obligatory)
c	call ets_set_femver	(not obligatory)
c	call ets_write_header
c	call ets_write_header2
c	call ets_write_record
c	...
c	call ets_close
c
c Usage reading:
c
c	open file
c	call ets_init
c	call ets_read_header
c	call dimets
c	call ets_get_title	(not obligatory)
c	call ets_get_date	(not obligatory)
c	call ets_get_femver	(not obligatory)
c	call ets_read_header2
c	call ets_read_record
c	...
c	call ets_close
c
c format of file:
c
c version 1
c
c	ftype,nvers
c	nkn,nlv,nvar
c	title
c	date,time				(version 4)
c	femver					(version 5)
c
c	(ilhkv(k),k=1,nkn)			empty if nlv <= 1
c	(hlv(k),k=1,nlv)			empty if nlv <= 1
c	(hkv(k),k=1,nkn)
c	
c	it,ivar,lmax
c	(c(1,k),k=1,nkn)			if nlv <= 1 or lmax <= 1
c	((c(l,k),l=1,ilhkv(k)),k=1,nkn)		if nlv > 1
c
c variable types:
c
c	 1	zeta
c	 6	x-vel
c	 7	y-vel
c	10	conz
c	11	salt
c	12	temp
c	15	oxygen
c	18	rms
c
c************************************************************
c************************************************************
c************************************************************
c internal routines
c************************************************************
c************************************************************
c************************************************************

	subroutine iniets

c sets up initial common block - internal routine

	implicit none

	include 'etsinf.h'

	integer i,n

	logical binit
	save binit
	data binit /.false./

	if( binit ) return

	binit = .true.

	etsitem = 0

	do n=1,ndim
	  do i=0,nitdim
	    etsvar(i,n) = 0
	  end do
	  do i=1,nchdim
	    etschar(i,n) = ' '
	  end do
	end do

	end

c************************************************************

	subroutine setets(iunit,nvers,nkn,nlv,nvar)

c sets up parameter common block - internal routine

	implicit none

	include 'etsinf.h'

	integer iunit,nvers,nkn,nlv,nvar

	integer n
	integer findets

c we do not check if unit has already been opened -> open with ifileo
c changed -> before calling this ets_init has to be called

	n = findets(iunit)

	if( n .eq. 0 ) then
	  call errets(iunit,'setets','Cannot find entry.')
	end if

	if( nvers .gt. 0 ) etsvar(1,n) = nvers
	if(   nkn .gt. 0 ) etsvar(2,n) = nkn
	if(   nlv .gt. 0 ) etsvar(4,n) = nlv
	if(  nvar .gt. 0 ) etsvar(5,n) = nvar

	end

c************************************************************

	subroutine getets(iunit,nvers,nkn,nlv,nvar)

c gets parameter common block - internal routine

	implicit none

	include 'etsinf.h'

	integer iunit,nvers,nkn,nlv,nvar

	integer n
	integer findets

	n = findets(iunit)
	if( n .eq. 0 ) then
	  call errets(iunit,'getets','File is not initialized.')
	end if

	nvers = etsvar(1,n)
	nkn   = etsvar(2,n)
	nlv   = etsvar(4,n)
	nvar  = etsvar(5,n)

	end

c************************************************************

	subroutine delets(iunit)

c closes ets file internal structure - internal routine
c
c please note that the file has still to be closed manually

	implicit none

	include 'etsinf.h'

	integer iunit

	integer n,i

	call findets_err(iunit,'delets'
     +			,'File is not open, cannot close.',n)

	do i=0,nitdim
	  etsvar(i,n) = 0
	end do
	do i=1,nchdim
	  etschar(i,n) = ' '
	end do

	end

c************************************************************

	subroutine dimets(iunit,nknddi,nlvddi)

c checks dimension of arrays - internal routine

	implicit none

	integer iunit,nknddi,nlvddi

	integer nvers,nkn,nlv,nvar

	call getets(iunit,nvers,nkn,nlv,nvar)

        if( nkn .gt. nknddi ) goto 99
        if( nlv .gt. nlvddi ) goto 99

	return
   99   continue
        write(6,*) 'nkn,nknddi : ',nkn,nknddi
        write(6,*) 'nlv,nlvddi : ',nlv,nlvddi
        stop 'error stop dimets: dimension error'
	end

c************************************************************
c************************************************************
c************************************************************

	subroutine errets(iunit,routine,text)

c error routine for ets - internal routine

	implicit none

	integer iunit
	character*(*) routine,text

	write(6,*) 'For unit ',iunit,' in routine ',routine
	write(6,*) text
	stop 'error stop errets'

	end

c************************************************************

	subroutine findets_err(iunit,routine,text,n)

c finds entry for iunit -> returns it in n or stops with error

	implicit none

	include 'etsinf.h'

	integer iunit
	character*(*) routine,text
	integer n

	integer findets

	n = findets(iunit)

	if( n .eq. 0 ) then
	  call errets(iunit,routine,text)
	end if

	end

c************************************************************

	function findets(iunit)

c finds entry for iunit - internal routine

	implicit none

	include 'etsinf.h'

	integer iunit
	integer findets

	integer n

	do n=1,min(etsitem+1,ndim)		!look at one entry more
	  if( etsvar(0,n) .eq. iunit ) goto 1
	end do
	n = 0
    1	continue

	if( n .gt. etsitem ) etsitem = n
	findets = n

	end

c************************************************************

	subroutine infoets(iunit,iout)

c writes info for unit - internal routine

	implicit none

	include 'etsinf.h'

	integer iunit,iout

	integer n,i

	call findets_err(iunit,'ets_info','Cannot find entry.',n)

	write(iout,*) 'iunit = ',iunit,' position = ',n
	
	do i=0,nitdim
	  write(iout,*) i,etsvar(i,n)
	end do

	end

c************************************************************
c************************************************************
c************************************************************
c public routines
c************************************************************
c************************************************************
c************************************************************

	subroutine ets_init(iunit,nversion)

	implicit none

	include 'etsinf.h'

	integer iunit
	integer nversion

	integer n,nvers
	integer findets

	call iniets

	if( iunit .le. 0 ) then
	  write(6,*) 'ets_init: Cannot initialize for this unit'
	  write(6,*) 'iunit = ',iunit
	  call errets(iunit,'ets_init','Impossible unit number.')
	end if

	nvers = nversion
	if( nvers .le. 0 ) nvers = maxvers

	if( nvers .gt. maxvers ) then
	  write(6,*) 'ets_init: Impossible version number'
	  write(6,*) 'nvers = ',nvers,'   maxvers = ',maxvers
	  call errets(iunit,'ets_init','Impossible version number.')
	end if

	if( nvers .lt. maxcomp ) then
	  write(6,*) 'ets_init: Old function call'
	  write(6,*) 'nvers = ',nvers,'   maxcomp = ',maxcomp
	  call errets(iunit,'ets_init','Old function call.')
	end if

	nvers = maxvers	!always write with highest version

	n = findets(iunit)
	if( n .ne. 0 ) then
	  call errets(iunit,'ets_init','Unit already open.')
	end if

	n = findets(0)
	if( n .eq. 0 ) then
	  call errets(iunit,'ets_init','No space left (ndim).')
	end if

	etsvar(0,n) = iunit
	etsvar(1,n) = nvers

	rewind(iunit)

	end

c************************************************************

	subroutine ets_close(iunit)

	implicit none

	integer iunit

	call delets(iunit)

	end

c************************************************************

	subroutine ets_check_dimension(iunit,nknddi,nlvddi)

	implicit none

	integer iunit,nknddi,nlvddi

	call dimets(iunit,nknddi,nlvddi)

	end

c************************************************************
c************************************************************
c************************************************************

	subroutine ets_get_date(iunit,date,time)

	implicit none

	include 'etsinf.h'

	integer iunit
	integer date,time

	integer n

	call findets_err(iunit,'ets_get_date','Cannot find entry.',n)

	date = etsvar(6,n)
	time = etsvar(7,n)

	end

c************************************************************

	subroutine ets_set_date(iunit,date,time)

	implicit none

	include 'etsinf.h'

	integer iunit
	integer date,time

	integer n

	call findets_err(iunit,'ets_set_date','Cannot find entry.',n)

	etsvar(6,n) = date
	etsvar(7,n) = time

	end

c************************************************************

	subroutine ets_get_title(iunit,title)

	implicit none

	include 'etsinf.h'

	integer iunit
	character*(*) title

	integer n

	call findets_err(iunit,'ets_get_title','Cannot find entry.',n)

	title = etschar(1,n)

	end

c************************************************************

	subroutine ets_set_title(iunit,title)

	implicit none

	include 'etsinf.h'

	integer iunit
	character*(*) title

	integer n

	call findets_err(iunit,'ets_set_title','Cannot find entry.',n)

	etschar(1,n) = title

	end

c************************************************************

	subroutine ets_get_femver(iunit,femver)

	implicit none

	include 'etsinf.h'

	integer iunit
	character*(*) femver

	integer n

	call findets_err(iunit,'ets_get_femver','Cannot find entry.',n)

	femver = etschar(2,n)

	end

c************************************************************

	subroutine ets_set_femver(iunit,femver)

	implicit none

	include 'etsinf.h'

	integer iunit
	character*(*) femver

	integer n

	call findets_err(iunit,'ets_set_femver','Cannot find entry.',n)

	etschar(2,n) = femver

	end

c************************************************************

	subroutine ets_get_params(iunit,nkn,nlv,nvar)

	implicit none

	include 'etsinf.h'

	integer iunit
	integer nkn,nlv,nvar

	integer nvers

	call getets(iunit,nvers,nkn,nlv,nvar)

	end

c************************************************************

	subroutine ets_set_params(iunit,nkn,nlv,nvar)

	implicit none

	include 'etsinf.h'

	integer iunit
	integer nkn,nlv,nvar

	call setets(iunit,0,nkn,nlv,nvar)

	end

c************************************************************

	subroutine ets_clone_params(iu_from,iu_to)

c clones data from one to other file 
c
c second file must have already been opened and initialized with ets_init
c should be only used to write file -> nvers should be max version

	implicit none

	include 'etsinf.h'

	integer iu_from
	integer iu_to

	integer i,nf,nt

	call findets_err(iu_from,'ets_clone_params'
     +				,'Cannot find entry.',nf)
	call findets_err(iu_to,'ets_clone_params'
     +				,'Cannot find entry.',nt)

	do i=2,nitdim		!unit and version are not cloned
	  etsvar(i,nt) = etsvar(i,nf)
	end do
	do i=1,nchdim
	  etschar(i,nt) = etschar(i,nf)
	end do

	end

c************************************************************
c************************************************************
c************************************************************

        subroutine ets_is_ets_file(iunit,nvers)

c checks if iunit is open on ets file - returns nvers
c
c nvers == 0    no nos file (ntype is different) or read error
c nvers < 0     version number is wrong
c nvers > 0     good nos file

        implicit none

        include 'etsinf.h'

        integer iunit,nvers

        integer ntype

        nvers = 0
        if( iunit .le. 0 ) return

        read(iunit,end=1,err=1) ntype,nvers

        if( ntype .ne. ftype ) nvers = 0
        if( nvers .le. 0 .or. nvers .gt. maxvers ) nvers = -abs(nvers)

    1   continue

	end

c************************************************************
c************************************************************
c************************************************************

	subroutine ets_peek_header(iunit,nkn,nlv,nvar,ierr)

	implicit none

	include 'etsinf.h'

	integer iunit
	integer nkn,nlv,nvar
	integer ierr

	integer n,nvers
	integer ntype,irec

	call iniets

c first record - find out what version

	irec = 1
	read(iunit,end=91,err=99) ntype,nvers

c control version number and type of file

	if( ntype .ne. ftype ) goto 97
	if( nvers .le. 0 .or. nvers .gt. maxvers ) goto 98

c next records

	irec = 2
	if( nvers .ge. 1 ) then
	  read(iunit,err=99)	 nkn,nlv,nvar
	else
	   stop 'error stop ets_peek_header: internal error (1)'
	end if

	rewind(iunit)
	rewind(iunit)

	ierr = 0

	return
   99	continue
	write(6,*) 'ets_peek_header: Error encountered while'
	write(6,*) 'reading record number ',irec
	write(6,*) 'of ETS file header'
	write(6,*) 'nvers = ',nvers
	ierr=99
	return
   98	continue
	write(6,*) 'ets_peek_header: Version not recognized : ',nvers
	ierr=98
	return
   97	continue
	write(6,*) 'ets_peek_header: Wrong type of file : ',ntype
	write(6,*) 'Expected ',ftype
	ierr=97
	return
   91	continue
	write(6,*) 'ets_peek_header: File is empty'
	backspace(iunit)
	ierr=91
	return
	end

c************************************************************

	subroutine ets_read_header(iunit,nkn,nlv,nvar,ierr)

c before this ets_init has to be called

	implicit none

	include 'etsinf.h'

	integer iunit
	integer nkn,nlv,nvar
	integer ierr

	integer n,nvers
	integer ntype,irec
	integer date,time
	character*80 line

	call iniets

	call findets_err(iunit,'ets_read_header','Cannot find entry.',n)

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
	if( nvers .ge. 1 ) then
	  read(iunit,err=99)	 nkn,nlv,nvar
	  read(iunit,err=99)	 line
	  call setets(iunit,nvers,nkn,nlv,nvar)
	  call ets_set_title(iunit,line)
	else
	   stop 'error stop ets_read_header: internal error (1)'
	end if

	irec = 3
	if( nvers .ge. 1 ) then
	  read(iunit,err=99)	 date,time
	  read(iunit,err=99)	 line
	  call ets_set_date(iunit,date,time)
	  call ets_set_femver(iunit,line)
	end if

	ierr=0

	return
   99	continue
	write(6,*) 'ets_read_header: Error encountered while'
	write(6,*) 'reading record number ',irec
	write(6,*) 'of ETS file header'
	write(6,*) 'nvers = ',nvers
	ierr=99
	return
   98	continue
	write(6,*) 'ets_read_header: Version not recognized : ',nvers
	ierr=98
	return
   97	continue
	write(6,*) 'ets_read_header: Wrong type of file : ',ntype
	write(6,*) 'Expected ',ftype
	ierr=97
	return
   91	continue
	write(6,*) 'ets_read_header: File is empty'
	ierr=91
	return
	end

c********************************************************************

	subroutine ets_write_header(iunit,nkn,nlv,nvar,ierr)

c writes first header of ETS file

	implicit none

	include 'etsinf.h'

	integer iunit
	integer nkn,nlv,nvar
	integer ierr

	integer n,nvers
	integer date,time
	character*80 title,femver

	call iniets

	call findets_err(iunit,'ets_write_header','Cannot find entry.',n)

	nvers = maxvers
	call setets(iunit,nvers,nkn,nlv,nvar)

	call ets_get_title(iunit,title)
	call ets_get_date(iunit,date,time)
	call ets_get_femver(iunit,femver)

	write(iunit)		ftype,maxvers
	write(iunit)		nkn,nlv,nvar
	write(iunit)		title
	write(iunit)		date,time
	write(iunit)		femver

	ierr=0

	end

c************************************************************

	subroutine ets_read_header2(iu,ilhkv,hlv,hkv
     +					,nodes,xg,yg,desc,ierr)

c reads second record of ETS file

	implicit none

	integer iu
	integer ilhkv(1)
	real hlv(1)
	real hkv(1)
	integer nodes(1)
	real xg(1),yg(1)
	character*80 desc(1)
	integer ierr

	logical bdata
	integer iunit
	integer k,l
	integer nvers,nkn,nlv,nvar
	integer nkni,irec

	bdata = iu .gt. 0	!with negative unit number skip arrays
	iunit = abs(iu)

	call getets(iunit,nvers,nkn,nlv,nvar)

	nkni = nkn
	if( .not .bdata ) then	!do not read arrays
	  nkn = 0
	  nlv = 0
	  nkni = 0
	else if( nlv .le. 1 ) then
	  do k=1,nkn
	    ilhkv(k) = 1
	  end do
	  hlv(1) = 10000.

	  nkni = 0
	  nlv = 0
	end if

c read records

	irec = 1
	if( nvers .ge. 1 ) then
	  read(iunit,err=99) (ilhkv(k),k=1,nkni)
	  read(iunit,err=99) (hlv(l),l=1,nlv)
	  read(iunit,err=99) (hkv(k),k=1,nkn)
	else
	   stop 'error stop ets_read_header2: internal error (1)'
	end if

	irec = 2
	if( nvers .ge. 1 ) then
	  read(iunit,err=99) (nodes(k),k=1,nkn)
	  read(iunit,err=99) (xg(k),k=1,nkn)
	  read(iunit,err=99) (yg(k),k=1,nkn)
	  read(iunit,err=99) (desc(k),k=1,nkn)
	end if

	ierr = 0

	return
   99	continue
	write(6,*) 'ets_read_header2: Error encountered while'
	write(6,*) 'reading record number ',irec
	write(6,*) 'of ETS file second header'
	write(6,*) 'nvers = ',nvers
	ierr=99
	return
	end

c************************************************************

	subroutine ets_write_header2(iunit,ilhkv,hlv,hkv
     +					,nodes,xg,yg,desc,ierr)

c writes second record of ETS file

	implicit none

	integer iunit
	integer ilhkv(*)
	real hlv(*)
	real hkv(*)
	integer nodes(*)
	real xg(*),yg(*)
	character*80 desc(*)
	integer ierr

	integer k,l
	integer nvers,nkn,nlv,nvar
	integer nkni

	call getets(iunit,nvers,nkn,nlv,nvar)

c only one layer

	nkni = nkn
	if( nlv .le. 1 ) then
	  nlv = 0
	  nkni = 0
	end if

c write records

	write(iunit) (ilhkv(k),k=1,nkni)
	write(iunit) (hlv(l),l=1,nlv)
	write(iunit) (hkv(k),k=1,nkn)
	write(iunit) (nodes(k),k=1,nkn)
	write(iunit) (xg(k),k=1,nkn)
	write(iunit) (yg(k),k=1,nkn)
	write(iunit) (desc(k),k=1,nkn)

	ierr = 0

	end

c************************************************************

	subroutine ets_read_record(iu,it,ivar,nlvddi,ilhkv,c,ierr)

c reads data record of ETS file

	implicit none

c arguments
	integer iu,it,ivar
	integer nlvddi
	integer ilhkv(1)
	real c(nlvddi,1)
	integer ierr
c local
	integer l,k,lmax
	integer nvers,nkn,nlv,nvar
	integer iunit
	logical bdata

	bdata = iu .gt. 0	!with negative unit number only time record
	iunit = abs(iu)

	call getets(iunit,nvers,nkn,nlv,nvar)

	lmax = min(nlv,nlvddi)

	if( nvers .ge. 1 ) then
	   ivar = 1
	   read(iunit,end=88,err=98) it,ivar,lmax
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
	   stop 'error stop ets_read_record: internal error (1)'
	end if

	ierr=0

	return
   88	continue
	backspace(iunit)
	ierr=-1
	return
   98	continue
	write(6,*) 'ets_read_record: Error while reading'
	write(6,*) 'time record of ETS file'
	ierr=98
	return
   99	continue
	write(6,*) 'ets_read_record: Error while reading'
	write(6,*) 'data record of ETS file'
	write(6,*) 'it = ',it,'  ivar = ',ivar
	ierr=99
	return
	end

c************************************************************

	subroutine ets_write_record(iunit,it,ivar,nlvddi,ilhkv,c,ierr)

c writes data record of ETS file

	implicit none

c arguments
	integer iunit,it,ivar
	integer nlvddi
	integer ilhkv(*)
	real c(nlvddi,*)
	integer ierr
c local
	integer l,k,lmax
	integer nvers,nkn,nlv,nvar

	call getets(iunit,nvers,nkn,nlv,nvar)

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
c************************************************************
c************************************************************

	subroutine ets_next_record(iunit,it,ivar,ierr)

c skips data record - only reads header of record

	implicit none

	integer iunit,it,ivar,ierr

	integer nlvddi
	integer ilhkv(1)
	real c(1,1)

	nlvddi = 1
	call ets_read_record(-iunit,it,ivar,nlvddi,ilhkv,c,ierr)

	end

c************************************************************

	subroutine ets_back_record(iunit)

c skips back one data record (contains two reads)

	implicit none

	integer iunit

	backspace(iunit)
	backspace(iunit)

	end

c************************************************************

	subroutine ets_skip_header(iunit,nvar,ierr)

	implicit none

	integer iunit,nvar,ierr

	integer nkn,nlv
	integer ilhkv(1)
	real hlv(1)
	real hkv(1)
	integer nodes(1)
	real xg(1),yg(1)
	character*80 desc(1)

	call ets_read_header(iunit,nkn,nlv,nvar,ierr)
	if( ierr .ne. 0 ) return
	call ets_read_header2(-iunit,ilhkv,hlv,hkv
     +					,nodes,xg,yg,desc,ierr)

	end

c************************************************************

	subroutine ets_skip_record(iunit,it,ivar,ierr)

	implicit none

	integer iunit,it,ivar,ierr

	call ets_next_record(iunit,it,ivar,ierr)

	end

c************************************************************
c************************************************************
c************************************************************

	subroutine infets(ivar,name)

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

