c
c $Id: subnos.f,v 1.15 2007-03-20 13:14:53 georg Exp $
c
c utility routines to read/write NOS file - file type 161
c
c contents :
c
c subroutine ininos
c subroutine setnos(iunit,nvers,nkn,nel,nlv,nvar)
c subroutine getnos(iunit,nvers,nkn,nel,nlv,nvar)
c
c subroutine dimnos(iunit,nkndim,neldim,nlvdim)
c
c subroutine rfnos(iunit,nvers,nkn,nel,nlv,nvar,title,ierr)
c subroutine wfnos(iunit,nvers,nkn,nel,nlv,nvar,title,ierr)
c subroutine rsnos(iunit,ilhkv,hlv,hev,ierr)
c subroutine wsnos(iunit,ilhkv,hlv,hev,ierr)
c subroutine rdnos(iunit,it,ivar,nlvdim,ilhkv,c,ierr)
c subroutine wrnos(iunit,it,ivar,nlvdim,ilhkv,c,ierr)
c
c revision log :
c
c 19.08.1998	ggu	version 3 implemented
c 18.11.1998	ggu	subroutines setnos, getnos, dimnos added
c 26.01.1999	ggu	some comments
c 10.03.2004	ggu	in wfnos make title 80 chars long
c 28.09.2006	ggu	new routine delnos to close nos file explicitly
c 28.09.2006	ggu	check if nos file is already open
c 07.03.2007	ggu	new routines rhnos, whnos
c 10.02.2012	ggu	new routines to skip records
c
c notes :
c
c format of file:
c
c version 3
c
c	mtype,nvers
c	nkn,nel,nlv,nvar
c	title
c
c	(ilhkv(k),k=1,nkn)			empty if nlv <= 1
c	(hlv(k),k=1,nlv)			empty if nlv <= 1
c	(hev(k),k=1,nel)
c	
c	it,ivar
c	((c(l,k),l=1,ilhkv(k)),k=1,nkn)		if nlv > 1
c	(c(1,k),k=1,nkn)			if nlv <= 1
c
c version 2
c
c	mtype,nvers
c	nkn,nlv,nvar
c	title
c
c	(ilhkv(k),k=1,nkn)			only if nlv > 1
c	
c	it,ivar
c	((c(l,k),l=1,ilhkv(k)),k=1,nkn)		if nlv > 1
c	(c(1,k),k=1,nkn)			if nlv <= 1
c
c version 1
c
c	mtype,nvers
c	title
c	nkn,nvar
c
c	it
c	(c(i),i=1,nkn*nvar)
c
c variable types:
c
c	10	conz
c	11	salt
c	12	temp
c	15	oxygen
c	18	rms
c
c todo :
c
c - routine to close file (read and write) -> clnos()
c - remove data structure from nosvar when closed -> iunit=0
c - in setnos look for first available iunit
c - in setnos check if iunit is already open
c
c************************************************************

	subroutine ininos

c sets up initial common block - internal routine

	implicit none

c parameters
	integer ftype,maxvers
	parameter(ftype=161,maxvers=3)
	integer ndim,nitdim
	parameter(ndim=30,nitdim=5)
c common
	integer mtype,maxver
	common /noscom/ mtype,maxver
	integer nositem,nosvar(0:nitdim,ndim)
	common /nosvar/nositem,nosvar
c local
	integer i,n
c save
	logical binit
	save binit
	save /noscom/
	save /nosvar/
c data
	data binit /.false./

	if( binit ) return

	binit = .true.

	mtype = ftype
	maxver = maxvers

	nositem = 0
	do n=1,ndim
	  do i=0,nitdim
	    nosvar(i,n) = 0
	  end do
	end do

	end

c************************************************************

	subroutine setnos(iunit,nvers,nkn,nel,nlv,nvar)

c sets up parameter common block - internal routine

	implicit none

c arguments
	integer iunit,nvers,nkn,nel,nlv,nvar
c parameters
	integer ndim,nitdim
	parameter(ndim=30,nitdim=5)
c common
	integer nositem,nosvar(0:nitdim,ndim)
	common /nosvar/nositem,nosvar
c local
	integer n

c we do not check if unit has already been opened -> open with ifileo

	do n=1,nositem
	  if( nosvar(0,n) .eq. 0 ) goto 1
	  if( nosvar(0,n) .eq. iunit ) goto 99
	end do
    1	continue
	if( n .gt. nositem ) nositem = n

	if( nositem .gt. ndim ) then
	   stop 'error stop setnos: ndim'
	end if

	nosvar(0,n) = iunit
	nosvar(1,n) = nvers
	nosvar(2,n) = nkn
	nosvar(3,n) = nel
	nosvar(4,n) = nlv
	nosvar(5,n) = nvar

	return
   99	continue
	write(6,*) 'unit = ',iunit
	stop 'error stop setnos: unit already open - please close first'
	end

c************************************************************

	subroutine getnos(iunit,nvers,nkn,nel,nlv,nvar)

c gets parameter common block - internal routine

	implicit none

c arguments
	integer iunit,nvers,nkn,nel,nlv,nvar
c parameters
	integer ndim,nitdim
	parameter(ndim=30,nitdim=5)
c common
	integer nositem,nosvar(0:nitdim,ndim)
	common /nosvar/nositem,nosvar

	integer n

	do n=1,nositem
	  if( nosvar(0,n) .eq. iunit ) goto 1
	end do

	write(6,*) 'Cannot read on unit ',iunit
	write(6,*) 'File is not initialized.'
	stop 'error stop getnos: no initialization'
    1	continue

	nvers = nosvar(1,n)
	nkn   = nosvar(2,n)
	nel   = nosvar(3,n)
	nlv   = nosvar(4,n)
	nvar  = nosvar(5,n)

	end

c************************************************************

	subroutine delnos(iunit)

c closes nos file internal structure - internal routine
c
c please note that the file has still to be closed manually

	implicit none

c arguments
	integer iunit
c parameters
	integer ndim,nitdim
	parameter(ndim=30,nitdim=5)
c common
	integer nositem,nosvar(0:nitdim,ndim)
	common /nosvar/nositem,nosvar

	integer n

	do n=1,nositem
	  if( nosvar(0,n) .eq. iunit ) goto 1
	end do

	write(6,*) 'Cannot close unit ',iunit
	write(6,*) 'File is not open.'
	stop 'error stop delnos: file not open'
    1	continue

	nosvar(0,n) = 0

	end

c************************************************************

	subroutine dimnos(iunit,nkndim,neldim,nlvdim)

c checks dimension of arrays

	implicit none

c arguments
	integer iunit,nkndim,neldim,nlvdim
c local
	integer nvers,nkn,nel,nlv,nvar

	call getnos(iunit,nvers,nkn,nel,nlv,nvar)

        if( nkn .gt. nkndim ) goto 99
        if( nel .gt. neldim ) goto 99
        if( nlv .gt. nlvdim ) goto 99

	return
   99   continue
        write(6,*) 'nkn,nkndim : ',nkn,nkndim
        write(6,*) 'nel,neldim : ',nel,neldim
        write(6,*) 'nlv,nlvdim : ',nlv,nlvdim
        stop 'error stop dimnos: dimension error'
	end

c************************************************************
c************************************************************
c************************************************************
c************************************************************
c************************************************************

	subroutine rfnos	(iunit,nvers
     +				,nkn,nel,nlv,nvar
     +				,title
     +				,ierr
     +				)

c reads first record of NOS file
c
c nvers		on entry maximal version that can be read
c		-> must be an input, used to check the corectness
c		.. of the call parameters
c		on return actual version read

	implicit none

c arguments
	integer iunit,nvers
	integer nkn,nel,nlv,nvar
	character*(*) title
	integer ierr
c common
	integer mtype,maxver
	common /noscom/ mtype,maxver
c local
	integer ntype,irec
	character*80 line

c initialize

	call ininos

	line = ' '		!must be 80 chars

c control newest version number for call

	if( maxver .ne. nvers ) goto 95

c rewind file

	rewind(iunit,err=96)

c first record - find out what version

	irec = 1
	read(iunit,end=91,err=99) ntype,nvers

c control version number and type of file

	if( ntype .ne. mtype ) goto 97
	if( nvers .le. 0 .or. nvers .gt. maxver ) goto 98

c next records

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
	   stop 'error stop rfnos: internal error (1)'
	end if

	call setnos(iunit,nvers,nkn,nel,nlv,nvar)
	title = line

	ierr=0

	return
   99	continue
	write(6,*) 'rfnos: Error encountered while'
	write(6,*) 'reading record number ',irec
	write(6,*) 'of NOS file header'
	write(6,*) 'nvers = ',nvers
	ierr=99
	return
   98	continue
	write(6,*) 'rfnos: Version not recognized : ',nvers
	ierr=98
	return
   97	continue
	write(6,*) 'rfnos: Wrong type of file : ',ntype
	write(6,*) 'Expected ',mtype
	ierr=97
	return
   96	continue
	write(6,*) 'rfnos: Cannot rewind file for unit : ',iunit
	ierr=96
	return
   95	continue
	write(6,*) 'rfnos: Old function call ',nvers
	write(6,*) 'nvers = ',nvers,'   maxver = ',maxver
	write(6,*) 'Please adjust call to rfnos and recompile'
	ierr=95
	return
   91	continue
	write(6,*) 'rfnos: File is empty'
	ierr=91
	return
	end

c********************************************************************

	subroutine wfnos	(iunit,nvers
     +				,nkn,nel,nlv,nvar
     +				,title
     +				,ierr
     +				)

c writes first record of NOS file
c
c nvers		on entry maximal version
c		-> must be an input, used to check the corectness
c		.. of the call parameters

	implicit none

c arguments
	integer iunit,nvers
	integer nkn,nel,nlv,nvar
	character*(*) title
	integer ierr
c common
	integer mtype,maxver
	common /noscom/ mtype,maxver

	character*80 line

        line = ' '
        line = title            !must be 80 chars

c initialize

	call ininos

c control newest version number for call

	if( nvers.ne.maxver ) goto 95

	rewind(iunit)

	write(iunit)		mtype,maxver
	write(iunit)		nkn,nel,nlv,nvar
	write(iunit)		line

	call setnos(iunit,nvers,nkn,nel,nlv,nvar)

	ierr=0

	return
   95	continue
	write(6,*) 'wfnos: Old function call'
	write(6,*) 'nvers = ',nvers,'   maxver = ',maxver
	write(6,*) 'Please adjust call to rfnos and recompile'
	ierr=95
	return
	end

c************************************************************

	subroutine rsnos(iunit,ilhkv,hlv,hev,ierr)

c reads second record of NOS file

	implicit none

c arguments
	integer iunit
	integer ilhkv(1)
	real hlv(1)
	real hev(1)
	integer ierr
c common
	integer mtype,maxver
	common /noscom/ mtype,maxver
c local
	integer k,l,ie
	integer nvers,nkn,nel,nlv,nvar

	call getnos(iunit,nvers,nkn,nel,nlv,nvar)

c only one layer

	if( nlv .le. 1 ) then
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
	   stop 'error stop rsnos: internal error (1)'
	end if

	ierr = 0

	return
   99	continue
	write(6,*) 'rsnos: Error encountered while'
	write(6,*) 'reading second part of NOS file header'
	ierr=99
	return
	end

c************************************************************

	subroutine wsnos(iunit,ilhkv,hlv,hev,ierr)

c writes second record of NOS file

	implicit none

c arguments
	integer iunit
	integer ilhkv(1)
	real hlv(1)
	real hev(1)
	integer ierr
c common
	integer mtype,maxver
	common /noscom/ mtype,maxver
c local
	integer k,l,ie
	integer nvers,nkn,nel,nlv,nvar

	call getnos(iunit,nvers,nkn,nel,nlv,nvar)

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

	return
	end

c************************************************************

	subroutine rdnos(iunit,it,ivar,nlvdim,ilhkv,c,ierr)

c reads data record of NOS file

	implicit none

c arguments
	integer iunit,it,ivar
	integer nlvdim
	integer ilhkv(1)
	real c(nlvdim,1)
	integer ierr
c common
	integer mtype,maxver
	common /noscom/ mtype,maxver
c local
	integer l,k
	integer nvers,nkn,nel,nlv,nvar

	call getnos(iunit,nvers,nkn,nel,nlv,nvar)

	if( nvers .eq. 1 ) then
	   ivar = 1
	   read(iunit,end=88,err=98) it
c	   read(iunit,end=99,err=99) ((c(l,k),l=1,nlv),k=1,nkn)
	   read(iunit,end=99,err=99) (c(1,k),k=1,nkn)
	else if( nvers .ge. 2 ) then
	   read(iunit,end=88,err=98) it,ivar
	   if( nlv .le. 1 ) then
             read(iunit,end=99,err=99) (c(1,k),k=1,nkn)
	   else
             read(iunit,end=99,err=99) ((c(l,k),l=1,ilhkv(k)),k=1,nkn)
	   end if
	else
	   stop 'error stop rdnos: internal error (1)'
	end if

	ierr=0

	return
   88	continue
	ierr=-1
	return
   98	continue
	write(6,*) 'rdnos: Error while reading'
	write(6,*) 'time record of NOS file'
	ierr=98
	return
   99	continue
	write(6,*) 'rdnos: Error while reading'
	write(6,*) 'data record of NOS file'
	write(6,*) 'it = ',it,'  ivar = ',ivar
	ierr=99
	return
	end

c************************************************************

	subroutine wrnos(iunit,it,ivar,nlvdim,ilhkv,c,ierr)

c writes data record of NOS file

	implicit none

c arguments
	integer iunit,it,ivar
	integer nlvdim
	integer ilhkv(1)
	real c(nlvdim,1)
	integer ierr
c common
	integer mtype,maxver
	common /noscom/ mtype,maxver
c local
	integer l,k
	integer nvers,nkn,nel,nlv,nvar

	call getnos(iunit,nvers,nkn,nel,nlv,nvar)

	write(iunit) it,ivar

	if( nlv .le. 1 ) then
	  write(iunit) (c(1,k),k=1,nkn)
	else
	  write(iunit) ((c(l,k),l=1,ilhkv(k)),k=1,nkn)
	end if

	ierr=0

	return
	end

c************************************************************

	subroutine rhnos	(iunit,nvers
     +				,nkndim,neldim,nlvdim
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
	integer nkndim,neldim,nlvdim
	integer nkn,nel,nlv,nvar
	integer ilhkv(1)
	real hlv(1)
	real hev(1)
	character*(*) title
c local
	integer ierr,l

        call rfnos(iunit,nvers,nkn,nel,nlv,nvar,title,ierr)
        if(ierr.ne.0) goto 99

        write(6,*) 'nvers    : ',nvers
        write(6,*) 'nkn,nel  : ',nkn,nel
        write(6,*) 'nlv,nvar : ',nlv,nvar
        write(6,*) 'title    : ',title

        call dimnos(iunit,nkndim,neldim,nlvdim)

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

	subroutine shnos	(iunit,nvers
     +				,title
     +				)

c skips all headers of NOS file
c
c nvers		on entry maximal version that can be read
c		-> must be an input, used to check the corectness
c		.. of the call parameters
c		on return actual version read

	implicit none

	integer iunit,nvers
        character*(*) title

	integer ierr

	call sfnos(iunit,nvers,title,ierr)
        if(ierr.ne.0) goto 99

        call ssnos(iunit,ierr)
        if(ierr.ne.0) goto 98

	return
   98	continue
	write(6,*) 'ierr: ',ierr
	stop 'error stop shnos: error reading second header'
   99	continue
	write(6,*) 'ierr: ',ierr
	stop 'error stop shnos: error reading first header'
	end

c************************************************************

	subroutine sfnos(iunit,nvers,title,ierr)

c reads first header of nos file and skips information

	implicit none

	integer iunit,nvers
        character*(*) title
	integer ierr

	integer nkn,nel,nlv,nvar

	call rfnos	(iunit,nvers
     +				,nkn,nel,nlv,nvar
     +				,title
     +				,ierr
     +				)

	call delnos(iunit)			!trick false parameters
	call setnos(iunit,nvers,0,0,0,0)

	end

c************************************************************

	subroutine ssnos(iunit,ierr)

c reads second header of nos file and skips information

	implicit none

	integer iunit
	integer ierr

	integer ilhkv(1)
	real hlv(1),hev(1)

	call rsnos(iunit,ilhkv,hlv,hev,ierr)

	end

c************************************************************

	subroutine sknos(iunit,it,ivar,ierr)

c reads data of nos file and skips information

	implicit none

	integer iunit
	integer it,ivar
	integer ierr

	integer nlvdim
	integer ilhkv(1)
	real c(1,1)

	call rdnos(iunit,it,ivar,nlvdim,ilhkv,c,ierr)

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

