c
c $Id: subfil.f,v 1.7 1998/07/13 10:04:40 georg Exp $
c
c file opening routines
c
c contents :
c
c function ifileq(iunit,quest,format,status)	asks for filename and opens file
c function ifileo(iunit,file,format,status)	opens file
c subroutine filna(iunit,name)			inquires file name of iunit
c function filex(name)				inquires if file name exists
c function ifilun(name)				inquires unit of file name
c subroutine filext(file,ext)			appends file extension
c
c revision log :
c
c 02.05.1998	ggu	new format binary for ifileo
c 07.05.1998	ggu	new subroutine filext
c 21.05.1998	ggu	unit 0 is ok in ifileo()
c 29.06.1998	ggu	filna revisited
c
c********************************************************
c
	function ifileq(iunit,quest,format,status)
c
c asks for filename and opens file
c
c iunit		unit number to be opened
c quest		text to be written to terminal
c format	f*ormatted or u*nformatted
c status	o*ld, n*ew, s*cratch  or u*nknown
c ifileq	effective unit number the file is opened
c		-1 if no file is opened, 0 if no name is given
c
	character*(*) quest,format,status
	character*80  name
c
	write(6,*) quest
c
	read(5,1000) name
c
	if(name.eq.' ') then
		ifileq=0
	else
		ifileq=ifileo(iunit,name,format,status)
	end if
c
	return
c
 1000	format(a)
c
	end
c
c***************************************************************
c
	function ifileo(iunit,file,format,status)
c
c opens file
c
c iunit		unit number to be opened
c file		file name to be opened
c format	f*ormatted, u*nformatted or b*inary
c status	o*ld, n*ew, s*cratch  or u*nknown
c ifileo	effective unit number the file is opened
c
c routine tries to open file at iunit. If this is not possible
c ...higher unit numbers are tested and the file is eventually
c ...opened at one of these units. The effective unit number
c ...is passed back in ifileo. -1 is returned if no file is opened.
c
	character*(*) file,format,status
c
	logical found,error,opened,ex
	character*1 cf,cs
	character*15 form,stat,access
c
	external ichanm
	integer ichanm
c
	ifileo=-1
	iu=iunit
	if( iu .le. 0 ) iu = 55
c
	cf=format(1:1)
	cs=status(1:1)
	call uplow(cf,'low')
	call uplow(cs,'low')
c
	access = 'sequential'
c
	if(cf.eq.'f') then
		form='formatted'
	else if(cf.eq.'u') then
		form='unformatted'
	else if(cf.eq.'b') then
		form='unformatted'
clahey#		access = 'transparent'
	else
		write(6,*) 'format keyword not recognized :',format
		return
	end if
c
	if(cs.eq.'o') then
		stat='old'
	else if(cs.eq.'n') then
c for VAX      change to	stat='new'
c for DOS/UNIX change to	stat='unknown'
c		stat='new'
		stat='unknown'
	else if(cs.eq.'s') then
		stat='scratch'
	else if(cs.eq.'u') then
		stat='unknown'
	else
		write(6,*) 'status keyword not recognized :',status
		return
	end if
c
	inquire(file=file,exist=ex)
	if(.not.ex.and.stat.eq.'old') then
		nfile=ichanm(file)
		if(nfile.le.0) nfile=1
		write(6,*) 'file does not exist : ',file(1:nfile)
		return
	end if
c
	found=.false.
	error=.false.
c
	do while(.not.found.and..not.error)
		inquire(iu,exist=ex)
c		error=.not.ex
		error=.false.	!?? for lahey
		if(error) then
			write(6,*) 'no unit available to open file'
			return
		else
			inquire(iu,opened=opened)
			found=.not.opened
		end if
		if(.not.found) iu=iu+1
	end do
c
	if(found) then
          open(		 unit=iu
     +			,file=file
     +			,form=form
     +			,status=stat
     +			,access=access
     +			,iostat=ios
     +	      )
		if(ios.ne.0) then
			nfile=ichanm(file)
			if(nfile.le.0) nfile=1
			write(6,*) 'error opening file : '
     +				,file(1:nfile)
			write(6,*) 'unit : ',iu,'  iostat : ',ios
			write(6,*) 'error : ',mod(ios,256)
		else
			rewind(iu)
			ifileo=iu
		end if
	end if
c
	return
	end
c
c*******************************************************

	subroutine filna(iunit,name)

c inquires file name of unit number iunit
c
c iunit		unit number
c name		file name (return value)

	implicit none

	integer iunit
	character*(*) name

	logical btest

	name = ' '

	inquire(iunit,opened=btest)
	if( .not. btest ) return

	inquire(iunit,named=btest)
	if( .not. btest ) return

	inquire(unit=iunit,name=name)

	end

c*******************************************************
c
	function filex(name)
c
c inquires if file name exists
c
c name		file name
c filex		.true. if file exists
c
	logical filex
	character*(*) name
c
	inquire(file=name,exist=filex)
c
	return
	end
c
c*******************************************************
c
	function ifilun(name)
c
c inquires unit of file name if opened
c
c name		file name
c ifilun	unit if file is opened, else 0
c
	integer ifilun
	character*(*) name
	logical open
	integer iunit
c
	inquire(file=name,opened=open)
	if(open) then
		inquire(file=name,number=iunit)
	else
		iunit=0
	end if
c
	ifilun=iunit
c
	return
	end

c*******************************************************

	subroutine filext(file,ext)

c appends file extension if not there

	implicit none

	character*(*) file,ext

	integer nf,ne,length
	integer nef,nel,nff,nfl
	integer ichanm

	length = len(file)

	nf = ichanm(file)
	ne = ichanm(ext)

	nef = 1
	nel = ne
	if( ext(1:1) .eq. '.' ) then	!in case dot is in extension
	  nef = 2
	  ne = ne - 1
	end if

	if( ne .le. 0 ) return
	if( nf+ne+1 .gt. length ) return

c	now we know that there is enough room for extension

	nff = nf - ne + 1
	nfl = nf

	if( nff .le. 0 .or. file(nff:nfl) .ne. ext(nef:nel) ) then
	  file(nfl+1:) = '.' // ext(nef:nel)
	end if

	end

c*******************************************************

