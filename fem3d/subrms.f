
!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2018  Georg Umgiesser
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

! revision log :
!
! 23.03.2010	ggu	changed v6.1.1
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62

!**************************************************************************

c file name routines (for VAX/VMS systems)
c
c contents :
c
c function namdef(...)			supplies default values for file name
c function iscnam(name,type,back)	scanns a file name
c
c***********************************************************
c
	function namdef(name,node,device,directory,file
     +			,extension,version,back)
c
c supplies default values for file name
c
c name		original file name
c node,...	defaults for different portions of file name
c back		complete file name (missing portions in 'name'
c		...have been substituted by their default values
c		...passed in node,...)
c namdef	length of file name in 'back'
c
	character*(*) name,node,device,directory,file
	character*(*) extension,version,back
	character*80 help1,help2
c
c test for error in name
c
	ichar1=iscnam(name,'node',help1)
c
	if(ichar1.lt.-1) goto 88	!error
c
	nextc=1
	nback=len(back)
c
c node
c
	ichar2=iscnam(node,'node',help2)
c
	if(ichar1.gt.0) then
		if(nextc.gt.nback) goto 90
		back(nextc:)=help1
		nextc=nextc+ichar1
	else if(ichar2.gt.0) then
		if(nextc.gt.nback) goto 90
		back(nextc:)=help2
		nextc=nextc+ichar2
	end if
c
c device
c
	ichar1=iscnam(name,'device',help1)
	ichar2=iscnam(device,'device',help2)
c
	if(ichar1.gt.0) then
		if(nextc.gt.nback) goto 90
		back(nextc:)=help1
		nextc=nextc+ichar1
	else if(ichar2.gt.0) then
		if(nextc.gt.nback) goto 90
		back(nextc:)=help2
		nextc=nextc+ichar2
	end if
c
c directory
c
	ichar1=iscnam(name,'directory',help1)
	ichar2=iscnam(directory,'directory',help2)
c
	if(ichar1.gt.0) then
		if(nextc.gt.nback) goto 90
		back(nextc:)=help1
		nextc=nextc+ichar1
	else if(ichar2.gt.0) then
		if(nextc.gt.nback) goto 90
		back(nextc:)=help2
		nextc=nextc+ichar2
	end if
c
c file
c
	ichar1=iscnam(name,'file',help1)
	ichar2=iscnam(file,'file',help2)
c
	if(ichar1.gt.0) then
		if(nextc.gt.nback) goto 90
		back(nextc:)=help1
		nextc=nextc+ichar1
	else if(ichar2.gt.0) then
		if(nextc.gt.nback) goto 90
		back(nextc:)=help2
		nextc=nextc+ichar2
	end if
c
c extension
c
	ichar1=iscnam(name,'extension',help1)
	ichar2=iscnam(extension,'extension',help2)
c
	if(ichar1.gt.0) then
		if(nextc.gt.nback) goto 90
		back(nextc:)=help1
		nextc=nextc+ichar1
	else if(ichar2.gt.0) then
		if(nextc.gt.nback) goto 90
		back(nextc:)=help2
		nextc=nextc+ichar2
	end if
c
c version
c
	ichar1=iscnam(name,'version',help1)
	ichar2=iscnam(version,'version',help2)
c
	if(ichar1.gt.0) then
		if(nextc.gt.nback) goto 90
		back(nextc:)=help1
		nextc=nextc+ichar1
	else if(ichar2.gt.0) then
		if(nextc.gt.nback) goto 90
		back(nextc:)=help2
		nextc=nextc+ichar2
	end if
c
	namdef=nextc-1
c
	goto 100
c
c error traps and end
c
   90	continue
	write(6,*) 'name too long. does not fit in back :'
	write(6,*) back
	namdef=-10
	goto 100
   88	continue
	namdef=ichar1
	goto 100
c
  100	continue
c
	return
	end
c
c****************************************************************
c
	function iscnam1(name,type,back)
c
c scanns a file name and gives back the desired component
c
c VMS file name format
c
c name		file name
c type		the desired portion of the file name
c		...type has to be one of these keywords :
c		... [ node  device  directory  file  extension  version ]
c		...non ambigous abreviations can be used
c back		the returned portion of the file name
c		...nothing is returned if the desired portion
c		...of the file is missing
c		...node and device are returned with their colons
c		...directory is returned with their brackets
c		...extension and version with their leading dot
c		...or semicolon
c iscnam1	length of the character string (return value)
c		...(0 if nothing is returned)
c		...if iscnam1 is negativ, an error has occured
c 		... error code :
c			-1: nothing in name
c			-2: error in file name
c			-3: ambigous type
c			-4: type not recognized
c			-5: dimension of back too small
c			-6: name too long (maximum 80 characters)
c
	character*(*) name,type,back
c
	character*80 node,device,directory,file,extension,version
	character*80 help
	character*1  c
	logical btest
c
	data btest /.false./	!.true. ==> test output on terminal
c
	nname=len(name)
	ntype=len(type)
	nback=len(back)
c
	node=' '
	device=' '
	directory=' '
	file=' '
	extension=' '
	version=' '
	back=' '
c	
	do i=1,nname
	if(name(i:i).ne.' ') goto 1
	end do
	goto 99
    1	ifirst=i
c
	i=ifirst
	istart=i
	indir=0
	lwhere=0
	c=name(i:i)
	do while(c.ne.' '.and.i.le.nname)
c
	if(c.eq.':') then
		if(i.lt.nname.and.name(i+1:i+1).eq.':') then	!node
			if(lwhere.gt.0) goto 98	!node must be first
			if(i.eq.istart) goto 98	!no node name
			if(i+1-istart+1.gt.80) goto 94 !name too long
			node=name(istart:i+1)
			istart=i+2
			i=i+1
			lwhere=1
		else				!device
			if(lwhere.gt.1) goto 98	!device must be second
			if(i.eq.istart) goto 98	!no device name
			if(i-istart+1.gt.80) goto 94 !name too long
			device=name(istart:i)
			istart=i+1
			lwhere=2
		end if
	else if(c.eq.'[') then			!start directory
		if(lwhere.gt.2) goto 98	!directory must be third
		if(i.ne.istart) goto 98 ![ must start a directory
		if(indir.eq.1) goto 98	!already in directory
		indir=1
	else if(c.eq.']') then			!end directory
		if(lwhere.gt.2) goto 98	!directory must be third
		if(i.eq.istart) goto 98	!] must end a directory
		if(indir.eq.0) goto 98	!not in directory
		if(i-istart+1.gt.80) goto 94	!name too long
		directory=name(istart:i)
		indir=0
		istart=i+1
		lwhere=3
	else if(c.eq.'.') then		!extension or version
		if(indir.eq.0) then	!no subdirectory
			if(lwhere.le.3) then		   !file
				if(i-istart.gt.80) goto 94 !name too long
				if(i.gt.istart) file=name(istart:i-1)
				istart=i
				lwhere=4
			else if(lwhere.eq.4) then	   !extension
				if(i-istart.gt.80) goto 94 !name too long
				extension=name(istart:i-1)
				istart=i
				lwhere=5
			else
				goto 98	!. not allowed here
			end if
		end if
	else if(c.eq.';') then				!version
		if(lwhere.le.3) then			!file
			if(i-istart.gt.80) goto 94	!name too long
			if(i.gt.istart) file=name(istart:i-1)
			istart=i
			lwhere=5
		else if(lwhere.eq.4) then		!extension
			if(i-istart.gt.80) goto 94	!name too long
			extension=name(istart:i-1)
			istart=i
			lwhere=5
		else
			goto 98		!; not allowed here
		end if
	end if
c
	i=i+1
	if(i.le.nname) c=name(i:i)
c
	end do
c
	if(indir.eq.1) goto 98		!directory not closed
c
	if(lwhere.le.3) then			!end file
		if(i-istart.gt.80) goto 94	!name too long
		if(i.gt.istart) file=name(istart:i-1)
		istart=i
		lwhere=4
	else if(lwhere.eq.4) then		!end extension
		if(i-istart.gt.80) goto 94	!name too long
		extension=name(istart:i-1)
		istart=i
		lwhere=5
	else if(lwhere.eq.5) then		!end version
		if(i-istart.gt.80) goto 94	!name too long
		version=name(istart:i-1)
		istart=i
		lwhere=6
	end if
c
	if(ntype.lt.1) goto 97
	c=type(1:1)
	if(c.eq.'n'.or.c.eq.'N') then
		help=node
	else if(c.eq.'d'.or.c.eq.'D') then
		if(ntype.lt.2) goto 97
		c=type(2:2)
		if(c.eq.'e'.or.c.eq.'E') then
			help=device
		else if(c.eq.'i'.or.c.eq.'I') then
			help=directory
		else if(c.eq.' ') then
			goto 97
		else
			goto 96
		end if
	else if(c.eq.'f'.or.c.eq.'F') then
		help=file
	else if(c.eq.'e'.or.c.eq.'E') then
		help=extension
	else if(c.eq.'v'.or.c.eq.'V') then
		help=version
	else if(c.eq.' ') then
		goto 97
	else
		goto 96
	end if
c
	do i=80,1,-1
	if(help(i:i).ne.' ') goto 3
	end do
	i=0
    3	continue
c
	if(i.gt.nback) goto 95	!back to small
c
	if(i.gt.0) back(1:i)=help(1:i)
c
	iscnam1=i
c
	goto 100
c
c +-----------------------------------------------------------+
c !	error traps
c +-----------------------------------------------------------+
c
   94	continue
	write(6,*) 'name too long (maximum 80 characters) :'
	write(6,*) name
	iscnam1=-6
	goto 100
   95	continue
	write(6,*) 'dimension of back too small for variable :'
	write(6,*) help
	iscnam1=-5
	goto 100
   96	continue
	write(6,*) 'type not recognized :'
	write(6,*) type
	iscnam1=-4
	goto 100
   97	continue
	write(6,*) 'ambigous type : ',type
	write(6,*) 'supply more characters'
	iscnam1=-3
	goto 100
   98	continue
	write(6,*) 'error in name :'
	write(6,*) name
	iscnam1=-2
	goto 100
   99	continue
	iscnam1=-1
	goto 100
c
c +--------------------------------------------------------------+
c !	test output and end
c +--------------------------------------------------------------+
c
  100	continue
c
	if(btest) then
		write(6,*) node
		write(6,*) device
		write(6,*) directory
		write(6,*) file
		write(6,*) extension
		write(6,*) version
		write(6,*) 'iscnam1,i,istart,indir,lwhere : '
		write(6,*)  iscnam1,i,istart,indir,lwhere
	end if
c
	return
	end
c
c****************************************************************
c
	function iscnam(name,type,back)
c
c scanns a file name and gives back the desired component
c
c version for UNIX file format; node, disk & version retained for
c    compatibility with VMS
c
c name		file name
c type		the desired portion of the file name
c		...type has to be one of these keywords :
c		... [ node  device  directory  file  extension  version ]
c		...non ambigous abreviations can be used
c back		the returned portion of the file name
c		...nothing is returned if the desired portion
c		...of the file is missing
c		...node and device are returned with their colons
c		...directory is returned with their brackets
c		...extension and version with their leading dot
c		...or semicolon
c iscnam	length of the character string (return value)
c		...(0 if nothing is returned)
c		...if iscnam is negativ, an error has occured
c 		... error code :
c			-1: nothing in name
c			-2: error in file name
c			-3: ambigous type
c			-4: type not recognized
c			-5: dimension of back too small
c			-6: name too long (maximum 80 characters)
c
	character*(*) name,type,back
c
	character*80 node,device,directory,file,extension,version
	character*80 help
	character*1  c
	logical btest
c
	data btest /.false./	!.true. ==> test output on terminal
c
	nname=len(name)
	ntype=len(type)
	nback=len(back)
c
	node=' '
	device=' '
	directory=' '
	file=' '
	extension=' '
	version=' '
	back=' '
c	
	do i=1,nname
	if(name(i:i).ne.' ') goto 1
	end do
	goto 99
    1	ifirst=i
c
	do i=nname,1,-1
	if(name(i:i).ne.' ') goto 2
	end do
	goto 99
    2	ilast=i
c
c find directory
c
	idir=ifirst-1
	do i=ilast,ifirst,-1
	  if(name(i:i).eq.'/') goto 5
	end do
    5	continue
	if(i.ge.ifirst) then
	  if(i-ifirst.ge.80) goto 94
	  idir=i
	  directory=name(ifirst:i)
	end if
c
c find extension
c
	iext=ilast+1
	do i=ilast,ifirst,-1
	  if(name(i:i).eq.'.') goto 6
	end do
    6	continue
c	if(i.gt.idir+1) then	!idir+1 to allow names like .login
	if(i.gt.idir) then
	  if(ilast-i.ge.80) goto 94
	  iext=i
	  extension=name(i:ilast)
	end if
c
c copy file name
c
c	write(6,*) idir,iext
c	write(6,*) directory
c	write(6,*) extension

	if(iext-idir-2.gt.80) then
	  goto 94
	else if(iext-idir-2.ge.0) then
	  file=name(idir+1:iext-1)
	end if
c
	if(ntype.lt.1) goto 97
	c=type(1:1)
	if(c.eq.'n'.or.c.eq.'N') then
		help=node
	else if(c.eq.'d'.or.c.eq.'D') then
		if(ntype.lt.2) goto 97
		c=type(2:2)
		if(c.eq.'e'.or.c.eq.'E') then
			help=device
		else if(c.eq.'i'.or.c.eq.'I') then
			help=directory
		else if(c.eq.' ') then
			goto 97
		else
			goto 96
		end if
	else if(c.eq.'f'.or.c.eq.'F') then
		help=file
	else if(c.eq.'e'.or.c.eq.'E') then
		help=extension
	else if(c.eq.'v'.or.c.eq.'V') then
		help=version
	else if(c.eq.' ') then
		goto 97
	else
		goto 96
	end if
c
	do i=80,1,-1
	if(help(i:i).ne.' ') goto 3
	end do
	i=0
    3	continue
c
	if(i.gt.nback) goto 95	!back to small
c
	if(i.gt.0) back(1:i)=help(1:i)
c
	iscnam=i
c
	goto 100
c
c +-----------------------------------------------------------+
c !	error traps
c +-----------------------------------------------------------+
c
   94	continue
	write(6,*) 'name too long (maximum 80 characters) :'
	write(6,*) name
	iscnam=-6
	goto 100
   95	continue
	write(6,*) 'dimension of back too small for variable :'
	write(6,*) help
	iscnam=-5
	goto 100
   96	continue
	write(6,*) 'type not recognized :'
	write(6,*) type
	iscnam=-4
	goto 100
   97	continue
	write(6,*) 'ambigous type : ',type
	write(6,*) 'supply more characters'
	iscnam=-3
	goto 100
c  98	continue	!$$ALPHA
c	write(6,*) 'error in name :'
c	write(6,*) name
c	iscnam=-2
c	goto 100
   99	continue
	iscnam=-1
	goto 100
c
c +--------------------------------------------------------------+
c !	test output and end
c +--------------------------------------------------------------+
c
  100	continue
c
	if(btest) then
		write(6,*) node
		write(6,*) device
		write(6,*) directory
		write(6,*) file
		write(6,*) extension
		write(6,*) version
		write(6,*) 'iscnam,i,istart,indir,lwhere : '
		write(6,*)  iscnam,i,istart,indir,lwhere
	end if
c
	return
	end
