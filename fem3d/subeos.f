
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

c utility routines to read/write EOS file - file type 167
c
c contents :
c
c subroutine inieos
c subroutine seteos(iunit,nvers,nkn,nel,nlv,nvar)
c subroutine geteos(iunit,nvers,nkn,nel,nlv,nvar)
c
c subroutine dimeos(iunit,nknddi,nelddi,nlvddi)
c
c subroutine rfeos(iunit,nvers,nkn,nel,nlv,nvar,title,ierr)
c subroutine wfeos(iunit,nvers,nkn,nel,nlv,nvar,title,ierr)
c subroutine rseos(iunit,ilhv,hlv,hev,ierr)
c subroutine wseos(iunit,ilhv,hlv,hev,ierr)
c subroutine rdeos(iunit,it,ivar,nlvddi,ilhv,c,ierr)
c subroutine wreos(iunit,it,ivar,nlvddi,ilhv,c,ierr)
c
c revision log :
c
c 31.08.2011	ggu	new routines EOS
c 01.09.2011	ggu	changed VERS_6_1_32
c 06.10.2011	ggu	bug fix for ilhv -> nel_lev
c 18.10.2011	ggu	changed VERS_6_1_33
c 10.07.2015	ggu	changed VERS_7_1_50
c 18.12.2018	ggu	changed VERS_7_5_52
c 16.02.2019	ggu	changed VERS_7_5_60
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
c	(ilhv(k),k=1,nel)			empty if nlv <= 1
c	(hlv(k),k=1,nlv)			empty if nlv <= 1
c	(hev(k),k=1,nel)
c	
c	it,ivar
c	((c(l,ie),l=1,ilhv(ie)),ie=1,nel)	if nlv > 1
c	(c(1,ie),ie=1,nel)			if nlv <= 1
c
c version 2	not existing
c
c version 1	not existing
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
c - routine to close file (read and write) -> cleos()
c - remove data structure from eosvar when closed -> iunit=0
c - in seteos look for first available iunit
c - in seteos check if iunit is already open
c
c************************************************************

	subroutine inieos

c sets up initial common block - internal routine

	implicit none

c parameters
	integer ftype,maxvers
	parameter(ftype=167,maxvers=3)
	integer ndim,nitdim
	parameter(ndim=30,nitdim=5)
c common
	integer mtype,maxver
	common /eoscom/ mtype,maxver
	integer eositem,eosvar(0:nitdim,ndim)
	common /eosvar/eositem,eosvar
c local
	integer i,n
c save
	logical binit
	save binit
	save /eoscom/
	save /eosvar/
c data
	data binit /.false./

	if( binit ) return

	binit = .true.

	mtype = ftype
	maxver = maxvers

	eositem = 0
	do n=1,ndim
	  do i=0,nitdim
	    eosvar(i,n) = 0
	  end do
	end do

	end

c************************************************************

	subroutine seteos(iunit,nvers,nkn,nel,nlv,nvar)

c sets up parameter common block - internal routine

	implicit none

c arguments
	integer iunit,nvers,nkn,nel,nlv,nvar
c parameters
	integer ndim,nitdim
	parameter(ndim=30,nitdim=5)
c common
	integer eositem,eosvar(0:nitdim,ndim)
	common /eosvar/eositem,eosvar
c local
	integer n

c we do not check if unit has already been opened -> open with ifileo

	do n=1,eositem
	  if( eosvar(0,n) .eq. 0 ) goto 1
	  if( eosvar(0,n) .eq. iunit ) goto 99
	end do
    1	continue
	if( n .gt. eositem ) eositem = n

	if( eositem .gt. ndim ) then
	   stop 'error stop seteos: ndim'
	end if

	eosvar(0,n) = iunit
	eosvar(1,n) = nvers
	eosvar(2,n) = nkn
	eosvar(3,n) = nel
	eosvar(4,n) = nlv
	eosvar(5,n) = nvar

	return
   99	continue
	write(6,*) 'unit = ',iunit
	stop 'error stop seteos: unit already open - please close first'
	end

c************************************************************

	subroutine geteos(iunit,nvers,nkn,nel,nlv,nvar)

c gets parameter common block - internal routine

	implicit none

c arguments
	integer iunit,nvers,nkn,nel,nlv,nvar
c parameters
	integer ndim,nitdim
	parameter(ndim=30,nitdim=5)
c common
	integer eositem,eosvar(0:nitdim,ndim)
	common /eosvar/eositem,eosvar

	integer n

	do n=1,eositem
	  if( eosvar(0,n) .eq. iunit ) goto 1
	end do

	write(6,*) 'Cannot read on unit ',iunit
	write(6,*) 'File is not initialized.'
	stop 'error stop geteos: no initialization'
    1	continue

	nvers = eosvar(1,n)
	nkn   = eosvar(2,n)
	nel   = eosvar(3,n)
	nlv   = eosvar(4,n)
	nvar  = eosvar(5,n)

	end

c************************************************************

	subroutine deleos(iunit)

c closes eos file internal structure - internal routine
c
c please note that the file has still to be closed manually

	implicit none

c arguments
	integer iunit
c parameters
	integer ndim,nitdim
	parameter(ndim=30,nitdim=5)
c common
	integer eositem,eosvar(0:nitdim,ndim)
	common /eosvar/eositem,eosvar

	integer n

	do n=1,eositem
	  if( eosvar(0,n) .eq. iunit ) goto 1
	end do

	write(6,*) 'Cannot close unit ',iunit
	write(6,*) 'File is not open.'
	stop 'error stop deleos: file not open'
    1	continue

	eosvar(0,n) = 0

	end

c************************************************************

	subroutine dimeos(iunit,nknddi,nelddi,nlvddi)

c checks dimension of arrays

	implicit none

c arguments
	integer iunit,nknddi,nelddi,nlvddi
c local
	integer nvers,nkn,nel,nlv,nvar

	call geteos(iunit,nvers,nkn,nel,nlv,nvar)

        if( nkn .gt. nknddi ) goto 99
        if( nel .gt. nelddi ) goto 99
        if( nlv .gt. nlvddi ) goto 99

	return
   99   continue
        write(6,*) 'nkn,nknddi : ',nkn,nknddi
        write(6,*) 'nel,nelddi : ',nel,nelddi
        write(6,*) 'nlv,nlvddi : ',nlv,nlvddi
        stop 'error stop dimeos: dimension error'
	end

c************************************************************
c************************************************************
c************************************************************
c************************************************************
c************************************************************

	subroutine rfeos	(iunit,nvers
     +				,nkn,nel,nlv,nvar
     +				,title
     +				,ierr
     +				)

c reads first record of EOS file
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
	common /eoscom/ mtype,maxver
c local
	integer ntype,irec
	character*80 line

c initialize

	call inieos

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
	if( nvers .eq. 3 ) then
	   read(iunit,err=99)	 nkn,nel,nlv,nvar
	   read(iunit,err=99)	 line
	else
	   write(6,*) 'version = ',nvers
	   stop 'error stop rfeos: internal error (1)'
	end if

	call seteos(iunit,nvers,nkn,nel,nlv,nvar)
	title = line

	ierr=0

	return
   99	continue
	write(6,*) 'rfeos: Error encountered while'
	write(6,*) 'reading record number ',irec
	write(6,*) 'of EOS file header'
	write(6,*) 'nvers = ',nvers
	ierr=99
	return
   98	continue
	write(6,*) 'rfeos: Version not recognized : ',nvers
	ierr=98
	return
   97	continue
	write(6,*) 'rfeos: Wrong type of file : ',ntype
	write(6,*) 'Expected ',mtype
	ierr=97
	return
   96	continue
	write(6,*) 'rfeos: Cannot rewind file for unit : ',iunit
	ierr=96
	return
   95	continue
	write(6,*) 'rfeos: Old function call ',nvers
	write(6,*) 'nvers = ',nvers,'   maxver = ',maxver
	write(6,*) 'Please adjust call to rfeos and recompile'
	ierr=95
	return
   91	continue
	write(6,*) 'rfeos: File is empty'
	ierr=91
	return
	end

c********************************************************************

	subroutine wfeos	(iunit,nvers
     +				,nkn,nel,nlv,nvar
     +				,title
     +				,ierr
     +				)

c writes first record of EOS file
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
	common /eoscom/ mtype,maxver

	character*80 line

        line = ' '
        line = title            !must be 80 chars

c initialize

	call inieos

c control newest version number for call

	if( nvers.ne.maxver ) goto 95

	rewind(iunit)

	write(iunit)		mtype,maxver
	write(iunit)		nkn,nel,nlv,nvar
	write(iunit)		line

	call seteos(iunit,nvers,nkn,nel,nlv,nvar)

	ierr=0

	return
   95	continue
	write(6,*) 'wfeos: Old function call'
	write(6,*) 'nvers = ',nvers,'   maxver = ',maxver
	write(6,*) 'Please adjust call to rfeos and recompile'
	ierr=95
	return
	end

c************************************************************

	subroutine rseos(iunit,ilhv,hlv,hev,ierr)

c reads second record of EOS file

	implicit none

c arguments
	integer iunit
	integer ilhv(1)
	real hlv(1)
	real hev(1)
	integer ierr
c common
	integer mtype,maxver
	common /eoscom/ mtype,maxver
c local
	integer l,ie
	integer nvers,nkn,nel,nlv,nvar
	integer nel_lev

	call geteos(iunit,nvers,nkn,nel,nlv,nvar)

c only one layer

	nel_lev = nel
	if( nlv .le. 1 ) then
	  do ie=1,nel
	    ilhv(ie) = 1
	  end do
	  hlv(1) = 10000.

	  nel_lev = 0
	  nlv = 0
	end if

c read records

	if( nvers .eq. 3 ) then
	  read(iunit,err=99) (ilhv(ie),ie=1,nel_lev)
	  read(iunit,err=99) (hlv(l),l=1,nlv)
	  read(iunit,err=99) (hev(ie),ie=1,nel)
	else
	   write(6,*) 'version = ',nvers
	   stop 'error stop rseos: internal error (1)'
	end if

	ierr = 0

	return
   99	continue
	write(6,*) 'rseos: Error encountered while'
	write(6,*) 'reading second part of EOS file header'
	ierr=99
	return
	end

c************************************************************

	subroutine wseos(iunit,ilhv,hlv,hev,ierr)

c writes second record of EOS file

	implicit none

c arguments
	integer iunit
	integer ilhv(1)
	real hlv(1)
	real hev(1)
	integer ierr
c common
	integer mtype,maxver
	common /eoscom/ mtype,maxver
c local
	integer l,ie
	integer nvers,nkn,nel,nlv,nvar
	integer nel_lev

	call geteos(iunit,nvers,nkn,nel,nlv,nvar)

c only one layer

	nel_lev = nel
	if( nlv .le. 1 ) then
	  nlv = 0
	  nel_lev = 0
	end if

c write records

	write(iunit) (ilhv(ie),ie=1,nel_lev)
	write(iunit) (hlv(l),l=1,nlv)
	write(iunit) (hev(ie),ie=1,nel)

	ierr = 0

	return
	end

c************************************************************

	subroutine rdeos(iunit,it,ivar,nlvddi,ilhv,c,ierr)

c reads data record of EOS file

	implicit none

c arguments
	integer iunit,it,ivar
	integer nlvddi
	integer ilhv(1)
	real c(nlvddi,1)
	integer ierr
c common
	integer mtype,maxver
	common /eoscom/ mtype,maxver
c local
	integer l,ie
	integer nvers,nkn,nel,nlv,nvar

	call geteos(iunit,nvers,nkn,nel,nlv,nvar)

	if( nvers .eq. 3 ) then
	   read(iunit,end=88,err=98) it,ivar
	   if( nlv .le. 1 ) then
             read(iunit,end=99,err=99) (c(1,ie),ie=1,nel)
	   else
             call linear2read(iunit,nlvddi,nel,ilhv,c,ierr)
             if( ierr /= 0 ) goto 99
             !read(iunit,end=99,err=99) ((c(l,ie),l=1,ilhv(ie)),ie=1,nel)
	   end if
	else
	   write(6,*) 'version = ',nvers
	   stop 'error stop rdeos: internal error (1)'
	end if

	ierr=0

	return
   88	continue
	ierr=-1
	return
   98	continue
	write(6,*) 'rdeos: Error while reading'
	write(6,*) 'time record of EOS file'
	ierr=98
	return
   99	continue
	write(6,*) 'rdeos: Error while reading'
	write(6,*) 'data record of EOS file'
	write(6,*) 'it = ',it,'  ivar = ',ivar
	ierr=99
	return
	end

c************************************************************

	subroutine wreos(iunit,it,ivar,nlvddi,ilhv,c,ierr)

c writes data record of EOS file

	implicit none

c arguments
	integer iunit,it,ivar
	integer nlvddi
	integer ilhv(1)
	real c(nlvddi,1)
	integer ierr
c common
	integer mtype,maxver
	common /eoscom/ mtype,maxver
c local
	integer l,ie
	integer nvers,nkn,nel,nlv,nvar

	call geteos(iunit,nvers,nkn,nel,nlv,nvar)

	write(iunit) it,ivar

	if( nlv .le. 1 ) then
	  write(iunit) (c(1,ie),ie=1,nel)
	else
          call linear2write(iunit,nlvddi,nel,ilhv,c,ierr)
	  if( ierr /= 0 ) stop 'error stop wreos'
	  !write(iunit) ((c(l,ie),l=1,ilhv(ie)),ie=1,nel)
	end if

	ierr=0

	return
	end

c************************************************************

	subroutine rheos	(iunit,nvers
     +				,nknddi,nelddi,nlvddi
     +				,nkn,nel,nlv,nvar
     +				,ilhv,hlv,hev
     +				,title
     +				)

c reads all headers of EOS file
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
	integer ilhv(1)
	real hlv(1)
	real hev(1)
	character*(*) title
c local
	integer ierr,l

        call rfeos(iunit,nvers,nkn,nel,nlv,nvar,title,ierr)
        if(ierr.ne.0) goto 99

        write(6,*) 'nvers    : ',nvers
        write(6,*) 'nkn,nel  : ',nkn,nel
        write(6,*) 'nlv,nvar : ',nlv,nvar
        write(6,*) 'title    : ',title

        call dimeos(iunit,nknddi,nelddi,nlvddi)

        call rseos(iunit,ilhv,hlv,hev,ierr)
        if(ierr.ne.0) goto 98

        write(6,*) 'Available levels: ',nlv
        write(6,*) (hlv(l),l=1,nlv)

	return
   98	continue
	write(6,*) 'ierr: ',ierr
	stop 'error stop rheos: error reading second header'
   99	continue
	write(6,*) 'ierr: ',ierr
	stop 'error stop rheos: error reading first header'
	end

c************************************************************

        subroutine wheos        (iunit,nvers
     +                          ,nkn,nel,nlv,nvar
     +				,ilhv,hlv,hev
     +                          ,title
     +                          )

c writes all headers of EOS file
c
c nvers         on entry maximal version
c               -> must be an input, used to check the corectness
c               .. of the call parameters

        implicit none

c arguments
        integer iunit,nvers
        integer nkn,nel,nlv,nvar
	integer ilhv(1)
	real hlv(1)
	real hev(1)
        character*(*) title

        integer ierr

        call wfeos(iunit,nvers,nkn,nel,nlv,nvar,title,ierr)
        if(ierr.ne.0) goto 1

	call wseos(iunit,ilhv,hlv,hev,ierr)
        if(ierr.ne.0) goto 1

	return
    1	continue
	write(6,*) 'ierr: ',ierr
	stop 'error stop wheos'
	end

c************************************************************

