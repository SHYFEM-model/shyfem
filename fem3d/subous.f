c
c $Id: subous.f,v 1.3 2004/09/28 13:31:55 georg Exp $
c
c utility routines to read/write OUS file - file type 161
c
c contents :
c
c subroutine inious
c subroutine setous(iunit,nvers,nkn,nel,nlv)
c subroutine getous(iunit,nvers,nkn,nel,nlv)
c
c subroutine dimous(iunit,nkndim,neldim,nlvdim)
c
c subroutine rfous(iunit,nvers,nkn,nel,nlv,href,hzmin,title,ierr)
c subroutine wfous(iunit,nvers,nkn,nel,nlv,href,hzmin,title,ierr)
c subroutine rsous(iunit,ilhv,hlv,hev,ierr)
c subroutine wsous(iunit,ilhv,hlv,hev,ierr)
c subroutine rdous(iunit,it,nlvdim,ilhv,z,ze,ut,vt,ierr)
c subroutine wrous(iunit,it,nlvdim,ilhv,z,ze,ut,vt,ierr)
c
c revision log :
c
c 21.08.2003	ggu	version 1 implemented
c 01.09.2003	ggu	first revision
c 02.09.2003	ggu	last bug fixes (nvers=3 -> nvers=1)
c 22.09.2004	ggu	bug fix in rdous/wrous -> ie instead of k
c 08.06.2011	ggu	new routine delous(), check for end in read
c
c notes :
c
c format of file:
c
c version 1
c
c	mtype,nvers
c	nkn,nel,nlv
c	href,hzmin
c	title
c
c	(ilhv(ie),ie=1,nel)
c	(hlv(l),l=1,nlv)
c	(hev(ie),ie=1,nel)
c	
c	it
c	(z(k),k=1,nkn)
c	(ze(i),i=1,3*nel)
c	((ut(l,ie),l=1,ilhv(k)),ie=1,nel)
c	((vt(l,ie),l=1,ilhv(k)),ie=1,nel)
c
c************************************************************

	subroutine inious

c sets up initial common block - internal routine

	implicit none

c parameters
	integer ftype,maxvers
	parameter(ftype=27,maxvers=1)
	integer ndim,nitdim
	parameter(ndim=10,nitdim=5)
c common
	integer mtype,maxver
	common /ouscom/ mtype,maxver
	integer ousitem,ousvar(0:nitdim,ndim)
	common /ousvar/ousitem,ousvar
c local
	integer i,n
c save
	logical binit
	save binit
	save /ouscom/
	save /ousvar/
c data
	data binit /.false./

	if( binit ) return

	binit = .true.

	mtype = ftype
	maxver = maxvers

	ousitem = 0
	do n=1,ndim
	  do i=0,nitdim
	    ousvar(i,n) = 0
	  end do
	end do

	end

c************************************************************

	subroutine setous(iunit,nvers,nkn,nel,nlv)

c sets up parameter common block - internal routine

	implicit none

c arguments
	integer iunit,nvers,nkn,nel,nlv
c parameters
	integer ndim,nitdim
	parameter(ndim=10,nitdim=5)
c common
	integer ousitem,ousvar(0:nitdim,ndim)
	common /ousvar/ousitem,ousvar

	integer n

c we do not check if unit has already been opened -> open with ifileo

        do n=1,ousitem
          if( ousvar(0,n) .eq. 0 ) goto 1
          if( ousvar(0,n) .eq. iunit ) goto 99
        end do
    1   continue
        if( n .gt. ousitem ) ousitem = n

	if( ousitem .gt. ndim ) then
	   stop 'error stop setous: ndim'
	end if

	ousvar(0,n) = iunit
	ousvar(1,n) = nvers
	ousvar(2,n) = nkn
	ousvar(3,n) = nel
	ousvar(4,n) = nlv

	return
   99   continue
        write(6,*) 'unit = ',iunit
        stop 'error stop setous: unit already open - please close first'
	end

c************************************************************

	subroutine getous(iunit,nvers,nkn,nel,nlv)

c gets parameter common block - internal routine

	implicit none

c arguments
	integer iunit,nvers,nkn,nel,nlv
c parameters
	integer ndim,nitdim
	parameter(ndim=10,nitdim=5)
c common
	integer ousitem,ousvar(0:nitdim,ndim)
	common /ousvar/ousitem,ousvar

	integer n

	do n=1,ousitem
	  if( ousvar(0,n) .eq. iunit ) goto 1
	end do

	write(6,*) 'Cannot read on unit ',iunit
	write(6,*) 'File is not initialized.'
	stop 'error stop getous: no initialization'
    1	continue

	nvers = ousvar(1,n)
	nkn   = ousvar(2,n)
	nel   = ousvar(3,n)
	nlv   = ousvar(4,n)

	end

c************************************************************

        subroutine delous(iunit)

c closes ous file internal structure - internal routine
c
c please note that the file has still to be closed manually

        implicit none

c arguments
        integer iunit
c parameters
        integer ndim,nitdim
        parameter(ndim=10,nitdim=5)
c common
        integer ousitem,ousvar(0:nitdim,ndim)
        common /ousvar/ousitem,ousvar

        integer n

        do n=1,ousitem
          if( ousvar(0,n) .eq. iunit ) goto 1
        end do

        write(6,*) 'Cannot close unit ',iunit
        write(6,*) 'File is not open.'
        stop 'error stop delous: file not open'
    1   continue

        ousvar(0,n) = 0

        end

c************************************************************

	subroutine dimous(iunit,nkndim,neldim,nlvdim)

c checks dimension of arrays

	implicit none

c arguments
	integer iunit,nkndim,neldim,nlvdim
c local
	integer nvers,nkn,nel,nlv

	call getous(iunit,nvers,nkn,nel,nlv)

        if( nkn .gt. nkndim ) goto 99
        if( nel .gt. neldim ) goto 99
        if( nlv .gt. nlvdim ) goto 99

	return
   99   continue
        write(6,*) 'nkn,nkndim : ',nkn,nkndim
        write(6,*) 'nel,neldim : ',nel,neldim
        write(6,*) 'nlv,nlvdim : ',nlv,nlvdim
        stop 'error stop dimous: dimension error'
	end

c************************************************************

	subroutine rfous	(iunit,nvers
     +				,nkn,nel,nlv
     +				,href,hzmin
     +				,title
     +				,ierr
     +				)

c reads first record of OUS file
c
c nvers		on entry maximal version that can be read
c		-> must be an input, used to check the corectness
c		.. of the call parameters
c		on return actual version read

	implicit none

c arguments
	integer iunit,nvers
	integer nkn,nel,nlv
	real href,hzmin
	character*80 title
	integer ierr
c common
	integer mtype,maxver
	common /ouscom/ mtype,maxver
c local
	integer ntype,irec

c initialize

	call inious

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
	   read(iunit,err=99)	 nkn,nel,nlv
	   read(iunit,err=99)	 href,hzmin
	   read(iunit,err=99)	 title
	else
	   stop 'error stop rfous: internal error (1)'
	end if

	call setous(iunit,nvers,nkn,nel,nlv)

	ierr=0

	return
   99	continue
	write(6,*) 'rfous: Error encountered while'
	write(6,*) 'reading record number ',irec
	write(6,*) 'of OUS file header'
	write(6,*) 'nvers = ',nvers
	ierr=99
	return
   98	continue
	write(6,*) 'rfous: Version not recognized : ',nvers
	ierr=98
	return
   97	continue
	write(6,*) 'rfous: Wrong type of file : ',ntype
	write(6,*) 'Expected ',mtype
	ierr=97
	return
   96	continue
	write(6,*) 'rfous: Cannot rewind file for unit : ',iunit
	ierr=96
	return
   95	continue
	write(6,*) 'rfous: Old function call ',nvers
	write(6,*) 'Please adjust call to rfous and recompile'
	ierr=95
	return
   91	continue
	write(6,*) 'rfous: File is empty'
	ierr=91
	return
	end

c********************************************************************

	subroutine wfous	(iunit,nvers
     +				,nkn,nel,nlv
     +				,href,hzmin
     +				,title
     +				,ierr
     +				)

c writes first record of OUS file
c
c nvers		on entry maximal version
c		-> must be an input, used to check the corectness
c		.. of the call parameters

	implicit none

c arguments
	integer iunit,nvers
	integer nkn,nel,nlv
	real href,hzmin
	character*80 title
	integer ierr
c common
	integer mtype,maxver
	common /ouscom/ mtype,maxver

c initialize

	call inious

c control newest version number for call

	if( nvers .ne. maxver ) goto 95

	rewind(iunit)

	write(iunit)		mtype,maxver
	write(iunit)		nkn,nel,nlv
	write(iunit)		href,hzmin
	write(iunit)		title

	call setous(iunit,nvers,nkn,nel,nlv)

	ierr=0

	return
   95	continue
	write(6,*) 'wfous: Old function call ',nvers
	write(6,*) 'wfous: should be maxver = ',maxver
	write(6,*) 'Please adjust call to wfous and recompile'
	ierr=95
	return
	end

c************************************************************

	subroutine rsous(iunit,ilhv,hlv,hev,ierr)

c reads second record of OUS file

	implicit none

c arguments
	integer iunit
	integer ilhv(1)
	real hlv(1)
	real hev(1)
	integer ierr
c common
	integer mtype,maxver
	common /ouscom/ mtype,maxver
c local
	integer l,ie
	integer nvers,nkn,nel,nlv

	call getous(iunit,nvers,nkn,nel,nlv)

c read records

	if( nvers .eq. 1 ) then
	  read(iunit,err=99) (ilhv(ie),ie=1,nel)
	  read(iunit,err=99) (hlv(l),l=1,nlv)
	  read(iunit,err=99) (hev(ie),ie=1,nel)
	else
	   stop 'error stop rsous: internal error (1)'
	end if

	ierr = 0

	return
   99	continue
	write(6,*) 'rsous: Error encountered while'
	write(6,*) 'reading second part of OUS file header'
	ierr=99
	return
	end

c************************************************************

	subroutine wsous(iunit,ilhv,hlv,hev,ierr)

c writes second record of OUS file

	implicit none

c arguments
	integer iunit
	integer ilhv(1)
	real hlv(1)
	real hev(1)
	integer ierr
c common
	integer mtype,maxver
	common /ouscom/ mtype,maxver
c local
	integer l,ie
	integer nvers,nkn,nel,nlv

	call getous(iunit,nvers,nkn,nel,nlv)

c write records

	write(iunit) (ilhv(ie),ie=1,nel)
	write(iunit) (hlv(l),l=1,nlv)
	write(iunit) (hev(ie),ie=1,nel)

	ierr = 0

	return
	end

c************************************************************

	subroutine rdous(iunit,it,nlvdim,ilhv,z,ze,ut,vt,ierr)

c reads data record of OUS file

	implicit none

c arguments
	integer iunit,it
	integer nlvdim
	integer ilhv(1)
	real z(1)
	real ze(1)
	real ut(nlvdim,1)
	real vt(nlvdim,1)
	integer ierr
c common
	integer mtype,maxver
	common /ouscom/ mtype,maxver
c local
	integer l,k,ie,i
	integer nvers,nkn,nel,nlv

	call getous(iunit,nvers,nkn,nel,nlv)

	if( nvers .eq. 1 ) then
	   read(iunit,end=88,err=98) it
	   read(iunit,end=99,err=99) (z(k),k=1,nkn)
	   read(iunit,end=99,err=99) (ze(i),i=1,3*nel)
	   read(iunit,end=99,err=99) ((ut(l,ie),l=1,ilhv(ie)),ie=1,nel)
	   read(iunit,end=99,err=99) ((vt(l,ie),l=1,ilhv(ie)),ie=1,nel)
	else
	   stop 'error stop rdous: internal error (1)'
	end if

	ierr=0

	return
   88	continue
	ierr=-1
	return
   98	continue
	write(6,*) 'rdous: Error while reading'
	write(6,*) 'time record of OUS file'
	ierr=98
	return
   99	continue
	write(6,*) 'rdous: Error while reading'
	write(6,*) 'data record of OUS file'
	write(6,*) 'it = ',it
	ierr=99
	return
	end

c************************************************************

	subroutine wrous(iunit,it,nlvdim,ilhv,z,ze,ut,vt,ierr)

c writes data record of OUS file

	implicit none

c arguments
	integer iunit,it
	integer nlvdim
	integer ilhv(1)
	real z(1)
	real ze(1)
	real ut(nlvdim,1)
	real vt(nlvdim,1)
	integer ierr
c common
	integer mtype,maxver
	common /ouscom/ mtype,maxver
c local
	integer l,k,ie,i
	integer nvers,nkn,nel,nlv

	call getous(iunit,nvers,nkn,nel,nlv)

	write(iunit) it
	write(iunit) (z(k),k=1,nkn)
	write(iunit) (ze(i),i=1,3*nel)
	write(iunit) ((ut(l,ie),l=1,ilhv(ie)),ie=1,nel)
	write(iunit) ((vt(l,ie),l=1,ilhv(ie)),ie=1,nel)

	ierr=0

	return
	end

c************************************************************

