c
c $Id: subext.f,v 1.7 2001/11/16 07:35:43 georg Exp $
c
c utility routines to read/write EXT file - file type 71
c
c contents :
c
c subroutine iniext
c subroutine rfext(iunit,nvmax,nvers,npoint,href,hzmin,title,ierr)
c subroutine wfext(iunit,nvmax,nvers,npoint,href,hzmin,title,ierr)
c subroutine rsext(iunit,kpoint,ipoint,hdep,x,y,ierr)
c subroutine wsext(iunit,kpoint,ipoint,hdep,x,y,ierr)
c subroutine rdext(iunit,it,u,v,z,ierr)
c subroutine wrext(iunit,it,u,v,z,ierr)
c
c revision log :
c
c 20.05.1998	ggu	cleaned up a bit
c 04.02.2000	ggu	wrrc77 from newpr to here
c 14.09.2015	ggu	some more helper routines
c 05.10.2015	ggu	handle error in backspace smoothly
c
c notes :
c
c variables used:
c
c	mtype,ftype	type of file (71)
c	maxver,maxvers	newest version
c
c	ierr		error code (0: ok,  <0: EOF,  >0: error)
c
c	iunit		file unit
c	nvmax		version of calling routine (must be maxver)
c	nvers		version of file
c	npoint		number of nodes written
c	href		reference level for water level
c	hzmin		minimum total water depth
c	title		title of run
c
c	kpoint(npoint)	external number of nodes
c	ipoint(npoint)	internal number of nodes
c	hdep(npoint)	depth at nodes
c	x(npoint)	x coordinate of nodes
c	y(npoint)	y coordinate of nodes
c
c	it		actual time of data record
c	u(npoint)	velocity in x direction
c	v(npoint)	velocity in y direction
c	z(npoint)	water level
c
c format of file:
c
c version 7
c
c	mtype,nvers
c	npoint
c	href,hzmin
c	title
c	(kpoint(i),i=1,npoint)
c	(ipoint(i),i=1,npoint)
c	(hdep(i),i=1,npoint)
c	(x(i),i=1,npoint)
c	(y(i),i=1,npoint)
c
c	it
c	(u(i),i=1,npoint)
c	(v(i),i=1,npoint)
c	(z(i),i=1,npoint)
c
c version 3-6
c
c 	nvers
c	npoint,(ipoint(i),i=1,npoint),(hdep(i),i=1,npoint),href,hzmin,title
c
c	it,(u(i),i=1,npoint),(v(i),i=1,npoint),(z(i),i=1,npoint)
c
c version 2
c
c 	nvers
c	npoint,(ipoint(i),i=1,npoint),(hdep(i),i=1,npoint),href,hzmin,title
c
c	float(it),(u(i),i=1,npoint),(v(i),i=1,npoint),(z(i),i=1,npoint)
c
c version 1
c
c	npoint,(ipoint(i),i=1,npoint)
c
c	float(it),(u(i),i=1,npoint),(v(i),i=1,npoint),(z(i),i=1,npoint)
c
c************************************************************
c
c old utility routines to read/write EXT file
c
c contents
c
c-----------------------------------------------------------------------
c function read7(iunit,ndim,nvers,knausm,knaus,hdep,href,hzmin,descrp)
c                       reads first record of file 7
c function writ7(iunit,ndim,nvers,knausm,knaus,hdep,href,hzmin,descrp)
c                       writes first record of file 7
c function rdrc7(iunit,nvers,it,knausm,xv)
c                       reads data record of file 7
c function skrc7(iunit,nvers,it,knausm,xv)
c                       skips one data record of file 7
c function wrrc7(iunit,nvers,it,knausm,knaus,xv)
c                       writes data record of file 7
c function wrrc77(iunit,nvers,it,knausm,knaus,u,v,z)
c			writes data record of extra point file
c-----------------------------------------------------------------------
c
c nvermx		maximum version recognized -> 6
c
c*********************************************************
c*********************************************************
c*********************************************************
c*********************************************************
c*********************************************************
c*********************************************************

	function check_ext_file(file)

	implicit none

        logical check_ext_file
        character*(*) file

        integer nb,nvers,knausm,ierr
        integer ifileo

        check_ext_file = .false.

        nb = ifileo(0,file,'unform','old')
        if( nb .le. 0 ) return
	call ext_check_header(nb,nvers,knausm,ierr)
        close(nb)

        check_ext_file = ierr == 0

	end

c*********************************************************

	subroutine ext_is_ext_file(iunit,nvers)

	implicit none

	integer iunit,nvers

	integer knausm,ierr

	call ext_check_header(iunit,nvers,knausm,ierr)

	if( ierr .ne. 0 ) nvers = 0

	rewind(iunit)

	end

c*********************************************************
c*********************************************************
c*********************************************************

	subroutine ext_peek_header(iunit,nvers,knausm)

	implicit none

	integer iunit,knausm

	integer nvers,ierr

	call ext_check_header(iunit,nvers,knausm,ierr)
	if( ierr .ne. 0 ) then
	  stop 'error stop ext_peek_header: error reading header'
	end if

	rewind(iunit)

	end

c*********************************************************

	subroutine ext_peek_record(iunit,nvers,it,ierr)

	implicit none

	integer iunit,nvers,it,ierr

	integer, parameter :: nvermx = 6
	real tt

	tt = 0

	if(nvers.ge.1.and.nvers.le.2) then
	  read(iunit,iostat=ierr)  tt
	  it = tt
	else if(nvers.ge.3.and.nvers.le.nvermx) then
	  read(iunit,iostat=ierr)  it
	else
	  stop 'error stop ext_peek_record: internal error (1)'
	end if

	if( ierr /= 0 ) return

	backspace(iunit,iostat=ierr)

	if( ierr /= 0 ) ierr = -1	!fake end of file

	end

c*********************************************************

	subroutine ext_check_header(iunit,nvers,knausm,ierr)

c checks version of ext file and returns number of points

	implicit none

	integer iunit,nvers,knausm,ierr

	integer, parameter :: nvermx = 6
	integer kaux,ios,i,it,j
	real haux,tt,xaux
	character*80 title

	nvers = 0
	knausm = 0
	ierr = -1

	read(iunit,iostat=ios) nvers
	if( ios .ne. 0 ) return

	if(nvers.ge.2.and.nvers.le.nvermx) then
	  read(iunit,iostat=ios)   knausm
     +                                  ,(kaux,j=1,knausm)
     +                                  ,(haux,j=1,knausm)
     +                                  ,haux
     +                                  ,haux
     +                                  ,title
	else
	  ios = 1
	end if
	if( ios .ne. 0 ) return

	if(nvers.ge.1.and.nvers.le.2) then
	  read(iunit,iostat=ios)  tt,(xaux,i=1,3*knausm)
	else if(nvers.ge.3.and.nvers.le.nvermx) then
	  read(iunit,iostat=ios)  it,(xaux,i=1,3*knausm)
	end if
	if( ios .ne. 0 ) return

	ierr = 0

	end

c*********************************************************
c*********************************************************
c*********************************************************

	subroutine ext_read_header(iunit,ndim,nvers,knausm,knaus,hdep
     +                          ,href,hzmin,descrp)
	implicit none
	integer iunit,ndim,nvers,knausm
	real href,hzmin
	character*80 descrp
	integer knaus(knausm)
	real hdep(knausm)

	integer ierr
	real read7

	ierr = read7(iunit,ndim,nvers,knausm,knaus,hdep
     +                          ,href,hzmin,descrp)

	if( ierr .ne. 0. ) then
	  stop 'error stop read_ext_header: read error'
	end if

	end

c*********************************************************

	subroutine ext_read_record(iunit,nvers,it,knausm,xv,ierr)
	implicit none
	integer iunit,nvers,it,knausm,ierr
	real xv(3*knausm)

	real rdrc7

	ierr = rdrc7(iunit,nvers,it,knausm,xv)

	end

c*********************************************************

	subroutine ext_write_header(iunit,ndim,nvers,knausm,knaus,hdep
     +                          ,href,hzmin,descrp)
	implicit none
	integer iunit,ndim,nvers,knausm
	real href,hzmin
	character*80 descrp
	integer knaus(knausm)
	real hdep(knausm)

	integer ierr
	real writ7

	ierr = writ7(iunit,ndim,nvers,knausm,knaus,hdep
     +                          ,href,hzmin,descrp)

	end

c*********************************************************

	subroutine ext_write_record(iunit,nvers,it,knausm,xv)
	implicit none
	integer iunit,nvers,it,knausm
	real xv(3*knausm)

	integer ierr,i
	integer knaus(knausm)
	real wrrc7

	do i=1,knausm
	  knaus(i) = i
	end do

	ierr = wrrc7(iunit,nvers,it,knausm,knaus,xv)

	end

c*********************************************************
c*********************************************************
c*********************************************************

	function read7(iunit,ndim,nvers,knausm,knaus,hdep
     +                          ,href,hzmin,descrp)
c
c reads first record of file 7
c
c error codes 11 21 31 35 41 61 71
c
	character*80 descrp
	integer knaus(ndim)
	real hdep(ndim)
c
	nvermx = 6
c
	rewind(iunit,iostat=ios)
c
	if(ios.ne.0) then
		write(6,*) 'Cannot rewind file for unit :',iunit
		read7=71.
		return
	end if
c
c first record
c
	read(iunit,iostat=ios) nvers
c
	if(ios.eq.0) then       ! no read error ==> first version
!
	else if(ios.lt.0) then  !eof
		write(6,*) 'EOF encountered while reading'
		write(6,*) 'first record of file 7 header'
		write(6,*) 'nvers =',nvers
		read7=21.
		return
	end if
c
c second record
c
	if(nvers.eq.1) then
		hzmin=0.05
		descrp=' '
		do j=1,knausm
		hdep(j)=100000.         !no dry areas
		end do
	else if(nvers.ge.2.and.nvers.le.nvermx) then
		read(iunit,iostat=ios)   knausm
     +                                  ,(knaus(j),j=1,knausm)
     +                                  ,(hdep(j),j=1,knausm)
     +                                  ,href
     +                                  ,hzmin
     +                                  ,descrp
	else
		write(6,*) 'version not recognized : ',nvers
		read7=11.
		return
	end if
c
	if(ios.gt.0) then       !error
		write(6,*) 'error while reading'
		write(6,*) 'second record of file 7 header'
		write(6,*) 'nvers =',nvers
		write(6,*) 'ios =',ios
		read7=35.
		return
	else if(ios.lt.0) then  !eof
		write(6,*) 'EOF encountered while reading'
		write(6,*) 'second record of file 7 header'
		write(6,*) 'nvers =',nvers
		read7=31.
		return
	else if(knausm.lt.0.or.knausm.gt.ndim) then     !knausm
		write(6,*) 'error reading knausm : ',knausm
		write(6,*) 'maximum dimension is : ',ndim
		read7=41.
		return
	end if
c
	read7=0.
c
	return
	end
c
c*********************************************************
c
	function writ7(iunit,ndim,nvers,knausm,knaus,hdep
     +                          ,href,hzmin,descrp)
c
c writes first record of file 7
c
c error codes 11
c ndim is dummy argument
c
	character*80 descrp
	integer knaus(ndim)
	real hdep(ndim)
c
	idummy=ndim        
	idummy=2*idummy
c
	nvermx = 6
c
	rewind iunit
c
	if(nvers.eq.1) then
		write(iunit)             knausm
     +                                  ,(knaus(j),j=1,knausm)
	else if(nvers.ge.2.and.nvers.le.nvermx) then
		write(iunit)             nvers
		write(iunit)             knausm
     +                                  ,(knaus(j),j=1,knausm)
     +                                  ,(hdep(j),j=1,knausm)
     +                                  ,href
     +                                  ,hzmin
     +                                  ,descrp
	else
		write(6,*) 'version not recognized : ',nvers
		writ7=11.
		return
	end if
c
	writ7=0.
c
	return
	end
c
c*********************************************************
c
	function rdrc7(iunit,nvers,it,knausm,xv)
c
c reads data record of file 7
c
c error codes 11 35
c EOF -1
c
	real xv(3*knausm)
c
	nvermx = 6
c
	if(nvers.ge.1.and.nvers.le.2) then
		read(iunit,iostat=ios)  tt
     +                                  ,(xv(i),i=1,3*knausm)
		it=iround(tt)
	else if(nvers.ge.3.and.nvers.le.nvermx) then
		read(iunit,iostat=ios)  it
     +                                  ,(xv(i),i=1,3*knausm)
	else
		write(6,*) 'version not recognized : ',nvers
		rdrc7=11.
		return
	end if

	!write(6,*) nvers,it,ios

	if(ios.gt.0) then       !error
		write(6,*) 'error while reading'
		write(6,*) 'data record of file 7'
		rdrc7=35.
		return
	else if(ios.lt.0) then  !eof
		rdrc7=-1.
		return
	end if
c
	rdrc7=0.
c
	return
	end
c
c************************************************************
c
	function skrc7(iunit,nvers,it,knausm,xv)
c
c skips one data record of file 7
c
c error codes 11 35
c EOF -1
c
	real xv(3*knausm)
c
	idummy=knausm        
	dummy=xv(1)
	dummy=idummy*dummy
c
	nvermx = 6
c
	if(nvers.ge.1.and.nvers.le.2) then
		read(iunit,iostat=ios)  tt
		it=iround(tt)
	else if(nvers.ge.3.and.nvers.le.nvermx) then
		read(iunit,iostat=ios)  it
	else
		write(6,*) 'version not recognized : ',nvers
		skrc7=11.
		return
	end if
c
	if(ios.gt.0) then       !error
		write(6,*) 'error while skipping'
		write(6,*) 'data record of file 7'
		skrc7=35.
		return
	else if(ios.lt.0) then  !eof
		skrc7=-1.
		return
	end if
c
	skrc7=0.
c
	return
	end
c
c************************************************************
c
	function wrrc7(iunit,nvers,it,knausm,knaus,xv)
c
c writes data record of float tracking file
c
c error codes 11
c
	integer knaus(knausm)
	real xv(*)
c
	nvermx = 6
c
	if(nvers.ge.1.and.nvers.le.2) then
		write(iunit)    float(it)
     +                          ,((xv(3*(knaus(j)-1)+i)
     +                          ,j=1,knausm)
     +                          ,i=1,3)
	else if(nvers.ge.3.and.nvers.le.nvermx) then
		write(iunit)    it
     +                          ,((xv(3*(knaus(j)-1)+i)
     +                          ,j=1,knausm)
     +                          ,i=1,3)
	else
		write(6,*) 'version not recognized : ',nvers
		wrrc7=11.
		return
	end if
c
	wrrc7=0.
c
	return
	end
c
c************************************************************

        function wrrc77(iunit,nvers,it,knausm,knaus,u,v,z)
c
c writes data record of extra point file
c
c error codes 11
c
        dimension knaus(*),u(*),v(*),z(knausm)
c
                write(iunit)    it
     +                          ,( u(knaus(j)),j=1,knausm )
     +                          ,( v(knaus(j)),j=1,knausm )
     +                          ,( z(knaus(j)),j=1,knausm )
c
        wrrc77=0.
c
        return
        end

c************************************************************
c************************************************************
c************************************************************
c************************************************************
c************************************************************
c************************************************************

	subroutine iniext

c sets up initial common block

	implicit none

c parameters
	integer ftype,maxvers
	parameter(ftype=71,maxvers=7)
c common
	integer mtype,maxver,nverso,npoext
	common /extcom/ mtype,maxver,nverso,npoext
c save
	logical binit
	save binit
	save /extcom/
c data
	data binit /.false./

	if( binit ) return

	binit = .true.

	mtype = ftype
	maxver = maxvers
	nverso = 0
	npoext = 0

	return
	end

c************************************************************

	subroutine rfext	(iunit,nvmax,nvers
     +				,npoint
     +				,href,hzmin
     +				,title
     +				,ierr
     +				)

c reads first record of EXT file

	implicit none

c arguments
	integer iunit,nvmax,nvers
	integer npoint
	real href,hzmin
	character*80 title
	integer ierr
c common
	integer mtype,maxver,nverso,npoext
	common /extcom/ mtype,maxver,nverso,npoext
c local
	integer ntype,irec

c initialize

	call iniext

c control newest version number for call

	if(maxver.ne.nvmax) goto 95

c rewind file

	rewind(iunit,err=96)

c first record - find out what version

	irec = 1
	read(iunit,err=99) ntype,nvers

c control version number and type of file

	if(ntype.ne.mtype) goto 97
	if(nvers.le.0.or.nvers.gt.maxver) goto 98

	if(nvers.lt.7) goto 91	!only type 7 or up

c next records

	irec = 2
	read(iunit,err=99)	 npoint
	read(iunit,err=99)	 href,hzmin
	read(iunit,err=99)	 title

	nverso=nvers
	npoext=npoint

	ierr=0

	return
   99	continue
	write(6,*) 'rfext: Error encountered while'
	write(6,*) 'reading record number ',irec
	write(6,*) 'of EXT file header'
	write(6,*) 'nvers = ',nvers
	ierr=99
	return
   98	continue
	write(6,*) 'rfext: Version not recognized : ',nvers
	ierr=98
	return
   97	continue
	write(6,*) 'rfext: Wrong type of file : ',ntype
	write(6,*) 'Expected ',mtype
	ierr=97
	return
   96	continue
	write(6,*) 'rfext: Cannot rewind file for unit : ',iunit
	ierr=96
	return
   95	continue
	write(6,*) 'rfext: Old function call ',nvmax
	write(6,*) 'Please adjust call to rfext and recompile'
	ierr=95
   91	continue
	write(6,*) 'rfext: Cannot read version ',nvers
	write(6,*) 'Convert to new version with EXTCONV'
	ierr=91
	return
	end

c************************************************************

	subroutine wfext	(iunit,nvmax,nvers
     +				,npoint
     +				,href,hzmin
     +				,title
     +				,ierr
     +				)

c writes first record of EXT file

	implicit none

c arguments
	integer iunit,nvmax,nvers
	integer npoint
	real href,hzmin
	character*80 title
	integer ierr
c common
	integer mtype,maxver,nverso,npoext
	common /extcom/ mtype,maxver,nverso,npoext

c initialize

	call iniext

c control newest version number for call

	if( nvmax.ne.maxver ) goto 95
	if( nvers.ne.maxver .and. nvers.ne.0 ) goto 98

	rewind(iunit)

	write(iunit)		mtype,maxver
	write(iunit)		npoint
	write(iunit)		href,hzmin
	write(iunit)		title

	nverso=maxver
	npoext=npoint

	ierr=0

	return
   98	continue
	write(6,*) 'wfext: Cannot write version ',nvers
	ierr=98
	return
   95	continue
	write(6,*) 'wfext: Old function call ',nvmax
	write(6,*) 'Please adjust call to wfext and recompile'
	ierr=95
	return
	end

c************************************************************

	subroutine rsext(iunit,kpoint,ipoint,hdep,x,y,ierr)

c reads second record of EXT file

	implicit none

c arguments
	integer iunit
	integer kpoint(npoext)
	integer ipoint(npoext)
	real hdep(npoext)
	real x(npoext),y(npoext)
	integer ierr
c common
	integer mtype,maxver,nverso,npoext
	common /extcom/ mtype,maxver,nverso,npoext
c local
	integer i,npoint

	npoint = npoext

	read(iunit,err=99) (kpoint(i),i=1,npoint)
	read(iunit,err=99) (ipoint(i),i=1,npoint)
	read(iunit,err=99) (hdep(i),i=1,npoint)
	read(iunit,err=99) (x(i),i=1,npoint)
	read(iunit,err=99) (y(i),i=1,npoint)

	ierr = 0

	return
   99	continue
	write(6,*) 'rsext: Error encountered while'
	write(6,*) 'reading second part of EXT file header'
	ierr=99
	return
	end

c************************************************************

	subroutine wsext(iunit,kpoint,ipoint,hdep,x,y,ierr)

c writes second record of EXT file

	implicit none

c arguments
	integer iunit
	integer kpoint(npoext)
	integer ipoint(npoext)
	real hdep(npoext)
	real x(npoext),y(npoext)
	integer ierr
c common
	integer mtype,maxver,nverso,npoext
	common /extcom/ mtype,maxver,nverso,npoext
c local
	integer i,npoint

	npoint = npoext

	write(iunit) (kpoint(i),i=1,npoint)
	write(iunit) (ipoint(i),i=1,npoint)
	write(iunit) (hdep(i),i=1,npoint)
	write(iunit) (x(i),i=1,npoint)
	write(iunit) (y(i),i=1,npoint)

	ierr = 0

	return
	end

c************************************************************

	subroutine rdext(iunit,it,u,v,z,ierr)

c reads data record of EXT file

	implicit none

c arguments
	integer iunit,it
	real u(npoext),v(npoext),z(npoext)
	integer ierr
c common
	integer mtype,maxver,nverso,npoext
	common /extcom/ mtype,maxver,nverso,npoext
c local
	integer i,npoint

	npoint = npoext

c time record

	read(iunit,end=88,err=98) it

c data record

	read(iunit,err=99) (u(i),i=1,npoint)
	read(iunit,err=99) (v(i),i=1,npoint)
	read(iunit,err=99) (z(i),i=1,npoint)

	ierr=0

	return
   88	continue
	ierr=-1
	return
   98	continue
	write(6,*) 'rdext: Error while reading'
	write(6,*) 'time record of EXT file'
	ierr=98
	return
   99	continue
	write(6,*) 'rdext: Error while reading'
	write(6,*) 'data record of EXT file'
	write(6,*) 'it = ',it
	ierr=99
	return
	end

c************************************************************

	subroutine wrext(iunit,it,u,v,z,ierr)

c writes data record of EXT file

	implicit none

c arguments
	integer iunit,it
	real u(npoext),v(npoext),z(npoext)
	integer ierr
c common
	integer mtype,maxver,nverso,npoext
	common /extcom/ mtype,maxver,nverso,npoext
c local
	integer i,npoint

	npoint = npoext

	write(iunit) it

	write(iunit) (u(i),i=1,npoint)
	write(iunit) (v(i),i=1,npoint)
	write(iunit) (z(i),i=1,npoint)

	ierr=0

	return
	end

c************************************************************

