!
! $Id: subext.f,v 1.7 2001/11/16 07:35:43 georg Exp $
!
! utility routines to read/write EXT file - file type 71
!
! contents :
!
! subroutine iniext
! subroutine rfext(iunit,nvmax,nvers,npoint,href,hzmin,title,ierr)
! subroutine wfext(iunit,nvmax,nvers,npoint,href,hzmin,title,ierr)
! subroutine rsext(iunit,kpoint,ipoint,hdep,x,y,ierr)
! subroutine wsext(iunit,kpoint,ipoint,hdep,x,y,ierr)
! subroutine rdext(iunit,it,u,v,z,ierr)
! subroutine wrext(iunit,it,u,v,z,ierr)
!
! revision log :
!
! 20.05.1998	ggu	cleaned up a bit
! 04.02.2000	ggu	wrrc77 from newpr to here
! 14.09.2015	ggu	some more helper routines
! 05.10.2015	ggu	handle error in backspace smoothly
!
! notes :
!
! variables used:
!
!	mtype,ftype	type of file (71)
!	maxver,maxvers	newest version
!
!	ierr		error code (0: ok,  <0: EOF,  >0: error)
!
!	iunit		file unit
!	nvmax		version of calling routine (must be maxver)
!	nvers		version of file
!	npoint		number of nodes written
!	href		reference level for water level
!	hzmin		minimum total water depth
!	title		title of run
!
!	kpoint(npoint)	external number of nodes
!	ipoint(npoint)	internal number of nodes
!	hdep(npoint)	depth at nodes
!	x(npoint)	x coordinate of nodes
!	y(npoint)	y coordinate of nodes
!
!	it		actual time of data record
!	u(npoint)	velocity in x direction
!	v(npoint)	velocity in y direction
!	z(npoint)	water level
!
! format of file:
!
! version 7
!
!	mtype,nvers
!	npoint
!	href,hzmin
!	title
!	(kpoint(i),i=1,npoint)
!	(ipoint(i),i=1,npoint)
!	(hdep(i),i=1,npoint)
!	(x(i),i=1,npoint)
!	(y(i),i=1,npoint)
!
!	it
!	(u(i),i=1,npoint)
!	(v(i),i=1,npoint)
!	(z(i),i=1,npoint)
!
! version 3-6
!
! 	nvers
!	npoint,(ipoint(i),i=1,npoint),(hdep(i),i=1,npoint),href,hzmin,title
!
!	it,(u(i),i=1,npoint),(v(i),i=1,npoint),(z(i),i=1,npoint)
!
! version 2
!
! 	nvers
!	npoint,(ipoint(i),i=1,npoint),(hdep(i),i=1,npoint),href,hzmin,title
!
!	float(it),(u(i),i=1,npoint),(v(i),i=1,npoint),(z(i),i=1,npoint)
!
! version 1
!
!	npoint,(ipoint(i),i=1,npoint)
!
!	float(it),(u(i),i=1,npoint),(v(i),i=1,npoint),(z(i),i=1,npoint)
!
!************************************************************
!
! old utility routines to read/write EXT file
!
! contents
!
!-----------------------------------------------------------------------
! function read7(iunit,ndim,nvers,knausm,knaus,hdep,href,hzmin,descrp)
!                       reads first record of file 7
! function writ7(iunit,ndim,nvers,knausm,knaus,hdep,href,hzmin,descrp)
!                       writes first record of file 7
! function rdrc7(iunit,nvers,it,knausm,xv)
!                       reads data record of file 7
! function skrc7(iunit,nvers,it,knausm,xv)
!                       skips one data record of file 7
! function wrrc7(iunit,nvers,it,knausm,knaus,xv)
!                       writes data record of file 7
! function wrrc77(iunit,nvers,it,knausm,knaus,u,v,z)
!			writes data record of extra point file
!-----------------------------------------------------------------------
!
! nvermx		maximum version recognized -> 6
!
!*********************************************************
!*********************************************************
!*********************************************************
!*********************************************************
!*********************************************************
!*********************************************************
!-----------------------------------------------------------------------
        module ioext
!-----------------------------------------------------------------------
        contains
!-----------------------------------------------------------------------

	function check_ext_file(file)

        use fil

	implicit none

        logical check_ext_file
        character*(*) file

        integer nb,nvers,knausm,ierr

        check_ext_file = .false.

        nb = ifileo(0,file,'unform','old')
        if( nb .le. 0 ) return
	call ext_check_header(nb,nvers,knausm,ierr)
        close(nb)

        check_ext_file = ierr == 0

	end

!*********************************************************

	subroutine ext_is_ext_file(iunit,nvers)

	implicit none

	integer iunit,nvers

	integer knausm,ierr

	call ext_check_header(iunit,nvers,knausm,ierr)

	if( ierr .ne. 0 ) nvers = 0

	rewind(iunit)

	end

!*********************************************************
!*********************************************************
!*********************************************************

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

!*********************************************************

	subroutine ext_peek_record(iunit,nvers,it,ierr)

	implicit none

	integer iunit,nvers,it,ierr

	integer, parameter :: nvermx = 6
	double precision tt

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

!*********************************************************

	subroutine ext_check_header(iunit,nvers,knausm,ierr)

! checks version of ext file and returns number of points

	implicit none

	integer iunit,nvers,knausm,ierr

	integer, parameter :: nvermx = 6
	integer kaux,ios,i,it,j
	double precision haux,tt,xaux
	character*80 title

	nvers = 0
	knausm = 0
	ierr = -1

	read(iunit,iostat=ios) nvers
	if( ios .ne. 0 ) return

	if(nvers.ge.2.and.nvers.le.nvermx) then
	  read(iunit,iostat=ios) knausm,(kaux,j=1,knausm),(haux,j=1,knausm),haux,haux,title
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

!*********************************************************
!*********************************************************
!*********************************************************

	subroutine ext_read_header(iunit,ndim,nvers,knausm,knaus,hdep,href,hzmin,descrp)
	implicit none
	integer iunit,ndim,nvers,knausm
	double precision href,hzmin
	character*80 descrp
	integer knaus(knausm)
	double precision hdep(knausm)

	integer ierr

	ierr = read7(iunit,ndim,nvers,knausm,knaus,hdep,href,hzmin,descrp)

	if( ierr .ne. 0. ) then
	  stop 'error stop read_ext_header: read error'
	end if

	end

!*********************************************************

	subroutine ext_read_record(iunit,nvers,it,knausm,xv,ierr)
	implicit none
	integer iunit,nvers,it,knausm,ierr
	double precision xv(3*knausm)

	ierr = rdrc7(iunit,nvers,it,knausm,xv)

	end

!*********************************************************

	subroutine ext_write_header(iunit,ndim,nvers,knausm,knaus,hdep,href,hzmin,descrp)
	implicit none
	integer iunit,ndim,nvers,knausm
	double precision href,hzmin
	character*80 descrp
	integer knaus(knausm)
	double precision hdep(knausm)

	integer ierr

	ierr = writ7(iunit,ndim,nvers,knausm,knaus,hdep,href,hzmin,descrp)

	end

!*********************************************************

	subroutine ext_write_record(iunit,nvers,it,knausm,xv)
	implicit none
	integer iunit,nvers,it,knausm
	double precision xv(3*knausm)

	integer ierr,i
	integer knaus(knausm)

	do i=1,knausm
	  knaus(i) = i
	end do

	ierr = wrrc7(iunit,nvers,it,knausm,knaus,xv)

	end

!*********************************************************
!*********************************************************
!*********************************************************

	function read7(iunit,ndim,nvers,knausm,knaus,hdep,href,hzmin,descrp)
!
! reads first record of file 7
!
! error codes 11 21 31 35 41 61 71
!
	character*80 descrp
	integer knaus(ndim)
	double precision hdep(ndim)
        double precision read7
        double precision href, hzmin
!
	nvermx = 6
!
	rewind(iunit,iostat=ios)
!
	if(ios.ne.0) then
		write(6,*) 'Cannot rewind file for unit :',iunit
		read7=71.
		return
	end if
!
! first record
!
	read(iunit,iostat=ios) nvers
!
	if(ios.eq.0) then       ! no read error ==> first version
!
	else if(ios.lt.0) then  !eof
		write(6,*) 'EOF encountered while reading'
		write(6,*) 'first record of file 7 header'
		write(6,*) 'nvers =',nvers
		read7=21.
		return
	end if
!
! second record
!
	if(nvers.eq.1) then
		hzmin=0.05
		descrp=' '
		do j=1,knausm
		hdep(j)=100000.         !no dry areas
		end do
	else if(nvers.ge.2.and.nvers.le.nvermx) then
	  read(iunit,iostat=ios) knausm,(knaus(j),j=1,knausm),(hdep(j),j=1,knausm)      &
     &                                  ,href,hzmin,descrp
	else
		write(6,*) 'version not recognized : ',nvers
		read7=11.
		return
	end if
!
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
!
	read7=0.
!
	return
	end
!
!*********************************************************
!
	function writ7(iunit,ndim,nvers,knausm,knaus,hdep,href,hzmin,descrp)
!
! writes first record of file 7
!
! error codes 11
! ndim is dummy argument
!
	character*80 descrp
	integer knaus(ndim)
	double precision hdep(ndim)
	double precision writ7
        double precision href, hzmin
!
	idummy=ndim        
	idummy=2*idummy
!
	nvermx = 6
!
	rewind iunit
!
	if(nvers.eq.1) then
		write(iunit) knausm,(knaus(j),j=1,knausm)
	else if(nvers.ge.2.and.nvers.le.nvermx) then
		write(iunit)             nvers
		write(iunit) knausm,(knaus(j),j=1,knausm),(hdep(j),j=1,knausm)  &
     &                                  ,href,hzmin,descrp
	else
		write(6,*) 'version not recognized : ',nvers
		writ7=11.
		return
	end if
!
	writ7=0.
!
	return
	end
!
!*********************************************************
!
	function rdrc7(iunit,nvers,it,knausm,xv)
!
! reads data record of file 7
!
! error codes 11 35
! EOF -1
!
        use utility

	double precision xv(3*knausm)
	double precision rdrc7,tt
!
	nvermx = 6
!
	if(nvers.ge.1.and.nvers.le.2) then
		read(iunit,iostat=ios)  tt,(xv(i),i=1,3*knausm)
		it=iround(tt)
	else if(nvers.ge.3.and.nvers.le.nvermx) then
		read(iunit,iostat=ios)  it,(xv(i),i=1,3*knausm)
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
!
	rdrc7=0.
!
	return
	end
!
!************************************************************
!
	function skrc7(iunit,nvers,it,knausm,xv)
!
! skips one data record of file 7
!
! error codes 11 35
! EOF -1
!
        use utility
	double precision xv(3*knausm)
	double precision tt
!
	idummy=knausm        
	dummy=xv(1)
	dummy=idummy*dummy
!
	nvermx = 6
!
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
!
	if(ios.gt.0) then       !error
		write(6,*) 'error while skipping'
		write(6,*) 'data record of file 7'
		skrc7=35.
		return
	else if(ios.lt.0) then  !eof
		skrc7=-1.
		return
	end if
!
	skrc7=0.
!
	return
	end
!
!************************************************************
!
	function wrrc7(iunit,nvers,it,knausm,knaus,xv)
!
! writes data record of float tracking file
!
! error codes 11
!
	integer knaus(knausm)
	double precision xv(*)
	double precision wrrc7
!
	nvermx = 6
!
	if(nvers.ge.1.and.nvers.le.2) then
	  write(iunit) float(it),((xv(3*(knaus(j)-1)+i),j=1,knausm),i=1,3)
	else if(nvers.ge.3.and.nvers.le.nvermx) then
	  write(iunit) it,((xv(3*(knaus(j)-1)+i),j=1,knausm),i=1,3)
	else
	  write(6,*) 'version not recognized : ',nvers
	  wrrc7=11.
	  return
	end if
!
	wrrc7=0.
!
	return
	end
!
!************************************************************

        function wrrc77(iunit,nvers,it,knausm,knaus,u,v,z)
!
! writes data record of extra point file
!
! error codes 11
!
        dimension knaus(knausm)
        double precision, dimension(:) :: u(knausm),v(knausm),z(knausm)
	double precision wrrc77
!
          write(iunit) it,( u(j),j=1,knausm ),( v(j),j=1,knausm ),( z(j),j=1,knausm )
!
        wrrc77=0.
!
        return
        end

!************************************************************
!************************************************************
!************************************************************
!************************************************************
!************************************************************
!************************************************************

	subroutine iniext

! sets up initial common block

	implicit none

! parameters
	integer ftype,maxvers
	parameter(ftype=71,maxvers=7)
! common
	integer mtype,maxver,nverso,npoext
	common /extcom/ mtype,maxver,nverso,npoext
! save
	logical binit
	save binit
	save /extcom/
! data
	data binit /.false./

	if( binit ) return

	binit = .true.

	mtype = ftype
	maxver = maxvers
	nverso = 0
	npoext = 0

	return
	end

!************************************************************

	subroutine rfext(iunit,nvmax,nvers,npoint,href,hzmin,title,ierr)

! reads first record of EXT file

	implicit none

! arguments
	integer iunit,nvmax,nvers
	integer npoint
	double precision href,hzmin
	character*80 title
	integer ierr
! common
	integer mtype,maxver,nverso,npoext
	common /extcom/ mtype,maxver,nverso,npoext
! local
	integer ntype,irec

! initialize

	call iniext

! control newest version number for call

	if(maxver.ne.nvmax) goto 95

! rewind file

	rewind(iunit,err=96)

! first record - find out what version

	irec = 1
	read(iunit,err=99) ntype,nvers

! control version number and type of file

	if(ntype.ne.mtype) goto 97
	if(nvers.le.0.or.nvers.gt.maxver) goto 98

	if(nvers.lt.7) goto 91	!only type 7 or up

! next records

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

!************************************************************

	subroutine wfext(iunit,nvmax,nvers,npoint,href,hzmin,title,ierr)

! writes first record of EXT file

	implicit none

! arguments
	integer iunit,nvmax,nvers
	integer npoint
	double precision href,hzmin
	character*80 title
	integer ierr
! common
	integer mtype,maxver,nverso,npoext
	common /extcom/ mtype,maxver,nverso,npoext

! initialize

	call iniext

! control newest version number for call

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

!************************************************************

	subroutine rsext(iunit,kpoint,ipoint,hdep,x,y,ierr)

! reads second record of EXT file

	implicit none

! arguments
	integer iunit
	integer kpoint(npoext)
	integer ipoint(npoext)
	double precision hdep(npoext)
	double precision x(npoext),y(npoext)
	integer ierr
! common
	integer mtype,maxver,nverso,npoext
	common /extcom/ mtype,maxver,nverso,npoext
! local
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

!************************************************************

	subroutine wsext(iunit,kpoint,ipoint,hdep,x,y,ierr)

! writes second record of EXT file

	implicit none

! arguments
	integer iunit
	integer kpoint(npoext)
	integer ipoint(npoext)
	double precision hdep(npoext)
	double precision x(npoext),y(npoext)
	integer ierr
! common
	integer mtype,maxver,nverso,npoext
	common /extcom/ mtype,maxver,nverso,npoext
! local
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

!************************************************************

	subroutine rdext(iunit,it,u,v,z,ierr)

! reads data record of EXT file

	implicit none

! arguments
	integer iunit,it
	double precision u(npoext),v(npoext),z(npoext)
	integer ierr
! common
	integer mtype,maxver,nverso,npoext
	common /extcom/ mtype,maxver,nverso,npoext
! local
	integer i,npoint

	npoint = npoext

! time record

	read(iunit,end=88,err=98) it

! data record

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

!************************************************************

	subroutine wrext(iunit,it,u,v,z,ierr)

! writes data record of EXT file

	implicit none

! arguments
	integer iunit,it
	double precision u(npoext),v(npoext),z(npoext)
	integer ierr
! common
	integer mtype,maxver,nverso,npoext
	common /extcom/ mtype,maxver,nverso,npoext
! local
	integer i,npoint

	npoint = npoext

	write(iunit) it

	write(iunit) (u(i),i=1,npoint)
	write(iunit) (v(i),i=1,npoint)
	write(iunit) (z(i),i=1,npoint)

	ierr=0

	return
	end

!************************************************************

!-----------------------------------------------------------------------
        end module ioext
!-----------------------------------------------------------------------
