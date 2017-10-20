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
c 05.10.2017	ggu	file 7 substituted with EXT file in output
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
c	knausm,lmax
c	atime0
c	href,hzmin
c	title,femver
c	(knaus(i),i=1,npoint)
c	(hdep(i),i=1,npoint)
c	(ilhkv(i),i=1,npoint)
c	(x(i),i=1,npoint)
c	(y(i),i=1,npoint)
c	(strings(i),i=1,npoint)
c	(hlv(i),i=1,lmax)
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

!==================================================================
        module extfile
!==================================================================

        implicit none

        integer, save :: ext_type = 947336
        integer, save :: ext_maxvers = 7

!==================================================================
        contains
!==================================================================

!==================================================================
        end module extfile
!==================================================================

	function check_ext_file(file)

	implicit none

        logical check_ext_file
        character*(*) file

        integer nb,nvers,knausm,lmax,nvar,ierr
        integer ifileo

        check_ext_file = .false.

        nb = ifileo(0,file,'unform','old')
        if( nb .le. 0 ) return
	call ext_check_header(nb,nvers,knausm,lmax,nvar,ierr)
        close(nb)

        check_ext_file = ( ierr == 0 )

	end

c*********************************************************

	subroutine ext_is_ext_file(iunit,nvers)

	implicit none

	integer iunit,nvers

	integer knausm,lmax,nvar,ierr

	call ext_check_header(iunit,nvers,knausm,lmax,nvar,ierr)

	if( ierr .ne. 0 ) nvers = 0

	end

c*********************************************************
c*********************************************************
c*********************************************************

	subroutine ext_peek_header(iunit,nvers,knausm,lmax,nvar,ierr)

	implicit none

	integer iunit,nvers,knausm,lmax,nvar,ierr

	call ext_check_header(iunit,nvers,knausm,lmax,nvar,ierr)

	end

c*********************************************************

	subroutine ext_peek_record(iunit,nvers,atime,ivar,ierr)

	implicit none

	integer iunit,nvers,ivar,ierr
	double precision atime

	integer it
	real tt

	atime = 0
	ivar = 0

	if(nvers.ge.1.and.nvers.le.2) then
	  read(iunit,iostat=ierr)  tt
	  atime = tt
	else if(nvers.ge.3.and.nvers.le.6) then
	  read(iunit,iostat=ierr)  it
	  atime = it
	else if(nvers.ge.7.and.nvers.le.7) then
	  read(iunit,iostat=ierr)  atime,ivar
	else
	  stop 'error stop ext_peek_record: internal error (1)'
	end if

	backspace(iunit)

	end

c*********************************************************

	subroutine ext_check_header(iunit,nvers,knausm,lmax,nvar,ierr)

c checks version of ext file and returns number of points

	implicit none

	integer iunit,nvers,knausm,lmax,nvar,ierr

	call ext_check_new_header(iunit,nvers,knausm,lmax,nvar,ierr)
	if( ierr == 0 ) return
	call ext_check_old_header(iunit,nvers,knausm,lmax,nvar,ierr)

	end

c*********************************************************

	subroutine ext_check_new_header(iunit,nvers,knausm,lmax,nvar,ierr)

c checks version of ext file and returns number of points ( nvers > 6 )

	use extfile

	implicit none

	integer iunit,nvers,knausm,lmax,nvar,ierr

	integer ios,ntype

	nvers = 0
	knausm = 0
	ierr = -1

	read(iunit,iostat=ios) ntype,nvers
	if( ios /= 0 ) goto 99
	if( ntype /= ext_type ) goto 99
	if( nvers <= 6 ) goto 98

	read(iunit,iostat=ios) knausm,lmax,nvar
	if( ios /= 0 ) goto 99

	ierr = 0

   99	continue
	rewind(iunit)		!we rewind in any case

	return
   98	continue
	write(6,*) 'cannot read old version with new routine...'
	stop 'error stop ext_check_new_header: internal error (1)'
	end

c*********************************************************

	subroutine ext_check_old_header(iunit,nvers,knausm,lmax,nvar,ierr)

c checks version of ext file and returns number of points ( nvers <= 6 )

	implicit none

	integer iunit,nvers,knausm,lmax,nvar,ierr

	integer kaux,ios,i,it,j
	real haux,tt,xaux
	character*80 title

	nvers = 0
	knausm = 0
	lmax = 1
	nvar = 1
	ierr = -1

	read(iunit,iostat=ios) nvers
	if( ios /= 0 ) goto 99
	if( nvers < 1 .or. nvers > 6 ) goto 99

	if(nvers.ge.2.and.nvers.le.6) then
	  read(iunit,iostat=ios)   knausm
     +                                  ,(kaux,j=1,knausm)
     +                                  ,(haux,j=1,knausm)
     +                                  ,haux
     +                                  ,haux
     +                                  ,title
	else
	  ios = 1
	end if
	if( ios /= 0 ) goto 99

	if(nvers.ge.1.and.nvers.le.2) then
	  read(iunit,iostat=ios)  tt,(xaux,i=1,3*knausm)
	else if(nvers.ge.3.and.nvers.le.6) then
	  read(iunit,iostat=ios)  it,(xaux,i=1,3*knausm)
	end if
	if( ios /= 0 ) goto 99

	ierr = 0

   99	continue
	rewind(iunit)		!we rewind in any case

	return
	end

c*********************************************************
c*********************************************************
c*********************************************************

	subroutine ext_read_header(iunit,nvers,knausm,lmax,nvar,ierr)

	implicit none

	integer iunit,nvers,knausm
	integer lmax,nvar
	integer ierr

	call ext_check_header(iunit,nvers,knausm,lmax,nvar,ierr) !this rewinds
	if( ierr /= 0 ) return

	if( nvers > 6 ) then
	  read(iunit)		!empty read - must succeed
	  read(iunit)		!empty read - must succeed
	end if

	end

c*********************************************************

	subroutine ext_read_header2(iunit,nvers,knausm,lmax
     +                          ,atime0
     +                          ,href,hzmin,descrp,femver
     +                          ,knaus,hdep,ilhkv,x,y,strings,hlv
     +				,ierr)

	implicit none

	integer iunit,nvers,knausm,lmax
	double precision atime0
	real href,hzmin
	character*80 descrp,femver
	integer knaus(knausm)
	real hdep(knausm)
	integer ilhkv(knausm)
	real x(knausm),y(knausm)
	character*80 strings(knausm)
	real hlv(lmax)
	integer ierr

	integer ndim
	real read7

	if( nvers <= 6 ) then
	  ndim = knausm
	  ierr = read7(iunit,ndim,nvers,knausm,knaus,hdep
     +                          ,href,hzmin,descrp)
	  if( ierr /= 0 ) return
	  femver = ' '
	  ilhkv = 1
	  x = 0.
	  y = 0.
	  strings = ' '
	  hlv(1) = 10000.
	else
	  read(iunit,iostat=ierr) atime0
	  if( ierr /= 0 ) return
	  read(iunit,iostat=ierr) href,hzmin
	  if( ierr /= 0 ) return
	  read(iunit,iostat=ierr) descrp,femver
	  if( ierr /= 0 ) return
	  read(iunit,iostat=ierr) knaus,hdep,ilhkv,x,y,strings
	  if( ierr /= 0 ) return
	  read(iunit,iostat=ierr) hlv
	  if( ierr /= 0 ) return
	end if

	end

c*********************************************************

	subroutine ext_read_record(iunit,nvers,atime,knausm,lmax
     +					,ivar,m,ilhkv,vals,ierr)

	implicit none

	integer, intent(in) :: iunit,nvers,knausm,lmax
	integer, intent(in) :: ilhkv(knausm)
	integer, intent(out) :: ivar,m,ierr
	double precision, intent(out) :: atime
	real, intent(out) :: vals(lmax,knausm,3)

	integer i,j,l,it,lm
	real xv(knausm,3)

	real rdrc7

	if( nvers <= 6 ) then
	  ierr = rdrc7(iunit,nvers,it,knausm,xv)
	  if( ierr /= 0 ) return
	  atime = it
	  ivar = 0
	  m = 3
	  vals(1,:,1) = xv(:,1)
	  vals(1,:,2) = xv(:,2)
	  vals(1,:,3) = xv(:,3)
	else
	  read(iunit,iostat=ierr) atime,ivar,m,lm
     +				,(((vals(l,j,i)
     +				,l=1,min(lm,ilhkv(j)))
     +				,j=1,knausm)
     +				,i=1,m)
	end if

	end

c*********************************************************

	subroutine ext_write_header(iunit,nvers,knausm,lmax,nvar,ierr)

	use extfile

	implicit none

	integer iunit,nvers,knausm,lmax,nvar,ierr

	if( nvers /= 0 .and. nvers /= ext_maxvers ) then
	  write(6,*) 'cannot write this version for EXT file: ',nvers
	  write(6,*) 'please either use 0 or ',ext_maxvers
	  ierr = 999
	  return
	end if

	rewind(iunit,iostat=ierr)
	if( ierr /= 0 ) return
	
	write(iunit,iostat=ierr) ext_type,ext_maxvers
	if( ierr /= 0 ) return

	write(iunit,iostat=ierr) knausm,lmax,nvar
	if( ierr /= 0 ) return

	end

c*********************************************************

	subroutine ext_write_header2(iunit,nvers,knausm,lmax
     +                          ,atime0
     +                          ,href,hzmin,descrp,femver
     +                          ,knaus,hdep,ilhkv,x,y,strings,hlv
     +				,ierr)

	implicit none

	integer iunit,nvers,knausm,lmax
	double precision atime0
	real href,hzmin
	character*80 descrp,femver
	integer knaus(knausm)
	real hdep(knausm)
	integer ilhkv(knausm)
	real x(knausm),y(knausm)
	character*80 strings(knausm)
	real hlv(lmax)
	integer ierr

	write(iunit,iostat=ierr) atime0
	if( ierr /= 0 ) return
	write(iunit,iostat=ierr) href,hzmin
	if( ierr /= 0 ) return
	write(iunit,iostat=ierr) descrp,femver
	if( ierr /= 0 ) return
	write(iunit,iostat=ierr) knaus,hdep,ilhkv,x,y,strings
	if( ierr /= 0 ) return
	write(iunit,iostat=ierr) hlv
	if( ierr /= 0 ) return

	end

c*********************************************************

	subroutine ext_write_record(iunit,nvers,atime,knausm,lmax
     +					,ivar,m,ilhkv,vals,ierr)

	implicit none

	integer, intent(in) :: iunit,nvers,knausm,lmax
	integer, intent(in) :: ilhkv(knausm)
	integer, intent(in) :: ivar,m
	double precision, intent(in) :: atime
	real, intent(in) :: vals(lmax,knausm,3)
	integer, intent(out) :: ierr

	integer i,j,l,lm

	if( ivar == 0 ) then
	  lm = 1
	  write(iunit,iostat=ierr) atime,ivar,m,lm
     +				,((vals(1,j,i)
     +				,j=1,knausm)
     +				,i=1,m)
	else
	  write(iunit,iostat=ierr) atime,ivar,m,lmax
     +				,(((vals(l,j,i)
     +				,l=1,min(lmax,ilhkv(j)))
     +				,j=1,knausm)
     +				,i=1,m)
	end if

	end

c*********************************************************
c*********************************************************
c*********************************************************
c old routines - needed to read/write until version 6
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
		write(6,*) 'first record of EXT header'
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
		write(6,*) 'second record of EXT header'
		write(6,*) 'nvers =',nvers
		write(6,*) 'ios =',ios
		read7=35.
		return
	else if(ios.lt.0) then  !eof
		write(6,*) 'EOF encountered while reading'
		write(6,*) 'second record of EXT header'
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
	!write(6,*) 'rdrc7: ',nvers,ios,it,knausm

	!write(6,*) nvers,it,ios

	if(ios.gt.0) then       !error
		write(6,*) 'error while reading'
		write(6,*) 'data record of EXT file'
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
		write(6,*) 'data record of EXT file'
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
c next routines to be deleted (never used or referenced)
c************************************************************
c************************************************************
c************************************************************

	subroutine rfext_000	(iunit,nvmax,nvers
     +				,npoint
     +				,href,hzmin
     +				,title
     +				,ierr
     +				)

c reads first record of EXT file

	use extfile

	implicit none

c arguments
	integer iunit,nvmax,nvers
	integer npoint
	real href,hzmin
	character*80 title
	integer ierr
c local
	integer ntype,irec

c control newest version number for call

	if(ext_maxvers.ne.nvmax) goto 95

c rewind file

	rewind(iunit,err=96)

c first record - find out what version

	irec = 1
	read(iunit,err=99) ntype,nvers

c control version number and type of file

	if(ntype.ne.ext_type) goto 97
	if(nvers.le.0.or.nvers.gt.ext_maxvers) goto 98

	if(nvers.lt.7) goto 91	!only type 7 or up

c next records

	irec = 2
	read(iunit,err=99)	 npoint
	read(iunit,err=99)	 href,hzmin
	read(iunit,err=99)	 title

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
	write(6,*) 'Expected ',ext_type
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

	subroutine wfext_000	(iunit,nvmax,nvers
     +				,npoint
     +				,href,hzmin
     +				,title
     +				,ierr
     +				)

c writes first record of EXT file

	use extfile

	implicit none

c arguments
	integer iunit,nvmax,nvers
	integer npoint
	real href,hzmin
	character*80 title
	integer ierr

c control newest version number for call

	if( nvmax.ne.ext_maxvers ) goto 95
	if( nvers.ne.ext_maxvers .and. nvers.ne.0 ) goto 98

	rewind(iunit)

	write(iunit)		ext_type,ext_maxvers
	write(iunit)		npoint
	write(iunit)		href,hzmin
	write(iunit)		title

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

	subroutine rsext_000(iunit,npoint,kpoint,ipoint,hdep,x,y,ierr)

c reads second record of EXT file

	implicit none

c arguments
	integer iunit
	integer npoint
	integer kpoint(npoint)
	integer ipoint(npoint)
	real hdep(npoint)
	real x(npoint),y(npoint)
	integer ierr
c local
	integer i

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

	subroutine wsext_000(iunit,npoint,kpoint,ipoint,hdep,x,y,ierr)

c writes second record of EXT file

	implicit none

c arguments
	integer iunit
	integer npoint
	integer kpoint(npoint)
	integer ipoint(npoint)
	real hdep(npoint)
	real x(npoint),y(npoint)
	integer ierr
c local
	integer i

	write(iunit) (kpoint(i),i=1,npoint)
	write(iunit) (ipoint(i),i=1,npoint)
	write(iunit) (hdep(i),i=1,npoint)
	write(iunit) (x(i),i=1,npoint)
	write(iunit) (y(i),i=1,npoint)

	ierr = 0

	return
	end

c************************************************************

	subroutine rdext_000(iunit,it,np,u,v,z,ierr)

c reads data record of EXT file

	implicit none

c arguments
	integer iunit,it,np
	real u(np),v(np),z(np)
	integer ierr
c local
	integer i

c time record

	read(iunit,end=88,err=98) it

c data record

	read(iunit,err=99) (u(i),i=1,np)
	read(iunit,err=99) (v(i),i=1,np)
	read(iunit,err=99) (z(i),i=1,np)

	ierr=0

	return
   88	continue
	backspace(iunit)
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

	subroutine wrext_000(iunit,it,np,u,v,z,ierr)

c writes data record of EXT file

	implicit none

c arguments
	integer iunit,it,np
	real u(np),v(np),z(np)
	integer ierr
c local
	integer i

	write(iunit) it

	write(iunit) (u(i),i=1,np)
	write(iunit) (v(i),i=1,np)
	write(iunit) (z(i),i=1,np)

	ierr=0

	return
	end

c************************************************************

