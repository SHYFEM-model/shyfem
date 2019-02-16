
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

c utility routines to read/write EXT file - file type 71
c
c contents :
c
c revision log :
c
c 20.05.1998	ggu	cleaned up a bit
c 04.02.2000	ggu	wrrc77 from newpr to here
c 14.09.2015	ggu	some more helper routines
c 05.10.2015	ggu	handle error in backspace smoothly
c 05.10.2017	ggu	file 7 substituted with EXT file in output
c 20.10.2017	ggu	completely restructured for version 7
c 30.08.2018	ggu	new version 8 to avoid implicit do
c 09.11.2018	ggu	more general version of linear routines
c
c notes :
c
c variables used:
c
c	ierr		error code (0: ok,  <0: EOF,  >0: error)
c	iunit		file unit
c
c	mtype		type of file (id)
c	nvers		version of file
c	knausm		number of nodes written
c	lmax		vertical dimension
c	nvar		number of rvariable records (+1)
c
c	atime0		reference time (absolute)
c	href		reference level for water level
c	hzmin		minimum total water depth
c	title		title of run
c	femver		version of shyfem
c
c	knaus(knausm)	external number of nodes
c	hdep(knausm)	depth at nodes
c	ilhkv(knausm)	vertical levels of node
c	x(knausm)	x coordinate of nodes
c	y(knausm)	y coordinate of nodes
c	strings(knausm)	description of node
c	hlv(lmax)	level depth
c
c	atime		time of record (absolute)
c	ivar		indicator of variable (0 for barotropic u,v,z)
c	m		number of variables 
c				(3 for ivar=0, 2 for ivar=2, else 1)
c	lm		vertical levels contained in record
c				(lm=1 for ivar=0, else lm=lmax)
c	vals(lmax,knausm,3)	value of variables
c
c	u,v,z		for record 0 barotropic velocity and water level
c
c format of file:
c
c version 8
c
c	as version 7, but records are written as
c
c	atime,ivar,m,lm,nlin,vals(1:nlin)
c
c version 7
c
c	mtype,nvers
c	knausm,lmax,nvar
c
c	atime0
c	href,hzmin
c	title,femver
c	knaus,hdep,ilhkv,x,y,strings
c	hlv
c
c	atime,ivar,m,lm,vals
c
c version 3-6
c
c 	nvers
c	knausm,(knaus(i),i=1,knausm),(hdep(i),i=1,knausm),href,hzmin,title
c
c	it,(u(i),i=1,knausm),(v(i),i=1,knausm),(z(i),i=1,knausm)
c
c version 2
c
c 	nvers
c	knausm,(knaus(i),i=1,knausm),(hdep(i),i=1,knausm),href,hzmin,title
c
c	float(it),(u(i),i=1,knausm),(v(i),i=1,knausm),(z(i),i=1,knausm)
c
c version 1
c
c	knausm,(knaus(i),i=1,knausm)
c
c	float(it),(u(i),i=1,knausm),(v(i),i=1,knausm),(z(i),i=1,knausm)
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
        integer, save :: ext_maxvers = 8

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
	else if(nvers.ge.7.and.nvers.le.8) then
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
     +                          ,href,hzmin,title,femver
     +                          ,knaus,hdep,ilhkv,x,y,strings,hlv
     +				,ierr)

	implicit none

	integer iunit,nvers,knausm,lmax
	double precision atime0
	real href,hzmin
	character*80 title,femver
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
     +                          ,href,hzmin,title)
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
	  read(iunit,iostat=ierr) title,femver
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
        integer nlin
        real, allocatable :: rlin(:)

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
        else if( nvers == 7 ) then
	  read(iunit,iostat=ierr) atime,ivar,m,lm
          if( ierr /= 0 ) return
	  if( m > 3 ) goto 99
          backspace(iunit)
	  call count_linear(lm,knausm,m,ilhkv,nlin)
	  allocate(rlin(nlin))
	  read(iunit,iostat=ierr) atime,ivar,m,lm,(rlin(i),i=1,nlin)
          call linear2vals(lm,knausm,m,ilhkv,vals,rlin,nlin)
!	  read(iunit,iostat=ierr) atime,ivar,m,lm
!     +				,(((vals(l,j,i)
!     +				,l=1,min(lm,ilhkv(j)))
!     +				,j=1,knausm)
!     +				,i=1,m)
        else
	  nlin = knausm*lmax*3
	  allocate(rlin(nlin))
	  read(iunit,iostat=ierr) atime,ivar,m,lm,nlin,(rlin(i),i=1,nlin)
          if( ierr /= 0 ) return
	  if( m > 3 ) goto 99
          call linear2vals(lm,knausm,m,ilhkv,vals,rlin,nlin)
	end if

	return
   99	continue
	write(6,*) 'm = ',m
	stop 'error stop ext_read_record: m > 3'
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
     +                          ,href,hzmin,title,femver
     +                          ,knaus,hdep,ilhkv,x,y,strings,hlv
     +				,ierr)

	implicit none

	integer iunit,nvers,knausm,lmax
	double precision atime0
	real href,hzmin
	character*80 title,femver
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
	write(iunit,iostat=ierr) title,femver
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
	real, intent(in) :: vals(lmax,knausm,m)
	integer, intent(out) :: ierr

	integer i,j,l,lm
        integer nlin
        real, allocatable :: rlin(:)

	if( ivar == 0 ) then
	  lm = 1
	  write(iunit,iostat=ierr) atime,ivar,m,lm
     +				,((vals(1,j,i)
     +				,j=1,knausm)
     +				,i=1,m)
	else
!	  write(iunit,iostat=ierr) atime,ivar,m,lmax
!     +				,(((vals(l,j,i)
!     +				,l=1,min(lmax,ilhkv(j)))
!     +				,j=1,knausm)
!     +				,i=1,m)
	  nlin = lmax*knausm*m
	  allocate(rlin(nlin))
          call vals2linear(lmax,knausm,m,ilhkv,vals,rlin,nlin)
	  write(iunit,iostat=ierr) atime,ivar,m,lmax,nlin
     +                                  ,(rlin(i),i=1,nlin)
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
     +                          ,href,hzmin,title)
c
c reads first record of file 7
c
c error codes 11 21 31 35 41 61 71
c
	character*80 title
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
		title=' '
		do j=1,knausm
		hdep(j)=100000.         !no dry areas
		end do
	else if(nvers.ge.2.and.nvers.le.nvermx) then
		read(iunit,iostat=ios)   knausm
     +                                  ,(knaus(j),j=1,knausm)
     +                                  ,(hdep(j),j=1,knausm)
     +                                  ,href
     +                                  ,hzmin
     +                                  ,title
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
     +                          ,href,hzmin,title)
c
c writes first record of file 7
c
c error codes 11
c ndim is dummy argument
c
	character*80 title
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
     +                                  ,title
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

