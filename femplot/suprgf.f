
!--------------------------------------------------------------------------
!
!    Copyright (C) 2009  Georg Umgiesser
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

c regular file time handling
c
c revision log :
c
c 10.03.2009	ggu	finished coding
c 27.03.2009	ggu	bug fix ITACT (itact was not adjourned)
c 31.03.2009	ggu	changed header of regular files (x0,.., description)
c 18.12.2018	ggu	changed VERS_7_5_52
c 21.05.2019	ggu	changed VERS_7_5_62
c
c*********************************************************************

	subroutine rgf_admin(it,file,nvar0,ndim,ifile,dfile)

c administers interpolation of regular fields

c can be called any time
c at first call ifile(1) must be -1 (to indicate need for initialization)
c ifile and dfile are returned and contain all needed information

	implicit none

	integer it		!time for interpolation
	character*(*) file	!name of file (only needed for first call)
	integer nvar0		!how many variables are expected
	integer ndim		!dimension of dfile
	integer ifile(7)	!info on regular field data
	real dfile(ndim,3)	!data of regular field
	
	integer iunit
	integer itact,itold,itnew
	integer i
	integer nvar
	integer nx,ny

	integer ifileo

	itact = it

c-------------------------------------------------------------
c initialize regular data files
c-------------------------------------------------------------

	iunit = ifile(1)
	if( iunit .eq. 0 ) return	!file has encountered EOF before

	if( iunit .lt. 0 ) then

	  do i=1,7
	    ifile(i) = 0
	  end do

	  if( file .eq. ' ' ) return	!no file given

	  write(6,*) 'initialization of regular file: ',file

	  iunit = 0
	  iunit = ifileo(0,file,'form','old')
	  if( iunit .le. 0 ) goto 99

	  iunit = -iunit

	  nvar = 0
	  nx = 0
	  ny = 0
	  call rgf_intp(iunit,ndim,itact,itold,itnew,nvar,nx,ny
     +				,dfile(1,1),dfile(1,2),dfile(1,3))

	  if( nvar .ne. nvar0 ) goto 97

	  write(6,*) 'regular file opened :',iunit,nvar,nx,ny

	  ifile(1) = iunit
	  ifile(2) = itold
	  ifile(3) = itnew
	  ifile(4) = itact
	  ifile(5) = nvar
	  ifile(6) = nx
	  ifile(7) = ny
	end if

c-------------------------------------------------------------
c interpolate in time
c-------------------------------------------------------------

	itold = ifile(2)
	itnew = ifile(3)
	itact = ifile(4)
	nvar  = ifile(5)
	nx    = ifile(6)
	ny    = ifile(7)

	if( it .ne. itact ) then
	  itact = it			!bugfix ITACT
	  call rgf_intp(iunit,ndim,itact,itold,itnew,nvar,nx,ny
     +				,dfile(1,1),dfile(1,2),dfile(1,3))
	end if

	ifile(1) = iunit
	ifile(2) = itold
	ifile(3) = itnew
	ifile(4) = itact

c-------------------------------------------------------------
c end of routine
c-------------------------------------------------------------

	return
   97	continue
	write(6,*) nvar,nvar0
	stop 'error stop rgf_admin: number of variables'
   99	continue
	write(6,*) iunit,file
	stop 'error stop rgf_admin: cannot open file'
	end

c*********************************************************************

	subroutine rgf_intp(iunit,ndim,it,itold,itnew
     +                          ,nvar,nx,ny,dold,dnew,dact)

c handles time interpolation of regular fields
c
c file must alread be open on unit iunit
c the first time must be called with negative iunit
c sets iunit = 0 if EOF has been found

	implicit none

	integer iunit
	integer ndim
	integer it,itold,itnew
	integer nvar
	integer nx,ny
	real dold(nx,ny,nvar)
	real dnew(nx,ny,nvar)
	real dact(nx,ny,nvar)

	logical bdata

c--------------------------------------------------------------
c read first time
c--------------------------------------------------------------

	if( iunit .eq. 0 ) return	!EOF found before

	if( iunit .lt. 0 ) then
	  iunit = -iunit
	  call rgf_read(iunit,ndim,itnew,nvar,nx,ny,dnew,bdata)
	  if( .not. bdata ) goto 99

	  itold = itnew
	  call rgf_copy(nvar,nx,ny,dnew,dold)
	end if

	if( nx*ny*nvar .le. 0 ) goto 98

c--------------------------------------------------------------
c read new record
c--------------------------------------------------------------

	do while( it .gt. itnew )
	  itold = itnew
	  call rgf_copy(nvar,nx,ny,dnew,dold)
	  call rgf_read(iunit,ndim,itnew,nvar,nx,ny,dnew,bdata)
	  if( .not. bdata ) then
	    itnew = it
	    iunit = 0
	  end if
	end do

c--------------------------------------------------------------
c interpolate
c--------------------------------------------------------------

	call rgf_time_interpolate(it,itold,itnew
     +				,nvar,nx,ny,dold,dnew,dact)

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	return
   98	continue
	write(6,*) 'nvar,nx,ny: ',nvar,nx,ny
	stop 'error stop rgf_intp: parameters are zero'
   99	continue
	stop 'error stop rgf_intp: no data'
	end

c*********************************************************************

	subroutine rgf_time_interpolate(it,itold,itnew
     +				,nvar,nx,ny,dold,dnew,dact)

c time interpolation of regular field data

	implicit none

	integer it,itold,itnew
	integer nvar
	integer nx,ny
	real dold(nx,ny,nvar)
	real dnew(nx,ny,nvar)
	real dact(nx,ny,nvar)

	integer i,ix,iy
	real rit,do,dn

        if( itnew .gt. itold ) then
          rit=float(it-itold)/float(itnew-itold)
        else
          rit = 0.
        end if

	do i=1,nvar
	  do iy=1,ny
	    do ix=1,nx
		do = dold(ix,iy,i)
		dn = dnew(ix,iy,i)
		dact(ix,iy,i) = rit*(dn-do) + do
	    end do
	  end do
	end do

	end

c*********************************************************************
c*********************************************************************
c*********************************************************************

	subroutine rgf_copy(nvar,nx,ny,dsource,dtarget)

c copy regular field record

	implicit none

	integer nvar
	integer nx,ny
	real dsource(nx,ny,nvar)
	real dtarget(nx,ny,nvar)

	integer i,ix,iy

	do i=1,nvar
	  do iy=1,ny
	    do ix=1,nx
		dtarget(ix,iy,i) = dsource(ix,iy,i)
	    end do
	  end do
	end do

	end

c*********************************************************************
c*********************************************************************
c*********************************************************************

	subroutine rgf_read(iunit,ndim,it,nvar,nx,ny,data,bdata)

c read complete record of regular field

	implicit none

	integer iunit
	integer ndim
	integer it
	integer nvar
	integer nx,ny
	real data(ndim)
	logical bdata

	integer nvar0,nx0,ny0

c------------------------------------------------------------
c read header
c------------------------------------------------------------

	call rgf_read_header(iunit,it,nvar0,nx0,ny0,bdata)

c------------------------------------------------------------
c check for errors
c------------------------------------------------------------

	if( .not. bdata ) then			!no more data
	  return
	else if( nvar*nx*ny .le. 0 ) then	!first call
	  nvar = nvar0
	  nx = nx0
	  ny = ny0
	else if( nvar.ne.nvar0 .or. nx.ne.nx0 .or. ny.ne.ny0 ) then
	  write(6,*) 'nvar,nx,ny:    ',nvar,nx,ny
	  write(6,*) 'nvar0,nx0,ny0: ',nvar0,nx0,ny0
	  stop 'error stop rgf_read: incompatible parameters'
	end if

	if( nx*ny*nvar .gt. ndim ) then
	  write(6,*) 'nx,ny,nvar:  ',nx,ny,nvar
	  write(6,*) 'ndim,needed: ',ndim,nx*ny*nvar
	  stop 'error stop rgf_read: ndim'
	end if

c------------------------------------------------------------
c read data
c------------------------------------------------------------

	call rgf_read_data(iunit,nvar,nx,ny,data)

c------------------------------------------------------------
c end of routine
c------------------------------------------------------------

	end

c*********************************************************************

	subroutine rgf_read_header(iunit,it,nvar,nx,ny,bdata)

c reads header of regular field

	implicit none

	integer iunit
	integer it
	integer nvar
	integer nx,ny
	logical bdata
	real x0,y0,dx,dy,flag

	!read(iunit,*,end=1) it,nvar,nx,ny
	read(iunit,*,end=1) it,nvar,nx,ny,x0,y0,dx,dy,flag
	bdata = .true.

	return
    1	continue
	bdata = .false.
	end

c*********************************************************************

	subroutine rgf_read_data(iunit,nvar,nx,ny,data)

c reads data of regular field

	implicit none

	integer iunit
	integer nvar
	integer nx,ny
	real data(nx,ny,nvar)

	integer i,ix,iy
	character*80 description

	do i=1,nvar
	  read(iunit,'(a)') description
	  read(iunit,*) ((data(ix,iy,i),ix=1,nx),iy=1,ny)
	end do

	end

c*********************************************************************
c*********************************************************************
c*********************************************************************

	subroutine rgf_info(nout,ifile)

c prints info of regular field to unit nout

	implicit none

	integer nout
	integer ifile(7)

	integer unit,i

	unit = nout
	if( nout .le. 0 ) unit = 6

	write(unit,*) 'rgf_info: ',(ifile(i),i=1,7)

	end

c*********************************************************************

	subroutine rgf_print(nout,ifile,dfile)

c prints data of regular field to unit nout

	implicit none

	integer nout
	integer ifile(7)
	real dfile(*)

	integer unit,i
	integer nvar,nx,ny,n

	unit = nout
	if( nout .le. 0 ) unit = 6

        nvar  = ifile(5)
        nx    = ifile(6)
        ny    = ifile(7)
        n = nx * ny * nvar

	write(unit,*) 'rgf_info: ',(ifile(i),i=1,7)
	write(unit,*) (dfile(i),i=1,n)

	end

c*********************************************************************

	subroutine rgf_get_val(ifile,dfile,ivar,i,j,value)

c gets value of regular field

	implicit none

	integer ifile(7)
	real dfile(*)
	integer ivar		!variable to get
	integer i,j		!coordinates of variable in grid
	real value		!value of point (i,j) (return)

	integer nvar,nx,ny,n,ip

        nvar  = ifile(5)
        nx    = ifile(6)
        ny    = ifile(7)
        n = nx * ny

	if( i .le. 0 .or. i .gt. nx ) goto 99
	if( j .le. 0 .or. j .gt. ny ) goto 99
	if( ivar .le. 0 .or. ivar .gt. nvar ) goto 99

	ip = (ivar-1)*n + (j-1)*nx + i
	value = dfile(ip)

	return
   99	continue
	write(6,*) 'nvar,nx,ny: ',nvar,nx,ny
	write(6,*) 'ivar,i,j:   ',ivar,i,j
	stop 'error stop rgf_get_val: input parameters'
	end

c*********************************************************************

