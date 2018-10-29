
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

	program tvis

c visualize moving vectors

	implicit none

	integer nzdim,nydim
c	parameter (nzdim=8,nydim=9)
	parameter (nzdim=8,nydim=25)

	character*60 file
	character*80 line
	real v(nydim,nzdim)
	real w(nydim,nzdim)
	real v1(nydim,nzdim)
	real w1(nydim,nzdim)

	logical bwait
	integer ixmin,iymin,ixmax,iymax
	integer ns,i,j,ierr,it,itnew
	integer nrec
	real ymin,ymax,zmin,zmax
	real xvmin,yvmin,xvmax,yvmax
	real yax1,zax1,yax2,zax2
	real fact

	nrec=0

c-------------------------------------------------------------

	ymin=0.		!min/max for plot
	ymax=10.
	ymax=26.
	zmin=0.
	zmax=9.

	yax1=0.5	!x/y for axis
	zax1=0.5
	yax2=9.5	!x/y for axis
	yax2=25.5	!x/y for axis
	zax2=8.5

	fact=0.01
	fact=0.1

c-------------------------------------------------------------

c getxwd returns size of window, that will be opened
c setxwd sets size of window to open (only befor qopen is called)
c (these units are in pixel, all subsequent units are real values)

	call getxwd(ixmin,iymin,ixmax,iymax)
	write(6,*) ixmin,iymin,ixmax,iymax
c	call setxwd(ixmin,iymin,ixmax+100,ixmax+100)

	write(6,*) 'Enter slowing parameter [0...9] : '
	read(5,'(i10)') ns

	if(ns.lt.0) then
	  ns=-ns
	  bwait=.true.
	else
	  ns=ns*100000
	  bwait=.false.
	end if
	write(6,*) ns

	write(6,*) 'Enter data file : '
	read(5,'(a)') file
	write(6,*) file
	if(file.eq.' ') stop
	open(66,err=98,file=file,form='formatted',status='old')

c qopen opens the window - has to be called once at the beginning

	call qopen

c qgetvp returns the dimension of the viewport (terminal) in viewport coord.
c (units of viewport coordinates are approx. in cm)

	call qgetvp(xvmin,yvmin,xvmax,yvmax)
	write(6,*) xvmin,yvmin,xvmax,yvmax

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	do j=1,nzdim
	  do i=1,nydim
	    v1(i,j)=0.
	    w1(i,j)=0.
	  end do
	end do

	call rdvect(66,v,w,nydim,nzdim,itnew,ierr)

	call qworld(ymin,zmin,ymax,zmax)

	call qstart

    2	continue

c qworld sets new world (window) coordinates - note that the corners
c ...of the viewport and window are coincident - if you want to avoid
c ...distortion, choose the dimensions of the world window and the
c ...viewport so that that the x/y is the same in both systems

c change color with qnewp [1-16]

c delete old line

	call qnewp(16)
	call plvect(v1,w1,nydim,nzdim,fact)

c plot axis

	call qnewp(1)
	call qline(yax1,zax1,yax2,zax1)
	call qline(yax2,zax1,yax2,zax2)
	call qline(yax2,zax2,yax1,zax2)
	call qline(yax1,zax2,yax1,zax1)

c plot new line

	call qnewp(2)
	call plvect(v,w,nydim,nzdim,fact)
	  
	call qout

	nrec=nrec+1
	it=itnew

	do j=1,nzdim
	  do i=1,nydim
	    v1(i,j)=v(i,j)
	    w1(i,j)=w(i,j)
	  end do
	end do

	call rdvect(66,v,w,nydim,nzdim,itnew,ierr)
	if(ierr.ne.0) goto 99

	if(bwait) then
		write(6,*) it,nrec
		if( mod(nrec,ns).eq.0 ) call wait
	else
		call slow(ns)
	end if

	goto 2

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   99	continue

	write(6,*) 'exit at time = ',it

	call qend

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c call qclose at the end of the session

	write(6,*) 'Press <CR> to exit ...'
	call wait

	call qclose

	stop
   98	write(6,*) 'No such file : '
	write(6,*) file
	stop 'error stop'
	end

c**************************************************************

	subroutine wait

	read(5,'(i10)') i

	end

	subroutine slow(ns)

	a=2.
	do i=1,ns
	  a=a*a
	  if(a.gt.10000.) a=2.
	end do

	return
	end

c**************************************************************

	subroutine rdvect(iunit,v,w,nydim,nzdim,it,ierr)

c reads in data

	implicit none

	integer iunit,ierr,it
	integer nydim,nzdim
	real v(nydim,nzdim)
	real w(nydim,nzdim)

	integer i,j,idum
	character*80 line

    1	continue

	read(iunit,'(a)',end=2) line

	if(line(1:4).eq.'vprv') then
	  do j=nzdim,1,-1
	    read(iunit,'(i2,2x,25f3.0)') idum,(v(i,j),i=1,nydim)
	  end do
	else if(line(1:4).eq.'wprv') then
	  read(line,'(10x,i7)') it
	  do j=nzdim,1,-1
	    read(iunit,'(i2,2x,25f3.0)') idum,(w(i,j),i=1,nydim)
	  end do
	  ierr=0
	  return
	end if

	goto 1
    2	continue

	ierr=-1

	return
	end

c*****************************************************

	subroutine plvect(v,w,nydim,nzdim,fact)

c plots vector

	implicit none

	integer nydim,nzdim
	real fact
	real v(nydim,nzdim)
	real w(nydim,nzdim)

	integer i,j
	real y,z,vv,ww
	real dd

	dd=0.02

	do j=1,nzdim
	  do i=1,nydim
	    y=i
	    z=j
	    vv=v(i,j)*fact
	    ww=w(i,j)*fact
	    call qline(y,z,y+vv,z+ww)
	    call qline(y-dd,z-dd,y+dd,z-dd)
	    call qline(y+dd,z-dd,y+dd,z+dd)
	    call qline(y+dd,z+dd,y-dd,z+dd)
	    call qline(y-dd,z+dd,y-dd,z-dd)
	  end do
	end do

	return
	end
