
!--------------------------------------------------------------------------
!
!    Copyright (C) 2018-2019  Georg Umgiesser
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
! 18.12.2018	ggu	changed VERS_7_5_52
! 14.02.2019	ggu	changed VERS_7_5_56
! 21.05.2019	ggu	changed VERS_7_5_62

	program bug

	implicit none

	integer ndim
	parameter (ndim=2000)

	integer icaver,icmax
	integer n,nn
	real xt(3,ndim)
	real yt(3,ndim)
	real xl(2,ndim)
	real yl(2,ndim)

	real xmin,ymin,xmax,ymax

	call qopen

	do
	  nn = 1000*(rand(0)-0.5)
	  read(5,*,end=9,err=9) icaver,icmax
	  write(6,*) icaver,icmax,nn
	  if( icaver .eq. 0 ) then
	    n = ndim
	    call read_iter(n,xt,yt,xl,yl)
	    if( n > ndim) stop 'error stop: ndim'
	    call get_minmax(n,xt,yt,xl,yl,xmin,ymin,xmax,ymax)
	    write(6,*) n,xmin,ymin,xmax,ymax
	    write(6,*) 'plotting...'
	    call plot_fetch(n,xt,yt,xl,yl,xmin,ymin,xmax,ymax)
	  end if
	end do

    9	continue

	call qclose

	end

c********************************************************

	subroutine plot_fetch(n,xt,yt,xl,yl,xmin,ymin,xmax,ymax)

	implicit none

	integer n
	real xt(3,n)
	real yt(3,n)
	real xl(2,n)
	real yl(2,n)
	real xmin,ymin,xmax,ymax
	real xmm,ymm

	integer i,ii
	real fact

	fact = 0.3

	call qstart

	call qsetvp(1.,1.,15.,15.)
	xmm = xmin + fact*(xmax-xmin)
	ymm = ymin + fact*(ymax-ymin)
	call qworld(xmin,ymin,xmm,ymm)

	do i=1,n
	  !write(6,*) i,xl(1,i),yl(1,i),xl(2,i),yl(2,i)
	  call qline(xl(1,i),yl(1,i),xl(2,i),yl(2,i))
	  call qmove(xt(3,i),yt(3,i))
	  do ii=1,3
	    call qplot(xt(ii,i),yt(ii,i))
	  end do
	end do

	call qend

	end

c********************************************************

	subroutine get_minmax(n,xt,yt,xl,yl,xmin,ymin,xmax,ymax)

	implicit none

	integer n
	real xt(3,n)
	real yt(3,n)
	real xl(2,n)
	real yl(2,n)
	real xmin,ymin,xmax,ymax

	integer i,ii

	xmin = xt(1,1)
	xmax = xt(1,1)
	ymin = yt(1,1)
	ymax = yt(1,1)

	do i=1,n
	  do ii=1,3
	    xmin=min(xmin,xt(ii,i))
	    ymin=min(ymin,yt(ii,i))
	    xmax=max(xmax,xt(ii,i))
	    ymax=max(ymax,yt(ii,i))
	  end do
	  do ii=1,2
	    xmin=min(xmin,xl(ii,i))
	    ymin=min(ymin,yl(ii,i))
	    xmax=max(xmax,xl(ii,i))
	    ymax=max(ymax,yl(ii,i))
	  end do
	end do

	end

c********************************************************

	subroutine read_iter(n,xt,yt,xl,yl)

	implicit none

	integer n
	real xt(3,n)
	real yt(3,n)
	real xl(2,n)
	real yl(2,n)

	integer ic,ie,iei,ien,ndim
	real wdir

	ndim = n

	read(5,*) ie,wdir
	read(5,*)
	read(5,*)
	write(6,*) 'starting from element ',ie,wdir

	n = 0
	do while( n < 1000 )
	  n = n + 1
	  call read_intersect(n,iei,xt,yt)
	  call read_fetch(n,ic,ie,ien,xl,yl)
	  if( ic .ne. n ) then
	    write(6,*) ic,n
	    stop 'error stop: icount ne n'
	  end if
	  !write(6,*) n,ic,iei,ie,ien
	end do

	read(5,*)
	read(5,*) 

	end

c********************************************************

	subroutine read_intersect(n,iei,xt,yt)

	implicit none

	integer n,iei
	real xt(3,n)
	real yt(3,n)

	integer ii
	integer iaux(5)
	character*30 text

	read(5,*)
	read(5,*) iei
	read(5,*) (xt(ii,n),ii=1,3)
	read(5,*) (yt(ii,n),ii=1,3)
	read(5,*)

	do
	  read(5,*) iaux
	  if( iaux(1) == 0 ) exit
	  if( iaux(1) == 9 ) then
	    read(5,'(a30)') text
	    !write(6,*) text
	  end if
	end do

	read(5,*)
	read(5,*)

	end

c********************************************************

	subroutine read_fetch(n,ic,ie,ien,xl,yl)

	implicit none

	integer n,ic,ie,ien
	real xl(2,n)
	real yl(2,n)

	read(5,*)
	read(5,*) ic
	read(5,*) ie,ien
	read(5,*) xl(1,n),yl(1,n),xl(2,n),yl(2,n)
	read(5,*)
	read(5,*)

	end

c********************************************************



