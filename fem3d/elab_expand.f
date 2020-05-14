
!--------------------------------------------------------------------------
!
!    Copyright (C) 2017,2019  Georg Umgiesser
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
! 17.11.2017	ggu	changed VERS_7_5_37
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62

!**************************************************************************

!***************************************************************

	subroutine adjust_reg_vertical(nlvddi,nx,ny,flag,am,il)

	implicit none

	integer nlvddi
	integer nx,ny
	real flag
        real am(nlvddi,nx,ny)
	integer il(nx,ny)

	integer ix,iy,l

	do iy=1,ny
	  do ix=1,nx
	    do l=1,nlvddi
	      if( am(l,ix,iy) == flag ) exit
	    end do
	    il(ix,iy) = max(1,l-1)
	  end do
	end do

	end

!***************************************************************

	subroutine reg_expand_3d(nlvddi,nx,ny,lmax,regexpand,flag,am)

	implicit none

	integer nlvddi,nx,ny,lmax
	integer regexpand
	real flag
        real am(nlvddi,nx*ny)

	logical bdebug
	integer l
	real am2d(nx*ny)

	bdebug = .true.
	bdebug = .false.

	if( regexpand >= 0 ) then
	  if( lmax <= 1 ) then
	    am2d = am(1,:)
	    call reg_expand_2d(nx,ny,regexpand,flag,am2d)
	    am(1,:) = am2d
	  else
	    do l=1,lmax
	      if( l > nlvddi ) exit
	      am2d = am(l,:)
	      if( bdebug ) write(6,*) 'expand level ',l
	      call reg_expand_2d(nx,ny,regexpand,flag,am2d)
	      am(l,:) = am2d
	    end do
	  end if
	end if

	!call reg_debug(nlvdi,nx,ny,flag,am)
	
	end

!***************************************************************

	subroutine reg_expand_2d(nx,ny,regexpand,flag,am)

	implicit none

	integer nx,ny
	integer regexpand
	real flag
        real am(nx,ny)

	logical bdebug
	integer iexp
	integer ix,iy,i,n
	integer iflag
	integer ixflag(nx*ny)
	integer iyflag(nx*ny)
	real val
	real amaux(nx,ny)

	bdebug = .true.
	bdebug = .false.

	iexp = regexpand
	if( iexp < 0 ) return

	if( iexp == 0 ) iexp = max(nx,ny)

!--------------------------------------------------------
! find nodes with flag
!--------------------------------------------------------

	iflag = 0
	do iy=1,ny
	  do ix=1,nx
	    if( am(ix,iy) == flag ) then
	      iflag = iflag + 1
	      ixflag(iflag) = ix
	      iyflag(iflag) = iy
	    end if
	  end do
	end do

!--------------------------------------------------------
! loop over nodes with flag
!--------------------------------------------------------

	if( bdebug ) write(6,*) 'iflag: ',iexp,iflag

	do 
	  if( iexp <= 0 ) exit
	  amaux = am
	  do i=1,iflag
	    ix = ixflag(i)
	    iy = iyflag(i)
	    if( amaux(ix,iy) /= flag ) then
	      write(6,*) ix,iy,amaux(ix,iy)
	      stop 'error stop reg_expand: not a flag...'
	    end if
	    call box_intp(ix,iy,nx,ny,flag,am,n,val)
	    if( n > 0 ) then
	      amaux(ix,iy) = val
	      ixflag(i) = -1
	      iyflag(i) = -1
	    end if
	  end do
	  am = amaux
	  call compress_flag(iflag,ixflag,iyflag)
	  iexp = iexp - 1
	  if( bdebug ) write(6,*) 'iflag: ',iexp,iflag
	end do

	end

!***************************************************************

	subroutine compress_flag(iflag,ixflag,iyflag)

	implicit none

	integer iflag
	integer ixflag(iflag)
	integer iyflag(iflag)

	integer ifree,icopy,i

        logical is_free,is_valid
        is_free(i) = ixflag(i) < 0
        is_valid(i) = ixflag(i) > 0

        ifree = 0
        icopy = iflag+1

        do while( ifree .lt. icopy )

          do i=ifree+1,icopy-1
            if( is_free(i) ) exit
          end do
          ifree = i
          if( ifree .eq. icopy ) exit

          do i=icopy-1,ifree+1,-1
            if( is_valid(i) ) exit
          end do
          icopy = i
          if( ifree .eq. icopy ) exit

	  ixflag(ifree) = ixflag(icopy)
	  iyflag(ifree) = iyflag(icopy)
	  ixflag(icopy) = -1
	  iyflag(icopy) = -1

        end do

	iflag = ifree - 1

	end

!***************************************************************

	subroutine box_intp(ix,iy,nx,ny,flag,am,n,val)

	implicit none

	integer ix,iy,nx,ny
	real flag
	real am(nx,ny)
	integer n
	real val

	integer ix0,ix1,iy0,iy1
	integer i,j

	ix0 = max(1,ix-1)
	ix1 = min(nx,ix+1)
	iy0 = max(1,iy-1)
	iy1 = min(ny,iy+1)

	n = 0
	val = 0.
	do j=iy0,iy1
	  do i=ix0,ix1
	    if( am(i,j) /= flag ) then
	      n = n + 1
	      val = val + am(i,j)
	    end if
	  end do
	end do

	if( n > 0 ) val = val / n

	end

!***************************************************************

