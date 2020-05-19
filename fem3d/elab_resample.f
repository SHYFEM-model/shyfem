
!--------------------------------------------------------------------------
!
!    Copyright (C) 2020  Georg Umgiesser
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

! routines to resample regular grid
!
! revision log :
!
! 21.04.2020	ggu	written from scratch
! 18.05.2020	ggu	some more improvements, but still not ready
!
!***************************************************************

	subroutine fem_resample_parse(string,regpar_in
     +				,regpar_out,nxn,nyn,id0x,id0y)

	use basin

	implicit none

	character*(*) string
	real regpar_in(7)
	integer id0x,id0y,id1x,id1y

	logical bdebug
	integer ianz
	real dx,dy,x0,y0,x1,y1
	real ddx,ddy
	double precision d(6)
	integer iscand
	real, parameter :: flag = -999.
	real regpar_out(7)
	character*80 s

	bdebug = .false.
	regpar = 0.

	s = string
	if( index(s,' ') > 0 ) then
	  write(6,*) 'no white space allowed in string: ',trim(s)
	  stop 'error stop fem_resample_parse: white space'
	end if
	call replace_with_flag(s,'-999.')

	nx = nint(regpar_in(1))
	ny = nint(regpar_in(2))
	x0 = regpar_in(3)
	y0 = regpar_in(4)
	dx = regpar_in(5)
	dy = regpar_in(6)
	flag = regpar_in(7)

	x1 = x0 + (nx-1)*dx
	y1 = y0 + (ny-1)*dy

	ianz = iscanf(s,f,4)

	if( ianz /= 4 ) then
	  write(6,*) 'resample string must have 4 values... cannot parse'
	  write(6,*) trim(string)
	  write(6,*) trim(s)
	  stop 'error stop fem_resample_parse: cannot parse'
	end if

	xn0 = f(1)
	call resample_bounds(x0,dx,xn0,id0x,.true.)
	yn0 = f(2)
	call resample_bounds(y0,dy,yn0,id0y,.true.)
	xn1 = f(3)
	call resample_bounds(x1,dx,xn1,id1x,.false.)
	yn1 = f(4)
	call resample_bounds(y1,dy,yn1,id1y,.false.)

	nxn = idx1 - idx0
	nyn = idy1 - idy0

	regpar_out = (/float(nxn),float(nyn),xn0,yn0,dx,dy,flag/)

	end

!***************************************************************

	subroutine fem_resample_setup1(nx,ny,regpar,ilhv
     +				,fmreg,fmextra
     +				,ilcoord,xcoord,ycoord,hcoord
     +				,xlon,ylat)

	use basin

	implicit none

	integer nx,ny
	real regpar(7)
	integer ilhv(nel)
	real fmreg(4,nx,ny)
	real fmextra(6,nkn)
	integer ilcoord(nx,ny)
	real xcoord(nx,ny)
	real ycoord(nx,ny)
	real hcoord(nx,ny)
	real xlon(nx)
	real ylat(ny)

	integer ix,iy,ie
	real dx,dy,x0,y0,x,y
	real hkv(nkn)
	real, save :: regpar_save(7) = 0.
	real, save :: regpar_out(7) = 0.

	if( any( regpar_save /= regpar_in ) ) then	!must initialize
	  call fem_resample_parse(string,regpar_in
     +				,regpar_out,nxn,nyn,id0x,id0y)
	  regpar_save = ragpar_in
	  if( allocated(hd_out) ) deallocate(hd_out)
	  if( allocated(il_out) ) deallocate(il_out)
	  if( allocated(data_out) ) deallocate(data_out)
	  allocate(hd_out(nxn,nyn))
	  allocate(il_out(nxn,nyn))
	  allocate(data_out(nlvdi,nxn,nyn))
	end if

	np = nxn*nyn

	end

!***************************************************************

	subroutine resample_bounds(xy,dxy,xyn,idxy,blow)

	implicit none

	real xy,dxy,xyn
	integer idxy
	logical blow

	real, parameter :: eps = 1.e-3

	if( xyn == flag .or. abs(xyn-xy)/dxy < eps ) then
	  xyn = xy
	  idxy = 0
	  return
	end if

	f = (xyn - xy) / dxy

	if( blow ) then
	  n = floor(f+eps)
	else
	  n = ceiling(f-eps)
	end if

	xyn = xy + n*dxy
	idxy = n

	end

!***************************************************************

	subroutine replace_with_flag(s,sflag)

	implicit none

	character*(*) s,sflag

	integer i

	do
	  i=index(s,',,')
	  if( i == 0 ) exit
	  s = s(1:i)//trim(sflag)//s(i+1:)
	end do

	if( s(1:1) == ',' ) then
	  s = trim(sflag)//s
	end if

	i = len_trim(s)
	if( s(i:i) == ',' ) then
	  s = s//trim(sflag)
	end if

	end

!***************************************************************

	subroutine compress_spaces(string)

	implicit none

	character*(*) string

	integer l,is,i
	character*1 c

	string = adjustl(string)
	if( index(string,'  ') <= 0 ) return	!no spaces to compress

	l = len_trim(string)

	is = 0
	do i=1,l
	  is = is + 1
	  c = string(is:is)
	  string(i:i) = c
	  if( c == ' ' ) then
	    do
	      c = string(is+1:is+1)
	      if( c /= ' ' ) exit
	      is = is + 1
	    end do
	  end if
	end do

	end

!***************************************************************

	subroutine resample_2d(nx,ny,data_in,nxn,nyn,idx,idy,data_out)

	implicit none

	integer nx,ny
	integer nxn,nyn
	integer idx,idy
	real data_in(nx*ny)
	real data_out(nxn*nyn)

	do ix=1,nxn
	  do iy=1,nyn
	    ito 
	    data_out(ito) = data_in(ifrom)
	  end do
	end do

	end

!***************************************************************

