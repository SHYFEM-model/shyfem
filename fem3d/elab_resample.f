
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
! 05.11.2021	ggu	finished for use in femelab
!
!***************************************************************

	subroutine fem_resample_parse(string,regpar_in
     +				,regpar_out,nxn,nyn,idx0,idy0)

	use basin

	implicit none

	character*(*) string
	real regpar_in(7)
	real regpar_out(7)
	integer idx0,idy0,idx1,idy1

	logical bdebug
	integer ianz
	integer nxn,nyn,nx,ny
	real dx,dy,x0,y0,x1,y1,xaux,yaux
	real ddx,ddy
	real xn0,yn0,xn1,yn1
	real flag
	real f(10)
	double precision d(6)
	integer iscand,iscanf
	character*80 s

	bdebug = .false.
	regpar_out = 0.

	s = string
	if( index(trim(s),' ') > 0 ) then
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
	call resample_bounds(x0,dx,xn0,idx0,.true.)
	yn0 = f(2)
	call resample_bounds(y0,dy,yn0,idy0,.true.)
	xn1 = f(3)
	call resample_bounds(x1,dx,xn1,idx1,.false.)
	yn1 = f(4)
	call resample_bounds(y1,dy,yn1,idy1,.false.)

	nxn = nx + idx1 - idx0
	nyn = ny + idy1 - idy0

	write(6,*) 'resample values: ----------'
	write(6,*) nxn,nyn
	write(6,*) x0,y0,x1,y1
	write(6,*) xn0,yn0,xn1,yn1
	xaux = xn0 + (nxn-1)*dx
	yaux = yn0 + (nyn-1)*dy
	write(6,*) xn0,yn0,xaux,yaux
	write(6,*) idx0,idx1,idy0,idy1
	write(6,*) 'resample end --------------'

	regpar_out = (/float(nxn),float(nyn),xn0,yn0,dx,dy,flag/)

	end

!***************************************************************

	subroutine resample_bounds(xy,dxy,xyn,idxy,blow)

	implicit none

	real xy,dxy,xyn
	integer idxy
	logical blow

	integer n
	real f
	real, parameter :: flag = -999.
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

	!write(6,*) 'orig ',trim(s)
	do
	  i=index(s,',,')
	  if( i == 0 ) exit
	  s = s(1:i)//trim(sflag)//s(i+1:)
	  !write(6,*) 'loop ',trim(s)
	end do

	if( s(1:1) == ',' ) then
	  s = trim(sflag)//s
	  !write(6,*) 'init ',trim(s)
	end if

	i = len_trim(s)
	if( s(i:i) == ',' ) then
	  s = s(1:i)//trim(sflag)
	  !write(6,*) 'end  ',trim(s)
	end if

	end

!***************************************************************

	subroutine fem_resample_setup(nn,regpar_in)

	use basin

	implicit none

	integer nn
	real regpar_in(7)

	real, save :: regpar_save(7) = 0.

	if( any( regpar_save /= regpar_in ) ) then
	  if( any( regpar_save /= 0. ) ) then
	    stop 'error stop: size of regular grid cannot change'
	  end if
	  regpar_save = regpar_in
	end if

	end

!***************************************************************

	subroutine resample_data(flag,nlvdi,nx,ny,il_in,hd_in,data_in
     +				,nxn,nyn,idx,idy
     +				,il_out,hd_out,data_out)

	implicit none

	integer nlvdi,nvar
	integer nx,ny
	integer nxn,nyn
	integer idx,idy
	real flag
	integer il_in(nx,ny)
	real hd_in(nx,ny)
	real data_in(nlvdi,nx,ny)
	integer il_out(nxn,nyn)
	real hd_out(nxn,nyn)
	real data_out(nlvdi,nxn,nyn)

	integer ix,iy,ixto,iyto
	integer ixs,ixe,iys,iye
	integer l,lmax

	il_out = 1
	hd_out = flag
	data_out = flag

        call resample_index(nx,nxn,idx,ixs,ixe)
        call resample_index(ny,nyn,idy,iys,iye)

	if( ixs < 1 .or. ixe > nx ) goto 97
	if( iys < 1 .or. iye > ny ) goto 96

	!write(6,*) 'resample data: '
	!write(6,*) nx,nxn,idx,ixs,ixe
	!write(6,*) ny,nyn,idy,iys,iye
	!write(6,*) ixs,ixe,iys,iye
	!write(6,*) ixs-idx,ixe-idx,iys-idy,iye-idy

	do ix=ixs,ixe
	  ixto = ix - idx
	  if( ixto < 1 .or. ixto > nxn ) goto 99
	  do iy=iys,iye
	    iyto = iy - idy
	    if( iyto < 1 .or. iyto > nyn ) goto 98
	    lmax = il_in(ix,iy)
	    il_out(ixto,iyto) = il_in(ix,iy)
	    hd_out(ixto,iyto) = hd_in(ix,iy)
	    do l=1,lmax
	      data_out(l,ixto,iyto) = data_in(l,ix,iy)
	    end do
	  end do
	end do

	return
   97	continue
	write(6,*) ixs,ixe,nx
	stop 'error stop resample_data: ixs,ixe out of bound'
   96	continue
	write(6,*) iys,iye,ny
	stop 'error stop resample_data: iys,iye out of bound'
   99	continue
	write(6,*) ixto,nxn
	stop 'error stop resample_data: ixto out of bound'
   98	continue
	write(6,*) iyto,nyn
	stop 'error stop resample_data: iyto out of bound'
	end

!***************************************************************

        subroutine resample_index(n,nn,id0,is,ie)

        implicit none

        integer n,nn,id0,is,ie

        integer id1,imin,imax

        id1 = nn - n + id0

        imin = 1 + id0
        imax = n + id1

        is = max(1,imin)
        ie = min(n,imax)

        end

!***************************************************************

