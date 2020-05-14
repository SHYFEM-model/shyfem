
!--------------------------------------------------------------------------
!
!    Copyright (C) 2017-2019  Georg Umgiesser
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

! convert nc files to fem files: domain handling
!
! contents :
!
!
! revision log :
!
! 25.05.2017	ggu	changed VERS_7_5_28
! 13.06.2017	ggu	changed VERS_7_5_29
! 11.07.2017	ggu	changed VERS_7_5_30
! 03.07.2018	ggu	revision control introduced
! 13.07.2018	ggu	changed VERS_7_4_1
! 16.02.2019	ggu	changed VERS_7_5_60
! 22.10.2019	ggu	check in check_regular_coords() for negative vals
!
!*****************************************************************
!*****************************************************************
!*****************************************************************

!=================================================================
	module nc_domain
!=================================================================

	implicit none

	integer, save :: ncx1,ncx2,ncy1,ncy2,ncz1,ncz2

!=================================================================
	end module nc_domain
!=================================================================

!=================================================================
	module nc_interpol
!=================================================================

	implicit none

	logical, save :: do_interpolation = .false.
	real, save :: nc_flag = -999.
	real, save, allocatable :: fm(:,:,:)

	logical, save :: do_single = .false.
	integer, save :: nsingle = 0
	real, save, allocatable :: xsingle(:)
	real, save, allocatable :: ysingle(:)

!=================================================================
	end module nc_interpol
!=================================================================

!*****************************************************************

	subroutine nc_set_domain(ix1,ix2,iy1,iy2,iz1,iz2)

	use nc_domain

	implicit none

	integer ix1,ix2,iy1,iy2,iz1,iz2

	ncx1 = ix1
	ncx2 = ix2
	ncy1 = iy1
	ncy2 = iy2
	ncz1 = iz1
	ncz2 = iz2

	end

!*****************************************************************

	subroutine nc_get_domain(ix1,ix2,iy1,iy2,iz1,iz2)

	use nc_domain

	implicit none

	integer ix1,ix2,iy1,iy2,iz1,iz2

	ix1 = ncx1
	ix2 = ncx2
	iy1 = ncy1
	iy2 = ncy2
	iz1 = ncz1
	iz2 = ncz2

	end

!*****************************************************************

	function nc_is_full_domain(nx,ny,nz)

	use nc_domain

	implicit none

	logical nc_is_full_domain
	integer nx,ny,nz

	nc_is_full_domain = .false.

	if( ncx1 /= 1 .or. nx /= ncx2 ) return
	if( ncy1 /= 1 .or. ny /= ncy2 ) return
	if( ncz1 /= 1 .or. nz /= ncz2) return

	nc_is_full_domain = .true.

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

        subroutine check_regular_coords(nx,ny,x,y
     +			,bregular,regpar)

! sets regpar and bregular
!
! if bregular is .true. than regpar describes fully the grid
! if it is .false. dx,dy are average values, x0,y0,x1,y1 are min/max values
! also for regular grid, dx,dy are average values

        implicit none

        integer nx,ny
        real x(nx,ny)
        real y(nx,ny)
        logical bregular		!return
        real regpar(9)			!return

	logical bdebug
        integer ix,iy
        real xtot,ytot,eps,dx,dy,dxx,dyy
	real xmin,xmax,ymin,ymax
	real dxmin,dxmax,dymin,dymax
	real flag,epsilon
	double precision dxxtot,dyytot

	bdebug = .true.
	bdebug = .false.

        bregular = .true.
	epsilon = 1.e-4
        regpar = 0.
        dx = 0.
        dy = 0.
	dxxtot = 0
	dyytot = 0

	xmin = minval(x)
	xmax = maxval(x)
        xtot = xmax - xmin
        eps = epsilon * xtot
        dx = x(2,1) - x(1,1)
	dxmin = dx
	dxmax = dx
        do iy=1,ny
          do ix=2,nx
            dxx = x(ix,iy) - x(ix-1,iy)
	    dxxtot = dxxtot + dxx
	    dxmin = min(dxx,dxmin)
	    dxmax = max(dxx,dxmax)
          end do
        end do
	dxx = dxxtot / (ny*(nx-1))

	if( dxmin <= 0 ) then
	  write(6,*) 'coords spacing is negative or zero ',dxmin
	  write(6,*) 'use the follwoing command to invert coords:'
	  write(6,*) 'ncpdq -a "-lon" in.nc out.nc'
	  write(6,*) '(this assumes that your x-coord is named lon)'
	  stop 'error stop check_regular_coords: dx<=0'
	else if( dxmax-dxmin < eps ) then
	  bregular = .true.
	else
	  write(6,*) 'coords not regular...' 
	  write(6,*) dxmin,dxmax,dxmax-dxmin,eps
	  bregular = .false.
	end if

	!call invert_coords(1,nx,ny,y)	!test for Nador

	ymin = minval(y)
	ymax = maxval(y)
        ytot = ymax - ymin
        eps = epsilon * ytot
        dy = y(1,2) - y(1,1)
	dymin = dy
	dymax = dy
        do iy=2,ny
          do ix=1,nx
            dyy = y(ix,iy) - y(ix,iy-1)
	    dyytot = dyytot + dyy
	    dymin = min(dyy,dymin)
	    dymax = max(dyy,dymax)
          end do
        end do
	dyy = dyytot / (nx*(ny-1))

	if( dymin <= 0 ) then
	  write(6,*) 'coords spacing is negative or zero ',dymin
	  write(6,*) 'use the follwoing command to invert coords:'
	  write(6,*) 'ncpdq -a "-lat" in.nc out.nc'
	  write(6,*) '(this assumes that your y-coord is named lat)'
	  stop 'error stop check_regular_coords: dy<=0'
	else if( dymax-dymin < eps ) then
	  bregular = .true.
	else
	  write(6,*) 'coords not regular...' 
	  write(6,*) dymin,dymax,dymax-dymin,eps
	  bregular = .false.
	end if
	  
	if( bdebug ) then
	  write(6,*) 'info check_regular_coords -----------------'
	  write(6,*) 'bregular: ',bregular
	  write(6,*) nx,ny
	  write(6,*) xmin,xmax
	  write(6,*) ymin,ymax
	  write(6,*) dx,dxx,dxxtot
	  write(6,*) dy,dyy,dyytot
	  write(6,*) dxmin,dxmax
	  write(6,*) dymin,dymax
	  write(6,*) '-------------------------------------------'
	end if

	flag = -999.
	call set_regpar(regpar,nx,ny,dxx,dyy,xmin,ymin,xmax,ymax,flag)

        end

!*****************************************************************

	subroutine invert_coords(ii,nx,ny,coord)

	implicit none

	integer ii
	integer nx,ny
	real coord(nx,ny)

	integer ix,iy,il
	real aux

	if( ii == 0 ) then	!x-direction
	  do iy=1,ny
	    do ix=1,nx/2
	      il = nx+1-ix
	      aux=coord(ix,iy)
	      coord(ix,iy) = coord(il,iy)
	      coord(il,iy) = aux
	    end do
	  end do
	else			!y-direction
	  do ix=1,nx
	    do iy=1,ny/2
	      il = ny+1-iy
	      aux=coord(ix,iy)
	      coord(ix,iy) = coord(ix,il)
	      coord(ix,il) = aux
	    end do
	  end do
	end if
	
	end
	
!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine handle_domain(bverb,dstring
     +				,bregular,regpar_data,regpar)

! uses bregular and regpar_data as input
! computes regpar which is the regular output domain
! uses dstring to decide how to build the output grid
!
! format for dstring:
! if regpar_data is regular: x0,y0,x1,y0
! if regpar_data is irregular: dx,dy,x0,y0,x1,y0

	implicit none

	logical bverb
	character*(*) dstring	!dx,dy,x0,y0,x1,y0
	logical bregular
	real regpar_data(9)
	real regpar(9)

	integer n,max
	integer nx,ny,i,nz
	integer ix1,ix2,iy1,iy2,iz1,iz2
	real xx0,yy0,xx1,yy1
	real x0,y0,x1,y1
	real dx,dy,x,y
	real f(7)
	real flag

	integer iscanf
	logical nc_is_full_domain

	max = 7
	f = 0.
	n = iscanf(dstring,f,max)

	if( n == max ) then
	  write(6,*) 'too many values in domain string'
	  write(6,*) trim(dstring)
	  stop 'error stop handle_domain: too many values'
	end if

	if( bregular ) then
	  call handle_regular_domain(n,f,regpar_data,regpar)
	else
	  call handle_irregular_domain(n,f,regpar_data,regpar)
	end if
	
	call get_regpar(regpar_data,nx,ny,dx,dy,x0,y0,x1,y1,flag)
	if( bverb ) then
          write(6,*) 'original domain: ',bregular
          write(6,*) 'nx,ny: ',nx,ny
          write(6,*) 'x0,y0: ',x0,y0
          write(6,*) 'x1,y1: ',x1,y1
          write(6,*) 'dx,dy: ',dx,dy
	end if

	call get_regpar(regpar,nx,ny,dx,dy,x0,y0,x1,y1,flag)
	call nc_get_domain(ix1,ix2,iy1,iy2,iz1,iz2)
	if( bverb ) then
	  write(6,*) 'final regular domain: '
          write(6,*) 'ix1,iy1: ',ix1,iy1
          write(6,*) 'ix2,iy2: ',ix2,iy2
          write(6,*) 'nx,ny: ',nx,ny
          write(6,*) 'x0,y0: ',x0,y0
          write(6,*) 'x1,y1: ',x1,y1
          write(6,*) 'dx,dy: ',dx,dy
	end if

	end

!*****************************************************************

	subroutine recompute_regular_domain(regpar)

! use ix1,ix2,iy1,iy2,iz1,iz2 to recompute regular domain
! dx,dy is not changed

	implicit none

	real regpar(9)

	integer ix1,ix2,iy1,iy2,iz1,iz2
	integer nx,ny
	real x0,y0,dx,dy,x1,y1

	call nc_get_domain(ix1,ix2,iy1,iy2,iz1,iz2)

	nx = nint(regpar(1))
	ny = nint(regpar(2))
	x0 = regpar(3)
	y0 = regpar(4)
	dx = regpar(5)
	dy = regpar(6)

	nx = ix2-ix1+1
	ny = iy2-iy1+1
	x0 = x0 + (ix1-1)*dx
	y0 = y0 + (iy1-1)*dy
	x1 = x0 + (nx-1)*dx
	y1 = y0 + (ny-1)*dy

	regpar(1) = nx
	regpar(2) = ny
	regpar(3) = x0
	regpar(4) = y0
	regpar(8) = x1
	regpar(9) = y1

	end

!*****************************************************************

	subroutine handle_regular_domain(n,f,regpar_data,regpar)

! recompute domain from regular grid
!
! only x0,y0,x1,y1 are used
! dx,dy are not used and not changed

	implicit none

	integer n
	real f(7)
	real regpar_data(9)
	real regpar(9)

	integer nx,ny,i
	integer ix1,ix2,iy1,iy2,iz1,iz2
	real dx,dy,x0,y0,x1,y1
	real flag
	real xt,yt,xtt,ytt
	real xx0,yy0,xx1,yy1
	real x,y

	call get_regpar(regpar_data,nx,ny,dx,dy,x0,y0,x1,y1,flag)

	if( n == 0 ) then
	  f = 0.
	else if( n == 6 ) then
	  write(6,*) 'coordinates are regular'
	  write(6,*) 'to specify new domain we need 4 values:'
	  write(6,*) 'x0,y0,x1,y1'
	  write(6,*) '6 values are given, assuming dx,dy,x0,y0,x1,y1'
	  write(6,*) 'therefore dx,dy are ignored'
	  f(1:4) = f(3:6)
	else if( n /= 4 ) then
	  write(6,*) 'coordinates are regular'
	  write(6,*) 'to specify new domain we need 4 values:'
	  write(6,*) 'x0,y0,x1,y1'
	  stop 'error stop handle_domain: need 4 values'
	end if
	xx0 = f(1)
	yy0 = f(2)
	xx1 = f(3)
	yy1 = f(4)

	call nc_get_domain(ix1,ix2,iy1,iy2,iz1,iz2)

	if( n > 0 ) then
	  do i=1,nx
	    x = x0 + (i-1)*dx
	    if( x <= xx0 ) ix1 = i
	    if( x > xx1 ) exit
	  end do
	  ix2 = i

	  do i=1,ny
	    y = y0 + (i-1)*dy
	    if( y <= yy0 ) iy1 = i
	    if( y > yy1 ) exit
	  end do
	  iy2 = i
	end if

	call nc_set_domain(ix1,ix2,iy1,iy2,iz1,iz2)

	regpar = regpar_data
	call recompute_regular_domain(regpar)

	end

!*****************************************************************

	subroutine handle_irregular_domain(n,f,regpar_data,regpar)

! make regular grid from an irregular one
!
! can be used also for unstructured domains
! only dx,dy,x0,y0,x1,y1 are used, where dx,dy are requested resolution

	implicit none

	integer n
	real f(7)
	real regpar_data(9)
	real regpar(9)

	integer nx,ny
	real dx,dy,x0,y0,x1,y1
	real flag
	real xt,yt,xtt,ytt
	real x00,y00,x11,y11
	real, save :: eps = 1.e-5
	integer nxx,nyy
	integer ix1,ix2,iy1,iy2,iz1,iz2

	call get_regpar(regpar_data,nx,ny,dx,dy,x0,y0,x1,y1,flag)

	if( n == 0 ) then
	  !nothing to be changed
	else if( n == 1 ) then
	  dx = f(1)
	  dy = dx
	else if( n == 2 ) then
	  dx = f(1)
	  dy = f(2)
	else if( n == 6 ) then
	  dx = f(1)
	  dy = f(2)
	  x0 = f(3)
	  y0 = f(4)
	  x1 = f(5)
	  y1 = f(6)
	else
	  write(6,*) 'total number given: ',n
	  write(6,*) 'possible values are: 1,2,6'
	  write(6,*) 'format of domain: dx[,dy[,x0,y0,x1,y0]]'
	  stop 'error stop handle_irregular_domain: numbers given'
	end if
	
	xt = x1 - x0
	yt = y1 - y0
	nxx = 1 + ceiling(xt/dx)
	nyy = 1 + ceiling(yt/dy)
	xtt = (nxx-1) * dx
	ytt = (nyy-1) * dy

	if( xtt < xt .or. ytt < yt .or. nxx < 2 .or. nyy < 2 ) then
	  write(6,*) nx,ny,nxx,nyy
	  write(6,*) dx,dy
	  write(6,*) x0,y0,x1,y1
	  write(6,*) xt,yt,xtt,ytt
	  stop 'error stop handle_irregular_domain: internal error (1)'
	end if

	x00 = x0 - 0.5*(xtt-xt)
	y00 = y0 - 0.5*(ytt-yt)
	x11 = x00 + (nxx-1) * dx
	y11 = y00 + (nyy-1) * dy

	if( x0 < x00 .or. y00 < y00 .or. x1 > x11 .or. y1 > y11 ) then
	  write(6,*) nx,ny,nxx,nyy
	  write(6,*) dx,dy
	  write(6,*) x0,y0,x1,y1
	  write(6,*) x00,y00,x11,y11
	  write(6,*) xt,yt,xtt,ytt
	  stop 'error stop handle_irregular_domain: internal error (2)'
	end if

	call set_regpar(regpar,nxx,nyy,dx,dy,x00,y00,x11,y11,flag)

	call nc_get_domain(ix1,ix2,iy1,iy2,iz1,iz2)
	call nc_set_domain(ix1,nxx,iy1,nyy,iz1,iz2)

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine set_regpar(regpar,nx,ny,dx,dy,x0,y0,x1,y1,flag)

	implicit none

	real regpar(9)
	integer nx,ny
	real dx,dy,x0,y0,x1,y1
	real flag

        regpar(1) = nx
        regpar(2) = ny
        regpar(3) = x0
        regpar(4) = y0
        regpar(5) = dx
        regpar(6) = dy
        regpar(7) = flag
        regpar(8) = x1
        regpar(9) = y1

	end

!*****************************************************************

	subroutine get_regpar(regpar,nx,ny,dx,dy,x0,y0,x1,y1,flag)

	implicit none

	real regpar(9)
	integer nx,ny
	real dx,dy,x0,y0,x1,y1
	real flag

        nx = nint(regpar(1))
        ny = nint(regpar(2))
        x0 = regpar(3)
        y0 = regpar(4)
        dx = regpar(5)
        dy = regpar(6)
        flag = regpar(7)
        x1 = regpar(8)
        y1 = regpar(9)

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

	function must_interpol()

	use nc_interpol

	implicit none

	logical must_interpol

	must_interpol = do_interpolation

	end

!*****************************************************************

	function is_single()

	use nc_interpol

	implicit none

	logical is_single

	is_single = do_single

	end

!*****************************************************************

	subroutine prepare_no_interpol

	use nc_interpol

	implicit none

	do_interpolation = .false.

	end

!*****************************************************************

	subroutine prepare_interpol(nx,ny,xx,yy,regpar)

	use nc_interpol

	implicit none

	integer nx,ny
	real xx(nx,ny)
	real yy(nx,ny)
	real regpar(9)

	integer ip,jp
	real dx,dy,x0,y0,x1,y1,flag

	call bas_insert_irregular(nx,ny,xx,yy)	!inserts data coords into basin

	call get_regpar(regpar,ip,jp,dx,dy,x0,y0,x1,y1,flag)
	call setgeo(x0,y0,dx,dy,flag)		!prepares for interp on reg

	nc_flag = flag
	do_interpolation = .true.

	allocate(fm(4,ip,jp))

	call av2fm(fm,ip,jp)			!computes interpolation matrix

	end

!*****************************************************************

	subroutine handle_interpol_2d(nx,ny,val,nxnew,nynew,valnew)

	use nc_interpol

	implicit none

	integer nx,ny
	real val(nx,ny)
	integer nxnew,nynew
	real valnew(nxnew,nynew)

	integer nxx,nyy,nz

	if( do_interpolation ) then
	  call do_interpol_2d(nx,ny,val,nxnew,nynew,valnew)
	else
	  nxx = nx
	  nyy = ny
	  nz = 1
	  call compress_data(nxx,nyy,nz,val,valnew)
	  if( nxx /= nxnew .or. nyy /= nynew ) then
	    write(6,*) 'nx,nxnew: ',nx,nxnew
	    write(6,*) 'ny,nynew: ',ny,nynew
	    stop 'error stop handle_interpol_2d: wrong final dimensions'
	  end if
	end if

	end

!*****************************************************************

	subroutine insert_interpol_2d(nx,ny,val,ndim,valaux,flag)

! inserts values val in data structure

	implicit none

	integer nx,ny
	real val(nx,ny)
	integer ndim
	real valaux(ndim)
	real flag

	integer k,ix,iy,n,i
	integer, save :: ixx(4) = (/0,1,1,0/)
	integer, save :: iyy(4) = (/0,0,1,1/)
	real v,vv

        k = 0

        do iy=1,ny
          do ix=1,nx
            k = k + 1
            valaux(k) = val(ix,iy)
          end do
        end do

        do iy=2,ny
          do ix=2,nx
            k = k + 1
	    if( k > ndim ) goto 99
	    v = 0
	    n = 0
	    do i=1,4
	      vv = val(ix-ixx(i),iy-iyy(i))
	      if( vv /= flag ) then
		n = n + 1
	        v = v + vv
	      end if
	    end do
	    if( n == 0 ) then
	      valaux(k) = flag
	    else
              valaux(k) = v/n
	    end if
          end do
        end do

	return
   99	continue
	stop 'error stop insert_interpol_2d: internal errror ndim'
	end

!*****************************************************************

	subroutine do_interpol_2d(nx,ny,val,nxnew,nynew,valnew)

	use nc_interpol

	implicit none

	integer nx,ny
	real val(nx,ny)
	integer nxnew,nynew
	real valnew(nxnew,nynew)

	integer ndim
	real flag
	real, allocatable :: valaux(:)

	flag = nc_flag
	ndim = nx*ny + (nx-1)*(ny-1)
	allocate(valaux(ndim))

	call insert_interpol_2d(nx,ny,val,ndim,valaux,flag)
	call fm2am2d(valaux,nxnew,nynew,fm,valnew)

	end

!*****************************************************************

	subroutine do_interpol_3d(nx,ny,nz,val,nxnew,nynew,valnew)

	implicit none

	integer nx,ny,nz
	real val(nx,ny,nz)
	integer nxnew,nynew
	real valnew(nxnew,nynew,nz)

	integer iz
	!real val2d(nx,ny)
	!real val2dnew(nxnew,nynew)
	real, allocatable :: val2d(:,:)
	real, allocatable :: val2dnew(:,:)

	allocate(val2d(nx,ny))
	allocate(val2dnew(nxnew,nynew))

	do iz=1,nz
	  val2d(:,:) = val(:,:,iz)
	  call do_interpol_2d(nx,ny,val2d,nxnew,nynew,val2dnew)
	  valnew(:,:,iz) = val2dnew(:,:)
	end do

	end

!*****************************************************************

	subroutine prepare_single(sfile,ns,nx,ny,xx,yy,regpar)

	use nc_interpol

	implicit none

	character*(*) sfile
	integer ns
        integer nx,ny
        real xx(nx,ny)
        real yy(nx,ny)
        real regpar(9)

	integer ndim,ip,jp

	ns = 0
	if( sfile == ' ' ) return

	ndim = 0
	call read_single_points(sfile,ndim,ns,xsingle,ysingle)
	allocate(xsingle(ns),ysingle(ns))
	ndim = ns
	call read_single_points(sfile,ndim,ns,xsingle,ysingle)
	nsingle = ns

	write(6,*) 'single coordinates read: ',ns,trim(sfile)

	call bas_insert_irregular(nx,ny,xx,yy)	!inserts data coords into basin

	nc_flag = regpar(7)
	do_interpolation = .true.
	do_single = .true.
	regpar(1) = 0.
	regpar(2) = 0.

	ip = ns
	jp = 1
	allocate(fm(4,ip,jp))
	call av2fm_single(fm,ns,xsingle,ysingle)

	end

!*****************************************************************

	subroutine read_single_points(sfile,ndim,ns,xs,ys)

	implicit none

	character*(*) sfile
	integer ndim
	integer ns
	real xs(ndim)
	real ys(ndim)

	integer ios
	real x,y

	ns = 0
	if( sfile == ' ' ) return

	open(1,file=sfile,form='formatted',status='old')

	do
	  read(1,*,iostat=ios) x,y
	  if( ios < 0 ) exit
	  if( ios > 0 ) goto 99
	  ns = ns + 1
	  if( ndim == 0 ) cycle		!only check size
	  if( ns > ndim ) goto 98
	  xs(ns) = x
	  ys(ns) = y
	end do

	close(1)

	return
   98	continue
	write(6,*) 'ndim = ',ndim
	stop 'error stop read_single_points: ndim'
   99	continue
	write(6,*) 'line = ',ns
	stop 'error stop read_single_points: read error'
	end

!*****************************************************************

	subroutine get_single_points(ndim,ns,xs,ys)

	use nc_interpol

	implicit none

	integer ndim
	integer ns
	real xs(ndim)
	real ys(ndim)

	integer ios
	real x,y

	ns = nsingle
	if( ndim == 0 ) return

	if( ndim < ns ) goto 98

	xs(1:ns) = xsingle
	ys(1:ns) = ysingle

	return
   98	continue
	write(6,*) 'ndim,ns = ',ndim,ns
	stop 'error stop get_single_points: ndim'
	end

!*****************************************************************

