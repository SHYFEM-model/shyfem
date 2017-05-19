
!=================================================================
	module nc_domain
!=================================================================

	implicit none

	integer, save :: ncx1,ncx2,ncy1,ncy2,ncz1,ncz2

!=================================================================
	end module nc_domain
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
     +			,bregular,regpar,iregpar)

        implicit none

        integer nx,ny
        real x(nx,ny)
        real y(nx,ny)
        logical bregular
        real regpar(7)
        real iregpar(9)

        integer ix,iy
        real xtot,ytot,eps,dx,dy,dxx,dyy
	real xmin,xmax,ymin,ymax
	double precision dxxtot,dyytot

        bregular = .true.
        regpar = 0.
        dx = 0.
        dy = 0.
	dxxtot = 0
	dyytot = 0

	xmin = minval(x)
	xmax = maxval(x)
        xtot = xmax - xmin
        eps = 1.e-5 * xtot
        dx = x(2,1) - x(1,1)
        do iy=1,ny
          do ix=2,nx
            dxx = x(ix,iy) - x(ix-1,iy)
            if( abs(dx-dxx) > eps ) bregular = .false.
	    dxxtot = dxxtot + dxx
          end do
        end do
	dxx = dxxtot / (ny*(nx-1))

	ymin = minval(y)
	ymax = maxval(y)
        ytot = ymax - ymin
        eps = 1.e-5 * ytot
        dy = y(1,2) - y(1,1)
        do iy=2,ny
          do ix=1,nx
            dyy = y(ix,iy) - y(ix,iy-1)
            if( abs(dy-dyy) > eps ) bregular = .false.
	    dyytot = dyytot + dyy
          end do
        end do
	dyy = dyytot / (nx*(ny-1))

        regpar(1) = nx
        regpar(2) = ny
        regpar(3) = x(1,1)
        regpar(4) = y(1,1)
        regpar(5) = dx
        regpar(6) = dy
        regpar(7) = -999.

        iregpar(1) = nx
        iregpar(2) = ny
        iregpar(3) = xmin
        iregpar(4) = ymin
        iregpar(5) = dxx
        iregpar(6) = dyy
        iregpar(7) = -999.
        iregpar(8) = xmax
        iregpar(9) = ymax

        end

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine handle_domain(dstring,bregular,regpar,iregpar)

	implicit none

	character*(*) dstring
	logical bregular
	real regpar(7)
	real iregpar(9)

	integer ianz,max
	integer nx,ny,i,nz
	integer ix1,ix2,iy1,iy2,iz1,iz2
	real xx0,yy0,xx1,yy1
	real x0,y0,x1,y1
	real dx,dy,x,y
	real f(7)

	integer iscanf
	logical nc_is_full_domain

	max = 7
	f = 0.
	ianz = iscanf(dstring,f,max)

	if( ianz == max ) then
	  write(6,*) 'too many values in domain string'
	  write(6,*) trim(dstring)
	  stop 'error stop handle_domain: too many values'
	end if

	if( bregular ) then
	  nx = nint(regpar(1))
	  ny = nint(regpar(2))
	  x0 = regpar(3)
	  y0 = regpar(4)
	  dx = regpar(5)
	  dy = regpar(6)
	  x1 = x0 + (nx-1)*dx
	  y1 = y0 + (ny-1)*dy

	  if( ianz == 0 ) then
	    f(1:2) = regpar(3:4)
	    f(3) = x0 + (nx-1)*dx
	    f(4) = y0 + (ny-1)*dy
	  else if( ianz /= 4 ) then
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

	  call nc_set_domain(ix1,ix2,iy1,iy2,iz1,iz2)

	end if
	
        write(6,*) 'regular domain: ',bregular
        write(6,*) 'nx,ny: ',nx,ny
        write(6,*) 'x0,y0: ',x0,y0
        write(6,*) 'x1,y1: ',x1,y1
        write(6,*) 'dx,dy: ',dx,dy

	nz = iz2
	if( .not. nc_is_full_domain(nx,ny,nz) ) then
	  write(6,*) 'changed regular domain: '
          write(6,*) 'ix1,iy1: ',ix1,iy1
          write(6,*) 'ix2,iy2: ',ix2,iy2
          write(6,*) 'nx,ny: ',ix2-ix1+1,iy2-iy1+1
          write(6,*) 'x0,y0: ',x0 + (ix1-1)*dx,y0 + (iy1-1)*dy
          write(6,*) 'x1,y1: ',x0 + (ix2-1)*dx,y0 + (iy2-1)*dy
          write(6,*) 'dx,dy: ',dx,dy
	end if

	end

!*****************************************************************

	subroutine recompute_regular_domain(regpar)

	implicit none

	real regpar(7)

	integer ix1,ix2,iy1,iy2,iz1,iz2
	integer nx,ny
	real x0,y0,dx,dy

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

	regpar(1) = nx
	regpar(2) = ny
	regpar(3) = x0
	regpar(4) = y0

	end

!*****************************************************************

