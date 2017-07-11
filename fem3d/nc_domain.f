
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
	real, save, allocatable :: valfem(:)

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

        implicit none

        integer nx,ny
        real x(nx,ny)
        real y(nx,ny)
        logical bregular		!return
        real regpar(9)			!return

        integer ix,iy
        real xtot,ytot,eps,dx,dy,dxx,dyy
	real xmin,xmax,ymin,ymax
	real flag
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

	if( bregular ) then
	  dxx = dx
	  dyy = dy
	end if

	flag = -999.
	call set_regpar(regpar,nx,ny,dxx,dyy,xmin,ymin,xmax,ymax,flag)

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

	subroutine prepare_no_interpol

	use nc_interpol

	implicit none

	do_interpolation = .false.

	end

!*****************************************************************

	subroutine prepare_interpol(nx,ny,xx,yy,regpar)

	use basin
	use nc_interpol

	implicit none

	integer nx,ny
	real xx(nx,ny)
	real yy(nx,ny)
	real regpar(9)

	integer ip,jp
	real dx,dy,x0,y0,x1,y1,flag

	call bas_insert_irregular(nx,ny,xx,yy)	!inserts data coordinates into basin

	call get_regpar(regpar,ip,jp,dx,dy,x0,y0,x1,y1,flag)
	call setgeo(x0,y0,dx,dy,flag)		!prepares for interpolation on reg

	nc_flag = flag
	do_interpolation = .true.

	allocate(fm(4,ip,jp))
	allocate(valfem(nkn))

	call av2fm(fm,ip,jp)			!computes interpolation matrix

	end

!*****************************************************************

	subroutine handle_interpol_2d(nx,ny,val,nxnew,nynew,valnew)

	use basin
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

	subroutine do_interpol_2d(nx,ny,val,nxnew,nynew,valnew)

	use basin
	use nc_interpol

	implicit none

	integer nx,ny
	real val(nx,ny)
	integer nxnew,nynew
	real valnew(nxnew,nynew)

	integer k,ix,iy,n,i
	integer, save :: ixx(4) = (/0,1,1,0/)
	integer, save :: iyy(4) = (/0,0,1,1/)
	real v,vv,flag

	flag = nc_flag

        k = 0
        do iy=1,ny
          do ix=1,nx
            k = k + 1
            valfem(k) = val(ix,iy)
          end do
        end do

        do iy=2,ny
          do ix=2,nx
            k = k + 1
	    v = 0
	    n = 0
	    do i=1,4
	      vv = val(ix-ixx(i),iy-iyy(i))
	      if( vv /= flag ) then
		n = n + 1
	        v = v + vv
	      end if
	    end do
            !v = val(ix,iy)+val(ix-1,iy)+val(ix,iy-1)+val(ix-1,iy-1)
	    if( n == 0 ) then
	      valfem(k) = flag
	    else
              valfem(k) = v/n
	    end if
          end do
        end do

	call fm2am2d(valfem,nxnew,nynew,fm,valnew)

	end

!*****************************************************************

	subroutine do_interpol_3d(nx,ny,nz,val,nxnew,nynew,valnew)

	implicit none

	integer nx,ny,nz
	real val(nx,ny,nz)
	integer nxnew,nynew
	real valnew(nxnew,nynew,nz)

	integer iz
	real val2d(nx,ny)
	real val2dnew(nxnew,nynew)

	do iz=1,nz
	  val2d(:,:) = val(:,:,iz)
	  call do_interpol_2d(nx,ny,val2d,nxnew,nynew,val2dnew)
	  valnew(:,:,iz) = val2dnew(:,:)
	end do

	end

!*****************************************************************

