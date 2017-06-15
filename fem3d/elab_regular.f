!
! routines for shy post processing output
!
! revision log :
!
! 23.09.2016    ggu     expand routine written
! 05.10.2016    ggu     routines dealing with regular field copied to own file
!
!***************************************************************

	subroutine fem_regular_parse(string,regpar,nx,ny)

	use basin

	implicit none

	character*(*) string
	real regpar(7)
	integer nx,ny

	logical bdebug
	integer ianz
	real dx,dy,x0,y0,x1,y1
	real ddx,ddy
	double precision d(6)
	integer iscand
	real, parameter :: flag = -999.

	bdebug = .false.
	regpar = 0.
	ianz = iscand(string,d,6)

	if( ianz == 0 ) then
	  nx = 0
	  ny = 0
	else if( ianz == 1 ) then
	  dx = d(1)
	  dy = dx
	else if( ianz == 2 ) then
	  dx = d(1)
	  dy = d(2)
	else if( ianz == 6 ) then
	  dx = d(1)
	  dy = d(2)
	else
	  write(6,*) ianz,'  ',trim(string)
	  write(6,*) 'read error or wrong number of parameters'
	  stop 'error stop fem_regular_parse: string'
	end if

	if( ianz == 0 ) then
	  return
	else if( ianz == 6 ) then
	  x0 = d(3)
	  y0 = d(4)
	  x1 = d(5)
	  y1 = d(6)
	else
	  x0 = minval(xgv)
	  y0 = minval(ygv)
	  x1 = maxval(xgv)
	  y1 = maxval(ygv)
	end if

        nx = 1 + ceiling((x1-x0)/dx)
        ny = 1 + ceiling((y1-y0)/dy)

	if( nx < 2 .or. ny < 2 ) then
          write(6,*) nx,ny
          write(6,*) dx,dy
          write(6,*) x0,y0,x1,y1
	  stop 'error stop fem_regular_parse: erroneous domain'
	end if

	if( ianz /= 6 ) then		!correct x0/y0
          x0 = dx * (int(x0/dx))
          y0 = dy * (int(y0/dy))
          x1 = dx * (int(x1/dx)+1)
          y1 = dy * (int(y1/dy)+1)
	end if

        x1 = x0 + (nx-1)*dx
        y1 = y0 + (ny-1)*dy

	if( bdebug ) then
	  write(6,*) 'fem_regular_parse: final domain'
          write(6,*) nx,ny
          write(6,*) dx,dy
          write(6,*) x0,y0,x1,y1
	end if

	regpar = (/float(nx),float(ny),x0,y0,dx,dy,flag/)

	call setgeo(x0,y0,dx,dy,flag)

	end

!***************************************************************

	subroutine fem_regular_setup(nx,ny,regpar,ilhv
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

	call makehkv_minmax(hkv,+1)
	call av2fm(fmreg,nx,ny)
	call fm_extra_setup(nx,ny,fmextra)

	ilcoord = 0
	hcoord = 0.
	
	x0 = regpar(3)
	y0 = regpar(4)
	dx = regpar(5)
	dy = regpar(6)

	do iy=1,ny
	  y = y0 + (iy-1)*dy
	  do ix=1,nx
	    x = x0 + (ix-1)*dx
	    xcoord(ix,iy) = x
	    ycoord(ix,iy) = y
	    ie = fmreg(4,ix,iy)
	    if( ie > 0 ) then
	      ilcoord(ix,iy) = ilhv(ie)
	      hcoord(ix,iy) = sum(hm3v(:,ie))/3.	!average in element
	    end if
	  end do
	end do

	! the next way to compute depth is compatible with nos/ous2nc

	call fm2am2d(hkv,nx,ny,fmreg,hcoord)	!other way to compute depth

	do ix=1,nx
	  x = x0 + (ix-1)*dx
	  xlon(ix) = x
	end do

	do iy=1,ny
	  y = y0 + (iy-1)*dy
	  ylat(iy) = y
	end do

	end

!***************************************************************

	subroutine fem_regular_interpolate(nx,ny,regexpand,lmax
     +			,fmreg,fmextra,il,cv3,am)

        use basin
	use levels
	!use shyelab_out

        implicit none

	integer nx,ny
	integer regexpand
        integer lmax                   !vertical dimension of regular array
	real fmreg(4,nx,ny)
	real fmextra(6,nkn)
	integer il(nx,ny)
        real cv3(nlvdi,nkn)            !values of fem array
        real am(nlvdi,nx*ny)     !interpolated values (return)

	logical binelem,bfromnode
	integer mode
	integer l
	real flag
	real fem2d(nkn)
	real am2d(nx*ny)

	mode = 3	!1: only in element, 3: only from nodes, 2: both
	mode = 2

	binelem = mode <= 2
	bfromnode = mode >= 2

	call getgeoflag(flag)

	am = flag

	if( binelem ) then
	  if( lmax <= 1 ) then
	    am2d = am(1,:)
	    fem2d = cv3(1,:)
	    call fm2am2d(fem2d,nx,ny,fmreg,am2d)
	    am(1,:) = am2d
	  else
	    call fm2am3d(nlvdi,ilhv,cv3,nlvdi,nx,ny,fmreg,am)
	  end if
	end if

	if( bfromnode ) then
	  if( lmax <= 1 ) then
	    am2d = am(1,:)
	    fem2d = cv3(1,:)
	    call fm_extra_2d(nx,ny,fmextra,fem2d,am2d)
	    am(1,:) = am2d
	  else
	    call fm_extra_3d(nlvdi,nlv,ilhkv,nx,ny,fmextra,cv3,am)
	  end if
	end if

	!call reg_debug(nlvdi,nx,ny,flag,am)

	if( regexpand >= 0 ) then
	  call reg_expand_3d(nlvdi,nx,ny,lmax,regexpand,flag,am)
	end if

	call adjust_reg_vertical(nlvdi,nx,ny,flag,am,il)
	
	end

!***************************************************************

	subroutine reg_debug(nlvddi,nx,ny,flag,am)

	implicit none

	integer nlvddi,nx,ny
	real flag
	real am(nlvddi,nx,ny)

	integer ix,iy,l,ix0,iy0
	real r

	write(6,*) 'reg_debug: ',nx,ny

	do iy=1,ny
	  if( am(1,23,iy) /= flag .and. am(1,22,iy) == flag ) then
	    write(6,*) 'ggguuu: ',iy,am(1,23,iy)
	  end if
	end do

	ix0 = 25
	iy0 = ny+1-51-5
	write(6,*) 'special: ',ix0,iy0,am(1:4,ix0,iy0)

	do ix=ix0-2,ix0+2
	  do iy=iy0-2,iy0+2
	    write(6,*) ix,iy,am(1:3,ix,iy)
	  end do
	end do

	ix=25
	iy=ny+1-51
	r=-5.
	write(6,*) 'setting: ',ix,iy,r
	am(1,ix,iy) = r
	am(1,ix,iy+1) = r
	am(1,ix+1,iy) = r
	am(1,ix+1,iy+1) = r

	end

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
!***************************************************************
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

