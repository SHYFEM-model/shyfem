c
c $Id: supiso.f,v 1.12 2010-02-26 15:29:19 georg Exp $
c
c routines for plotting isoline and colering area
c
c contents :
c
c subroutine isoline(val,nval,dis,mode)
c subroutine pliso(f,x,y,dist,isoanz,fnull,fiso)
c subroutine plcol(x,y,z,color,rlev,ncol,fnull)
c 
c revision log :
c
c 17.10.2001    ggu     use parameter isolin in isoline()
c 05.12.2003    ggu     filling is done in fill_area()
c 05.12.2003    ggu     is_r_nan() introduced for BUG search
c 11.03.2004    ggu     debugging in fill_area()
c 02.03.2005    ggu     in isoline: only plot if not flag
c 06.12.2008    ggu     plot isolines only where area is existing (bplot)
c 09.01.2009    ggu     deleted plcol0 (not used), new set_fxy_vals()
c 09.01.2009    ggu     plot only part of isolines using make_single_isolines()
c 23.02.2010    ggu     restructured and commented, use generic color table
c 18.08.2011    ggu     use isoinp to decide if interpolate or not in element
c 29.10.2013    ggu     new routine pldif (plotting isolines between areas)
c 18.04.2018    ggu     use bplot to decide what element to plot
c
c****************************************************

	subroutine isoline(val,nval,dis,mode)

c plots color and isolines in elements
c
c if isoanz in /isolin/ is 0 dis is used to determine isolines
c ...else isoanz gives number of isolines on fiso
c
c val           array with values
c nval          dimension of val
c dis           distance of isolines (dis>0)
c mode		0: isolines with dis  1: isolines in /isolin/
c		2: color  3: color on element values

	use basin
	use color
	use mod_hydro_plot

	implicit none

c argument
	integer nval
	real val(nval)
	real dis
	integer mode
c local
	character*80 line
	real f(3),x(3),y(3)
	real dist,flag
	integer ie,ii,kn
	integer inull
	integer isolin
	integer isoinp
	integer icsave

        real getpar
c       logical is_r_nan

c--------------------------------------------------------------------
c initialization
c--------------------------------------------------------------------

        isolin = nint(getpar('isolin'))	!plot isoline also for color
        isoinp = nint(getpar('isoinp'))	!interpolate in element?

	if( mode .le. 2 .and. nkn .ne. nval ) then
		write(6,*) 'nval must be nkn :',nval,nkn
		write(6,*) 'Cannot execute routine isoline'
		stop 'error stop isoline: nval /= nkn'
	end if

	call get_flag(flag)

c--------------------------------------------------------------------
c find isolines to be plotted
c--------------------------------------------------------------------

	if( mode .eq. 0 ) then
		dist = dis
	else
		dist = 0.
	end if

c--------------------------------------------------------------------
c loop over elements
c--------------------------------------------------------------------

	!write(6,*) 'isoline: isoanz... ',isoanz,mode

	if( mode .ge. 2 ) then

c	  -----------------------------------------
c	  color plot
c	  -----------------------------------------

	  call qcomm('plotting elements')
	  call get_color_table(icsave)
	  !call set_color_table(-1)
	  do ie=1,nel
	    if( .not. bplot(ie) ) cycle
	    call set_fxy_vals(ie,flag,val,f,x,y,inull)
	    if( mode .eq. 3 ) then				!element values
	      f = val(ie)
	      call plcol(x,y,f,ciso,fiso,isoanz+1,fnull)
	    else if( isoinp .gt. 0 ) then			!interpolate
	      call plcol(x,y,f,ciso,fiso,isoanz+1,fnull)
	    else
	      call plnode(x,y,f,ciso,fiso,isoanz+1,fnull)	!plot node
	    end if
	  end do
	  call set_color_table(icsave)
	  call qcomm('finished plotting elements')

	else

c	  -----------------------------------------
c	  isoline plot
c	  -----------------------------------------

	  do ie=1,nel
	    if( .not. bplot(ie) ) cycle
	    call set_fxy_vals(ie,flag,val,f,x,y,inull)
	    if( inull .eq. 0 ) then
	      call pliso(f,x,y,dist,isoanz,flag,fiso)
	    end if
	  end do

	end if

c--------------------------------------------------------------------
c plot isolines
c--------------------------------------------------------------------

        if( mode .eq. 2 .and. isolin .gt. 0 ) then	!plot single isolines
	  call qcomm('plotting isolines')
          call qgray(0.)
	  call make_single_isolines(isolin)		!set nriso,riso
          do ie=1,nel
	    if( .not. bplot(ie) ) cycle
	    call set_fxy_vals(ie,flag,val,f,x,y,inull)
	    if( inull .eq. 0 ) then
	      if( isoinp .gt. 0 ) then
		call pliso(f,x,y,0.,nriso,flag,riso)
	      else
		call pldif(f,x,y,0.5,flag)
	      end if
	    end if
          end do
	  call qcomm('finished plotting isolines')
        end if
	
c--------------------------------------------------------------------
c	end of routine
c--------------------------------------------------------------------

	end

c****************************************************

	subroutine isoreg(regval,nval,regpar,dis,mode)

c plots color and isolines of regular grid
c
c if isoanz in /isolin/ is 0 dis is used to determine isolines
c ...else isoanz gives number of isolines on fiso
c
c val           array with values
c nval          dimension of val
c dis           distance of isolines (dis>0)
c mode		0: isolines with dis  1: isolines in /isolin/
c		2: color  3: color on element values

	use basin
	use color

	implicit none

c argument
	integer nval
	real regval(nval)
	real regpar(7)
	real dis
	integer mode
c local
	character*80 line
	real x(4),y(4),z(4)
	real dist,flag
	logical bexreg
	integer ie,ii,kn
	integer inull
	integer icsave
	integer nx,ny,ix,iy
	integer ierr
	real x0,y0,dx,dy
	real, allocatable :: femval(:)

        real getpar
c       logical is_r_nan

c--------------------------------------------------------------------
c initialization
c--------------------------------------------------------------------

	bexreg = nint(getpar('iexreg')) > 0	!plot in half valid box

c--------------------------------------------------------------------
c find isolines to be plotted
c--------------------------------------------------------------------

	if( mode .eq. 0 ) then
		dist = dis
	else
		dist = 0.
	end if

	nx = nint(regpar(1))
	ny = nint(regpar(2))
	x0 = regpar(3)
	y0 = regpar(4)
	dx = regpar(5)
	dy = regpar(6)
	flag = regpar(7)

	fnull = flag

c--------------------------------------------------------------------
c loop over elements
c--------------------------------------------------------------------

	!write(6,*) 'isoline: isoanz... ',isoanz,mode

	if( mode .eq. 2 ) then

c	  -----------------------------------------
c	  color plot
c	  -----------------------------------------

	  call qcomm('plotting regular grid')
	  call get_color_table(icsave)
	  !call set_color_table(-1)

	  do iy=2,ny
	    do ix=2,nx
	      call set_box_val(nx,ny,ix,iy,x0,y0,dx,dy,regval,x,y,z)
	      if( bexreg ) call extend_box_val(z,flag)	!plot in half valid box
	      call plot_box_val(x,y,z,ciso,fiso,isoanz+1,fnull)
	    end do
	  end do

	  call set_color_table(icsave)
	  call qcomm('finished plotting regular grid')

	else if( mode .eq. 4 ) then		!interpolate to fem grid


	end if

	end

c***************************************************************

	subroutine extend_box_val(z,flag)

	implicit none

	real z(4)
	real flag

	integer ic,j
	real tot

	ic = count( z == flag )

	if( ic == 0 .or. ic == 4 ) return

	tot = 0.
	do j=1,4
	  if( z(j) /= flag ) tot = tot + z(j)
	end do
	tot = tot / (4-ic)

	where( z == flag ) z = tot

	end

c***************************************************************

	subroutine set_box_val(nx,ny,ix,iy,x0,y0,dx,dy,val,x,y,z)

	implicit none

	integer nx,ny,ix,iy
	real x0,y0,dx,dy
	real val(nx,ny)
	real x(4),y(4),z(4)

	real xl,yl,xh,yh

	xl = x0 + (ix-2)*dx
	yl = y0 + (iy-2)*dy
	xh = x0 + (ix-1)*dx
	yh = y0 + (iy-1)*dy

	x(1) = xl
	x(2) = xh
	x(3) = xh
	x(4) = xl
	y(1) = yl
	y(2) = yl
	y(3) = yh
	y(4) = yh

	z(1) = val(ix-1,iy-1)
	z(2) = val(ix,iy-1)
	z(3) = val(ix,iy)
	z(4) = val(ix-1,iy)

	end

c***************************************************************

	subroutine  pldif(f,x,y,dist,flag)

c plots isolines between different regions

	implicit none

c arguments
	real f(3)		!values at vertices
	real x(3)		!x-coordinates of vertices
	real y(3)		!y-coordinates of vertices
	real dist		!distance between isolines
	real flag		!null value -> do not plot

c local
	integer ii,i1
	real xm,ym,xc,yc
c save

	xm = 0.
	ym = 0.
	do ii=1,3
	  xm = xm + x(ii)
	  ym = ym + y(ii)
	end do
	xm = xm / 3.
	ym = ym / 3.

	do ii=1,3
	  i1 = mod(ii,3) + 1
	  if( f(ii) .ne. flag .and. f(i1) .ne. flag ) then
	    if( abs(f(ii)-f(i1)) .gt. dist ) then
	      xc = 0.5*(x(ii)+x(i1))
	      yc = 0.5*(y(ii)+y(i1))
	      call qmove(xc,yc)
	      call qplot(xm,ym)
	    end if
	  end if
	end do

	end

c***************************************************************

	subroutine pliso(f,x,y,dist,isoanz,fnull,fiso)

c plots isolines

	implicit none

c arguments
	real f(3)		!values at vertices
	real x(3)		!x-coordinates of vertices
	real y(3)		!y-coordinates of vertices
	real dist		!distance between isolines
	integer isoanz		!number of isolines in fiso
	real fnull		!null value -> do not plot
	real fiso(isoanz)	!values of isolines
c local
	real fmin,fmax
	real fisol
	real xp(2),yp(2)
	real df,fac
	integer i,iisol,ip
	integer is,iss
	integer i1,i2
	integer inull,isum
c save
	real epsp
	save epsp
	data epsp /1.e-10/

c get min/max

	fmin = min(f(1),f(2),f(3))
	fmax = max(f(1),f(2),f(3))

c look for null value (if just one null value, we plot border)

	inull = 0
	isum = 0
	do i=1,3
	  if( f(i) .eq. fnull ) then
		inull = inull + 1
		isum = isum + i
	  end if
	end do

	if( inull .eq. 1 ) then
	  i1 = mod(isum,3) + 1
	  i2 = mod(i1,3) + 1
	  call qmove(x(i1),y(i1))
	  call qplot(x(i2),y(i2))
	end if

	if( inull .gt. 0 ) return
	  
c if constant no processing

	if( fmin .eq. fmax ) return

c starting point

	iisol = 0
	if( dist .gt. 0. ) then
		fisol = int(fmin/dist)*dist
	else
		i = 1
		do while( i .le. isoanz .and. fiso(i) .lt. fmin )
			i = i + 1
		end do
		if( i .gt. isoanz ) return
		iisol = i
		fisol = fiso(i)
	end if

c loop over isolines

	do while( fisol .lt. fmax )

	is=0
	iss=0
	do i=1,3
	  if( fisol .eq. f(i) ) then
	    is = is + 1
	    iss = iss + i
	  end if
	end do

	if( is .eq. 3 ) then		!constant
		return
	else if( is .eq. 2 ) then	!two vertices equal 
		i = 6 - iss
		i1 = mod(i,3) + 1
		i2 = mod(i1,3) + 1
		call qmove(x(i1),y(i1))
		call qplot(x(i2),y(i2))
	else
		ip = 1
		if( is .eq. 1 ) then	!this is first point
			xp(1) = x(iss)
			yp(1) = y(iss)
			ip = 2
		end if

		i1 = 3
		i2 = 1
		do while( ip .le. 2 .and. i2 .le. 3 )
		  df = f(i2) - f(i1)
		  if( abs(df) .gt. epsp ) then
		    fac = (fisol-f(i1)) / df
		    if( fac .gt. 0. .and. fac .lt. 1. ) then
        	      xp(ip) = (x(i2)-x(i1))*fac + x(i1)
        	      yp(ip) = (y(i2)-y(i1))*fac + y(i1)
		      ip = ip + 1
		    end if
		  end if
		  i1 = i2
		  i2 = i2 + 1
		end do

		if( ip .eq. 3 ) then
		  call qmove(xp(1),yp(1))
		  call qplot(xp(2),yp(2))
		end if
	end if

c new isoline

	if( dist .gt. 0. ) then
	  fisol = fisol + dist
	else
	  iisol = iisol + 1
	  if( iisol .gt. isoanz ) then
	    fisol = fmax + 1.
	  else
	    fisol = fiso(iisol)
	  end if
	end if

	end do	!do while over isolines

c end of loop

	return
	end

c*********************************************************************

	subroutine plcol(x,y,z,color,rlev,ncol,fnull)

c interpolates and plots colors in triangle
c
c x,y,z		coordinates and values of vertices in triangle
c color		array of colors to use (in total ncol colors)
c rlev		levels to use (in total ncol-1 levels)
c ncol		total number of colors in color
c fnull		null value -> do not interpolate

	implicit none

	integer ncol
	real fnull
	real x(3),y(3),z(3)
	real color(ncol),rlev(ncol-1)

	logical btypa,bdiff,bnull
	real xp(7),yp(7)
        integer idebug
	integer i,ii,i1,i2,imin
	integer icont
	integer ieck,icol
	integer ipf,ipb,ipa,ip
	real zmin,zmax
	real xa,ya,za
	real xb,yb,zb
	real x1,y1,z1
	real x2,y2,z2
	real zr,dz

c       logical is_r_nan

	real eps
	save eps
	data eps /1.e-5/	!befor eps=0	!$$ALPHA - ERROR

c find min/max %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	call mima(z,3,zmin,zmax)

        idebug = 0

c see if null value %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	bnull = .false.
	do i=1,3
	  if( z(i) .eq. fnull ) bnull = .true.
	end do

	if( bnull ) return

c put nodes in order %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c if there is more then one node with z == zmin ==> btypa = .true.
c (corresponds to one free (high) node)
c the first low node is called (xa,ya), the second node (xb,yb)
c the free node is called (x1,y1)
c
c if there is only one node with z = zmin ==> btypa = .false.
c this node is called (xa,ya)
c a dummy node (xb,yb) is created with (xa,ya) = (xb,yb)
c the free nodes are called (x1,y1) and (x2,y2)
c
c in all cases the following is true:
c
c za == zb <= z1 <= z2
c
c for btypa == .true.  => (xa,ya) != (xb,yb) and there is no (x2,y2)
c for btypa == .false. => (xa,ya) == (xb,yb) and (x2,y2) exists
c
c avoid compiler warnings

	xa = 0.
	ya = 0.
	za = 0.
	xb = 0.
	yb = 0.
	zb = 0.
	x1 = 0.
	y1 = 0.
	z1 = 0.
	x2 = 0.
	y2 = 0.
	z2 = 0.
	btypa = .false.

c find first and last node

	bdiff = .false.
	icont=0
	ieck=0
	do ii=3,1,-1
	if(z(ii).eq.zmin) then
	    if(icont.eq.0) then
		xa=x(ii)
		ya=y(ii)
		za=z(ii)
		icont=icont+1
		ieck=ieck+ii
	    else if(icont.eq.1) then
		xb=x(ii)
		yb=y(ii)
		zb=z(ii)
		icont=icont+1
		ieck=ieck+ii
	    else if(icont.eq.2) then
		icont=icont+1
	    end if
	end if
	end do

c find free nodes

	if(icont.eq.1) then
		xb=xa
		yb=ya
		zb=za
		i1=mod(ieck,3)+1
		i2=mod(i1,3)+1
		if( z(i1) .gt. z(i2) ) then
			ii = i1
			i1 = i2
			i2 = ii
		end if
		x1=x(i1)
		y1=y(i1)
		z1=z(i1)
		x2=x(i2)
		y2=y(i2)
		z2=z(i2)
		btypa=.false.
	else if(icont.eq.2) then
		ii=6-ieck
		x1=x(ii)
		y1=y(ii)
		z1=z(ii)
		btypa=.true.
	else if(icont.eq.3) then
		ii=1
		x1=x(ii)
		y1=y(ii)
		z1=z(ii)
		btypa=.true.
	end if

c find first level .ge. nodes a and b %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	do icol=1,ncol-1
c		if(rlev(icol).ge.zmin) goto 2
		if(rlev(icol).gt.zmin) goto 2
	end do
    2	continue
	imin = icol

c start plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	do icol = imin,ncol

c zr is actual level (last level must be always true)

	if(icol.lt.ncol) then
		zr=rlev(icol)
	else
		zr=zmax+1.
	end if

c two cases

	if(btypa) then

c case a (one free node)

		if(z1.le.zr) then	!case a1, fill remaining triangle
			xp(1)=xa
			yp(1)=ya
			xp(2)=x1
			yp(2)=y1
			xp(3)=xb
			yp(3)=yb
			ip=3
                        call fill_area(idebug,"case a1"
     +                                  ,ip,xp,yp,color(icol))
			goto 99
		else			!case a2, one free node
			xp(1)=xa
			yp(1)=ya
			dz=z1-za
			if(dz.gt.eps) then
				xp(2)=xa+(x1-xa)*(zr-za)/dz
				yp(2)=ya+(y1-ya)*(zr-za)/dz
			else
				xp(2)=(xa+x1)*0.5
				yp(2)=(ya+y1)*0.5
			end if
			dz=z1-zb
			if(dz.gt.eps) then
				xp(3)=xb+(x1-xb)*(zr-zb)/dz
				yp(3)=yb+(y1-yb)*(zr-zb)/dz
			else
				xp(3)=(xb+x1)*0.5
				yp(3)=(yb+y1)*0.5
			end if
			xp(4)=xb
			yp(4)=yb
			ip=4
                        call fill_area(idebug,"case a2"
     +                                  ,ip,xp,yp,color(icol))
			xa=xp(2)
			ya=yp(2)
			za=zr
			xb=xp(3)
			yb=yp(3)
			zb=zr
		end if
	else

c case b (two free nodes)

		if( z2 .le. zr ) then	!case b1, fill remaining triangle

			xp(1)=xa
			yp(1)=ya
			xp(2)=x1
			yp(2)=y1
			xp(3)=x2
			yp(3)=y2
			ip = 3
			if( bdiff ) then
				xp(4)=xb
				yp(4)=yb
				ip = 4
			end if
                        call fill_area(idebug,"case b1"
     +                                  ,ip,xp,yp,color(icol))

			goto 99

		else if( z1 .gt. zr ) then !case b2, still two free nodes

			ip = 1
			xp(ip)=xa
			yp(ip)=ya

			ip=ip+1
			dz=z1-za
			if(dz.gt.eps) then
				xp(ip)=xa+(x1-xa)*(zr-za)/dz
				yp(ip)=ya+(y1-ya)*(zr-za)/dz
			else
				xp(ip)=(xa+x1)*0.5
				yp(ip)=(ya+y1)*0.5
			end if

			ip=ip+1
			dz=z2-zb
			if(dz.gt.eps) then
				xp(ip)=xb+(x2-xb)*(zr-zb)/dz
				yp(ip)=yb+(y2-yb)*(zr-zb)/dz
			else
				xp(ip)=(xb+x2)*0.5
				yp(ip)=(yb+y2)*0.5
			end if

			if( bdiff ) then
				ip = ip + 1
				xp(ip)=xb
				yp(ip)=yb
			end if

                        call fill_area(idebug,"case b2"
     +                                  ,ip,xp,yp,color(icol))

			xa = xp(2)
			ya = yp(2)
			za = zr
			xb = xp(3)
			yb = yp(3)
			zb = zr
			bdiff = .true.

		else				!case b3, now one free node

			ip = 1
			xp(ip)=xa
			yp(ip)=ya

			ip=ip+1
			xp(ip)=x1
			yp(ip)=y1

			ip=ip+1
			dz=z2-z1
			if(dz.gt.eps) then
				xp(ip)=x1+(x2-x1)*(zr-z1)/dz
				yp(ip)=y1+(y2-y1)*(zr-z1)/dz
			else
				xp(ip)=(x2+x1)*0.5
				yp(ip)=(y2+y1)*0.5
			end if

			ip=ip+1
			dz=z2-zb
			if(dz.gt.eps) then
				xp(ip)=xb+(x2-xb)*(zr-zb)/dz
				yp(ip)=yb+(y2-yb)*(zr-zb)/dz
			else
				xp(ip)=(xb+x2)*0.5
				yp(ip)=(yb+y2)*0.5
			end if

			if( bdiff ) then
				ip = ip + 1
				xp(ip)=xb
				yp(ip)=yb
			end if

                        call fill_area(idebug,"case b3"
     +                                  ,ip,xp,yp,color(icol))

			xa = xp(3)
			ya = yp(3)
			za = zr
			xb = xp(4)
			yb = yp(4)
			zb = zr
			x1=x2
			y1=y2
			z1=z2
			bdiff = .true.
			btypa = .true.

		end if
	end if

	end do

   99	continue

	return
	end

c***************************************************************

	subroutine plnode(x,y,z,color,rlev,ncol,fnull)

c plots colors of nodal values without interpolation
c
c x,y,z		coordinates and values of vertices in triangle
c color		array of colors to use (in total ncol colors)
c rlev		levels to use (in total ncol-1 levels)
c ncol		total number of colors in color
c fnull		null value -> do not interpolate

	implicit none

	integer ncol
	real fnull
	real x(3),y(3),z(3)
	real color(1),rlev(1)

	integer idebug,ip,ii
	real xp(4),yp(4)
	real col

	real get_color

	idebug = 0
	ip = 4		!always 4 nodes

	do ii=1,3
	  if( z(ii) .ne. fnull ) then
	    call make_xy(ii,x,y,xp,yp)
	    col = get_color(z(ii),ncol,color,rlev)
            call fill_area(idebug,"no_intp",ip,xp,yp,col)
	  end if
	end do

	end

c***************************************************************

	subroutine plot_box_val(x,y,z,ciso,fiso,ncol,fnull)

	implicit none

	integer ncol
	real fnull
	real x(4),y(4),z(4)
	real ciso(ncol),fiso(ncol-1)

	real xx(5),yy(5),zz(5)
	real xxx(3),yyy(3),zzz(3)

	xx(1:4) = x
	yy(1:4) = y
	zz(1:4) = z
	xx(5) = sum(x)/4.
	yy(5) = sum(y)/4.
	zz(5) = sum(z)/4.

	if( count( zz == fnull ) > 0 ) zz = fnull

	call set_tria_val(xx,yy,zz,5,1,2,xxx,yyy,zzz)
	call plcol(xxx,yyy,zzz,ciso,fiso,ncol,fnull)
	call set_tria_val(xx,yy,zz,5,2,3,xxx,yyy,zzz)
	call plcol(xxx,yyy,zzz,ciso,fiso,ncol,fnull)
	call set_tria_val(xx,yy,zz,5,3,4,xxx,yyy,zzz)
	call plcol(xxx,yyy,zzz,ciso,fiso,ncol,fnull)
	call set_tria_val(xx,yy,zz,5,4,1,xxx,yyy,zzz)
	call plcol(xxx,yyy,zzz,ciso,fiso,ncol,fnull)

	end

c***************************************************************

	subroutine set_tria_val(xx,yy,zz,i1,i2,i3,xxx,yyy,zzz)

	implicit none

	integer i1,i2,i3
	real xx(5),yy(5),zz(5)
	real xxx(3),yyy(3),zzz(3)

	xxx(1) = xx(i1)
	xxx(2) = xx(i2)
	xxx(3) = xx(i3)
	yyy(1) = yy(i1)
	yyy(2) = yy(i2)
	yyy(3) = yy(i3)
	zzz(1) = zz(i1)
	zzz(2) = zz(i2)
	zzz(3) = zz(i3)

	end

c***************************************************************

	function get_color(z,ncol,color,rlev)

c returns color for value z

	implicit none

	real get_color
	real z
	integer ncol
	real color(ncol),rlev(ncol-1)

	integer i

	do i=1,ncol-1
	  if( rlev(i) .gt. z ) goto 1
	end do
    1	continue

	get_color = color(i)

	end

c***************************************************************

        subroutine make_xy(in,x,y,xp,yp)

c makes x/y of finite volume of local node ii

	implicit none

	integer in
	real x(3),y(3)
	real xp(4),yp(4)

	integer ii,ib,ia
	real xc,yc

	xc = 0.
	yc = 0.
	do ii=1,3
	  xc = xc + x(ii)
	  yc = yc + y(ii)
	end do
	xc = xc / 3.
	yc = yc / 3.

	ib = mod(in+1,3)+1
	ia = mod(in,3)+1
	xp(1) = x(in)
	yp(1) = y(in)
	xp(2) = 0.5*(x(in)+x(ia))
	yp(2) = 0.5*(y(in)+y(ia))
	xp(3) = xc
	yp(3) = yc
	xp(4) = 0.5*(x(in)+x(ib))
	yp(4) = 0.5*(y(in)+y(ib))

	end

c***************************************************************

        subroutine fill_area(idebug,text,ip,xp,yp,col)

c fills area of polygon given by xp,yp with color

        implicit none

        integer idebug
        character*(*) text
        integer ip
        real xp(ip)
        real yp(ip)
        real col

        if( idebug .gt. 999 ) then
          write(6,*) 'idebug ---------------------------------'
          write(6,*) idebug,ip,xp,yp,col,text
          write(6,*) '----------------------------------------'
	  call qcomm(text)
        end if

        if( col .lt. 0. ) return

	call qsetc(col)
	call qafill(ip,xp,yp)

        end

c***************************************************************

	subroutine set_fxy_vals(ie,flag,val,f,x,y,inull)

c returns x,y,val of element, and indication of flag values

	use basin

        implicit none

        integer ie              !element for which info is needed
        real flag               !value of flag
        real val(nkn)             !nodal value to be extracted
        real f(3),x(3),y(3)     !return values of val,x,y
        integer inull           !total number of flags found

	integer ii,kn

	inull = 0

	do ii=1,3
	    kn=nen3v(ii,ie)
	    f(ii)=val(kn)
	    if( f(ii) .eq. flag ) inull = inull + 1
	    x(ii)=xgv(kn)
	    y(ii)=ygv(kn)
	end do

	end

c***************************************************************

	subroutine set_xy(ie,x,y)

c returns x,y of element

	use basin

        implicit none

        integer ie              !element for which info is needed
        real x(3),y(3)          !return values of x,y

	integer ii,kn

	do ii=1,3
	    kn=nen3v(ii,ie)
	    x(ii)=xgv(kn)
	    y(ii)=ygv(kn)
	end do

	end

c***************************************************************

