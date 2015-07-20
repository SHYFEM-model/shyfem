c
c $Id: subproj.f,v 1.5 2010-03-11 15:36:38 georg Exp $
c
c subroutines for handling projection
c
c revision log :
c
c 17.11.2011    ggu     written from scratch
c 10.01.2012    ggu     bug fix: c_param was real
c 21.01.2015    ggu     code to handle projection in both directions
c
c****************************************************************            

	subroutine handle_projection

c handles projection

	implicit none

	logical bspheric

	logical is_spherical

	integer icall
	save icall
        data icall / 0 /
        
        if( icall .ne. 0 ) return

	icall = 1
       
c---------------------------------------------------------------
c initialization
c---------------------------------------------------------------

	bspheric = is_spherical()

	if( bspheric ) then	!lat/lon -> cartesian
	  call proj_geo2cart
	else			!cartesian -> lat/lon
	  call proj_cart2geo
	end if

	end

c****************************************************************            

	subroutine proj_cart2geo

c handles projection - converts x/y to lat/lon

	use mod_tides
	use basin

	implicit none

        include 'param.h'

	integer mode,iproj,i
	double precision c_param(9)

        real getpar

	integer icall
	save icall
        data icall / 0 /
        
        if( icall .ne. 0 ) return

	icall = 1
       
c---------------------------------------------------------------
c initialization
c---------------------------------------------------------------

	mode = 1	!from cartesian to lat/lon

        iproj = nint(getpar('iproj'))

	if( iproj .eq. 0 ) then
	  !nothing
	else if( iproj .eq. 1 ) then
	  c_param(1) = getpar('c_fuse')
	  c_param(2) = getpar('c_x0')
	  c_param(3) = getpar('c_y0')
	else if( iproj .eq. 2 ) then
	  c_param(1) = getpar('c_zone')
	  c_param(2) = getpar('c_x0')
	  c_param(3) = getpar('c_y0')
	else if( iproj .eq. 3 ) then
	  c_param(1) = getpar('c_phi')
	  c_param(2) = getpar('c_lon0')
	  c_param(3) = getpar('c_lat0')
	else if( iproj .eq. 4 ) then
	  c_param(1) = getpar('c_lamb')
	  c_param(2) = getpar('c_x0')
	  c_param(3) = getpar('c_y0')
	  c_param(4) = getpar('c_skal')
	else
	  write(6,*) 'iproj = ',iproj
	  stop 'error stop proj_cart2geo: value for iproj not allowed'
	end if

	call init_coords(iproj,c_param)
	call convert_coords(mode,nkn,xgv,ygv,xgeov,ygeov)
	xcartv = xgv
	ycartv = ygv

	write(6,*) 'start of proj_cart2geo'
	write(6,*) 'mode  = ',mode
	write(6,*) 'iproj = ',iproj

	write(6,*) (xgv(i),i=1,5)
	write(6,*) (ygv(i),i=1,5)
	write(6,*) (xgeov(i),i=1,5)
	write(6,*) (ygeov(i),i=1,5)
	write(6,*) (xcartv(i),i=1,5)
	write(6,*) (ycartv(i),i=1,5)

	write(6,*) 'end of proj_cart2geo'

	end

c****************************************************************            

	subroutine proj_geo2cart

c handles projection - converts lat/lon to x/y

	use mod_tides
	use basin

	implicit none

        include 'param.h'

	integer mode,iproj,i,k
	double precision c_param(9)
	double precision c_lat0,c_lon0,c_phi
	real xmin,ymin,xmax,ymax

	integer icall
	save icall
        data icall / 0 /
        
        if( icall .ne. 0 ) return

	icall = 1

	mode = -1		!lat/lon to x/y
        iproj = 3		!always use cpp

	xmin = xgv(1)
	ymin = ygv(1)
	xmax = xgv(1)
	ymax = ygv(1)
	do k=1,nkn
	  xmin = min(xmin,xgv(k))
	  ymin = min(ymin,ygv(k))
	  xmax = max(xmax,xgv(k))
	  ymax = max(ymax,ygv(k))
	end do

	c_phi  = 0.5*(ymax-ymin)
	c_lat0 = 0.5*(ymax-ymin)
	c_lon0 = 0.5*(xmax-xmin)

	c_param(1) = c_phi
	c_param(2) = c_lon0
	c_param(3) = c_lat0

	call init_coords(iproj,c_param)
	call convert_coords(mode,nkn,xcartv,ycartv,xgv,ygv)
	xgeov = xgv
	ygeov = ygv

	write(6,*) 'start of proj_geo2cart'
	write(6,*) 'mode  = ',mode
	write(6,*) 'iproj = ',iproj

	write(6,*) (xgv(i),i=1,5)
	write(6,*) (ygv(i),i=1,5)
	write(6,*) (xgeov(i),i=1,5)
	write(6,*) (ygeov(i),i=1,5)
	write(6,*) (xcartv(i),i=1,5)
	write(6,*) (ycartv(i),i=1,5)

	write(6,*) 'end of proj_geo2cart'

	end

c****************************************************************            

        subroutine baric_cart(ie,x,y)

c finds baricentre of element
c
c ie            number of element
c x,y           coordinates of baricentre (return value)

	use mod_tides
	use basin

        implicit none

c arguments
        integer ie
        real x,y
c common blocks
c local variables
        integer i,kkk
        real xb,yb

        include 'param.h'

        xb=0.
        yb=0.
        do i=1,3
           kkk=nen3v(i,ie)
           xb=xb+xcartv(kkk)
           yb=yb+ycartv(kkk)
        end do

        x=xb/3.
        y=yb/3.

        end

c****************************************************************            

        subroutine getexy_cart(ie,x,y)

c gets coordinates x/y for element ie

	use mod_tides
	use basin

        implicit none

        integer ie
        real x(3), y(3)

        integer k,ii

        include 'param.h'

        do ii=1,3
          k = nen3v(ii,ie)
          x(ii) = xcartv(k)
          y(ii) = ycartv(k)
        end do

        end

c****************************************************************            


