c
c $Id: supbas.f,v 1.12 2010-02-26 15:29:19 georg Exp $
c
c routines for plotting basin
c
c contents :
c
c subroutine basinit			internal initialization
c
c subroutine bash(mode)			plots basin (shell)
c
c subroutine basin(mode)		plots basin
c subroutine bpixel			writes the pixel size as a comment 
c subroutine frame(mode)		plots frame
c subroutine basmima(mode)		computes min/max of basin
c subroutine setbas(x0,y0,x1,y1)	sets min/max of basin to plot
c subroutine getbas(x0,y0,x1,y1)	gets min/max of basin to plot
c subroutine boundline			plots boundary line for lagoon
c
c subroutine reggrid(ngrid,dist,gray)	plots regular grid
c
c subroutine basscale(dist)		computes typical length scale for basin
c
c function roundm(r,mode)		rounds r to the closest value
c function rround(r,rmaster,mode)	rounds r to next rmaster value
c function rdist(xmin,ymin,xmax,ymax)	computes gridspacing 
c function divdist(x,n,mode)		divides x into n equal pieces
c
c subroutine handle_spherical		handles spherical coordinates
c subroutine plot_reg_grid		handles plotting of regular grid
c subroutine label_reg_grid		handles labeling of regular grid
c
c revision log :
c
c 16.02.1999  ggu     bpixel: write bounding box to ps file (as comment)
c 09.02.2000  ggu     use inboxdim to compute box to plot
c 12.06.2000  ggu     get gray value for basin mode 3 from str file
c 11.02.2001  ggu     routine to compute typical length scale
c 21.08.2003  ggu     occupy is called elsewhere
c 16.12.2004  ggu     changed reggrid to plot regular grid
c 02.03.2005  ggu     in bash: bug fix -> get size of grid if not given
c 12.06.2009  ggu     new routines to handle spherical coords. & regular grid
c 15.06.2009  ggu     call to reggrid() changeed -> pass in gray value
c 14.09.2009  ggu     new routine divdist()
c 22.02.2010  ggu     new routine bw_frame() to plot bw scale around plot
c 09.04.2010  ggu     bug fix in frac_pos() -> maybe compiler error
c 17.05.2011  ggu     new routine basin_number()
c 30.08.2012  ggu     new routines to automatically label spherical grid
c 24.10.2012  ggu     bug in labelling non spherical grid (returned -1)
c 02.05.2013  ggu     handle fact in spherical coords
c 02.05.2013  ggu     meteo point plotting (plot_meteo_points())
c 13.06.2013  ggu     bug fix in spherical_fact() -> set fact to 1
c 13.12.2013  ggu     new mode=4 for plotting gray grid over scalar variable
c 30.05.2014  ggu     new metpnt for meteo points, imicro computed
c 10.02.2015  ggu     also plot other points, also regular points
c 12.10.2015  ggu     fix in rround() to handle rmaster==0 case
c
c notes:
c
c rules for regular frame around plot:
c
c if reggrd is given		-> use this value
c if reggrd == 0		-> do not plot frame
c else (reggrd==-1)
c if no legend is requested	-> do not plot frame
c if legend is requested
c	if spherical		-> plot frame and no scale/north
c	else			-> plot scale/north and no frame
c
c*************************************************************

	subroutine basinit

c internal initialization

	implicit none

	real xmin,ymin,xmax,ymax
	common /bamima/ xmin,ymin,xmax,ymax
	save /bamima/

	logical binit
	save binit
	data binit /.false./

	if( binit ) return

	binit = .true.

	call basmima(0)		!compute exact dimensions as default
	call bastlscale		!computes typical length scale

	end

c**************************************************************

	subroutine bash(mode)

! plots basin (shell)
!
! mode	
!	0: only scaling  
!	1: net  
!	2: boundary and legend
!	3: net in gray (for bathymetry - use bgray)
!	4: net in gray (for scalar and velocities - use bsgray)
!
! call first with mode == 0 (scaling9
! then call with desired mode
! must be called with mode == 2 before closing the plot

	implicit none

	integer mode

	integer ifreg
	real x0,y0,x1,y1
	real x0leg,y0leg,x1leg,y1leg
	real dxygrd,x,y
	real cgray
	character*80 bndlin,file
	real getpar
	logical inboxdim
        logical is_spherical
	logical is_box_given

!--------------------------------------
! initializing
!--------------------------------------

	cgray = getpar('cislnd')	!gray scale for islands

	call basinit

!--------------------------------------
! if mode == 0 -> do scaling
!--------------------------------------

	if( mode .eq. 0 ) then

	  call basmima(0)
	  if( inboxdim(' ',x0,y0,x1,y1) ) then	!size of plotting is given
	    call setbas(x0,y0,x1,y1)
	  end if

	  call getbas(x0,y0,x1,y1)		!bug fix 2.3.2005

	  !--------------------------------------
	  ! prepare regular grid
	  !--------------------------------------

	  dxygrd = getpar('dxygrd')
	  x = x0 + dxygrd/2.
	  y = y0 + dxygrd/2.
	  call setgeo(x,y,dxygrd,dxygrd,-999.)

	  call annote		!annotation
	  call plot_basin(0)	!scaling
	  call label_reg_grid
	  call plot_islands(cgray)

	  return

	end if

	! frame -> old position, now moved down, see if it workes, 19.11.2014
	!call frame(0)
	!call plot_reg_grid

!--------------------------------------
! plot
!--------------------------------------

        call qcomm('Plotting basin')
        call plot_basin(mode)

!--------------------------------------
! end of plot - now only labeling
!--------------------------------------

	if( mode .ne. 2 ) return

	!--------------------------------------
	! boundary line and islands
	!--------------------------------------

	call boundline

	!--------------------------------------
	! frame
	!--------------------------------------

	call frame(0)
	call plot_reg_grid

	!--------------------------------------
	! user defined legend
	!--------------------------------------

	call legplo

	!--------------------------------------
	! legend (north and scale)
	!--------------------------------------

	!if( inboxdim('leg',x0,y0,x1,y1) ) then	!write legend
	if( is_box_given('leg') ) then	!write legend
          if( is_spherical() ) then
            write(6,*) 'coordinates are spherical'
            write(6,*) 'no north and scale written...'
            return
	  else
	    if( inboxdim('leg',x0leg,y0leg,x1leg,y1leg) ) then
	      call scale_legend(x0leg,y0leg,x1leg,y1leg)
	    end if
          end if
	end if

	!--------------------------------------
	! special output
	!--------------------------------------

	call getfnm('metpnt',file)
	if( file .ne. ' ' ) then
	  call plot_meteo_points(1,file)
	end if

	call getfnm('obspnt',file)
	if( file .ne. ' ' ) then
	  call plot_meteo_points(2,file)
	end if

	ifreg = nint(getpar('ifreg'))
	if( ifreg .gt. 0 ) then
	  call plot_regular_points
	end if

!--------------------------------------
! end of routine
!--------------------------------------

	end

c*************************************************************

	subroutine plot_basin(mode)

c plots basin
c
c mode	
c	0: only scaling  
c	1: net  
c	2: boundary  
c	3: net in gray (for bathymetry - use bgray)
c	4: net in gray (for scalar and velocities - use bsgray)

	use mod_geom
	use basin

	implicit none

	integer mode

	include 'param.h'

	real xmin,ymin,xmax,ymax
	common /bamima/ xmin,ymin,xmax,ymax

	integer ie,kn,ii,k
	real gray

	real getpar

	call basinit

	if( mode .eq. 0 ) then				!scaling
	  call qworld(xmin,ymin,xmax,ymax)
	  call qrcfy
	  call handle_spherical
	  call bpixel
	else if( mode .eq. 2 ) then			!plot boundary
	  call qgray(0.)
	  do k=1,nkn
	    if( kantv(1,k) .ne. 0 ) then
	      kn=kantv(1,k)
	      if( kn .gt. k ) call qline(xgv(k),ygv(k),xgv(kn),ygv(kn))
	      kn=kantv(2,k)
	      if( kn .gt. k ) call qline(xgv(k),ygv(k),xgv(kn),ygv(kn))
	    end if
	  end do
	else if( mode .gt. 0 .and. mode .le. 4 ) then	!plot grid
	  call qgray(0.)
	  if( mode .eq. 3 ) then
	    gray = getpar('bgray')
	    call qgray(gray)
	  end if
	  if( mode .eq. 4 ) then
	    gray = getpar('bsgray')
	    if( gray .lt. 0. ) return	!do not plot net
	    call qgray(gray)
	  end if

	  do ie=1,nel
	    kn = nen3v(3,ie)
	    call qmove(xgv(kn),ygv(kn))
	    do ii=1,3
	      kn = nen3v(ii,ie)
	      call qplot(xgv(kn),ygv(kn))
	    end do
	  end do
	else
	  write(6,*) 'mode = ',mode
	  stop 'error stop basin: internal error'
	end if

	call qgray(0.)

	end

c*************************************************************

	subroutine basin_number(mode)

c plots basin with node and element numbers
c
c mode	1: node number   2: element number
c	positive: external    negative: internal

	use mod_geom
	use basin

	implicit none

	integer mode

	include 'param.h'

	real xmin,ymin,xmax,ymax
	common /bamima/ xmin,ymin,xmax,ymax

	integer ie,kn,ii,k,ke,iee,i
	real gray,x,y
	character*10 s

	integer ialfa,ipext,ieext

	write(6,*) '======================================='
	write(6,*) '======================================='
	write(6,*) '======================================='
	write(6,*) 'debug print of node and element numbers'
	write(6,*) 'please remove in production code'
	write(6,*) '======================================='
	write(6,*) '======================================='
	write(6,*) '======================================='

	call qgray(0.)
	call qtxts(12)

	if( abs(mode) .eq. 1 ) then

	  do k=1,nkn
	    ke = k
	    if( mode .gt. 0 ) ke = ipext(k)
	    x = xgv(k)
	    y = ygv(k)
	    i = ialfa(float(ke),s,-1,-1)
	    call qtext(x,y,s)
	  end do

	else if( abs(mode) .eq. 2 ) then

	  call qtxtcc(0,0)
	  do ie=1,nel
	    x = 0.
	    y = 0.
	    do ii=1,3
	      k = nen3v(ii,ie)
	      x = x + xgv(k)
	      y = y + ygv(k)
	    end do
	    x = x / 3.
	    y = y / 3.
	    iee = ie
	    if( mode .gt. 0 ) iee = ieext(ie)
	    i = ialfa(float(iee),s,-1,-1)
	    call qtext(x,y,s)
	  end do
	  call qtxtcc(-1,-1)
	else

	  write(6,*) mode
	  stop 'error stop basin_number: mode'

	end if

	end

c*****************************************************************

	subroutine bpixel

c writes the pixel size as a comment to PS file

	implicit none

	character*80 line
	real fact
	real xvmin,yvmin,xvmax,yvmax
	real dxy
	integer idxy
	integer ixmin,iymin,ixmax,iymax

	call qgetvp(xvmin,yvmin,xvmax,yvmax)

	fact = 28.3	!factor to translate from PS dimensions to pixel
	dxy = 1.	!translation of plot
	idxy = 5	!leave some space around plot

	ixmin = fact * ( xvmin + dxy ) - idxy
	iymin = fact * ( yvmin + dxy ) - idxy
	ixmax = fact * ( xvmax + dxy ) + idxy
	iymax = fact * ( yvmax + dxy ) + idxy

c	write(6,*) 'message from basin() ...'
c	write(6,*) '****** ',xvmin,yvmin,xvmax,yvmax
c	write(6,*) '****** ',ixmin,iymin,ixmax,iymax

	write(line,'(a,4i8)') 'InternalBoundingBox'
     +				,ixmin,iymin,ixmax,iymax
	call qcomm(line)

	end

c*****************************************************************

	subroutine frame(mode)

c plots frame
c
c 0: just box around plot (only mode allowed)

	implicit none

	integer mode
	integer iaux
	real pxareg,pyareg,pxereg,pyereg
	real plon0,plat0,dlon,dlat
	real alonmin,alatmin,alonmax,alatmax
	real alon,alat
	real x,y

	real xmin,ymin,xmax,ymax
	common /bamima/ xmin,ymin,xmax,ymax

	call basinit

	pxareg = xmin
	pyareg = ymin
	pxereg = xmax
	pyereg = ymax

	call qcomm('Plotting frame around basin')
	call qgray(0.)
	call qmove(pxareg,pyareg)
	call qplot(pxereg,pyareg)
	call qplot(pxereg,pyereg)
	call qplot(pxareg,pyereg)
	call qplot(pxareg,pyareg)

	if( mode .eq. 0 ) return

	stop 'error stop frame: only mode == 0 is allowed'
	end

c*************************************************************

	subroutine basmima(mode)

c computes min/max of basin
c
c mode	0: exact dimensions  1: larger dimensions

	use basin

	implicit none

	integer mode

	include 'param.h'

	real xmin,ymin,xmax,ymax
	common /bamima/ xmin,ymin,xmax,ymax

	integer k
	real x,y
	real dist

	real rdist,rround

	call basinit

	xmin = xgv(1)
	xmax = xgv(1)
	ymin = ygv(1)
	ymax = ygv(1)
	do k=2,nkn
	  x = xgv(k)
	  y = ygv(k)
	  if( x .gt. xmax ) xmax = x
	  if( x .lt. xmin ) xmin = x
	  if( y .gt. ymax ) ymax = y
	  if( y .lt. ymin ) ymin = y
	end do

	if( mode .eq. 0 ) return

	dist = rdist(xmin,ymin,xmax,ymax)

	xmin = rround(xmin,dist,-1)
	ymin = rround(ymin,dist,-1)
	xmax = rround(xmax,dist,+1)
	ymax = rround(ymax,dist,+1)

	end

c*************************************************************

	subroutine setbas(x0,y0,x1,y1)

c sets min/max of basin to plot

	implicit none

	real x0,y0,x1,y1

	real xmin,ymin,xmax,ymax
	common /bamima/ xmin,ymin,xmax,ymax

	call basinit

	xmin = x0
	ymin = y0
	xmax = x1
	ymax = y1

	end

c*************************************************************

	subroutine getbas(x0,y0,x1,y1)

c gets min/max of basin to plot

	implicit none

	real x0,y0,x1,y1

	real xmin,ymin,xmax,ymax
	common /bamima/ xmin,ymin,xmax,ymax

	x0 = xmin
	y0 = ymin
	x1 = xmax
	y1 = ymax

	end

c*************************************************************

	subroutine boundline

c plots boundary line for lagoon

	implicit none

	character*80 bndlin

	integer n
	real x,y

	call basinit

	call getfnm('bndlin',bndlin)
	if( bndlin .eq. " " ) return

	call qcomm('plotting boundary line')
	call qgray(0.)

	open(1,file=bndlin,status='old',err=99)
    1	continue
	read(1,*,end=2) x,y,n
	if( n .eq. 1 ) then
		call qmove(x,y)
	else
		call qplot(x,y)
	end if
	goto 1
    2	continue
	close(1)

	return
   99	continue
	write(6,*) 'error opening boundary file:'
	write(6,'(a60)') bndlin
	stop 'error stop boundline'
	end

c*************************************************************

	subroutine reggrid(ngrid,dist,gray)

c plots regular grid
c
c if ngrid > 0                          use ngrid
c if ngrid <= 0 and dist > 0.           use dist
c otherwise                             do nothing

	implicit none

        integer ngrid           !divide into ngrid parts
	real dist               !regular distance of grid
	real gray		!gray value

	real xdmin,ydmin,xdmax,ydmax
	real x,y,x0,y0,dx,dy
	integer n,nx,ny

	real xmin,ymin,xmax,ymax
	common /bamima/ xmin,ymin,xmax,ymax

	integer iround
	real rround

	call basinit

	xdmin = rround(xmin,dist,-1)
	xdmax = rround(xmax,dist,+1)
	ydmin = rround(ymin,dist,-1)
	ydmax = rround(ymax,dist,+1)

        if( ngrid .gt. 0 ) then
          nx = ngrid - 1
          ny = ngrid - 1
          x0 = xmin
          y0 = ymin
          dx = (xmax-xmin)/ngrid
          dy = (ymax-ymin)/ngrid
        else if( dist .gt. 0. ) then
	  nx = iround((xdmax-xdmin)/dist)
	  ny = iround((ydmax-ydmin)/dist)
          x0 = xdmin
          y0 = ydmin
          dx = dist
          dy = dist
        else
          return
        end if

	call qcomm('regular grid')
	call qgray(gray)

	do n=0,nx
	  x = x0 + n * dx
	  if( x .gt. xmin .and. x .lt. xmax ) then
	    call qline(x,ymin,x,ymax)
	  end if
	end do

	do n=0,ny
	  y = y0 + n * dy
	  if( y .gt. ymin .and. y .lt. ymax ) then
	    call qline(xmin,y,xmax,y)
	  end if
	end do

	call qgray(0.)

	end

c*************************************************************

	subroutine bastlscale

c computes typical length scale for basin
c
c dxygrd is used first
c then typls is used if given
c else it is computed from grid

	use basin, only : nkn,nel,ngr,mbw

	implicit none


	integer ie
	real area,ao
	real dist,typls			!typical length scale
	real fact,afact
	real dxygrd
	double precision acu

	real getpar,aomega_elem

	acu = 0.
	do ie=1,nel
	  ao = aomega_elem(ie)
	  acu = acu + ao
	end do
	area = 24. * acu / nel		!a little bit bigger than one element

	call spherical_fact(fact,afact)	!correct for spherical coordinates
	dist = sqrt(area/afact)

	typls = getpar('typls')
	dxygrd = getpar('dxygrd')
	if( dxygrd .gt. 0. ) then
	  typls = dxygrd
	else if( typls .gt. 0. ) then
	  typls = typls
	else
	  typls = dist
	end if
	call putpar('typls',typls)

	write(6,*) 'typical length scale for basin: ',dist
	write(6,*) 'typical length scale used     : ',typls

	end

c***************************************************
c***************************************************
c***************************************************

	function roundm(r,mode)

c rounds r to the closest given value
c
c r		value to round
c mode
c		 0: do not round, return r
c		 1: round to higher value
c		-1: round to lower value

	implicit none

	real roundm
	real r
	integer mode

	double precision raux,fact,sign,rr
	integer i

	integer ndim
	parameter (ndim = 4)
	real rmaster(ndim),eps
	save rmaster,eps
	data eps /1.e-5/
	data rmaster /1.,2.,5.,10./
c	data rmaster /1.,2.,4.,5.,8.,10./

	roundm = r
	raux = 0
	if( mode .eq. 0 ) return

	fact = 1.0d+0
	if( r .lt. 0. ) then
	  sign = -1.
	  rr = -r
	else
	  sign = 1.
	  rr = r
	end if

	do while( rr*fact .gt. 10. )
	  fact = 0.1d+0 * fact
	end do
	do while( rr*fact .lt. 1. )
	  fact = 10.d+0 * fact
	end do

	rr = rr * fact		!rr is between 1. and 10.

	if( mode .gt. 0 ) then
	  do i=ndim,1,-1
	    if( rr .le. rmaster(i) + eps ) raux = rmaster(i)
	  end do
	else if( mode .lt. 0 ) then
	  do i=1,ndim
	    if( rr .ge. rmaster(i) - eps ) raux = rmaster(i)
	  end do
	end if

	roundm = sign * raux / fact

	end

c****************************************************

	function rround(r,rmaster,mode)

c rounds r to next rmaster value
c
c r		value to round
c rmaster	value to which r is rounded (must be positive)
c mode
c		 0: do not round, return r
c		 1: round to higher value
c		-1: round to lower value
c
c negative values are respected

	implicit none

	real rround
	real r,rmaster
	integer mode

	integer iaux
	real raux

	if( rmaster == 0. ) then	!FIXME
	  rround = r
	  return
	end if

	iaux = r/rmaster
	raux = iaux * rmaster

	if( mode .gt. 0 ) then
	  if( raux .lt. r ) raux = raux + rmaster
	else if( mode .lt. 0 ) then
	  if( raux .gt. r ) raux = raux - rmaster
	else
	  raux = r
	end if

	rround = raux

	return
	end

c****************************************************

        function rdist(xmin,ymin,xmax,ymax)

c computes gridspacing for frame (4-7 grid lines)
c
c xmin,ymin     coordinates of lower left point
c xmax,ymax     coordinates of upper rigth point

	implicit none

	real rdist
	real xmin,ymin,xmax,ymax

	real xdist,ydist,dist,fdist
	integer istell,lines

	integer iround

        xdist=xmax-xmin
        ydist=ymax-ymin

        if( xdist .gt. ydist ) then
		dist = xdist
	else
		dist = ydist
	end if

        dist=iround(dist)

        istell=log10(dist)
        fdist=10**istell
        lines=dist/fdist

        if(lines.le.3) fdist=fdist*0.5
        if(lines.ge.8) fdist=fdist*2.

        rdist=fdist

        return
	end

c**************************************************************

        function divdist(x,n,mode)

c divides x into n equal pieces (with rounding)

        implicit none

	real divdist	!length part (rounded)
        real x          !length to be divided
        integer n       !number of pieces to create
        integer mode    !0: get closest to n -1:get less  +1: get more

	logical debug
        integer lines,lines_high,lines_low
        real dist_high,dist_low
        real dist,fdist

        real roundm

	debug = .false.

        dist=nint(x)

        dist_high = roundm(dist/n,+1)
        dist_low  = roundm(dist/n,-1)

        lines_high = dist / dist_low
        lines_low  = dist / dist_high

        if( mode .lt. 0 ) then
          fdist = dist_high
        else if( mode .gt. 0 ) then
          fdist = dist_low
        else
          if( lines_high-n .lt. n-lines_low ) then
            fdist = dist_low
          else
            fdist = dist_high
          end if
        end if

	if( debug ) then
        write(6,*) '-------------------'
        write(6,*) x,dist,n,mode
        write(6,*) lines_low,lines_high
        write(6,*) dist_high,dist_low
        write(6,*) fdist,int(dist/fdist)
        write(6,*) '-------------------'
	end if

        divdist = fdist

        end

c**************************************************************

	subroutine frac_pos(r,np)

c computes number of fractional digits of real r

	implicit none

	real r
	integer np

	real eps
	integer ieps,ir

	eps = 1.e-5
	ieps = 100000

	np = 0
	ir = nint(ieps*abs(r))
	if( ir .eq. 0 ) return
	np = 5

	do while( 10*(ir/10) .eq. ir )
	  ir = ir/10
	  np = np - 1
	  !write(6,*) 'new pos: ',ir,np
	end do

	end

c**************************************************************

	subroutine frac_pos1(r,np)

c computes number of fractional digits of real r

	implicit none

	real r
	integer np

	real rr,ri
	real eps
	integer ieps

	eps = 1.e-5
	ieps = nint(1./eps)

	rr = r

	rr = abs(r)
	rr = eps*nint(rr/eps)
	ri = float(int(rr))
	np = 0
	write(6,*) 'in frac_pos: ',r,rr,ri,np

	do while( abs(rr-ri) .gt. eps )
	  !write(6,*) np,r,rr
	  np = np + 1
	  if( np .gt. 5 ) goto 99
	  rr = rr * 10.
	  ri = float(int(rr))
	  write(6,*) 'in frac_pos: ',r,rr,ri,np
	end do

	return
   99	continue
	write(6,*) np,r,rr,ri,rr-ri
	stop 'error stop frac_pos: internal error'
	end

c**************************************************************
c**************************************************************
c**************************************************************

	subroutine spherical_fact(fact,afact)

c computes factors for for spherical coordinates

	implicit none

	real fact	!factor in x direction
	real afact	!area factor between geo and cartesian coordinates

	real geomile,onedeg
	parameter( geomile = 1855.325 , onedeg = 60.*geomile )

	real y,pi,rad
	real x0,y0,x1,y1

	real getpar
	logical is_spherical

	fact = 1.
	afact = 1.

        if( .not. is_spherical() ) return	!only for spherical

	call getbas(x0,y0,x1,y1)
	y = 0.5*(y0+y1)
	pi = 4.*atan(1.)
	rad = pi/180.

	fact = cos(y*rad)
	afact = fact * onedeg * onedeg

	end

c**************************************************************

	subroutine handle_spherical

c handles spherical coordinates

	implicit none

	real fact,afact
	logical is_spherical

	call spherical_fact(fact,afact)

	call qfact(fact,1.0)

        if( is_spherical() ) then
	  write(6,*) 'Using factor for spherical coordinates: ',fact
	else
	  write(6,*) 'Using factor for coordinates: ',fact
	end if

	end

c**************************************************************
c**************************************************************
c**************************************************************

	subroutine adjust_reg_grid_spacing(dreg,imicro)

c checks if regular grid should be written

	implicit none

	real dreg
	integer imicro

	logical is_spherical,is_box_given

	if( dreg .ge. 0. ) return		!already given
        if( .not. is_spherical() ) then		!only for spherical
	  if( dreg .lt. 0. ) dreg = 0.
	  return
	end if
	!if( .not. is_box_given('leg') ) return	!no legend was requested

	call compute_reg_grid_spacing(dreg,imicro)

	end

c**************************************************************

	subroutine compute_reg_grid_spacing(dreg,imicro)

c tries to find best regular grid spacing value

	implicit none

	real dreg
	integer imicro

	real x0,y0,x1,y1
	real dx,dy,dxy
	real dsreg

	real rnext,rnextsub

	call getbas(x0,y0,x1,y1)

	dx = x1 - x0
	dy = y1 - y0
	!dxy = max(dx,dy)
	dxy = min(dx,dy)	!use smaller side

	dxy = dxy/4.		!around 4 grid lines

	dreg = rnext(dxy,1)
	dsreg = rnextsub(dreg)
	imicro = nint(dreg/dsreg)

	write(6,*) 'reg,micro: ',dreg,dsreg,imicro

	write(6,*) 'new reggrd = ',dreg,x0,x1,y0,y1

	end

c**************************************************************

	subroutine plot_reg_grid

c handles plotting of regular grid

	implicit none

	integer ngrid
	real reggrd
	real reggry
	real dreg

	real getpar

	write(6,*) 'starting plot_reg_grid...'

	reggrd = getpar('reggrd')
	reggry = getpar('reggry')

	!call adjust_reg_grid_spacing(reggrd)	!check if plotted automatically
	! we do not need the above call, because we automatically
	! only label the plot, but do not plot lines
	! if you want lines you have to explicitly set reggrd

	if( reggrd .le. 0. ) goto 1
	if( reggry .ge. 1. ) goto 1	!no white painting

	ngrid = 0
	call reggrid(ngrid,reggrd,reggry)

    1	continue
	write(6,*) 'ending plot_reg_grid...'

	end

c**************************************************************

	subroutine label_reg_grid

c handles labeling of regular grid

	implicit none

	real xmin,ymin,xmax,ymax
	common /bamima/ xmin,ymin,xmax,ymax

	integer nx,ny,n,i,nc
	integer imicro
	integer nxymax
	real reggrd
	real xvmin,yvmin,xvmax,yvmax
	real xdmin,ydmin,xdmax,ydmax
	real x0,y0,dx,dy
	real dist,x,y
	real size,ftext
	character*10 string
	logical bdebug

	real getpar
	real rround
	integer ialfa

	bdebug = .false.
	bdebug = .true.

	call basinit

	write(6,*) 'starting label_reg_grid...'

	size = 0.5	!leave this space around plot for labeling
	ftext = 2.8	!factor to shift text vertically
	nxymax = 50	!not more than these number of reg grids

	reggrd = getpar('reggrd')
	imicro = nint(getpar('regdst'))

	call adjust_reg_grid_spacing(reggrd,imicro)	!check if automatic

	if( bdebug ) write(6,*) 'reggrd: ',reggrd

	if( reggrd .eq. 0. ) return

	call frame(0)

	call qgetvp(xvmin,yvmin,xvmax,yvmax)

	xvmin = xvmin + size
	yvmin = yvmin + size
	xvmax = xvmax - size
	yvmax = yvmax - size

	call qsetvpnc(xvmin,yvmin,xvmax,yvmax)

c here labeling

	dist = reggrd
	call frac_pos(dist,nc)
	if( nc .eq. 0 ) nc = -1
	write(6,*) 'label_reg_grid: ',dist,nc,imicro

	xdmin = rround(xmin,dist,-1)
	xdmax = rround(xmax,dist,+1)
	ydmin = rround(ymin,dist,-1)
	ydmax = rround(ymax,dist,+1)

	nx = nint((xdmax-xdmin)/dist)
	ny = nint((ydmax-ydmin)/dist)
        x0 = xdmin
        y0 = ydmin
        dx = dist
        dy = dist

	if( nx .gt. nxymax .or. ny .gt. nxymax ) goto 99

	call qcomm('labeling regular grid')
	call qfont('Times-Roman')
	call qgray(0.0)

	call qtxts(9)

	do n=0,nx
	 x = x0 + n * dx
	 if( x .gt. xmin .and. x .lt. xmax ) then
	  i = ialfa(x,string,nc,-1)
	  call qtxtcr(0.,-ftext)
	  call qtext(x,ymax,string(1:i))
	  call qtxtcr(0.,+ftext)
	  call qtext(x,ymin,string(1:i))
	 end if
	end do

	do n=0,ny
	 y = y0 + n * dy
	 if( y .gt. ymin .and. y .lt. ymax ) then
	  i = ialfa(y,string,nc,-1)
	  call qtxtr(90.)
	  call qtxtcr(0.,-ftext)
	  call qtext(xmin,y,string(1:i))
	  call qtxtr(-90.)
	  call qtxtcr(0.,-ftext)
	  call qtext(xmax,y,string(1:i))
	 end if
	end do

	call qtxtcr(-1.,-1.)
	call qtxtr(0.)
	call qgray(0.)

	call bw_frame(imicro,x0,y0,dx,dy,xmin,xmax,ymin,ymax)

	call qsetvp(xvmin,yvmin,xvmax,yvmax)

	write(6,*) 'ending label_reg_grid...'

	return
   99	continue
	write(6,*) 'nx,ny: ',nx,ny
	stop 'error stop label_reg_grid: nx,ny too high'
	end

c**************************************************************

	subroutine bw_frame(imicro,x0,y0,dx,dy,xmin,xmax,ymin,ymax)

c plot black/white frame around geographical grid

	implicit none

	integer imicro
	real x0,y0,dx,dy,xmin,xmax,ymin,ymax

	integer i
	real xcm,ycm
	real hx,hy
	real ddx,ddy
	real x,y,xx0,yy0
	real xn,yn
	real fact

c------------------------------------------------------------
c set parameters
c------------------------------------------------------------

	if( dx .le. 0. .or. dy .le. 0. ) goto 97

	call qcm(xcm,ycm)
	hx = 0.1*xcm
	hy = 0.1*ycm

	fact = 1.5

c	here we should define how many boxes inbetween labels -> empirical

	if( imicro .le. 0 ) then
	  ddx = dx/2.
	  if( ddx .gt. fact*xcm ) ddx = dx/4.
	  ddy = dy/2.
	  if( ddy .gt. fact*ycm ) ddy = dy/4.
	else
	  ddx = dx/imicro
	  ddy = dy/imicro
	end if

c------------------------------------------------------------
c plot along x-axis
c------------------------------------------------------------

	x = x0
	do while( x .lt. xmin )
	  x = x + dx
	end do
	xx0 = x

	x = x0
	if( imicro .ge. 0 ) x = 2.*xmax
	do while( x .le. xmax )		!plot ticks on x axis
	  call qline(x,ymax,x,ymax+hx)
	  call qline(x,ymin,x,ymin-hx)
	  x = x + dx
	end do

	x = xx0
	if( imicro .lt. 0 ) x = xmax
	i = 0
	do while( x .lt. xmax )
	  i = mod(i+1,2)
	  xn = min(x+ddx,xmax)
	  call make_box(i,x,ymax,xn,ymax+hy)
	  call make_box(i,x,ymin,xn,ymin-hy)
	  x = x + ddx
	end do

	x = xx0
	if( imicro .lt. 0 ) x = xmin
	i = 1
	do while( x .gt. xmin )
	  i = mod(i+1,2)
	  xn = max(x-ddx,xmin)
	  call make_box(i,xn,ymax,x,ymax+hy)
	  call make_box(i,xn,ymin,x,ymin-hy)
	  x = x - ddx
	end do

c------------------------------------------------------------
c plot small corner boxes
c------------------------------------------------------------

	if( imicro .ge. 0 ) then
	  call make_box(0,xmin-hx,ymax,xmin,ymax+hy)
	  call make_box(0,xmax,ymax,xmax+hx,ymax+hy)
	  call make_box(0,xmin-hx,ymin-hy,xmin,ymin)
	  call make_box(0,xmax,ymin-hy,xmax+hx,ymin)
	end if

c------------------------------------------------------------
c plot along y-axis
c------------------------------------------------------------

	y = y0
	do while( y .lt. ymin )
	  y = y + dy
	end do
	yy0 = y

	y = y0
	if( imicro .ge. 0 ) y = 2.*ymax
	do while( y .le. ymax )		!plot ticks on y axis
	  call qline(xmin,y,xmin-hx,y)
	  call qline(xmax,y,xmax+hx,y)
	  y = y + dy
	end do

	y = yy0
	if( imicro .lt. 0 ) y = ymax
	i = 0
	do while( y .lt. ymax )
	  i = mod(i+1,2)
	  yn = min(y+ddy,ymax)
	  call make_box(i,xmin-hx,y,xmin,yn)
	  call make_box(i,xmax,y,xmax+hx,yn)
	  y = y + ddy
	end do

	y = yy0
	if( imicro .lt. 0 ) y = ymin
	i = 1
	do while( y .gt. ymin )
	  i = mod(i+1,2)
	  yn = max(y-ddy,ymin)
	  call make_box(i,xmin-hx,yn,xmin,y)
	  call make_box(i,xmax,yn,xmax+hx,y)
	  y = y - ddy
	end do

c------------------------------------------------------------
c end of routine
c------------------------------------------------------------

	return
   97	continue
	write(6,*) 'dx,dy: ',dx,dy
	stop 'error stop bw_frame: dx,dy'
	end

c**************************************************************

	subroutine make_box(ic,x1,y1,x2,y2)

	implicit none

	integer ic
	real x1,y1,x2,y2

	if( ic .eq. 1 ) then
	  call qrfill(x1,y1,x2,y2)
	end if

	call qline(x1,y1,x2,y1)
	call qline(x2,y1,x2,y2)
	call qline(x2,y2,x1,y2)
	call qline(x1,y2,x1,y1)

	end

c**************************************************************
c**************************************************************
c**************************************************************

	subroutine plot_regular_points

c plots regular points

	implicit none

	include 'supout.h'

        real xmin,ymin,xmax,ymax
        common /bamima/ xmin,ymin,xmax,ymax

	integer nx,ny
	real ddx,ddy,dxy
	real x0,y0,dx,dy,flag
	real fact,afact
	real x,y

        nx = nint(regp(1))
        ny = nint(regp(2))
        x0 = regp(3)
        y0 = regp(4)
        dx = regp(5)
        dy = regp(6)
        flag = regp(7)

	if( dx <= 0. .or. dy <= 0. ) return
	if( nx <= 0. .or. ny <= 0. ) return

	write(6,*) 'plotting regular points:'
	write(6,*) x0,y0,dx,dy,flag
	write(6,*) xmin,ymin,xmax,ymax

	do while( x0 + dx < xmin )
	  x0 = x0 + dx
	end do
	do while( y0 + dy < ymin )
	  y0 = y0 + dy
	end do

	nx = (xmax-xmin)/dx
	ny = (ymax-ymin)/dy
	nx = max(nx,20)
	ny = max(ny,20)
	ddx = 0.3*(xmax-xmin)/(nx-1)
	ddy = 0.3*(ymax-ymin)/(ny-1)
	dxy = min(ddx,ddy)
	call spherical_fact(fact,afact)
	ddx = dxy/fact
	ddy = dxy

	write(6,*) nx,ny,ddx,ddy

	call qgray(0.)

	y = y0
	do while( y <= ymax )
	  x = x0
	  do while( x <= xmax )
	    call plot_plus(x,y,ddx,ddy)
	    x = x + dx
	  end do
	  y = y + dy
	end do

	end

c**************************************************************

	subroutine plot_meteo_points(mode,file)

c plots special points from meteo file
c
c name of coords.dat has to passed into the subroutine
c if sea_land.dat exists it is used, otherwise we can do without

	implicit none

	integer mode
	character*(*) file

	integer ndim
	parameter (ndim=30000)

	integer nx,ny,nz,n
	integer idum,i,ios
	real xx,yy
	real fact,afact
	real xmin,xmax,ymin,ymax
	real dx,dy,dxy
	real x(ndim),y(ndim),rf(ndim)
	save n,x,y,rf,dx,dy,dxy

	integer icall
	save icall
	data icall / 0 /

	if( icall .eq. -1 ) return

	if( icall .eq. 0 ) then
	  open(1,file=file,status='old',form='formatted',err=88)

	  if( mode .eq. 1 ) then
	    read(1,*) nx,ny,nz
	    n = nx*ny
	    if( n .gt. ndim ) goto 99
	    do i=1,n
	      read(1,*) idum,idum,x(i),y(i)
	    end do
	  else if( mode .eq. 2 ) then
	    i = 0
	    do
	      !read(1,*,iostat=ios) xx,yy
	      read(1,*,iostat=ios) idum,idum,xx,yy
	      if( ios < 0 ) exit
	      if( ios > 0 ) goto 98
	      i = i + 1
	      if( i .gt. ndim ) goto 99
	      x(i) = xx
	      y(i) = yy
	    end do
	    n = i
	    nx = 10
	    ny = 10
	  else
	    write(6,*) 'mode = ',mode
	    stop 'error stop plot_meteo_points: mode not recognized'
	  end if

	  close(1)

	  call mima(x,n,xmin,xmax)
	  call mima(y,n,ymin,ymax)

	  dx = 0.3*(xmax-xmin)/(nx-1)
	  dy = 0.3*(ymax-ymin)/(ny-1)
	  dxy = min(dx,dy)
	  call spherical_fact(fact,afact)
	  dx = dxy/fact
	  dy = dxy

	  do i=1,nx*ny
	    rf(i) = 1.
	  end do
	  open(1,file='sea_land.dat',status='old',form='formatted',err=77)
	  read(1,*)
	  read(1,*) nx,ny,nz
	  if( n .ne. nx*ny ) goto 97
	  read(1,*) (rf(i),i=1,n)
	  close(1)
   77	  continue
	  icall = 1
	end if

	write(6,*) 'plotting meteo points'
	call qcomm('plotting meteo points')
	call qgray(0.)
	!dx = 0.02
	!dy = dx
	do i=1,n
	  if( rf(i) .gt. 0 ) then
	    call qgray(0.)
	    call plot_plus(x(i),y(i),dx,dy)
	  else
	    call qgray(0.5)
	    call plot_cross(x(i),y(i),dx,dy)
	  end if
	end do

	return
   88	continue
	icall = -1
	write(6,*) 'message from plot_meteo_points:'
	write(6,*) 'no coords.dat file ... cannot plot coordinates'
	write(6,*) 'file: ',file
	stop 'error stop plot_meteo_points: file'
   98	continue
	write(6,*) 'file: ',trim(file)
	write(6,*) 'mode,irec: ',mode,i
	stop 'error stop plot_meteo_points: read error'
   99	continue
	write(6,*) n,ndim
	stop 'error stop plot_meteo_points: ndim'
   97	continue
	write(6,*) nx,ny,n,ndim
	stop 'error stop plot_meteo_points: different length'
	end

c**************************************************************

	subroutine plot_plus(x,y,dx,dy)

	implicit none

	real x,y,dx,dy

	call qline(x,y-dy,x,y+dy)
	call qline(x-dx,y,x+dx,y)

	end

	subroutine plot_cross(x,y,dx,dy)

	implicit none

	real x,y,dx,dy

	call qline(x-dx,y-dy,x+dx,y+dy)
	call qline(x-dx,y+dx,x+dx,y-dx)
	
	end

c**************************************************************

	subroutine plot_islands(cgray)

c plots islands gray

	use mod_geom
	use basin

	implicit none

	real cgray	!color
	

	include 'param.h'

	real xmin,ymin,xmax,ymax
	common /bamima/ xmin,ymin,xmax,ymax

	logical bouter
	integer k,kstart,kk,kn,ko
	integer n,nis
	real x0,y0,area,cg
	real xa(nkn),ya(nkn)

	real areapoly

	if( cgray < 0 ) return

	cg = cgray
	bouter = .false.		! plot external island
	if( cg >= 2 ) then
	  cg = cg - 2
	  bouter = .true.
	end if

	call qgray(cg)

	x0 = 0.5*(xmax-xmin)
	y0 = 0.5*(ymax-ymin)

	nis = 0

	do k=1,nkn
	  if( kantv(1,k) > 0 ) then
	    kstart = k
	    kk = k
	    n = 0
	    nis = nis + 1
	    do
	      kn = kantv(1,kk)
	      kantv(1,kk) = -kantv(1,kk)
	      ko = kk
	      kk = kn
	      if( kn == 0 ) stop 'error stop plot_islands: internal error'
	      if( kn < 0 ) exit
	      n = n + 1
	      xa(n) = xgv(kk)
	      ya(n) = ygv(kk)
	    end do
	    area = areapoly(n,xa,ya)
	    if( kstart /= ko ) then
		stop 'error stop plot_islands: internal error'
	    end if
	    !write(69,*) 'debug: ',kk,kn,ko,kstart
	    !write(69,*) 'island found: ',nis,n,area
	    write(6,*) 'island found: ',nis,n,area
	    if( area < 0 ) then	!real island
	      call qafill(n,xa,ya)
	    else if( bouter ) then
	      call plot_outer(n,xa,ya)
	    end if
	  end if
	end do

	do k=1,nkn
	  if( kantv(1,k) < 0 ) kantv(1,k) = -kantv(1,k)
	  if( kantv(2,k) < 0 ) kantv(2,k) = -kantv(2,k)
	end do

	call qgray(0.)

	end

c**************************************************************

	subroutine plot_outer(n,xa,ya)

c plots outer island

	implicit none

	integer n
	real xa(n),ya(n)

	real xmin,ymin,xmax,ymax
	common /bamima/ xmin,ymin,xmax,ymax

	integer i,ilow
	real dx,dy,xlow,ylow
	real xn(n+7),yn(n+7)

	dx = xmax - xmin
	dy = ymax - ymin

	ilow = 0
	ylow = ya(1) + 1.
	do i=1,n
	  if( ya(i) < ylow ) then
	    ylow = ya(i)
	    ilow = i
	  end if
	end do
	xlow = xa(ilow)

	do i=1,ilow
	  xn(i) = xa(i)
	  yn(i) = ya(i)
	end do
	
	i = ilow

	i = i + 1
	xn(i) = xlow
	yn(i) = ymin - dy

	i = i + 1
	xn(i) = xmin - dx
	yn(i) = ymin - dy

	i = i + 1
	xn(i) = xmin - dx
	yn(i) = ymax + dy

	i = i + 1
	xn(i) = xmax + dx
	yn(i) = ymax + dy

	i = i + 1
	xn(i) = xmax + dx
	yn(i) = ymin - dy

	i = i + 1
	xn(i) = xlow
	yn(i) = ymin - dy

	do i=ilow,n
	  xn(i+7) = xa(i)
	  yn(i+7) = ya(i)
	end do

	call qafill(n+6,xn,yn)

	end

c**************************************************************

