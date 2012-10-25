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
c subroutine sbline( bndfile )		stores name of boundary file
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

	real xmin,ymin,xmax,ymax
	common /bamima/ xmin,ymin,xmax,ymax
	character*80 bndfil
	common /bndfil/bndfil
	save /bamima/
	save /bndfil/

	logical binit
	save binit
	data binit /.false./

	if( binit ) return

	binit = .true.

	call basmima(0)		!compute exact dimensions as default
	call bastlscale		!computes typical length scale
	bndfil = " "

	end

c**************************************************************

	subroutine bash(mode)

c plots basin (shell)
c
c mode	0: only scaling  1: net  2: boundary  3: net in gray
c
c bash MUST be called first with mode == 0, and then
c with the desired mode

	implicit none

	integer mode

	real x0,y0,x1,y1
	real dxygrd,x,y
	character*80 bndlin
	real getpar
	logical inboxdim
        logical is_spherical
	logical is_box_given

	call basinit

c if mode == 0 -> do scaling

	if( mode .eq. 0 ) then

	  call basmima(0)
	  if( inboxdim(' ',x0,y0,x1,y1) ) then	!size of plotting is given
	    call setbas(x0,y0,x1,y1)
	  end if

	  call getbas(x0,y0,x1,y1)		!bug fix 2.3.2005

c	  prepare regular grid

	  dxygrd = getpar('dxygrd')
	  x = x0 + dxygrd/2.
	  y = y0 + dxygrd/2.
	  call setgeo(x,y,dxygrd,dxygrd,-999.)

	  call annote		!annotation
	  call basin(0)		!scaling
	  call label_reg_grid

	  return

	end if

c frame

	call frame(0)
	call plot_reg_grid

c plot

        call qcomm('Plotting basin')
        call basin(mode)

c boundary line

	call getfnm('bndlin',bndlin)
	if( bndlin .ne. " " ) call sbline(bndlin)
	call boundline

c user defined legend

	call legplo

c legend (north and scale)

	!if( inboxdim('leg',x0,y0,x1,y1) ) then	!write legend
	if( is_box_given('leg') ) then	!write legend
          if( is_spherical() ) then
            write(6,*) 'coordinates are spherical'
            write(6,*) 'no north and scale written...'
            return
	  else
	    call legend(x0,y0,x1,y1)
          end if
	end if

c end of routine

	end

c*************************************************************

	subroutine basin(mode)

c plots basin
c
c mode	0: only scaling  1: net  2: boundary  3: net in gray

	implicit none

	integer mode

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nen3v(3,1)
	integer ipev(1), ipv(1)
	real xgv(1), ygv(1)
	integer kantv(2,1)
	common /xgv/xgv, /ygv/ygv
	common /nen3v/nen3v
	common /ipev/ipev, /ipv/ipv
	common /kantv/kantv

	real xmin,ymin,xmax,ymax
	common /bamima/ xmin,ymin,xmax,ymax

	integer ie,kn,ii,k
	real gray

	real getpar

	call basinit

	if( mode .eq. 0 ) then
	  call qworld(xmin,ymin,xmax,ymax)
	  call qrcfy
	  call handle_spherical
	  call bpixel
	end if

	if( mode .eq. 0 ) return

	if( mode .eq. 1 .or. mode .eq. 3 ) then		!net

	call qgray(0.)
	if( mode .eq. 3 ) then
	  gray = getpar('bgray')
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

	call qgray(0.)

	end if						!boundary

	if( mode .eq. 2 .or. mode .eq. 3 ) then

	call qgray(0.)

	do k=1,nkn
	  if( kantv(1,k) .ne. 0 ) then
	    kn=kantv(1,k)
	    if( kn .gt. k ) call qline(xgv(k),ygv(k),xgv(kn),ygv(kn))
	    kn=kantv(2,k)
	    if( kn .gt. k ) call qline(xgv(k),ygv(k),xgv(kn),ygv(kn))
	  end if
	end do

	end if

c	call boundline

	end

c*************************************************************

	subroutine basin_number(mode)

c plots basin with node and element numbers
c
c mode	1: node number   2: element number
c	positive: external    negative: internal

	implicit none

	integer mode

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nen3v(3,1)
	integer ipev(1), ipv(1)
	real xgv(1), ygv(1)
	integer kantv(2,1)
	common /xgv/xgv, /ygv/ygv
	common /nen3v/nen3v
	common /ipev/ipev, /ipv/ipv
	common /kantv/kantv

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

	character*80 line
	real fact
	real xvmin,yvmin,xvmax,yvmax
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
c 0: just box around plot    1: geographical coordinates

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

c========================
	plon0 = 16.
	plat0 = 41.
	dlon = 1.
	dlat = 1.
c========================

	call mercoo(1,alonmin,alatmin,pxareg,pyareg,1,plat0,plon0,100.)
	call mercoo(1,alonmax,alatmax,pxereg,pyereg,1,plat0,plon0,100.)

	iaux=alonmin
	alonmin=iaux
	iaux=alatmin
	alatmin=iaux

	call qgray(0.5)

	alon=alonmin
	alat=alatmin

	do while( alon .le. alonmax )
	  call mercoo(0,alon,alat,x,y,1,plat0,plon0,100.)
	  call qline(x,pyareg,x,pyereg)
	  alon = alon + dlon
	end do
	  
	do while( alat .le. alatmax )
	  call mercoo(0,alon,alat,x,y,1,plat0,plon0,100.)
	  call qline(pxareg,y,pxereg,y)
	  alat = alat + dlat
	end do

	call qgray(0.)

	return
	end

c*************************************************************

	subroutine basmima(mode)

c computes min/max of basin
c
c mode	0: exact dimensions  1: larger dimensions

	implicit none

	integer mode

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nen3v(3,1)
	real xgv(1), ygv(1)
	common /nen3v/nen3v
	common /xgv/xgv, /ygv/ygv

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

	subroutine sbline( bndfile )

c stores name of boundary file

	implicit none

	character*(*) bndfile

	character*80 bndfil
	common /bndfil/bndfil

	bndfil = bndfile

	end

c*************************************************************

	subroutine boundline

c plots boundary line for lagoon

	implicit none

	character*80 bndfil
	common /bndfil/bndfil

	integer n
	real x,y

	call basinit

	if( bndfil .eq. " " ) return

	call qcomm('plotting boundary line')
	call qgray(0.)

	open(1,file=bndfil,status='old',err=99)
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
	write(6,'(a60)') bndfil
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

	implicit none

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	integer ie
	real area,ao
	real dist,typls			!typical length scale
	real dxygrd
	double precision acu

	real getpar,aomega_elem

	acu = 0.
	do ie=1,nel
	  ao = aomega_elem(ie)
	  acu = acu + ao
	end do
	area = 24. * acu / nel

	dist = sqrt(area)

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

	subroutine handle_spherical

c handles spherical coordinates

	implicit none

	real fact,y,pi,rad
	real x0,y0,x1,y1

	real getpar
	logical is_spherical

        if( .not. is_spherical() ) return	!only for spherical

	call getbas(x0,y0,x1,y1)
	y = 0.5*(y0+y1)
	pi = 4.*atan(1.)
	rad = pi/180.

	fact = cos(y*rad)

	call qfact(fact,1.0)

	write(6,*) 'Using factor for spherical coordinates: ',y,fact

	end

c**************************************************************

	subroutine adjust_reg_grid_spacing(dreg)

c checks if regular grid should be written

	implicit none

	real dreg

	logical is_spherical,is_box_given

	if( dreg .ge. 0. ) return		!already given
        if( .not. is_spherical() ) then		!only for spherical
	  if( dreg .lt. 0. ) dreg = 0.
	  return
	end if
	if( .not. is_box_given('leg') ) return	!no legend was requested

	call compute_reg_grid_spacing(dreg)

	end

c**************************************************************

	subroutine compute_reg_grid_spacing(dreg)

c tries to find best regular grid spacing value

	implicit none

	real dreg

	real x0,y0,x1,y1
	real dx,dy,dxy

	real rnext

	call getbas(x0,y0,x1,y1)

	dx = x1 - x0
	dy = y1 - y0
	dxy = max(dx,dy)

	dxy = dxy/4.		!around 4 grid lines

	dreg = rnext(dxy,1)

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

	if( reggrd .le. 0. ) return
	if( reggry .ge. 1. ) return	!no white painting

	ngrid = 0
	call reggrid(ngrid,reggrd,reggry)

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
	real reggrd
	real xvmin,yvmin,xvmax,yvmax
	real xdmin,ydmin,xdmax,ydmax
	real x0,y0,dx,dy
	real dist,x,y
	real size,ftext
	character*10 string

	real getpar
	real rround
	integer ialfa

	call basinit

	write(6,*) 'starting label_reg_grid...'

	size = 0.5	!leave this space around plot for labeling
	ftext = 2.8	!factor to shift text vertically

	reggrd = getpar('reggrd')
	imicro = nint(getpar('regdst'))

	call adjust_reg_grid_spacing(reggrd)	!check if plotted automatically

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
	write(6,*) dist,nc

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

	end

c**************************************************************

	subroutine bw_frame(imicro,x0,y0,dx,dy,xmin,xmax,ymin,ymax)

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

