c
c $Id: subcoord.f,v 1.5 2010-02-26 17:35:06 georg Exp $
c
c coordinate handling routines
c
c revision log :
c
c 07.05.2009    ggu     new framework (to be finished), new basic projection
c
c********************************************************************

	subroutine init_coords(iproj,c_param)

c initialization of transformation routines
c
c please note: the coordinates (xtrans,ytrans) are subtracted from the 
c	computed cartesian grid and are in meters

	implicit none

	integer iproj			!1: gauss-boaga  2: UTM
        double precision c_param(9)     !parameters for projection

	integer zone
	double precision xtrans,ytrans
	double precision lon0,lat0,phi

c if GB, zone (fuse) must be 1 or 2
c if UTM, zone must be 1-60
c
c GB: Venice is in fuse 2, central meridian is 15
c UTM: Venice is in zone 33, central meridian is 15

	integer proj
	common /coords1/ proj
	save /coords1/

	if( iproj .eq. 1 ) then		!GB
          zone = nint( c_param(1) )     !fuse for gauss-boaga, zone for UTM
          xtrans = c_param(2)           !extra shift in x [m]
          ytrans = c_param(3)           !extra shift in y [m]
	  call gb_init(zone)
	  call gb_trans(xtrans,ytrans)
	else if( iproj .eq. 2 ) then	!UTM
	  stop 'error stop init_coords: UTM not yet implemented'
	else if( iproj .eq. 3 ) then	!basic projection
          lon0 = c_param(1)             !longitude of origin
          lat0 = c_param(2)             !latitude of origin
          phi  = c_param(3)             !central latitude
	  call bs_init(lon0,lat0,phi)
	else
	  write(6,*) 'iproj = ',iproj
	  stop 'error stop init_coords: unknown projection'
	end if

	proj = iproj

	end

c********************************************************************

	subroutine convert_coords(mode,n,xc,yc,xg,yg)

c conversion of coordinates

	implicit none

	integer mode		! +1: cart to geo  -1: geo to cart
	integer n		! number of points
	real xc(1),yc(1)	! cartesian coordinates 
	real xg(1),yg(1)	! geographical coordinates 

        external gb_c2g,gb_g2c
        external bs_c2g,bs_g2c

	integer proj
	common /coords1/ proj

	integer i
	double precision east,north,lon,lat

	if( mode .eq. 1 ) then
	  write(6,*) 'using projection = ',proj
	  if( proj .eq. 1 ) then
            call apply_coords(n,xc,yc,xg,yg,gb_c2g)
	  else if( proj .eq. 3 ) then
            call apply_coords(n,xc,yc,xg,yg,bs_c2g)
	  end if
	else if( mode .eq. -1 ) then
	  if( proj .eq. 1 ) then
            call apply_coords(n,xg,yg,xc,yc,gb_g2c)
	  else if( proj .eq. 3 ) then
            call apply_coords(n,xc,yc,xg,yg,bs_g2c)
	  end if
	else
	  stop 'error stop convert_coords: unknown mode'
	end if

	end

c********************************************************************

        subroutine apply_coords(n,xs,ys,xt,yt,sub)
        
        implicit none

        external sub
        integer n
        real xs(1), ys(1)
        real xt(1), yt(1)

        integer i
        double precision xsource,ysource,xtarget,ytarget

	do i=1,n
	  xsource = xs(i)
	  ysource = ys(i)
	  call sub(xsource,ysource,xtarget,ytarget)
	  xt(i) = xtarget
	  yt(i) = ytarget
        end do

        end

c********************************************************************

	subroutine genproj_all_test
	end

c********************************************************************

c         program genproj_main
c         call genproj_all_test
c         end

c********************************************************************

	subroutine bs_init(lon0,lat0,phi)

        double precision lon0,lat0,phi

	double precision lon_0,lat_0,phi_0,xfact,yfact
	common /bs_param1/ lon_0,lat_0,phi_0,xfact,yfact
	save /bs_param1/

        double precision pi, rad

        lon_0 = lon0
        lat_0 = lat0
        phi_0 = phi

        pi = 4. * atan(1.d+0)
        rad = pi / 180.

        yfact = 60 * 1852               !distance of 1 min = 1 nautical mile
        xfact = yfact * cos(phi*rad)

        end

c********************************************************************

	subroutine bs_g2c(lon,lat,x,y)

c transformation from (lon/lat) to cartesian (x,y)

	implicit none

	double precision lon,lat !geographical coordinates
	double precision x,y  !cartesian coordinates (return)

	double precision lon_0,lat_0,phi_0,xfact,yfact
	common /bs_param1/ lon_0,lat_0,phi_0,xfact,yfact
	save /bs_param1/

        x = (lon-lon_0) * xfact
        y = (lat-lat_0) * yfact

        end

c********************************************************************

	subroutine bs_c2g(x,y,lon,lat)

c transformation from to cartesian (x,y) to (lon/lat)

	implicit none

	double precision x,y  !cartesian coordinates
	double precision lon,lat !geographical coordinates (return)

	double precision lon_0,lat_0,phi_0,xfact,yfact
	common /bs_param1/ lon_0,lat_0,phi_0,xfact,yfact
	save /bs_param1/

        lon =  lon_0 + x / xfact
        lat =  lat_0 + y / yfact

        end

c********************************************************************

c
c gauss-boaga transformations
c
c contents :
c
c gauss-boaga core routines
c
c subroutine gb_init(fuse)
c		intializes gaus-boaga routines
c subroutine gb_trans(xtrans,ytrans)
c		intializes translation of gauss-boaga coordinates
c subroutine gb_g2c(lon,lat,east,north)
c		transformation from (lon/lat) to gauss-boaga (east,north)
c subroutine gb_c2g(east,north,lon,lat)
c		transformation from gauss-boaga (east,north) to (lon/lat)
c
c helper routines
c
c subroutine lag2gb(xlag,ylag,e,n)
c		transforms from lagoon coordinates to gauss-boaga
c subroutine gb2lag(e,n,xlag,ylag)
c		transforms from gauss-boaga to lagoon coordinates
c subroutine lag2geo(xlag,ylag,lon,lat)
c		transforms from lagoon coordinates to lon/lat
c subroutine geo2lag(lon,lat,xlag,ylag)
c		transforms from lon/lat to lagoon coordinates
c
c ...other test routines...
c
c usage :
c
c gb_init(fuse) must be called at beginning to set fuse
c gb_trans(xtrans,ytrans) may be called if different origin is desired
c	if not called no translation will be performed
c
c example: in order to pass from lat/lon to lagoon coordinates use
c
c	call gb_init(2)
c	call gb_trans(2280000.D+0,5000000.D+0)
c
c and then
c
c	call gb_g2c(lon,lat,east,north)
c
c notes :
c
c gauss-boaga transformations
c
c fuso I  = W = 32 UTM   meridiano centrale  9 deg
c fuso II = E = 33 UTM   meridiano centrale 15 deg
c
c scale 0.9996
c
c UTM  uses ED50 (European Data 1950)
c GB (IGM) uses Roma40
c
c both use international ellipsoid (Hayford)
c
c ellipsoid:  c = 6397376.633  es = e2 = 0.0067681702
c
c Monte Mario:
c
c Roma40: 41 55 25.510 12 27 08.400
c Ed50:  41 55 31.487 12 27 08.933
c
c------------------------------------------------------------
c
c for Venice lagoon:venlag62_new2.grd
c
c fuso = 2

c to get GB coordinates (x_gb,y_gb) from venlag62 (x_lag,y_lag) compute:
c
c x_gb = x_lag + x0
c y_gb = y_lag + y0
c
c with
c
c x0 = 2330000. - 50000. = 2280000;
c y0 = 5000000.;
c
c revision log :
c
c 16.02.2006    ggu     separated from other routines and commented
c 11.03.2009    ggu     cleaned, easier initialization
c
c**********************************************************************

	subroutine gb_init(fuse)

c intializes gaus-boaga routines
c
c to obtain E/N use E = x + x0    N = y
c with x0 = 1500 km (fuse=1) and x0 = 2520 km (fuse=2)

	implicit none

	integer fuse			!1 or 2


c---------------------------------------------------------------- 
c proj.h
c----------------------------------------------------------------

	double precision aearth,bearth,flat,rflat,e,es
	common /proj_param01/ aearth,bearth,flat,rflat,e,es
	save /proj_param01/

	double precision k0
	common /proj_param02/ k0
	save /proj_param02/

	double precision zero,one,two,four,half
	common /proj_param11/ zero,one,two,four,half
	save /proj_param11/

	double precision pi,half_pi,quarter_pi,rad,rrad
	common /proj_param12/ pi,half_pi,quarter_pi,rad,rrad
	save /proj_param12/

	double precision eps,tol
	common /proj_param13/ eps,tol
	save /proj_param13/

	logical debug
	common /proj_param21/ debug
	save /proj_param21/

c---------------------------------------------------------------- 

c	include 'proj.h'

	integer fuso			!1 or 2
	common /gb_param1/ fuso
	double precision lambda0,x0,xtrans0,ytrans0
	common /gb_param2/ lambda0,x0,xtrans0,ytrans0
	save /gb_param1/,/gb_param2/

	call proj_init

	fuso = fuse

	if( fuso .eq. 1 ) then
	  lambda0 = 9.
	  x0 = 1500000.
	else if( fuso .eq. 2 ) then
	  lambda0 = 15.
	  x0 = 2520000.
	else
	  stop 'error stop gb_init: fuse'
	end if

	xtrans0 = 0.		! change origin of gb coordinates
	ytrans0 = 0.

c use own constants for international ellipsoid (Hayford)

	aearth	= 6397376.633
	k0	= 0.9996
	es	= 0.0067681702
	e	= sqrt(es)

        if( .not. debug ) return
        write(6,*) 'gb_init:'
        write(6,*) fuso,lambda0,x0
        write(6,*) aearth,k0
        write(6,*) e,es

	end

c**********************************************************************

	subroutine gb_trans(xtrans,ytrans)

c intializes translation of gauss-boaga coordinates
c
c (xtrans,ytrans) are the coordinates of the GB system that
c will be the new origin of the translation
c -> so these coordinates will be subtracted from the GB coordinates

	implicit none

	double precision xtrans,ytrans

	double precision lambda0,x0,xtrans0,ytrans0
	common /gb_param2/ lambda0,x0,xtrans0,ytrans0

	xtrans0 = xtrans
	ytrans0 = ytrans

	end

c**********************************************************************

	 subroutine gb_g2c(lon,lat,east,north)

c transformation from (lon/lat) to gauss-boaga (east,north)

	 implicit none

	 double precision lon,lat !geographical coordinates
	 double precision east,north  !coordinates gauss-boaga (return)


c---------------------------------------------------------------- 
c proj.h
c----------------------------------------------------------------

	double precision aearth,bearth,flat,rflat,e,es
	common /proj_param01/ aearth,bearth,flat,rflat,e,es
	save /proj_param01/

	double precision k0
	common /proj_param02/ k0
	save /proj_param02/

	double precision zero,one,two,four,half
	common /proj_param11/ zero,one,two,four,half
	save /proj_param11/

	double precision pi,half_pi,quarter_pi,rad,rrad
	common /proj_param12/ pi,half_pi,quarter_pi,rad,rrad
	save /proj_param12/

	double precision eps,tol
	common /proj_param13/ eps,tol
	save /proj_param13/

	logical debug
	common /proj_param21/ debug
	save /proj_param21/

c---------------------------------------------------------------- 

c	 include 'proj.h'

	 double precision lambda0,x0,xtrans0,ytrans0
	 common /gb_param2/ lambda0,x0,xtrans0,ytrans0

	 double precision a1,a2,a4,a6
	 parameter ( a1 = 111092.08210 , a2 = 16100.59187 )
	 parameter ( a4 = 16.96942 , a6 = 0.02226 )
	
	 double precision lambdap,phi,nu,nu1,xi,aux
	 !double precision c,e2
	 double precision x,y
	 
	 !c = aearth
	 !e2 = es

	 lambdap = rad * (lon - lambda0)
	 phi = rad * lat
	
	 nu1 = sqrt( 1.+es*(cos(phi)**2) )
	 xi = atan( tan(phi) / cos(nu1*lambdap) )
	 nu = sqrt( 1.+es*(cos(xi)**2) )
	
	 aux = cos(xi) * tan(lambdap) / nu
	 x = aearth * log(aux+sqrt(aux*aux+1.))
	 y = a1*xi/rad - a2*sin(2.*xi) + a4*sin(4.*xi) - a6*sin(6.*xi)
	
	 east = x + x0
	 north = y
	
	 east = east - xtrans0
	 north = north - ytrans0

	 if( .not. debug ) return
	 write(6,*) 'gb_g2c:'
	 write(6,*) lambdap/rad,lon
	 write(6,*) nu1,xi/rad,nu
	 write(6,*) aux,x,y

	 end
	
c**********************************************************************
	
	 subroutine gb_c2g(east,north,lon,lat)
	
c transformation from gauss-boaga (east,north) to (lon/lat)
	
	 implicit none
	
	 double precision east,north  !coordinates gauss-boaga
	 double precision lon,lat !geographical coordinates (return)
	

c---------------------------------------------------------------- 
c proj.h
c----------------------------------------------------------------

	double precision aearth,bearth,flat,rflat,e,es
	common /proj_param01/ aearth,bearth,flat,rflat,e,es
	save /proj_param01/

	double precision k0
	common /proj_param02/ k0
	save /proj_param02/

	double precision zero,one,two,four,half
	common /proj_param11/ zero,one,two,four,half
	save /proj_param11/

	double precision pi,half_pi,quarter_pi,rad,rrad
	common /proj_param12/ pi,half_pi,quarter_pi,rad,rrad
	save /proj_param12/

	double precision eps,tol
	common /proj_param13/ eps,tol
	save /proj_param13/

	logical debug
	common /proj_param21/ debug
	save /proj_param21/

c---------------------------------------------------------------- 

c	 include 'proj.h'

	 double precision lambda0,x0,xtrans0,ytrans0
	 common /gb_param2/ lambda0,x0,xtrans0,ytrans0

	 double precision a1,a2,a4,a6
	 parameter ( a1 = 111092.08210 , a2 = 16100.59187 )
	 parameter ( a4 = 16.96942 , a6 = 0.02226 )
	 double precision b1,b2,b4,b6
	 parameter ( b1 = 1. , b2 = 0.1449300705 )
	 parameter ( b4 = 0.0002138508 , b6 = 0.0000004322 )
	
	 double precision lambdap,phi,nu,nu1,xi,aux
	 double precision x,y,auxr
	 
	 x = east - x0 + xtrans0
	 y = north + ytrans0
	
	 !write(6,*) x,y,x0,xtrans0,ytrans0,east,north

	 aux = y / a1
	 auxr = aux * rad
	 xi = aux + b2*sin(2.*auxr) + b4*sin(4.*auxr) + b6*sin(6.*auxr)
	 xi = xi * rad
	 nu = sqrt( 1.+es*(cos(xi)**2) )
	
	 lambdap = atan( nu*sinh(x/aearth) / cos(xi) )
	 phi = atan( tan(xi) * cos(nu*lambdap) )
	
	 lon = lambdap/rad + lambda0
	 lat = phi/rad
	
	 if( .not. debug ) return
	 write(6,*) 'gb_c2g:'
	 write(6,*) aux,xi/rad,nu,x
	 write(6,*) lambdap/rad,phi/rad

	 end
	
c**********************************************************************
c**********************************************************************
c**********************************************************************
c**********************************************************************
c**********************************************************************

	subroutine lag2gb(xlag,ylag,e,n)

c transforms from lagoon coordinates to gauss-boaga

	implicit none

	double precision xlag,ylag	!lagoon coordinates
	double precision e,n		!gauss boaga

	double precision x0,y0
	parameter ( x0 = 2330000.-50000., y0 = 5000000. )

        e = xlag + x0
        n = ylag + y0

	end

c**********************************************************************

	subroutine gb2lag(e,n,xlag,ylag)

c transforms from gauss-boaga to lagoon coordinates

	implicit none

	double precision e,n		!gauss boaga
	double precision xlag,ylag	!lagoon coordinates

	double precision x0,y0
	parameter ( x0 = 2330000.-50000., y0 = 5000000. )

        xlag = e - x0
        ylag = n - y0

	end

c**********************************************************************

	subroutine lag2geo(xlag,ylag,lon,lat)

c transforms from lagoon coordinates to lon/lat

	implicit none

	double precision xlag,ylag	!lagoon coordinates
	double precision lon,lat	!lon/lat

	double precision e,n		!gauss boaga

	call lag2gb(xlag,ylag,e,n)
	call gb_c2g(e,n,lon,lat)

	end

c**********************************************************************

	subroutine geo2lag(lon,lat,xlag,ylag)

c transforms from lon/lat to lagoon coordinates

	implicit none

	double precision lon,lat	!lon/lat
	double precision xlag,ylag	!lagoon coordinates

	double precision e,n		!gauss boaga

	call gb_g2c(lon,lat,e,n)
	call gb2lag(e,n,xlag,ylag)

	end

c**********************************************************************
c**********************************************************************
c**********************************************************************
c**********************************************************************
c**********************************************************************
	
	 subroutine gb_sl_test
	
c test somma lombardo

c FIXME -> find expected values and insert them
	
	 implicit none
	
	 integer fuso
	 double precision deg,min,sec,sign
	 double precision lon_mm,phi_sl,lamb_sl,lon,lat,x,y

	 write(6,*) 'test somma lombardo...'

	 call dms2dec(12.d0,27.d0,8.4d0,1.d0,lon_mm)

	 fuso = 1
	 call gb_init(fuso)

	 call dms2dec(45.d0,41.d0,19.d0,1.d0,phi_sl)
	 call dms2dec(3.d0,43.d0,47.d0,-1.d0,lamb_sl)

	 lon = lamb_sl + lon_mm
	 lat = phi_sl
	
	 call gb_g2c(lon,lat,x,y)
	
	 write(6,*) lon,lat,x,y

	 x = 1477277.49
	 y = 5059027.45
	
	 call gb_c2g(x,y,lon,lat)

	 phi_sl = lat
	 lamb_sl = lon - lon_mm

	 write(6,*) x,y,lamb_sl,phi_sl
	 call dec2dms(lamb_sl,deg,min,sec,sign)
	 write(6,*) sign,deg,min,sec
	 call dec2dms(phi_sl,deg,min,sec,sign)
	 write(6,*) sign,deg,min,sec

	 write(6,*) 'end of test somma lombardo...'

	 end 
 
c****************************************

	subroutine deg_test

	implicit none

	integer i,idum,ndim
	double precision dc,dc2,geo2cent,eps
	double precision deg,min,sec,sign
	real proj_ran1
	
	write(6,*) 'test degrees...'

	idum = 98765
	eps = 1.e-8
	ndim = 10000
	
	do i=1,ndim
	   dc = 360. * proj_ran1(idum)
	   dc = dc - 180.
	   call dec2dms(dc,deg,min,sec,sign)
	   call dms2dec(deg,min,sec,sign,dc2)
	   if( abs(dc-dc2) .gt. eps ) then
	     write(6,*) '*** ',i,deg,min,sec,dc,dc2
	   end if
	   !write(6,*) i,deg,min,sec,dc,dc2
	end do

	write(6,*) 'finished testing ',ndim,' loops'
	write(6,*) 'end of test degrees...'

	 end
 
c****************************************

	 subroutine gb_rand_test

	 implicit none

	 integer i,idum,nerr,ndim,fuso
	 double precision eps,res,resmax
	 double precision x,y,lon,lat,lon1,lat1
	 real r
	 real proj_ran1

	 write(6,*) 'test random gauss-boaga...'

	 idum = 98765
	 eps = 2.e-7
	 resmax = 0.
	 nerr = 0
	 ndim = 10000
	 fuso = 0

	 do i=1,ndim
	   fuso = mod(fuso,2) + 1
	   r = 6. * proj_ran1(idum)
	   lon = 6.*fuso + r
	   r = 12. * proj_ran1(idum)
	   lat = 36. + r

	   call gb_init(fuso)
	   call gb_g2c(lon,lat,x,y)
	   call gb_c2g(x,y,lon1,lat1)

	   res = sqrt( (lat-lat1)**2 + (lon-lon1)**2 )
	   resmax = max(resmax,res)
	   if( res .gt. eps ) then
	     write(6,*) '*** ',i,fuso,lon,lon1,lat,lat1,res
	     write(6,*) '***   ',x,y
	     nerr = nerr + 1
	   end if
	   !write(6,*) i,fuso,lon,lon1,lat,lat1,res
	   !write(6,*) '      ',x,y
	 end do

	 write(6,*) 'finished testing ',ndim,' loops ',resmax
	 write(6,*) 'number of errors: ',nerr
	 write(6,*) 'end of test random gauss-boaga...'

	 end
 
c****************************************

	 subroutine gb_all_test
	 call gb_sl_test
	 call deg_test
	 call gb_rand_test
	 end

c****************************************

c	 program gb_main
c	 call gb_all_test
c	 end

c****************************************

c*******************************************************************
c
c Ellipsoid:
c
c Ellipsoid					major axis	1/f
c---------------------------------------------------------------------------
c Airy 1830   					6,377,563   	299.33
c Everest 1830  				6,377,276.3  	300.80
c Bessel 1841 					6,377,397.2 	299.15
c Clarke 1866 					6,378,206.4 	294.98
c Clarke 1880 					6,378,249.2 	293.47
c International 1924 				6,378,388 	297
c Krasovsky 1940 				6,378,245 	298.3
c International Astronomical Union 1968		6,378,160 	298.25
c WGS 72 (1972) 				6,378,135 	298.26
c GRS 80 (1980) 				6,378,137 	298.26
c WGS 84 (1984) 				6,378,137 	298.25722 
c
c Datum:
c
c Datum		Area			Origin			Ellipsoid
c----------------------------------------------------------------------------
c WGS 1984	Global			Earth center of mass	WGS 84
c NAD 1983	North America, Carib.	Earth center of mass	GRS 80
c NAD 1927	North America		Meades Ranch		Clarke 1866
c European 1950	Europe, N.Africa	Potsdam			International 
c
c*******************************************************************

	subroutine proj_set_ellipse(a,aux)

c sets parameters for ellipse
c
c as second parameter any of e, e2=es, f, 1/f, b can be given

	implicit none

	double precision a,aux


c---------------------------------------------------------------- 
c proj.h
c----------------------------------------------------------------

	double precision aearth,bearth,flat,rflat,e,es
	common /proj_param01/ aearth,bearth,flat,rflat,e,es
	save /proj_param01/

	double precision k0
	common /proj_param02/ k0
	save /proj_param02/

	double precision zero,one,two,four,half
	common /proj_param11/ zero,one,two,four,half
	save /proj_param11/

	double precision pi,half_pi,quarter_pi,rad,rrad
	common /proj_param12/ pi,half_pi,quarter_pi,rad,rrad
	save /proj_param12/

	double precision eps,tol
	common /proj_param13/ eps,tol
	save /proj_param13/

	logical debug
	common /proj_param21/ debug
	save /proj_param21/

c---------------------------------------------------------------- 

c	include 'proj.h'

	double precision phi

c first compute phi = b / a

	if( aux .gt. 6.10E6 ) then	!b
	  phi = aux/a
	else if( aux .ge. 290. ) then	!1/f
	  phi = 1. - 1./aux
	else if( aux .le. 0.004 ) then	!f
	  phi = 1. - aux
	else if( aux .le. 0.007 ) then	!e2 = es
	  phi = sqrt( 1. - aux )
	else if( aux .le. 0.09 ) then	!e
	  phi = sqrt( 1. - aux*aux )
	else
	  write(6,*) 'second variable: ',aux
	  stop 'error stop proj_set_ellipse: variable out of range'
	end if

c now compute all parameters

	aearth = a
	bearth = a * phi
	flat = 1. - phi
	rflat = 1. / flat
	es = 1. - phi*phi
	e = sqrt(es)

c final check

	if( rflat .lt. 290. .or. rflat .gt. 310. ) then
	  write(6,*) 'inverse flattening: ',rflat
	  stop 'error stop proj_set_ellipse: variable out of range'
	end if
	
	end

c*******************************************************************

	subroutine proj_set_scale_factor(k)

	implicit none

	double precision k


c---------------------------------------------------------------- 
c proj.h
c----------------------------------------------------------------

	double precision aearth,bearth,flat,rflat,e,es
	common /proj_param01/ aearth,bearth,flat,rflat,e,es
	save /proj_param01/

	double precision k0
	common /proj_param02/ k0
	save /proj_param02/

	double precision zero,one,two,four,half
	common /proj_param11/ zero,one,two,four,half
	save /proj_param11/

	double precision pi,half_pi,quarter_pi,rad,rrad
	common /proj_param12/ pi,half_pi,quarter_pi,rad,rrad
	save /proj_param12/

	double precision eps,tol
	common /proj_param13/ eps,tol
	save /proj_param13/

	logical debug
	common /proj_param21/ debug
	save /proj_param21/

c---------------------------------------------------------------- 

c	include 'proj.h'

	k0 = k

	end

c*******************************************************************

	subroutine proj_set_debug(bdebug)

	implicit none

	logical bdebug


c---------------------------------------------------------------- 
c proj.h
c----------------------------------------------------------------

	double precision aearth,bearth,flat,rflat,e,es
	common /proj_param01/ aearth,bearth,flat,rflat,e,es
	save /proj_param01/

	double precision k0
	common /proj_param02/ k0
	save /proj_param02/

	double precision zero,one,two,four,half
	common /proj_param11/ zero,one,two,four,half
	save /proj_param11/

	double precision pi,half_pi,quarter_pi,rad,rrad
	common /proj_param12/ pi,half_pi,quarter_pi,rad,rrad
	save /proj_param12/

	double precision eps,tol
	common /proj_param13/ eps,tol
	save /proj_param13/

	logical debug
	common /proj_param21/ debug
	save /proj_param21/

c---------------------------------------------------------------- 

c	include 'proj.h'

	debug = bdebug

	end

c*******************************************************************

	subroutine proj_print

	implicit none


c---------------------------------------------------------------- 
c proj.h
c----------------------------------------------------------------

	double precision aearth,bearth,flat,rflat,e,es
	common /proj_param01/ aearth,bearth,flat,rflat,e,es
	save /proj_param01/

	double precision k0
	common /proj_param02/ k0
	save /proj_param02/

	double precision zero,one,two,four,half
	common /proj_param11/ zero,one,two,four,half
	save /proj_param11/

	double precision pi,half_pi,quarter_pi,rad,rrad
	common /proj_param12/ pi,half_pi,quarter_pi,rad,rrad
	save /proj_param12/

	double precision eps,tol
	common /proj_param13/ eps,tol
	save /proj_param13/

	logical debug
	common /proj_param21/ debug
	save /proj_param21/

c---------------------------------------------------------------- 

c	include 'proj.h'

	write(6,*) 'proj_init:'
	write(6,*) 'aearth,bearth: ',aearth,bearth
	write(6,*) 'flat,rflat: ',flat,rflat
	write(6,*) 'e,es: ',e,es
	write(6,*) 'k0: ',k0

	end

c*******************************************************************

	subroutine proj_init

	implicit none


c---------------------------------------------------------------- 
c proj.h
c----------------------------------------------------------------

	double precision aearth,bearth,flat,rflat,e,es
	common /proj_param01/ aearth,bearth,flat,rflat,e,es
	save /proj_param01/

	double precision k0
	common /proj_param02/ k0
	save /proj_param02/

	double precision zero,one,two,four,half
	common /proj_param11/ zero,one,two,four,half
	save /proj_param11/

	double precision pi,half_pi,quarter_pi,rad,rrad
	common /proj_param12/ pi,half_pi,quarter_pi,rad,rrad
	save /proj_param12/

	double precision eps,tol
	common /proj_param13/ eps,tol
	save /proj_param13/

	logical debug
	common /proj_param21/ debug
	save /proj_param21/

c---------------------------------------------------------------- 

c	include 'proj.h'

	integer icall
	save icall
	data icall /0/

	if( icall .ne. 0 ) return
	icall = 1

	call proj_set_ellipse( 6378137. d0 , 298.25722 d0 )	!WGS 84 (1984)
	call proj_set_scale_factor( 1. d0 )

	zero	= 0.
	one	= 1.
	two	= 2.
	four	= 4.
	half	= 1./2.

	pi		= 4.*atan(one)
	half_pi		= pi/two
	quarter_pi	= half_pi/two
	rad		= pi / 180.
	rrad		= 1. / rad

	eps	= 1.d-10
	tol	= 1.d-7

	debug	= .true.
	debug	= .false.

	if( debug ) call proj_print

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************
c*******************************************************************
c*******************************************************************

	function phi2(ts)

	implicit none

	double precision phi2
	double precision ts


c---------------------------------------------------------------- 
c proj.h
c----------------------------------------------------------------

	double precision aearth,bearth,flat,rflat,e,es
	common /proj_param01/ aearth,bearth,flat,rflat,e,es
	save /proj_param01/

	double precision k0
	common /proj_param02/ k0
	save /proj_param02/

	double precision zero,one,two,four,half
	common /proj_param11/ zero,one,two,four,half
	save /proj_param11/

	double precision pi,half_pi,quarter_pi,rad,rrad
	common /proj_param12/ pi,half_pi,quarter_pi,rad,rrad
	save /proj_param12/

	double precision eps,tol
	common /proj_param13/ eps,tol
	save /proj_param13/

	logical debug
	common /proj_param21/ debug
	save /proj_param21/

c---------------------------------------------------------------- 

c	include 'proj.h'

	integer i,niter
	parameter(niter = 15)
	double precision phi,eccnth,con,aux,dphi

	eccnth = half * e
	phi = half_pi - two * atan(ts)

	do i=1,niter
	  con = e * sin(phi)
	  aux = ( (one-con)/(one+con) ) ** eccnth
	  dphi = half_pi - two * atan(ts*aux) - phi
	  phi = phi + dphi
	  if( debug ) write(6,*) 'iterating: ',phi,dphi,i
	  if( abs(dphi) .lt. eps ) then
	    phi2 = phi
	    return
	  end if
	end do

	stop 'error stop phi2: could not converge'
	end

c*******************************************************************

	function tsfn(phi)

	implicit none

	double precision tsfn
	double precision phi


c---------------------------------------------------------------- 
c proj.h
c----------------------------------------------------------------

	double precision aearth,bearth,flat,rflat,e,es
	common /proj_param01/ aearth,bearth,flat,rflat,e,es
	save /proj_param01/

	double precision k0
	common /proj_param02/ k0
	save /proj_param02/

	double precision zero,one,two,four,half
	common /proj_param11/ zero,one,two,four,half
	save /proj_param11/

	double precision pi,half_pi,quarter_pi,rad,rrad
	common /proj_param12/ pi,half_pi,quarter_pi,rad,rrad
	save /proj_param12/

	double precision eps,tol
	common /proj_param13/ eps,tol
	save /proj_param13/

	logical debug
	common /proj_param21/ debug
	save /proj_param21/

c---------------------------------------------------------------- 

c	include 'proj.h'

	double precision aux1,aux2

	aux1 = e * sin(phi)
	aux2 = ( one - aux1 ) / ( one + aux1 )
	if( e .le. zero ) aux2 = one

	tsfn = tan( half * (half_pi - phi) ) / aux2 ** (half*e)

	end

c*******************************************************************

	function msfn(phi)

	implicit none

	double precision msfn
	double precision phi


c---------------------------------------------------------------- 
c proj.h
c----------------------------------------------------------------

	double precision aearth,bearth,flat,rflat,e,es
	common /proj_param01/ aearth,bearth,flat,rflat,e,es
	save /proj_param01/

	double precision k0
	common /proj_param02/ k0
	save /proj_param02/

	double precision zero,one,two,four,half
	common /proj_param11/ zero,one,two,four,half
	save /proj_param11/

	double precision pi,half_pi,quarter_pi,rad,rrad
	common /proj_param12/ pi,half_pi,quarter_pi,rad,rrad
	save /proj_param12/

	double precision eps,tol
	common /proj_param13/ eps,tol
	save /proj_param13/

	logical debug
	common /proj_param21/ debug
	save /proj_param21/

c---------------------------------------------------------------- 

c	include 'proj.h'

	double precision sinphi, cosphi

	sinphi = sin(phi)
	cosphi = cos(phi)

	msfn = cosphi / sqrt (one - es * sinphi * sinphi)

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine dec2dms(a,d,m,s,sign)

	implicit none

	double precision a,d,m,s,sign

	double precision aa,aux
	double precision zero,one,sixty
	parameter(zero=0.,one=1.,sixty=60.)

	aa = abs(a)
	d = int(aa)
	aux = sixty*(aa-d)
	m = int(aux)
	s = sixty*(aux-m)

	sign = one
	if( a .lt. zero ) sign = -sign

	end

c*******************************************************************

	subroutine dms2dec(d,m,s,sign,a)

c transforms deg min sec to decimal degrees

	double precision d,m,s,sign,a

	double precision sixty,sixty2
	parameter(sixty=60.,sixty2=3600.)

	a = d + m/sixty + s/sixty2
	a = a * sign

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************
c*******************************************************************
c*******************************************************************
	
	function proj_ran1(idum)

	implicit none

	real proj_ran1
	integer idum	! idum < 0 => reset

	integer m1,m2,m3,ia1,ia2,ia3,ic1,ic2,ic3
	real rm1,rm2
	parameter (m1=259200,ia1=7141,ic1=54773,rm1=3.8580247e-6)
	parameter (m2=134456,ia2=8121,ic2=28411,rm2=7.4373773e-6)
	parameter (m3=243000,ia3=4561,ic3=51349)
	integer ndim
	parameter (ndim=97)

	integer ix1,ix2,ix3
	real r(ndim)
	save ix1,ix2,ix3,r
	integer j

	integer iff
	save iff
	data iff /0/

	if ( idum.lt.0 .or. iff.eq.0 ) then
	  iff=1
	  ix1=mod(abs(ic1-idum),m1)
	  ix1=mod(ia1*ix1+ic1,m1)
	  ix2=mod(ix1,m2)
	  ix1=mod(ia1*ix1+ic1,m1)
	  ix3=mod(ix1,m3)
	  do j=1,ndim
	    ix1=mod(ia1*ix1+ic1,m1)
            ix2=mod(ia2*ix2+ic2,m2)
            r(j)=(float(ix1)+float(ix2)*rm2)*rm1
	  end do
	endif

	ix1=mod(ia1*ix1+ic1,m1)
	ix2=mod(ia2*ix2+ic2,m2)
	ix3=mod(ia3*ix3+ic3,m3)
	j=1+(ndim*ix3)/m3
	if(j.gt.ndim.or.j.lt.1) then
	  write(6,*) j
	  stop 'error stop proj_ran1: internal error'
	end if

	proj_ran1=r(j)
	r(j)=(float(ix1)+float(ix2)*rm2)*rm1

	end

c*******************************************************************

	function proj_ran_intv(idum,rmin,rmax)

	implicit none

	real proj_ran_intv
	integer idum
	real rmin,rmax

	real proj_ran1

	proj_ran_intv = rmin + proj_ran1(idum) * (rmax-rmin)

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************
c*******************************************************************
c*******************************************************************

c        PIN->e = sqrt(PIN->es);
c        PIN->ra = 1. / PIN->a;
c        PIN->one_es = 1. - PIN->es;
c        if (PIN->one_es == 0.) { pj_errno = -6; goto bum_call; }
c        PIN->rone_es = 1./PIN->one_es;

c pj_tsfn(double phi, double sinphi, double e) {
c         sinphi *= e;
c         return (tan (.5 * (HALFPI - phi)) /
c            pow((1. - sinphi) / (1. + sinphi), .5 * e));
c }

c pj_msfn(double sinphi, double cosphi, double es) {
c         return (cosphi / sqrt (1. - es * sinphi * sinphi));
c }

c #define TOL 1.0e-10
c #define N_ITER 15
c 
c         double
c pj_phi2(double ts, double e) {
c         double eccnth, Phi, con, dphi;
c         int i;
c 
c         eccnth = .5 * e;
c         Phi = HALFPI - 2. * atan (ts);
c         i = N_ITER;
c         do {
c                 con = e * sin (Phi);
c                 dphi = HALFPI - 2. * atan (ts * pow((1. - con) /
c                    (1. + con), eccnth)) - Phi;
c                 Phi += dphi;
c         } while ( fabs(dphi) > TOL && --i);
c         if (i <= 0)
c                 pj_errno = -18;
c         return Phi;
c }

c*******************************************************************

	subroutine proj_all_test
	end

c*******************************************************************

c	program proj_main
c	call proj_all_test
c	end

c*******************************************************************


	subroutine all_test
	write(6,*) '================================================='
	call genproj_all_test
	write(6,*) '================================================='
	call gb_all_test
	write(6,*) '================================================='
	call proj_all_test
	write(6,*) '================================================='
	end

c	program main_test
c	call all_test
c	end
