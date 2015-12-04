!
! $Id: subcoord.f,v 1.5 2010-02-26 17:35:06 georg Exp $
!
! coordinate handling routines
!
! revision log :
!
! 07.05.2009    ggu     new framework (to be finished), new basic projection
! 26.05.2010    ggu     new utm, bs renamed to cpp
! 07.12.2010    ggu     bug fix in convert_coords for mode=-1
! 16.02.2011    ggu     new proj=4 (non standard UTM) implemented
! 18.11.2011    ggu     in CCP changed order of params in c_param
! 30.09.2015    ccf     introduced module projection
! 03.12.2015    ggu     bug fix because arrays were defined length 1
!
! usage :
!
!	example is for Gauss-Boaga, other projections work similar
!	example only converts one point
!
!	params(1) = 2
!	params(2) = 0.
!	params(3) = 0.
!	call init_coords(1,params)
!
!	mode = -1			!from geo to GB
!	n = 1
!	xg = 12.
!	yg = 45.
!	call convert_coords(mode,n,xc,yc,xg,yg)
!
!********************************************************************

!==================================================================
        module projection
!==================================================================

        integer, save		     	     :: proj = 0

        double precision, private, parameter :: zero  = 0.d0
        double precision, private, parameter :: one   = 1.d0
        double precision, private, parameter :: two   = 2.d0
        double precision, private, parameter :: four  = 4.d0
        double precision, private, parameter :: half  = one/two
        double precision, private, parameter :: pi    = four*atan(one)
        double precision, private, parameter :: half_pi = pi/two
        double precision, private, parameter :: rad   = pi / 180.d0
        logical, private		     :: debug = .FALSE. 

        double precision, private	     :: k0
        double precision, private	     :: aearth
        double precision, private	     :: bearth
        double precision, private     	     :: flat
        double precision, private	     :: rflat
        double precision, private	     :: es
        double precision, private	     :: e
!GB
        integer, private, save  	     :: fuso		!(1-2)
        double precision, private	     :: lambda0
        double precision, private	     :: x0
        double precision, private	     :: xtrans0
        double precision, private	     :: ytrans0
!UTM
        integer, private, save               :: utm_zone	!(1-60)
        double precision, private	     :: ep2,dn,es_4,es_6
        double precision, private	     :: z1,z2,z3,z4
        double precision, private	     :: j1,j2,j3,j4
!CPP
        double precision, private	     :: lon_0,lat_0
        double precision, private	     :: phi_0
        double precision, private	     :: xfact,yfact
!LAG
        double precision,private,parameter :: x0lag = 2330000.-50000.
        double precision,private,parameter :: y0lag = 5000000.

!==================================================================
        contains
!==================================================================

!********************************************************************
!********************************************************************
!
! ROUTINES FOR CPP TRANSFORMATIONS
!
! equidistant cylindrical projection
!
! carte parallelogrammatique projection or CPP (Marinus of Tyre)
!
!********************************************************************

        subroutine cpp_init(lon0,lat0,phi)

        implicit none

        double precision :: lon0,lat0,phi

        lon_0 = lon0
        lat_0 = lat0
        phi_0 = phi

        yfact = 60 * 1852               !distance of 1 min = 1 nautical mile
        xfact = yfact * cos(phi*rad)

        end subroutine cpp_init

!********************************************************************

        subroutine cpp_g2c(lon,lat,x,y)

! transformation from (lon/lat) to cartesian (x,y)

        implicit none

        double precision :: lon,lat !geographical coordinates
        double precision :: x,y     !cartesian coordinates (return)

        x = (lon-lon_0) * xfact
        y = (lat-lat_0) * yfact

        end subroutine cpp_g2c

!********************************************************************

        subroutine cpp_c2g(x,y,lon,lat)

! transformation from to cartesian (x,y) to (lon/lat)

        implicit none

        double precision :: x,y  !cartesian coordinates
        double precision :: lon,lat !geographical coordinates (return)

        lon =  lon_0 + x / xfact
        lat =  lat_0 + y / yfact

        end subroutine cpp_c2g

!********************************************************************
!********************************************************************
!
! ROUTINES FOR GAUSS-BOAGA TRANSFORMATIONS
!
! contents :
!
! gauss-boaga core routines
!
! subroutine gb_init(fuse)
!		intializes gaus-boaga routines
! subroutine gb_trans(xtrans,ytrans)
!		intializes translation of gauss-boaga coordinates
! subroutine gb_g2c(lon,lat,east,north)
!		transformation from (lon/lat) to gauss-boaga (east,north)
! subroutine gb_c2g(east,north,lon,lat)
!		transformation from gauss-boaga (east,north) to (lon/lat)
!
! helper routines
!
! subroutine lag2gb(xlag,ylag,e,n)
!		transforms from lagoon coordinates to gauss-boaga
! subroutine gb2lag(e,n,xlag,ylag)
!		transforms from gauss-boaga to lagoon coordinates
! subroutine lag2geo(xlag,ylag,lon,lat)
!		transforms from lagoon coordinates to lon/lat
! subroutine geo2lag(lon,lat,xlag,ylag)
!		transforms from lon/lat to lagoon coordinates
!
! ...other test routines...
!
! usage :
!
! gb_init(fuse) must be called at beginning to set fuse
! gb_trans(xtrans,ytrans) may be called if different origin is desired
!	if not called no translation will be performed
!
! example: in order to pass from lat/lon to lagoon coordinates use
!
!	call gb_init(2)
!	call gb_trans(2280000.D+0,5000000.D+0)
!
! and then
!
!	call gb_g2c(lon,lat,east,north)
!
! notes :
!
! gauss-boaga transformations
!
! fuso I  = W = 32 UTM   meridiano centrale  9 deg
! fuso II = E = 33 UTM   meridiano centrale 15 deg
!
! scale 0.9996
!
! UTM  uses ED50 (European Data 1950)
! GB (IGM) uses Roma40
!
! both use international ellipsoid (Hayford)
!
! ellipsoid:  c = 6397376.633  es = e2 = 0.0067681702
!
! Monte Mario:
!
! Roma40: 41 55 25.510 12 27 08.400
! Ed50:  41 55 31.487 12 27 08.933
!
!------------------------------------------------------------
!
! for Venice lagoon:venlag62_new2.grd
!
! fuso = 2

! to get GB coordinates (x_gb,y_gb) from venlag62 (x_lag,y_lag) compute:
!
! x_gb = x_lag + x0lag
! y_gb = y_lag + y0lag
!
! with
!
! x0 = 2330000. - 50000. = 2280000;
! y0 = 5000000.;
!
! revision log :
!
! 16.02.2006    ggu     separated from other routines and commented
! 11.03.2009    ggu     cleaned, easier initialization
!
!**********************************************************************

        subroutine gb_init(fuse)

! intializes gaus-boaga routines
!
! to obtain E/N use E = x + x0    N = y
! with x0 = 1500 km (fuse=1) and x0 = 2520 km (fuse=2)

        implicit none

        integer :: fuse			!1 or 2

        call proj_init

        fuso = fuse

        select case (fuso)
          case (1)
            lambda0 = 9.
            x0 = 1500000.
          case(2)
            lambda0 = 15.
            x0 = 2520000.
          case default
            stop 'error stop gb_init: fuse'
        end select

        xtrans0 = 0.		! change origin of gb coordinates
        ytrans0 = 0.

! use own constants for international ellipsoid (Hayford)

        aearth	= 6397376.633
        k0	= 0.9996
        es	= 0.0067681702
        e	= sqrt(es)

        if( .NOT. debug ) return
        write(6,*) 'gb_init:'
        write(6,*) fuso,lambda0,x0
        write(6,*) aearth,k0
        write(6,*) e,es

        end subroutine gb_init

!**********************************************************************

        subroutine gb_trans(xtrans,ytrans)

! intializes translation of gauss-boaga coordinates
!
! (xtrans,ytrans) are the coordinates of the GB system that
! will be the new origin of the translation
! -> so these coordinates will be subtracted from the GB coordinates

        implicit none

        double precision :: xtrans,ytrans

        xtrans0 = xtrans
        ytrans0 = ytrans

        end subroutine gb_trans

!**********************************************************************

        subroutine gb_g2c(lon,lat,east,north)

! transformation from (lon/lat) to gauss-boaga (east,north)

        implicit none

        double precision :: lon,lat !geographical coordinates
        double precision :: east,north  !coordinates gauss-boaga (return)

        double precision :: a1,a2,a4,a6
        parameter ( a1 = 111092.08210 , a2 = 16100.59187 )
        parameter ( a4 = 16.96942 , a6 = 0.02226 )
        	
        double precision :: lambdap,phi,nu,nu1,xi,aux
!double precision c,e2
        double precision :: x,y
        	 
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

        if( .NOT. debug ) return
        write(6,*) 'gb_g2c:'
        write(6,*) lambdap/rad,lon
        write(6,*) nu1,xi/rad,nu
        write(6,*) aux,x,y

        end subroutine gb_g2c
        	
!**********************************************************************
        	
        subroutine gb_c2g(east,north,lon,lat)
        	
! transformation from gauss-boaga (east,north) to (lon/lat)
        	
        implicit none
        	
        double precision :: east,north  !coordinates gauss-boaga
        double precision :: lon,lat !geographical coordinates (return)
        	
        double precision :: a1,a2,a4,a6
        parameter ( a1 = 111092.08210 , a2 = 16100.59187 )
        parameter ( a4 = 16.96942 ,     a6 = 0.02226 )
        double precision :: b1,b2,b4,b6
        parameter ( b1 = 1. ,           b2 = 0.1449300705 )
        parameter ( b4 = 0.0002138508 , b6 = 0.0000004322 )
        	
        double precision :: lambdap,phi,nu,nu1,xi,aux
        double precision :: x,y,auxr
        	 
        x = east - x0 + xtrans0
        y = north + ytrans0
        	
        aux = y / a1
        auxr = aux * rad
        xi = aux + b2*sin(2.*auxr) + b4*sin(4.*auxr) + b6*sin(6.*auxr)
        xi = xi * rad
        nu = sqrt( 1.+es*(cos(xi)**2) )
        	
        lambdap = atan( nu*sinh(x/aearth) / cos(xi) )
        phi = atan( tan(xi) * cos(nu*lambdap) )
        	
        lon = lambdap/rad + lambda0
        lat = phi/rad
        	
        if( .NOT. debug ) return
        write(6,*) 'gb_c2g:'
        write(6,*) aux,xi/rad,nu,x
        write(6,*) lambdap/rad,phi/rad

        end subroutine gb_c2g
        	
!**********************************************************************

        subroutine lag2gb(xlag,ylag,xgb,ygb)

! transforms from lagoon coordinates to gauss-boaga

        implicit none

        double precision :: xlag,ylag	!lagoon coordinates
        double precision :: xgb,ygb	!gauss boaga

        xgb = xlag + x0lag
        ygb = ylag + y0lag

        end subroutine lag2gb

!**********************************************************************

        subroutine gb2lag(xgb,ygb,xlag,ylag)

! transforms from gauss-boaga to lagoon coordinates

        implicit none

        double precision :: xgb,ygb	!gauss boaga
        double precision :: xlag,ylag	!lagoon coordinates

        xlag = xgb - x0lag
        ylag = ygb - y0lag

        end subroutine gb2lag

!**********************************************************************

        subroutine lag2geo(xlag,ylag,lon,lat)

! transforms from lagoon coordinates to lon/lat

        implicit none

        double precision :: xlag,ylag	!lagoon coordinates
        double precision :: lon,lat	!lon/lat

        double precision :: xgb,ygb	!gauss boaga

        call lag2gb(xlag,ylag,xgb,ygb)
        call gb_c2g(xgb,ygb,lon,lat)

        end subroutine lag2geo

!**********************************************************************

        subroutine geo2lag(lon,lat,xlag,ylag)

! transforms from lon/lat to lagoon coordinates

        implicit none

        double precision :: lon,lat	!lon/lat
        double precision :: xlag,ylag	!lagoon coordinates

        double precision :: xgb,ygb	!gauss boaga

        call gb_g2c(lon,lat,xgb,ygb)
        call gb2lag(xgb,ygb,xlag,ylag)

        end subroutine geo2lag

!**********************************************************************
!**********************************************************************
!
! ROUTINES FOR UTM TRANSFORMATIONS
!
! contents :
!
! utm core routines
!
! subroutine utm_init(zone)
!		intializes gaus-boaga routines
!
! subroutine utm_trans(xtrans,ytrans)
!		intializes translation of utm coordinates
! subroutine utm_set_zone(zone)
!		changes zone of UTM - must have already been initialized
!
! subroutine utm_compute_zone(lon,zone)
!		computes and sets zone of UTM given longitude
! subroutine utm_compute_band(lat,band)
!		computes band of UTM given latitude
! subroutine utm_east_north(bsouth)
!		use false easting/northing (bsouth true for south)
!
! subroutine utm_g2c(lon,lat,x,y)
!		transformation from (lon/lat) to utm (east,north)
! subroutine utm_c2g(x,y,lon,lat)
!		transformation from utm (east,north) to (lon/lat)
!
! ...other test routines...
!
! usage :
!
! utm_init(zone) must be called at beginning to set zone and initialize
!	auxiliary variables
! utm_trans(xtrans,ytrans) may be called if different origin is desired
!	if not called no translation will be performed
!
! no false easting and northing is performed in order to non decrease
! accuracy. if you want to have traditional UTM coordinates you have to
! do one of the following:
!
! 1)
!	add 500000 to x
!	if in southern hemisphere add 10000000 to y
! 2)
!	call utm_trans(-500000.d+0,+0.d+0)	!northern hemisphere
!	call utm_trans(-500000.d+0,-1.d+7)	!southern hemisphere
! 3)
!	call utm_east_north(.false.)		!northern hemisphere
!	call utm_east_north(.true.)		!southern hemisphere
!
! example: in order to pass from lat/lon to x/y close to Venice lagoon:
!
!
! Venice is in zone 33, band T
!
!	call utm_init(33)
!	call utm_east_north(.false.)	!if false easting/northing is desired
!
!	call utm_g2c(lon,lat,x,y)
!
! revision log :
!
! 22.04.2010    ggu     written starting from GB routines
! 25.04.2010    ggu     finished and tested
!
!**********************************************************************

        subroutine utm_init(zone)

! standard call to initialize UTM routines

        implicit none

        integer :: zone			!UTM zone [1-60]

        call utm_set_zone(zone)
        call utm_init_internal

        end subroutine utm_init

!**********************************************************************

        subroutine utm_init_nonstd(lambda)

! non-standard call to initialize UTM routines -> give central meridian

        implicit none

        double precision :: lambda		!central meridian

        call utm_set_central_meridian(lambda)
        call utm_init_internal

        end subroutine utm_init_nonstd

!**********************************************************************

        subroutine utm_init_internal

! intializes utm routines

        implicit none

        double precision :: e1

        call proj_init
!call utm_set_zone(zone)	!done outside

        k0 = 0.9996

        xtrans0 = 0.		! change origin of utm coordinates
        ytrans0 = 0.

        ep2 = es / (1.-es)
        dn = (aearth-bearth)/(aearth+bearth)

        es_4 = es*es
        es_6 = es_4*es

        z1 = 1. - es/4. - 3.*es_4/64. - 5.*es_6/256.
        z2 = 3.*es/8. + 3.*es_4/32. + 45.*es_6/1024.
        z3 = 15.*es_4/256. + 45.*es_6/1024.
        z4 = 35.*es_6/3072

        e1 = (1.-sqrt(1.-es)) / (1.+sqrt(1.-es))

        j1 = 3.*e1/2. - 27.*(e1**3)/32.
        j2 = 21.*(e1**2)/16. - 55.*(e1**4)/32.
        j3 = 151.*(e1**3)/96.
        j4 = 1097.*(e1**4)/512.

        if( .NOT. debug ) return
        write(6,*) 'utm_init:'
        write(6,*) lambda0
        write(6,*) aearth,k0
        write(6,*) e,es

        end subroutine utm_init_internal

!**********************************************************************

        subroutine utm_set_scale_factor(k)

! sets scale factor in UTM (for non standard applications)

        implicit none

        double precision :: k

        call proj_set_scale_factor(k)

        end subroutine utm_set_scale_factor

!**********************************************************************

        subroutine utm_trans(xtrans,ytrans)

! intializes translation of utm coordinates
!
! (xtrans,ytrans) are the coordinates of the GB system that
! will be the new origin of the translation
! -> so these coordinates will be subtracted from the GB coordinates

        implicit none

        double precision :: xtrans,ytrans

        xtrans0 = xtrans
        ytrans0 = ytrans

        end subroutine utm_trans

!**********************************************************************

        subroutine utm_set_zone(zone)

! changes zone of UTM - must have already been initialized

        implicit none

        integer :: zone

        if( zone < 1 .OR. zone > 60 ) then
            write(6,*) 'zone = ',zone
            stop 'error stop utm_set_zone: zone out of bounds [1-60]'
        end if

        utm_zone = zone
        lambda0 = (zone-1)*6 - 180 + 3	!3 centers inside zone

        end subroutine utm_set_zone

!**********************************************************************

        subroutine utm_set_central_meridian(lambda)

! sets central meridian directly (for non standard UTM projections)

        implicit none

        double precision :: lambda

        utm_zone = 0		! not valid
        lambda0 = lambda

        end subroutine utm_set_central_meridian

!**********************************************************************

        subroutine utm_compute_zone(lon,zone)

! computes and sets zone of UTM given longitude

        implicit none

        double precision :: lon
        integer :: zone

        integer :: long

        if( lon < -180. .OR. lon > 180. ) then
            write(6,*) 'lon = ',lon
            stop 'error stop utm_compute_zone: lon out of bounds'
        end if

        long = 180 + lon
        zone = min(60,1+long/6)

        call utm_set_zone(zone)

        end subroutine utm_compute_zone

!**********************************************************************

        subroutine utm_compute_band(lat,band)

! computes band of UTM given latitude (no exceptions)

        implicit none

        double precision :: lat
        character(1) :: band

        integer :: lati,index

        character(20) :: bands
        save bands
        data bands /'CDEFGHJKLMNPQRSTUVWX'/

        if( lat < -80. .OR. lat > 84. ) then
            write(6,*) 'lat = ',lat
            stop 'error stop utm_compute_band: lat out of bounds'
        end if

        lati = 80 + lat
        index = min(20,1+lati/8)

        band = bands(index:index)

        end subroutine utm_compute_band

!**********************************************************************

        subroutine utm_east_north(bsouth)

! use false easting/northing (bsouth must be true for southern hemisphere)

        implicit none

        logical :: bsouth

        if( bsouth ) then
            call utm_trans( -500000.d+0 , -1.d+7 )	!southern hemisphere
        else
            call utm_trans( -500000.d+0 , 0.d+0 )		!northern hemisphere
        end if

        end subroutine utm_east_north

!**********************************************************************

        subroutine utm_g2c(lon,lat,x,y)

! transformation from (lon/lat) to utm (x,y)

        implicit none

        double precision :: lon,lat	!geographical coordinates
        double precision :: x,y		!coordinates UTM (return)

        double precision :: long,phi,p,nu,M
        double precision :: k1,k2,k3,k4,k5
        	 
        long = lon
        if( long > 180. ) long = long - 180.
        phi = rad * lat
        p = rad * (long - lambda0)
        nu = aearth / sqrt( 1.-es*(sin(phi)**2) )
        	
        M = z1*phi - z2*sin(2.*phi) + z3*sin(4.*phi) - z4*sin(6.*phi)
        M = aearth * M

        k1 = M * k0
        k2 = k0 * nu * sin(2.*phi)/4.
        k3 = k0 * nu * sin(phi)*((cos(phi)**3)/24.) 
     +       * (5. - tan(phi)**2 + 9.*ep2*cos(phi)**2 
     +       + 4.*ep2*ep2*cos(phi)**4)
        k4 = k0 * nu * cos(phi)
        k5 = k0 * nu * ((cos(phi)**3)/6.) 
     +       * (1. - tan(phi)**2 + ep2*cos(phi)**2)

        y = k1 + k2*p**2 + k3*p**4
        x = k4*p + k5*p**3

        x = x - xtrans0
        y = y - ytrans0

        if( .NOT. debug ) return
        write(6,*) '======================='
        write(6,*) 'utm_g2c:'
        write(6,*) 'lambda0: ',lambda0
        write(6,*) 'lat: ',lat
        write(6,*) 'lon: ',lon
        write(6,*) 'phi: ',phi
        write(6,*) 'ep2: ',ep2
        write(6,*) 'es: ',es
        write(6,*) 'nu: ',nu
        write(6,*) 'M: ',M
        write(6,*) '======================='

        end subroutine utm_g2c
        	
!**********************************************************************
        	
        subroutine utm_c2g(x,y,lon,lat)
        	
! transformation from utm (x,y) to (lon/lat)
        	
        implicit none
        	
        double precision :: x,y  !coordinates utm
        double precision :: lon,lat !geographical coordinates (return)
        	
        double precision :: xx,yy,M,mu,fp
        double precision :: c1,t1,r1,n1
        double precision :: d,d2,d4,d6
        double precision :: q1,q2,q3,q4,q5,q6,q7
        double precision :: rrad 

        rrad = 1.d0 / rad

        xx = x  + xtrans0
        yy = y  + ytrans0
        	
        M = yy / k0
        mu = M / (aearth * z1)

        fp = mu + j1*sin(2.*mu) + j2*sin(4.*mu) +
     +       j3*sin(6.*mu) + j4*sin(8.*mu)

        c1 = ep2 * cos(fp)**2
        t1 = tan(fp)**2
        r1 = aearth*(1.-es)/(1.-es*sin(fp)**2)**(1.5)
        n1 = aearth/sqrt(1.-es*sin(fp)**2)

        d  = xx / (n1*k0)
        d2 = d*d
        d4 = d2*d2
        d6 = d4*d2
        	
        q1 = n1 * tan(fp) / r1
        q2 = d2/2.
        q3 = (5. + 3.*t1 + 10.*c1 - 4.*c1*c1 - 9.*ep2)*d4/24.
        q4 = (61. + 90.*t1 + 298.*c1 + 45.*t1*t1 -
     +       3.*c1*c1 - 252.*ep2) * d6 / 720.
        q5 = d
        q6 = (1. + 2.*t1 + c1) * d2*d / 6.
        q7 = (5. - 2.*c1 + 28.*t1 - 3.*c1*c1 +
     +       8.*ep2 + 24.*t1*t1) * d4 * d / 120.

        lat = rrad * (fp - q1*(q2-q3+q4))
        lon = lambda0 + rrad * (q5-q6+q7)/cos(fp)
        	
        if( .NOT. debug ) return
        write(6,*) '======================='
        write(6,*) 'utm_c2g:'
        write(6,*) 'M: ',M
        write(6,*) 'mu: ',mu
        write(6,*) 'fp: ',fp
        write(6,*) 'n1: ',n1
        write(6,*) 't1: ',t1
        write(6,*) 'r1: ',r1
        write(6,*) 'c1: ',c1
        write(6,*) '======================='

        end subroutine utm_c2g
        	
!**********************************************************************
!**********************************************************************
!**********************************************************************
!*******************************************************************
!
! Ellipsoid:
!
! Ellipsoid					major axis	1/f
!---------------------------------------------------------------------------
! Airy 1830   					6,377,563   	299.33
! Everest 1830  				6,377,276.3  	300.80
! Bessel 1841 					6,377,397.2 	299.15
! Clarke 1866 					6,378,206.4 	294.98
! Clarke 1880 					6,378,249.2 	293.47
! International 1924 				6,378,388 	297
! Krasovsky 1940 				6,378,245 	298.3
! International Astronomical Union 1968		6,378,160 	298.25
! WGS 72 (1972) 				6,378,135 	298.26
! GRS 80 (1980) 				6,378,137 	298.25722
! WGS 84 (1984) 				6,378,137 	298.25722
!
! Datum:
!
! Datum		Area			Origin			Ellipsoid
!----------------------------------------------------------------------------
! WGS 1984	Global			Earth center of mass	WGS 84
! NAD 1983	North America, Carib.	Earth center of mass	GRS 80
! NAD 1927	North America		Meades Ranch		Clarke 1866
! European 1950	Europe, N.Africa	Potsdam			International
!
!*******************************************************************

        subroutine proj_set_ellipse(a,aux)

! sets parameters for ellipse
!
! as second parameter any of e, e2=es, f, 1/f, b can be given
!
! a		equatorial radius [m]
! b		polar radius [m]
! phi		phi = b/a
! f		flattening, f = (a-b)/a = 1 - phi
! e2 (=es)	eccentricity squared, e2 = 1 - (b/a)**2
! e		eccentricity, e = sqrt(e2)

        implicit none

        double precision :: a,aux
        double precision :: phi

! first compute phi = b / a

        if( aux > 6.10E6 ) then	!b
            phi = aux/a
        else if( aux >= 290. ) then	!1/f
            phi = 1. - 1./aux
        else if( aux <= 0.004 ) then	!f
            phi = 1. - aux
        else if( aux <= 0.007 ) then	!e2 = es
            phi = sqrt( 1. - aux )
        else if( aux <= 0.09 ) then	!e
            phi = sqrt( 1. - aux*aux )
        else
            write(6,*) 'second variable: ',aux
            stop 'error stop proj_set_ellipse: 2. var out of range'
        end if

! now compute all parameters

        aearth = a
        bearth = a * phi
        flat = 1. - phi
        rflat = 1. / flat
        es = 1. - phi*phi
        e = sqrt(es)

! final check

        if( rflat < 290. .OR. rflat > 310. ) then
            write(6,*) 'inverse flattening: ',rflat
            stop 'error stop proj_set_ellipse: variable out of range'
        end if
        	
        end subroutine proj_set_ellipse

!*******************************************************************

        subroutine proj_set_scale_factor(k)

        implicit none

        double precision :: k

        k0 = k

        end subroutine proj_set_scale_factor

!*******************************************************************

        subroutine proj_set_debug(bdebug)

        implicit none

        logical :: bdebug

        debug = bdebug

        end subroutine proj_set_debug

!*******************************************************************

        subroutine proj_print

        implicit none

        write(6,*) 'proj_init:'
        write(6,*) 'aearth,bearth: ',aearth,bearth
        write(6,*) 'flat,rflat: ',flat,rflat
        write(6,*) 'e,es: ',e,es
        write(6,*) 'k0: ',k0

        end subroutine proj_print

!*******************************************************************

        subroutine proj_init

        implicit none

        integer, save :: icall
        data icall /0/

        if( icall /= 0 ) return
        icall = 1

        call proj_set_ellipse( 6378137.d0 , 298.25722d0 )	!WGS 84 (1984)
        call proj_set_scale_factor( 1.d0 )

        if( debug ) call proj_print

        end subroutine proj_init

!*******************************************************************
!*******************************************************************
!*******************************************************************
!*******************************************************************
!*******************************************************************

        function phi2(ts)

        implicit none

        double precision :: phi2
        double precision :: ts

        integer :: i,niter
        parameter(niter = 15)
        double precision :: phi,eccnth,con,aux,dphi
        double precision, parameter :: eps   = 1.d-10

        eccnth = half * e
        phi = half_pi - two * atan(ts)

        do i=1,niter
            con = e * sin(phi)
            aux = ( (one-con)/(one+con) ) ** eccnth
            dphi = half_pi - two * atan(ts*aux) - phi
            phi = phi + dphi
            if( debug ) write(6,*) 'iterating: ',phi,dphi,i
            if( abs(dphi) < eps ) then
                phi2 = phi
                return
            end if
        end do

        stop 'error stop phi2: could not converge'

        end function phi2

!*******************************************************************

        function tsfn(phi)

        implicit none

        double precision :: tsfn
        double precision :: phi

        double precision :: aux1,aux2

        aux1 = e * sin(phi)
        aux2 = ( one - aux1 ) / ( one + aux1 )
        if( e <= zero ) aux2 = one

        tsfn = tan( half * (half_pi - phi) ) / aux2 ** (half*e)

        end function tsfn

!*******************************************************************

        function msfn(phi)

        implicit none

        double precision :: msfn
        double precision :: phi
        double precision :: sinphi, cosphi

        sinphi = sin(phi)
        cosphi = cos(phi)

        msfn = cosphi / sqrt (one - es * sinphi * sinphi)

        end function msfn

!*******************************************************************
!*******************************************************************
!*******************************************************************
!*******************************************************************
!*******************************************************************

        subroutine dec2dms(a,d,m,s,sign)

        implicit none

        double precision :: a,d,m,s,sign

        double precision :: aa,aux
        double precision :: sixty
        parameter(sixty=60.)

        aa = abs(a)
        d = int(aa)
        aux = sixty*(aa-d)
        m = int(aux)
        s = sixty*(aux-m)

        sign = one
        if( a < zero ) sign = -sign

        end subroutine dec2dms

!*******************************************************************

        subroutine dms2dec(d,m,s,sign,a)

! transforms deg min sec to decimal degrees

        double precision :: d,m,s,sign,a

        double precision :: sixty,sixty2
        parameter(sixty=60.,sixty2=3600.)

        a = d + m/sixty + s/sixty2
        a = a * sign

        end subroutine dms2dec

!==================================================================
        end module projection
!==================================================================

        subroutine init_coords(iproj,c_param)

! initialization of transformation routines
!
! please note: the coordinates (xtrans,ytrans) are subtracted from the
!	computed cartesian grid and are in meters
!
! iproj:
!
!	1	gauss-boaga
!	2	UTM
!	3	equidistant cylindrical, carte parallelogrammatique (CPP)
!	4	UTM non standard
!
! c_param:
!	meaning of c_param changes depending on projection
!	please see below the exact meaning of the parameters
!
! if GB, zone (fuse) must be 1 or 2
! if UTM, zone must be 1-60
!
! GB: Venice is in fuse 2, central meridian is 15
! UTM: Venice is in zone 33, central meridian is 15

        use projection

        implicit none

        integer :: iproj 		  	  !type of projection
        double precision, dimension(9) :: c_param !parameters for projection

        integer :: fuse,zone
        double precision :: xtrans,ytrans
        double precision :: lon0,lat0,phi
        double precision :: lambda,k

        select case (iproj)
          case ( 0 )                     !no projection
          case ( 1 )                     !Guass-Boaga
            fuse = nint( c_param(1) )     !fuse for gauss-boaga
            xtrans = c_param(2)           !extra shift in x [m]
            ytrans = c_param(3)           !extra shift in y [m]
            call gb_init(fuse)
            call gb_trans(xtrans,ytrans)
          case ( 2 )                     !UTM
            zone = nint( c_param(1) )     !zone for UTM
            xtrans = c_param(2)           !extra shift in x [m]
            ytrans = c_param(3)           !extra shift in y [m]
            call utm_init(zone)
            call utm_trans(xtrans,ytrans)
          case ( 3 )               	!equidistant cylindrical
            phi  = c_param(1)             !central latitude
            lon0 = c_param(2)             !longitude of origin
            lat0 = c_param(3)             !latitude of origin
            call cpp_init(lon0,lat0,phi)
          case ( 4 )              	!UTM (non standard)
            lambda = c_param(1)           !central meridian for UTM
            xtrans = c_param(2)           !extra shift in x [m]
            ytrans = c_param(3)           !extra shift in y [m]
            k = c_param(4)                !scale factor
            call utm_init_nonstd(lambda)
            call utm_set_scale_factor(k)
            call utm_trans(xtrans,ytrans)
          case default
            write(6,*) 'iproj = ',iproj
            stop 'error stop init_coords: unknown projection'
        end select

        proj = iproj

        end subroutine init_coords

!********************************************************************

        subroutine convert_coords(mode,n,xc,yc,xg,yg)

! conversion of coordinates

        use projection

        implicit none

        integer, intent(in)  :: mode		! +1: cart to geo  -1: geo to cart
        integer, intent(in)  :: n		! number of points
        real, intent(inout)  :: xc(n),yc(n)	! cartesian coordinates
        real, intent(inout)  :: xg(n),yg(n)	! geographical coordinates

        write(6,*) 'using projection = ',proj

        select case ( mode )
          case ( 1 )        		!cart to geo

            select case (proj)
              case ( 0 )    		!no projection
                xg = xc
                yg = yc
              case ( 1 )                !Guass-Boaga
                call apply_coords(n,xc,yc,xg,yg,gb_c2g)
              case ( 2 )                !UTM
                call apply_coords(n,xc,yc,xg,yg,utm_c2g)
              case ( 3 )                !CPP
                call apply_coords(n,xc,yc,xg,yg,cpp_c2g)
              case ( 4 )                !UTM non standard
                call apply_coords(n,xc,yc,xg,yg,utm_c2g)
              case default
                stop 'error stop convert_coords: internal error'
            end select

          case ( -1 )              	!geo to cart

            select case (proj)
              case ( 0 )                !no projection
                xc = xg
                yc = yg
              case ( 1 )                !Guass-Boaga
                call apply_coords(n,xg,yg,xc,yc,gb_g2c)
              case ( 2 )                !UTM
                call apply_coords(n,xg,yg,xc,yc,utm_g2c)
              case ( 3 )                !CPP
                call apply_coords(n,xg,yg,xc,yc,cpp_g2c)
              case ( 4 )                !UTM non standard
                call apply_coords(n,xg,yg,xc,yc,utm_g2c)
              case default
                stop 'error stop convert_coords: internal error'
            end select

          case default
            stop 'error stop convert_coords: unknown mode'
        end select

        end subroutine convert_coords

!********************************************************************

        subroutine apply_coords(n,xs,ys,xt,yt,sub)
                
        use projection

        implicit none

        integer, intent(in) :: n
        real, intent(in)    :: xs(n), ys(n)
        real, intent(out)   :: xt(n), yt(n)
        external sub

        integer :: i
        double precision :: xsource,ysource,xtarget,ytarget

        do i=1,n
            xsource = xs(i)
            ysource = ys(i)
            call sub(xsource,ysource,xtarget,ytarget)
            xt(i) = xtarget
            yt(i) = ytarget
        end do

        end subroutine apply_coords

!*******************************************************************
!*******************************************************************
!********** AUXILIARY ROUTINES FOR TESTING *************************
!*******************************************************************
!*******************************************************************
        	
        function proj_ran1(idum)

        implicit none

        real :: proj_ran1
        integer :: idum	! idum < 0 => reset

        integer :: m1,m2,m3,ia1,ia2,ia3,ic1,ic2,ic3
        real :: rm1,rm2
        parameter (m1=259200,ia1=7141,ic1=54773,rm1=3.8580247e-6)
        parameter (m2=134456,ia2=8121,ic2=28411,rm2=7.4373773e-6)
        parameter (m3=243000,ia3=4561,ic3=51349)
        integer :: ndim
        parameter (ndim=97)

        integer :: ix1,ix2,ix3
        real :: r(ndim)
        save ix1,ix2,ix3,r
        integer :: j

        integer :: iff
        save iff
        data iff /0/

        if ( idum < 0 .OR. iff == 0 ) then
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
        if(j > ndim .OR. j < 1) then
            write(6,*) j
            stop 'error stop proj_ran1: internal error'
        end if

        proj_ran1=r(j)
        r(j)=(float(ix1)+float(ix2)*rm2)*rm1

        end function proj_ran1

!*******************************************************************

        function proj_ran_intv(idum,rmin,rmax)

        implicit none

        real :: proj_ran_intv
        integer :: idum
        real :: rmin,rmax

        real :: proj_ran1

        proj_ran_intv = rmin + proj_ran1(idum) * (rmax-rmin)

        end function proj_ran_intv

!********************************************************************
!********************************************************************
!********************************************************************


!********************************************************************

        subroutine genproj_all_test
        end subroutine genproj_all_test

!********************************************************************

!         program genproj_main
!         call genproj_all_test
!         end program genproj_main

!*******************************************************************

        subroutine proj_all_test
        end subroutine proj_all_test

!*******************************************************************

!       program proj_main
!       call proj_all_test
!       end

!*******************************************************************

        subroutine all_test
        write(6,*) '================================================='
        call genproj_all_test
        write(6,*) '================================================='
        call cpp_all_test
        write(6,*) '================================================='
        call gb_all_test
        write(6,*) '================================================='
        call utm_all_test
        write(6,*) '================================================='
        call proj_all_test
        write(6,*) '================================================='
        end subroutine all_test

!********************************************************************
!	program main_test
!	call all_test
!	end program main_test

!********************************************************************

!********************************************************************

        subroutine cpp_all_test
        end subroutine cpp_all_test

!********************************************************************

!	program cpp_main
!	call cpp_all_test
!	end program cpp_main

!**********************************************************************
!**********************************************************************
!**********************************************************************
        	
        subroutine gb_sl_test
        	
! test somma lombardo

! FIXME -> find expected values and insert them
        	
        use projection

        implicit none
        	
        integer :: fuso
        double precision :: deg,min,sec,sign
        double precision :: lon_mm,phi_sl,lamb_sl,lon,lat,x,y

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

        end subroutine gb_sl_test
         
!****************************************

        subroutine deg_test

        use projection

        implicit none

        integer :: i,idum,ndim
        double precision :: dc,dc2,geo2cent,eps
        double precision :: deg,min,sec,sign
        real :: proj_ran1
        	
        write(6,*) 'test degrees...'

        idum = 98765
        eps = 1.e-8
        ndim = 10000
        	
        do i=1,ndim
            dc = 360. * proj_ran1(idum)
            dc = dc - 180.
            call dec2dms(dc,deg,min,sec,sign)
            call dms2dec(deg,min,sec,sign,dc2)
            if( abs(dc-dc2) > eps ) then
                write(6,*) '*** ',i,deg,min,sec,dc,dc2
            end if
!write(6,*) i,deg,min,sec,dc,dc2
        end do

        write(6,*) 'finished testing ',ndim,' loops'
        write(6,*) 'end of test degrees...'

        end subroutine deg_test
         
!****************************************

        subroutine gb_rand_test

        use projection

        implicit none

        integer :: i,idum,nerr,ndim,fuso
        double precision :: eps,res,resmax
        double precision :: x,y,lon,lat,lon1,lat1
        real :: r
        real :: proj_ran1

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
            if( res > eps ) then
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

        end subroutine gb_rand_test
         
!****************************************

        subroutine gb_all_test
        call gb_sl_test
        call deg_test
        call gb_rand_test
        end subroutine gb_all_test

!****************************************

!	 program gb_main
!	 call gb_all_test
!	 end program gb_main

!**********************************************************************
!**********************************************************************

        subroutine utm_sl_test
        	
! test somma lombardo and other

! FIXME -> find expected values and insert them
        	
        use projection

        implicit none
        	
        integer :: fuso
        double precision :: deg,min,sec,sign
        double precision :: phi_sl,lamb_sl,lon,lat,x,y
        double precision :: xa,xe,ye

        write(6,*) 'test somma lombardo...'

        call utm_init(32)

        call dms2dec(45.d0,41.d0,0.d0,1.d0,phi_sl)
        call dms2dec(8.d0,42.d0,0.d0,1.d0,lamb_sl)

        lon = lamb_sl
        lat = phi_sl
        	
        call utm_g2c(lon,lat,x,y)
        write(6,*) 'somma lombardo:'
        write(6,*) 'lat/lon: ',lon,lat
        write(6,*) 'raw utm: ',x,y
        write(6,*) 'adjusted: ',x+500000.,y
        write(6,*) 'expected: ',476638.,5058908.
        call utm_c2g(x,y,lon,lat)
        write(6,*) 'inverse: ',lon,lat

        lat = 47.3782
        lon = 8.2325

        call utm_g2c(lon,lat,x,y)
        write(6,*) 'other point:'
        write(6,*) 'lat/lon: ',lon,lat
        write(6,*) 'raw utm: ',x,y
        write(6,*) 'adjusted: ',x+500000.,y
        write(6,*) 'expected: ',442063.49206,5247475.33300
        call utm_c2g(x,y,lon,lat)
        write(6,*) 'inverse: ',lon,lat

        write(6,*) 'end of test somma lombardo...'

        end subroutine utm_sl_test
         
!****************************************

        subroutine utm_rand_test

        use projection

        implicit none

        integer :: i,idum,nerr,ndim,zone
        double precision :: eps,res,resmax
        double precision :: x,y,lon,lat,lon1,lat1
        real :: r
        real :: proj_ran1

        write(6,*) 'test random utm...'

        idum = 98765
        eps = 1.e-6
        resmax = 0.
        nerr = 0
        ndim = 100000

        call utm_init(1)	!initialize with any zone

        do i=1,ndim
            lon = 360.*proj_ran1(idum) - 180.
            lat = 180.*proj_ran1(idum) - 90

            call utm_compute_zone(lon,zone)
            call utm_g2c(lon,lat,x,y)
            call utm_c2g(x,y,lon1,lat1)

            res = sqrt( (lat-lat1)**2 + (lon-lon1)**2 )
            resmax = max(resmax,res)
            if( res > eps ) then
                write(6,*) '*** ',i,zone,lon,lon1,lat,lat1,res
                write(6,*) '***   ',x,y
                nerr = nerr + 1
            end if
!write(6,*) i,zone,lon,lon1,lat,lat1,res
!write(6,*) '      ',x,y
        end do

        write(6,*) 'finished testing ',ndim,' loops '
        write(6,*) 'number of errors:   ',nerr
        write(6,*) 'maximum difference: ',resmax
        write(6,*) 'end of test random utm...'

        end subroutine utm_rand_test
         
!****************************************

        subroutine utm_band_test

        use projection

        implicit none

        double precision :: lat
        character(1) :: band

        lat = -79.5

        do while( lat <= 84. )
            call utm_compute_band(lat,band)
            write(6,*) lat,'  ',band
            lat = lat + 1.
        end do

        end subroutine utm_band_test

!****************************************

        subroutine utm_all_test
!call utm_band_test
        call utm_sl_test
        call utm_rand_test
        end subroutine utm_all_test

!****************************************

!	 program utm_main
!	 call utm_all_test
!	 end program utm_main

!****************************************

!*******************************************************************
!*******************************************************************
!*******************************************************************
!*******************************************************************
!*******************************************************************
!
!        PIN->e = sqrt(PIN->es);
!        PIN->ra = 1. / PIN->a;
!        PIN->one_es = 1. - PIN->es;
!        if (PIN->one_es == 0.) { pj_errno = -6; goto bum_call; }
!        PIN->rone_es = 1./PIN->one_es;
!
! pj_tsfn(double phi, double sinphi, double e) {
!         sinphi *= e;
!         return (tan (.5 * (HALFPI - phi)) /
!            pow((1. - sinphi) / (1. + sinphi), .5 * e));
! }
!
! pj_msfn(double sinphi, double cosphi, double es) {
!         return (cosphi / sqrt (1. - es * sinphi * sinphi));
! }
!
! #define TOL 1.0e-10
! #define N_ITER 15
!
!         double
! pj_phi2(double ts, double e) {
!         double eccnth, Phi, con, dphi;
!         int i;
!
!         eccnth = .5 * e;
!         Phi = HALFPI - 2. * atan (ts);
!         i = N_ITER;
!         do {
!                 con = e * sin (Phi);
!                 dphi = HALFPI - 2. * atan (ts * pow((1. - con) /
!                    (1. + con), eccnth)) - Phi;
!                 Phi += dphi;
!         } while ( fabs(dphi) > TOL && --i);
!         if (i <= 0)
!                 pj_errno = -18;
!         return Phi;
! }
!
!==================================================================

