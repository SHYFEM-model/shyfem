c
c $Id: tideg.f,v 1.7 2010-02-22 15:38:36 georg Exp $
c
c revision log :
c
c 08.04.1999	ggu	written from scratch
c 13.04.1999	ggu	itide introduced
c 19.04.1999	ggu	itide changed to rtide
c 27.08.2007    ccf     isphe for spherical coordinate system
c 11.03.2009    ggu     not adjusted for new x/ygeov -> still to be handled
c 05.02.2010    ccf     include new tideforc routine for tidal potential
c 18.11.2011    ggu     excluded projection code from tide
c 12.12.2011    ccf     new component (MSm) in tidal computation
c 14.01.2015    ccf     small bug fix rounding date and time
c
c********************************************************************

	subroutine tideini

	use mod_tides
	use basin

	implicit none

	include 'param.h'
	include 'pkonst.h'


	real rtide
	common /tidcom/ rtide
	save /tidcom/

	real xmin,xmax,ymin,ymax
	integer iproj
	integer k

	double precision dgetpar
        integer isphe          !if = 1  coordinates are in spherical system
        integer date,time

	rtide = dgetpar('rtide')

	if( rtide .le. 0. ) return

	iproj = nint(dgetpar('iproj'))
	call get_coords_ev(isphe)

	if( isphe .le. 0 .and. iproj .eq. 0 ) then
	  write(6,*) 'isphe,iproj: ',isphe,iproj
	  write(6,*) 'for tidal potential either '
	  write(6,*) 'isphe or iproj must be set'
	  stop 'error stop tideini: no lat/lon coordinates'
	end if

	do k=1,nkn
	  zeqv(k) = 0.
	end do

	date = nint(dgetpar('date'))
	time = nint(dgetpar('time'))
	call dtsini(date,time)

	write(6,*) 'Tidal potential active'

	end

c********************************************************************

	subroutine tidenew(it)

	use mod_tides
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer it

	include 'param.h'
	include 'pkonst.h'


	real rtide
	common /tidcom/ rtide

	integer iy,im,id,ih,imn,isec
	integer k
	real efact		!reduction due to earth deformation
	parameter( efact=0.693 )


	if( rtide .le. 0 ) return
	
	iy = 1990	!year
	im = 1		!month
	id = 1		!day
	iy = 2003	!year
	im = 8		!month
	id = 1		!day
	iy = 2004	!year
	im = 1		!month
	id = 1		!day


	ih = 0		!hour
	imn = 0		!minute

	call dts2dt(it,iy,im,id,ih,imn,isec)

	if(iy .lt. 50)then
          iy = iy + 2000
	else if(iy .lt. 100)then
          iy = iy + 1900
	endif

	call tidesc(iy,im,id,ih,imn,it)

	do k=1,nkn
	  call tidelc(xgeov(k),ygeov(k),zeqv(k))
	  zeqv(k) = rtide * zeqv(k)*efact
	end do

	end

!********************************************************************
!*     Astronomical tidal force as body + earth + load               *
!*********************************************************************
! INPUT ARGUMENT LIST:
!     time    - time from beginning of simulation (days)
!     iday    - day of the month
!     imonth  - month of the year
!     iyear   - year (4-digit)
!     nkn     - total node number
!
! OUTPUT ARGUMENT LIST:
!     eeq     - array with tidal forcing
!
! The model cosiders these tidal costituents:
!  - Semidiurnal Species:
!    - M2  principal lunar
!    - S2  principal solar
!    - N2  elliptical lunar
!    - K2  declination luni-solar
!  - Diurnal Species:
!    - K1  declination luni-solar
!    - O1  principal lunar
!    - P1  principal solar
!    - Q1  elliptical lunar
! - Long-Period Species:
!    - Mf  fortnightly lunar
!    - Mm  monthly lunar
!    - Ssa semiannual solar
!    - MSm S0-semiannual solar
!
!*********************************************************************

	subroutine tideforc(it)

	use mod_tides
	use mod_depth
	use mod_hydro
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer it

! common
	include 'param.h'
	include 'pkonst.h'


	real rtide
	common /tidcom/ rtide

! local
	real lat,lon,eqt
	integer iy,im,id,ih,imn,isec
	real time
	integer k,jd
	real hour,tis
	integer ntd		!number of tidal costituents
	parameter (ntd=12)
	real chi(ntd)		!astronomical arguments [rad]
	real loadb		!loading tide factor [0.054,0.047,= lbd*depth]
	real lbd		!coefficient for load
	parameter (lbd = 7e-5)

!--------------------------------------------------------
!------ compute tidal potential? ------------------------
!--------------------------------------------------------

	if( rtide .le. 0 ) return
	
!--------------------------------------------------------
!------ computes julian day - 1 -------------------------
!--------------------------------------------------------

	call dts2dt(it,iy,im,id,ih,imn,isec)

	if(iy .lt. 50)then
          iy = iy + 2000
	else if(iy .lt. 100)then
          iy = iy + 1900
	endif

	time = it / 86400.
	hour = time - int(time)
        call date2j(iy,im,id,jd)
	tis = jd - 1.

!--------------------------------------------------------
!------ computes astronomical arguments -----------------
!--------------------------------------------------------

	call get_chi(tis,iy,chi)

!--------------------------------------------------------
!------ computes eq. tide and account for load tide -----
!--------------------------------------------------------

	do k = 1,nkn
	  lat = ygeov(k)
	  lon = xgeov(k)
	  call equilibrium_tide(lat,lon,hour,chi,eqt)
	  loadb = lbd * (hkv(k) + zov(k))
	  zeqv(k) = eqt + loadb*zov(k) 
	end do

!--------------------------------------------------------
!------ end of routine ----------------------------------
!--------------------------------------------------------

	end

!********************************************************************
! This subroutine computes astronomical arguments chi (in degrees)

        subroutine get_chi(tis,iyear,chi)

        implicit none

        include 'param.h'
  
        real tis		!julian day - 1
        integer iyear		!year
        integer ntd		!number of tidal costituents
        parameter (ntd=12)
        real chi(ntd)		!astronomical arguments [rad]

        integer i
        real degrad
        parameter( degrad = 57.2957795131 ) ! 360/(2*pi)

        real h01,h02,h03          !  h0 = h01 + h02*T + h03*t**2
        parameter(h01 = 279.69668, h02 = 36000.768930485, h03=3.03E-04)
        real s01,s02,s03,s04      !  s0 = s01 + s02*t + s03*t**2 + s04*t**3
        parameter(s01 = 270.434358, s02 = 481267.88314137, 
     &		  s03 = -0.001133,  s04 = 1.9E-06 )
        real p01,p02,p03,p04      !  p0 = p01 + p02*t + p03*t**2 + p04*t**3
        parameter(p01 = 334.329653, p02 = 4069.0340329575, 
     &		  p03 = -0.010325,  p04 = -1.2E-05 )
        real t1,t2                !  t = t1 + t2*d
        parameter( t1 = 0.74996579132, t2 = 2.7378508846E-05 )

        real day0,day,d,t,h0,s0,p0

        day0  = 1.0

        day = day0 + tis
        d  = day + 365.0*(iyear - 1975.0) + int((iyear-1973.0)/4.0)
        t  = t1 + t2*d
        h0 = h01 + h02*t + h03*(t**2.0)
        s0 = s01 + s02*t + s03*(t**2.0) + s04*(t**3.0)
        p0 = p01 + p02*t + p03*(t**2.0) + p04*(t**3.0)

        chi(1) = 2*h0 - 2*s0           !M2
        chi(2) = 0.0                   !S2
        chi(3) = 2*h0 - 3*s0 + p0      !N2
        chi(4) = 2*h0                  !K2

        chi(5) = h0        + 90.0      !K1
        chi(6) = h0 - 2*s0 - 90.0      !O1
        chi(7) = -h0       - 90.0      !P1
        chi(8) = h0 - 3*s0 + p0 -90.0  !Q1

        chi(9) =      2*s0             !Mf
        chi(10) =       s0 - p0        !Mm
        chi(11) = 2*h0                 !Ssa
        chi(12) =-2*h0 + s0 +p0	       !MSm

! ----- Convert to Radians
        do i=1,ntd
          chi(i) = chi(i)/degrad
        enddo

        end

!*********************************************************************
! This subroutine computes the equilibrium tide as a sum of the single
! costituents

        subroutine equilibrium_tide(lat,lon,hour,chi,eqtide)

	implicit none

	integer ntd			!Number of tidal costituents
        parameter (ntd=12)
        real lat,lon			!Latitude and longitude [degree]
	real hour			!Universal standard time [days]
        real chi(ntd)			!Astronomic arguments
        real eqtide			!Equilibrium tide [m]

	real ome(ntd)			!Tidal constituent frequency [2Pi/s]
	real amp(ntd)			!Tidal constituent amplitude [m]
	real loven(ntd)			!Elasticy factor 
	real mask_t(ntd)		!Mask: =1, tide included

        data ome/1.40519e-4,1.45444e-4,1.37880e-4,1.458426e-4,
     &           0.72921e-4,0.67598e-4,0.72523e-4,0.649584e-4,
     &		 0.053234e-04,0.026392e-04,0.003982e-04,0.049252e-04/
	data amp/0.242334,0.112841,0.046398,0.030704,
     &           0.141465,0.100514,0.046843,0.019256,
     &		 0.041742,0.022026,0.019446,0.004239/
        data loven/0.6930,0.6930,0.6930,0.6930,
     &             0.7364,0.6950,0.7059,0.6946,  
     &             0.6920,0.6920,0.6920,0.6920/  
	data mask_t/1., 1., 1., 1., 
     &		    1., 1., 1., 1.,
     &		    1., 1., 1., 1./

	integer n,m,l
        real pi2,rad
        parameter(pi2=6.28318530717796d0,rad=pi2/360.d0)

        real colat,twolon,time
        real sins,sind,sinl

        eqtide = 0.
        time   = hour * 86400.
	if (lon .lt. 0.0) lon = 360. + lon
        lon    = rad * lon
        twolon = 2.0*lon
        colat  = rad * (90. - lat)
        sins   = sin(colat)**2.0
	sind   = sin(2.0*colat)
	sinl   = 3.0*sins - 2.0

        do n = 1,4
! ------- Semidiurnal species -------
          eqtide = eqtide + mask_t(n) * loven(n) * 
     & 		   amp(n) * sins *cos(ome(n)*time + chi(n) + twolon)
! ------- Diurnal species -------
          m = n+4
          eqtide = eqtide + mask_t(m) * loven(m) * 
     & 		   amp(m) * sind *cos(ome(m)*time + chi(m) + lon)
! ------- Long-period species -------
          l = m+4
          eqtide = eqtide + mask_t(l) * loven(l) * 
     & 		   amp(l) * sinl *cos(ome(l)*time + chi(l))
        enddo

	return
	end

!*********************************************************************
