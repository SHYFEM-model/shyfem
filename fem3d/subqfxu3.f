
!--------------------------------------------------------------------------
!
!    Copyright (C) 2002-2004,2008,2010-2011,2019  Georg Umgiesser
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

! heat flux module (utility 3, radiation)
!
! contents :
!
! revision log :
!
! 09.12.2002	ggu	new import into file
! 09.04.2003	ggu	new routines for day light
! 17.08.2004	ggu	bises and jdmon deleted (now in newdat.f)
! 22.09.2004	ggu	bug fix for idmon/jdmon
! 08.10.2008	ggu	jdmon0 and bises0 re-introduced
! 23.03.2010	ggu	changed v6.1.1
! 16.02.2011	ggu	decl() commented
! 16.02.2019	ggu	changed VERS_7_5_60
! 08.03.2023	ggu	commented and added test routines
!
!*****************************************************************************

	subroutine qsun(year,month,day,hour,alat,cc,qs)

! computes solar radiation incident to sea water
!
! uses only date, time, latitude, and cloud cover
!
! Gill, p.34

	implicit none

	integer year	!year (fully qualified, e.g., 1995)
	integer month	!month [1,12]
	integer day	!day [1,28-31]
	real hour	!hour [0-24)		!ATTENTION!! -> REAL
	real alat	!latitude [deg]
	real cc		!cload cover [fractional] (0-1)
	real qs		!solar radiation into water [W/m**2] (return)

	real qtopat	!solar radiation at top of atmosphere [W/m**2]
	real albedo	!mean albedo [fractional]
	real absorb	!average absorption of atmosphere [fractional]

	real ad		!declination of sun [rad]
	real rr		!relative distance from sun [fraction]
	real ai		!angle of incidence (with respect to zenith) [rad]
	real cosi	!cos(ai)

	albedo = 0.06
	absorb = 0.25

	call decl(year,month,day,ad,rr)

	call incid(hour,alat,ad,ai,cosi)
	call radmax(rr,cosi,qtopat)

	qs = qtopat * (1.-absorb) * (1.-albedo) * (1.-0.7*cc)

	qs = max(qs,0.)

!	write(6,*) month,day,hour,qtopat,qs

	end

!**************************************************************************

      subroutine srad(ia,im,id,ah,qrt,qrdif,qstot,cosi,qsmax,cc)

! still to comment and check - not used

! input:   ia, im, id, ah, qrt, qrdif
! output:  qstot, qsmax, cc
!
!    ia,im,id,ah = year, month, day, hour
!    qrt = total solar radiation flux [W/m**2]
!    qrdif = diffuse solar radiation flux [W/m**2]
!    qstot = net solar radiation flux (reflection already subtracted) [W/m**2]
!    qsmax = maximum solar radiation flux, theoretically calculated [W/m**2]
!    cc = cloud cover
!
!    cosi = cos angle of incidence
!
! cloud cover cc: it is computed only during the day, when cosi>0
!                 during the night its value is the same as the last 
!                 hour (with cosi > 0) of that day.

	real, save :: ccold = 1.
	real, save :: alat = 45.5	!average latitude

	stop 'do not use this routine: srad'

        qidir = qrt - qrdif

        call decl(ia,im,id,ad,rr)
        call incid(ah,alat,ad,ai,cosi)

        if(cosi.ge.0.) then			!day
          call rifl(ai,ro)
          call totq(qidir,qrdif,ro,qstot)
          call radmax(rr,cosi,qsmax)
          call cloud(qstot,qsmax,cc)
	  ccold = cc
        else					!night
          qstot=0.
          qsmax=0.
	  cc = ccold
        endif

        end

!**************************************************************************

      subroutine decl(iy,im,id,ad,rr)

! computes declination and relative distance from sun

	implicit none

	integer iy		!year
	integer im		!month
	integer id		!day
	real ad			!declination [rad] (return)
	real rr			!relative distance from sun (return)

	integer nd
	real pi
	real b

        integer jdmon0		! bug fix 22.09.2004

        pi=4*atan(1.)

        nd=id+jdmon0(iy,im-1)	! Julian day (1 Jan => nd = 1)
	b = 2.*pi*(284.+nd)/365.

        ad=(23.45*pi/180.)*sin(2.*pi*(284.+nd)/365.)
        rr = 1. + 0.033*cos(b)

        end

!**************************************************************************

        subroutine incid(ah,alat,ad,ai,cosi)

! computes angle of incidence

        implicit none

	real ah		!hour
	real alat	!latitude [deg]
	real ad		!solar declination
	real ai		!angle of incidence respect to zenith [rad] (return)
	real cosi	!cos(ai) (return)
      
	real, save :: pi=4*atan(1.)
	real omega,t12,alatr

        pi=4*atan(1.)
        omega=2*pi/24.
        alatr=alat*pi/180.
        t12=ah-12.

        cosi=sin(ad)*sin(alatr)+cos(ad)*cos(alatr)*cos(omega*t12)
        ai=acos(cosi)
     
        end

!**************************************************************************

        subroutine rifl(ai,ro)

! computes fraction of reflected radiation

	implicit none

	real ai		!angle of incidence of solar radiation [rad]
	real ro		!fraction of reflected radiation [0-1]

	real an2,aj,ap,am

        an2=1.33

        aj = asin(sin(ai)/an2)
        ap = aj + ai
        am = aj - ai

        ro = 0.5*( (sin(am)/sin(ap))**2 + (tan(am)/tan(ap))**2 ) 

        end

!**************************************************************************

        subroutine totq(qidir,qidif,ro,qitot)

! computes total ratiation from direct and diffuse radiation

	implicit none

	real qidir	!direct solar radiation [W/m**2]
	real qidif	!diffuse solar radiation [W/m**2]
	real ro		!fraction of reflected radiation [0-1]
	real qitot	!total solar radiation [W/m**2] (return)

	real, parameter ::  rodif = 0.0589

        qitot = (1-ro)*qidir + (1-rodif)*qidif

        end

!**************************************************************************

        subroutine radmax(rr,cosi,qsmax)

! computes maximum solar radiation
!
! Lazzarin R., 1981: Sistemi solari attivi

	implicit none

	real rr		!relative distance from sun [fraction]
	real cosi	!cos(ai) with ai angle of incidence
	real qsmax	!maxium solar radiation [ W/m**2] (return)

	real, parameter :: sc = 1353.	!solar constant [W/m**2]

        qsmax = rr*sc*cosi

	end

!**************************************************************************

         subroutine cloud(qstot,qsmax,cc)

! computes cloud cover from solar radiation and theoretical maximum
!
! if negative leaves old value

	implicit none

	real qstot	!total solar radiation [W/m**2]
	real qsmax	!theoretical maximum solar radiation [W/m**2]
	real cc		!cloud cover [0-1] (return)

	real cctmp

        cctmp = 1. - qstot/qsmax
        if (cctmp.ge.0.) cc=cctmp

	end

!**************************************************************************
!**************************************************************************
!**************************************************************************
! jdmon0 and bises0 introduced to be independent from other files
!**************************************************************************
!**************************************************************************
!**************************************************************************

        function jdmon0(year,month)

! returns total days from January to end of month

        implicit none

        integer jdmon0
        integer year,month

        logical bises0

        integer jmm(0:12)
        save jmm
        data jmm /0,31,59,90,120,151,181,212,243,273,304,334,365/

        jdmon0=jmm(month)

        if( month .ge. 2 .and. bises0(year) ) jdmon0 = jdmon0 + 1

        end

!**************************************************************************

        function bises0(year)

! true if year is bisestial

	implicit none

        logical bises0
        integer year

        if(  ( mod(year,4) .eq. 0 .and. mod(year,100) .ne. 0 )
     +                  .or. mod(year,400) .eq. 0 ) then
                bises0 = .true.
        else
                bises0 = .false.
        end if

        end

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

	subroutine daylight(lat,jday,rnh)

! computes length of daylight in hours

	implicit none

	real lat	!latitude [deg]
	integer jday	!julian day [1-365]
	real rnh	!daylight hours [h] (return)

	real pi,rad,phi,delta,omegas
	real aux,lataux

	pi = 4.*atan(1.)
	rad = pi/180.

	lataux = lat
	if( lataux .ge. 90. ) lataux = 89.999
	if( lataux .le. -90. ) lataux = -89.999

	phi = lataux * rad
	delta = 0.4093 * sin( (2.*pi*jday)/365 - 1.39 )
	aux = -tan(phi) * tan(delta)

	if( aux .ge. 1. ) then
	  rnh = 0.
	else if( aux .le. -1 ) then
	  rnh = 24.
	else
	  omegas = acos( aux )
	  rnh = (24./pi) * omegas
	end if

	if( lat .eq. 90. ) then
	  !write(6,*) jday,delta,-tan(phi),aux,rnh
	end if

	end

!**************************************************************************
!**************************************************************************
!**************************************************************************

	subroutine test_solar_radiation

! computes theoretical solar radiation

	implicit none

	integer year,month,days,day,hour,secs
	integer dayslast,daysend
	real rhour,alat,cc,qs

	integer jdmon0

	alat = 45.5		!latitude

	cc = 0.
	dayslast = 0
	secs = 0

	write(6,*) 'creating solar radiation for latitude',alat
	write(66,*) '#   seconds   solar_rad'

	do year=2012,2012
	  do month=1,12
            days = jdmon0(year,month)
	    daysend = days - dayslast
	    !write(6,*) 'month,days: ',month,days,daysend
	    dayslast = days
	    do day=1,daysend
	      do hour=1,24
	        rhour = hour
		secs = secs + 3600
	        call qsun(year,month,day,rhour,alat,cc,qs)
		write(66,'(i12,f12.2)') secs,qs
	      end do
	    end do
	  end do
	end do


	end

!**************************************************************************

	subroutine test_daylight

	implicit none

	integer ndim
	parameter (ndim=9)

	real alat,rnh
	integer j,l

	real rl(ndim),rn(ndim)
	data rl / 0.,30.,45.,60.,70.,75.,80.,85.,90. /

	write(6,*) 'creating day length for different latitudes'
	write(67,'(a,10f7.2)') '#  day',rl

	do j=1,365
	  do l=1,ndim
	    alat = rl(l)
	    call daylight(alat,j,rnh)
	    rn(l) = rnh
	  end do
	  write(67,1000) j,(rn(l),l=1,ndim)
	end do

	return
 1000	format(i6,10f7.2)
	end

!**************************************************************************
!	program test_main_radiation
!	call test_solar_radiation
!	call test_daylight
!	end
!**************************************************************************

