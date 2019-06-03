
!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2018  Georg Umgiesser
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

c heat flux module (utility 3, radiation)
c
c contents :
c
c revision log :
c
c 09.12.2002	ggu	new import into file
c 09.04.2003	ggu	new routines for day light
c 17.08.2004	ggu	bises and jdmon deleted (now in newdat.f)
c 22.09.2004	ggu	bug fix for idmon/jdmon
c 08.10.2008	ggu	jdmon0 and bises0 re-introduced
c 23.03.2010	ggu	changed v6.1.1
c 16.02.2011	ggu	decl() commented
c 16.02.2019	ggu	changed VERS_7_5_60
c
c*****************************************************************************

	subroutine qsun(year,month,day,hour,cc,qs)

c computes solar radiation incident to sea water
c
c uses only date, time and cloud cover
c
c Gill, p.34

	implicit none

	integer year	!year (fully qualified, e.g., 1995)
	integer month	!month [1,12]
	integer day	!day [1,28-31]
	real hour	!hour [0-24)		!ATTENTION!! -> REAL
	real cc		!cload cover [fractional] (0-1)
	real qs		!solar radiation into water [W/m**2]

	real qtopat	!solar radiation at top of atmosphere [W/m**2]
	real albedo	!mean albedo [fractional]
	real absorb	!average absorption of atmosphere [fractional]

	real ad		!declination of sun [rad]
	real rr		!?? (correction??)
	real ai		!angle of incidence (with respect to zenith) [rad]
	real cosi	!cos(ai)

	albedo = 0.06
	absorb = 0.25

	call decl(year,month,day,ad,rr)

	call incid(hour,ad,ai,cosi)
	call radmax(rr,cosi,qtopat)

	qs = qtopat * (1.-absorb) * (1.-albedo) * (1.-0.7*cc)

	qs = max(qs,0.)

c	write(6,*) month,day,hour,qtopat,qs

	end

c**************************************************************************

      subroutine srad(ia,im,id,ah,qrt,qrdif,qstot,cosi,qsmax,cc)

c  input:   ia, im, id, ah, qrt, qrdif
c  output:  qstot, qsmax, cc
c
c    ia,im,id,ah = year, month, day, hour
c    qrt = total solar radiation flux [W/m**2]
c    qrdif = diffuse solar radiation flux [W/m**2]
c    qstot = net solar radiation flux (reflection already subtracted) [W/m**2]
c    qsmax = maximum solar radiation flux, theoretically calculated [W/m**2]
c    cc = cloud cover
c
c    cosi = cos angle of incidence
c
c cloud cover cc: it is computed only during the day, when cosi>0
c                 during the night its value is the same as the last 
c                 hour (with cosi > 0) of that day.

	real ccold
	save ccold
	data ccold /1./

        qidir = qrt - qrdif

        call decl(ia,im,id,ad,rr)
        call incid(ah,ad,ai,cosi)

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

c**************************************************************************

      subroutine decl(iy,im,id,ad,rr)

c computes declination and relative distance from sun

	implicit none

	integer iy		!year
	integer im		!month
	integer id		!day
	real ad			!declination [rad] (return)
	real rr			!relative distance (return)

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

c**************************************************************************

      subroutine incid(ah,ad,ai,cosi)

c input:
c  ah = hour
c  ad = declinazione solare
c output:
c  ai = angolo di incidenza della luce solare (in rad rispetto allo zenith)
c  cosi = cos(ai)

      pi=4*atan(1.)
      omega=2*pi/24.
      alat=45.5
      alatr=alat*pi/180.
      t12=ah-12.

      cosi=sin(ad)*sin(alatr)+cos(ad)*cos(alatr)*cos(omega*t12)
      ai=acos(cosi)
     
      end

c**************************************************************************

      subroutine rifl(ai,ro)

c input:  ai = angolo di incidenza della luce solare [rad]
c output: ro = percentage of reflected radiation [0-1]

      an2=1.33

      aj = asin(sin(ai)/an2)
      ap = aj + ai
      am = aj - ai

      ro = 0.5*( (sin(am)/sin(ap))**2 + (tan(am)/tan(ap))**2 ) 

      end

c**************************************************************************

      subroutine totq(qidir,qidif,ro,qitot)

      rodif=0.0589
      qitot = (1-ro)*qidir + (1-rodif)*qidif

      end

c**************************************************************************

         subroutine radmax(rr,cosi,qsmax)

c input:
c   rr = ??
c   cosi = cos(ai)
c output:
c   qsmax = theoretical maximum radiation [W/m**2]
c
c  Lazzarin R., 1981: Sistemi solari attivi
c  sc = costante solare [W/m**2]

      sc = 1353.
      qsmax = rr*sc*cosi

      end

c**************************************************************************

         subroutine cloud(qstot,qsmax,cc)

         cctmp = 1. - qstot/qsmax
         if (cctmp.ge.0.) cc=cctmp

      end

c**************************************************************************
c**************************************************************************
c**************************************************************************
c jdmon0 and bises0 introduced to be independent from other files
c**************************************************************************
c**************************************************************************
c**************************************************************************

        function jdmon0(year,month)

c returns total days from January to end of month

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

c**************************************************************************

        function bises0(year)

c true if year is bisestial

        logical bises0
        integer year

        if(  ( mod(year,4) .eq. 0 .and. mod(year,100) .ne. 0 )
     +                  .or. mod(year,400) .eq. 0 ) then
                bises0 = .true.
        else
                bises0 = .false.
        end if

        end

c**************************************************************************
c**************************************************************************
c**************************************************************************
c**************************************************************************
c**************************************************************************

	subroutine daylight(lat,jday,rnh)

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

c**************************************************************************

	subroutine test_daylight

	implicit none

	integer ndim
	parameter (ndim=9)

	real lat,rnh
	integer j,l

	real rl(ndim),rn(ndim)
c	data rl / 0.,30.,45.,60.,89. /
c	data rl / 60.,65.,70.,75.,80. /
	data rl / 0.,30.,45.,60.,70.,75.,80.,85.,90. /

	lat = 45.

	do j=1,365
	  do l=1,ndim
	    lat = rl(l)
	    call daylight(lat,j,rnh)
	    rn(l) = rnh
	  end do
	  write(6,1000) j,(rn(l),l=1,ndim)
	end do

 1000	format(i4,10f6.2)
	end

c**************************************************************************

c	call test_daylight
c	end

c**************************************************************************

