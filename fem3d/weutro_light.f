
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

c support routines for solar radiation
c
c contents :
c
c revision log :
c
c 20.08.2004    ggu	written from scratch
c 26.08.2004    ggu	new linstant (same as rintens from weutro.f)
c
c notes :
c
c***********************************************************

	function declin(jd)

c computes solar declination

	implicit none

	real declin	!solar declination [degrees]
	integer jd	!julian day [1-365]

	real pi
	parameter(pi=3.14159)

	declin = 23.45 * cos( 2. * pi * (jd-172) / 365. )

	end

c***********************************************************

	function dayl(jd,rlat)

c computes day length (in hours)

	implicit none

	real dayl	!daylength [hours]
	integer jd	!julian day [1-365]
	real rlat	!latitude [degrees]

	real pi,rad
	parameter(pi=3.14159,rad=pi/180.)

	real dec,aux
	real declin

	dec = declin(jd)
	aux = tan(rlat*rad) * tan(dec*rad)
	dayl = (24./pi) * acos( -aux )

	end

c***********************************************************

	function dayl_max(rlat)

c computes maximum day length during year (in hours)

	implicit none

	real dayl_max	!maximum daylength [hours]
	real rlat	!latitude [degrees]

	real pi,rad
	parameter(pi=3.14159,rad=pi/180.)

	real dec,aux

	dec = 23.45
	aux = tan(rlat*rad) * tan(dec*rad)
	dayl_max = (24./pi) * acos( -aux )

	end

c***********************************************************

	function cd_solrad(jd,rlat)

c computes average clear day solar radiation

	implicit none

	real cd_solrad	!average clear day solar radiation [ly/day]
	integer jd	!julian day [1-365]
	real rlat	!latitude [degrees]

	real pi,rad,factw
	parameter(pi=3.14159,rad=pi/180.,factw=0.4843)

	real ld,sinday,cosphi,sinphi
	real aux,a,b
	real dayl_max

	ld = dayl_max(rlat)
	sinday = sin(pi*ld/24.)
	cosphi = cos(rlat*rad)
	sinphi = sin(rlat*rad)

	aux = 0.29 * cosphi + 0.52
	a = sinphi*(46.355*ld-574.3885) + 816.41*cosphi*sinday
	b = sinphi*(574.3885-1.509*ld) - 26.59*cosphi*sinday

	a = a * aux
	b = b * aux

	cd_solrad = a + b * cos( 2.*pi*(jd-172)/365. )

c	cd_solrad = cd_solrad * factw	!radiation in W/m**2

	end

c***********************************************************

	subroutine get_radiation(jd,rlat,fday,itot,imax)

c computes parameters for solar radiation
c
c itot is average solar radiation over one whole day

	implicit none

	integer jd	!julian day [1-365]
	real rlat	!latitude [degrees]
	real fday	!fractional day length [0-1]
	real itot	!average clear day solar radiation [ly/day]
	real imax	!maximum solar radiation at noon [ly/day]

	real pi
	parameter(pi=3.14159)

	real dayl, cd_solrad

	fday = dayl(jd,rlat)/24.
	itot = cd_solrad(jd,rlat)	!itot
	imax = itot*pi/(2.*fday)	!max at noon

	end

c***********************************************************

        subroutine linstant(t,itot,fday,iinst)

c computes light intensity during a day given itot and fday

        implicit none

        real t          !day [0-365]
        real itot       !average incident light intensity over one whole day
        real fday       !fraction of day length [0-1]
        real iinst      !instantaneous light intensity at time t (return)

        real pi
        parameter( pi = 3.14159 )

        real tday,aux

        tday = t - int(t)                       !fraction in day
        if( tday .lt. 0. ) tday = tday + 1.     !negative days
        tday = tday - 0.5                       !maximum at noon

        aux = pi / fday

        if( abs(tday) .le. fday/2 ) then
          iinst = 0.5 * itot * aux * cos( tday * aux )
        else
          iinst = 0.
        end if

        end

c***********************************************************

	subroutine test_solrad

c tests solar radiation routines

	implicit none

	integer jd
	real phi
	real cd,day,dec,daym
	real pi

	real dayl,declin,dayl_max,cd_solrad,cdmax,fday

	pi = 4. * atan(1.)
	phi = 45.
	phi = 50.

	do jd=1,365
	  cd = cd_solrad(jd,phi)	!itot
	  day = dayl(jd,phi)
	  dec = declin(jd)
	  daym = dayl_max(phi)
	  fday = day/24.
	  cdmax = cd*pi/(2.*fday)	!max at noon
	  write(6,*) jd,day,fday,daym,dec,cd,cdmax
	end do

	end

c***********************************************************

	subroutine luxlen_init(file)

c initializes luxlen (reads file)

	implicit none

	character*(*) file

	real luxv(365)
	real fdayv(365)
	common /luxlen_c/ luxv,fdayv
	save /luxlen_c/

	integer n
	real tdummy

c--------------------------------------------------------------
c read in file
c--------------------------------------------------------------

        open(2,file='lux.dat',status='old',form='formatted')
	do n=1,365
          read(2,*) tdummy,luxv(n),fdayv(n)
	end do
	close(2)

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	end

c***********************************************************

	subroutine luxlen(t,lux,fday)

c given a day t returns light intensity lux and day length fday

	implicit none

	real t		!day [0-365]
	real lux	!light intensity [ly/day]
	real fday	!fractional day length [0-1]

	real luxv(365)
	real fdayv(365)
	common /luxlen_c/ luxv,fdayv
	save /luxlen_c/

	integer n
	integer it2n

c--------------------------------------------------------------
c compute values
c--------------------------------------------------------------

	n = it2n(t)

	lux = luxv(n)
	fday = fdayv(n)

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	end

c***********************************************************

	function it2n(t)

c converts time [day] to pointer into array [1-365]

	implicit none

	integer it2n
	real t

	integer n

	n = t
	if( t .lt. 0. ) n = n - 1
	n = mod(n,365)
	if( n .lt. 0 ) n = n + 365	!handle negative days
	n = n + 1			!n is in [1,365]

	it2n = n

	end

c***********************************************************

c	call test_solrad
c	end

c***********************************************************

