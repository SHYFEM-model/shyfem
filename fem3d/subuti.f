
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

c utility routines for 2d/3d model
c
c contents :
c
c getxy(k,x,y)					coordinates x/y for node k
c getexy(ie,x,y)				coordinates x/y for element ie
c
c real function areael(ie)			area for element ie
c real function areavl(k)			area for finite volume k
c
c function flxnod(k)            		discharge around node k
c
c subroutine energ(ielem,kenerg,penerg)		kinetic & potential energy
c
c subroutine stmima(a,nkn,nlvdi,ilhkv,amin,amax)
c                                       computes min/max of 3d field
c subroutine n2ebar(cn,ce)
c		copies concentrations from node value to element value
c
c revision log :
c
c 19.08.1998	ggu	new routines volco0, volno0
c 26.08.1998	ggu	subroutine stmima transferred from newbcl0
c 18.12.1999	ggu	/iweich/ -> /iwegv/ (bug)
c 29.03.2000	ggu	new routine getxy and getexy
c 16.05.2000	ggu	routine volel removed
c 28.04.2009	ggu	links re-structured
c 23.03.2010	ggu	changed v6.1.1
c 08.06.2010	ggu	new routine for computing 3D kin/pot energy
c 22.07.2010	ggu	changed VERS_6_1_9
c 15.12.2010	ggu	changed VERS_6_1_14
c 07.07.2011	ggu	bug fix in areael (unstable computation)
c 14.07.2011	ggu	changed VERS_6_1_27
c 23.12.2014	ggu	changed VERS_7_0_11
c 19.01.2015	ggu	changed VERS_7_1_2
c 19.01.2015	ggu	changed VERS_7_1_3
c 01.04.2015	ggu	changed VERS_7_1_7
c 29.04.2015	ggu	energy now in Joule
c 21.05.2015	ggu	changed VERS_7_1_11
c 05.06.2015	ggu	changed VERS_7_1_12
c 10.07.2015	ggu	changed VERS_7_1_50
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 24.07.2015	ggu	changed VERS_7_1_82
c 16.12.2015	ggu	changed VERS_7_3_16
c 18.12.2015	ggu	changed VERS_7_3_17
c 05.12.2017	ggu	changed VERS_7_5_39
c 22.02.2018	ggu	changed VERS_7_5_42
c 19.04.2018	ggu	changed VERS_7_5_45
c 16.02.2019	ggu	changed VERS_7_5_60
c
c******************************************

	subroutine getxy(k,x,y)

c gets coordinates x/y for node k

	use basin

	implicit none

	integer k
	real x,y

	x = xgv(k)
	y = ygv(k)

	end

c******************************************

	subroutine getexy(ie,x,y)

c gets coordinates x/y for element ie

	use basin

	implicit none

	integer ie
	real x(3), y(3)

	integer k,ii

	do ii=1,3
	  k = nen3v(ii,ie)
	  x(ii) = xgv(k)
	  y(ii) = ygv(k)
	end do

	end

c******************************************

	function areael(ie)

c area for element ie
c
c double precision version - bug fix 07.07.2011

	use basin

	implicit none

c arguments
	real areael
	integer ie
c local
	integer kn1,kn2,kn3
	real*8 x1,x2,x3,y1,y2,y3
	real*8 a1,a2,a3
	real*8 half

	half = 0.5

	kn1=nen3v(1,ie)
	kn2=nen3v(2,ie)
	kn3=nen3v(3,ie)

	x1=xgv(kn1)
	y1=ygv(kn1)
	x2=xgv(kn2)
	y2=ygv(kn2)
	x3=xgv(kn3)
	y3=ygv(kn3)

	!a1=x2*y3-x3*y2
	!a2=x3*y1-x1*y3
	!a3=x1*y2-x2*y1
	!areael = half*(a1+a2+a3)

	areael = half * ( (x2-x1) * (y3-y1) - (x3-x1) * (y2-y1) )

	end

c******************************************

	function areavl(k)

c area for finite volume k

	use mod_geom
	use evgeom
	use basin

	implicit none

c arguments
	real areavl
	integer k
c local
	logical blink
	integer ie,ii,i,nl
	integer elems(maxlnk)
	real area

	blink = .true.
	area=0.

	if( blink ) then

	  call get_elems_around(k,maxlnk,nl,elems)

          do i=1,nl
            ie=elems(i)
            area=area+ev(10,ie)
          end do

	else

	  do ie=1,nel
	    do ii=1,3
	      if(nen3v(ii,ie).eq.k) then
		area=area+ev(10,ie)
	      end if
	    end do
	  end do

	end if

	areavl = 4. * area

	end

c*****************************************************************

        function flxnod(k)

c computes discharge (m**3/sec) into finite volume around node k
c ... value flxnod has to be multiplied by dt to obtain total volume
c
c depending on value of blink uses link structure or not
c
c discharge into node n:     Q = 12 * aj * ( b(n)*U + c(n)*V )
c volume difference:         dV = dt * Q

	use mod_geom
	use mod_hydro_baro
	use evgeom
	use basin

        implicit none

c arguments
        real flxnod
        integer k
c local
        real flux
        integer i,nl,ie,ii
	integer elems(maxlnk)
        logical blink
c function
        integer ithis
c statement function
        real fflux
        fflux(ii,ie) = ev(10,ie) *
     +                  ( unv(ie)*ev(3+ii,ie) + vnv(ie)*ev(6+ii,ie) )

        blink = .true.
        flux=0.

        if( blink ) then

	  call get_elems_around(k,maxlnk,nl,elems)

          do i=1,nl
            ie=elems(i)
            ii=ithis(k,ie)
            flux = flux + fflux(ii,ie)
          end do

        else

          do ie=1,nel
            do ii=1,3
              if( nen3v(ii,ie) .eq. k ) then
                flux = flux + fflux(ii,ie)
              end if
            end do
          end do

        end if

        flxnod = 12. * flux

        end

c*****************************************************************

	subroutine energ(ielem,kenerg,penerg)

c computation of kinetic & potential energy [Joule]
c
c	pot = (g/2) * rho * area * z*z
c	kin = (1/2) * rho * area * (U*U+V*V)/H

	use mod_depth
	use mod_hydro_baro
	use mod_hydro
	use evgeom
	use basin

	implicit none

	integer ielem		!element (0 for total energy - all elements)
	real kenerg		!kinetic energy (return)
	real penerg		!potential energy relative to z=0 (return)

	include 'pkonst.h'

	integer ie,ii,ie1,ie2
	double precision aj,pot,kin,z,zz

	if(ielem.gt.0.and.ielem.le.nel) then
	  ie1 = ielem
	  ie2 = ielem
	else
	  ie1 = 1
	  ie2 = nel
	end if

	kin=0.
	pot=0.
	do ie=ie1,ie2
	  zz=0.
	  aj=ev(10,ie)
	  do ii=1,3
	    z = zenv(ii,ie)
	    zz = zz + z*z
	  end do
	  pot=pot+aj*zz
	  kin=kin+aj*(unv(ie)**2+vnv(ie)**2)/hev(ie)
	end do

	penerg=2.*grav*pot	! 2 = 12 / 3 / 2
	kenerg=6.*kin		! 6 = 12 / 2

	end

c***************************************************************

	subroutine energ3d(kenergy,penergy,ksurf,ia_ignore)

c computation of kinetic & potential energy [Joule]
c
c	pot = (g/2) * rho * area * z*z
c	kin = (1/2) * rho * area * (U*U+V*V)/H

	use mod_layer_thickness
	use mod_ts
	use mod_hydro
	use evgeom
	use levels
	use basin
	use shympi

	implicit none

	real kenergy		!kinetic energy (return)
	real penergy		!potential energy relative to z=0 (return)
	real ksurf		!kinetic energy of surface layer (return)
	integer ia_ignore	!area code to be ignored

	include 'pkonst.h'

	integer ie,ii,l,lmax,ia,k,ntot
	double precision area,pot,kin,kinsurf,z,zz,ke
	double precision h,uu,vv,rho

	pot=0.
	kin=0.
	kinsurf=0.

	ntot = nel
	ntot = nel_unique

	do ie=1,ntot

	  area = 12. * ev(10,ie)
	  ia = iarv(ie)
	  if( ia .eq. ia_ignore ) cycle

	  zz=0.
	  rho = 0.
	  do ii=1,3
	    k = nen3v(ii,ie)
	    rho = rho + rhov(1,k)
	    z = zenv(ii,ie)
	    zz = zz + z*z
	  end do
	  rho = rowass + rho/3.
	  zz = zz / 3.
          pot = pot + area * rho * zz

	  lmax = ilhv(ie)
	  do l=1,lmax
	    rho = 0.
	    do ii=1,3
	      k = nen3v(ii,ie)
	      rho = rho + rhov(l,k)
	    end do
	    rho = rowass + rho/3.
	    h = hdenv(l,ie)
	    uu = utlnv(l,ie)
	    vv = vtlnv(l,ie)
	    ke = area * rho * (uu*uu + vv*vv) / h
	    kin = kin + ke
	    if( l == 1 ) kinsurf = kinsurf + ke
	  end do

	end do

        penergy = 0.5*grav*pot
        kenergy = 0.5*kin
	ksurf = 0.5*kinsurf

	end

c***************************************************************

        subroutine stmima(a,nkn,nlvddi,ilhkv,amin,amax)

c computes min/max of 3d field

        implicit none

c arguments
        integer nkn,nlvddi
        real a(nlvddi,nkn)
        integer ilhkv(nkn)
        real amin,amax
c local
        integer lmax,k,l

        amin=a(1,1)
        amax=amin

        do k=1,nkn
          lmax=ilhkv(k)
          do l=1,lmax
            if(a(l,k).gt.amax) amax=a(l,k)
            if(a(l,k).lt.amin) amin=a(l,k)
          end do
        end do

        end

c**********************************************************************

        subroutine n2ebar(cn,ce)

c copies concentrations from node value to element value (only wet areas)

	use mod_geom_dynamic
	use evgeom
	use basin

        implicit none

        real cn(nkn)
        real ce(3,nel)

        integer ie,ii,k

        do ie=1,nel
         if( iwegv(ie) .eq. 0 ) then   !wet
          do ii=1,3
            k = nen3v(ii,ie)
            ce(ii,ie) = cn(k)
          end do
         end if
        end do

        end

c******************************************************

