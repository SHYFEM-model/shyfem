c
c $Id: subuti.f,v 1.11 2009-05-21 09:24:00 georg Exp $
c
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
c 28.04.2009    ggu     links re-structured
c 08.06.2010    ggu     new routine for computing 3D kin/pot energy
c 07.07.2011    ggu     bug fix in areael (unstable computation)
c 29.04.2015    ggu     energy now in Joule
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

	subroutine energ3d(kenergy,penergy,ia_ignore)

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

	implicit none

	real kenergy		!kinetic energy (return)
	real penergy		!potential energy relative to z=0 (return)
	integer ia_ignore	!area code to be ignored

	include 'pkonst.h'

	integer ie,ii,l,lmax,ia,k
	double precision area,pot,kin,z,zz
	double precision h,uu,vv,rho

	kin=0.
	pot=0.

	do ie=1,nel

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
	    kin = kin + area * rho * (uu*uu + vv*vv) / h
	  end do

	end do

        penergy = 0.5*grav*pot
        kenergy = 0.5*kin

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

