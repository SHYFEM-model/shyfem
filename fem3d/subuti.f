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
c
c******************************************

	subroutine getxy(k,x,y)

c gets coordinates x/y for node k

	implicit none

	integer k
	real x,y

	real xgv(1), ygv(1)
	common /xgv/xgv, /ygv/ygv

	x = xgv(k)
	y = ygv(k)

	end

c******************************************

	subroutine getexy(ie,x,y)

c gets coordinates x/y for element ie

	implicit none

	integer ie
	real x(3), y(3)

	integer k,ii

	integer nen3v(3,1)
	common /nen3v/nen3v
	real xgv(1), ygv(1)
	common /xgv/xgv, /ygv/ygv

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

	implicit none

c arguments
	real areael
	integer ie
c common
	integer nen3v(3,1)
	real xgv(1),ygv(1)
	common /nen3v/nen3v
	common /xgv/xgv,/ygv/ygv
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

	implicit none

c arguments
	real areavl
	integer k
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	include 'links.h'
	integer nen3v(3,1)
	common /nen3v/nen3v
	include 'ev.h'
c local
	logical blink
	integer ie,ii,i,nl
	real area

	blink = .true.
	area=0.

	if( blink ) then

	  call set_elem_links(k,nl)

          do i=1,nl
            ie=lnk_elems(i)
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

        implicit none

c arguments
        real flxnod
        integer k
c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real unv(1),vnv(1)
        common /unv/unv, /vnv/vnv
        integer nen3v(3,1)
        common /nen3v/nen3v
	include 'ev.h'
	include 'links.h'
c local
        real flux
        integer i,nl,ie,ii
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

	  call set_elem_links(k,nl)

          do i=1,nl
            ie=lnk_elems(i)
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
c
c computation of kinetic & potential energy (per unit mass) in element ielem
c
c for ielem=0 --> total kinetic & potential energy
c
c	pot = (g/2) * aj * 4  * z(m)*z(m)
c	kin = (1/2) * aj * 12 * u*u/h
c
	implicit none

	include 'param.h'

c arguments
	integer ielem
	real kenerg,penerg
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real grav,fcor,dcor,dirn,rowass,roluft
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft
	integer nen3v(3,1)
	real hev(1),znv(1)
	real unv(1),vnv(1)
	common /nen3v/nen3v
	common /hev/hev, /znv/znv
	common /unv/unv ,/vnv/vnv
	include 'ev.h'

        real zeov(3,neldim),zenv(3,neldim)
        common /zeov/zeov, /zenv/zenv

c local
	integer ie,ii
	real*8 aj,pot,kin,z,zz
c
	ie=ielem
c
	if(ie.gt.0.and.ie.le.nel) then
	  z=0.
	  aj=ev(10,ie)
	  do ii=1,3
	    zz = zenv(ii,ie)
	    z=z+zz*zz
	  end do
	  pot=aj*z
	  kin=aj*(unv(ie)**2+vnv(ie)**2)/hev(ie)
	else
	  kin=0.
	  pot=0.
	  do ie=1,nel
	    z=0.
	    aj=ev(10,ie)
	    do ii=1,3
	      zz = zenv(ii,ie)
	      z=z+zz*zz
	    end do
	    pot=pot+aj*z
	    kin=kin+aj*(unv(ie)**2+vnv(ie)**2)/hev(ie)
	  end do
	end if
c
	penerg=2.*grav*pot	! 2 = 12 / 3 / 2
	kenerg=6.*kin		! 6 = 12 / 2
c
	return
	end

c***************************************************************

	subroutine energ3d(kenergy,penergy)

c computation of kinetic & potential energy (per unit mass) in element ielem
c
c for ielem=0 --> total kinetic & potential energy
c
c	pot = (g/2) * aj * 4  * z(m)*z(m)
c	kin = (1/2) * aj * 12 * u*u/h

	implicit none

	real kenergy,penergy

	include 'param.h'
	include 'ev.h'

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real grav,fcor,dcor,dirn,rowass,roluft
	common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft

	integer ilhv(1)
	common /ilhv/ilhv
        real hdenv(nlvdim,neldim)
        common /hdenv/hdenv

        real zeov(3,neldim),zenv(3,neldim)
        common /zeov/zeov, /zenv/zenv
        real utlnv(nlvdim,neldim)
        common /utlnv/utlnv
        real vtlnv(nlvdim,neldim)
        common /vtlnv/vtlnv

c local
	integer ie,ii,l,lmax
	double precision area,pot,kin,z,zz
	double precision h,uu,vv

	kin=0.
	pot=0.

	do ie=1,nel

	  area = 12. * ev(10,ie)

	  z=0.
	  do ii=1,3
	    zz = zenv(ii,ie)
	    z = z + zz*zz
	  end do
          pot=pot+area*z/3.

	  lmax = ilhv(ie)
	  do l=1,lmax
	    h = hdenv(l,ie)
	    uu = utlnv(l,ie)
	    vv = vtlnv(l,ie)
	    kin = kin + area * (uu*uu + vv*vv) / h
	    !write(13,*) ie,l,uu,vv,h,kin
	  end do

	end do

        penergy = 0.5*grav*pot      ! 0.5 = 1 / 2
        kenergy = 0.5*kin           ! 0.5 = 1 / 2

	!write(6,*) kenergy,penergy

	end

c***************************************************************

        subroutine stmima(a,nkn,nlvdim,ilhkv,amin,amax)

c computes min/max of 3d field

        implicit none

c arguments
        integer nkn,nlvdim
        real a(nlvdim,1)
        integer ilhkv(1)
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

        implicit none

        real cn(1)
        real ce(3,1)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        integer nen3v(3,1)
        common /nen3v/nen3v
        integer iwegv(1)
        common /iwegv/iwegv
	include 'ev.h'

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

