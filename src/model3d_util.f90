!
! $Id: subuti.f,v 1.11 2009-05-21 09:24:00 georg Exp $
!
! utility routines for 2d/3d model
!
! contents :
!
! getxy(k,x,y)					coordinates x/y for node k
! getexy(ie,x,y)				coordinates x/y for element ie
!
! double precision function areael(ie)			area for element ie
! double precision function areavl(k)			area for finite volume k
!
! function flxnod(k)            		discharge around node k
!
! subroutine energ(ielem,kenerg,penerg)		kinetic & potential energy
!
! subroutine stmima(a,nkn,nlvdi,ilhkv,amin,amax)
!                                       computes min/max of 3d field
! subroutine n2ebar(cn,ce)
!		copies concentrations from node value to element value
!
! revision log :
!
! 19.08.1998	ggu	new routines volco0, volno0
! 26.08.1998	ggu	subroutine stmima transferred from newbcl0
! 18.12.1999	ggu	/iweich/ -> /iwegv/ (bug)
! 29.03.2000	ggu	new routine getxy and getexy
! 16.05.2000	ggu	routine volel removed
! 28.04.2009    ggu     links re-structured
! 08.06.2010    ggu     new routine for computing 3D kin/pot energy
! 07.07.2011    ggu     bug fix in areael (unstable computation)
! 29.04.2015    ggu     energy now in Joule
!
!******************************************
!-----------------------------------------------------------------------
        module model3d_util
!-----------------------------------------------------------------------
        contains
!-----------------------------------------------------------------------

	subroutine getxy(k,x,y)

! gets coordinates x/y for node k

	use basin

	implicit none

	integer k
	double precision x,y

	x = xgv(k)
	y = ygv(k)

	end

!******************************************

	subroutine getexy(ie,x,y)

! gets coordinates x/y for element ie

	use basin

	implicit none

	integer ie
	double precision x(3), y(3)

	integer k,ii

	do ii=1,3
	  k = nen3v(ii,ie)
	  x(ii) = xgv(k)
	  y(ii) = ygv(k)
	end do

	end

!******************************************

	function areael(ie)

! area for element ie
!
! double precision version - bug fix 07.07.2011

	use basin

	implicit none

! arguments
	double precision areael
	integer ie
! local
	integer kn1,kn2,kn3
	double precision x1,x2,x3,y1,y2,y3
	double precision a1,a2,a3
	double precision half

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

!******************************************

	function areavl(k)

! area for finite volume k

	use geom
	use evgeom
	use basin
        use lnku

	implicit none

! arguments
	double precision areavl
	integer k
! local
	logical blink
	integer ie,ii,i,nl
	integer elems(maxlnk)
	double precision area

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

!*****************************************************************

        function flxnod(k)

! computes discharge (m**3/sec) into finite volume around node k
! ... value flxnod has to be multiplied by dt to obtain total volume
!
! depending on value of blink uses link structure or not
!
! discharge into node n:     Q = 12 * aj * ( b(n)*U + c(n)*V )
! volume difference:         dV = dt * Q

	use geom
	use hydro_baro
	use evgeom
	use basin
        use lnku

        implicit none

! arguments
        double precision flxnod
        integer k
! local
        double precision flux
        integer i,nl,ie,ii
	integer elems(maxlnk)
        logical blink
! statement function
        double precision fflux
        fflux(ii,ie) = ev(10,ie) * ( unv(ie)*ev(3+ii,ie) + vnv(ie)*ev(6+ii,ie) )

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

!*****************************************************************

	subroutine energ(ielem,kenerg,penerg)

! computation of kinetic & potential energy [Joule]
!
!	pot = (g/2) * rho * area * z*z
!	kin = (1/2) * rho * area * (U*U+V*V)/H

	use depth
	use hydro_baro
	use hydro_admin
	use evgeom
	use basin

	implicit none

	integer ielem		!element (0 for total energy - all elements)
	double precision kenerg		!kinetic energy (return)
	double precision penerg		!potential energy relative to z=0 (return)

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

!***************************************************************

	subroutine energ3d(kenergy,penergy,ia_ignore)

! computation of kinetic & potential energy [Joule]
!
!	pot = (g/2) * rho * area * z*z
!	kin = (1/2) * rho * area * (U*U+V*V)/H

	use layer_thickness
	use ts
	use hydro_admin
	use evgeom
	use levels
	use basin
        use shympi      !ivb

	implicit none

	double precision kenergy		!kinetic energy (return)
	double precision penergy		!potential energy relative to z=0 (return)
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

        penergy = shympi_sum(penergy)
        kenergy = shympi_sum(kenergy)

	end

!***************************************************************

        subroutine stmima(a,nkn,nlvddi,ilhkv,amin,amax)

! computes min/max of 3d field

        use shympi

        implicit none

! arguments
        integer nkn,nlvddi
        double precision a(nlvddi,nkn)
        integer ilhkv(nkn)
        double precision amin,amax
! local
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

        amin = shympi_min(amin)
        amax = shympi_max(amax)

        end

!**********************************************************************

        subroutine n2ebar(cn,ce)

! copies concentrations from node value to element value (only wet areas)

	use geom_dynamic
	use evgeom
	use basin

        implicit none

        double precision cn(nkn)
        double precision ce(3,nel)

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

!******************************************************

        !ivb
!        subroutine energ3d2(kenergyx,kenergyy,kenergy,penergy)
!
!c computation of kinetic & potential energy [Joule]
!c
!c potential = int_V g * rho * z * dV
!c kinetic   = int_V (1/2) * rho * dV * (u*u+u*u)
!
!        use layer_thickness
!        use ts
!        use hydro_admin
!        use evgeom
!        use levels
!        use basin
!        use shympi
!        !use elems_dealing
!
!        implicit none
!
!        double precision kenergy            !kinetic energy (return)
!        double precision kenergyx           !kinetic energy x direction (return)
!        double precision kenergyy           !kinetic energy y direction (return)
!        double precision penergy            !potential energy relative to z=0 (return)
!
!        include 'pkonst.h'
!
!        integer ie,ii,l,lmax,ia,k
!        double precision area,pot,kin,z,zz,kinx,kiny
!        double precision h,uu,vv,rho
!        double precision vol
!        double precision areanode
!
!        kin=0.
!        pot=0.
!
!        ! I use volumes around density points
!        do k=1,nkn
!          zz = 20.0 + zov(k) ! the bottom (z=-20) is the reference
!          lmax = ilhkv(k)
!          do l=1,lmax
!            zz   = zz - 0.5*hdkov(l,k)
!            area = areanode(l,k,+1)
!            rho  = rowass + rhov(l,k)
!            vol  = area * hdknv(l,k)
!            pot  = pot + rho * zz * vol
!            zz   = zz - 0.5*hdknv(l,k)
!          end do
!        end do
!
!        penergy = grav*pot
!
!        do ie=1,nel
!          area = 12. * ev(10,ie)
!          lmax = ilhv(ie)
!          do l=1,lmax
!            rho = 0.
!            do ii=1,3
!              k = nen3v(ii,ie)
!              rho = rho + rhov(l,k)
!            end do
!            rho = rowass + rho/3.
!            h   = hdenv(l,ie)
!            uu  = utlnv(l,ie)
!            vv  = vtlnv(l,ie)
!            kinx= kinx + area * rho * (uu*uu) / h
!            kiny= kiny + area * rho * (vv*vv) / h
!            kin = kin + area * rho * (uu*uu + vv*vv) / h
!          end do
!        end do
!
!        kenergyx= 0.5*kinx
!        kenergyy= 0.5*kiny
!        kenergy = 0.5*kin
!
!        penergy = shympi_sum(penergy)
!        kenergy = shympi_sum(kenergy)
!        kenergyx= shympi_sum(kenergyx)
!        kenergyy= shympi_sum(kenergyy)
!
!        end
!-----------------------------------------------------------------------
        end module model3d_util
!-----------------------------------------------------------------------
