c
c $Id: newtra.f,v 1.22 2010-03-11 15:36:38 georg Exp $
c
c routines for various transformations
c
c contents :
c
c subroutine ttov			transforms transports to velocities
c subroutine vtot			transforms velocities to transports
c subroutine uvtop0			transforms bar. transp. to nod. veloc.
c subroutine uvtopr			transforms velocities to nodal values
c subroutine prtouv			transforms nodal values to element vel.
c subroutine uvint			computation of barotropic transports
c subroutine austau(vv)			computes aux vectors for austausch term
c subroutine baro2l			distribute barotropic velocities 
c subroutine setxv			sets obsolete data structure xv
c subroutine getuv(l,k,u,v)             accessor routine to get velocities u/v
c subroutine copy_uvz                   copies u/v/z to old time step
c subroutine make_prvel                 makes print velocities and xv
c subroutine init_uv                    initializes uvz values
c subroutine e2n2d(elv,nov,aux)         transforms element to nodal values
c subroutine n2e2d(nov,elv)             transforms nodal to element values 2D
c subroutine n2e3d(nlvdi,nov,elv)       transforms nodal to element values 3D
c
c revision log :
c
c 10.08.2003    ggu     new routines copy_uvz, make_prvel, init_uvz
c 18.09.2003    ggu     new routine e2n2d
c 04.12.2003    ggu     new routine n2e2d
c 30.03.2004	ccf	new routine n2e2d and n2e3d
c 15.10.2004	ggu	new routine smagorinsky started
c 14.01.2005	ggu	bug (nlvdi) in n2e3d fixed
c 24.02.2005	ggu	smagorinsky into subdif.f
c 15.03.2005	ggu	austv eliminated, austau() into subdif.f
c 27.06.2005	ggu	bug in vtot corrected
c 10.04.2008	ggu	copy velocities at nodes in copy_uvz()
c 01.03.2010	ggu	new version of n2e3d()
c 11.03.2010	ggu	new routine check_volume(); init w only if no restart
c 16.02.2011	ggu	new routine e2n3d() and e2n3d_minmax()
c 27.01.2012	deb&ggu	routines adapted for sigma levels
c 03.12.2015	ccf&ggu	code optimized
c 07.04.2016	ggu	new routine aver_nodal()
c 19.05.2016	ggu	use where construct where possible
c
c****************************************************************************

	subroutine ttov

c transforms transports to velocities

	use mod_layer_thickness
	use mod_hydro_vel
	use mod_hydro
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	where( hdenv > 0 )
	  ulnv = utlnv / hdenv
	  vlnv = vtlnv / hdenv
	end where

	end

c******************************************************************

	subroutine vtot

c transforms velocities to transports

	use mod_layer_thickness
	use mod_hydro_vel
	use mod_hydro
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	utlnv = ulnv * hdenv
	vtlnv = vlnv * hdenv

	end

c******************************************************************
c
	subroutine uvtop0
c
c transforms barotropic transports to nodal velocities
c
	use mod_geom_dynamic
	use mod_depth
	use mod_hydro_baro
	use mod_hydro_print
	use mod_hydro
	use evgeom
	use basin

	implicit none

c local
	logical bcolin
	integer ie,k,ii
	real aj,zm,hm
	real vv(nkn)
c function
	real getpar
	integer iround
c
	bcolin=iround(getpar('iclin')).ne.0
c
	up0v = 0.
	vp0v = 0.
	vv   = 0.
c
	do ie=1,nel
	  if( iwegv(ie) /= 0 ) cycle
	  aj=ev(10,ie)
	  zm=0.
	  do ii=1,3
	    zm=zm+zenv(ii,ie)
	  end do
	  zm=zm/3.
	  hm=hev(ie)
	  if(.not.bcolin) hm=hm+zm
	  do ii=1,3
	    k=nen3v(ii,ie)
	    vv(k)=vv(k)+aj
	    up0v(k)=up0v(k)+aj*unv(ie)/hm
	    vp0v(k)=vp0v(k)+aj*vnv(ie)/hm
	  end do
	end do

	where ( vv > 0. ) 
          up0v = up0v / vv
          vp0v = vp0v / vv
	end where

	return
	end

c******************************************************************

	subroutine uvtopr

c transforms velocities to nodal values

	use mod_geom_dynamic
	use mod_hydro_print
	use mod_hydro_vel
	use evgeom
	use levels
	use basin

	implicit none

	integer ie,l,k,ii
	integer lmax
	real aj
	!real vv(nlvdi,nkn)
	real, allocatable :: vv(:,:)

	allocate(vv(nlvdi,nkn))
	uprv = 0.
	vprv = 0.
	vv   = 0.

c baroclinic part

	do ie=1,nel
	  if ( iwegv(ie) /= 0 ) cycle
          lmax = ilhv(ie)
	  aj=ev(10,ie)
	  do l=1,lmax
	    do ii=1,3
	      k=nen3v(ii,ie)
	      vv(l,k)=vv(l,k)+aj
	      uprv(l,k)=uprv(l,k)+aj*ulnv(l,ie)
	      vprv(l,k)=vprv(l,k)+aj*vlnv(l,ie)
	    end do
	  end do
	end do

	where ( vv > 0. ) 
	  uprv = uprv / vv
	  vprv = vprv / vv
	end where

c vertical velocities -> we compute average over one layer

	do l=1,nlv
	  wprv(l,:)=0.5*(wlnv(l,:)+wlnv(l-1,:))
	end do

	wprv(0,:) = 0.

	deallocate(vv)

	end

c******************************************************************
c
	subroutine prtouv
c
c transforms nodal values to element values (velocities)
c
	use mod_geom_dynamic
	use mod_hydro_print
	use mod_hydro_vel
	use levels
	use basin

	implicit none

	integer ie,l,k,ii
	real u,v
c
c baroclinic part
c
	do ie=1,nel
	 if( iwegv(ie) .eq. 0 ) then
	  do l=1,ilhv(ie)
	    u=0.
	    v=0.
	    do ii=1,3
	      k=nen3v(ii,ie)
	      u=u+uprv(l,k)
	      v=v+vprv(l,k)
	    end do
	    ulnv(l,ie)=u/3.
	    vlnv(l,ie)=v/3.
	  end do
	 else
	    ulnv(:,ie)=0.
	    vlnv(:,ie)=0.
	 end if
	end do
c
c vertical velocities -> from layer average to interface values
c
	wlnv(nlv,:) = 0.
	do l=nlv-1,0,-1
	  wlnv(l,:)=2.*wprv(l,:)-wlnv(l+1,:)
	end do
c
	return
	end
c
c*****************************************************************
c
	subroutine uvint
c
c computation of barotropic part of transports
c
	use mod_hydro_baro
	use mod_hydro
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer ie,l
	real u,v
c
	do ie=1,nel
	  u=0.
	  v=0.
	  do l=1,ilhv(ie)
	    u=u+utlnv(l,ie)
	    v=v+vtlnv(l,ie)
	  end do
	  unv(ie)=u
	  vnv(ie)=v
	end do
c
	return
	end
c
c*****************************************************************

	subroutine check_volume

c checks for negative volume (depth)
c
c only first layer has to be checked

	use mod_layer_thickness
	use mod_hydro
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	logical bstop,bsigma
	integer nsigma
	integer k,ie,ke,iee,ii
	real h,z
	real hsigma

	integer ipext,ieext

	bstop = .false.

	call get_sigma(nsigma,hsigma)
	bsigma = nsigma .gt. 0

	if( .not. bsigma ) then		!not sure we need this - FIXME
	 do k=1,nkn
	  h = hdknv(1,k)
	  if( h .le. 0. ) then
	    z = znv(k)
	    ke = ipext(k)
	    write(6,*) 'negative depth in node (layer 1): '
	    write(6,*) '   node (int/ext): ',k,ke
	    write(6,*) '   depth,zeta    : ',h,z
	    bstop = .true.
	  end if
	 end do
	end if

	do ie=1,nel
	  h = hdenv(1,ie)
	  if( h .le. 0. ) then
	    iee = ieext(ie)
	    write(6,*) 'negative depth in elem (layer 1): '
	    write(6,*) '   element (int/ext): ',ie,iee
	    write(6,*) '   depth: ',h
	    write(6,*) '   water levels: ',(zenv(ii,ie),ii=1,3)
	    bstop = .true.
	  end if
	end do

	if( bstop ) then
	  write(6,*) 'Negative volumes in first layer found.'
	  write(6,*) 'Drying of an element with more than one layer'
	  write(6,*) 'is not allowed.'
	  write(6,*) 'This can be due to either a computed water level'
	  write(6,*) 'which is too low or because of instability.'
	  write(6,*) 'In the first case check the boundary data and/or'
	  write(6,*) 'increase the layer thickness of the first layer.'
	  write(6,*) 'In the second case look for the causes of the'
	  write(6,*) 'instability and maybe reduce the time step.'
	  stop 'error stop check_volume: negative volume'
	end if

	end

c*****************************************************************
c
	subroutine baro2l

c distribute barotropic velocities onto layers (only in dry elements)

	use mod_geom_dynamic
	use mod_hydro_baro
	use mod_hydro_vel
	use mod_hydro
	use levels
	use basin

	implicit none

	logical bsigma
	integer nsigma
	integer ie,ilevel,ii,l
	real hsigma
c functions
        integer ieext

	call get_sigma(nsigma,hsigma)
	bsigma = nsigma .gt. 0

	if( bsigma ) return

	do ie=1,nel
	  if( iwegv(ie) .gt. 0 ) then	!dry
	    ilevel=ilhv(ie)
	    if( ilevel .gt. 1 ) then
		write(6,*) 'drying in more than one layer'
		write(6,*) ie,ieext(ie),ilevel
		do ii=1,3
		  write(6,*) zeov(ii,ie),zenv(ii,ie),hm3v(ii,ie)
		end do
                do l=1,ilevel
                  write(6,*) utlnv(l,ie),vtlnv(l,ie)
                end do
                write(6,*) unv(ie),vnv(ie)
		stop 'error stop baro2l: drying in more than one layer'
	    end if
	    utlnv(1,ie) = unv(ie)
	    vtlnv(1,ie) = vnv(ie)
	  end if
	end do

	end

c******************************************************************

	subroutine setxv

c sets obsolete data structure xv

	use mod_hydro_print
	use mod_hydro
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	xv(1,:) = up0v(:)
	xv(2,:) = vp0v(:)
	xv(3,:) = znv(:)

	end

c******************************************************************

	subroutine getuv(l,k,u,v)

c accessor routine to get velocities u/v

	use mod_hydro_print

	implicit none

        integer l,k
        real u,v

        u = uprv(l,k)
        v = vprv(l,k)

        end

c******************************************************************

	subroutine copy_uvz

c copies u/v/z to old time step

	use mod_hydro_baro
	use mod_hydro_print
	use mod_hydro_vel
	use mod_hydro

	implicit none

	zov   = znv
        upro  = uprv
        vpro  = vprv
	zeov  = zenv
	uov   = unv			!$$UVBARO
	vov   = vnv
	ulov  = ulnv
	vlov  = vlnv
	utlov = utlnv
	vtlov = vtlnv
	wlov  = wlnv

	end

c******************************************************************

	subroutine make_prvel

c makes print velocities and xv from new level arrays

	use basin

	implicit none

	call uvtopr
	call uvtop0
	call setxv

	end

c******************************************************************

	subroutine init_uv

c initializes uvz values from zenv, utlnv, vtlnv, hdenv

	use basin

	implicit none

	logical rst_use_restart
	real dzeta(nkn)

	call ttov			!velocities from layer transports
	call uvint			!barotropic transports

	if( .not. rst_use_restart(5) ) then
	  call hydro_vertical(dzeta)	!vertical velocities
	end if

	call make_prvel			!nodal values

	call copy_uvz			!copy to old time step

	end

c******************************************************************
c******************************************************************
c******************************************************************

	subroutine e2n2d(elv,nov,aux)

c transforms element values to nodal values (weights are area)
c
c (2D version)

	use evgeom
	use basin

	implicit none

c arguments
        real elv(nel)     !array with element values (in)
        real nov(nkn)     !array with nodal values (out)
        real aux(nkn)     !aux array (nkn)

c local
        integer k,ie,ii
        real area,value

c-----------------------------------------------------------
c initialize arrays
c-----------------------------------------------------------

        nov = 0.
        aux = 0.

c-----------------------------------------------------------
c accumulate values
c-----------------------------------------------------------

        do ie=1,nel
          area = 4.*ev(10,ie)
          value = elv(ie)
          do ii=1,3
            k = nen3v(ii,ie)
            nov(k) = nov(k) + area*value
            aux(k) = aux(k) + area
          end do
        end do

c-----------------------------------------------------------
c compute final value
c-----------------------------------------------------------

        nov = nov / aux

c-----------------------------------------------------------
c end of routine
c-----------------------------------------------------------

        end

c******************************************************************

	subroutine e2n3d(nlvddi,elv,nov,aux)

c transforms element values to nodal values (weights are area)
c
c (3D version)

	use evgeom
	use levels
	use basin

	implicit none

c arguments
	integer nlvddi
        real elv(nlvddi,nel)     !array with element values (in)
        real nov(nlvddi,nkn)     !array with nodal values (out)
        real aux(nlvddi,nkn)     !aux array (nkn)

c local
        integer k,ie,ii,l,lmax
        real area,value

c-----------------------------------------------------------
c initialize arrays
c-----------------------------------------------------------

        nov = 0.
        aux = 0.

c-----------------------------------------------------------
c accumulate values
c-----------------------------------------------------------

        do ie=1,nel
          area = 4.*ev(10,ie)
	  lmax = ilhv(ie)
	  do l=1,lmax
            value = elv(l,ie)
            do ii=1,3
              k = nen3v(ii,ie)
              nov(l,k) = nov(l,k) + area*value
              aux(l,k) = aux(l,k) + area
	    end do
          end do
        end do

c-----------------------------------------------------------
c compute final value
c-----------------------------------------------------------

	where ( aux > 0. ) 
	  nov = nov / aux
	end where

c-----------------------------------------------------------
c end of routine
c-----------------------------------------------------------

        end

c******************************************************************

	subroutine e2n3d_minmax(mode,nlvddi,elv,nov)

c transforms element values to nodal values (no weights - use min/max)
c
c (3D version)

	use levels
	use basin

	implicit none

c arguments
	integer mode		!min (-1) or max (+1)
	integer nlvddi		!vertical dimension
        real elv(nlvddi,nel)      !array with element values (in)
        real nov(nlvddi,nkn)      !array with nodal values (out)

c local
        integer k,ie,ii,l,lmax
        real rinit,value

c-----------------------------------------------------------
c initialize arrays
c-----------------------------------------------------------

	if( mode .eq. -1) then
	  rinit = 1.e+30
	else if( mode .eq. 1 ) then
	  rinit = -1.e+30
	else
	  stop 'error stop e2n3d_minmax: unknown mode'
	end if

        nov = rinit

c-----------------------------------------------------------
c accumulate values
c-----------------------------------------------------------

        do ie=1,nel
	  lmax = ilhv(ie)
	  do l=1,lmax
            value = elv(l,ie)
            do ii=1,3
              k = nen3v(ii,ie)
	      if( mode .eq. 1 ) then
                nov(l,k) = max(nov(l,k),value)
	      else
                nov(l,k) = min(nov(l,k),value)
	      end if
	    end do
          end do
        end do

c-----------------------------------------------------------
c end of routine
c-----------------------------------------------------------

        end

c******************************************************************
c******************************************************************
c******************************************************************

        subroutine n2e2d(nov,elv)

c transforms nodal values to element values
c
c (2D version)

	use basin

        implicit none

        real nov(nkn)     !array with nodal values (in)
        real elv(nel)     !array with element values (out)

        integer k,ie,ii
        real acu,value

c-----------------------------------------------------------
c convert values
c-----------------------------------------------------------

        do ie=1,nel
          acu = 0.
          do ii=1,3
            k = nen3v(ii,ie)
            value = nov(k)
            acu = acu + value
          end do
          elv(ie) = acu / 3.
        end do

c-----------------------------------------------------------
c end of routine
c-----------------------------------------------------------

        end

c******************************************************************

        subroutine n2e3d(nlvddi,nov,elv)

c transforms nodal values to element values
c
c (3D version)

	use levels
	use basin

        implicit none

c arguments
        integer nlvddi		!vertical dimension of arrays
        real nov(nlvddi,nkn)	!array with nodal values (in)
        real elv(nlvddi,nel)	!array with element values (out)

c local
        integer k,ie,ii,l,lmax
        real acu,value

c-----------------------------------------------------------
c convert values
c-----------------------------------------------------------

        do ie=1,nel
          lmax = ilhv(ie)
          do l = 1,lmax
	    acu = 0.
            do ii=1,3
              k = nen3v(ii,ie)
              value = nov(l,k)
              acu = acu + value
	    end do
            elv(l,ie) = acu / 3.
          end do
	end do

c-----------------------------------------------------------
c end of routine
c-----------------------------------------------------------

        end

c******************************************************************

	subroutine aver_nodal(val,aver)

c computes average of val (defined on nodes) over total basin

	use evgeom
	use basin

	implicit none

        real val(nkn)     !array with element values (in)
	real aver	  !average (out)

        integer k,ie,ii
        double precision area,value
        double precision accum,area_tot

c-----------------------------------------------------------
c initialize arrays
c-----------------------------------------------------------

	accum = 0.
	area_tot = 0.

c-----------------------------------------------------------
c accumulate values
c-----------------------------------------------------------

        do ie=1,nel
          area = 4.*ev(10,ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
            value = val(k)
	    accum = accum + value * area
	    area_tot = area_tot + area
          end do
        end do

c-----------------------------------------------------------
c compute final value
c-----------------------------------------------------------

	if ( area_tot > 0. ) accum = accum / area_tot

	aver = accum

c-----------------------------------------------------------
c end of routine
c-----------------------------------------------------------

        end

c******************************************************************

