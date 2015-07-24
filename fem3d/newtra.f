c
c $Id: newtra.f,v 1.22 2010-03-11 15:36:38 georg Exp $
c
c routines for various transformations
c
c contents :
c
c subroutine ttov			transforms transports to velocities
c subroutine vtot			transforms velocities to transports
c subroutine uvtop0(vv)			transforms bar. transp. to nod. veloc.
c subroutine uvtopr(vv)			transforms velocities to nodal values
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
c
c****************************************************************************
c
	subroutine ttov
c
c transforms transports to velocities
c
	use mod_layer_thickness
	use mod_hydro_vel
	use mod_hydro
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none
c
c parameters
	include 'param.h'
c common
c local
	integer ie,l,ilevel
	real h,rh

	do ie=1,nel

	  ilevel=ilhv(ie)

	  do l=1,ilevel
	    h = hdenv(l,ie)
	    rh = 1. / h
	    ulnv(l,ie)=utlnv(l,ie)*rh
	    vlnv(l,ie)=vtlnv(l,ie)*rh
	  end do

	end do

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

c parameters
	include 'param.h'
c common
c local
	integer ie,l,ilevel
	real h

	do ie=1,nel

	  ilevel=ilhv(ie)

	  do l=1,ilevel
	    h = hdenv(l,ie)
	    utlnv(l,ie)=ulnv(l,ie)*h
	    vtlnv(l,ie)=vlnv(l,ie)*h
	  end do

	end do

	end
c
c******************************************************************
c
	subroutine uvtop0(vv)
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
c
c arguments
	real vv(1)
c common
	include 'param.h'

c local
	logical bcolin
	integer ie,k,ii
	real aj,zm,hm
c function
	real getpar
	integer iround
c
	bcolin=iround(getpar('iclin')).ne.0
c
	do k=1,nkn
	  up0v(k)=0.
	  vp0v(k)=0.
	  vv(k)=0.
	end do
c
	do ie=1,nel
	 if( iwegv(ie) .eq. 0 ) then
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
	 end if
	end do
c
	do k=1,nkn
	  if(vv(k).gt.0.) then
	    up0v(k)=up0v(k)/vv(k)
	    vp0v(k)=vp0v(k)/vv(k)
	  end if
	end do
c
	return
	end
c
c******************************************************************
c
	subroutine uvtopr(vv)
c
c transforms velocities to nodal values
c
	use mod_geom_dynamic
	use mod_hydro_print
	use mod_hydro_vel
	use evgeom
	use levels
	use basin

	implicit none
c
c parameters
	include 'param.h'
c arguments
	real vv(1)
c common

c local
	integer ie,l,k,ii
	real aj
c
	do k=1,nkn
	  do l=1,nlv
	    uprv(l,k)=0.
	    vprv(l,k)=0.
	  end do
	end do
c
c baroclinic part
c
	do l=1,nlv
c
	  do k=1,nkn
	    vv(k)=0.
	  end do
c
	  do ie=1,nel
	   if( iwegv(ie) .eq. 0 ) then
	    if(l.le.ilhv(ie)) then
	      aj=ev(10,ie)
	      do ii=1,3
	        k=nen3v(ii,ie)
	        vv(k)=vv(k)+aj
	        uprv(l,k)=uprv(l,k)+aj*ulnv(l,ie)
	        vprv(l,k)=vprv(l,k)+aj*vlnv(l,ie)
	      end do
	    end if
	   end if
	  end do
c
	  do k=1,nkn
	    if(vv(k).gt.0.) then
	      uprv(l,k)=uprv(l,k)/vv(k)
	      vprv(l,k)=vprv(l,k)/vv(k)
	    end if
	  end do
c
	end do
c
c vertical velocities -> we compute average over one layer
c
	do k=1,nkn
	  do l=nlv,1,-1
	    wprv(l,k)=0.5*(wlnv(l,k)+wlnv(l-1,k))
	  end do
	  wprv(0,k)=0.
	end do
c
	return
	end
c
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
c
c parameters
	include 'param.h'
c common
c local
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
	  do l=1,ilhv(ie)
	    ulnv(l,ie)=0.
	    vlnv(l,ie)=0.
	  end do
	 end if
	end do
c
c vertical velocities -> from layer average to interface values
c
	do k=1,nkn
	  wlnv(nlv,k)=0.
	  do l=nlv-1,0,-1
	    wlnv(l,k)=2.*wprv(l,k)-wlnv(l+1,k)
	  end do
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
c
c parameter
	include 'param.h'
c common
c local
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

c parameters
	include 'param.h'
c common


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

c parameters
	include 'param.h'
c common
c local
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

c common
	include 'param.h'
c local
	integer k

	do k=1,nkn
	  xv(1,k) = up0v(k)
	  xv(2,k) = vp0v(k)
	  xv(3,k) = znv(k)
	end do

	end

c******************************************************************

	subroutine getuv(l,k,u,v)

c accessor routine to get velocities u/v

	use mod_hydro_print

	implicit none

        integer l,k
        real u,v

c parameters
	include 'param.h'
c common

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
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'









	integer k,ie,ii,l

	do k=1,nkn
	  zov(k)=znv(k)
          do l = 1,nlv
            upro(l,k) = uprv(l,k)
            vpro(l,k) = vprv(l,k)
	  end do
	end do

	do ie=1,nel				!ZEONV
	  do ii=1,3
	    zeov(ii,ie) = zenv(ii,ie)
	  end do
	end do

	do ie=1,nel
	  uov(ie)=unv(ie)			!$$UVBARO
	  vov(ie)=vnv(ie)
	end do

	do ie=1,nel
	 do l=1,nlv
	  ulov(l,ie)=ulnv(l,ie)
	  vlov(l,ie)=vlnv(l,ie)
	  utlov(l,ie)=utlnv(l,ie)
	  vtlov(l,ie)=vtlnv(l,ie)
	 end do
	end do

	do k=1,nkn
	 do l=0,nlv
	  wlov(l,k)=wlnv(l,k)
	 end do
	end do

	end

c******************************************************************

	subroutine make_prvel

c makes print velocities and xv from new level arrays

	use mod_aux_array

	implicit none

	include 'param.h'


	call uvtopr(v1v)
	call uvtop0(v1v)
	call setxv

	end

c******************************************************************

	subroutine init_uv

c initializes uvz values from zenv, utlnv, vtlnv, hdenv

	use mod_aux_array

	implicit none

	include 'param.h'


	logical has_restart

	call ttov			!velocities from layer transports
	call uvint			!barotropic transports

	if( .not. has_restart(5) ) then
	  call sp256w(v1v,saux1,saux2)	!vertical velocities
	end if

	call make_prvel		!nodal values

	call copy_uvz		!copy to old time step

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
        real elv(1)     !array with element values (in)
        real nov(1)     !array with nodal values (out)
        real aux(1)     !aux array (nkn)

c common
	include 'param.h'

c local
        integer k,ie,ii
        real area,value

c-----------------------------------------------------------
c initialize arrays
c-----------------------------------------------------------

        do k=1,nkn
          nov(k) = 0.
          aux(k) = 0.
        end do

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

        do k=1,nkn
          if( aux(k) .gt. 0. ) then
            nov(k) = nov(k) / aux(k)
          end if
        end do

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
        real elv(nlvddi,1)     !array with element values (in)
        real nov(nlvddi,1)     !array with nodal values (out)
        real aux(nlvddi,1)     !aux array (nkn)

c common
	include 'param.h'

c local
        integer k,ie,ii,l,lmax
        real area,value

c-----------------------------------------------------------
c initialize arrays
c-----------------------------------------------------------

        do k=1,nkn
	  lmax = ilhkv(k)
	  do l=1,lmax
            nov(l,k) = 0.
            aux(l,k) = 0.
	  end do
        end do

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

        do k=1,nkn
	  lmax = ilhkv(k)
	  do l=1,lmax
            if( aux(l,k) .gt. 0. ) then
              nov(l,k) = nov(l,k) / aux(l,k)
            end if
	  end do
        end do

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
        real elv(nlvddi,1)      !array with element values (in)
        real nov(nlvddi,1)      !array with nodal values (out)

c common
	include 'param.h'

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

        do k=1,nkn
	  lmax = ilhkv(k)
	  do l=1,lmax
            nov(l,k) = rinit
	  end do
        end do

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

        real nov(1)     !array with nodal values (in)
        real elv(1)     !array with element values (out)


	include 'param.h'

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
        real nov(nlvddi,1)	!array with nodal values (in)
        real elv(nlvddi,1)	!array with element values (out)

c common
	include 'param.h'

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

