!
! $Id: newtra.f,v 1.22 2010-03-11 15:36:38 georg Exp $
!
! routines for various transformations
!
! contents :
!
! subroutine ttov			transforms transports to velocities
! subroutine vtot			transforms velocities to transports
! subroutine uvtop0			transforms bar. transp. to nod. veloc.
! subroutine uvtopr			transforms velocities to nodal values
! subroutine prtouv			transforms nodal values to element vel.
! subroutine uvint			computation of barotropic transports
! subroutine austau(vv)			computes aux vectors for austausch term
! subroutine baro2l			distribute barotropic velocities 
! subroutine setxv			sets obsolete data structure xv
! subroutine getuv(l,k,u,v)             accessor routine to get velocities u/v
! subroutine copy_uvz                   copies u/v/z to old time step
! subroutine make_prvel                 makes print velocities and xv
! subroutine e2n2d(elv,nov,aux)         transforms element to nodal values
! subroutine n2e2d(nov,elv)             transforms nodal to element values 2D
! subroutine n2e3d(nlvdi,nov,elv)       transforms nodal to element values 3D
!
! revision log :
!
! 10.08.2003    ggu     new routines copy_uvz, make_prvel, init_uvz
! 18.09.2003    ggu     new routine e2n2d
! 04.12.2003    ggu     new routine n2e2d
! 30.03.2004	ccf	new routine n2e2d and n2e3d
! 15.10.2004	ggu	new routine smagorinsky started
! 14.01.2005	ggu	bug (nlvdi) in n2e3d fixed
! 24.02.2005	ggu	smagorinsky into subdif.f
! 15.03.2005	ggu	austv eliminated, austau() into subdif.f
! 27.06.2005	ggu	bug in vtot corrected
! 10.04.2008	ggu	copy velocities at nodes in copy_uvz()
! 01.03.2010	ggu	new version of n2e3d()
! 11.03.2010	ggu	new routine check_volume(); init w only if no restart
! 16.02.2011	ggu	new routine e2n3d() and e2n3d_minmax()
! 27.01.2012	deb&ggu	routines adapted for sigma levels
! 03.12.2015	ccf&ggu	code optimized
!
!****************************************************************************
!----------------------------------------------------------------------------
        module transforms
!----------------------------------------------------------------------------
        contains
!----------------------------------------------------------------------------

	subroutine ttov
!
! transforms transports to velocities
!
	use layer_thickness
	use hydro_vel
	use hydro_admin
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer ie,l,ilevel
	double precision h,rh

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

!******************************************************************

	subroutine vtot

! transforms velocities to transports

	use layer_thickness
	use hydro_vel
	use hydro_admin
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer ie,l,ilevel
	double precision h

	do ie=1,nel

	  ilevel=ilhv(ie)

	  do l=1,ilevel
	    h = hdenv(l,ie)
	    utlnv(l,ie)=ulnv(l,ie)*h
	    vtlnv(l,ie)=vlnv(l,ie)*h
	  end do

	end do

	end
!
!******************************************************************
!
	subroutine uvtop0
!
! transforms barotropic transports to nodal velocities
!
	use geom_dynamic
	use depth
	use hydro_baro
	use hydro_print
	use hydro_admin
	use evgeom
	use basin
	use shympi
	use area
        use utility
        use para
        use mpi_io_admin

	implicit none

! local
	logical bcolin
	integer ie,k,ii,e
	double precision aj,zm,hm
!
	bcolin=iround(getpar('iclin')).ne.0
!
	up0v = 0.d0
	vp0v = 0.d0

#ifdef DEBUGON
        call shympi_exchange_halo_2d_elems(iwegv)
        call shympi_exchange_halo_3d_elems3(zenv)
        call shympi_exchange_halo_2d_elems(hev)
        call shympi_exchange_halo_2d_elems(unv)
        call shympi_exchange_halo_2d_elems(vnv)
!
	do e=1,nel_local
          if(bmpi) then
            ie=domain%elems%mapID(e)
          else
            ie=e
          end if
#else
	do ie=1,nel
#endif
	  if( iwegv(ie) /= 0 ) cycle
	  aj=4.*ev(10,ie)
	  zm=0.
	  do ii=1,3
	    zm=zm+zenv(ii,ie)
	  end do
	  zm=zm/3.
	  hm=hev(ie)
	  if(.not.bcolin) hm=hm+zm
	  do ii=1,3
	    k=nen3v(ii,ie)
	    up0v(k)=up0v(k)+aj*unv(ie)/hm
	    vp0v(k)=vp0v(k)+aj*vnv(ie)/hm
	  end do
	end do

        !       use only areakv(1,k) for weighting (surface layer)
        !call shympi_comment('shympi_elem: exchange up0v, vp0v')
#ifndef DEBUGON
        call shympi_exchange_and_sum_2D_nodes(up0v)
        call shympi_exchange_and_sum_2D_nodes(vp0v)
#endif

	where (areakv(1,:) > 0.) 
          up0v = up0v / areakv(1,:)
          vp0v = vp0v / areakv(1,:)
	end where

	return
	end
!
!******************************************************************
!
        subroutine uvtopr
!
! transforms velocities to nodal values
!
        use geom_dynamic
        use hydro_print
        use hydro_vel
        use evgeom
        use levels
        use basin
        use shympi
        use area

        implicit none

        include 'femtime.h'

        integer ie,l,k,ii,e
        integer lmax
        double precision aj

        uprv = 0.d0
        vprv = 0.d0
!
! baroclinic part
!

#ifdef DEBUGON
        call shympi_exchange_halo_2d_elems(ilhv)
        call shympi_exchange_halo_2d_elems(iwegv)
        call shympi_exchange_halo_3d_elems(ulnv)
        call shympi_exchange_halo_3d_elems(vlnv)

        do e=1,nel_local
          if(bmpi) then
            ie=domain%elems%mapID(e)
          else
            ie=e
          end if
#else
        do ie=1,nel
#endif
          if ( iwegv(ie) /= 0 ) cycle
          lmax = ilhv(ie)
          aj=4.*ev(10,ie)
          do l=1,lmax
            do ii=1,3
              k=nen3v(ii,ie)
              uprv(l,k)=uprv(l,k)+aj*ulnv(l,ie)
              vprv(l,k)=vprv(l,k)+aj*vlnv(l,ie)
            end do
          end do
        end do

        !call shympi_comment('shympi_elem: exchange uprv, vprv')
#ifndef DEBUGON
        call shympi_exchange_and_sum_3D_nodes(uprv)
        call shympi_exchange_and_sum_3D_nodes(vprv)
#endif

        where ( areakv > 0. )
          uprv = uprv / areakv
          vprv = vprv / areakv
        end where

! vertical velocities -> we compute average over one layer
!
        do l=1,nlv
          wprv(l,:)=0.5*(wlnv(l,:)+wlnv(l-1,:))
        end do

        wprv(0,:) = 0.
!
        return
        end
!
!******************************************************************
!
	subroutine prtouv
!
! transforms nodal values to element values (velocities)
!
	use geom_dynamic
	use hydro_print
	use hydro_vel
	use levels
	use basin

	implicit none

	integer ie,l,k,ii
	double precision u,v
!
! baroclinic part
!
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
!
! vertical velocities -> from layer average to interface values
!
	wlnv(nlv,:) = 0.
	do l=nlv-1,0,-1
	  wlnv(l,:)=2.*wprv(l,:)-wlnv(l+1,:)
	end do
!
	return
	end
!
!*****************************************************************
!
	subroutine uvint
!
! computation of barotropic part of transports
!
	use hydro_baro
	use hydro_admin
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer ie,l
	double precision u,v
!
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
!
	return
	end
!
!*****************************************************************

	subroutine check_volume

! checks for negative volume (depth)
!
! only first layer has to be checked

	use layer_thickness
	use hydro_admin
        use fem_util
        use sigma_admin
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	logical bstop,bsigma
	integer nsigma
	integer k,ie,ke,iee,ii
	double precision h,z
	double precision hsigma
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

!*****************************************************************
!
	subroutine baro2l

! distribute barotropic velocities onto layers (only in dry elements)

	use geom_dynamic
	use hydro_baro
	use hydro_vel
	use hydro_admin
	use levels
	use basin
        use fem_util
        use sigma_admin

	implicit none

	logical bsigma
	integer nsigma
	integer ie,ilevel,ii,l
	double precision hsigma
! functions

	call get_sigma(nsigma,hsigma)
	bsigma = nsigma .gt. 0

	if( bsigma ) return

	do ie=1,nel
	  if( iwegv(ie) .gt. 0 ) then	!dry
	    ilevel=ilhv(ie)
	    if( ilevel .gt. 1 ) then
		write(6,*) 'drying in more than one layer'
		write(6,*) 'ie,iext,ilevel:',ie,ieext(ie),ilevel
	        write(6,*) 'zeov(ii,ie),zenv(ii,ie),hm3v(ii,ie)'
		do ii=1,3
		  write(6,*) zeov(ii,ie),zenv(ii,ie),hm3v(ii,ie)
		end do
                write(6,*) 'utlnv(l,ie),vtlnv(l,ie)'
                do l=1,ilevel
                  write(6,*) utlnv(l,ie),vtlnv(l,ie)
                end do
                write(6,*) 'unv(ie),vnv(ie)'
                write(6,*) unv(ie),vnv(ie)
		stop 'error stop baro2l: drying in more than one layer'
	    end if
	    utlnv(1,ie) = unv(ie)
	    vtlnv(1,ie) = vnv(ie)
	  end if
	end do

	end

!******************************************************************

	subroutine setxv

! sets obsolete data structure xv

	use hydro_print
	use hydro_admin
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	xv(1,:) = up0v(:)
	xv(2,:) = vp0v(:)
	xv(3,:) = znv(:)

	end

!******************************************************************

	subroutine copy_uvz

! copies u/v/z to old time step

	use hydro_baro
	use hydro_print
	use hydro_vel
	use hydro_admin

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

!******************************************************************

        subroutine make_prvel

! makes print velocities and xv from new level arrays

        use basin
        use hydro_print
        use shympi

        implicit none

        call uvtopr     !computes uprv,vprv
        call uvtop0     !computes up0v,vp0v
        call setxv      !sets xv from up0v,vp0v,znv

        !call shympi_comment('exchanging uprv, vprv, up0v, vp0v')
        call shympi_exchange_3d_node(uprv)
        call shympi_exchange_3d_node(vprv)
        call shympi_exchange_2d_node(up0v)
        call shympi_exchange_2d_node(vp0v)
        !call shympi_barrier

        end

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine e2n2d(elv,nov,aux)

! transforms element values to nodal values (weights are area)
!
! (2D version)

	use evgeom
	use basin

	implicit none

! arguments
        double precision elv(nel)     !array with element values (in)
        double precision nov(nkn)     !array with nodal values (out)
        double precision aux(nkn)     !aux array (nkn)

! local
        integer k,ie,ii
        double precision area,value

!-----------------------------------------------------------
! initialize arrays
!-----------------------------------------------------------

        nov = 0.
        aux = 0.

!-----------------------------------------------------------
! accumulate values
!-----------------------------------------------------------

        do ie=1,nel
          area = 4.*ev(10,ie)
          value = elv(ie)
          do ii=1,3
            k = nen3v(ii,ie)
            nov(k) = nov(k) + area*value
            aux(k) = aux(k) + area
          end do
        end do

!-----------------------------------------------------------
! compute final value
!-----------------------------------------------------------

        nov = nov / aux

!-----------------------------------------------------------
! end of routine
!-----------------------------------------------------------

        end

!******************************************************************

	subroutine e2n3d(nlvddi,elv,nov,aux)

! transforms element values to nodal values (weights are area)
!
! (3D version)

	use evgeom
	use levels
	use basin

	implicit none

! arguments
	integer nlvddi
        double precision elv(nlvddi,nel)     !array with element values (in)
        double precision nov(nlvddi,nkn)     !array with nodal values (out)
        double precision aux(nlvddi,nkn)     !aux array (nkn)

! local
        integer k,ie,ii,l,lmax
        double precision area,value

!-----------------------------------------------------------
! initialize arrays
!-----------------------------------------------------------

        nov = 0.
        aux = 0.

!-----------------------------------------------------------
! accumulate values
!-----------------------------------------------------------

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

!-----------------------------------------------------------
! compute final value
!-----------------------------------------------------------

	where ( aux > 0. ) 
	  nov = nov / aux
	end where

!-----------------------------------------------------------
! end of routine
!-----------------------------------------------------------

        end

!******************************************************************

	subroutine e2n3d_minmax(mode,nlvddi,elv,nov)

! transforms element values to nodal values (no weights - use min/max)
!
! (3D version)

	use levels
	use basin

	implicit none

! arguments
	integer mode		!min (-1) or max (+1)
	integer nlvddi		!vertical dimension
        double precision elv(nlvddi,nel)      !array with element values (in)
        double precision nov(nlvddi,nkn)      !array with nodal values (out)

! local
        integer k,ie,ii,l,lmax
        double precision rinit,value

!-----------------------------------------------------------
! initialize arrays
!-----------------------------------------------------------

	if( mode .eq. -1) then
	  rinit = 1.e+30
	else if( mode .eq. 1 ) then
	  rinit = -1.e+30
	else
	  stop 'error stop e2n3d_minmax: unknown mode'
	end if

        nov = rinit

!-----------------------------------------------------------
! accumulate values
!-----------------------------------------------------------

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

!-----------------------------------------------------------
! end of routine
!-----------------------------------------------------------

        end

!******************************************************************
!******************************************************************
!******************************************************************

        subroutine n2e2d(nov,elv)

! transforms nodal values to element values
!
! (2D version)

	use basin

        implicit none

        double precision nov(nkn)     !array with nodal values (in)
        double precision elv(nel)     !array with element values (out)

        integer k,ie,ii
        double precision acu,value

!-----------------------------------------------------------
! convert values
!-----------------------------------------------------------

        do ie=1,nel
          acu = 0.
          do ii=1,3
            k = nen3v(ii,ie)
            value = nov(k)
            acu = acu + value
          end do
          elv(ie) = acu / 3.
        end do

!-----------------------------------------------------------
! end of routine
!-----------------------------------------------------------

        end

!******************************************************************

        subroutine n2e3d(nlvddi,nov,elv)

! transforms nodal values to element values
!
! (3D version)

	use levels
	use basin

        implicit none

! arguments
        integer nlvddi		!vertical dimension of arrays
        double precision nov(nlvddi,nkn)	!array with nodal values (in)
        double precision elv(nlvddi,nel)	!array with element values (out)

! local
        integer k,ie,ii,l,lmax
        double precision acu,value

!-----------------------------------------------------------
! convert values
!-----------------------------------------------------------

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

!-----------------------------------------------------------
! end of routine
!-----------------------------------------------------------

        end

!******************************************************************

!----------------------------------------------------------------------------
        end module transforms
!----------------------------------------------------------------------------
