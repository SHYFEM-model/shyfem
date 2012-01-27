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
c subroutine n2e3d(nlvdim,nov,elv)      transforms nodal to element values 3D
c
c revision log :
c
c 10.08.2003    ggu     new routines copy_uvz, make_prvel, init_uvz
c 18.09.2003    ggu     new routine e2n2d
c 04.12.2003    ggu     new routine n2e2d
c 30.03.2004	ccf	new routine n2e2d and n2e3d
c 15.10.2004	ggu	new routine smagorinsky started
c 14.01.2005	ggu	bug (nlvdim) in n2e3d fixed
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
	implicit none
c
c parameters
	include 'param.h'
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nlvdi,nlv
	common /level/ nlvdi,nlv
	integer ilhv(1)
	common /ilhv/ilhv
	real utlnv(nlvdim,1),vtlnv(nlvdim,1)
	common /utlnv/utlnv, /vtlnv/vtlnv
	real ulnv(nlvdim,1),vlnv(nlvdim,1)
	common /ulnv/ulnv, /vlnv/vlnv
        real hdenv(nlvdim,1)
        common /hdenv/hdenv
c local
	integer ie,l,ilevel
	real h,rh

	if(nlvdim.ne.nlvdi) stop 'error stop : level dimension in ttov'

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

	implicit none

c parameters
	include 'param.h'
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nlvdi,nlv
	common /level/ nlvdi,nlv
	integer ilhv(1)
	common /ilhv/ilhv
	real utlnv(nlvdim,1),vtlnv(nlvdim,1)
	common /utlnv/utlnv, /vtlnv/vtlnv
	real ulnv(nlvdim,1),vlnv(nlvdim,1)
	common /ulnv/ulnv, /vlnv/vlnv
        real hdenv(nlvdim,1)
        common /hdenv/hdenv
c local
	integer ie,l,ilevel
	real h

	if(nlvdim.ne.nlvdi) stop 'error stop : level dimension in vtot'

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
	implicit none
c
c arguments
	real vv(1)
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nen3v(3,1)
	real up0v(1),vp0v(1)
	real unv(1),vnv(1)
	real znv(1),hev(1)
	common /nen3v/nen3v
	common /up0v/up0v, /vp0v/vp0v
	common /unv/unv, /vnv/vnv
	common /znv/znv, /hev/hev
	include 'ev.h'

        real zeov(3,1),zenv(3,1)
        common /zeov/zeov, /zenv/zenv
	integer iwegv(1)
	common /iwegv/iwegv
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
	implicit none
c
c parameters
	include 'param.h'
c arguments
	real vv(1)
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nlvdi,nlv
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /level/ nlvdi,nlv
	integer ilhv(1)
	integer nen3v(3,1)
	real ulnv(nlvdim,1),vlnv(nlvdim,1),wlnv(0:nlvdim,1)
	real uprv(nlvdim,1),vprv(nlvdim,1),wprv(0:nlvdim,1)
	common /ilhv/ilhv
	common /nen3v/nen3v
	common /ulnv/ulnv, /vlnv/vlnv, /wlnv/wlnv
	common /uprv/uprv, /vprv/vprv, /wprv/wprv
	include 'ev.h'

	integer iwegv(1)
	common /iwegv/iwegv
c local
	integer ie,l,k,ii
	real aj
c
	if(nlvdim.ne.nlvdi) stop 'error stop : level dimension in uvtopr'
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
	implicit none
c
c parameters
	include 'param.h'
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nlvdi,nlv
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /level/ nlvdi,nlv
	integer ilhv(1)
	integer nen3v(3,1)
	real ulnv(nlvdim,1),vlnv(nlvdim,1),wlnv(0:nlvdim,1)
	real uprv(nlvdim,1),vprv(nlvdim,1),wprv(0:nlvdim,1)
	common /ilhv/ilhv
	common /nen3v/nen3v
	common /ulnv/ulnv, /vlnv/vlnv, /wlnv/wlnv
	common /uprv/uprv, /vprv/vprv, /wprv/wprv
	integer iwegv(1)
	common /iwegv/iwegv
c local
	integer ie,l,k,ii
	real u,v
c
	if(nlvdim.ne.nlvdi) stop 'error stop : level dimension in prtouv'
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
	implicit none
c
c parameter
	include 'param.h'
c common
	integer nlvdi,nlv
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /level/ nlvdi,nlv
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer ilhv(1)
	real unv(1),vnv(1)
	real utlnv(nlvdim,1),vtlnv(nlvdim,1)
	common /ilhv/ilhv
	common /unv/unv, /vnv/vnv
	common /utlnv/utlnv, /vtlnv/vtlnv
c local
	integer ie,l
	real u,v
c
	if(nlvdim.ne.nlvdi) stop 'error stop : level dimension uvint'
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

	implicit none

c parameters
	include 'param.h'
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nlvdi,nlv
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /level/ nlvdi,nlv

        real znv(1)
        common /znv/znv
        real zenv(3,1)
        common /zenv/zenv
        real hdknv(nlvdim,1)
        common /hdknv/hdknv
        real hdenv(nlvdim,1)
        common /hdenv/hdenv

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

	implicit none

c parameters
	include 'param.h'
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nlvdi,nlv
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /level/ nlvdi,nlv
	integer ilhv(1)
	real utlnv(nlvdim,1),vtlnv(nlvdim,1)
	real ulnv(nlvdim,1),vlnv(nlvdim,1)
	common /utlnv/utlnv, /vtlnv/vtlnv
	common /ulnv/ulnv, /vlnv/vlnv
	common /ilhv/ilhv
	integer iwegv(1)
	common /iwegv/iwegv
	real unv(1),vnv(1)
	common /unv/unv, /vnv/vnv
        real zeov(3,1),zenv(3,1)
        common /zeov/zeov, /zenv/zenv
	real hm3v(3,1)
	common /hm3v/hm3v
c local
	logical bsigma
	integer nsigma
	integer ie,ilevel,ii,l
	real hsigma
c functions
        integer ieext

	if(nlvdim.ne.nlvdi) stop 'error stop : level dimension in ttov'

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

	implicit none

c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real up0v(1),vp0v(1)
	common /up0v/up0v, /vp0v/vp0v
	real znv(1)
	common /znv/znv
	real xv(3,1)
	common /xv/xv
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

	implicit none

        integer l,k
        real u,v

c parameters
	include 'param.h'
c common
	real uprv(nlvdim,1),vprv(nlvdim,1),wprv(0:nlvdim,1)
	common /uprv/uprv, /vprv/vprv, /wprv/wprv

        u = uprv(l,k)
        v = vprv(l,k)

        end

c******************************************************************

	subroutine copy_uvz

c copies u/v/z to old time step

	implicit none

	include 'param.h'

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nlvdi,nlv
	common /level/ nlvdi,nlv

	real utlov(nlvdim,1), vtlov(nlvdim,1)
	common /utlov/utlov, /vtlov/vtlov
	real utlnv(nlvdim,1), vtlnv(nlvdim,1)
	common /utlnv/utlnv, /vtlnv/vtlnv

	real ulov(nlvdim,1), vlov(nlvdim,1)
	common /ulov/ulov,  /vlov/vlov
	real ulnv(nlvdim,1), vlnv(nlvdim,1)
	common /ulnv/ulnv,  /vlnv/vlnv

        real uprv(nlvdim,1)
        common /uprv/uprv
        real vprv(nlvdim,1)
        common /vprv/vprv

        real upro(nlvdim,1)
        common /upro/upro
        real vpro(nlvdim,1)
        common /vpro/vpro

	real wlov(0:nlvdim,1)
	common /wlov/wlov
	real wlnv(0:nlvdim,1)
	common /wlnv/wlnv

	real uov(1),vov(1)
	common /uov/uov, /vov/vov
	real unv(1),vnv(1)
	common /unv/unv, /vnv/vnv

	real zov(1)
	common /zov/zov
	real znv(1)
	common /znv/znv
        real zeov(3,1),zenv(3,1)
        common /zeov/zeov, /zenv/zenv

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

	implicit none

	include 'param.h'

	real v1v(1)
	common /v1v/v1v

	call uvtopr(v1v)
	call uvtop0(v1v)
	call setxv

	end

c******************************************************************

	subroutine init_uv

c initializes uvz values from zenv, utlnv, vtlnv, hdenv

	implicit none

	include 'param.h'

	real saux1(nlvdim,1)
	common /saux1/saux1
	real saux2(nlvdim,1)
	common /saux2/saux2

	logical has_restart

	call ttov			!velocities from layer transports
	call uvint			!barotropic transports

	if( .not. has_restart(5) ) then
	  call sp256w(saux1,saux2)	!vertical velocities
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

	implicit none

c arguments
        real elv(1)     !array with element values (in)
        real nov(1)     !array with nodal values (out)
        real aux(1)     !aux array (nkndim)

c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nen3v(3,1)
	common /nen3v/nen3v
	include 'ev.h'

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

	subroutine e2n3d(nlvdim,elv,nov,aux)

c transforms element values to nodal values (weights are area)
c
c (3D version)

	implicit none

c arguments
	integer nlvdim
        real elv(nlvdim,1)     !array with element values (in)
        real nov(nlvdim,1)     !array with nodal values (out)
        real aux(nlvdim,1)     !aux array (nkndim)

c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nen3v(3,1)
	common /nen3v/nen3v
        integer ilhv(1),ilhkv(1)
        common /ilhv/ilhv, /ilhkv/ilhkv
	include 'ev.h'

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

	subroutine e2n3d_minmax(mode,nlvdim,elv,nov)

c transforms element values to nodal values (no weights - use min/max)
c
c (3D version)

	implicit none

c arguments
	integer mode		!min (-1) or max (+1)
	integer nlvdim		!vertical dimension
        real elv(nlvdim,1)      !array with element values (in)
        real nov(nlvdim,1)      !array with nodal values (out)

c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nen3v(3,1)
	common /nen3v/nen3v
        integer ilhv(1),ilhkv(1)
        common /ilhv/ilhv, /ilhkv/ilhkv

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

        implicit none

        real nov(1)     !array with nodal values (in)
        real elv(1)     !array with element values (out)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        integer nen3v(3,1)
        common /nen3v/nen3v

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

        subroutine n2e3d(nlvdim,nov,elv)

c transforms nodal values to element values
c
c (3D version)

        implicit none

c arguments
        integer nlvdim		!vertical dimension of arrays
        real nov(nlvdim,1)	!array with nodal values (in)
        real elv(nlvdim,1)	!array with element values (out)

c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer ilhv(1),ilhkv(1)
        common /ilhv/ilhv, /ilhkv/ilhkv
        integer nen3v(3,1)
        common /nen3v/nen3v

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

