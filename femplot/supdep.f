c
c $Id: supdep.f,v 1.5 2009-11-18 17:16:00 georg Exp $
c
c routines for averaging depth
c
c contents :
c
c subroutine mkht(hetv,href)		makes hetv (elem depth of actual layer)
c subroutine mkht3(nlvddi,het3v,href)	makes het3v (3D depth structure)
c function hlthick(l,lmax,hl)		layer thickness
c
c revision log :
c
c 26.05.2000    ggu     routines written from scratch
c 17.09.2008    ggu     routine mkht changed for layer = -1
c 13.10.2009    ggu     new routine mkht3
c 13.10.2009    ggu     new routine mkht3
c 17.12.2010    ggu     substituted hv with hkv, new routine hlthick()
c 30.03.2011    ggu     new routine mkareafvl()
c 14.11.2011    ggu     use get_layer_thickness() for layer structure
c 23.02.2012    ccf     bug fix in mkht (get_layer_thicknes, get_sigma_info)
c 16.03.2012    deb     bug fix in mkht3 (get_layer_thicknes, get_sigma_info)
c 10.06.2013    ggu     bug fix in mkht,mkht3 (get_layer_thicknes_e)
c 05.09.2013    ggu     adapt to new get_layer_thickness()
c 27.05.2016    ggu     mkhkv,mkhev deleted
c
c******************************************************************

	subroutine mkht(hetv,href)

c makes hetv (elementwise depth of actual layer)
c
c uses level to decide what to do

	use mod_depth
	use mod_hydro
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	real hetv(nel)
	real href

	logical bdebug
	integer ie,ii
	integer level,lmax
	real zeta
	integer nsigma,nlvaux
	real hsigma
        real hl(nlv)

	integer getlev
	real hlthick

c-------------------------------------------------------------------
c initialization
c-------------------------------------------------------------------

        bdebug = .true.
        bdebug = .false.

	level = getlev()
        call get_sigma_info(nlvaux,nsigma,hsigma)

c-------------------------------------------------------------------
c handle different kind of levels
c-------------------------------------------------------------------

	do ie=1,nel
	  lmax = ilhv(ie)
	  call compute_levels_on_element(ie,zenv,zeta)
	  call get_layer_thickness(lmax,nsigma,hsigma,zeta,hev(ie),hlv,hl)
	  hetv(ie) = hlthick(level,lmax,hl)
	end do

c-------------------------------------------------------------------
c debug output
c-------------------------------------------------------------------

	if( bdebug ) then
	  ie = 1644
	  write(6,*) 'debugging mkht...'
	  write(6,*) ilhv(ie),hetv(ie),hev(ie)
	  do ii=1,ilhv(ie)
	    write(6,*) hlv(ii)
	  end do
	end if

c-------------------------------------------------------------------
c end of routine
c-------------------------------------------------------------------

	end

c******************************************************************

	subroutine mkht3(nlvddi,het3v,href)

c makes het3v (3D depth structure)

	use mod_depth
	use mod_hydro
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

c arguments
	integer nlvddi
	real het3v(nlvddi,1)
	real href
c local
	logical bdebug
	integer ie,ii
	integer l,lmax
	integer nlvaux,nsigma
	real hsigma
	real zeta

c-------------------------------------------------------------------
c initialization
c-------------------------------------------------------------------

        bdebug = .true.
        bdebug = .false.

        call get_sigma_info(nlvaux,nsigma,hsigma)

c-------------------------------------------------------------------
c compute layer thickness
c-------------------------------------------------------------------

	do ie=1,nel
	  lmax = ilhv(ie)
	  call compute_levels_on_element(ie,zenv,zeta)
	  call get_layer_thickness(lmax,nsigma,hsigma,zeta,hev(ie),hlv
     +				,het3v(1,ie))
!	  call get_layer_thickness_e(ie,lmax,bzeta,nsigma,hsigma
!     +				,het3v(1,ie))
	end do

c-------------------------------------------------------------------
c end of routine
c-------------------------------------------------------------------

	end

c******************************************************************

	function hlthick(l,lmax,hl)

c computes thickness of layer l
c
c works also for sigma layers

	implicit none

	real hlthick		! layer thickness (return)
	integer l		! layer to compute thickness
	integer lmax		! maximum layers
	real hl(1)		! level thickness

	integer ll
	real hm

	if( l .eq. 0 ) then			! total layer
	  hm = 0.
	  do ll=1,lmax
	    hm = hm + hl(ll)
	  end do
	  hlthick = hm
	else if( l .ge. 1 .and. l .le. lmax ) then
	  hlthick = hl(l)
	else if( l .eq. -1 ) then		! bottom layer
	  hlthick = hl(lmax)
	else
	  hlthick = 0.
	end if

	end

c******************************************************************

	subroutine mkareafvl

c makes area of finite volume

	use mod_hydro_plot
	use evgeom
	use basin

	implicit none

	integer ie,ii,k
	real afvl

	do k=1,nkn
	  arfvlv(k) = 0.
	end do

	do ie=1,nel
	  afvl = 4. * ev(10,ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    arfvlv(k) = arfvlv(k) + afvl
	  end do
	end do

	end

c******************************************************************

