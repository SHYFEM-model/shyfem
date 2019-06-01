
!--------------------------------------------------------------------------
!
!    Copyright (C) 2000,2008-2011,2013,2016  Georg Umgiesser
!    Copyright (C) 2012  Christian Ferrarin
!    Copyright (C) 2012  Debora Bellafiore
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
c 26.05.2000	ggu	routines written from scratch
c 17.09.2008	ggu	routine mkht changed for layer = -1
c 13.10.2009	ggu	new routine mkht3
c 13.10.2009	ggu	new routine mkht3
c 23.03.2010	ggu	changed v6.1.1
c 17.12.2010	ggu	substituted hv with hkv, new routine hlthick()
c 30.03.2011	ggu	new routine mkareafvl()
c 14.04.2011	ggu	changed VERS_6_1_22
c 14.11.2011	ggu	use get_layer_thickness() for layer structure
c 22.11.2011	ggu	changed VERS_6_1_37
c 23.02.2012	ccf	bug fix in mkht (get_layer_thicknes, get_sigma_info)
c 16.03.2012	dbf	bug fix in mkht3 (get_layer_thicknes, get_sigma_info)
c 10.06.2013	ggu	bug fix in mkht,mkht3 (get_layer_thicknes_e)
c 05.09.2013	ggu	adapt to new get_layer_thickness()
c 12.09.2013	ggu	changed VERS_6_1_67
c 23.12.2014	ggu	changed VERS_7_0_11
c 19.01.2015	ggu	changed VERS_7_1_2
c 19.01.2015	ggu	changed VERS_7_1_3
c 05.05.2015	ggu	changed VERS_7_1_10
c 10.07.2015	ggu	changed VERS_7_1_50
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 25.05.2016	ggu	changed VERS_7_5_10
c 27.05.2016	ggu	mkhkv,mkhev deleted
c 12.01.2017	ggu	changed VERS_7_5_21
c 25.10.2018	ggu	changed VERS_7_5_51
c 18.12.2018	ggu	changed VERS_7_5_52
c 13.03.2019	ggu	changed VERS_7_5_61
c 21.05.2019	ggu	changed VERS_7_5_62
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
	real het3v(nlvddi,nel)
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
	real hl(lmax)		! level thickness

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

