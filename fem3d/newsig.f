
!--------------------------------------------------------------------------
!
!    Copyright (C) 2010-2013,2015,2019  Georg Umgiesser
!    Copyright (C) 2012  Debora Bellafiore
!    Copyright (C) 2012  Christian Ferrarin
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

c routines for sigma levels
c
c revision log :
c
c 23.03.2010	ggu	changed v6.1.1
c 15.12.2010	ggu	changed VERS_6_1_14
c 16.12.2010	ggu	program partially finished
c 23.03.2011	ggu	changed VERS_6_1_21
c 19.09.2011	ggu	new routine set_bsigma()
c 18.10.2011	ggu	changed VERS_6_1_33
c 04.11.2011	ggu	new routines for hybrid levels
c 10.11.2011	ggu	adjust depth for hybrid levels
c 11.11.2011	ggu	error check in set_hkv_and_hev()
c 11.11.2011	ggu	in check_hsigma_crossing set zeta levels to const depth
c 18.11.2011	ggu	restructured hybrid - adjustment to bashsigma
c 09.12.2011	ggu	changed VERS_6_1_38
c 12.12.2011	ggu	eliminated (stupid) compiler bug (getpar)
c 27.01.2012	dbf&ggu	adapted for hybrid levels
c 23.02.2012	ccf	bug fix in set_hybrid_depth (no call to get_sigma)
c 30.03.2012	ggu	changed VERS_6_1_51
c 05.11.2012	ggu	changed VERS_6_1_60
c 05.09.2013	ggu	no set_sigma_hkv_and_hev()
c 12.09.2013	ggu	changed VERS_6_1_67
c 19.01.2015	ggu	changed VERS_7_1_3
c 05.05.2015	ggu	changed VERS_7_1_10
c 05.06.2015	ggu	changed VERS_7_1_12
c 16.02.2019	ggu	changed VERS_7_5_60
c 13.03.2019	ggu	changed VERS_7_5_61
c
c notes :
c
c important files where sigma levels are explicitly needed:
c
c	newini.f		set up of structure
c	subele.f		set new layer thickness
c
c	newbcl.f		for computation of rho
c	newexpl.f		for baroclinic term
c
c	lagrange_flux.f		limit zeta layers to surface layer
c
c********************************************************************
c********************************************************************
c********************************************************************

	subroutine get_bsigma(bsigma)

c returns bsigma which is true if sigma layers are used

	implicit none

	logical bsigma

	real getpar

	bsigma = nint(getpar('nsigma')) .gt. 0

	end

c********************************************************************

	subroutine get_sigma(nsigma,hsigma)

	implicit none

	integer nsigma
	real hsigma

	real getpar

	nsigma = nint(getpar('nsigma'))
	hsigma = getpar('hsigma')

	end

c********************************************************************

	subroutine set_sigma(nsigma,hsigma)

	implicit none

	integer nsigma
	real hsigma

	real getpar

	call putpar('nsigma',float(nsigma))
	call putpar('hsigma',hsigma)

	end 

c********************************************************************
c********************************************************************
c********************************************************************

	subroutine make_sigma_levels(nsigma,hlv)

	implicit none

	integer nsigma
	real hlv(nsigma)

	integer l
	real hl

	if( nsigma .le. 0 ) stop 'error stop make_sigma_levels: nsigma'

        hl = -1. / nsigma
        do l=1,nsigma
          hlv(l) = l * hl
        end do

	end

c********************************************************************

	subroutine make_zeta_levels(lmin,hmin,dzreg,nlv,hlv)

	implicit none

	integer lmin
	real hmin,dzreg
	integer nlv
	real hlv(nlv)

	integer l
	real hbot

	if( dzreg .le. 0. ) stop 'error stop make_zeta_levels: dzreg'

        hbot = hmin
	if( lmin .gt. 0 ) hlv(lmin) = hbot

        do l=lmin+1,nlv
          hbot = hbot + dzreg
          hlv(l) = hbot
        end do

	end

c********************************************************************

	subroutine set_hybrid_depth(lmax,zeta,htot
     +					,hlv,nsigma,hsigma,hlfem)

c sets depth structure and passes it back in hlfem

	implicit none

	integer lmax		!total number of layers
	real zeta		!water level
	real htot		!total depth (without water level)
	real hlv(1)		!depth structure (zeta, sigma or hybrid)
	integer nsigma		!number of sigma levels
	real hsigma		!depth of hybrid closure
	real hlfem(1)		!converted depth values (return)

	logical bsigma
	integer l,i
	real hsig

	bsigma = nsigma .gt. 0

	if( nsigma .gt. 0 ) then
          hsig = min(htot,hsigma) + zeta

	  do l=1,nsigma-1
            hlfem(l) = -zeta - hsig * hlv(l)
	  end do

	  hlfem(nsigma) = -zeta + hsig
	end if

        do l=nsigma+1,lmax
          hlfem(l) = hlv(l)
        end do

	if( nsigma .lt. lmax ) hlfem(lmax) = htot	!zeta or hybrid

c check ... may be deleted

	do l=2,lmax
	  if( hlfem(l) - hlfem(l-1) .le. 0. ) then
	    write(6,*) (hlfem(i),i=1,lmax)
	    stop 'error stop set_hybrid_depth: hlfem'
	  end if
	end do

	end

c********************************************************************






