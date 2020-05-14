
!--------------------------------------------------------------------------
!
!    Copyright (C) 1999,2003-2005,2009-2017,2019  Georg Umgiesser
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

c revision log :
c
c 06.04.1999	ggu	completely restructured
c 04.06.1999	ggu	new statistics are computed
c 08.09.2003	ggu	mode 5 -> write depth values from elements
c 23.09.2004	ggu	interpolq() changed for bathy interpolation
c 02.10.2004	ggu	interpole() for exponential interpolation
c 12.05.2005	ggu	pass hmin to interpolation functions
c 06.04.2009	ggu	read param.h
c 24.04.2009	ggu	new call to rdgrd()
c 21.05.2009	ggu	restructured to allow for nodal interpolation
c 16.12.2010	ggu	bug fix in transfer_depth()
c 02.12.2011	ggu	introduction of nminimum - hardcoded for now
c 16.03.2012	ggu	autoregression introduced (make_auto_corr,interpola)
c 16.03.2012	ggu	default value for umfact set to 3, new mode = 3
c 01.06.2012	ggu	some more changes
c 13.06.2013	ggu	copy_depth() renamed to transfer_depth()
c 13.02.2014	ggu	new data written, can read also bas file
c 05.03.2014	ggu	subroutines copied to other routine
c 10.10.2015	ggu	changed VERS_7_3_2
c 16.12.2015	ggu	depth ht is now passed in for square interpol
c 11.04.2016	ggu	meaning of ufact has changed
c 25.05.2017	ggu	changed VERS_7_5_28
c 16.02.2019	ggu	changed VERS_7_5_60
c
c****************************************************************

        subroutine basbathy

c performs bathymetry interpolation in basin
c
c takes care of lat/lon coordinates

	use mod_depth
	use evgeom
	use basin
	use grd
	use basutil

	implicit none

	integer np
	real, allocatable :: xp(:)
	real, allocatable :: yp(:)
	real, allocatable :: dp(:)
	real, allocatable :: ap(:)

	integer node,nit
	integer n,i,nn
	integer nk,ne,nl,nne,nnl
        integer ner,nco,nknh,nelh,nli
	integer isphe

	integer mode,ike,idepth
	integer nminimum
	real ufact,umfact
	real :: flag = -999

	integer nt
        real, allocatable :: xt(:)
        real, allocatable :: yt(:)
        real, allocatable :: at(:)
        real, allocatable :: ht(:)

c-----------------------------------------------------------------
c what to do
c-----------------------------------------------------------------

	mode = bmode
	idepth = 1
	if( ball ) idepth = 2
	if( btype >= 0 ) idepth = 3
	ike = 1
	if( bnode ) ike = 2
	ufact = usfact		!factor for standard deviation (abs if negative)
	umfact = uxfact		!factor for maximum radius

	nminimum = 1	!minimum number of points to be used for interpolation

c-----------------------------------------------------------------
c read in bathymetry file
c-----------------------------------------------------------------

	write(6,*) 'reading bathymetry file : ',trim(bfile)

	call grd_read(bfile)
	call grd_get_params(nk,ne,nl,nne,nnl)

	np = nk
	allocate(xp(np),yp(np),dp(np),ap(np))

	call grd_get_nodes(np,xp,yp,dp)
	call grd_close

	write(6,*) 'finished reading bathymetry file : ',trim(bfile)

c-----------------------------------------------------------------
c allocate arrays for basin
c-----------------------------------------------------------------

	nn = max(nkn,nel)
	allocate(xt(nn),yt(nn),at(nn),ht(nn))

c-----------------------------------------------------------------
c handling of depth and coordinates
c-----------------------------------------------------------------

	call get_coords_ev(isphe)
	call set_dist(isphe)

	hkv = flag
	call set_depth_i(idepth,btype,nknh,nelh)

c-----------------------------------------------------------------
c node_test
c-----------------------------------------------------------------

	call node_test

        if( ike .eq. 1 ) then                           !elementwise
          call prepare_on_elem(nt,xt,yt,at,ht,ufact)
        else                                            !nodewise
          call prepare_on_node(nt,xt,yt,at,ht,ufact)
        end if

c-----------------------------------------------------------------
c interpolate
c-----------------------------------------------------------------

	write(6,*) 'starting interpolation...'

	if( mode .eq. 1 ) then
	  call interpole(np,xp,yp,dp,nt,xt,yt,at,ht,umfact,nminimum)
        else if( mode .eq. 2 ) then
	  call interpolq(np,xp,yp,dp,ht)
	else if( mode .eq. 3 ) then
	  call make_auto_corr(np,xp,yp,dp,ap,ufact)
	  call interpola(np,xp,yp,dp,ap,nt,xt,yt,at,ht)
        else
          write(6,*) 'wrong choice for mode : ',mode
          stop 'error stop'
	end if

	write(6,*) 'finished interpolation...'

c-----------------------------------------------------------------
c transfer depth
c-----------------------------------------------------------------

        if( ike .eq. 1 ) then                           !elementwise
	  hev(1:nel) = ht(1:nel)
        else                                            !nodewise
	  hkv(1:nkn) = ht(1:nkn)
        end if

	call transfer_depth(ike)	!copy to nodes/elems (also sets hm3v)

c-----------------------------------------------------------------
c write
c-----------------------------------------------------------------

        call basin_to_grd

        call grd_write('basbathy.grd')
        write(6,*) 'The basin has been written to basbathy.grd'

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	return
	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

        subroutine set_depth_i(idepth,btype,nknh,nelh)

c handles depth values

        use mod_depth
        use basin

        implicit none

        integer idepth          !1: only at missing points 2: everywhere
				!3: only at type=btype
	integer btype		!type where to interpolate (if > 0)
        integer nknh,nelh       !return - depth values found

        integer k,ie,ii,ia
        real flag

        flag = -999.

        if( idepth .eq. 2 ) then
	  hkv = flag
	  hev = flag
        end if

        if( idepth .eq. 3 ) then
          do ie=1,nel
            ia = iarv(ie)
            if( ia == btype ) then
		hev(ie) = flag
                do ii=1,3
                  k = nen3v(ii,ie)
		  hkv(k) = flag
                end do
            end if
          end do
	end if

        nknh = 0
        nelh = 0

        do k=1,nkn
          if( hkv(k) .le. flag ) nknh = nknh + 1
        end do
        do ie=1,nel
          if( hev(ie) .le. flag ) nelh = nelh + 1
        end do

        end

c*******************************************************************

