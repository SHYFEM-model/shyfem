
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
c 14.02.2022	ggu	revisted, some bug fixes
c
c****************************************************************

        subroutine basbathy

c performs bathymetry interpolation in basin
c
c takes care of lat/lon coordinates
c grd has already been copied to basin

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
	integer nhe,nhk,nh
        integer ner,nco,nknh,nelh,nli
	integer isphe

	integer mode
	integer nminimum
	real ufact,umfact
	real :: flag = -999

	integer nt
        real, allocatable :: xt(:)
        real, allocatable :: yt(:)
        real, allocatable :: at(:)
        real, allocatable :: ht(:)

        real, allocatable :: hk(:)
        real, allocatable :: he(:)

c-----------------------------------------------------------------
c general check
c-----------------------------------------------------------------

	allocate(hk(nkn),he(nel))

	if( breadbas ) then
	  write(6,*) 'BAS file has been read...'
	  write(6,*) 'this branch still has to be finished...'
	  stop 'error stop basbathy: interpolating bas file not ready'
	else
	  write(6,*) 'GRD file has been read...'
	  call grd_get_params(nk,ne,nl,nne,nnl)
	  call grd_get_element_depth(he)
	  call grd_get_nodal_depth(hk)
	  if( ball ) then	!interpolate in all items, also with depth
	    he = flag
	    hk = flag
	  end if
	  nhe = count( he /= flag )
	  nhk = count( hk /= flag )
	  write(6,*) 'nkn = ',nk,'  nel = ',ne
	  write(6,*) 'nhk = ',nhk,'  nhe = ',nhe
	  if( nhe == 0 .and. nhk == 0 ) then
	    !ok
	  else if( nhe > 0 .and. nhk > 0 ) then
	    write(6,*) 'existing depth values on nodes and elements'
	    write(6,*) 'before interpolating first delete one of these'
	    stop 'error stop basbathy: nhe>0 and nhk>0'
	  else if( nhe > 0 .and. bnode ) then
	    write(6,*) 'interpolation on nodes requested but nhe>0'
	    write(6,*) 'depth values on elements not allowed'
	    stop 'error stop basbathy: nhe>0 and nodal interpolation'
	  else if( nhk > 0 .and. .not. bnode ) then
	    write(6,*) 'interpolation on elements requested but nhk>0'
	    write(6,*) 'depth values on nodes not allowed'
	    stop 'error stop basbathy: nhk>0 and element interpolation'
	  end if
	end if

	if( bnode .and. bmode == 2 ) then
	  write(6,*) 'uniform interpolation only allowed on elements'
	  stop 'error stop basbathy: incompatible parameters'
	end if

c-----------------------------------------------------------------
c what to do
c-----------------------------------------------------------------

	mode = bmode
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

	nn = max(nk,ne)
	allocate(xt(nn),yt(nn),at(nn),ht(nn))

c-----------------------------------------------------------------
c handling of depth and coordinates
c-----------------------------------------------------------------

	call get_coords_ev(isphe)
	call set_dist(isphe)

	where( iarv == btype ) he = flag

	hev = he
	hkv = hk

c-----------------------------------------------------------------
c node_test
c-----------------------------------------------------------------

	call node_test	!just for coherence

        if( bnode ) then
	  nt = nkn
          call prepare_on_node(nt,xt,yt,at,ht,ufact)
        else
	  nt = nel
          call prepare_on_elem(nt,xt,yt,at,ht,ufact)
        end if

	nh = count( ht(1:nt) /= flag )
	write(6,*) 'items with/without depth: ',nh,nt-nh

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

	if( bnode ) then
	  hkv(1:nkn) = ht(1:nkn)
        else
	  hev(1:nel) = ht(1:nel)
        end if

	call transfer_depth(bnode)	!copy to nodes/elems (also sets hm3v)

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

