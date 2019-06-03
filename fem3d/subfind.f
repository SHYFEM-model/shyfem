
!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2018  Georg Umgiesser
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

c routines for interpolation onto regular grid
c
c contents :
c
c subroutine setgeo(x0,y0,dx,dy,flag)
c		sets grid values for interpolation
c subroutine getgeo(x0,y0,dx,dy,flag)
c		gets grid values for interpolation
c subroutine getgeoflag(flag)
c		gets flag value for interpolation
c
c subroutine av2am(av,am,ip,jp)
c		interpolation of av onto a regular net
c subroutine av2amk(bwater,av,am,ip,jp)
c		interpolation of av onto a regular net
c function intri(x,y,xp,yp)
c		point in triangle or not
c function intrid(x,y,xp,yp)
c		point in triangle or not (double precision version)
c subroutine am2av(am,av,ip,jp)
c		interpolation of am onto finite element mesh
c function am2val(am,ip,jp,xx,yy)
c		interpolation of am onto finite element mesh
c subroutine ave2am(av,am,ip,jp)
c		interpolation of av (elementwise) onto a regular net
c
c subroutine mkmask(bwater,zv,href,hzoff)
c		makes mask in element
c
c subroutine mimareg(am,ip,jp,amin,amax)
c		computes min/max of regular matrix (without flag values)
c subroutine a2char(am,ac,ip,jp)
c		creates 1 char representation of matrix
c subroutine prchar(ac,ip,jp)
c		prints 1 char representation of matrix
c
c subroutine femintp(ie,z,xp,yp,zp)
c               interpolation in element (with ev)
c subroutine elemintp(x,y,z,xp,yp,zp)
c               interpolation in element (no ev)
c
c subroutine find_elem_from_old(ieold,xp,yp,ielem)
c		finds element for point (xp,yp) starting from ieold
c subroutine find_element(xp,yp,ielem)
c		finds element for point (xp,yp)
c function in_element(ie,xp,yp)
c		checks if point (xp,yp) is in element ie
c subroutine get_xy_elem(ie,x,y)
c		returns x,y of vertices of element ie
c
c revision log :
c
c 18.11.1998	ggu	routine commented
c 18.11.1998	ggu	routine setgeo introduced
c 19.11.1998	ggu	routines a2char, prchar added
c 19.10.1999	ggu	routine mkmask added from subutl
c 25.11.2004	ggu	new routines femintp and elemintp for interpolation
c 14.03.2005	ggu	new routines for interpolation in element
c 11.03.2009	ggu	new helper routine getgeoflag()
c 12.06.2009	ggu	passing to double precision, intrid, bug bug_f_64bit
c 26.01.2011	ggu&mbj	handling extrapolation in am2av()
c 27.01.2011	ggu&ccf	bug fix in find_elem_from_old() BUG_27.01.2011
c 31.03.2011	ggu	new routine elemmask()
c 24.11.2011	ggu	new routine find_close_elem()
c 20.06.2012	ggu	new routine get_scal_elem()
c 07.10.2012	ggu	new routine av2fm()
c 10.10.2012	ggu	new routine fm2am2d() and fm2am3d()
c 26.10.2012	ggu	bug fix: do not access not existing storage
c 30.05.2014	ggu	in av2amk() do not interpolate for flag values
c 07.07.2014	ggu	new routine intp_reg()
c 25.09.2015	ggu	new routines intp_reg_nodes(), intp_reg_elems()
c 05.05.2016	ggu	file restructured (module)
c 14.05.2016	ggu	allow for extension of grid -> bregextend
c 23.06.2016	ggu	allow for eps in computing box
c 23.09.2016	ggu	allow for eps in computing box and reg intp
c 23.04.2017	ggu	new routine intp_reg_single_nodes()
c 25.05.2017	ggu	changed VERS_7_5_28
c 16.02.2019	ggu	changed VERS_7_5_60
c 13.03.2019	ggu	changed VERS_7_5_61
c
c notes :
c
c pxareg,pyareg         coordinates of lower left point of matrix
c pxdreg,pydreg         grid size of matrix
c pzlreg                value of z for land points
c
c******************************************************
c******************************************************
c******************************************************

	subroutine find_close_elem(ieold,xp,yp,ielem)

c finds element for point (xp,yp) starting from ieold
c
c uses data structure ev and ieltv

	use mod_geom
	use evgeom
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer ieold
	real xp,yp
	integer ielem	!element number on return

	logical binit,bdebug
	integer ie,ii,iside,lmax,loop
	real xi,ximin
	double precision a(3),b(3),c(3)

	logical in_element

	bdebug = .false.
	lmax = 10

	call is_init_ev(binit)

c-------------------------------------------------------------
c check if old element is given -> if not test all elements
c-------------------------------------------------------------

	if( bdebug ) write(6,*) 'ggu_xi (1) ',ieold
	if( ieold .le. 0 .or. ieold .gt. nel ) then
	  call find_element(xp,yp,ielem)
	  return
	end if

	if( bdebug ) write(6,*) 'ggu_xi (2) ',ieold
	if( .not. binit ) then
	  call find_elem_from_old(ieold,xp,yp,ielem)
	  return
	end if

c-------------------------------------------------------------
c start from old element
c-------------------------------------------------------------

	if( bdebug ) write(6,*) 'ggu_xi (3) ',ieold

	loop = 0
	ie = ieold
	do while( ie .gt. 0 )
	  iside = 0
	  ximin = 1.e+30
	  call xi_abc(ie,a,b,c)
	  do ii=1,3
	    xi = a(ii) + b(ii)*xp + c(ii)*yp
	    if( bdebug ) write(6,*) 'ggu_xiii ',ie,ii,xi
	    if( xi .lt. ximin ) then
	      ximin = xi
	      iside = ii
	    end if
	  end do
	  if( ximin .ge. 0. ) then
	    ielem = ie
	    return
	  end if
	  if( iside .le. 0 .or. iside .gt. 3 ) then
	    if( bdebug ) write(6,*) '******** ',iside,ie,ximin,xi
	    ie = 0
	  else
	    ie = ieltv(iside,ie)
	    if( bdebug ) write(6,*) 'ggu_xiii iterate',ie,iside,ximin
	  end if
	  loop = loop + 1
	  if( loop .gt. lmax ) ie = 0
	end do

	ielem = 0

	end

c******************************************************

	subroutine find_elem_from_old(ieold,xp,yp,ielem)

c finds element for point (xp,yp) starting from ieold

	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer ieold
	real xp,yp
	integer ielem	!element number on return

	logical in_element
	integer iem,iep

c-------------------------------------------------------------
c check if old element is given -> if not test all elements
c-------------------------------------------------------------

	if( ieold .le. 0 .or. ieold .gt. nel ) then
	  call find_element(xp,yp,ielem)
	  return
	end if

c-------------------------------------------------------------
c check if in old element
c-------------------------------------------------------------

	if( in_element(ieold,xp,yp) ) then
	  ielem = ieold
	  return
	end if

c-------------------------------------------------------------
c start from old element going upwards and downwards
c-------------------------------------------------------------

	iem = ieold-1
	if( iem .lt. 1 ) iem = nel		!BUG_27.01.2011
	iep = ieold+1
	if( iep .gt. nel ) iep = 1		!BUG_27.01.2011

	do while( iem .ne. ieold .and. iep .ne. ieold )
	  if( in_element(iem,xp,yp) ) then
	    ielem = iem
	    return
	  end if
	  iem = iem - 1
	  if( iem .lt. 1 ) iem = nel

	  if( in_element(iep,xp,yp) ) then
	    ielem = iep
	    return
	  end if
	  iep = iep + 1
	  if( iep .gt. nel ) iep = 1
	end do

c-------------------------------------------------------------
c no element found
c-------------------------------------------------------------

	ielem = 0

c-------------------------------------------------------------
c end of routine
c-------------------------------------------------------------

	end

c******************************************************

	subroutine find_element(xp,yp,ielem)

c finds element for point (xp,yp)

	use basin, only : nkn,nel,ngr,mbw

	implicit none

	real xp,yp
	integer ielem	!element number on return


	integer ie
	logical in_element

	do ie=1,nel
	  if( in_element(ie,xp,yp) ) then
		  ielem = ie
		  return
	  end if
	end do

	ielem = 0

	end

c******************************************************

	function in_element(ie,xp,yp)

c checks if point (xp,yp) is in element ie

	use basin

	implicit none

	logical in_element
	integer ie
	real xp,yp

	integer ii,k,in
	real xmin,ymin,xmax,ymax
	real x(3),y(3)

	integer intri

	in_element = .false.

	do ii=1,3
	    k = nen3v(ii,ie)
	    x(ii) = xgv(k)
	    y(ii) = ygv(k)
	end do

	xmin = min(x(1),x(2),x(3))
	ymin = min(y(1),y(2),y(3))
	xmax = max(x(1),x(2),x(3))
	ymax = max(y(1),y(2),y(3))

	if( xp .ge. xmin .and. xp .le. xmax ) then
	  if( yp .ge. ymin .and. yp .le. ymax ) then
		in = intri(x,y,xp,yp)
		if( in .gt. 0 ) in_element = .true.
	  end if
	end if

	end

c******************************************************

