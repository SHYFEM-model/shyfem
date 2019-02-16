
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

c revision log :
c
c 06.04.1999    ggu     completely restructured
c 04.06.1999    ggu     new statistics are computed
c 08.09.2003    ggu     mode 5 -> write depth values from elements
c 23.09.2004    ggu     interpolq() changed for bathy interpolation
c 02.10.2004    ggu     interpole() for exponential interpolation
c 01.11.2004    ggu     whole program simplyfied
c 06.12.2008    ggu     smoothing introduced
c 06.04.2009    ggu     read param.h
c 29.05.2009    ggu     does only depth limiting and smoothing
c 20.11.2009    ggu     possibility to smooth only on specific areas
c 30.03.2011    ggu     new routines to delete elements
c 13.06.2013    ggu     copy_depth() renamed to transfer_depth()
c
c****************************************************************

        subroutine bas_smooth

c performs modifications on basin
c
c takes care of lat/lon coordinates

	use mod_geom
	use mod_depth
	use evgeom
	use basin
	use grd
	use basutil

	implicit none

	integer niter
	integer ike
	real alpha
	real f(4)

	logical is_depth_unique
	integer iscanf

c-----------------------------------------------------------------
c what to do
c-----------------------------------------------------------------

        !write(6,*)
        !write(6,*) 'Enter parameters for smoothing:'
        !write(6,*)
        !write(6,*) 'Either alpha or a1,h1,a2,h2'
        !write(6,*)
	!write(6,*) 'Enter parameters: '
	!read(5,'(a)') line
	!n = iscanf(line,f,4)
	!if( n .ne. 1 .and. n .ne. 4 ) goto 96
	!if( n .eq. 1 ) then
	!  f(3) = f(1)
	!  f(2) = 0.
	!  f(4) = 10000.
	!end if
	!write(6,*) 'parameters used :',(f(i),i=1,4)
        !write(6,*)

	niter = iter
	alpha = asmooth

	f(1) = alpha
	f(2) = 0.
	f(3) = alpha
	f(4) = 10000.

c-----------------------------------------------------------------
c handling depth
c-----------------------------------------------------------------

        ike = 1
	if( is_depth_unique() ) ike = 2

c-----------------------------------------------------------------
c smooth
c-----------------------------------------------------------------

	call limit_depth(ike,hmin,hmax)

	if( niter .gt. 0 ) then
	  call smooth_bathy(ike,niter,f)
	end if

	call transfer_depth(ike)	!copy to nodes/elements

c-----------------------------------------------------------------
c special
c-----------------------------------------------------------------

	!call delete_elements(0.)

c-----------------------------------------------------------------
c write
c-----------------------------------------------------------------

        call basin_to_grd
        call grd_write('bassmooth.grd')
        write(6,*) 'The basin has been written to bassmooth.grd'

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	return
	end

c*******************************************************************

	subroutine limit_depth(ike,hmin,hmax)

	use mod_depth
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer ike
	real hmin,hmax

	integer ie,k

	if( ike .eq. 1 ) then		!elementwise

	  do ie=1,nel
	    if( hev(ie) .gt. hmax ) hev(ie) = hmax
	    if( hev(ie) .lt. hmin ) hev(ie) = hmin
	  end do

	else				!nodewise

	  do k=1,nkn
	    if( hkv(k) .gt. hmax ) hkv(k) = hmax
	    if( hkv(k) .lt. hmin ) hkv(k) = hmin
	  end do

	end if

	end

c*******************************************************************

	subroutine smooth_bathy(ike,niter,f)

c smoothes depth values

	implicit none

	integer ike
	integer niter
	real f(4)

	if( ike .eq. 1 ) then
	  call smooth_bathy_elem(niter,f)
	else
	  stop 'error stop: cannot yet smooth on nodes'
	end if

	end

c*******************************************************************

	subroutine smooth_bathy_elem(niter,f)

c smoothes depth values

	use mod_depth
	use evgeom
	use basin

	implicit none

	integer niter
	real f(4)

	integer ie,ii,k,i,nok
	real x,y,d,h,hold,hnew,ao
	real hmin
	real alpha,beta
        real h1,h2,a1,a2
	integer iaux,inum,ityp
	real xt(3), yt(3)
	real v1v(nkn)
	integer ihev(nel)
	logical inconvex

c--------------------------------------------------------------
c set parameters
c--------------------------------------------------------------

	a1 = f(1)
	h1 = f(2)
	a2 = f(3)
	h2 = f(4)

	write(6,*) 'smoothing bathymetry: ',niter
	write(6,*) '  params: ',(f(i),i=1,4)

c--------------------------------------------------------------
c iterate over smoother
c--------------------------------------------------------------

	do i=1,niter

c	  -----------------------------------------------
c	  compute values at nodes (averages of element)
c	  -----------------------------------------------

	  do k=1,nkn
	    hkv(k) = 0.
	    v1v(k) = 0.
	  end do

	  do ie=1,nel
	    ao = ev(10,ie)
	    do ii=1,3
	      k = nen3v(ii,ie)
	      hkv(k) = hkv(k) + ao * hev(ie)
	      v1v(k) = v1v(k) + ao
	    end do
	  end do

	  do k=1,nkn
	    hkv(k) = hkv(k) / v1v(k)
	  end do

c	  -----------------------------------------------
c	  average to element
c	  -----------------------------------------------

	  nok = 0

	  do ie=1,nel
	    h = 0.
	    do ii=1,3
	      k = nen3v(ii,ie)
	      h = h + hkv(k)
	    end do
	    hnew = h / 3.
	    hold = hev(ie)

            beta = (hold-h1)/(h2-h1)
            beta = max(beta,0.)
            beta = min(beta,1.)

            alpha = a1 + beta * (a2-a1)
            alpha = a1

	    !call coords_ok(ie,alpha)	!customize to smooth on specific areas
	    if( alpha .gt. 0. ) nok = nok + 1
            !write(6,*) ie,hold,alpha

	    hev(ie) = (1.-alpha) * hold + alpha * hnew
	  end do

c	  -----------------------------------------------
c	  write to terminal
c	  -----------------------------------------------

	  write(6,*) 'pass ',i,' of ',niter,'  elements smoothed ',nok

	end do

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	end

c*******************************************************************

	subroutine coords_ok(ie,alpha)

	implicit none

	integer ie
	real alpha

	real x,y,x1,y1,x2,y2

	x1 = 103242.
	y1 = 62770.
	x2 = 104133.
	y2 = 64304.

	x1 = 0.
	y1 = 0.
	x2 = 10000.
	y2 = 10000.

	call baric(ie,x,y)

	if( x .ge. x1 .and. x .le. x2 ) then
	  if( y .ge. y1 .and. y .le. y2 ) then
	    return				!ok - keep alpha
	  end if
	end if

	alpha = 0.0

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine delete_elements(hmin)

c deletes elements with depth lower then hmin

	use mod_geom
	use mod_depth
	use basin

	implicit none

	integer nlkdi
	real hmin

	integer icon(nel)

	integer icol,ibig,ie

	write(6,*)
	write(6,*) 'deleting elements with depth <= ',hmin
	write(6,*)

	call delete_elements_depth(hmin)

	call set_ev
	call set_geom

	write(6,*)
	write(6,*) 'checking connectivity of basin'
	write(6,*)

	call check_connection(icon,icol,ibig)

	do ie=1,nel
	  if( icon(ie) .ne. ibig ) then
	    hev(ie) = hmin - 1.
	  end if
	end do

	write(6,*)
	write(6,*) 'deleting not connected areas'
	write(6,*)

	call delete_elements_depth(hmin)

	call set_ev
	call set_geom

	write(6,*)
	write(6,*) 'final check'
	write(6,*)

	call check_connection(icon,icol,ibig)

	write(6,*)
	write(6,*) 'end deleting elements'
	write(6,*)

	end

c*******************************************************************

	subroutine delete_elements_depth(hmin)

c deletes elements with depth lower then hmin

	use mod_depth
	use basin

	implicit none

	real hmin

	integer ie,ii,k
	integer n,ieh,kh

	integer ind(nel)	!index for nodes/elements to substitute
	integer rind(nel)	!reverse index

c-----------------------------------------
c delete elements
c-----------------------------------------

	n = nel
	call determine_shift(n,ind,hev,hmin)

	do ie=1,nel
	  ieh = ind(ie)
	  if( ieh .gt. 0 ) then
	    ipev(ie) = ipev(ieh)
	    iarv(ie) = iarv(ieh)
	    hev(ie) = hev(ieh)
	    do ii=1,3
	      nen3v(ii,ie) = nen3v(ii,ieh)
	    end do
	  end if
	end do

	write(6,*) 'new elements: ',nel,n
	nel = n

c-----------------------------------------
c flag unused nodes and delete (using hkv)
c-----------------------------------------

	do k=1,nkn
	  hkv(k) = hmin - 1
	  rind(k) = 0
	end do

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    hkv(k) = 1
	  end do
	end do

	n = nkn
	call determine_shift(n,ind,hkv,hmin)

	do k=1,nkn
	  kh = ind(k)
	  if( kh .gt. 0 ) then
	    ipv(k) = ipv(kh)
	    xgv(k) = xgv(kh)
	    ygv(k) = ygv(kh)
	  end if
	  rind(kh) = k
	end do

	write(6,*) 'new nodes: ',nkn,n
	nkn = n

c-----------------------------------------
c adjust element index
c-----------------------------------------

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    kh = rind(k)
	    if( kh .gt. 0 ) nen3v(ii,ie) = kh
	  end do
	end do

c-----------------------------------------
c end of routine
c-----------------------------------------

	end

c*******************************************************************

	subroutine determine_shift(n,ind,h,hmin)

c determines what element/nodes must be shifted

	implicit none

	integer n
	integer ind(n)
	real h(n)
	real hmin

	integer i,ih
	logical bfound

	do i=1,n
	  ind(i) = 0
	end do

	bfound = .true.
	i = 0
	ih = n+1
	do while( i+1 .lt. ih )
	  i = i + 1
	  if( h(i) .le. hmin ) then
	    bfound = .false.
	    do while( ih-1 .gt. i .and. .not. bfound )
	      ih = ih - 1
	      if( h(ih) .gt. hmin ) then
	        bfound = .true.
		ind(i) = ih
	      end if
	    end do
	    !write(6,*) i,ih
	  end if
	end do

	if( .not. bfound ) i = i - 1	!last element to be excluded
	n = i

	end

c*******************************************************************

	subroutine check_connection(icon,icol,ibig)

	use mod_geom
	use basin

	implicit none

	integer icon(nel)
	integer icol,ibig

	integer ie
	integer i,nc,ic
	integer icolor(nel)

	icol = 0
	icon = 0
	icolor = 0

	do ie=1,nel
	  if( icon(ie) .eq. 0 ) then
	    icol = icol + 1
	    call color_area(ie,icol,icon)
	  end if
	end do

	do ie=1,nel
	  i = icon(ie)
	  icolor(i) = icolor(i) + 1
	end do
	
	nc = 0
	do i=1,icol
	  ic = icolor(i)
	  if( ic .gt. nc ) then
	    ibig = i
	    nc = ic
	  end if
	end do

	write(6,*) 'number of connected areas: ',icol
	write(6,*) 'biggest area: ',ibig,nc

	end

c*******************************************************************

	subroutine color_area(iestart,icol,icon)

	use basin
	use mod_geom

	implicit none

	integer iestart,icol
	integer icon(nel)

	integer ip,ien,ii,ie
	integer list(nel)

	ip = 1
	list(ip) = iestart

	do while( ip .gt. 0 ) 
	  ie = list(ip)
	  icon(ie) = icol
	  !write(6,*) icol,ip,ie
	  ip = ip -1
	  do ii=1,3
	    ien = ieltv(ii,ie)
	    if( ien .gt. 0 ) then
	      if( icon(ien) .eq. 0 ) then
	        ip = ip + 1
	        list(ip) = ien
	      else if( icon(ien) .ne. icol ) then
	        goto 99
	      end if
	    end if
	  end do
	end do

	return
   99	continue
	write(6,*) ip,ie,ien,icol,icon(ien)
	stop 'error stop color_area: internal error (1)'
	end

c*******************************************************************

