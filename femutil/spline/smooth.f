
!--------------------------------------------------------------------------
!
!    Copyright (C) 2002,2004-2005,2011-2012  Georg Umgiesser
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

! smoothing program
!
! contents :
!
! revision log :
! 
! 01.02.2002    ggu     handle periodic line
! 06.08.2002    ggu     read also depth value
! 23.09.2004    ggu     adapted for malta, bug fix
! 19.10.2005    ggu     documentation and description
! 02.12.2011    ggu     bug fix in intpdep() and reduce()
! 02.12.2011    ggu     use depth also for smoothing (change in distxy())
! 27.05.2012    ggu     renamed ndim to nsdim in smooth()
!
! description :
!
! This program smothes and reduces the number of points on a line
!
! The program ask for the GRD file that contains the lines
! It also asks for the smoothing parameter (sigma) and reduction
!	parameter (reduct)
! The units of these values must be in the units of the line coordinates.
! 
! Smoothing:
!
! The program uses a gaussian curve as the smoothing operator, with
! its standard deviation given by sigma. Therefore, if (x,y) coordinates
! of the lines are in meters and sigma=100, the lengthscale for the 
! smoothing is about 100 meters, which means that points closer than 100
! meters to the point will have much more effect than points further away
! than 100 m.
!
! Reduction:
!
! The program uses the parameter reduct to decide if a point should
! be kept or not. If reduct=100, than all points closer than 100 m
! to the reference point are discarded. At the end a line is left with
! a distance between the points of minimum 100 m.
!
! Adaption of resolution:
!
! If the points of the line have no depth values, then the 
! smoothing and reduction is working as described above. If, however,
! depth values do exists, then they are used to change the resolution
! of the line points. For a value of 2, twice the points in the reduction
! algorithm are retained. In other words, the value for reduct is changed
! locally to a value of 100/2=50, and only points closer as 50 are
! eliminated. Line points without depth receive an implicit value of 1.
! The value -1 indicates that this point should never be eliminated 
! nor should it be smoothed.
!
! To do:
!
! the fixed points should be handled better (stop smoothing, reducing)
! adaption of reduce on short lines, look at line globally
!
!***********************************************************************

	program smooth

! smoothing program for lines

	use basin
	use grd
	!use mod_depth

	implicit none

	integer ner
	integer nco,nk,ne,nl,nne,nnl
	integer iabg
	logical bstop
	character*80 file

	logical bperiod
	integer nt,nll

	real, allocatable :: xt(:)
	real, allocatable :: yt(:)
	real, allocatable :: ht(:)
	real, allocatable :: wt(:)
	integer, allocatable :: kt(:)
	real, allocatable :: dist(:)

	integer l,lmax
	integer nline
	integer nnode
	real area,areamax
	integer inodes
	real sigma
	real reduct

	call shyfem_copyright('gridr - smoothing of lines')

!-----------------------------------
! sigma		smoothing parameter (size of gaussian kernel)
! reduct	reduction of points in line (in meters)
!-----------------------------------

	call handle_command_line(file,sigma,reduct,inodes,iabg)

!------------------------------------------------------

	nline = 0
	nnode = 0

        call grd_read(file)
        call grd_get_params(nk,ne,nl,nne,nnl)

        write(6,'(a,5i7)') 'grid parameters: ',nk,ne,nl
        if( nk .le. 0 ) stop

        call grd_to_basin
	!call mod_depth_init(nkn,nel)

	write(6,*) 'nodes    : ',nk
	write(6,*) 'elements : ',ne
	write(6,*) 'lines    : ',nl
	write(6,*) 'nodes in element list : ',nne
	write(6,*) 'nodes in line list    : ',nnl

	allocate(xt(nnl),yt(nnl),ht(nnl),wt(nnl),kt(nnl))
	allocate(dist(nnl))

	call insnod(nkn,ipv)

!--------------------------------------------------------
! find max area (outer line)
!--------------------------------------------------------

	lmax = 0
	areamax = 0.
	do l=1,nl
	  nll = nnl
	  call extrli(l,nl,ipplv,ialv,ipntlv,inodlv,xgv,ygv,hhnv
     +				,xt,yt,ht,kt,nll,nt)
	  call area_of_polygon(nll,xt,yt,area)
	  area = abs(area)
	  if( area > areamax ) then
	    areamax = area
	    lmax = l
	  end if
	end do
	write(6,*) 'outer line: ',lmax,ipplv(lmax),area

!--------------------------------------------------------
! set up standard resolution
!--------------------------------------------------------

	call init_resolution(iabg,nkn,xgv,ygv,nel,nen3v,iarv,hm3v)
	call setup_resolution(inodes,area,sigma,reduct)

	open(99,file='smooth.grd',status='unknown',form='formatted')
	open(98,file='reduce.grd',status='unknown',form='formatted')

	do l=1,nl
	  nll = nnl
	  nline = ipplv(l)
	  call extrli(l,nl,ipplv,ialv,ipntlv,inodlv,xgv,ygv,hhnv
     +				,xt,yt,ht,kt,nll,nt)
	  !call wrline(95,nline,nnode,nll,xt,yt,ht,nt,bperiod)
	  call mkperiod(xt,yt,nll,bperiod)
	  call make_dist(nll,xt,yt,kt,dist,bperiod)
	  !call wrline(96,nline,nnode,nll,xt,yt,ht,nt,bperiod)
	  call intpdep(nline,ht,nll,bperiod)
	  call compute_weights(nll,xt,yt,ht,wt)
	  !call wrline(97,nline,nnode,nll,xt,yt,wt,nt,bperiod)
	  call smooth_line(sigma,xt,yt,wt,nll,bperiod)
	  call wrline(99,nline,nnode,nll,xt,yt,wt,nt,bperiod)
	  call reduce_line(reduct,xt,yt,wt,kt,nll,bperiod)
	  call wrline(98,nline,nnode,nll,xt,yt,wt,nt,bperiod)
	end do

	write(6,*) 'routine finished...'

	close(99)
	close(98)

	write(6,*) 'files smooth.grd and reduce.grd written'

	end

!********************************************************
!********************************************************
!********************************************************

	subroutine smooth_line(sigma,xt,yt,ht,nl,bperiod)

! smoothing

	implicit none

	real sigma
	real xt(nl)
	real yt(nl)
	real ht(nl)
	integer nl
	logical bperiod

	integer i
	real distot
	real, allocatable :: gk(:)
	real, allocatable :: raux(:)
	real, allocatable :: dxy(:)

        if( sigma .le. 0. ) return

! set up dxy

	allocate(gk(-nl:nl))
	allocate(raux(-nl:2*nl))
	allocate(dxy(-nl:2*nl))

	call distxy(nl,xt,yt,ht,dxy,bperiod)

	distot = 0.
	do i=1,nl-1
	  distot = distot + dxy(i)
	end do
	if( bperiod ) distot = distot + dxy(nl)

	write(6,*) 'distance : ',nl,distot,bperiod

! smooth x and y

	call grsmooth(nl,sigma,xt,raux,dxy,ht,bperiod)
	call grsmooth(nl,sigma,yt,raux,dxy,ht,bperiod)

	end

!********************************************************

	  subroutine reduce_line(reduct,xt,yt,ht,kt,nl,bperiod)

! reduces points in line

	implicit none

	real reduct
	real xt(nl)
	real yt(nl)
	real ht(nl)
	integer kt(nl)
	integer nl
	logical bperiod

	integer i,nnew
	real rr,rtot
	real dist,dx,dy
	real h

! very simplicistic approach

	if( reduct <= 0 ) return

	nnew = 0
	rtot = 0.
	if( ht(1) .lt. 0. .or. ht(1) .ge. 0. ) then	!bug fix 23.09.2004
	    nnew = nnew + 1
	    xt(nnew) = xt(1)
	    yt(nnew) = yt(1)
	    ht(nnew) = ht(1)
	end if
	do i=2,nl
	  dx = xt(i) - xt(i-1)
	  dy = yt(i) - yt(i-1)
	  h = ht(i)
	  dist = sqrt( dx*dx + dy*dy )
	  !write(6,*) dist,h,dist/h,rtot,reduct
	  if( h .gt. 0. ) dist = dist * h
	  rtot = rtot + dist
	  if( rtot .gt. reduct .or. h .lt. 0. ) then
	    nnew = nnew + 1
	    xt(nnew) = xt(i)
	    yt(nnew) = yt(i)
	    ht(nnew) = ht(i)
	    rtot = 0.
	  end if
	end do

	if( nnew .lt. 3 ) then
	  write(6,*) 'line too short -> eliminated'
	  nnew = 0
	end if

	write(6,*) 'reduce_line: new nodes = ',nnew
	nl = nnew

	end

!********************************************************

	  subroutine intpdep(nline,ht,nl,bperiod)

! interpolates depth values for points in line
!
! negative values are not interpolated and are left alone
! after calling this subroutine all ht are set (maybe to -1)

	implicit none

        integer nline
	real ht(1)
	integer nl
	logical bperiod

	integer i
	integer ifirst,ilast,inext,nonzero
	integer nval
	real value,vfirst,vlast
	real aux,hflag
	logical bdebug

	bdebug = .true.
	bdebug = nline .eq. 9
	bdebug = nline .eq. 6
	bdebug = .false.

        hflag = -990.

!------------------------------------------------------------
! look for values greater than 0
!------------------------------------------------------------

	nonzero = 0
	ifirst = 0
	ilast = 0

	do i=1,nl
	  if( ht(i) .gt. 0. ) then
	    nonzero = nonzero + 1
	    if( nonzero .eq. 1 ) ifirst = i
	    ilast = i
	  end if
	end do

      if( bdebug ) write(6,*) 'nz: ',nl,nonzero,ifirst,ilast

!------------------------------------------------------------
! if no value or 1 value found -> set everything constant
!------------------------------------------------------------

	if( nonzero .le. 1 ) then
	  if( nonzero .eq. 0 ) then
	    value = 1.
	  else
	    value = ht(ifirst)
	  end if
	  do i=1,nl
	    if( ht(i) /= -1. ) ht(i) = value
	  end do
	  return
	end if

!------------------------------------------------------------
! interpolate extreme values
!------------------------------------------------------------

	vfirst = ht(ifirst)
	vlast = ht(ilast)
	nval = nl - ilast + ifirst
	aux = (vfirst-vlast)/nval
	inext = 0

	if( bdebug ) write(6,*) bperiod
	if( bdebug ) write(6,*) ifirst,ilast,nval,vfirst,vlast

	i = 1
	do while( ht(i) == -1. )
	  i = i + 1
	end do
	if( bperiod ) then
	  inext = nl - ilast + i
	  value = vlast + inext * aux
	else
	  value = ht(ifirst)
	end if
	if( bdebug ) write(6,*) 'first: ',i,inext,value,ht(i)
	ht(i) = value
	ifirst = i

	i = nl
	do while( ht(i) == -1. )
	  i = i - 1
	end do
	if( bperiod ) then
	  inext = i - ilast
	  value = vlast + inext * aux
	else
	  value = ht(ilast)
	end if
	if( bdebug ) write(6,*) 'last: ',i,inext,value,ht(i)
	ht(i) = value
	ilast = i

!------------------------------------------------------------
! interpolate central values
!------------------------------------------------------------

	do while( ifirst .lt. ilast )
	  do i=ifirst+1,nl
	    if( ht(i) .gt. 0. ) goto 1
	  end do
	  stop 'error stop intpdep: impossible branch (1)'
    1	  continue
	  inext = i
	  nval = inext - ifirst
	  if( bdebug ) write(6,*) ifirst,inext,ilast,ht(ifirst),ht(inext)
	  do i=ifirst+1,inext-1
	    value = ht(ifirst) + (i-ifirst) * (ht(inext)-ht(ifirst))/nval
	    if( ht(i) .ge. 0. .or. ht(i) .lt. hflag ) ht(i) = value
	    if( bdebug ) write(6,*) i,value,ht(i)
	  end do
	  ifirst = inext
	end do

!------------------------------------------------------------
! do some sanity check
!------------------------------------------------------------

	do i=1,nl
	  if( ht(i) .eq. 0 ) then
	    write(6,*) i,nl,ht(i)
	    stop 'error stop intpdep: values not interpolated'
	  end if
	  if( bdebug ) write(6,*) i,ht(i)
	end do

!------------------------------------------------------------
! end of routine
!------------------------------------------------------------

	end

!********************************************************

	subroutine handle_command_line(file,sigma,reduct,inodes,iabg)

	use clo

	implicit none

	character*(*) file
	real sigma
	real reduct
	integer inodes
	integer iabg

	integer nfile

        call clo_init('smooth','grd-file','1.0')

        call clo_add_info('smoothes line and reduces points')
        call clo_add_option('sigma',0.
     +            ,'standard deviation for smoothing')
        call clo_add_option('reduct',0.
     +            ,'reduction of points with smaller distance')
        call clo_add_option('inodes',0
     +            ,'number of desired nodes (inside outer line)')
        call clo_add_option('bg',0
     +            ,'element type of background grid')

        call clo_parse_options(1)       !expecting 1 file

        call clo_get_option('sigma',sigma)
        call clo_get_option('reduct',reduct)
        call clo_get_option('inodes',inodes)
        call clo_get_option('bg',iabg)

        nfile = clo_number_of_files()
        if( nfile > 0 ) call clo_get_file(1,file)

	write(6,*) 'file name: ',trim(file)
	write(6,*) 'sigma: ',sigma
	write(6,*) 'reduct: ',reduct
	write(6,*) 'inodes: ',inodes
	write(6,*) 'iabg: ',iabg

	end

!********************************************************

