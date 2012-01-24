c
c $Id: lagrange_init.f,v 1.8 2010-03-11 15:36:38 georg Exp $
c
c initialization routines
c
c revision log :
c
c 12.12.2007    ggu	written from scratch
c 11.09.2009    ggu	more checks for initial seeding (linepoint_in_element)
c 11.03.2010    ggu	small bug fix in lagbound_seed_particles() -> init np
c 11.03.2010    ggu	in lgr_init_line() check if at least one line is read
c 06.10.2011    ggu	in check_elements() initialize iout, icheck (bug)
c
c*******************************************************************

	subroutine lgr_init_shell

c manages release of particles

	implicit none

	character*80 line
	integer npoints
	real getpar

	npoints = nint(getpar('nbdy'))
        call getfnm('lagra',line)

        if( line .eq. ' ' ) then      !total lagoon
            call lgr_init_total(npoints)
        else                          !use boundary
            call lgr_init_line(npoints)
        end if

	end

c*******************************************************************

	subroutine lgr_init_total(npoints)

c manages release of particles over whole domain

	implicit none

	integer npoints

	integer n
	real tot_area,dxy
	real x(1),y(1)

	call compute_total_area(tot_area)
	dxy = sqrt( tot_area / npoints )

	n = -1
	!call lagbound_seed_particles(dxy,n,x,y)	!should be the same
	call lgr_init_all(dxy)

	end

c*******************************************************************

	subroutine lgr_init_line(npoints)

c manages release of particles inside polygon(s)

	implicit none

	integer npoints

	integer ndim
	parameter (ndim=100)

	real x(ndim),y(ndim)

	integer iunit,n,nline
	real tot_area,dxy
	logical lagbound_read_next

	call lagbound_get_area(tot_area)
        dxy = sqrt( tot_area / npoints )

	write(6,*) 'lgr_init_line: ',npoints,tot_area,dxy

	call lagbound_open_file(iunit)
	nline = 0
	do while( lagbound_read_next(iunit,ndim,n,x,y) )
	  nline = nline + 1
	  write(6,*) 'new line read: ',nline,n
	  call lagbound_seed_particles(dxy,n,x,y)
	end do
	call lagbound_close_file(iunit)

	if( nline .eq. 0 ) then
	  write(6,*) 'no line found in file'
	  write(6,*) 'cannot initialize particles'
	  stop 'error stop lgr_init_line: no line'
	end if

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine lagbound_get_area(arealine)

c computes area of polygon(s)

	implicit none

	real arealine

	integer ndim
	parameter (ndim=100)

	real x(ndim),y(ndim)

	logical lagbound_read_next
	logical debug
	integer iunit,n,i
	real area,areatot
	real areapoly

	debug = .false.
	areatot = 0.

	call lagbound_open_file(iunit)

	do while( lagbound_read_next(iunit,ndim,n,x,y) )
	  area = areapoly(n,x,y)
	  areatot = areatot + area
	  if( debug ) then
	    write(6,*) 'new line ',n,area
	    do i=1,n
	      write(6,*) i,x(i),y(i)
	    end do
	  end if
	end do

	call lagbound_close_file(iunit)

	write(6,*) 'total area: ',areatot

	arealine = areatot

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine laggrd_read_file

	implicit none

	include 'param.h'

c	--------------------------------------
	integer nlidim,nlndim
	parameter(nlidim=100,nlndim=nkndim)
	real xbnd(nlndim), ybnd(nlndim)
	common /xybnd/xbnd,ybnd
	integer nline,iline
	integer ibnd(0:nlidim)
	common /ibnd/nline,iline,ibnd
	save /xybnd/,/ibnd/
c	--------------------------------------

	character*80 file
	logical bstop
	integer nco,nkn,nel,nli
	integer ipnv(nkndim)
	integer ipev(neldim)
	integer iplv(nlidim)
	integer ianv(nkndim)
	integer iaev(neldim)
	integer ialv(nlidim)
	real hnv(nkndim)
	real hev(neldim)
	real hlv(nlidim)
	real xgv(nkndim),ygv(nkndim)
	integer nen3v(3,neldim)
	integer ipntlv(0:nlidim)
	integer inodlv(nlndim)

	integer ipaux(nkndim)
	integer i,j,k,low,high

	call getfnm('lagra',file)

        call rdgrd(
     +                   file
     +                  ,bstop
     +                  ,nco,nkn,nel,nli
     +                  ,nkndim,neldim,nlidim,nlndim
     +                  ,ipnv,ipev,iplv
     +                  ,ianv,iaev,ialv
     +                  ,hnv,hev,hlv
     +                  ,xgv,ygv
     +                  ,nen3v
     +                  ,ipntlv,inodlv
     +                  )

	if( bstop ) goto 99

	call ex2in(nkn,nel,nli,ipnv,ipaux,nen3v,inodlv,bstop)
	if( bstop ) goto 98

	nline = nli
	ibnd(i) = 0
	do i=1,nline
	  ibnd(i) = ipntlv(i)
	  low = ibnd(i-1) + 1
	  high = ibnd(i)
	  do j=low,high
	    k = inodlv(j)
	    xbnd(j) = xgv(k)
	    ybnd(j) = ygv(k)
	  end do
	end do

	iline = 0
	
	return
   98	continue
	stop 'error stop laggrd_read_file: error converting node numbers'
   99	continue
	stop 'error stop laggrd_read_file: error reading grd file'
	end

c*******************************************************************

	function laggrd_read_next(ndim,n,x,y)

c this routine is never used

	implicit none

	logical laggrd_read_next
	integer ndim,n
	real x(1),y(1)

	include 'param.h'

c	--------------------------------------
	integer nlidim,nlndim
	parameter(nlidim=100,nlndim=nkndim)
	real xbnd(nlndim), ybnd(nlndim)
	common /xybnd/xbnd,ybnd
	integer nline,iline
	integer ibnd(0:nlidim)
	common /ibnd/nline,iline,ibnd
	save /xybnd/,/ibnd/
c	--------------------------------------

	integer i,j,low,high

	iline = iline + 1

	low = ibnd(iline-1) + 1
	high = ibnd(iline)
	n = high - low + 1
	if( n .gt. ndim ) stop 'error stop laggrd_read_next: ndim'

	i = 0
	do j=low,high
	  i = i + 1
	  x(i) = xbnd(j)
	  y(i) = xbnd(j)
	end do

	if( x(1) .eq. x(n) .and. y(1) .eq. y(n) ) n = n - 1	!closed line

	laggrd_read_next = .true.

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine lagbound_open_file(iunit)

	implicit none

	integer iunit
	character*80 line
	integer ifileo

	call getfnm('lagra',line)
        iunit=ifileo(0,line,'form','old')

	if( iunit .le. 0 ) stop 'error stop open_lagbound_file'

	end

c*******************************************************************

	subroutine lagbound_close_file(iunit)

	implicit none

	integer iunit

	close(iunit)

	end

c*******************************************************************

	function lagbound_read_next(iunit,ndim,n,x,y)

c reads polygons from file

	implicit none

	logical lagbound_read_next
	integer iunit
	integer ndim,n
	real x(1),y(1)

	integer iflag
	real xa,ya

	integer isave
	real xsave,ysave
	save isave,xsave,ysave
	data isave /0/

	n = 0
	iflag = 0
	lagbound_read_next = .true.

	if( isave .eq. -1 ) then	!EOF already encountered
	  lagbound_read_next = .false.
	  isave = 0			!next read would succeed
	  return
	else if( isave .eq. 1 ) then	!first point already read
	  n = 1
	  x(n) = xsave
	  y(n) = ysave
	end if

    1	continue
	  read(iunit,*,end=2) xa,ya,iflag
	  if( iflag .ne. 0 .and. isave .gt. 0 ) goto 2	!not 1.point on 1.line
	  n = n + 1
	  if( n .gt. ndim ) stop 'error stop lagbound_read_next: ndim'
	  x(n) = xa
	  y(n) = ya
	  goto 1
    2	continue

	if( iflag .eq. 1 ) then	!start of new line -> save for next call
	  xsave = xa
	  ysave = ya
	  isave = 1
	else			!EOF
	  isave = -1
	end if

	if( n .gt. 0 ) then
	  if( x(1) .eq. x(n) .and. y(1) .eq. y(n) ) n = n - 1	!closed line
	end if

	if( n .lt. 3 ) then
	  lagbound_read_next = .false.
	end if

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine check_elements(n,x,y,iflag)

c flags elements which are in/out-side or at border of line

	implicit none

	integer n
	real x(1),y(1)
	integer iflag(1)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nen3v(3,1)
	common /nen3v/nen3v
	real xgv(1)
	common /xgv/xgv
	real ygv(1)
	common /ygv/ygv
	real v1v(1)
	common /v1v/v1v

	logical inpoly
	logical debug
	integer ie,k,ii,iin
	integer iout,icheck
	real xlmin,xlmax,ylmin,ylmax
	real xmin,xmax,ymin,ymax
	real xa,ya,xe(3),ye(3)

	debug = .false.

c------------------------------------------------------------
c if n < 0 flag all nodes as inside
c------------------------------------------------------------

	if( n .lt. 0 ) then
	  do ie=1,nel
	    iflag(ie) = 1
	  end do
	  return
	end if

c------------------------------------------------------------
c compute min/max of line
c------------------------------------------------------------

	call xy_minmax(n,x,y,xlmin,xlmax,ylmin,ylmax)

	if( debug ) then
	  write(6,*) 'in check_elements ',n
	  do ii=1,n
	    write(6,*) ii,x(ii),y(ii)
	  end do
	end if

	write(6,*) 'min/max of line: ',xlmin,xlmax,ylmin,ylmax

c------------------------------------------------------------
c flag nodes that are inside line (v1v(k) = 1)
c------------------------------------------------------------

	do k=1,nkn
	  v1v(k) = 0.

	  xa = xgv(k)
	  ya = ygv(k)

	  if( xa .lt. xlmin ) goto 1
	  if( xa .gt. xlmax ) goto 1
	  if( ya .lt. ylmin ) goto 1
	  if( ya .gt. ylmax ) goto 1

	  if( inpoly(n,x,y,xa,ya) ) v1v(k) = 1.

    1	  continue
	end do

	iin = 0
	do k=1,nkn
	  if( v1v(k) .ne. 0. ) iin = iin + 1
	end do
	write(6,*) 'nodes inside: ',iin

c------------------------------------------------------------
c flag elements that are inside (1), outside (-1), or on border (0)
c------------------------------------------------------------

	do ie=1,nel

	  call xy_and_minmax(ie,xe,ye,xmin,xmax,ymin,ymax)	!xe,ye not used

	  iflag(ie) = -1

	  if( xmax .lt. xlmin ) goto 2
	  if( xmin .gt. xlmax ) goto 2
	  if( ymax .lt. ylmin ) goto 2
	  if( ymin .gt. ylmax ) goto 2

	  iin = 0
	  do ii=1,3
	    k = nen3v(ii,ie)
	    if( v1v(k) .ne. 0. ) iin = iin + 1
	  end do

	  !if( iin .ne. 0 ) write(6,*) 'element: ',ie,iin

	  if( iin .eq. 0 ) then			!no vertex inside line
	    iflag(ie) = -1
	  else if( iin .eq. 3 ) then		!all vertices inside line
	    iflag(ie) = 1
	  else
	    iflag(ie) = 0
	  end if

c         ---------------------------------------------------
c	  if iflag != 0 we should still check if no line point is in element
c         ---------------------------------------------------

	  if( iflag(ie) .ne. 0 ) then
	    call linepoint_in_element(n,x,y,xe,ye,iin)
	    if( iin .gt. 0 ) iflag(ie) = 0
	  end if

    2	  continue
	end do

c------------------------------------------------------------
c final statistics
c------------------------------------------------------------

	iin = 0
	iout = 0
	icheck = 0
	do ie=1,nel
	  if( iflag(ie) .eq. -1 ) iout = iout + 1
	  if( iflag(ie) .eq.  0 ) icheck = icheck + 1
	  if( iflag(ie) .eq. +1 ) iin = iin + 1
	end do

	write(6,*) 'total number of elements:           ',nel
	write(6,*) 'elements containing no line points: ',iout
	write(6,*) 'elements fully in line:             ',iin
	write(6,*) 'elements to be checked:             ',icheck

c------------------------------------------------------------
c end of routine
c------------------------------------------------------------

	end

c*******************************************************************

	subroutine linepoint_in_element(n,x,y,xe,ye,inside)

c checks if one of the line points is in element given by xe,ye

	implicit none

	integer n
	real x(1),y(1)
	real xe(3),ye(3)
	integer inside		!0: no point inside (return)

	integer i
	integer intri

	inside = 1

	do i=1,n
          if( intri(xe,ye,x(i),y(i)) .gt. 0 ) return
	end do

	inside = 0	!no point inside found

	end

c*******************************************************************

	subroutine lagbound_seed_particles(dxy,n,x,y)

c release in partial area

	implicit none

	include 'param.h'

	real dxy
	integer n
	real x(1),y(1)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	integer iflag(neldim)

	integer ie,i,j
	integer nin,iin,ifl,np
	integer imin,imax,jmin,jmax
	real xmin,xmax,ymin,ymax
	real x0,y0
	real xp,yp
	real xe(3),ye(3)
	integer nc(-1:1)

	integer intri
	logical inpoly

c-----------------------------------------------------------------
c get central node (x0,y0)
c-----------------------------------------------------------------

	if( n .lt. 0 ) then	!total domain
	  call basin_center(x0,y0)
	else
	  call xy_center(n,x,y,x0,y0)
	end if

c-----------------------------------------------------------------
c flag elements as internal (1), external (-1) and on border (0)
c-----------------------------------------------------------------

	call check_elements(n,x,y,iflag)

c-----------------------------------------------------------------
c some statistics
c-----------------------------------------------------------------

	do i=-1,1
	  nc(i) = 0
	end do

	do ie=1,nel
	  nin = iflag(ie)
	  if( nin .lt. -1 .or. nin .gt. 1 ) then
	    stop 'error stop lagbound_seed_particles'
	  end if
	  nc(nin) = nc(nin) + 1
	end do

	write(6,*) 'lagbound_seed_particles: ',nel,n,nc

c-----------------------------------------------------------------
c loop on elements
c-----------------------------------------------------------------

	np = 0
        do ie=1,nel

	  ifl = iflag(ie)
	  if( ifl .eq. -1 ) goto 1

          call xy_and_minmax(ie,xe,ye,xmin,xmax,ymin,ymax)

          imin=(xmin-x0)/dxy + 0
          imax=(xmax-x0)/dxy + 2
          jmin=(ymin-y0)/dxy + 0
          jmax=(ymax-y0)/dxy + 2

          do i=imin,imax
            do j=jmin,jmax
              xp=i*dxy+x0
              yp=j*dxy+y0

              iin=intri(xe,ye,xp,yp)

              if( iin .ne. 0 ) then
		if( ifl .eq. 1 .or. inpoly(n,x,y,xp,yp) ) then
                  np = np + 1
                  call insert_particle(ie,0.,xp,yp)
		end if
              end if
            end do
          end do

    1	  continue
        end do

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

        write(6,*) 'lagr_init_all: ',np,' particles inserted'

	end

c*******************************************************************

	subroutine lgr_init_all(dxy)

c release in total area

	implicit none

	real dxy

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	integer ie,i,j
	integer nin,iin
	integer imin,imax,jmin,jmax
	real xmin,xmax,ymin,ymax
	real x0,y0
	real xp,yp
	real xe(3),ye(3)

	integer intri

c------------------------------------------------------------------
c compute regular mesh and insert particles
c------------------------------------------------------------------

	nin = 0
	call basin_center(x0,y0)

	do ie=1,nel

	  call xy_and_minmax(ie,xe,ye,xmin,xmax,ymin,ymax)

          imin=(xmin-x0)/dxy + 0
          imax=(xmax-x0)/dxy + 2
          jmin=(ymin-y0)/dxy + 0
          jmax=(ymax-y0)/dxy + 2

	  do i=imin,imax
            do j=jmin,jmax
              xp=i*dxy+x0
              yp=j*dxy+y0

              iin=intri(xe,ye,xp,yp)

	      if( iin .ne. 0 ) then
		nin = nin + 1
	  	call insert_particle(ie,0.,xp,yp)
	      end if
	    end do
	  end do

	end do

	write(6,*) 'lagr_init_all: ',nin,' particles inserted'

c------------------------------------------------------------------
c end of routine
c------------------------------------------------------------------

	end

c*******************************************************************

