!
! $Id: lagrange_init.f,v 1.8 2010-03-11 15:36:38 georg Exp $
!
! initialization routines
!
! revision log :
!
! 12.12.2007    ggu	written from scratch
! 11.09.2009    ggu	more checks for initial seeding (linepoint_in_element)
! 11.03.2010    ggu	small bug fix in lagbound_seed_particles() -> init np
! 11.03.2010    ggu	in lgr_init_line() check if at least one line is read
! 06.10.2011    ggu	in check_elements() initialize iout, icheck (bug)
! 16.02.2012    ggu	bug fix reading lines
! 20.05.2015    ccf	give type to released particles
!
!*******************************************************************
!----------------------------------------------------------------------------
        module lagrange_init
!----------------------------------------------------------------------------
        contains
!----------------------------------------------------------------------------

	subroutine lgr_init_shell

! manages release of particles

        use para

	implicit none

	character*80 line
	integer npoints

	npoints = nint(getpar('nbdy'))
	if( npoints <= 0 ) return

        call getfnm('lagra',line)

        if( line .eq. ' ' ) then      !total lagoon
            call lgr_init_total(npoints)
        else                          !use boundary
            call lgr_init_line(npoints)
        end if

	end

!*******************************************************************

	subroutine lgr_init_total(npoints)

! manages release of particles over whole domain

	use lagrange_util

	implicit none

	integer npoints

	integer n
	double precision tot_area,dxy
	double precision x(1),y(1)

	call compute_total_area(tot_area)
	dxy = sqrt( tot_area / npoints )

	n = -1
	!call lagbound_seed_particles(dxy,n,x,y)	!should be the same
	call lgr_init_all(dxy)

	end

!*******************************************************************

	subroutine lgr_init_line(npoints)

! manages release of particles inside polygon(s)
! particles in different polygons have different types

	implicit none

	integer npoints

	integer ndim
	parameter (ndim=100)

	double precision x(ndim),y(ndim)

	integer iunit,n,nline
	double precision tot_area,dxy

	call lagbound_get_area(tot_area)
        dxy = sqrt( tot_area / npoints )

	write(6,*) 'lgr_init_line: ',npoints,tot_area,dxy

	call lagbound_open_file(iunit)
	nline = 0
	do while( lagbound_read_next(iunit,ndim,n,x,y) )
	  nline = nline + 1
	  write(6,*) 'new line read: ',nline,n
	  call lagbound_seed_particles(dxy,nline,n,x,y)
	end do
	call lagbound_close_file(iunit)
	write(6,*) 'total number of lines read: ',nline

	if( nline .eq. 0 ) then
	  write(6,*) 'no line found in file'
	  write(6,*) 'cannot initialize particles'
	  stop 'error stop lgr_init_line: no line'
	end if

	end

!*******************************************************************
!*******************************************************************
!*******************************************************************

	subroutine lagbound_get_area(arealine)

! computes area of polygon(s)

        use geometrical

	implicit none

	double precision arealine

	integer ndim
	parameter (ndim=100)

	double precision x(ndim),y(ndim)

	logical debug
	integer iunit,n,i
	double precision area,areatot

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

!*******************************************************************
!*******************************************************************
!*******************************************************************

	subroutine lagbound_open_file(iunit)

        use fil
        use para

	implicit none

	integer iunit
	character*80 line

	call getfnm('lagra',line)
        iunit=ifileo(0,line,'form','old')

	if( iunit .le. 0 ) stop 'error stop open_lagbound_file'

	end

!*******************************************************************

	subroutine lagbound_close_file(iunit)

	implicit none

	integer iunit

	close(iunit)

	end

!*******************************************************************

	function lagbound_read_next(iunit,ndim,n,x,y)

! reads polygons from file

	implicit none

	logical lagbound_read_next
	integer iunit
	integer ndim,n
	double precision x(1),y(1)

	integer iflag
	double precision xa,ya

	integer isave,iline
	double precision xsave,ysave
	save isave,iline,xsave,ysave
	data isave,iline /0,0/

	n = 0
	iflag = 0
	lagbound_read_next = .true.

	if( isave .eq. -1 ) then	!EOF already encountered
	  lagbound_read_next = .false.
	  isave = 0			!next read would succeed
	  iline = 0
	  return
	else if( isave .eq. 1 ) then	!first point already read
	  n = 1
	  x(n) = xsave
	  y(n) = ysave
	end if

    1	continue
	  read(iunit,*,end=2) xa,ya,iflag
	  if( iflag .ne. 0 .and. n .gt. 0 ) goto 2	!not 1.point on 1.line
	  n = n + 1
	  if( n .gt. ndim ) stop 'error stop lagbound_read_next: ndim'
	  x(n) = xa
	  y(n) = ya
	  goto 1
    2	continue

	if( iflag .ne. 0 ) then	!start of new line -> save for next call
	  xsave = xa
	  ysave = ya
	  isave = 1
	else			!EOF
	  isave = -1
	end if

	if( n .gt. 1 ) then
	  if( x(1) .eq. x(n) .and. y(1) .eq. y(n) ) n = n - 1	!closed line
	end if

	iline = iline + 1

	if( n .gt. 0 .and. n .lt. 3 ) then
	  write(6,*) 'line has less then 3 vertices: ',iline,n
	  stop 'error stop lagbound_read_next: degenerated line'
	else if( n .eq. 0 ) then
	  lagbound_read_next = .false.
	end if

	end

!*******************************************************************
!*******************************************************************
!*******************************************************************

	subroutine check_elements(n,x,y,iflag)

! flags elements which are in/out-side or at border of line

	use basin
        use geometrical
	use lagrange_util

	implicit none

	integer n
	double precision x(1),y(1)
	integer iflag(1)

	include 'param.h'

	logical debug
	integer ie,k,ii,iin
	integer iout,icheck
	double precision xlmin,xlmax,ylmin,ylmax
	double precision xmin,xmax,ymin,ymax
	double precision xa,ya,xe(3),ye(3)
	double precision v1v(nkn)

	debug = .false.

!------------------------------------------------------------
! if n < 0 flag all nodes as inside
!------------------------------------------------------------

	if( n .lt. 0 ) then
	  do ie=1,nel
	    iflag(ie) = 1
	  end do
	  return
	end if

!------------------------------------------------------------
! compute min/max of line
!------------------------------------------------------------

	call xy_minmax(n,x,y,xlmin,xlmax,ylmin,ylmax)

	if( debug ) then
	  write(6,*) 'in check_elements ',n
	  do ii=1,n
	    write(6,*) ii,x(ii),y(ii)
	  end do
	end if

	write(6,*) 'min/max of line: ',xlmin,xlmax,ylmin,ylmax

!------------------------------------------------------------
! flag nodes that are inside line (v1v(k) = 1)
!------------------------------------------------------------

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

!------------------------------------------------------------
! flag elements that are inside (1), outside (-1), or on border (0)
!------------------------------------------------------------

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

!         ---------------------------------------------------
!	  if iflag != 0 we should still check if no line point is in element
!         ---------------------------------------------------

	  if( iflag(ie) .ne. 0 ) then
	    call linepoint_in_element(n,x,y,xe,ye,iin)
	    if( iin .gt. 0 ) iflag(ie) = 0
	  end if

    2	  continue
	end do

!------------------------------------------------------------
! final statistics
!------------------------------------------------------------

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

!------------------------------------------------------------
! end of routine
!------------------------------------------------------------

	end

!*******************************************************************

	subroutine linepoint_in_element(n,x,y,xe,ye,inside)

! checks if one of the line points is in element given by xe,ye

        use regular

	implicit none

	integer n
	double precision x(1),y(1)
	double precision xe(3),ye(3)
	integer inside		!0: no point inside (return)

	integer i

	inside = 1

	do i=1,n
          if( intri(xe,ye,x(i),y(i)) .gt. 0 ) return
	end do

	inside = 0	!no point inside found

	end

!*******************************************************************

	subroutine lagbound_seed_particles(dxy,iln,n,x,y)

! release in partial area

	use lagrange_data
	use lagrange_util
        use lagrange_inout
        use regular
        use geometrical
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	double precision dxy
	integer iln		!line number for particle type
	integer n
	double precision x(1),y(1)

	include 'param.h'

	integer iflag(nel)

	integer ie,i,j
	integer nin,iin,ifl,np
	integer imin,imax,jmin,jmax
	double precision xmin,xmax,ymin,ymax
	double precision x0,y0
	double precision xp,yp
	double precision xe(3),ye(3)
	integer nc(-1:1)

!-----------------------------------------------------------------
! get central node (x0,y0)
!-----------------------------------------------------------------

	if( n .lt. 0 ) then	!total domain
	  call basin_center(x0,y0)
	else
	  call xy_center(n,x,y,x0,y0)
	end if

!-----------------------------------------------------------------
! flag elements as internal (1), external (-1) and on border (0)
!-----------------------------------------------------------------

	call check_elements(n,x,y,iflag)

!-----------------------------------------------------------------
! some statistics
!-----------------------------------------------------------------

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

!-----------------------------------------------------------------
! loop on elements
!-----------------------------------------------------------------

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
                  call insert_particle_3d(ie,iln,0.d0,xp,yp)
		end if
              end if
            end do
          end do

    1	  continue
        end do

!-----------------------------------------------------------------
! end of routine
!-----------------------------------------------------------------

        write(6,*) 'lagr_init_all: ',np,' particles inserted'

	end

!*******************************************************************

	subroutine lgr_init_all(dxy)

! release in total area

        use regular
	use lagrange_util
        use lagrange_inout
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	double precision dxy

	!include 'param.h'
	!include 'lagrange.h'

	integer ie,i,j
	integer nin,iin
	integer imin,imax,jmin,jmax
	double precision xmin,xmax,ymin,ymax
	double precision x0,y0
	double precision xp,yp
	double precision xe(3),ye(3)

	integer, parameter :: ity=0	!particle type = 0 everywere

!------------------------------------------------------------------
! compute regular mesh and insert particles
!------------------------------------------------------------------

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
	  	call insert_particle_3d(ie,ity,0.d0,xp,yp)
	      end if
	    end do
	  end do

	end do

	write(6,*) 'lagr_init_all: ',nin,' particles inserted'

!------------------------------------------------------------------
! end of routine
!------------------------------------------------------------------

	end

!*******************************************************************

!----------------------------------------------------------------------------
        end module lagrange_init
!----------------------------------------------------------------------------
