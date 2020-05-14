
!--------------------------------------------------------------------------
!
!    Copyright (C) 2007,2009-2012,2014-2019  Georg Umgiesser
!    Copyright (C) 2015  Christian Ferrarin
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

c initialization routines
c
c revision log :
c
c 12.12.2007	ggu	written from scratch
c 11.09.2009	ggu	more checks for initial seeding (linepoint_in_element)
c 11.03.2010	ggu	small bug fix in lagbound_seed_particles() -> init np
c 11.03.2010	ggu	in lgr_init_line() check if at least one line is read
c 23.03.2010	ggu	changed v6.1.1
c 06.10.2011	ggu	in check_elements() initialize iout, icheck (bug)
c 18.10.2011	ggu	changed VERS_6_1_33
c 24.01.2012	ggu	changed VERS_6_1_41
c 16.02.2012	ggu	bug fix reading lines
c 23.12.2014	ggu	changed VERS_7_0_11
c 19.01.2015	ggu	changed VERS_7_1_3
c 23.04.2015	ggu	changed VERS_7_1_8
c 20.05.2015	ccf	give type to released particles
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 23.09.2015	ggu	changed VERS_7_2_4
c 16.11.2015	ggu	changed VERS_7_3_14
c 15.02.2016	ggu	handle negative areas (line given clockwise)
c 01.04.2016	ggu	changed VERS_7_5_7
c 31.03.2017	ggu	changed VERS_7_5_24
c 25.10.2018	ggu	changed VERS_7_5_51
c 16.02.2019	ggu	changed VERS_7_5_60
c
c*******************************************************************

	subroutine lgr_init_shell

c manages release of particles

	implicit none

	character*80 line,traj
	integer npoints
	real getpar

	npoints = nint(getpar('nbdy'))
	if( npoints <= 0 ) return

        call getfnm('lgrlin',line)
        call getfnm('lgrtrj',traj)

        if( traj .ne. ' ' .and. line .ne. ' ') then
          write(*,*)'lgrlin and lgrtrj not allowed'
          stop 'error stop lgr_init_shell'
        end if

        if( traj .ne. ' ' ) return

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
c particles in different polygons have different types

	implicit none

	integer npoints

	integer ndim
	parameter (ndim=100)

	real x(ndim),y(ndim)

	integer iunit,n,nline
	real tot_area,dxy
	logical lagbound_read_next

	call lagbound_get_area(tot_area)
	if( tot_area < 0. ) then
	  write(6,*) 'area inside line: ',tot_area
	  stop 'error stop lgr_init_line: negative area'
	end if
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

c*******************************************************************

	subroutine lgr_init_traj(dtime)

c manages release of particles at specific time and positions
c particles in different points have different types

        use mod_lagrange

	implicit none

	double precision, intent(in)	:: dtime

        character*20 			:: line
        integer				:: date,time
	real				:: x,y
	integer				:: iunit,ierr,ierrp
	integer				:: ie,n,i
        double precision, save		:: atime,atime0

        integer, save                   :: icall = 0

	integer, save        		:: npoi,nt,ntot
	double precision, allocatable, save	:: dtimet(:)
        real, allocatable, save		:: xt(:),yt(:)
	double precision, save		:: dtlast
        character*80 traj
        integer ifileo
        real getpar

        if( icall == -1 ) return

!---------------------------------------------------------------
! initialization: read file and store time and position
!---------------------------------------------------------------

        if( icall == 0 ) then
          call getfnm('lgrtrj',traj)
          if( traj .eq. ' ' ) icall = -1
          if( icall .eq. -1 ) return

          npoi = nint(getpar('nbdy'))
          ntot = npoi*max(1,abs(ipvert))

          nt = 0
          ierr = 0
          iunit=ifileo(0,traj,'form','old')
          do while( ierr .eq. 0 )
            read(iunit,*,iostat=ierr,end=98)
            if( ierr /= 0 ) exit
            nt = nt + 1
          end do
 98       close(iunit)
          allocate(dtimet(nt))
          allocate(xt(nt))
          allocate(yt(nt))

          call get_absolute_ref_time(atime0)
          nt = 0
          ierr = 0
          iunit=ifileo(0,traj,'form','old')
          do while( ierr .eq. 0 )
            read(iunit,*,iostat=ierr,end=99) line,x,y
            if( ierr /= 0 ) exit
            nt = nt + 1
            call dts_string2time(line,atime,ierrp)
            dtimet(nt) = atime - atime0
            xt(nt) = x
            yt(nt) = y
          end do
 99       close(iunit)

          dtlast = dtimet(nt)

          icall = 1
        end if

!---------------------------------------------------------------
! normal call
!---------------------------------------------------------------

        if ( dtime > dtlast ) return

        do i = 1,nt
          if ( dtimet(i) == dtime ) then
            write(6,*) 'release of particles for lagrangian model'
            x = xt(i)
            y = yt(i)
            call find_element(x,y,ie)
            do n = 1,npoi
              call insert_particle_3d(ie,i,0.,x,y)
            end do
            call get_timeline(dtime,line)
            write(6,*) 'new particles released: ',ntot,'at time: ',line
          end if
        end do

        end subroutine lgr_init_traj	

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
	  if( area < 0. ) area = -area	!in case line is given clockwise
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

	subroutine lagbound_open_file(iunit)

	implicit none

	integer iunit
	character*80 line
	integer ifileo

	call getfnm('lgrlin',line)
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
	real x(ndim),y(ndim)

	integer iflag
	real xa,ya

	integer isave,iline
	real xsave,ysave
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

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine lagbound_seed_particles(dxy,iln,n,x,y)

c release in partial area

	use mod_lagrange
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	real dxy
	integer iln		!line number for particle type
	integer n
	real x(n),y(n)

	integer ieflag(nel)
	integer ikflag(nkn)

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

	call check_elements(n,x,y,ieflag,ikflag)

c-----------------------------------------------------------------
c some statistics
c-----------------------------------------------------------------

	do i=-1,1
	  nc(i) = 0
	end do

	do ie=1,nel
	  nin = ieflag(ie)
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

	  ifl = ieflag(ie)
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
                  call insert_particle_3d(ie,iln,0.,xp,yp)
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

	use basin, only : nkn,nel,ngr,mbw

	implicit none

	real dxy

	integer ie,i,j
	integer nin,iin
	integer imin,imax,jmin,jmax
	real xmin,xmax,ymin,ymax
	real x0,y0
	real xp,yp
	real xe(3),ye(3)

	integer intri
	integer, parameter :: ity=0	!particle type = 0 everywere

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
	  	call insert_particle_3d(ie,ity,0.,xp,yp)
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

