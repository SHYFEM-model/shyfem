
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

! subroutines for handling projection
!
! revision log :
!
! 23.03.2010	ggu	changed v6.1.1
! 26.04.2010	ggu	changed VERS_6_1_6
! 26.04.2010	ggu	changed VERS_6_1_7
! 15.12.2010	ggu	changed VERS_6_1_14
! 17.02.2011	ggu	changed VERS_6_1_18
! 17.11.2011	ggu	written from scratch
! 10.01.2012	ggu	bug fix: c_param was real
! 13.10.2014	ggu	changed VERS_7_0_2
! 19.01.2015	ggu	changed VERS_7_1_3
! 21.01.2015	ggu	code to handle projection in both directions
! 30.09.2015	ccf	introduced module coordinates
! 10.10.2015	ggu	changed VERS_7_3_2
! 13.10.2015	ggu	changed VERS_7_3_5
! 05.11.2015	ggu	changed VERS_7_3_12
! 05.12.2017	ggu	changed VERS_7_5_39
! 16.02.2019	ggu	changed VERS_7_5_60
!
!****************************************************************

!==================================================================
        module coordinates
!==================================================================

        implicit none

        integer, private, save  :: nkn_proj = 0

        integer, save  :: iproj    ! call parameter to activate proj

        real, allocatable, save :: xgeov(:)
        real, allocatable, save :: ygeov(:)
        real, allocatable, save :: xcartv(:)
        real, allocatable, save :: ycartv(:)

!==================================================================
        contains
!==================================================================

        subroutine coordinates_init(nkn)

        integer :: nkn

        if( nkn == nkn_proj ) return

        if( nkn_proj > 0 ) then
            deallocate(xgeov)
            deallocate(ygeov)
            deallocate(xcartv)
            deallocate(ycartv)
        end if

        nkn_proj = nkn

        if( nkn == 0 ) return

        allocate(xgeov(nkn))
        allocate(ygeov(nkn))
        allocate(xcartv(nkn))
        allocate(ycartv(nkn))

        end subroutine coordinates_init

!************************************************************

        subroutine proj_cart2geo

! handles projection - converts x/y to lat/lon

        use basin

        implicit none

        integer 			:: mode,i
        double precision, dimension(9)	:: c_param
        real				:: getpar

        mode = 1	!from cartesian to lat/lon

        select case (iproj)
          case ( 0 )				!no projection
          case ( 1 )				!Guass-Boaga
            c_param(1) = getpar('c_fuse')
            c_param(2) = getpar('c_x0')
            c_param(3) = getpar('c_y0')
          case ( 2 )				!UTM
            c_param(1) = getpar('c_zone')
            c_param(2) = getpar('c_x0')
            c_param(3) = getpar('c_y0')
          case ( 3 )				!CPP
            c_param(1) = getpar('c_phi')
            c_param(2) = getpar('c_lon0')
            c_param(3) = getpar('c_lat0')
          case ( 4 )				!UTM non standard
            c_param(1) = getpar('c_lamb')
            c_param(2) = getpar('c_x0')
            c_param(3) = getpar('c_y0')
            c_param(4) = getpar('c_skal')
          case default
            write(6,*) 'iproj = ',iproj
            stop 'error stop proj_cart2geo: value for iproj not allowed'
        end select

        call init_coords(iproj,c_param)
        call convert_coords(mode,nkn,xgv,ygv,xgeov,ygeov)

        xcartv = xgv
        ycartv = ygv

        write(6,*) 'start of proj_cart2geo'
        write(6,*) 'mode  = ',mode
        write(6,*) 'iproj = ',iproj

        write(6,1000) (xgv(i),i=1,5)
        write(6,1000) (ygv(i),i=1,5)
        write(6,1000) (xgeov(i),i=1,5)
        write(6,1000) (ygeov(i),i=1,5)
        write(6,1000) (xcartv(i),i=1,5)
        write(6,1000) (ycartv(i),i=1,5)

        write(6,*) 'end of proj_cart2geo'

	return
1000	format(5g14.6)
        end subroutine proj_cart2geo

!****************************************************************

        subroutine proj_geo2cart

! handles projection - converts lat/lon to x/y

        use basin
        use shympi

        implicit none

        integer 			:: mode,i
        double precision, dimension(9)	:: c_param
        double precision 		:: c_lat0,c_lon0,c_phi
        real 				:: xmin,ymin,xmax,ymax

        mode = -1		!from lat/lon to cartesian
        iproj = 3		!always use cpp

	xmin = shympi_min(xgv)
	ymin = shympi_min(ygv)
	xmax = shympi_max(xgv)
	ymax = shympi_max(ygv)

	write(6,*) 'min/max: ',xmin,ymin,xmax,ymax

        c_phi  = 0.5*(ymax-ymin)
        c_lat0 = 0.5*(ymax-ymin)
        c_lon0 = 0.5*(xmax-xmin)

        c_param(1) = c_phi
        c_param(2) = c_lon0
        c_param(3) = c_lat0

        call init_coords(iproj,c_param)
        call convert_coords(mode,nkn,xcartv,ycartv,xgv,ygv)

        xgeov = xgv
        ygeov = ygv

        write(6,*) 'start of proj_geo2cart'
        write(6,*) 'mode  = ',mode
        write(6,*) 'iproj = ',iproj

        write(6,1000) (xgv(i),i=1,5)
        write(6,1000) (ygv(i),i=1,5)
        write(6,1000) (xgeov(i),i=1,5)
        write(6,1000) (ygeov(i),i=1,5)
        write(6,1000) (xcartv(i),i=1,5)
        write(6,1000) (ycartv(i),i=1,5)

        write(6,*) 'end of proj_geo2cart'

	return
1000	format(5g14.6)
        end subroutine proj_geo2cart

!****************************************************************

        subroutine baric_cart(ie,x,y)

! finds baricentre of element

! ie            number of element
! x,y           coordinates of baricentre (return value)

        use basin

        implicit none

        integer, intent(in)	:: ie
        real, intent(out)	:: x,y

        integer			:: ii,k
        real			:: xb,yb

        xb = 0.
        yb = 0.
        do ii=1,3
            k  = nen3v(ii,ie)
            xb = xb + xcartv(k)
            yb = yb + ycartv(k)
        end do

        x = xb/3.
        y = yb/3.

        end subroutine baric_cart

!****************************************************************

        subroutine getexy_cart(ie,x,y)

! gets coordinates x/y for element ie

        use basin

        implicit none

        integer, intent(in)	 	:: ie
        real, intent(out), dimension(3) :: x, y

        integer :: k,ii

        do ii=1,3
            k = nen3v(ii,ie)
            x(ii) = xcartv(k)
            y(ii) = ycartv(k)
        end do

        end subroutine getexy_cart

!==================================================================
        end module coordinates
!==================================================================

        subroutine handle_projection

! handles projection

        use coordinates
	use shympi

        implicit none

        logical 	:: bspheric
        logical 	:: is_spherical
	real		:: getpar

        bspheric = is_spherical()
        iproj = nint(getpar('iproj'))

        write(6,*) 'start of handle_projection'
	write(6,*) 'bspheric,iproj: ',bspheric,iproj

        if( bspheric ) then	!lat/lon -> cartesian
            call proj_geo2cart
        else			!cartesian -> lat/lon
            call proj_cart2geo
        end if

! check results - values must be equal over domains

	call shympi_check_2d_node(xgeov,'xgeov')
	call shympi_check_2d_node(ygeov,'ygeov')
	call shympi_check_2d_node(xcartv,'xcartv')
	call shympi_check_2d_node(ycartv,'ycartv')

        write(6,*) 'end of handle_projection'

        end subroutine handle_projection

!****************************************************************

