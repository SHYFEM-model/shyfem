!
! $Id: subproj.f,v 1.5 2010-03-11 15:36:38 georg Exp $
!
! subroutines for handling projection
!
! revision log :
!
! 17.11.2011    ggu     written from scratch
! 10.01.2012    ggu     bug fix: c_param was real
! 21.01.2015    ggu     code to handle projection in both directions
! 30.09.2015    ccf     introduced module coordinates
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

        implicit none

        integer 			:: mode,i
        double precision, dimension(9)	:: c_param
        double precision 		:: c_lat0,c_lon0,c_phi
        real 				:: xmin,ymin,xmax,ymax

        mode = -1		!from lat/lon to cartesian
        iproj = 3		!always use cpp

        xmin = MINVAL(xgv)
        ymin = MINVAL(ygv)
        xmax = MAXVAL(xgv)
        ymax = MAXVAL(ygv)

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

        write(6,*) 'end of handle_projection'

        end subroutine handle_projection

!****************************************************************
