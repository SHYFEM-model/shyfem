c
c $Id: suputi.f,v 1.5 2004/03/09 17:28:09 georg Exp $
c
c general utilities for plot routines
c revision log :
c
c 16.10.2001    ggu     new routine fpfeil()
c 03.09.2003    ggu     bug fix XBXT in fpfeil
c 05.10.2003    ggu     bug fix S_IS_0 in pfeil
c
c**************************************************

	subroutine pbox(xmin,ymin,xmax,ymax)

c plots box around rectangle with actual color

	implicit none

	real xmin,ymin,xmax,ymax

        call qline(xmin,ymin,xmin,ymax)
        call qline(xmin,ymax,xmax,ymax)
        call qline(xmax,ymax,xmax,ymin)
        call qline(xmax,ymin,xmin,ymin)

	end

c*************************************************************

        subroutine pwbox(mode,xmin,ymin,xmax,ymax)

c plots box with white background (erasing)
c
c mode != 0 -> plot outline

        implicit none

        integer mode
        real xmin,ymin,xmax,ymax

        call qwhite(1)
        call qgray(1.)
        call qrfill(xmin,ymin,xmax,ymax)
        call qwhite(0)
        call qgray(0.)

        if( mode .ne. 0 ) call pbox(xmin,ymin,xmax,ymax)

        end

c**************************************************

        subroutine pfeil(x,y,u,v,s)

c draws arrow
c
c x,y           coordinates of arrow base
c u,v           size of arrow
c s             scaling of vector

        implicit none

        real dl,cos,sin
        parameter (dl=0.1)
        parameter (cos=.866025,sin=.5)

        real x,y,u,v,s
        real xu,yv,d
        real x1,y1,x2,y2

        if(u.eq.0..and.v.eq.0.) return
        if(s.eq.0.) return                      !S_IS_0

        xu=x+u*s
        yv=y+v*s

        d=0.2*s

        x1=xu-d*(cos*u+sin*v)
        y1=yv-d*(-sin*u+cos*v)
        x2=xu-d*(cos*u-sin*v)
        y2=yv-d*(sin*u+cos*v)

        call qmove(x,y)
        call qplot(xu,yv)
        call qplot(x1,y1)
        call qplot(x2,y2)
        call qplot(xu,yv)

        end

c**************************************************

        subroutine fpfeil(xb,yb,xt,yt,s)

c draws full arrow
c
c xb,yb         coordinates of arrow base
c xt,yt         coordinates of arrow tip
c s             size of tip (in fraction, 1 -> whole length)

        implicit none

        real xb,yb,xt,yt
        real s

        real cos,sin
        parameter (cos=.866025,sin=.5)  !cos(30), sin(30)

        real u,v
        real d
        real xs,ys
        real x1,y1,x2,y2
        real x3(3),y3(3)

        if( xb .eq. xt .and. yb .eq. yt ) return	!bug fix XBXT

        u = xt - xb
        v = yt - yb

        d=s

        x1=xt-d*(cos*u+sin*v)
        y1=yt-d*(-sin*u+cos*v)
        x2=xt-d*(cos*u-sin*v)
        y2=yt-d*(sin*u+cos*v)

        xs = 0.5 * ( x1 + x2 )   !this is the point on the bottom of the
        ys = 0.5 * ( y1 + y2 )   !...arrow, lying on the line

        call qmove(xb,yb)
        !call qplot(xt,yt)      !if line is too thick, this is not nice
        call qplot(xs,ys)

        x3(1) = xt
        y3(1) = yt
        x3(2) = x1
        y3(2) = y1
        x3(3) = x2
        y3(3) = y2

        call qafill(3,x3,y3)

        !call qplot(x1,y1)
        !call qplot(x2,y2)
        !call qplot(xt,yt)

        end

c**************************************************

