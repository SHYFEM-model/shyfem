c
c $Id: testdash.f,v 1.1 2004/09/21 08:53:24 georg Exp $
c
c revision log :
c
c 13.04.1999	ggu	written from scratch
c
c*************************************************************

	module mod_testdash

	implicit none

        real, save :: sigma_tot = 1.3
        real, save :: sigma_r = 1.0
        real, save :: d_sigma = 0.2
        integer, save :: ndec = 1

	end module mod_testdash

c*************************************************************

	program testdash

c plots simulation

	implicit none

c parameters
	integer ndim
	parameter (ndim=40)

	real wx(ndim)
	real wy(ndim)

	integer n
	real rmax

	n = 1
	rmax = 20.

c-----------------------
	call qopen
c-----------------------

        call taylor_start
	call taylor_set(9.,5.5,2.,-1)
        call taylor_full
        call taylor_data(8.,0.9,'A')
        call taylor_data(6.,0.6,'B')
        call taylor_end

c-----------------------
	call qclose
c-----------------------

	end

c*****************************************************************

        subroutine taylor_start

        implicit none

        call qstart

        call qfont('Times-Roman')
        call qgray(0.)
        call qtxts(10)

        end

c*****************************************************************

        subroutine taylor_end

        implicit none

        call qend

        end

c*****************************************************************

        subroutine taylor_set(s_tot,s_r,d_s,nd)

	use mod_testdash

        implicit none

        real s_tot      !maximum of sigma to plot
        real s_r        !sigma of reference point
        real d_s        !plot marker every d_s sigma
        integer nd      !decimals for sigma value

        sigma_tot = s_tot
          sigma_r = s_r
          d_sigma = d_s
             ndec = nd

        end

c*****************************************************************

        subroutine taylor_get(s_tot,s_r,d_s,nd)

	use mod_testdash

        implicit none

        real s_tot      !maximum of sigma to plot
        real s_r        !sigma of reference point
        real d_s        !plot marker every d_s sigma
        integer nd      !decimals for sigma value

        s_tot = sigma_tot
          s_r = sigma_r
          d_s = d_sigma
           nd = ndec

        end

c*****************************************************************

        subroutine taylor_data(s_f,r,text)

	use mod_testdash

        implicit none

        real s_f                !sigma of time series
        real r                  !correlation
        character*(*) text      !text to write

        real beta,x,y
        real width,height

        beta = acos(r)
        x = s_f * cos(beta)
        y = s_f * sin(beta)

        call qtsize(text,width,height)
        x = x - 0.5*width
        y = y - 0.5*height

        call qtext(x,y,text)
        
        end

c*****************************************************************

	subroutine taylor_full

	implicit none

        real rad,dr,sigmar
        integer ndec

        real rmax,rmin
        real v0,dvx,dvy
        real alpha
        real r
        real pi,aux,cor,beta,x,y
        real cosphi,phi
        real width,height,dbeta
        real aleg,phileg
        real tick,x1,y1
        integer n,j,ndig
        character*30 text

        integer ialfa

        integer ncor,nticks
        parameter(ncor=8,nticks=11)
        real corcof(ncor)
        real cortik(nticks)
        data corcof /0.2,0.4,0.6,0.8,0.9,0.95,0.99,1.0/
        data cortik /0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.99/

c-------------------------------------------------------
c set up
c-------------------------------------------------------

        pi = 4. * atan(1.)
        aux = 180./pi
        alpha = 1.3
        aleg = 3./4.
        tick = 0.01 * rad

        call taylor_get(rad,sigmar,dr,ndec)

        rmax = alpha * rad
        rmin = rmax - rad

        v0 = 1.
        dvx = 18.
        dvy = dvx * (rmax+rmin) / (2.*rmax)
	call qsetvp(v0,v0,v0+dvx,v0+dvy)
	call qworld(-rmax,-rmin,rmax,rmax)

c-------------------------------------------------------
c domain
c-------------------------------------------------------

	call qline(-rad,0.,rad,0.)
	call qline(0.,0.,0.,rad)
	call qarc(0.,0.,rad,0.,180.);

c-------------------------------------------------------
c standard deviation
c-------------------------------------------------------

        call setdash(4,0.02)
        n = 0.1 + rad/dr

        do j=1,n
          r = j * dr
	  call qarc(0.,0.,r,0.,180.);

          call qline(r,0.,r,tick)
          call qline(-r,0.,-r,tick)

          ndig = ialfa(r,text,ndec,-1)
          call qtsize(text,width,height)
          y = 0. - 1.8 * height
          x = r - 0.5 * width
          call qtext(x,y,text)
          x = - r - 0.5 * width
          call qtext(x,y,text)
        end do

        r = 0.
        ndig = ialfa(r,text,ndec,-1)
        call qtsize(text,width,height)
        y = 0. - 1.8 * height
        x = r - 0.5 * width
        call qtext(x,y,text)

c       ticks on axis

        call setdash(0,0.)
        do j=1,n
          r = j * dr
          call qline(r,0.,r,tick)
          call qline(-r,0.,-r,tick)
        end do

c-------------------------------------------------------
c correlation coefficient
c-------------------------------------------------------

        dbeta = 1./aux              !for starting point
        call setdash(4,0.02)
        do j=1,ncor
          cor = corcof(j)
          beta = acos(cor)
          x = rad*cos(beta)
          y = rad*sin(beta)
          write(6,*) cor,beta,beta*aux,x,y
	  call qline(0.,0.,x,y)
	  call qline(0.,0.,-x,y)

          call make_text(cor,text)
          call qtxtr(aux*beta)
          x = 1.03*rad*cos(beta-dbeta)
          y = 1.03*rad*sin(beta-dbeta)
          call qtext(x,y,text)

          call make_text(-cor,text)
          call qtsize(text,width,height)
          call qtxtr(-aux*beta)
          x = (1.03*rad+width)*cos(beta-dbeta)
          y = (1.03*rad+width)*sin(beta-dbeta)
          call qtext(-x,y,text)
        end do

        call make_text(0.0,text)
        call qtxtr(90.)
        beta = 0.5 * pi
        x = 1.03*rad*cos(beta-dbeta)
        y = 1.03*rad*sin(beta-dbeta)
        call qtext(x,y,text)

        call qtxtr(0.)

c       ticks on axis

        call setdash(0,0.)
        do j=1,nticks
          cor = cortik(j)
          r = j * dr
          beta = acos(cor)
          x = rad*cos(beta)
          y = rad*sin(beta)
          x1 = (rad-tick)*cos(beta)
          y1 = (rad-tick)*sin(beta)
          call qline(x,y,x1,y1)
          call qline(-x,y,-x1,y1)
        end do

c-------------------------------------------------------
c rms pattern
c-------------------------------------------------------

        call setdash(5,0.02)
        n = 0.1 + (rad+sigmar)/dr
        do j=1,n
          r = j * dr
          cosphi = (rad**2 - r**2 - sigmar**2) / (2.*r*sigmar)
          if( cosphi .ge. 1. ) then
            phi = 0.
          else
            phi = aux * acos(cosphi)
          end if
	  call qarc(sigmar,0.,r,phi,180.);

          phileg = aleg*180. + (1.-aleg)*phi
          x = (0.05*dr+r)*cos(phileg/aux) + sigmar
          y = (0.05*dr+r)*sin(phileg/aux)
          ndig = ialfa(r,text,ndec,-1)
          call qtxtr(phileg-90.)
          call qtext(x,y,text)
          write(6,*) j,cosphi,phi,phileg
        end do

        call qtxtr(0.)

c-------------------------------------------------------
c sigma r (reference point)
c-------------------------------------------------------

        call setdash(0,0.)
	call qarc(sigmar,0.,rad/100.,0.,360.);

c-------------------------------------------------------
c title
c-------------------------------------------------------

        text = 'Correlation Coefficient'
        call qtsize(text,width,height)
        x = -0.5 * width
        y = 1.15 * rad
        call qtext(x,y,text)

        text = 'Standard Deviation'
        call qtsize(text,width,height)
        x = -0.5 * width
        y = -4.0 * height
        call qtext(x,y,text)

c-------------------------------------------------------
c end of routine
c-------------------------------------------------------

	end

c*****************************************************************

        subroutine make_text(val,text)

        implicit none

        real val
        character*(*) text

        integer space,comma
        character*6 format

        space=3
        comma=1
        if( val .lt. 0. ) space = space + 1
        if( abs(val) .gt. 0.9 ) then
          space = space + 1
          comma = comma + 1
        end if

        write(format,'(a2,i1,a1,i1,a1)') '(f',space,'.',comma,')'
        write(text,format) val

        end

c*****************************************************************

	subroutine setdash(i,fact)

	implicit none

        integer i
        real fact

        real dash1(2)
        real dash2(2)
        real dash3(2)
        real dash4(2)
        real dash5(2)

        data dash1 /3.,3./
        data dash2 /1.,1./
        data dash3 /3.,1./
        data dash4 /1.,4./
        data dash5 /5.,5./

        if( i .eq. 0 ) then
          call qdash0
        else if( i .eq. 1 ) then
          call qdash(fact,0.,2,dash1)
        else if( i .eq. 2 ) then
          call qdash(fact,0.,2,dash2)
        else if( i .eq. 3 ) then
          call qdash(fact,0.,2,dash3)
        else if( i .eq. 4 ) then
          call qdash(fact,0.,2,dash4)
        else if( i .eq. 5 ) then
          call qdash(fact,0.,2,dash5)
        else
          call qdash0
        end if

        end

c*****************************************************************

