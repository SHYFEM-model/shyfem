
!--------------------------------------------------------------------------
!
!    Copyright (C) 2004,2010-2011,2019  Georg Umgiesser
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
c 05.04.2004	ggu	written from scratch
c 23.03.2010    ggu     changed v6.1.1
c 15.11.2011	ggu	new routine taylor_set_config() and taylor_set_param()
c 23.11.2011	ggu	rcol for colored lines
c 09.12.2011    ggu     changed VERS_6_1_38
c 14.02.2019    ggu     changed VERS_7_5_56
c
c notes :
c
c*************************************************************
c*****************************************************************
c*****************************************************************

        subroutine taylor_start

c this routine must be called at the beginning of a taylor plot

        implicit none

        call taylor_init

        call qstart

        call qfont('Times-Roman')
        call qgray(0.)
        call qtxts(10)

        end

c*****************************************************************

        subroutine taylor_end

c this routine must be called at the end of a taylor plot

        implicit none

        call qend

        end

c*****************************************************************
c*****************************************************************
c*****************************************************************

        subroutine taylor_init

c initializes taylor routines (internal routine, do not call explicitly)

        implicit none

        call taylor_set_param(1.3,1.0,0.2,1)  !initialize normalized
        call taylor_set_config(.false.,.false.,.false.)	!no lines
        call taylor_set_text('Standard Deviation','Correlation')

        end

c*****************************************************************

        subroutine taylor_set_text(std_text,cor_text)

c sets parameters for plot

        implicit none

        character*(*) std_text         !text to be written for std
        character*(*) cor_text         !text to be written for cor

        character*60 std,cor
        common /taylor_t/std,cor
        save /taylor_t/

        std = std_text
        cor = cor_text

        end

c*****************************************************************

        subroutine taylor_get_text(std_text,cor_text)

c sets parameters for plot

        implicit none

        character*(*) std_text         !text to be written for std
        character*(*) cor_text         !text to be written for cor

        character*60 std,cor
        common /taylor_t/std,cor
        save /taylor_t/

        std_text = std
        cor_text = cor

        end

c*****************************************************************
c next only for compatibility ... do not use anymore
c*****************************************************************

        subroutine taylor_set(sigma_tot,sigma_r,d_sigma,ndec)
	call taylor_set_param(sigma_tot,sigma_r,d_sigma,ndec)
	end

        subroutine taylor_get(sigma_tot,sigma_r,d_sigma,ndec)
	call taylor_get_param(sigma_tot,sigma_r,d_sigma,ndec)
	end

c*****************************************************************
c finished compatibility section
c*****************************************************************

        subroutine taylor_set_param(sigma_tot,sigma_r,d_sigma,ndec)

c sets parameters for plot

        implicit none

        real sigma_tot      !maximum of sigma to plot
        real sigma_r        !sigma of reference point
        real d_sigma        !plot marker every d_s sigma
        integer ndec        !decimals for sigma value

        real sigma_tot_c,sigma_r_c,d_sigma_c
        integer ndec_c
        common /taylor/sigma_tot_c,sigma_r_c,d_sigma_c,ndec_c
        save /taylor/

        sigma_tot_c = sigma_tot
          sigma_r_c = sigma_r
          d_sigma_c = d_sigma
             ndec_c = ndec

        end

c*****************************************************************

        subroutine taylor_get_param(sigma_tot,sigma_r,d_sigma,ndec)

c returns parameters for plot

        implicit none

        real sigma_tot      !maximum of sigma to plot
        real sigma_r        !sigma of reference point
        real d_sigma        !plot marker every d_s sigma
        integer ndec        !decimals for sigma value

        real sigma_tot_c,sigma_r_c,d_sigma_c
        integer ndec_c
        common /taylor/sigma_tot_c,sigma_r_c,d_sigma_c,ndec_c
        save /taylor/

        sigma_tot = sigma_tot_c
          sigma_r = sigma_r_c
          d_sigma = d_sigma_c
             ndec = ndec_c

        end

c*****************************************************************

        subroutine taylor_set_config(bstd,bcor,brms)

c sets configuration parameters for plot

        implicit none

	logical bstd		!plot lines for standard deviation
	logical bcor		!plot lines for correlation coefficient
	logical brms		!plot lines for RMS error

        logical bstd_c,bcor_c,brms_c
        common /taylor_log/bstd_c,bcor_c,brms_c
        save /taylor_log/

	bstd_c = bstd
	bcor_c = bcor
	brms_c = brms

        end

c*****************************************************************

        subroutine taylor_get_config(bstd,bcor,brms)

c gets configuration parameters for plot

        implicit none

	logical bstd		!plot lines for standard deviation
	logical bcor		!plot lines for correlation coefficient
	logical brms		!plot lines for RMS error

        logical bstd_c,bcor_c,brms_c
        common /taylor_log/bstd_c,bcor_c,brms_c
        save /taylor_log/

	bstd = bstd_c
	bcor = bcor_c
	brms = brms_c

        end

c*****************************************************************
c*****************************************************************
c*****************************************************************

        subroutine taylor_data(sigma_f,r_corr,text)

c writes text (centered) at data point

        implicit none

        real sigma_f                !sigma of time series
        real r_corr                 !correlation coefficient
        character*(*) text          !text to write

        real beta,x,y
        real width,height

        beta = acos(r_corr)
        x = sigma_f * cos(beta)
        y = sigma_f * sin(beta)

        call qtsize(text,width,height)
        x = x - 0.5*width
        y = y - 0.5*height

        call qtext(x,y,text)
        
        end

c*****************************************************************

        subroutine taylor_point(sigma_f,r_corr,text)

c writes text with point at data point

        implicit none

        real sigma_f                !sigma of time series
        real r_corr                 !correlation coefficient
        character*(*) text          !text to write

        real beta,x,y
        character*50 line           !aux for text to write

        beta = acos(r_corr)
        x = sigma_f * cos(beta)
        y = sigma_f * sin(beta)

        line = '.'//text

        call qtext(x,y,line)
        
        end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine taylor_full

	implicit none

        real rad,dr,sigmar
        integer ndec

        real x,y
        real width,height
        real xmin,phimax,dc
        character*60 text
        character*60 std_text,cor_text
        logical brms,bcor,bquart,bstd

        integer ialfa

c-------------------------------------------------------
c set up
c-------------------------------------------------------

        bquart = .false.

        call taylor_get_param(rad,sigmar,dr,ndec)
        call taylor_get_config(bstd,bcor,brms)
        call taylor_get_text(std_text,cor_text)

        xmin = -rad
        phimax = 180.
        dc = 0.2

        if( bquart ) then
          xmin = 0.
          phimax = 90.
          dc = 0.1
        end if

c-------------------------------------------------------
c set up domain
c-------------------------------------------------------

        call set_domain(rad,bquart)

c-------------------------------------------------------
c plot domain
c-------------------------------------------------------

	call qline(xmin,0.,rad,0.)
	call qline(0.,0.,0.,rad)
	call qarc(0.,0.,rad,0.,phimax);

c-------------------------------------------------------
c standard deviation
c-------------------------------------------------------

        call std(rad,sigmar,dr,ndec,bquart,bstd)

c-------------------------------------------------------
c correlation coefficient
c-------------------------------------------------------

        call corr_coff(rad,dc,1.,bcor)
        if( .not. bquart ) call corr_coff(rad,dc,-1.,bcor)

c-------------------------------------------------------
c rms pattern
c-------------------------------------------------------

        if( brms ) then
          call rms_pattern(rad,sigmar,dr,ndec,bquart)
        end if

c-------------------------------------------------------
c sigma r (reference point)
c-------------------------------------------------------

        call setdash(0,0.)
	call qarc(sigmar,0.,rad/100.,0.,360.);

c-------------------------------------------------------
c title
c-------------------------------------------------------

        text = cor_text
        call qtsize(text,width,height)
        x = -0.5 * width
        y = 1.15 * rad
        call qtext(x,y,text)

        text = std_text
        call qtsize(text,width,height)
        x = -0.5 * width
        y = -4.0 * height
        call qtext(x,y,text)

c-------------------------------------------------------
c end of routine
c-------------------------------------------------------

	end

c*****************************************************************

	subroutine taylor_quart

	implicit none

        real rad,dr,sigmar
        integer ndec

        real x,y
        real width,height
        real xmin,phimax,dc
        real aux
        character*60 text
        character*60 std_text,cor_text
        logical brms,bcor,bquart,bstd

        integer ialfa

c-------------------------------------------------------
c set up
c-------------------------------------------------------

        bquart = .true.

        call taylor_get_param(rad,sigmar,dr,ndec)
        call taylor_get_config(bstd,bcor,brms)
        call taylor_get_text(std_text,cor_text)

        xmin = -rad
        phimax = 180.
        dc = 0.2

        if( bquart ) then
          xmin = 0.
          phimax = 90.
          dc = 0.1
        end if

c-------------------------------------------------------
c set up domain
c-------------------------------------------------------

        call set_domain(rad,bquart)

c-------------------------------------------------------
c plot domain
c-------------------------------------------------------

	call qline(xmin,0.,rad,0.)
	call qline(0.,0.,0.,rad)
	call qarc(0.,0.,rad,0.,phimax);

c-------------------------------------------------------
c standard deviation
c-------------------------------------------------------

        call std(rad,sigmar,dr,ndec,bquart,bstd)

c-------------------------------------------------------
c correlation coefficient
c-------------------------------------------------------

        call corr_coff(rad,dc,1.,bcor)
        if( .not. bquart ) call corr_coff(rad,dc,-1.,bcor)

c-------------------------------------------------------
c rms pattern
c-------------------------------------------------------

        if( brms ) then
          call rms_pattern(rad,sigmar,dr,ndec,bquart)
        end if

c-------------------------------------------------------
c sigma r (reference point)
c-------------------------------------------------------

        call setdash(0,0.)
	call qarc(sigmar,0.,rad/100.,0.,360.);

c-------------------------------------------------------
c title
c-------------------------------------------------------

        text = cor_text
        call qtsize(text,width,height)
        aux = 1. / sqrt(8.)
        x = 1.15 * rad / sqrt(2.)
        y = 1.15 * rad / sqrt(2.)
        x = x - aux * width
        y = y + aux * height
        call qtxtr(-45.)
        call qtext(x,y,text)

        text = std_text
        call qtsize(text,width,height)

        x = -6.0 * height
        y = 0.5 * (rad-width)
        call qtxtr(90.)
        call qtext(x,y,text)

        x = 0.5 * (rad-width)
        y = -6.0 * height
        call qtxtr(0.)
        call qtext(x,y,text)

c-------------------------------------------------------
c end of routine
c-------------------------------------------------------

	end

c*****************************************************************
c*****************************************************************
c*****************************************************************
c
c internal routines
c
c*****************************************************************
c*****************************************************************
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

        subroutine corr_coff(rad,dcor,corfac,bcor)

c plots correlation coefficient (lines and ticks)

        implicit none

        real rad                !radius of arc
        real dcor               !delta of corr to write
        real corfac             !+/- 1 (to plot pos or neg circle)
        logical bcor            !plot lines ?

        real pi,aux,dbeta,fact,tick
        real cor,beta,x,y,x1,y1
        real width,height
        real r,rcol
        integer j
        logical bpos
        character*40 text

        integer ncor,nticks
        !parameter(ncor=8,nticks=11)
        !real corcof(ncor)
        !real cortik(nticks)
        !data corcof /0.2,0.4,0.6,0.8,0.9,0.95,0.99,1.0/
        !data cortik /0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.99/
        real corcof(50)
        real cortik(50)

c-------------------------------------------------------
c initialize constants
c-------------------------------------------------------

	rcol = 0.3
	rcol = -1.

        pi = 4. * atan(1.)
        aux = 180./pi
        dbeta = 1./aux              !for starting point
        fact = 0.002 * rad
        tick = 0.01 * rad

        bpos = corfac .gt. 0.

        call setdash(4,fact)
	if( rcol .ge. 0. ) call qhue(rcol)

        call make_cor_intv(dcor,ncor,nticks,corcof,cortik)
        write(6,*) 'corcof: ',ncor,(corcof(j),j=1,ncor)
        write(6,*) 'cortik: ',nticks,(cortik(j),j=1,nticks)

c-------------------------------------------------------
c plot dashed lines and label
c-------------------------------------------------------

        do j=1,ncor
          cor = corcof(j)
          beta = acos(cor)
          x = rad*cos(beta)
          y = rad*sin(beta)
          x = corfac * x
          if( bcor ) call qline(0.,0.,x,y)

          cor = corfac * cor
          call make_text(cor,text)
          call qtsize(text,width,height)
          call qtxtr(corfac*aux*beta)

          write(6,*) 'cor: ',cor,beta,beta*aux,x,y

          r = 1.03 * rad
          if( .not. bpos ) r = r + width

          x = r*cos(beta-dbeta)
          y = r*sin(beta-dbeta)
          x = corfac * x
          call qtext(x,y,text)

        end do

c-------------------------------------------------------
c label zero
c-------------------------------------------------------

        if( bpos ) then
          call make_text(0.0,text)
          call qtxtr(90.)
          beta = 0.5 * pi
          x = 1.03*rad*cos(beta-dbeta)
          y = 1.03*rad*sin(beta-dbeta)
          call qtext(x,y,text)
        end if

        call qtxtr(0.)

c-------------------------------------------------------
c plot ticks on axis
c-------------------------------------------------------

        call setdash(0,0.)
        do j=1,nticks
          cor = cortik(j)
          beta = acos(cor)
          x = rad*cos(beta)
          y = rad*sin(beta)
          x1 = (rad-tick)*cos(beta)
          y1 = (rad-tick)*sin(beta)
          x = corfac * x
          x1 = corfac * x1
          call qline(x,y,x1,y1)
          !call qline(-x,y,-x1,y1)
        end do

	call qgray(0.)

c-------------------------------------------------------
c end of routine
c-------------------------------------------------------

        end

c*****************************************************************

        subroutine make_cor_intv(dcor,ncor,nticks,corcof,cortik)

c makes intervals for corr coef

        implicit none

        real dcor               !line density (only 0.1 and 0.2)
        integer ncor            !total number of entries in corcof
        integer nticks          !total number of entries in cortik
        real corcof(1)          !labels for corr coef
        real cortik(1)          !ticks for corr coef

        integer nd,i,j

        nd = nint(10.*dcor)

        if( nd .ne. 1 .and .nd .ne. 2 ) then
          write(6,*) 'wrong value for dcor: ',dcor
          write(6,*) 'only 0.1 or 0.2 allowed'
          stop 'error stop make_cor_intv'
        end if

        j = 0
        do i=nd,8,nd
          j = j + 1
          corcof(j) = 0.1 * i
        end do

        corcof(j+1) = 0.9
        corcof(j+2) = 0.95
        corcof(j+3) = 0.99
        corcof(j+4) = 1.0
        ncor = j+4

        j = 0
        do i=nd,18,nd
          j = j + 1
          cortik(j) = 0.05 * i
        end do

        if( nd .eq. 2 ) then
          cortik(j+1) = 0.95
          cortik(j+2) = 0.99
          nticks = j+2
        else
          do i=1,9
            j = j + 1
            cortik(j) = 0.9 + 0.01 * i
          end do
          nticks = j
        end if

        end

c*****************************************************************

        subroutine rms_pattern(rad,sigmar,dr,ndec,bquart)

c plots rms pattern

        implicit none

        real rad,sigmar,dr
        integer ndec
        logical bquart          !want only quart of disk ?

        integer n,j,ndig
        real pi,aux,phimin
        real aleg,fact,r,phimax
        real cosphi,phi
        real phileg,x,y
	real rcol
        character*40 text

        integer ialfa

	rcol = 0.6
	rcol = -1.

        pi = 4. * atan(1.)
        aux = 180./pi
        phimin = 10.            !minimum arc to be plotted
        aleg = 3./4.
        fact = 0.002 * rad

        call setdash(5,fact)
	if( rcol .ge. 0. ) call qhue(rcol)
        n = 0.1 + (rad+sigmar)/dr

        do j=1,n
          r = j * dr
          cosphi = (rad**2 - r**2 - sigmar**2) / (2.*r*sigmar)
          if( cosphi .ge. 1. ) then
            phi = 0.
          else
            phi = aux * acos(cosphi)
          end if
          phimax = 180.
          if( bquart .and. sigmar .lt. r ) then
            phimax = aux * acos(-sigmar/r)
          end if

          phileg = aleg*phimax + (1.-aleg)*phi
          x = (0.05*dr+r)*cos(phileg/aux) + sigmar
          y = (0.05*dr+r)*sin(phileg/aux)
          ndig = ialfa(r,text,ndec,-1)

          if( phimax-phi .gt. phimin ) then
	    call qarc(sigmar,0.,r,phi,phimax);
            call qtxtr(phileg-90.)
            call qtext(x,y,text)
          end if
          write(6,*) 'rms: ',j,cosphi,phi,phileg,phimax
        end do

        call qtxtr(0.)
        call qgray(0.)

        end

c*****************************************************************

        subroutine std(rad,sigmar,dr,ndec,bquart,bstd)

c plots std

        implicit none

        real rad,sigmar,dr
        integer ndec
        logical bquart          !want only quart of disk ?
        logical bstd            !want all lines (.true.) or only at sigmar ?

        integer n,j,ndig
        real pi,aux,phimax
        real fact,tick
        real x,y,r,rcol
        real width,height
        character*40 text

        integer ialfa

	rcol = 0.9
	rcol = -1.

        pi = 4. * atan(1.)
        aux = 180./pi
        phimax = 180.
        if( bquart ) phimax = 90.
        fact = 0.002 * rad
        tick = 0.01 * rad

c-------------------------------------------------------
c make circles
c-------------------------------------------------------

        call setdash(4,fact)
	if( rcol .ge. 0. ) call qhue(rcol)
        n = 0.1 + rad/dr

        if( bstd ) then
          do j=1,n
            r = j * dr
	    call qarc(0.,0.,r,0.,phimax);
          end do
        else
	  call qarc(0.,0.,sigmar,0.,phimax);
        end if

c-------------------------------------------------------
c label zero
c-------------------------------------------------------

        call setdash(0,0.)

        r = 0.
        ndig = ialfa(r,text,ndec,-1)
        call qtsize(text,width,height)

        x = r - 0.5 * width
        y = 0. - 1.8 * height
        call qtext(x,y,text)

        if( bquart ) then
          x = 0. - 1.8 * width
          y = r - 0.5 * height
          call qtext(x,y,text)
        end if

c-------------------------------------------------------
c make ticks and label
c-------------------------------------------------------

        do j=1,n
          r = j * dr

          ndig = ialfa(r,text,ndec,-1)
          call qtsize(text,width,height)

          x = r - 0.5 * width
          y = 0. - 1.8 * height
          call qtext(x,y,text)
          call qline(r,0.,r,tick)

          if( .not. bquart ) then
            x = - r - 0.5 * width
            call qtext(x,y,text)
            call qline(-r,0.,-r,tick)
          else
            x = 0. - 2.0 * width
            y = r - 0.5 * height
            call qtext(x,y,text)
            call qline(0.,r,tick,r)
          end if

          write(6,*) 'std: ',j,r,ndec,ndig,x,y,width
     +                  ,' |'//text(1:ndig)//'|'
        end do

	call qgray(0.)

c-------------------------------------------------------
c end of routine
c-------------------------------------------------------

        end

c*****************************************************************

        subroutine set_domain(rad,bquart)

c sets up domain

        implicit none

        real rad
        logical bquart

        real v0,dvx,dvy
        real rmax,rmin
        real xmin,xmax,ymin,ymax
        real alpha

        alpha = 1.3
        rmax = alpha * rad
        rmin = (alpha - 1.) * rad

        xmin = -rmax
        if( bquart ) xmin = -rmin
        xmax = rmax
        ymin = -rmin
        ymax = rmax

        v0 = 1.
        dvx = 18.
        dvy = dvx * (ymax-ymin) / (xmax-xmin)
	call qsetvp(v0,v0,v0+dvx,v0+dvy)
	call qworld(xmin,ymin,xmax,ymax)

        end

c*****************************************************************
c*****************************************************************
c*****************************************************************
c
c taylor example
c
c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine taylor_example

c shows usage of taylor routines

	implicit none

c-----------------------
	call qopen
c-----------------------

        call taylor_start
	call taylor_set_param(9.,5.5,2.,-1)
        call taylor_full
        call taylor_data(8.,0.9,'A')
        call taylor_data(6.,0.6,'B')
        call taylor_point(8.,0.4,'C')
        call taylor_end

c-----------------------

        call taylor_start
	call taylor_set_param(9.,5.5,2.,-1)
        call taylor_quart
        call taylor_data(8.,0.9,'A')
        call taylor_data(6.,0.6,'B')
        call taylor_point(8.,0.4,'C')
        call taylor_end

c-----------------------

        call taylor_start
	call taylor_set_param(1.3,1.0,0.2,1)  !normalized
	call taylor_set_config(.true.,.true.,.false.)
	call taylor_set_text('Standard Deviation (Normalized)'
     +                  ,'Correlation Coefficient')
        call taylor_quart
        call taylor_data(1.2,0.9,'A')
        call taylor_data(0.9,0.6,'B')
        call taylor_point(1.0,0.4,'C')
        call taylor_end

c-----------------------
	call qclose
c-----------------------

	end

c*****************************************************************
c
c to test the taylor routines uncomment the following three lines
c
c*****************************************************************

c        program taylor_test
c        call taylor_example
c        end

c*****************************************************************

