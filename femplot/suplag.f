
!--------------------------------------------------------------------------
!
!    Copyright (C) 2005,2007-2009,2012,2015  Georg Umgiesser
!    Copyright (C) 2018  Christian Ferrarin
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

c utility routines for plotsim
c
c contents :
c
c subroutine plolagr
c
c revision log :
c
c 03.03.2005    ggu     revised lagrangian plotting routines
c 05.06.2007    ggu     plot lagrangian particles in color (LAG_COLOR)
c 13.06.2008    ggu     use also z for plotting
c 17.09.2008    ggu     new version for plotting lagrangian particles in color
c 27.01.2009    ggu     better error check reading particles, routines deleted
c 23.01.2012    ggu     use nbdydim in plolagr, plot z with color in plo_xy
c 01.10.2012    ggu     use station to get color in plot (ip_station)
c 23.04.2015    ggu     plotting for 3d particles
c 19.09.2018    ccf     plotting lagrangian trajectories
c
c**********************************************************


c**********************************************************
c**********************************************************
c**********************************************************

        subroutine plo_traj(n,nt,nrec,lgmean,xall,yall,rall,aplot,
     +                       nm,xmll,ymll,rmll,mplot,title)

        use mod_hydro_plot
	use plotutil
        use basin

	implicit none

	integer n,nt,nrec,lgmean
	real xall(n,0:nt)
	real yall(n,0:nt)
	real rall(n,0:nt)
	integer aplot(n)
	integer nm
	real xmll(nm,0:nt)
	real ymll(nm,0:nt)
	real rmll(nm,0:nt)
	integer mplot(nm)

        character*(*) title

        integer icsave
	real lw
	real pmin,pmax,flag
	real getpar

	integer ir

	call qstart
        call annotes(title)
        call bash(0)

	if( n > 0 ) then
	  call get_flag(flag)
          call get_minmax_flag(rall(:,50),n,pmin,pmax,flag)
          if (bverb) write(6,*) 'min/max: ',n,pmin,pmax 
          call colauto(pmin,pmax)
          call qcomm('Plotting lgr trajectories')

          call get_color_table(icsave)
	  ! all particles
          if ( lgmean < 2 ) then
            lw = 0.001
            call qlwidth(lw)		!line thickness
            if ( lgmean == 1 ) aplot = 2
  	    do ir = 1,n
              call plo_traj_line(nrec,xall(ir,:),yall(ir,:),
     +              rall(ir,:),aplot(ir),flag)
            end do

!	    do ir = 1,nrec
!              call plo_traj_line(n,xall(:,ir),yall(:,ir),
!     +	  	    xall(:,ir-1),yall(:,ir-1),
!     +              rall(:,ir),aplot,flag)
!            end do
          end if

	  ! mean particles
          if ( lgmean >= 2 ) then
            lw = 0.001
          else
            lw = 0.1
          end if
          call qlwidth(lw)		!line thickness
	  do ir = 1,nm
            call plo_traj_line(nrec,xmll(ir,:),ymll(ir,:),
     +              rmll(ir,:),mplot(ir),flag)
          end do

	  ! plot 1 mean trajectory on top with tichk line
          if ( lgmean == 3 ) then
            call qcomm('Plotting first lgr trajectory')
            lw = 0.1
            call qlwidth(lw)		!line thickness
            call plo_traj_line(nrec,xmll(1,:),ymll(1,:),
     +              rmll(1,:),mplot(1),flag)
          end if

          call set_color_table(icsave)
          call qlwidth(0.010)
	  call colsh
	  call qgray(0.)
	end if

        call bash(2)
	call qend

	end

c**********************************************************

        subroutine plo_part(n,xlag,ylag,rlag,iplot,title)

	use plotutil

        implicit none

        integer n
        real xlag(n)
        real ylag(n)
        real rlag(n)            !variable to plot
        integer iplot(n)
        character*(*) title

        real pmin,pmax,flag
        real getpar

        call qstart
        call annotes(title)
        call bash(0)

        if( n > 0 ) then
          call get_flag(flag)
          call get_minmax_flag(rlag,n,pmin,pmax,flag)
          if (bverb) write(6,*) 'min/max: ',n,pmin,pmax 
          call colauto(pmin,pmax)
          call qcomm('Plotting particles')
          call plo_xy_new(n,xlag,ylag,rlag,iplot)
          call colsh
        end if

        call bash(2)
        call qend

        end

c**********************************************************
c**********************************************************
c**********************************************************

        subroutine write_xy(n,xlag,ylag,zlag,llag,lmax)

        implicit none

        integer n
        real xlag(n)
        real ylag(n)
        real zlag(n)
        integer llag(n)
        integer lmax

        integer i,l
        real x,y,z

        do i=1,n
           x = xlag(i)
           y = ylag(i)
           z = zlag(i)
           l = llag(i)
           write(6,*) i,l,x,y,z
        end do

        end

c**********************************************************

        subroutine plo_xy_new(n,xlag,ylag,rlag,iplot)

        use color

        implicit none

        integer n
        real xlag(n)
        real ylag(n)
        real rlag(n)
        integer iplot(n)

        integer i
        real x,y,r
        logical bplot
        integer icsave
        real col
        real get_color

        call get_color_table(icsave)
        !call set_color_table(-1)
        call qlwidth(0.065)             !particle size

!todo assign color defined in the apn and in function of variable to plot zlag

        do i=1,n
          bplot = iplot(i) .eq. 1
          if( .not. bplot ) cycle

          x = xlag(i)
          y = ylag(i)
          r = rlag(i)

          col = get_color(r,isoanz+1,ciso,fiso)

          call qsetc(col)
          call plot_single_particle(x,y)
        end do
        call set_color_table(icsave)

        call qlwidth(0.010)
        call qgray(0.)

        end

c**********************************************************

        subroutine plo_traj_line(n,xlag,ylag,rlag,iplot,flag)
     
	use color

        implicit none
      
	integer n
	real xlag(n)
	real ylag(n)
	real xold(n)
	real yold(n)
	real rlag(n)
	integer iplot
        real flag

	integer i
        real x,y,r,xo,yo
	real xmean,ymean,rmean,xmeao,ymeao
	logical bplot
        real col
        real get_color

!todo assign color defined in the apn and in function of variable to plot zlag

	bplot = iplot > 0
	do i=2,n
	  if( .not. bplot ) cycle

	  x = xlag(i)
	  y = ylag(i)
	  r = rlag(i)
	  xo = xlag(i-1)
	  yo = ylag(i-1)
          if ( x == flag .or. xo == flag ) cycle

          if ( iplot == 1 ) then
            col = get_color(r,isoanz+1,ciso,fiso)
            call qsetc(col)
          else
            call qgray(0.5)
	  end if
          
          call qline(xo,yo,x,y)
	end do

       end

c**********************************************************
! CCF TO BE DELETED

        subroutine plo_xy(n,xlag,ylag,zlag,llag,hlag,tlag,lmax)
     
        implicit none
       
	integer n
	real xlag(n)
	real ylag(n)
	real zlag(n)
	integer llag(n)
	real hlag(n)
	integer tlag(n)
	integer lmax

	logical bcolor,bbottom,bzeta,bstation,bplot,blayer
	logical bdepth, btype
	integer i,is,ip_station,l
        real x,y,z,x1,y1,x2,y2
	real h,r,c
	real cwater,cbottom,col
	integer t,tm

c--------------------------------------
	bcolor = .false.	!use color to plot lagrangian particles
	bcolor = .true.		!use color to plot lagrangian particles
	bbottom = .true.	!use color to indicate p. on bottom
	bbottom = .false.	!use color to indicate p. on bottom
	blayer = .true.		!use layer to decide on color to use
	blayer = .false.	!use layer to decide on color to use
	bzeta = .true.		!use zeta to decide on color to use
	bzeta = .false.		!use zeta to decide on color to use
	bdepth = .false.	!use relative depth to decide on color to use
	btype  = .true.		!use type to decide on color to use
	bstation = .false.	!use zeta to get station number (and color)
	ip_station = 0		!if different from 0 -> plot only this station

	cwater = 1.0		!color to use for particles in water
	cbottom = 0.3		!color to use for particles on bottom
	cwater = 0.45		!color to use for particles in water
	cbottom = 0.9		!color to use for particles on bottom
c--------------------------------------

        call qlwidth(0.065)
	if( bcolor ) call qhue(cwater)		!LAG_COLOR
	is = 0
	bplot = .true.

	tm = -999
	do i=1,n
	  t = tlag(i)
	  tm = max(tm,t)
	end do

	do i=1,n

	   x = xlag(i)
	   y = ylag(i)
	   z = zlag(i)
	   l = llag(i)
	   h = hlag(i)
	   t = tlag(i)

	   if( bcolor ) then
	     if( bstation ) then
               call get_station_color(z,is,col)
               call qhue(col)
	       bplot = ip_station .eq. 0 .or. ip_station .eq. is
	     else if( blayer .and. lmax > 1 ) then
	       col = float(l)/lmax
	       !write(6,*) i,l,x,y,z,col
	       call qhue(col)
	     else if( bzeta ) then
	       call qhue(z)
	     else if( bbottom ) then
	       if( z .lt. 1.0 ) then		! in water column
	  	call qhue(cwater)
	       !else
	  !	call qhue(cbottom)
	       elseif ( z.eq.1. .and. l.eq.lmax) then	! on bottom
	  	call qhue(cbottom)
	       else 
	  	call qhue(cwater)
	       end if
	     else if( bdepth ) then
	       if ( h .gt. 0.95 .or. l.eq.0 ) then
	         call qgray(0.)
	       else
	         call qhue(h)
	       end if
	     else if( btype ) then
	         r = float(t)/float(tm)
	         call qhue(r)
		 !c = r*h
	         !if ( h .gt. 0.95 .or. l.eq.0 ) then
		 !  c = r
	         !  call qhue(c)
	 	 !else
	         !  call qhue(c)
	         !end if
	     end if
	   end if

	   if( bplot ) call plot_single_particle(x,y)

	end do
        call qlwidth(0.025)

	call qgray(0.)

       end

c**********************************************************

	subroutine plot_single_particle(x,y)

	implicit none

	real x,y

	real dr
	real x1,y1,x2,y2

	dr = 0.001

	!write(6,*) 'plot_single_particel: ',x,y

        x1=x+dr
        y1=y+dr
        x2=x-dr
        y2=y+dr        
        call qline(x1,y1,x2,y2)
        x1=x-dr
        y1=y-dr
        x2=x-dr
        y2=y+dr
        call qline(x1,y1,x2,y2)
        x1=x-dr
        y1=y-dr
        x2=x+dr
        y2=y-dr
        call qline(x1,y1,x2,y2)
        x1=x+dr
        y1=y+dr
        x2=x+dr
        y2=y-dr
        call qline(x1,y1,x2,y2)           

	end

c**********************************************************

        subroutine get_station_color(z,is,col)

        implicit none

        real z		!depth released (not changed in 2D)
        integer is	!station number (return)
        real col	!color to be used (return)

        integer n,next
        real colintern,colextern,deltac

        !n = 28
        n = 39
        next = 4
        colintern = 0.35
        colextern = 0.55
        deltac=(0.9-colextern)/next
        is = nint(n*z)  !station number

       !if( is.eq.27) then
       !if( is.eq.39.or.is.eq.37) then
       if( is.eq.1.or.is.eq.2) then
           col = colextern
       !else if( is .eq. 26 ) then
       !else if(is.eq.34.or.is.eq.31) then
       else if(is.eq.3.or.is.eq.4) then
           !col = colextern+0.05
           col = colextern+deltac
       !else if( is .eq. 22 ) then
       !else if(is.eq.27.or.is.eq.26 ) then
       else if(is.eq.5.or.is.eq.6 ) then
           !col = colextern+0.1
           col = colextern+deltac*3
       !else if( is .eq. 4 ) then
       !else if(is.eq.23 .or.is.eq.8 ) then
       else if(is.eq.7 .or.is.eq.8 ) then
           !col = colextern+0.15
           col = colextern+deltac*2
       !else if( is .eq. 1 ) then
       !else if( is.eq.2.or. is.eq.4 ) then
       else if( is.eq.9 ) then
           !col = colextern+0.20
           col = colextern+deltac*4


!if( (is.ge.29).and.(is.le.32)) then
!       col = colextern
        else
                !col = colintern
                col = 0. + (colintern)*z
        end if

!       if( col .eq. colextern ) then
                !col = colextern + (0.85-colextern)*(z)
!       else
!               col = 0. + (colintern)*z
!       end if

	end

c**********************************************************
