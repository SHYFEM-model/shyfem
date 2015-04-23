c
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
c
c**********************************************************


c**********************************************************
c**********************************************************
c**********************************************************

        subroutine plo_part(n,xlag,ylag,zlag,llag,lmax,title)

	implicit none

	integer n
	real xlag(1)
	real ylag(1)
	real zlag(1)
	integer llag(1)
	integer lmax
        character*(*) title

	real pmin,pmax
	real getpar

	call qstart
        call annotes(title)
        call bash(0)

        call qcomm('Plotting particles')
        call plo_xy(n,xlag,ylag,zlag,llag,lmax)

        call bash(2)
	call qend

	end

c**********************************************************
c**********************************************************
c**********************************************************

	subroutine lag_get_header_new(iunit,nvers,lmax)

	implicit none

	integer iunit,nvers,lmax

	integer mtype

	lmax = 1

	read(iunit,end=98,err=99) mtype,nvers

	if( mtype .ne. 367265 ) goto 97
	if( nvers .le. 2 ) goto 95
	if( nvers .gt. 5 ) goto 96

	if( nvers .ge. 5 ) read(iunit) lmax

	return
   95	continue
	write(6,*) mtype,nvers
	stop 'error stop lag_get_header: cannot read this version'
   96	continue
	write(6,*) mtype,nvers
	stop 'error stop lag_get_header: unknown nvers'
   97	continue
	write(6,*) mtype,nvers
	stop 'error stop lag_get_header: unknown mtype'
   98	continue
	stop 'error stop lag_get_header: file is empty'
   99	continue
	stop 'error stop lag_get_header: error reading header'
	end

c**********************************************************

	subroutine lag_get_xy_new(iunit,nvers,ndim,it,n
     +			,xlag,ylag,zlag,llag)

	implicit none

	integer iunit
	integer nvers
	integer ndim
	integer it
	integer n
	real xlag(1),ylag(1),zlag(1)
	integer llag(1)

	integer nbdy,nn,nout
	integer id,ie,ies,i,lb
	real x,y,z,xs,ys,zs,ts

c-----------------------------------------------------
c read in particles
c-----------------------------------------------------

	read(iunit,end=100,err=99) it,nbdy,nn,nout

	!write(6,*) 'reading lgr: ',iunit,ndim

	n = 0

        do i=1,nn
	  if( nvers .eq. 3 ) then
            read(iunit,err=98) id,x,y,z,ie,xs,ys,zs,ies
	  else if( nvers .eq. 4 ) then
            read(iunit,err=98) id,x,y,z,ie,xs,ys,zs,ies,ts
	  else if( nvers .eq. 5 ) then
            read(iunit,err=98) id,x,y,z,ie,lb,xs,ys,zs,ies,ts
	  else
	    write(6,*) 'nvers = ',nvers
	    stop 'error stop lag_get_xy_new: internal error (1)'
	  end if
          if( ie .gt. 0 ) then
	    n = n + 1
	    if( n .gt. ndim ) goto 97
	    xlag(n) = x
	    ylag(n) = y
	    zlag(n) = z
	    llag(n) = lb
	  end if
        end do

	write(6,*) 'particles read ',n

	return
  100	continue
	n = -1
	return
   97	continue
	stop 'error stop lag_get_xy_new: ndim too small'
   98	continue
	stop 'error stop lag_get_xy_new: error reading data record'
   99	continue
	stop 'error stop lag_get_xy_new: error reading time record'
	end

c**********************************************************
c**********************************************************
c**********************************************************

	subroutine lag_get_header(iunit,nvers)

	implicit none

	integer iunit,nvers

	nvers = 1

	end

c**********************************************************

	subroutine lag_get_xy(iunit,ndim,it,n,xlag,ylag,zlag)

c old subroutine - do not use anymore

	implicit none

	integer iunit
	integer ndim
	integer it
	integer n
	real xlag(1),ylag(1),zlag(1)

	integer nout
	real rnout
	real x,y
	real wint,wx,wy

c-----------------------------------------------------
c read in particles
c-----------------------------------------------------

	!write(6,*) 'reading lgr: ',iunit,ndim

	n = 0

    1	continue

          read(iunit,*,end=2,err=99) x,y

          if(x.eq.0 .and. y.eq.0) goto 3

	  n = n + 1
	  if( n .gt. ndim ) goto 97
	  xlag(n) = x
	  ylag(n) = y
	  zlag(n) = 0.5

	goto 1

c-----------------------------------------------------
c read in header
c-----------------------------------------------------

    3	continue

	!write(6,*) 'finished reading particles ',n

	read(iunit,*) rnout
	read(iunit,*) it
	read(iunit,*) wint
	read(iunit,*) wx,wy

	return

c-----------------------------------------------------
c end of read
c-----------------------------------------------------

    2	continue

	n = -1

	return
   97	continue
	stop 'error stop lag_get_xy: ndim too small'
   99	continue
	write(6,*) n,x,y
	stop 'error stop lag_get_xy: read error'
	end

c**********************************************************

	subroutine plolagr

	implicit none

	include 'param.h'

	include 'basin.h'


	integer ndim
	parameter(ndim=nbdydim)

	integer iunit,it,n,nvers,lmax
	real xlag(ndim)
	real ylag(ndim)
	real zlag(ndim)
	integer llag(ndim)

	logical ptime_ok,ptime_end
	integer nrec
        integer ifemop

        iunit = ifemop('.lgr','unform','unknown')
        if( iunit .le. 0 ) stop
	!call lag_get_header(iunit,nvers)	!old
	call lag_get_header_new(iunit,nvers,lmax)
	write(6,*) 'lagrangian unit opened: ',iunit,nvers

	nrec = 0
                        
    1   continue

	!call lag_get_xy(iunit,ndim,it,n,xlag,ylag,zlag)	!old
	call lag_get_xy_new(iunit,nvers,ndim,it,n,xlag,ylag,zlag,llag)
	if( n .lt. 0 ) goto 2
        nrec = nrec + 1
        write(6,*) nrec,it,n

	call ptime_set_itime(it)

	if( ptime_end() ) goto 2

	if( ptime_ok() ) then
	  write(6,*) 'plotting particles ',it,n
          call plo_part(n,xlag,ylag,zlag,llag,lmax,'particles')
	end if

	goto 1
    2   continue

	end

c**********************************************************

        subroutine plo_xy(n,xlag,ylag,zlag,llag,lmax)
     
        implicit none
       
	integer n
	real xlag(1)
	real ylag(1)
	real zlag(1)
	integer llag(1)
	integer lmax

	logical bcolor,bbottom,bzeta,bstation,bplot,blayer
	integer i,is,ip_station,l
        real x,y,z,x1,y1,x2,y2
	real cwater,cbottom,col

c--------------------------------------
	bcolor = .false.	!use color to plot lagrangian particles
	bcolor = .true.		!use color to plot lagrangian particles
	bbottom = .false.	!use color to indicate p. on bottom
	blayer = .true.		!use layer to decide on color to use
	bzeta = .true.		!use zeta to decide on color to use
	bzeta = .false.		!use zeta to decide on color to use
	bstation = .false.	!use zeta to get station number (and color)
	ip_station = 0		!if different from 0 -> plot only this station

	cwater = 1.0		!color to use for particles in water
	cbottom = 0.3		!color to use for particles on bottom
c--------------------------------------

        call qlwidth(0.035)
	if( bcolor ) call qhue(cwater)		!LAG_COLOR
	is = 0
	bplot = .true.

	do i=1,n

	   x = xlag(i)
	   y = ylag(i)
	   z = zlag(i)
	   l = llag(i)

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
	       else				! on bottom
	  	call qhue(cbottom)
	       end if
	     end if
	   end if

	   if( bplot ) call plot_single_particle(x,y)

	end do

	call qgray(0.)

       end

c**********************************************************

	subroutine plot_single_particle(x,y)

	implicit none

	real x,y

	real dr
	real x1,y1,x2,y2

	dr = 0.1

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

