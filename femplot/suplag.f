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

        subroutine plo_part(n,xlag,ylag,rlag,iplot,title)

	implicit none

	integer n
	real xlag(1)
	real ylag(1)
	real rlag(1)		!variable to plot
	integer iplot(1)
        character*(*) title

	real pmin,pmax
	real getpar

	call qstart
        call annotes(title)
        call bash(0)

        call get_minmax_flag(rlag,n,pmin,pmax)
        write(6,*) 'min/max: ',n,pmin,pmax 
        call colauto(pmin,pmax)

!todohere: convert rlag to colors to plot

        call qcomm('Plotting particles')
        call plo_xy_new(n,xlag,ylag,rlag,iplot)

	call colsh

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

	subroutine lag_peek_xy_header(iunit,it,nb,nn,nout)

	read(iunit,end=100,err=99) it,nbdy,nn,nout
	backspace(iunit)

	return
  100	continue
	backspace(iunit)
	n = -1
	return
   99	continue
	stop 'error stop lag_get_xy_new: error reading time record'
	end

c**********************************************************

	subroutine lag_alloc(nn
     +		,xlag,ylag,zlag,llag,hlag,tlag,alag,clag,rlag,iplot)

	implicit none

	integer nn

	real, allocatable :: xlag(:),ylag(:)	
	real, allocatable :: zlag(:)		
	real, allocatable :: hlag(:)		
	integer, allocatable :: llag(:)		
	integer, allocatable :: tlag(:)		
	real, allocatable :: alag(:) 		
        real, allocatable :: clag(:)            
        real, allocatable :: rlag(:)            
        integer, allocatable :: iplot(:)            

	if( allocated(zlag) ) then
	  deallocate(xlag)
	  deallocate(ylag)
	  deallocate(zlag)
	  deallocate(hlag)
	  deallocate(llag)
	  deallocate(tlag)
	  deallocate(alag)
	  deallocate(clag)
	  deallocate(rlag)
	  deallocate(iplot)
	end if

	allocate(xlag(nn))
	allocate(ylag(nn))
	allocate(zlag(nn))
	allocate(hlag(nn))
	allocate(llag(nn))
	allocate(tlag(nn))
	allocate(alag(nn))
	allocate(clag(nn))
	allocate(rlag(nn))
	allocate(iplot(nn))

	end

c**********************************************************

	subroutine lag_get_xy_new(iunit,nvers,ndim,it,n
     +			,xlag,ylag,zlag,llag,hlag,tlag,alag,clag)

	implicit none

	integer iunit
	integer nvers
	integer ndim
	integer it
	integer n
	real xlag(1),ylag(1)	!particle position
	real zlag(1)		!relative depth of particle on level l [0-1]
	real hlag(1)		!absolute relative depth of particle [0-1]
	integer llag(1)		!verical laver of particle
	integer tlag(1)		!particle type
	real alag(1) 		!age of particle [s]
        real clag(1)            !custom property of particle

	integer nbdy,nn,nout
	integer id,ie,ies,i,lb,t
	real x,y,z,h,xs,ys,zs,ts,c

c-----------------------------------------------------
c read in particles
c-----------------------------------------------------

	read(iunit,end=100,err=99) it,nbdy,nn,nout

	n = 0

        do i=1,nn
	  if( nvers .eq. 3 ) then
            read(iunit,err=98) id,x,y,z,ie,xs,ys,zs,ies
	  else if( nvers .eq. 4 ) then
            read(iunit,err=98) id,x,y,z,ie,xs,ys,zs,ies,ts
	  else if( nvers .eq. 5 ) then
            read(iunit,err=98) id,x,y,z,ie,lb,h,t,c,xs,ys,zs,ies,ts
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
	    hlag(n) = h
	    tlag(n) = t
	    clag(n) = c
	    alag(n) = it - ts
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

	subroutine plolagr

	use mod_plot2d
	use basin

	implicit none

	include 'param.h'

	integer iunit,it,i,n,nvers,lmax,l
	real, allocatable :: xlag(:),ylag(:)	
	real, allocatable :: zlag(:)		
	real, allocatable :: hlag(:)		
	integer, allocatable :: llag(:)		
	integer, allocatable :: tlag(:)		
	real, allocatable :: alag(:) 		
        real, allocatable :: clag(:)            
        real, allocatable :: rlag(:)            
        integer, allocatable :: iplot(:)            

        character*80 name
	logical ptime_ok,ptime_end
	integer nrec
	integer nb,nout
        integer level
        integer ifemop
        integer getlev,getvar

	INTERFACE
	subroutine lag_alloc(nn
     +		,xlag,ylag,zlag,llag,hlag,tlag,alag,clag,rlag,iplot)
	integer nn
	real, allocatable :: xlag(:),ylag(:)	
	real, allocatable :: zlag(:)		
	real, allocatable :: hlag(:)		
	integer, allocatable :: llag(:)		
	integer, allocatable :: tlag(:)		
	real, allocatable :: alag(:) 		
        real, allocatable :: clag(:)            
        real, allocatable :: rlag(:)            
        integer, allocatable :: iplot(:)            
	end subroutine
	END INTERFACE

c----------------------------------------------------------------
c open lgr file and read header
c----------------------------------------------------------------

        iunit = ifemop('.lgr','unform','unknown')
        if( iunit .le. 0 ) stop
	call lag_get_header_new(iunit,nvers,lmax)
	write(6,*) 'lagrangian unit opened: ',iunit,nvers

	nrec = 0
                        
    1   continue

c----------------------------------------------------------------
c read lgr data
c----------------------------------------------------------------

	call lag_peek_xy_header(iunit,it,nb,n,nout)
	if( n .lt. 0 ) goto 2
	call lag_alloc(n
     +		,xlag,ylag,zlag,llag,hlag,tlag,alag,clag,rlag,iplot)

	call lag_get_xy_new(iunit,nvers,n,it,n
     +                     ,xlag,ylag,zlag,llag,hlag,tlag,alag,clag)

	if( n .lt. 0 ) goto 2
        nrec = nrec + 1
        write(6,*) nrec,it,n

c----------------------------------------------------------------
c set what to plot with color with option varnam in plots
c----------------------------------------------------------------

	call getfnm('varnam',name)

	if ( name .eq. '' ) then	    !nothing
	   rlag = 0.
	else if( name .eq. 'type' ) then    !type of particle
	   rlag = tlag
	else if( name .eq. 'depth' ) then   !absolute realtive depth
	   rlag = hlag
	else if( name .eq. 'age' ) then	    !age [s]
	   !rlag = alag
	   rlag = alag / 86400.		    !age [d]
	else if( name .eq. 'custom' ) then  !custom
	   rlag = clag
        else
            goto 99
	end if

c----------------------------------------------------------------
c set vertical level to plot with option lev3d in plots
c----------------------------------------------------------------
!todohere: not working with lev3d = -1

	iplot = 1
        level = getlev()

	if ( level .ne. 0 ) then
	  iplot = 0
	  do i = 1,n
	    l = llag(i)
	    if ( l .eq. level ) iplot(i) = 1
  	  end do
	end if

c----------------------------------------------------------------
c plot particles 
c----------------------------------------------------------------
	
	call ptime_set_itime(it)

	if( ptime_end() ) goto 2

	if( ptime_ok() ) then
	  write(6,*) 'plotting particles ',it,n
          call plo_part(n,xlag,ylag,rlag,iplot,'particles')
	end if

	goto 1
    2   continue

        return
   99   continue
        write(6,*) 'Unknown varnam ',name
        stop 'error stop plolagr: varnam error'

c----------------------------------------------------------------
c end of routine
c----------------------------------------------------------------

	end

c**********************************************************

        subroutine write_xy(n,xlag,ylag,zlag,llag,lmax)

        implicit none

        integer n
        real xlag(1)
        real ylag(1)
        real zlag(1)
        integer llag(1)
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
     
        implicit none
      
        include 'color.h'
 
	integer n
	real xlag(1)
	real ylag(1)
	real rlag(1)
	integer iplot(1)

	integer i
        real x,y,r
	logical bplot
        integer icsave
        real col
        real get_color

        call get_color_table(icsave)
        call set_color_table(-1)
        call qlwidth(0.065)		!particle size

!todo assign color defined in the apn and in function of variable to plot zlag

	do i=1,n
	  bplot = iplot(i) .eq. 1
	  if( .not. bplot ) cycle

	  x = xlag(i)
	  y = ylag(i)
	  r = rlag(i)

          col = get_color(r,isoanz+1,ciso,fiso)

          call qhue(col)
	  call plot_single_particle(x,y)
	end do
        call set_color_table(icsave)

        call qlwidth(0.010)
	call qgray(0.)

       end

c**********************************************************
! CCF TO BE DELETED

        subroutine plo_xy(n,xlag,ylag,zlag,llag,hlag,tlag,lmax)
     
        implicit none
       
	integer n
	real xlag(1)
	real ylag(1)
	real zlag(1)
	integer llag(1)
	real hlag(1)
	integer tlag(1)
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
