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
c
c**********************************************************


c**********************************************************
c**********************************************************
c**********************************************************

        subroutine plo_part(n,xlag,ylag,zlag,title)

	implicit none

	integer n
	real xlag(1)
	real ylag(1)
	real zlag(1)
        character*(*) title

	real pmin,pmax
	real getpar
        integer gettime

	call qstart
        call annotes(gettime(),title)
        call bash(0)

        call qcomm('Plotting particles')
        call plo_xy(n,xlag,ylag,zlag)

        call bash(2)
	call qend

	end

c**********************************************************
c**********************************************************
c**********************************************************

	subroutine lag_get_header_new(iunit,nvers)

	implicit none

	integer iunit,nvers

	integer mtype

	!read(iunit,err=99) mtype,nvers
	read(iunit,end=98,err=99) mtype,nvers

	if( mtype .ne. 367265 ) goto 97
	if( nvers .gt. 3 ) goto 96

	return
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

	subroutine lag_get_xy_new(iunit,ndim,it,n,xlag,ylag,zlag)

	implicit none

	integer iunit
	integer ndim
	integer it
	integer n
	real xlag(1),ylag(1),zlag(1)

	integer nbdy,nn,nout
	integer id,ie,ies,i
	real x,y,z,xs,ys,zs,ts

c-----------------------------------------------------
c read in particles
c-----------------------------------------------------

	read(iunit,end=100,err=99) it,nbdy,nn,nout

	!write(6,*) 'reading lgr: ',iunit,ndim

	n = 0

        do i=1,nn
          read(iunit,err=98) id,x,y,z,ie,xs,ys,zs,ies,ts
          if( ie .gt. 0 ) then
	    n = n + 1
	    if( n .gt. ndim ) goto 97
	    xlag(n) = x
	    ylag(n) = y
	    zlag(n) = z
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

	character*80 descrr
	common /descrr/ descrr

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	integer ndim
	parameter(ndim=nbdydim)

	integer iunit,it,n,nvers
	real xlag(ndim)
	real ylag(ndim)
	real zlag(ndim)

	logical oktime
	integer nrec
        integer ifemop

        iunit = ifemop('.lgr','unform','unknown')
        if( iunit .le. 0 ) stop
	!call lag_get_header(iunit,nvers)	!old
	call lag_get_header_new(iunit,nvers)
	write(6,*) 'lagrangian unit opened: ',iunit,nvers

	nrec = 0
	call timeset(0,0,0)
        call timeask
                        
    1   continue

	!call lag_get_xy(iunit,ndim,it,n,xlag,ylag,zlag)	!old
	call lag_get_xy_new(iunit,ndim,it,n,xlag,ylag,zlag)
	if( n .lt. 0 ) goto 2
        nrec = nrec + 1
        write(6,*) nrec,it,n

	if( oktime(it) ) then
	  write(6,*) 'plotting particles ',it,n
          call plo_part(n,xlag,ylag,zlag,'particles')
	end if

	goto 1
    2   continue

	end

c**********************************************************

        subroutine plo_xy(n,xlag,ylag,zlag)
     
        implicit none
       
	integer n
	real xlag(1)
	real ylag(1)
	real zlag(1)

	logical bcolor,bbottom,bzeta
	integer i
        real x,y,z,x1,y1,x2,y2
	real cwater,cbottom

c--------------------------------------
	bcolor = .true.		!use color to plot lagrangian particles
	bbottom = .false.	!use color to indicate p. on bottom
	bzeta = .true.		!use zeta to decide on color to use

	cwater = 1.0		!color to use for particles in water
	cbottom = 0.3		!color to use for particles on bottom
c--------------------------------------

        call qlwidth(0.035)
	if( bcolor ) call qhue(cwater)		!LAG_COLOR

	do i=1,n

	   x = xlag(i)
	   y = ylag(i)
	   z = zlag(i)

	   if( bcolor ) then
	     if( bzeta ) then
	       call qhue(z)
	     else if( bbottom ) then
	       if( z .lt. 1.0 ) then		! in water column
	  	call qhue(cwater)
	       else				! on bottom
	  	call qhue(cbottom)
	       end if
	     end if
	   end if

           x1=x+0.1
           y1=y+0.1
           x2=x-0.1
           y2=y+0.1        
           call qline(x1,y1,x2,y2)
           x1=x-0.1
           y1=y-0.1
           x2=x-0.1
           y2=y+0.1
           call qline(x1,y1,x2,y2)
           x1=x-0.1
           y1=y-0.1
           x2=x+0.1
           y2=y-0.1
           call qline(x1,y1,x2,y2)
           x1=x+0.1
           y1=y+0.1
           x2=x+0.1
           y2=y-0.1
           call qline(x1,y1,x2,y2)           

	end do

	call qgray(0.)

       end

c**********************************************************

