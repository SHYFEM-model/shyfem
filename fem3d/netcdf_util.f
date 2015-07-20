c
c $Id: netcdf_util.f,v 1.15 2009-11-18 16:50:37 georg Exp $
c
c utilities for netcdf conversion
c
c revision log :
c
c 21.01.2013    ggu     routines transfered from ous2nc.f
c 20.02.2013    ggu     new routines get_period() and check_period()
c
c******************************************************************

        subroutine write_time(it)

        implicit none

        integer it

        character*40 line

        call dtsgf(it,line)
        write(6,*) 'time: ',it,'   ',line

        end

c******************************************************************

	subroutine read_date_and_time(date,time)

	implicit none

	integer date,time

	integer ierr
	character*20 line

	write(6,*)
	write(6,*) 'You can specify a date for fem-time 0'
	write(6,*) 'This is only used if the simulation does not'
	write(6,*) 'contain a date and time stamp.'
	if( date > 0 ) then
	  call dtsunform_pack(date,time,line,ierr)
	  write(6,*) 'Default for date: ',line
	else
	  write(6,*) 'There is no default date... using 0'
	end if
	write(6,*) 'Format: YYYY-MM-DD[::dd:mm:ss]'
	write(6,*) 'Enter date: (return for default)'

	read(5,'(a)') line

	if( line == ' ' ) return

	call dtsunform_pack(date,time,line,ierr)
	if( ierr .ne. 0 ) goto 99

	call dtsform_pack(date,time,line)
	write(6,*) 'Chosen date: ',line

	return
   99	continue
	write(6,*) 'cannot parse date format: ',line
	stop 'error stop read_date_and_time: date format'
	end

c******************************************************************

	subroutine get_period(iperiod,its,ite,nfreq)

	implicit none

	integer iperiod		! type of period (0 for none) (return)
	integer its,ite		! time limit (start, end) (return)
	integer nfreq		! frequency of output (return)

	integer n
	integer dates,datee
	integer year,month,day,hour,min,sec
	integer ierr
	character*80 line
	real f(10)
	double precision d(10)

	integer iscanf,iscand

	iperiod = 0
	its = -1
	ite = -1

	write(6,*) 'Do you want to specify period of extraction?'
	write(6,*) '  0 or return  all of file'
	write(6,*) '  1            give start, end and frequency'
	read(5,'(a)') line
	n = iscanf(line,f,1)

	if( n .le. 0 ) f(1) = 0.
	iperiod = f(1)
	if( iperiod .le. 0 ) return

	if( n .eq. 1 ) then
	  write(6,*) '  format for date is: YYYY-MM-DD[::dd:mm:ss]'

	  write(6,*) '  Enter start date (return for start of data):'
	  read(5,'(a)') line
	  if( line .ne. ' ' ) then
	    call dtsunform(year,month,day,hour,min,sec,line,ierr)
	    if( ierr .ne. 0 ) goto 99
	    call dts2it(its,year,month,day,hour,min,sec)
	  end if

	  write(6,*) '  Enter end date (return for end of data):'
	  read(5,'(a)') line
	  if( line .ne. ' ' ) then
	    call dtsunform(year,month,day,hour,min,sec,line,ierr)
	    if( ierr .ne. 0 ) goto 99
	    call dts2it(ite,year,month,day,hour,min,sec)
	  end if

	  write(6,*) '  Enter frequency (return for every record):'
	  read(5,'(i10)') nfreq
	end if

	return
   99	continue
	write(6,*) 'cannot parse date format: ',line
	stop 'error stop get_period: date format'
	end

c******************************************************************

	subroutine check_period(it,iperiod,its,ite,nfreq,bwrite)

	implicit none

	integer it		! fem time
	integer iperiod		! type of period (0 for none)
	integer its,ite		! time limit (start, end)
	integer nfreq		! frequency of output
	logical bwrite		! write output? (return)

	integer icall
	save icall
	data icall / 0 /

	icall = icall + 1

	bwrite = .true.
	if( iperiod .le. 0 ) return

	bwrite = .false.
	if( its .ne. -1 .and. it .lt. its ) return
	if( ite .ne. -1 .and. it .gt. ite ) return

	if( mod(icall,nfreq) .ne. 0 ) return

	bwrite = .true.

	end

c******************************************************************

	subroutine get_lmax_reg(nx,ny,fm,ilhv,lmax)

c computes max lmax for regular domain

	implicit none

	integer nx,ny
	real fm(4,nx,ny)
	integer ilhv(1)
	integer lmax		!max level (return)

	integer i,j,ie

	lmax = 0

	do j=1,ny
	  do i=1,nx
	    ie = nint(fm(4,i,j))
	    lmax = max(lmax,ilhv(ie))
	  end do
	end do

	end

c******************************************************************
c******************************************************************
c******************************************************************

	subroutine get_dimensions(nxdim,nydim,nx,ny,x0,y0,dx,dy,xlon,ylat)

c gets dimensions for reguar grid

	use basin

	implicit none

	include 'param.h'

	integer nxdim,nydim,nx,ny
	real x0,y0,dx,dy
	real xlon(nxdim)
	real ylat(nydim)



	integer i,n
	real x1,y1,dxy
	real xmin,ymin,xmax,ymax
	real f(10)
	character*80 line

	integer iscanf
	real rnext

	call mima(xgv,nkn,xmin,xmax)
	call mima(ygv,nkn,ymin,ymax)

	write(6,*) 'min/max of domain:'
	write(6,*) 'xmin/xmax: ',xmin,xmax
	write(6,*) 'ymin/ymax: ',ymin,ymax

	nx = min(100,nxdim)
	ny = min(100,nydim)
	dx = (xmax - xmin) / nx
	dy = (ymax - ymin) / ny
	dxy = max(dx,dy)
	dxy = rnext(dxy,3)
	x0 = xmin
	y0 = ymin
	x1 = xmax
	y1 = ymax

	!write(6,*) 'Default for dx/dy: ',dx,dy,dxy
	!write(6,*) 'Default for dx/dy: ',dxy

	write(6,*)
	write(6,*) 'Enter dx[,dy]: (return for unstructured output)'
	read(5,'(a)') line
	n = iscanf(line,f,2)
	if( n .le. 0 ) then
	  dx = 0.
	  dy = 0.
	else if( n .eq. 1 ) then
	  dx = f(1)
	  dy = dx
	else
	  dx = f(1)
	  dy = f(2)
	end if
	write(6,*) 'dx,dy: ',dx,dy

	if( dx .le. 0. .or. dy .le. 0. ) then
	  nx = 0
	  ny = 0
	  return
	end if

	x0 = dx * (int(xmin/dx))
	y0 = dy * (int(ymin/dy))
	x1 = dx * (int(xmax/dx)+1)
	y1 = dy * (int(ymax/dy)+1)
	nx = 1 + nint((x1-x0)/dx)
	ny = 1 + nint((y1-y0)/dy)

	write(6,*) 'limits of domain:'
	write(6,*) 'dx,dy      : ',dx,dy
	write(6,*) 'x0,y0,x1,y1: ',x0,y0,x1,y1
	write(6,*) 'nx,ny      : ',nx,ny

	do while( .true. )
	  write(6,*)
	  write(6,*) 'Enter x0,y0,x1,y1: (return for default)'
	  read(5,'(a)') line
	  n = iscanf(line,f,4)
	  if( n.eq. 0 .or. n .eq. 4 ) exit
	  write(6,*) 'Please either accept default or enter 4 values.'
	end do

	if( n .eq. 4 ) then
	  x0 = f(1)
	  y0 = f(2)
	  x1 = f(3)
	  y1 = f(4)
	end if

	nx = 1 + nint((x1-x0)/dx)
	ny = 1 + nint((y1-y0)/dy)
	x1 = x0 + (nx-1)*dx
	y1 = y0 + (ny-1)*dy

	write(6,*) 'Final parameters: '
	write(6,*) 'dx,dy: ',dx,dy
	write(6,*) 'x0,y0,x1,y1: ',x0,y0,x1,y1
	write(6,*) 'nx,ny: ',nx,ny

	if( nx .le. 0 ) goto 98
	if( ny .le. 0 ) goto 98
	if( nx .gt. nxdim ) goto 99
	if( ny .gt. nydim ) goto 99

	do i=1,nx
	  xlon(i) = x0 + (i-1)*dx
	end do

	do i=1,ny
	  ylat(i) = y0 + (i-1)*dy
	end do

	return
   98	continue
	write(6,*) 'nx,ny: ',nx,ny
	stop 'error stop get_dimensions: error in nx,ny'
   99	continue
	write(6,*) 'nx,nxdim: ',nx,nxdim
	write(6,*) 'ny,nydim: ',ny,nydim
	write(6,*) 'please increase nxdim,nydim'
	stop 'error stop get_dimensions: nxdim/nydim'
	end

c******************************************************************

	subroutine write_dimensions(nx,ny,x0,y0,dx,dy)

	implicit none

	integer nx,ny
	real x0,y0,dx,dy

	real x1,y1

	x1 = x0 + nx*dx
	y1 = y0 + ny*dy

	write(6,*) 'Final parameters used: '
	write(6,*) 'dx,dy: ',dx,dy
	write(6,*) 'x0,y0,x1,y1: ',x0,y0,x1,y1
	write(6,*) 'nx,ny: ',nx,ny

	end

c******************************************************************

