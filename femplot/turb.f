
c*************************************************************

	programm plotturb

c plots turbidity

	implicit none

c parameters
	integer ntdim,nlvdim
	parameter (ntdim=50000,nlvdim=20)

	integer nt,nlv
	integer iunit
	real a(nlvdim,ntdim)

	iunit = 5

	call rdfile(iunit,ntdim,nlvdim,nt,nlv,a)

	write(6,*) 'nt,nlv : ',nt,nlv

	call qopen
	call turb(nt,nlv,nlvdim,a)
	call qclose

c----------------------------------------------
c end of routine
c----------------------------------------------

	end

c*****************************************************************

	subroutine turb(nt,nlv,nlvdim,a)

c plots turbulence

	implicit none

	integer nt,nlv,nlvdim
	real a(nlvdim,1)

	integer it,l
	real am,at,col,acaux
	real colmin,colmax
	real atmin,atmax
	real xmin,ymin,xmax,ymax
	real rx1,ry1,rx2,ry2
	real amin,amax,dm,dmax
	integer itmax

	colmin = 0.3
	colmax = 0.7

	am = 0.
	dmax = 0.

	do it=1,nt
	  amin = a(1,it)
	  amax = a(1,it)
	  do l=1,nlv
	    at = a(l,it)
	    am = max(am,at)
	    amax = max(amax,at)
	    amin = min(amin,at)
	  end do
	  dm = (amax-amin) / (amax+amin)
	  if( dm .gt. dmax ) then
	    dmax = dm
	    itmax = it
	  end if
	end do

	write(6,*) 'Maximum: ',am,dmax,itmax

	atmin = 0.
	atmax = am

	acaux = (colmax-colmin)/(atmax-atmin)

	call qstart()
	xmin = 0.
	ymin = 0.
	xmax = nt
	ymax = nlv
	call qworld(xmin,ymin,xmax,ymax)

	do it=1,nt
	  do l=1,nlv
	    at = a(l,it)
	    col = colmin + (at-atmin)*acaux

	    call qhue(col)
	    rx1 = it-1
	    rx2 = it
	    ry1 = l-1
	    ry2 = l
	    call qrfill(rx1,ry1,rx2,ry2)
	  end do
	end do

	call qend()

	end

c*****************************************************************

	subroutine rdrec(iunit,nlvdim,nlv,a,ier)

c reads one record

	implicit none

	integer iunit
	integer nlvdim,nlv
	real a(1)
	integer ier

	integer it,i
	real rit

	read(iunit,*,end=99) rit,nlv
	if( nlv .gt. nlvdim ) stop 'error stop rdrec: nlvdim'
	read(iunit,*) (a(i),i=1,nlv)

	ier = 0

	return
   99	continue
	ier = -1
	return
	end

c*****************************************************************

	subroutine rdfile(iunit,ntdim,nlvdim,nt,nlv,a)

c reads file

	implicit none

	integer iunit
	integer ntdim,nlvdim,nt,nlv
	real a(nlvdim,ntdim)

	integer it,ier

	it = 0

    1	continue
	  it = it + 1
	  if( it .gt. ntdim) stop 'error stop rdfile: ntdim'
	  call rdrec(iunit,nlvdim,nlv,a(1,it),ier)
	  if( ier .ne. 0 ) goto 2
	goto 1
    2	continue

	nt = it - 1

	return
   99	continue
	ier = -1
	return
	end

c*****************************************************************

