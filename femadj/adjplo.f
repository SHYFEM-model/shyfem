c
c $Id: adjplo.f,v 1.6 2007-03-20 13:19:42 georg Exp $
c
c description :
c
c plotting routines
c
c contents :
c
c subroutine plobas(nkn,nel,nen3v,ngrade,xgv,ygv)
c			plots basin
c subroutine pltgrd(n,k,k1,k2,xgv,ygv)
c			plots grade in color
c subroutine pltsgrd(igr,nkn,nel,nen3v,ngrade,xgv,ygv)
c			plots only special grade
c subroutine write_grid(file,nkn,nel,xgv,ygv,nen3v)
c			writes quick and dirty results to file
c
c revision log :
c
c 19.05.2003    ggu     plot some more info on debug plot
c
c***********************************************************

	subroutine plobas

c plots basin

	use mod_adj_grade
	use basin

	implicit none

	include 'param.h'

	integer ie,ii,k,n,i1,i2
	real xmin,xmax,ymin,ymax

	call mima(xgv,nkn,xmin,xmax)
	call mima(ygv,nkn,ymin,ymax)

	call qstart

	call qworld(xmin,ymin,xmax,ymax)
	call qrcfy

	do ie=1,nel
	  k = nen3v(3,ie)
	  call qmove(xgv(k),ygv(k))
	  do ii=1,3
	    k = nen3v(ii,ie)
c	    if( k .gt. 0 ) call qplot(xgv(k),ygv(k))	!HACK
	    call qplot(xgv(k),ygv(k))
	  end do
	end do

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    n = ngrade(k)
	    if( n .ne. 6 .and. n .gt. 0 ) then
		i1 = mod(ii,3) + 1
		i2 = mod(i1,3) + 1
		call pltgrd(n,k,nen3v(i1,ie),nen3v(i2,ie),xgv,ygv)
	    end if
	  end do
	end do

	call qend

	end

c***********************************************************

	subroutine pltgrd(n,k,k1,k2,xgv,ygv)

c plots grade in color

	implicit none

	integer n,k,k1,k2
	real xgv(1),ygv(1)
	real x(3),y(3)
	real col

	if( n .lt. 6 ) col = 0.35
	if( n .lt. 5 ) col = 0.15
	if( n .gt. 6 ) col = 0.65
	if( n .gt. 7 ) col = 0.85
	if( n .gt. 8 ) col = 0.95

	x(1) = xgv(k)
	y(1) = ygv(k)
	x(2) = xgv(k) + 0.3 * (xgv(k1)-xgv(k))
	y(2) = ygv(k) + 0.3 * (ygv(k1)-ygv(k))
	x(3) = xgv(k) + 0.3 * (xgv(k2)-xgv(k))
	y(3) = ygv(k) + 0.3 * (ygv(k2)-ygv(k))

	call qhue(col)
	call qafill(3,x,y)

	end

c***********************************************************

	subroutine pltsgrd(igr,nkn,nel,nen3v,ngrade,xgv,ygv)

c plots only special grade

	implicit none

	integer igr
	integer nkn,nel
	integer nen3v(3,1)
	integer ngrade(1)
	real xgv(1),ygv(1)

	integer k,n

	do k=1,nkn
	  n = ngrade(k)
	  if( n .ne. igr ) then
	    ngrade(k) = ngrade(k) - 1000
	  end if
	end do
	call plobas
	do k=1,nkn
	  n = ngrade(k)
	  if( n .ne. igr ) then
	    ngrade(k) = ngrade(k) + 1000
	  end if
	end do

	end

c***********************************************************

	subroutine wr0grd

c writes quick and dirty results from basin to file

	implicit none

	call write_grid('new0.grd')

	end

c***********************************************************

	subroutine write_grid(file)

c writes quick and dirty results to file

	implicit none

	character*(*) file

	call basin_to_grd
	call grd_write(file)

	end

c***********************************************************
c***********************************************************
c***********************************************************

	subroutine plosel2(ie1,ie2)

c plots element

	use mod_adj_grade
	use basin

	implicit none

	include 'param.h'

	integer ie1,ie2

	character*11 line
	integer ie,ii,k,n
	integer iev(2),iee,i,kk,id
	real xmin,xmax,ymin,ymax
	real dx,dy
	real x,y,xc,yc

	integer ialfa
	iev(1) = ie1
	iev(2) = ie2

	xmin = 1.e+35
	xmax = -xmin
	ymin = xmin
	ymax = -ymin
	x = 0.
	y = 0.

	do iee=1,2
	  ie = iev(iee)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    n = ngrade(k)
	    do i=1,n
	      kk = ngri(i,k)
	      x = xgv(kk)
	      y = ygv(kk)
	      if( x .gt. xmax ) xmax = x
	      if( x .lt. xmin ) xmin = x
	      if( y .gt. ymax ) ymax = y
	      if( y .lt. ymin ) ymin = y
	    end do
	  end do
	end do

	dx = xmax - xmin
	dy = ymax - ymin
	xmin = xmin - 0.2*dx
	xmax = xmax + 0.2*dx
	ymin = ymin - 0.2*dy
	ymax = ymax + 0.2*dy

	call qstart

	call qworld(xmin,ymin,xmax,ymax)
	call qrcfy

	do iee=1,2
	  ie = iev(iee)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    xc = xgv(k)
	    yc = ygv(k)
	    id = ialfa(float(k),line,-1,-1)
	    call qtext(xc,yc,line)
	    n = ngrade(k)
	    do i=1,n
	      kk = ngri(i,k)
	      x = xgv(kk)
	      y = ygv(kk)
	      call qline(xc,yc,x,y)
	      id = ialfa(float(kk),line,-1,-1)
	      call qtext(x,y,line)
	    end do
	  end do
	end do

	call qend

	end

c***********************************************************

	subroutine plosno(k)

c plots node and neighborhood

	use mod_adj_grade
	use basin

	implicit none

	include 'param.h'

	integer k

	integer n,ip,i,kk
	real xmin,xmax,ymin,ymax

	call mnmx0(k,xmin,xmax,ymin,ymax)

	call qstart

	call qworld(xmin,ymin,xmax,ymax)
	call qrcfy

        n = ngrade(k)

	call plonn(k)

	do i=1,n
	   kk = ngri(i,k)
	   call plosno0(kk)
	   call plonn(kk)
	end do

	call qend

	end

c***********************************************************

	subroutine plosno0(k)

c plots one node

	use mod_adj_grade
	use basin

	implicit none

	include 'param.h'

	integer k

	integer n,ip,i,kk
	real x,y,xc,yc

	xc = xgv(k)
	yc = ygv(k)
        n = ngrade(k)

	do i=1,n
	   kk = ngri(i,k)
	   x = xgv(kk)
	   y = ygv(kk)
	   call qline(xc,yc,x,y)
	end do

	end

c***********************************************************

	subroutine plonn(k)

c plots node number

	use mod_adj_grade
	use basin

	implicit none

	include 'param.h'

	integer k

	character*11 line
	logical isbound
	integer id,idnew,idtot
	integer ng
	real x,y
	integer ialfa

	line = ' '

	isbound = nbound(k) .eq. 1
        ng = ngrade(k)
	id = ialfa(float(k),line,-1,-1)
	idtot = id

	if( isbound ) then
	  line(id+1:id+4) = ' (B)'
	  idtot = id + 4
	else if( ng .ne. 6 ) then
	  line(id+1:id+2) = ' ('
	  idnew = ialfa(float(ng),line(id+3:),-1,-1)
	  idtot = id + 3 + idnew
	  line(idtot:idtot) = ')'
	end if

	x = xgv(k)
	y = ygv(k)
	call qtext(x,y,line(1:idtot))

	end

c******************************************************

	subroutine mnmx0(k,xmin,xmax,ymin,ymax)

c computes min/max of (x,y) of node k and neighbors

	use mod_adj_grade

	implicit none

	include 'param.h'

	integer k
	real xmin,xmax,ymin,ymax

	integer n,ip,i,kk
	real dx,dy

	xmin = 1.e+35
	xmax = -xmin
	ymin = xmin
	ymax = -ymin

        n = ngrade(k)

	call mnmx(k,xmin,xmax,ymin,ymax)

	do i=1,n
	   kk = ngri(i,k)
	   call mnmx(kk,xmin,xmax,ymin,ymax)
	end do

	dx = xmax - xmin
	dy = ymax - ymin
	xmin = xmin - 0.2*dx
	xmax = xmax + 0.2*dx
	ymin = ymin - 0.2*dy
	ymax = ymax + 0.2*dy

	end

c******************************************************

	subroutine mnmx(k,xmin,xmax,ymin,ymax)

c computes min/max of (x,y) of node k
c xmin... must be already initialized

	use mod_adj_grade
	use basin

	implicit none

	include 'param.h'

	integer k
	real xmin,xmax,ymin,ymax

	integer n,ip,i,kk
	real x,y

        n = ngrade(k)

	do i=1,n
	   kk = ngri(i,k)
	   x = xgv(kk)
	   y = ygv(kk)
	   if( x .gt. xmax ) xmax = x
	   if( x .lt. xmin ) xmin = x
	   if( y .gt. ymax ) ymax = y
	   if( y .lt. ymin ) ymin = y
	end do

	end

c******************************************************

