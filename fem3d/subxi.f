
!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2018  Georg Umgiesser
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

c xi routines
c
c contents :
c
c
c revision log :
c
c 21.01.2015	ggu	new routine compute_cartesian_coords()
c 23.04.2015	ggu	compute s value for particle path
c 16.05.2015	ccf	add some more check in xit_start_end for small numbers
c
c***********************************************************
c***********************************************************
c***********************************************************
c***********************************************************

	subroutine xy2xit(xx,yy,x,y,xi)

c given x/y returns internal coordinates xi

	implicit none

	double precision xx(3),yy(3)
	double precision x,y
	double precision xi(3)

	integer ii
	double precision a(3),b(3),c(3)

	call xit_abc(xx,yy,a,b,c)

	do ii=1,3
	  xi(ii) = a(ii) + b(ii)*x + c(ii)*y
	end do

	end
	
c***********************************************************

	subroutine xit2xy(xx,yy,x,y,xi)

c given internal coordinates xi returns x/y

	implicit none

	double precision xx(3),yy(3)
	double precision x,y
	double precision xi(3)

	double precision eps
	parameter (eps = 1.e-8)

	integer ii,i1,i2,iimax
	double precision det,detmax
	double precision a(3),b(3),c(3)
	double precision aa(3)

	call xit_abc(xx,yy,a,b,c)

	iimax = 0
	detmax = 0.
	do ii=1,3
	  aa(ii) = (xi(ii)-a(ii))
	  i1 = mod(ii,3) + 1
	  i2 = mod(i1,3) + 1
	  det = b(i2)*c(i1) - b(i1)*c(i2)
	  if( abs(det) > detmax ) then
	    detmax = det
	    iimax = ii
	  end if
	end do

	!det = b(2)*c(1) - b(1)*c(2)
	!if( abs(det) < eps ) goto 99
	!y = (aa(1)*b(2) - aa(2)*b(1)) / det
	!x = -(aa(1)*c(2) - aa(2)*c(1)) / det

	ii = iimax
	i1 = mod(ii,3) + 1
	i2 = mod(i1,3) + 1
	det = b(i2)*c(i1) - b(i1)*c(i2)
	if( abs(det) < eps ) goto 99
	y = (aa(i1)*b(i2) - aa(i2)*b(i1)) / det
	x = -(aa(i1)*c(i2) - aa(i2)*c(i1)) / det

	return
   99	continue
	write(6,*) 'det too small: ',det
	write(6,*) x,y
	write(6,*) xi
	write(6,*) a
	write(6,*) b
	write(6,*) c
	stop 'error stop xi2xy: det == 0'
	end
	
c***********************************************************

	subroutine xit_abc(xx,yy,a,b,c)

c returns a,b,c to compute xi
c
c natural coordinates in triangle:   xi(i) = a(i) + b(i)*x + c(i)*y    i=1,3

	implicit none

	double precision xx(3),yy(3)
	double precision a(3),b(3),c(3)

	double precision x1,y1,x2,y2,x3,y3
	double precision a1,a2,a3,aj,aji

  	x1=xx(1)
	y1=yy(1)
	x2=xx(2)
	y2=yy(2)
	x3=xx(3)
	y3=yy(3)

	a1=x2*y3-x3*y2
	a2=x3*y1-x1*y3
	a3=x1*y2-x2*y1
	aj = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)
	aji=1./aj

	a(1)=a1*aji
	a(2)=a2*aji
	a(3)=a3*aji
	b(1)=(y2-y3)*aji
	c(1)=(x3-x2)*aji
	b(2)=(y3-y1)*aji
	c(2)=(x1-x3)*aji
	b(3)=(y1-y2)*aji
	c(3)=(x2-x1)*aji

	end

c***********************************************************

	subroutine xit_check(xx,yy,area,reg)

c checks triangle for regularity
c
c natural coordinates in triangle:   xi(i) = a(i) + b(i)*x + c(i)*y    i=1,3

	implicit none

	double precision xx(3),yy(3)
	double precision area,reg

	integer ii,iii
	double precision x1,y1,x2,y2,x3,y3
	double precision aj
	double precision circ,creg

  	x1=xx(1)
	y1=yy(1)
	x2=xx(2)
	y2=yy(2)
	x3=xx(3)
	y3=yy(3)

	circ = 0.
	do ii=1,3
	  iii = mod(ii,3) + 1
	  circ = circ + sqrt((xx(ii)-xx(iii))**2+(yy(ii)-yy(iii))**2)
	end do

	aj = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)
	area = 0.5*aj

	creg = 0.
	if( area > 0. ) then
	  creg = 6. * sqrt( area/sqrt(3.d0) )	!lenght of regular triangle
	end if
	
	reg = creg/circ

	end

c***********************************************************

	subroutine xit_start_end(iflux,alpha,xip,xis,xie,s)

c find start and end point of line through point xip
c
c natural coordinates in triangle:   xi(i) = a(i) + b(i)*x + c(i)*y    i=1,3

	implicit none

	integer iflux			!index of flux point [1-3]
	double precision alpha		!ratio of flux exit at side iflux [0-1]
	double precision xip(3)		!coordinates of point
	double precision xis(3)		!start coordinate (return)
	double precision xie(3)		!end coordinate (return)
	double precision s		!point parameter [0-1] (return)

	integer i1,i2,i3,ii
	integer iflag,ielem
	double precision as,ae
	double precision beta,xi3,a,gamma
	double precision, parameter  :: small=1.d-12

	xis = 0.d0
	xie = 0.d0

	gamma = -1.d0			!not yet computed
	i1 = iflux
	xi3 = alpha*(1.d0-xip(i1))
	i3 = mod(i1+1,3) + 1

	if( xi3 > xip(i3) ) then	!right sub-triangle
	  ielem = 1
	  a = alpha
	  i2 = mod(i1,3) + 1
	  i3 = mod(i2,3) + 1
	  if( xip(i3) == 0.d0 ) gamma = 1.d0
	  if( xip(i1) == 0.d0 ) gamma = 0.d0
	else				!left sub-triangle
	  ielem = 2
	  a = 1.d0 - alpha
	  i3 = mod(i1,3) + 1
	  i2 = mod(i3,3) + 1
	  if( xip(i3) == 0.d0 ) gamma = 1.d0
	  if( xip(i1) == 0.d0 ) gamma = 0.d0
	end if

	if( a .eq. 0.d0 ) then
	  beta = 0.d0
	else if( xip(i3) > a ) then
	  write(6,*) iflux,alpha
	  write(6,*) i1,i2,i3,a
	  write(6,*) xip
	  stop 'error stop xit_start_end: internal error (2)'
	else
	  beta = 1.d0 - xip(i1) - xip(i3)/a
	end if

	beta = dmin1(beta,1.0d0)
        if ( beta < small ) beta = 0.d0

	as = beta
	ae = a*(1.d0-beta)

        xis(i1) = 1.d0-as
        xis(i2) = as
        xis(i3) = 0.d0

        xie(i1) = 0.d0
        xie(i2) = 1.d0-ae
        xie(i3) = ae

	if( gamma >= 0.d0 ) then	!already computed
	  !nothing
	  iflag = 0
	else if( beta == 1.d0 ) then
	  gamma = 0.d0
	  iflag = 1
	else if( xip(i1) == 0.d0 ) then
	  gamma = 0.d0
	  iflag = 2
	else
	  gamma = xip(i1)/(1.d0-beta)
	  iflag = 3
	end if

	gamma = dmin1(gamma,1.0d0)
        if ( gamma < small ) gamma = 0.d0

        s = 1.d0 - gamma
	if( gamma .eq. 1.d0 ) s = 0.d0

	if( s < 0.d0 .or. s > 1.d0 ) then
	  write(6,*) 'error computing s: ',s,gamma
	  write(6,*) alpha,a,beta
	  write(6,*) i1,iflag,ielem
	  write(6,*) xip
	  write(6,*) xis
	  write(6,*) xie
	  stop 'error stop xit_start_end: internal error (1)'
	end if

	end

c***********************************************************

