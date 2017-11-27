c
c $Id: zinit.f,v 1.5 2009-04-07 10:43:57 georg Exp $
c
	program zinit

c prepares initial z condition

	use basin

        include 'param.h'

c
	integer iapini
c
	include 'simul.h'

	real, allocatable :: hv(:)
c
c----------------------------------------------------------------
c k1,k2		nodes that define x/y of node 1/2
c z1,z2		initial z values at nodes 1/2
c mode		0: nothing  1: linear  2: cosine
c
c	data k1,k2 / 5,221 /
c	data k1,k2 / 5,167 /
c	data k1,k2 / 5,317 /
c	data z1,z2 / -.30,+.30 /
c	data z1,z2 / -2.30,+2.30 /
c	data z1,z2 / -1.30,+1.30 /
c	data mode /2/
c	data z1,z2 / .10,0. /
c	data mode /1/
c	data z1,z2 / -0.0475,+0.0475 /
c	data mode /1/
c
c bas092 -----------------------------
	data k1,k2 / 5,221 /
	data z1,z2 / -.30,+.30 /
	data mode /2/
c venlag61 -----------------------------
c	data k1,k2 / 4289,4340 /
c	data z1,z2 / -.40,+.40 /
c	data mode /2/
c----------------------------------------------------------------
c
	if(iapini(1,0,0,0).eq.0) stop

	allocate(hv(nkn))
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	pi=acos(-1.)
c
	kk1 = ipext(k1)
	kk2 = ipext(k2)
	if( kk1 .le. 0 ) stop 'error stop: no node'
	if( kk2 .le. 0 ) stop 'error stop: no node'

	x1=xgv(kk1)
	y1=ygv(kk1)
	x2=xgv(kk2)
	y2=ygv(kk2)
	rl = (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1)
	zm = (z1+z2)/2.
	za = (z1-z2)/2.
c
	write(6,*) '------------ (kext,kint,z,x,y) -------------'
	write(6,1000) 'first  node : ',k1,kk1,z1,x1,y1
	write(6,1000) 'second node : ',k2,kk2,z2,x2,y2
 1000	format(1x,a,2i7,3f12.3)
c
	do i=1,nkn
	  xi=xgv(i)
	  yi=ygv(i)
c	  r varies linearily from 0 to 1 (r is distance from (x1,y1))
	  r = ( (xi-x1)*(x2-x1) + (yi-y1)*(y2-y1) ) / rl
	  if(mode.eq.1) then
c	    hv is linear -- hv(r=0)=z1 , hv(r=1)=z2
	    hv(i) = z1 + (z2-z1)*r
	  else if(mode.eq.2) then
c	    hv is cosinus -- hv(r=0)=z1 , hv(r=1)=z2 , hv(r=1/2)=(z1+z2)/2
	    hv(i) = zm  +  za * cos( r*pi )
	  else
c	    hv is average 
	    hv(i)=(z1+z2)*0.5
	  end if
	end do

	write(6,*) (hv(i),i=1,nkn)

	open(55,file='new.ini',status='unknown',form='unformatted')
	write(55) (hv(i),i=1,nkn)
	close(55)

	write(6,*) 'data written to file new.ini'

	end
