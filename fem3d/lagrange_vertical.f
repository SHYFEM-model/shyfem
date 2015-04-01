c
c deal with vertical velocities
c
c******************************************************

	subroutine getzvel(ie,z0,l0,is_in,is_out,a_in,a_out,w)

c returns vertical velocity to be used in lagrangian model

	implicit none

	integer ie	!element number
	real z0		!depth of particle
	integer l0	!layer in which particle is in [1-nlv]
	integer is_in	!side from which particle comes [1-3]
	integer is_out	!side to which particle goes [1-3]
	real a_in	!relative postion of side is_in [0-1]
	real a_out	!relative postion of side is_out [0-1]
	real w		!computed vertical velocity for ie and l0 (return)

	integer ii
	real win,wout

!	-----------------------------------------------
!	velocity at entering and exiting point
!	-----------------------------------------------

	if( ie .le. 0 ) then
	  w = 0.
	  return
	end if

	ii = mod(is_in,3) + 1
	write(6,*) 'gguuut: ',ie,l0,ii,a_in
	call getzvel_point(ie,l0,ii,a_in,win)

	ii = mod(is_out,3) + 1
	call getzvel_point(ie,l0,ii,a_out,wout)

!	-----------------------------------------------
!	average velocity in element
!	-----------------------------------------------

	w = 0.5 * ( win + wout )

!	-----------------------------------------------
!	end of routine
!	-----------------------------------------------

	end

c******************************************************

	subroutine getzvel_point(ie,l0,is,a,w)

c returns vertical velocity at point given by ie,l0,is,a

	implicit none

	integer ie	!element number
	integer l0	!layer in which particle is in
	integer is	!side where particle is [1-3]
	real a		!relative postion of side
	real w		!computed vertical velocity for ie and l0 (return)

	include 'param.h'
	include 'basin.h'
	include 'hydro_vel.h'

	integer ii,k1,k2
	real wo1,wo2,wn1,wn2
	real w1,w2

!	-----------------------------------------------
!	velocity at first node
!	-----------------------------------------------

	ii = is

	k1 = nen3v(ii,ie)
	wo1 = wlov(l0,k1) + wlov(l0-1,k1)
	wn1 = wlnv(l0,k1) + wlnv(l0-1,k1)
	w1 = 0.25 * ( wo1 + wn1 )

!	-----------------------------------------------
!	velocity at second node
!	-----------------------------------------------

	ii = mod(ii,3) + 1

	k2 = nen3v(ii,ie)
	wo2 = wlov(l0,k2) + wlov(l0-1,k2)
	wn2 = wlnv(l0,k2) + wlnv(l0-1,k2)
	w2 = 0.25 * ( wo2 + wn2 )

!	-----------------------------------------------
!	interpolate velocity
!	-----------------------------------------------

	w = (1.-a)*w1 + a*w2

!	-----------------------------------------------
!	end of routine
!	-----------------------------------------------

	end

c******************************************************

	subroutine lagr_layer_thickness(ie,nlv,hl)

c computes layer thickness for element ie

	implicit none

	integer ie		!element number
	integer nlv		!max number of layers
	real hl(nlv)		!layer thickness (return)

	include 'param.h'
	include 'levels.h'
	include 'depth.h'
	include 'hydro.h'

	integer lmax,nsigma,ii
	real hsigma
	real z,h

        !call compute_sigma_info(lmax,hlv,nsigma,hsigma)
	call get_sigma_info(lmax,nsigma,hsigma)
	if( lmax > nlv ) goto 99

	h = hev(ie)
	z = 0.
	do ii=1,3
	  z = z + zenv(ii,ie) + zeov(ii,ie)
	end do
	z = z / 6.

        call get_layer_thickness(lmax,nsigma,hsigma,z,h,hlv,hl)

	return
   99	continue
	write(6,*) 'nlv,lmax: ',nlv,lmax
	stop 'error stop lagr_layer_thickness: incompatible layers'
	end

c*******************************************************************

        subroutine vertpos(zn0,deltat,layd,w,zn1,ztime,addl)

	implicit none

	real zn0		!initial relative vertical position in layer
	real deltat		!time spent in element(horizontal)
	real w			!average vertical velocity
	real zn1		!final relative vertical position in layer
	real ztime		!time spent in element (vertical)
        integer addl		!relative movement to next layer [-1,0,+1]
        real layd		!layer thickness

c return are zn1,ztime,addl

c lstd		max vertical distance to be traveled
c cpstd		computed distance to be traveled

	real cpdst,lstd

        if(w.gt.0) lstd=layd*(1-zn0)
        if(w.lt.0) lstd=-layd*zn0
        if(w.eq.0) lstd=0

        cpdst=deltat*w
        zn1=zn0+(cpdst/layd)
        
	if( w == 0 ) then
	  ztime = 2*deltat
	else
          ztime=lstd/w
	end if

        if(ztime.le.deltat)then
          if(w.gt.0)addl=1
          if(w.lt.0)addl=-1
        else
          addl=0
        end if
	
	if(zn1.gt.1)zn1=0
        if(zn1.lt.0)zn1=1

        end

c************************************************************

	subroutine getalfa(side,d,near,far,re)

	implicit none

        integer side		!number of side of element [1-3]
	real d			!relative distance from closest vertex
	integer near		!number of closest vertex [1-3]
	integer far		!to be deleted...
	real re			!relative distance to first point of side

	integer ind

        ind=mod(side,3)+1
	
	if(ind.eq.near)then
	  re=d
	else
	  re=1-d
	end if
	
	end 

c************************************************************

