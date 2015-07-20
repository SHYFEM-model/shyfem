c
c $Id: newcra.f,v 1.7 2008-07-16 15:41:39 georg Exp $
c
c routines for radiation condition
c
c contents :
c
c subroutine crador(v1,v2,mode,res)	computes wave celerity in triangle
c subroutine gwi(hia,hik,ibrad,b,c,zold,dt,h,area)
c					implements radiation condition (GWI)
c
c revision log :
c
c 23.03.2006    ggu     changed time step to real
c
c notes :
c
c can use only for non dry elements
c
c*****************************************************************
c
	subroutine crador(v1,v2,mode,res)
c
c computes wave celerity in triangle
c
	use mod_depth
	use mod_bnd_aux
	use mod_bound_dynamic
	use mod_hydro_baro
	use mod_hydro
	use evgeom
	use basin

	implicit none
c
c arguments
	integer mode
	real v1(1),v2(1)
	real res
c parameters
	integer ndim
	parameter (ndim=55)
c common
	include 'param.h'
	include 'pkonst.h'
	include 'femtime.h'
	include 'mkonst.h'
c local
	integer ie,ii,k,irand,n
	real dt,area,bcz,c,cmax
	integer ier,i,ix,iy,iok
	real x,y
	real z0,z2,z1n,z1o
	real cl,rbold,rbaux
	real*8 zaux1,zaux2
c functions
	integer iround,locsps
c
c imin,imax	first and last entry in arrays nb,zb,rb,iipp
c nb		internal node number of first three rows/columns
c		...0 is at boundary line, ...
c zb		for ORE is value of water level on first inner line
c		...at time level n-1 (new time level is n+1)
c		...for other rad.b.c. is used as an aux vector
c rb		for ORE contains the new water levels that have to
c		...be imposed and are copied to rzv with mode=2
c		...for other rad.b.c. is used as an aux vector
c iipp		for every boundary node gives internal element number
c		...of elements that have one side on the boundary line
c		...and contain the boundary node
c
c mode:
c	0	ORE with FEM, crad (wave celerity) is computed in every element
c	1	ORE with FD, sets up rb and saves old time level to zb (sp159)
c	2	ORE with FD, copies rb to rzv (spb11)
c	3	Schrimpf, use only elements from iipp
c	4	Schrimpf, compute transports for points on first inner line and
c		...use all elements that have the node in common, use this
c		...transport value to compute zeta at boundary
c	5	Schrimpf, as 4, but average zeta values on boundary
c	6	Schrimpf, use all elements that have one noe in common
c		...with boundary line, 4 elements are used to compute the
c		...transport at the boundary line (including one element
c		...that is not including the boundary node where zeta
c		...is to be computed)
c	7	Schrimpf implicit, elements are chosen as under 6
c	8	Schrimpf FEM, only elements are used that have node in common
c	9	Schrimpf as 6, but iterated to real time balance
c
	logical bfirst
	save bfirst
	integer imin,imax
	integer nb(0:2,ndim)
	real*8 epsa
	real zb(ndim),rb(ndim),dd(2*ndim)
	integer iipp(2,ndim)
	save nb,zb,rb,iipp,imin,imax,epsa
	data bfirst /.true./
	data epsa / 1.d-15 /
c
	if(bfirst) then
		bfirst=.false.
		write(6,*) '**** setup of crador ****'
		imin=10000
		imax=0
		do i=1,ndim
		  do ii=0,2
		    nb(ii,i)=0
		  end do
		end do
		do k=1,nkn
		  x=xgv(k)
		  y=ygv(k)
		  ix=iround(xgv(k)/3000.)+1
		  iy=iround(ygv(k)/3000.)
		  if(iy.le.2) then
		    if(ix.gt.ndim) then
			stop 'error stop : crador : ndim'
		    end if
		    if(ix.lt.imin) imin=ix
		    if(ix.gt.imax) imax=ix
		    nb(iy,ix) = k
		  end if
		end do
		write(6,*)' crador : ',imin,imax
		if(imax+1.gt.ndim) then	!for averaging
		  stop 'error stop : crador : imax+1'
		end if
		if(imin-1.lt.1) then	!for averaging
		  stop 'error stop : crador : imin-1'
		end if
		ier=0
		do i=imin,imax
		  do ii=0,2
		    if(nb(ii,i).eq.0) then
			ier=1
			write(6,*) 'error ',ii,i
		    end if
		  end do
		end do
		if(ier.ne.0) stop 'error stop : crador'
		do i=imin,imax
		  write(6,*) i,(nb(ii,i),ii=0,2)
		end do
		do i=imin,imax
		  zb(i)=zov(nb(1,i))
		  rb(i)=zb(i)
		end do

c		Schrimpf boundary condition --> set up iipp

		do ie=1,nel
		  i=0
		  do ii=1,3
		    y=ygv(nen3v(ii,ie))
		    iy=iround(y/3000.)
		    if(iy.eq.0) i=i+1
		  end do
		  if(i.eq.2) then	!two boundary nodes in element
		   do ii=1,3
		    x=xgv(nen3v(ii,ie))
		    y=ygv(nen3v(ii,ie))
		    ix=iround(x/3000.)+1
		    iy=iround(y/3000.)
		    if(iy.eq.0) then
		      if(iipp(1,ix).eq.0) then
			iipp(1,ix)=ie
		      else if(iipp(2,ix).eq.0) then
			iipp(2,ix)=ie
		      else
			write(6,*) iipp(1,ix),iipp(2,ix),ix,iy,ie
			stop 'error iipp'
		      end if
		    end if
		   end do
		  end if
		end do

c		double first and last pointer to avoid exeption

	        if(iipp(2,imin).ne.0) then
		  write(6,*) 'error iipp(2,imin) : ',imin,iipp(2,imin)
		  stop 'error stop : crador'
		end if
	        if(iipp(2,imax).ne.0) then
		  write(6,*) 'error iipp(2,imax) : ',imax,iipp(2,imax)
		  stop 'error stop : crador'
		end if
		iipp(2,imin)=iipp(1,imin)
		iipp(2,imax)=iipp(1,imax)

c		control

		do i=imin,imax
		  do ii=1,2
		    if(iipp(ii,i).eq.0) then
			ier=1
			write(6,*) 'error ',ii,i
		    end if
		  end do
		end do
		if(ier.ne.0) stop 'error stop : crador'

		do i=imin,imax
		  write(6,*) i,iipp(1,i),iipp(2,i)
		end do
	end if
c		
	if(mode.eq.0) then
c
c compute ORE with FEM - crad -> wave celerity
c
	cmax=3000./90.
c
	call get_timestep(dt)
c
	do k=1,nkn
	  v1(k)=0.
	  v2(k)=0.
	end do
c
	do ie=1,nel
	  bcz=0.
	  irand=0
	  do ii=1,3
	    k=nen3v(ii,ie)
	    if(rzv(k).ne.flag) irand=1	!could jump to end of outer do loop
	    bcz=bcz+ev(6+ii,ie)*zov(k)		!top or bottom (c) !OOO
c	    bcz=bcz+ev(6+ii,ie)*(znv(k)+zov(k))	!top or bottom (c)
c	    bcz=bcz+ev(3+ii,ie)*(znv(k)+zov(k))	!right or left (b)
	  end do
	  if(irand.eq.0) then
	    area=ev(10,ie)
	    do ii=1,3
	      k=nen3v(ii,ie)
	      v1(k)=v1(k)+area
	      v2(k)=v2(k)+area*bcz	!bottom or left
c	      v2(k)=v2(k)-area*bcz	!right or top
	    end do
	  end if
	end do
c
	do k=1,nkn
	  if(v2(k).ne.0) then
	    v1(k)=v1(k)*(znv(k)-zov(k))/(dt*v2(k))	!OOO
c	    v1(k)=v1(k)*(znv(k)-zov(k))/(0.5*dt*v2(k))
	  end if
	end do
c
	do ie=1,nel
	  n=0
	  c=0.
	  do ii=1,3
	    k=nen3v(ii,ie)
	    if(rzv(k).eq.flag) then
	      n=n+1
	      c=c+v1(k)
	    end if
	  end do
	  if(n.gt.0) c=c/n
	  if(c.lt.0.) then
	    c=0.
	  else if(c.gt.cmax) then
	    c=cmax
	  end if
	  crad(ie)=c
	end do
c
	else if(mode.eq.1) then
c
c compute ORE via FD
c
	do i=imin,imax        
	  z0=zov(nb(0,i))
	  z1n=znv(nb(1,i))
	  z1o=zb(i)
	  z2=zov(nb(2,i))
c                             
	  zaux1 = z1n+z1o-2.*z2
	  zaux2 = z1o-z1n
	  if(abs(zaux1).gt.epsa) then
	    cl = zaux2/zaux1
	  else
	    cl = 0.
	  end if
c
c	if(mod(i,20).eq.0) write(6,*) 'cl : ',cl,zaux1,zaux2
	  if( cl.lt.0. ) cl=0.
	  if( cl.gt.1. ) cl=1.
c
	  rb(i) = ( z0*(1.-cl) + 2.*cl*z1n ) / ( 1.+cl )
	end do
c
c	do i=imin,imax,10
c	 write(6,*) rb(i)
c	end do
c
c averaging
c
	rb(imin-1)=rb(imin+1)
	rb(imax+1)=rb(imax-1)
	do i=imin,imax
c	  rb(i) = 0.25*( rb(i-1) + 2.*rb(i) + rb(i+1) )
	end do
c
c old time level
c
	do i=imin,imax
	  zb(i) = zov(nb(1,i))
	end do
c
	else if(mode.eq.2) then
c
c copy results to rzv
c
	  do i=imin,imax
	    rzv(nb(0,i)) = rb(i)
	  end do
c
	else if(mode.eq.3) then
c
c Schrimpf radiation condition
c
	  if(niter.le.3) write(6,*) 'SSSSSSSSSSSSSchrimpf'

	  do i=imin,imax
	    rb(i)=0.
	    do ii=1,2
	      ie=iipp(ii,i)
	      rb(i)=rb(i)-vnv(ie)/sqrt(grav*hev(ie))
	    end do
	    rb(i)=0.5*rb(i)
	  end do 
	  do i=imin,imax
	    rzv(nb(0,i)) = rb(i)
	  end do
c
	else if(mode.eq.4.or.mode.eq.5) then
c
	  if(niter.le.3) write(6,*) 'SsSsSsSschrimpffffffff...'

	  do i=imin,imax
	    rb(i)=0.
	    zb(i)=0.
	  end do
	  do ie=1,nel
	    area=ev(10,ie)
	    do ii=1,3
		    x=xgv(nen3v(ii,ie))
		    y=ygv(nen3v(ii,ie))
		    ix=iround(x/3000.)+1
		    iy=iround(y/3000.)
		    if(iy.eq.1) then
			k=nen3v(ii,ie)
	      		rb(ix)=rb(ix)-area*vnv(ie)/sqrt(grav*hev(ie))
			zb(ix)=zb(ix)+area
		    end if
	    end do
	  end do
c
	  do i=imin,imax
	    rb(i) = rb(i)/zb(i)
	  end do
c
c averaging
c
	if(mode.eq.5) then

	  if(niter.le.3) write(6,*) 'averaging...'

	rb(imin-1)=rb(imin+1)
	rb(imax+1)=rb(imax-1)
	rbold=rb(imin-1)
	do i=imin,imax
	  rbaux = rb(i)
	  rb(i) = 0.25*( rbold + 2.*rb(i) + rb(i+1) )
	  rbold = rbaux
	end do

	end if
c
	  do i=imin,imax
	    rzv(nb(0,i)) = rb(i)
	  end do

	else if(mode.eq.6) then

	  if(niter.le.3) write(6,*) 'Schrimpf equilibrated...'

	  do i=imin,imax
	    rb(i)=0.
	    zb(i)=0.
	  end do

	  do ie=1,nel
	    area=ev(10,ie)
	    x=0.
	    iok=0
	    do ii=1,3
		    x=x+xgv(nen3v(ii,ie))
		    y=ygv(nen3v(ii,ie))
		    iy=iround(y/3000.)
		    if(iy.eq.0) iok=1
	    end do
	    if(iok.eq.1) then
	    	x=x/3.
		ix=iround(x/3000.)+1
	  	rb(ix)=rb(ix)-area*vnv(ie)/sqrt(grav*hev(ie))
		zb(ix)=zb(ix)+area
		ix=ix+1
	  	rb(ix)=rb(ix)-area*vnv(ie)/sqrt(grav*hev(ie))
		zb(ix)=zb(ix)+area
	    end if
	  end do

	  do i=imin,imax
	    rzv(nb(0,i)) = rb(i)/zb(i)
	  end do

	else if(mode.eq.7) then

	  if(niter.le.3) write(6,*) 'Schrimpf spatial implicit...'

c all elements with same area

	  do i=imin,imax
	    rb(i)=0.
	    zb(i)=0.
	  end do

	  do ie=1,nel
	    x=0.
	    iok=0
	    area=1.
	    do ii=1,3
		    x=x+xgv(nen3v(ii,ie))
		    y=ygv(nen3v(ii,ie))
		    iy=iround(y/3000.)
		    if(iy.eq.0) iok=1
	    end do
	    if(iok.eq.1) then
	    	x=x/3.
		ix=iround(x/3000.)+1
	  	rb(ix)=rb(ix)-area*vnv(ie)/sqrt(grav*hev(ie))
		zb(ix)=zb(ix)+area
		ix=ix+1
	  	rb(ix)=rb(ix)-area*vnv(ie)/sqrt(grav*hev(ie))
		zb(ix)=zb(ix)+area
	    end if
	  end do

	  do i=imin,imax
	    rb(i-imin+1) = rb(i)
	  end do

	  n=imax-imin+1

c fill in main and upper diagonal

	  dd(locsps(1,1,n,1))=2.
	  dd(locsps(1,1+1,n,1))=1.
	  do i=2,n-1
	    dd(locsps(i,i,n,1))=4.
	    dd(locsps(i,i+1,n,1))=1.
	  end do
	  dd(locsps(n,n,n,1))=2.

	  call mchb(rb,dd,n,1,1,1,1.e-6,ier)
	  if(ier.ne.0) stop 'crador : inverting matrix'

	  do i=imin,imax
	    rzv(nb(0,i)) = rb(i-imin+1)
	  end do

	else if(mode.eq.8) then

	  if(niter.le.3) write(6,*) 'Schrimpf with finite elements...'

	  do i=imin,imax
	    rb(i)=0.
	    zb(i)=0.
	  end do

	  do ie=1,nel
	    area=ev(10,ie)
	    x=0.
	    iok=0
	    do ii=1,3
		    x=xgv(nen3v(ii,ie))
		    y=ygv(nen3v(ii,ie))
		    ix=iround(x/3000.)+1
		    iy=iround(y/3000.)
		    if(iy.eq.0) then
	  		rb(ix)=rb(ix)-area*vnv(ie)/sqrt(grav*hev(ie))
			zb(ix)=zb(ix)+area
		    end if
	    end do
	  end do

	  do i=imin,imax
	    rzv(nb(0,i)) = rb(i)/zb(i)
	  end do

	else if(mode.eq.9) then

	  if(niter.le.3) write(6,*) 'Schrimpf iterated...'

	  do i=imin,imax
	    rb(i)=0.
	    zb(i)=0.
	  end do

	  do ie=1,nel
	    area=ev(10,ie)
	    x=0.
	    iok=0
	    do ii=1,3
		    x=x+xgv(nen3v(ii,ie))
		    y=ygv(nen3v(ii,ie))
		    iy=iround(y/3000.)
		    if(iy.eq.0) iok=1
	    end do
	    if(iok.eq.1) then
	    	x=x/3.
		ix=iround(x/3000.)+1
	  	rb(ix)=rb(ix)-area*vnv(ie)/sqrt(grav*hev(ie))
		zb(ix)=zb(ix)+area
		ix=ix+1
	  	rb(ix)=rb(ix)-area*vnv(ie)/sqrt(grav*hev(ie))
		zb(ix)=zb(ix)+area
	    end if
	  end do

	  res=0.
	  do i=imin,imax
	    rb(i)=rb(i)/zb(i)
	    res=res+abs(rb(i)-rzv(nb(0,i)))
	    rzv(nb(0,i)) = rb(i)
	  end do
	  write(6,*) (rb(i),i=imin,imax,10)

	end if
c
	return
	end
c
c*****************************************************************
c
	subroutine gwi(hia,hik,ibrad,b,c,zold,dt,h,area)
c
c implements radiation condition (GWI)
c
	implicit none
c
c arguments
	real hia(3,3),hik(3)
	integer ibrad(3)
	real b(3),c(3),zold(3)
	real dt,h,area
c local
	logical btest,bltest,bgtest
	integer i,ib,ib1,ib2,ii,ii1,ii2
	integer j
	integer ibtot,itest
	real e(3)
	real right,eps,aux
	real auxib,auxib1,auxib2	
	real auxm(3,3),rm(3),rmaux(3)
	integer iaux(3)
c save
	save bgtest,itest
	data bgtest,itest /.true. , 0 /
c	
c find out how many nodes are on boundary
c
	ib = 0
	ib1 = 0
	ib2 = 0
	ii = 0
	ii1 = 0
	ii2 = 0

	ibtot=0
	do i=1,3
	  if(ibrad(i).eq.1) ibtot=ibtot+1
	end do
c
	if(ibtot.eq.0) then
	  return
	else if(ibtot.eq.1) then
	  if(ibrad(1).eq.1) then
	    ib=1
	    ii1=2
	    ii2=3
	  else if(ibrad(2).eq.1) then
	    ib=2
	    ii1=1
	    ii2=3
	  else
	    ib=3
	    ii1=1
	    ii2=2
	  end if
	else if(ibtot.eq.2) then
	  if(ibrad(1).eq.0) then
	    ii=1
	    ib1=2
	    ib2=3
	  else if(ibrad(2).eq.0) then
	    ii=2
	    ib1=1
	    ib2=3
	  else
	    ii=3
	    ib1=1
	    ib2=2
	  end if
	else
	  write(6,*) 'three boundary nodes not possible for GWI'
	  stop 'error stop : gwi'
	end if
c
c ibtot is total number of nodes on boundary
c
c ib,ib1,ib2 are boundary nodes
c ii,ii1,ii2 are internal nodes
c
c compute constants
c
c in computation of e : 
c		use b for right/left  boundary
c		use c for lower/upper boundary
c		use + for upper/right boundary
c		use - for lower/left  boundary
c
	aux=2.*dt*sqrt(9.81*h)
 	right=0.
	do i=1,3
	  e(i) = +aux*c(i)
	  right = right + e(i)*zold(i)
	end do
	right=-2.*right
c
	do i=1,3
	  do j=1,3
	    auxm(i,j)=hia(i,j)
	  end do
	  rm(i)=hik(i)
	end do	    
c
c decide if write to terminal
c
	if(bgtest) then	!once bgtest is wrong it is always wrong
	  bltest=.false.
	  do i=1,3
	    if(abs(zold(i)).ge.0.1) bltest=.true.
	  end do
c	  bltest=.true.
	  if(bltest) itest=itest+1
	  bgtest=itest.le.200
	end if
	btest=bgtest.and.bltest
c
c start writing
c	
	if(btest) then
	  do i=1,3
	    write(6,*) (hia(i,j),j=1,3),hik(i),ibrad(i)
	  end do
	end if
c
	call matinv(auxm,iaux,rmaux,3,3)
c
	do i=1,3
	  aux=0.
	  do j=1,3
	    aux=aux+auxm(i,j)*rm(j)
	  end do
	  rmaux(i)=aux
	end do
c
	if(btest) then
	  write(6,*)
	  write(6,*) (rmaux(i),i=1,3)
	  write(6,*)
	end if
c
c assimilate into matrix
c
	if(ibtot.eq.1) then
c
c ib      is boundary node
c ii1,ii2 are internal nodes
c
	  do i=1,3
	    hia(ib,i)=area*e(i)
	  end do
	  hia(ib,ib)=hia(ib,ib)+area*1.
	  hik(ib)=area*right
	  auxib=1./hia(ib,ib)
c
	  aux=hia(ii1,ib)*auxib
	  hia(ii1,ib)=0.
	  hia(ii1,ii1)=hia(ii1,ii1)-hia(ib,ii1)*aux
	  hia(ii1,ii2)=hia(ii1,ii2)-hia(ib,ii2)*aux
	  hik(ii1)=hik(ii1)-hik(ib)*aux
c
	  aux=hia(ii2,ib)*auxib
	  hia(ii2,ib)=0.
	  hia(ii2,ii2)=hia(ii2,ii2)-hia(ib,ii2)*aux
	  hia(ii2,ii1)=hia(ii2,ii1)-hia(ib,ii1)*aux
	  hik(ii2)=hik(ii2)-hik(ib)*aux
c
	  eps=(hia(ii2,ii1)-hia(ii1,ii2))/(hia(ii2,ii1)+hia(ii1,ii2))
c
	  aux=hia(ib,ii1)*auxib
	  hia(ii1,ib)=hia(ib,ii1)
	  hia(ii1,ii1)=hia(ii1,ii1)*(1.+eps)+hia(ib,ii1)*aux
	  hia(ii1,ii2)=hia(ii1,ii2)*(1.+eps)+hia(ib,ii2)*aux
	  hik(ii1)=hik(ii1)*(1.+eps)+hik(ib)*aux
c
	  aux=hia(ib,ii2)*auxib
	  hia(ii2,ib)=hia(ib,ii2)
	  hia(ii2,ii2)=hia(ii2,ii2)*(1.-eps)+hia(ib,ii2)*aux
	  hia(ii2,ii1)=hia(ii2,ii1)*(1.-eps)+hia(ib,ii1)*aux
	  hik(ii2)=hik(ii2)*(1.-eps)+hik(ib)*aux
c
	else if(ibtot.eq.2) then
c
c ib1,ib2 are boundary nodes
c ii      is internal node
c
	  do i=1,3
	    hia(ib1,i)=area*e(i)
	    hia(ib2,i)=area*e(i)
	  end do
	  hia(ib1,ib1)=hia(ib1,ib1)+area*1.
	  hia(ib2,ib2)=hia(ib2,ib2)+area*1.
	  hik(ib1)=area*right
	  hik(ib2)=area*right
	  auxib1=1./hia(ib1,ib1)
	  auxib2=1./hia(ib2,ib2)
c
	  aux=hia(ib1,ib2)*auxib2
	  hia(ib1,ib1)=hia(ib1,ib1)-hia(ib2,ib1)*aux
	  hia(ib1,ib2)=0.
	  hia(ib1,ii)=hia(ib1,ii)*(1.-aux)
	  hik(ib1)=hik(ib1)*(1.-aux)
c
	  aux=hia(ib2,ib1)*auxib1
	  hia(ib2,ib2)=hia(ib2,ib2)-hia(ib1,ib2)*aux
	  hia(ib2,ib1)=0.
	  hia(ib2,ii)=hia(ib2,ii)*(1.-aux)
	  hik(ib2)=hik(ib2)*(1.-aux)
c
	  auxib1=1./hia(ib1,ib1)
	  auxib2=1./hia(ib2,ib2)
c
	  hia(ii,ii)=hia(ii,ii)-hia(ii,ib1)*auxib1*hia(ib1,ii)
	  hia(ii,ii)=hia(ii,ii)-hia(ii,ib2)*auxib2*hia(ib2,ii)
	  hik(ii)=hik(ii)-hia(ii,ib1)*auxib1*hik(ib1)
	  hik(ii)=hik(ii)-hia(ii,ib2)*auxib2*hik(ib2)
c
	  hia(ii,ii)=hia(ii,ii)+auxib1*hia(ib1,ii)**2
	  hia(ii,ii)=hia(ii,ii)+auxib2*hia(ib2,ii)**2
	  hik(ii)=hik(ii)-auxib1*hia(ib1,ii)*hik(ib1)
	  hik(ii)=hik(ii)-auxib2*hia(ib2,ii)*hik(ib2)
	  hia(ii,ib1)=hia(ib1,ii)
	  hia(ii,ib2)=hia(ib2,ii)
c
	end if
c
	do i=1,3
	  do j=1,3
	    auxm(i,j)=hia(i,j)
	  end do
	  rm(i)=hik(i)
	end do	    
c
	if(btest) then
	  do i=1,3
	    write(6,*) (hia(i,j),j=1,3),hik(i),ibrad(i)
	  end do
	end if
c
	call matinv(auxm,iaux,rmaux,3,3)
c
	do i=1,3
	  aux=0.
	  do j=1,3
	    aux=aux+auxm(i,j)*rm(j)
	  end do
	  rmaux(i)=aux
	end do
c
	if(btest) then
	  write(6,*)
	  write(6,*) (rmaux(i),i=1,3)
	  write(6,*)
	end if
c
	return
	end
c
c*************************************************************
c
