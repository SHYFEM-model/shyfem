c
c $Id: subzzz.f,v 1.6 2003/09/03 14:31:39 georg Exp $
c
c mercator projection routines
c
c contents :
c
c subroutine mercin(phi,omega,scale)
c		initializes common block for mercator transformation
c subroutine g2c(n,alon,alat,x,y)
c		transforms coordinates from geographical to cartesian
c subroutine c2g(n,alon,alat,x,y)
c		transforms coordinates from cartesian to geographical
c
c subroutine xy2geo(alon,alat,x,y)
c subroutine geo2xy(alon,alat,x,y)
c
c function phidis(phi)
c		computes distance between phi and phi0 in km
c
c subroutine uvcomp(...)	transformation coord. cartesian <--> polar
c subroutine mercoo(...)	transformation coord. geographical <--> spatial
c function fifi(fi)		compute distance between fi and fi0
c subroutine zgrid0o(...)	reset matrices and vectors
c
c revision log :
c
c 19.11.1998	ggu	completly re-written -> user friendly
c 24.05.1999	ggu	new routines xy2geo and geo2xy
c 03.09.2003	ggu	changed for compiler warnings on origin
c
c**********************************************************

	subroutine mercin(phi,omega,scale)

c initializes common block for mercator transformation
c
c must be called befor using the other routines
c
c phi		latitude of the projection in degrees
c omega		longitude of reference in degrees
c scale		scale of the map used to digitize
c		e.g., 10000. for scale 1:10000 and units are in cm. 
c		To transform to/from meters use scale = 100.
c
c phi0          parallelo di proiezione in gradi centesimali
c omega0        meridiano di riferimento in gradi centesimali
c scale0        scala della carta geografica digitalizzata
c               ...(10000. per 1:10000, se x,y e' gia in metri
c               ...o se metri sono desiderati : sca=100.)
c
c e             eccentricita' della terra secondo Hayford (1950)
c mu            modulo dei logaritmi di Briggs
c rea           raggio della terra sferica all'equatore in km
c nu            raggio della terra ellissoidica in km

	implicit none

	real phi,omega,scale

	real e,mu,rea,pi,rad,co,n0,ta0log,phi0,omega0,scale0
	common /mercom/ e,mu,rea,pi,rad,co,n0,ta0log,phi0,omega0,scale0
	save /mercom/

        e = 0.0819918904
        mu = 0.4342944819
        rea = 6370.0

        pi = 4.0 * atan( 1. )
        rad = pi/180.

        co = cos(phi*rad)
        n0 = rea/(1.-e*e*sin(phi*rad)*sin(phi*rad))

        ta0log = tan((45.+phi/2.)*rad)
	ta0log = alog10(ta0log)

        phi0 = phi
        omega0 = omega
        scale0 = scale * 1.e-5	! km -> cm

	end
	
c**********************************************************

	subroutine g2c(n,alon,alat,x,y)

c transforms coordinates from geographical to cartesian

	implicit none

	integer n		!length of array to transform
	real alon(1),alat(1)	!geographical coordinates
	real x(1),y(1)		!cartesian coordinates

	real e,mu,rea,pi,rad,co,n0,ta0log,phi0,omega0,scale0
	common /mercom/ e,mu,rea,pi,rad,co,n0,ta0log,phi0,omega0,scale0

	integer i
	real xx,yy
	real delta,rsca,aux
	real phidis

	rsca = 1. / scale0
	aux = n0 * co * rad

        do i=1,n
           delta = alon(i) - omega0
           xx = aux * delta
           yy = phidis( alat(i) )
           x(i) = xx * rsca
           y(i) = yy * rsca
        end do

	end

c**********************************************************

	subroutine c2g(n,alon,alat,x,y)

c transforms coordinates from cartesian to geographical

	implicit none

	integer n		!length of array to transform
	real alon(1),alat(1)	!geographical coordinates
	real x(1),y(1)		!cartesian coordinates

	real e,mu,rea,pi,rad,co,n0,ta0log,phi0,omega0,scale0
	common /mercom/ e,mu,rea,pi,rad,co,n0,ta0log,phi0,omega0,scale0

	integer i
	real xx,yy,c,f1
	real rrad,aux1,aux2

	aux1 = 1. / ( n0 * co * rad )
	aux2 = mu / ( n0 * co )
	rrad = 1. / rad

        do i=1,n
           xx = x(i)*scale0
           yy = y(i)*scale0
           alon(i) = omega0 + xx * aux1
           c = yy * aux2
           f1 = 10**(c+ta0log)
           f1 = atan(f1)
           f1 = (f1*rrad-45.)*2
           alat(i) = f1
        end do

	end

c**********************************************************

	subroutine xy2geo(alon,alat,x,y)

	implicit none

	real alon,alat,x,y
	real alonv(1),alatv(1),xv(1),yv(1)

	alonv(1) = alon
	alatv(1) = alat
	xv(1)    = x
	yv(1)    = y

	call c2g(1,alonv,alatv,xv,yv)

	alon = alonv(1)
	alat = alatv(1)
	x    = xv(1)
	y    = yv(1)

	end
	
c**********************************************************

	subroutine geo2xy(alon,alat,x,y)

	implicit none

	real alon,alat,x,y
	real alonv(1),alatv(1),xv(1),yv(1)

	alonv(1) = alon
	alatv(1) = alat
	xv(1)    = x
	yv(1)    = y

	call g2c(1,alonv,alatv,xv,yv)

	alon = alonv(1)
	alat = alatv(1)
	x    = xv(1)
	y    = yv(1)

	end
	
c**********************************************************

	function phidis(phi)

c computes distance between phi and phi0 in km
c
c phi		latitude
c phidis	distance in km

	real phi
	real phidis

	real e,mu,rea,pi,rad,co,n0,ta0log,phi0,omega0,scale0
	common /mercom/ e,mu,rea,pi,rad,co,n0,ta0log,phi0,omega0,scale0

	real ta,si,si0,x1,x2

	ta = tan((45.+phi/2)*rad)
	si = sin(phi*rad)
	si0 = sin(phi0*rad)
	x1 = alog10(ta)-ta0log
	x2 = e*e*(si-si0)+(e**4)*(si**3-si0**3)/3.

	phidis = n0*co*(x1/mu-x2)

	end

c**********************************************************
c**********************************************************
c********************************************* old programs
c**********************************************************
c**********************************************************

	subroutine uvcomp(in,u,v,dir,sp,n)
c
c trasformazione coordinate cartesiane in coordinate polari
c ...e viceversa
c
c  STEFANO 4 NOVEMBRE 1985  VERSIONE 1
c  GEORG   4 FEBBRAIO 1989
c
c  in		modo di operazione
c		0: calcola le componenti di un vettore date 
c		...direzione e velocita'
c		1:  calcola direzione e velocita' date le componenti 
c		...u (asse x) e v (asse y)
c  u,v		coordinate cartesiane in direzione x,y
c  dir,sp	coordinate polari (direzione e modulo)
c  n		numero delle coordinate nei vettori
c
	dimension u(1),v(1),dir(1),sp(1)
c
	rad=90./asin(1.)
c
	if(in.eq.0) then	!polari --> cartesiane
	   do i=1,n
	   u(i)=dir(i)		!u utilizzato come vettore ausiliare
	   if(u(i).le.90.and.u(i).ge.0)   u(i)=90.-u(i)
	   if(u(i).le.360.and.u(i).ge.90) u(i)=450.-u(i)
	   u(i)=u(i)/rad
	   u(i)=sp(i)*cos(u(i))
	   v(i)=sp(i)*sin(u(i))
	   end do
	else			!cartesiane --> polari
	   do i=1,n
	   if(u(i).eq.0.) u(i)=0.001
	   sp(i)=sqrt(u(i)**2+v(i)**2)
	   alfa=atan(v(i)/u(i))*rad
	   if(u(i).gt.0.and.v(i).gt.0.or.u(i).gt.0.and.v(i).lt.0)
     +		alfa=90.-alfa
	   if(u(i).lt.0.and.v(i).gt.0.or.u(i).lt.0.and.v(i).lt.0)
     +		alfa=270.-alfa
	   dir(i)=alfa
	   end do
	end if
c
	return
	end
c
c**********************************************************
c
	subroutine mercoo(id,alon,alat,x,y,n,fih0,omeh0,scal)
c
c  per trasformare coordinate geografiche in coordinate spaziali e
c  viceversa secondo la proiezione del Mercatore secante
c
c  STEFANO  giugno 1985
c  GEORG    gennaio 1989
c
c  id		modo di utilizzo della routine
c		...0 : trasformazione di coordinate da 
c		...      geografiche a spaziali in cm
c		...1 : viceversa
c  alon,alat	vettori delle coordinate geografiche in gradi centesimali
c  x,y		vettori delle coordinate spaziali in cm
c  n		numero di dati in vettore
c
c  fih0		parallelo di proiezione in gradi centesimali
c  omeh0	meridiano di riferimento in gradi centesimali
c  scal		scala della carta geografica digitalizzata
c		...(10000. per 1:10000, se x,y e' gia in metri
c		...o se metri sono desiderati : sca=100.)
c
c  e		eccentricita' della terra secondo Hayford (1950)
c  mu		modulo dei logaritmi di Briggs
c  rea		raggio della terra sferica all'equatore in km
c  nu		raggio della terra ellissoidica in km
c
c  ATTENZIONE: le coordinate geografiche devono essere espresse in gradi
c  centesimali
c
c  gli assi di riferimento dai quali x e y vengono calcolate sono il
c  il parallelo di proiezione (fi0) e il meridiano di riferimento (ome0)
c
	real n0,mu
	dimension alon(1),alat(1),x(1),y(1)
	common /mcoo/ fi0,ome0,n0,sca
c
	fi0=fih0
	ome0=omeh0
	sca=scal
c
	e=0.0819918904
	rea=6370.0
	pi=2.0*asin(1.)
	rad=pi/180.
	mu=0.4342944819
c
	co=cos(fi0*rad)
	n0=rea/(1-e*e*sin(fi0*rad)*sin(fi0*rad))
c
	if(id.eq.0) then	!geographic --> spacial
	     do i=1,n
		delta=alon(i)-ome0
		x(i)=n0*co*delta*rad
     		y(i)=fifi(alat(i))
		x(i)=x(i)/(sca*1.e-5)	! trasforma da km in cm
		y(i)=y(i)/(sca*1.e-5)	! trasforma da km in cm
	     end do
	else			!spacial --> geographic
	     do i=1,n
		xx=x(i)*sca*1.e-5	! trasforma da cm in km
		alon(i)=xx/(n0*co*rad)+ome0
  		yy=y(i)*sca*1.e-5	! trasforma da cm in km
		ta0=tan((45.+fi0/2)*rad)
		c=yy*mu/(n0*co)
		f1=10**(c+alog10(ta0))
		f1=atan(f1)
		f1=(f1/rad-45.)*2
     		alat(i)=f1
	     end do
	end if
c
	return
	end
c
c**********************************************************************
c
	function fifi(fi)
c
c  per calcolare la distanza in km tra fi e fi0
c
c  fi		lattitudine del punto
c  fifi		distanza fra fi e fi0
c
	real n0,mu
	common /mcoo/ fi0,ome0,n0,sca
c
	rad=asin(1.)/90.
	e=.0819918904
	mu=.4342944819
	co=cos(fi0*rad)
c
	ta=tan((45.+fi/2)*rad)
	ta0=tan((45.+fi0/2)*rad)
	si=sin(fi*rad)
	si0=sin(fi0*rad)
	x1=alog10(ta)-alog10(ta0)
	x2=e*e*(si-si0)+(e**4)*(si**3-si0**3)/3
	fifi=n0*co*(x1/mu-x2)
c
	return
	end
c
c**************************************************************
c
	subroutine zgrid0o(nlat,nlon,nll,imnew,pout,knxt,zpij)
c
c reset matrices and vectors
c
c old version, use version in library of stefano --> zgrid0
c
	dimension pout(nlon,nlat)
	dimension imnew(nlat),knxt(nll),zpij(nll)
c
	do j=1,nlat
	  imnew(j)=0
	  do i=1,nlon
	    pout(i,j)=1.e35
	  end do
	end do
c
	do i=1,nll
	  knxt(i)=0
	  zpij(i)=0.
	end do
c
	return
	end
