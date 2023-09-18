!-------------------------------------------------------------------
        module heat_util2
!-------------------------------------------------------------------
        contains
!-------------------------------------------------------------------
!*******************************************************************

	subroutine satur(t,p,ew,rw,qw)

! computes saturation values for water vapor
!
! e/p = r / ( eps + r )	=>  r = eps / ( p/e - 1 )
! r   = q / ( 1 - q )	=>  q = r / ( r + 1 )

	implicit none

	double precision t		!temperature [C]			- in
	double precision p		!pressure [mb]				- in
	double precision ew		!saturation vapor pressure [mb]		- out
	double precision rw		!saturation mixing ratio [0-1]		- out
	double precision qw		!saturation specific humidity [0-1]	- out

	double precision eps				!molecular weight ratio
	parameter(eps=0.62197)
	double precision fice
!	parameter (fice=0.00422)		!over ice
	parameter (fice=0.0)
	double precision fsalt
!	parameter (fsalt=0.98)			!over salt
	parameter (fsalt=1.0)

	double precision aux,fw

	aux = (0.7859+0.03477*t)/(1.+0.00412*t)
	aux = aux + fice * t
	fw = 1. + 1.e-6 * p * (4.5+0.0006*t*t)

	ew = fsalt * fw * 10.**aux

	rw = eps / ( p/ew - 1. )
	qw = rw / ( rw + 1. )

	end

!*******************************************************************

	subroutine vapor(t,p,u,e,r,q)

! computes values for water vapor given relative humidity
!
! e/p = r / ( eps + r )	=>  r = eps / ( p/e - 1 )
! r   = q / ( 1 - q )	=>  q = r / ( r + 1 )
! u = r / rw

	implicit none

	double precision t		!temperature [C]			- in
	double precision p		!pressure [mb]				- in
	double precision u		!relative humidity [%] ([0-100])	- in
	double precision e		!vapor pressure [mb]			- out
	double precision r		!mixing ratio [0-1]			- out
	double precision q		!specific humidity [0-1]		- out

	double precision eps				!molecular weight ratio
	parameter(eps=0.62197)

	double precision ew,rw,qw

	call satur(t,p,ew,rw,qw)	!compute saturation values

	r = 0.01 * u * rw
	q = r / ( r + 1. )
	e = p * r / ( eps + r )

	end

!*******************************************************************

	subroutine smithb(w,cd)

	double precision w		!wind speed [m/s]	- in
	double precision cd		!drag coefficient [-]	- out

	if( w .le. 6. ) then
	  cd = 1.1e-3
	else
	  cd = 0.61e-3 + 0.063e-3 * w
	end if

	end

!*******************************************************************

	subroutine charnock(w,cd)

	double precision w		!wind speed [m/s]	- in
	double precision cd		!drag coefficient [-]	- out

	double precision rho	!density of air [kg/m**3]
	double precision g		!gravity [m/s**2]
	double precision z		!anemometer height [m]
	double precision a		!Charnock constant [-]
	double precision k		!von Karman constant [-]

	parameter(rho=1.25,g=9.81,z=10.,a=0.0185,k=0.4)

	double precision cdold,aux,tau

	cd = 1.5e-3
	cdold = 0.
	
	do while( abs(cd-cdold) .gt. 1.e-5 )
	  cdold = cd
	  tau = cd * rho * w * w
	  aux = rho * g * z / ( a * tau )
	  cd = ( k / log(aux) )**2
	end do

	end

!*******************************************************************

	subroutine chartest

        double precision w,cd1,cd2

	do i=1,20
	  w = i
	  call charnock(w,cd1)
	  call smithb(w,cd2)
	  write(6,*) w,cd1,cd2
	end do

	end

!*******************************************************************
!*******************************************************************
!*******************************************************************

!	call chartest
!	end

!*******************************************************************

!-------------------------------------------------------------------
        end module heat_util2
!-------------------------------------------------------------------
