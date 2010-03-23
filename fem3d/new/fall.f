
c********************************************************************

	subroutine fall_vel(gd,fall)

	implicit none

	real gd,fall

	real rhos,rhow,visk

	rhos = 2500.
	rhow = 1025
	visk = 1.3e-3 / rhow

	call comp_fall(rhos,rhow,gd,visk,fall)

	end

c********************************************************************

	subroutine comp_fall(rhos,rhow,gd,visk,fall)

	implicit none

	real rhos,rhow,gd,visk,fall

	real g,drho,aux,dx
	
      G = 9.81
      DRHO = RHOS-RHOW
      AUX = (DRHO*G*GD)/RHOW

C COMPUTE THE DIMENSIONLESS PARTICLE DIAMETER
      DX = GD * (((DRHO/RHOW)*G)/VISK**2.)**(1./3.)

C COMPUTE FALL VELOCITY (Soulsby, 1997)
      FALL = (VISK/GD)*((10.36**2. + 1.049*DX**3.)**0.5 - 10.36)

	end

c********************************************************************

	subroutine fall_test

	implicit none

	integer i
	real gd,fact,fall

	gd = 1000.	!in micron
	fact = 1.

	gd = gd * 1.e-6

	do i=1,10
	  gd = gd / fact
	  call fall_vel(gd,fall)
	  write(6,*) i,1.e+6*gd,fall
	  fact = fact*1.2
	end do

	end

c********************************************************************

	program main
	call fall_test
	end

c********************************************************************

