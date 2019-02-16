
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

c routines for k-epsilon model
c
c revision log :
c
c 25.02.1999	ggu	finished writing routines
c 10.04.2008	ggu&ccf	new gotm routine call integrated
c 01.06.2011	ggu&cpb	call to do_gotm_turb was wrong (BUG)
c 10.02.2015	ggu	in gotm call to do_gotm_turb() with ddt and not dt
c
c notes :
c
c u,v,rho are evaluated at the center of the layer
c other variables (k,eps,rkm,rkt) at the top and bottom of the layers
c
c********************************************************************

	subroutine gotm(lmax,dt,rho0,taus,taub,dl,u,v,rho
     +				,rkm,rkt,k,eps,len)

c computes turbulent diffusion coefficients through gotm

	implicit none

c arguments
	integer lmax		!number of levels [1...lmax]
	real dt			!time step
	real rho0		!reference density
	real taus		!surface stress at new time level
	real taub		!bottom stress at new time level
	real dl(1:lmax)		!bottom depth of single levels
	real u(1:lmax)		!velocity of single levels in x direction
	real v(1:lmax)		!velocity of single levels in y direction
	real rho(1:lmax)	!density of single levels
	real rkm(0:lmax)	!eddy coefficient for momentum		(out)
	real rkt(0:lmax)	!eddy coefficient for scalar quantities (out)
	real k(0:lmax)		!turbulent kinetic energy		(out)
	real eps(0:lmax)	!turbulent dissipation			(out)
	real len(0:lmax)	!length scale				(out)

	integer nlevdi
	parameter(nlevdi=200)
	real g,karman
	parameter(g=9.81,karman=0.4)

	real*8 dzk(nlevdi)		!z difference between k (eps) levels
	real*8 dzr(0:nlevdi)		!z difference between rho (u,v) levels
	real*8 dzl(0:nlevdi)		!needed for gotm
	real*8 b(nlevdi)			!buoyancy
	real*8 n2(0:nlevdi)		!Brunt-V. frequency
	real*8 m2(0:nlevdi)		!shear production

	real*8 hh(0:nlevdi)
	real*8 nn(0:nlevdi)
	real*8 ss(0:nlevdi)
	real*8 num(0:nlevdi)
	real*8 nuh(0:nlevdi)
	real*8 rkin(0:nlevdi)
	real*8 reps(0:nlevdi)
	real*8 rlen(0:nlevdi)

	integer l,lmax1,lp
	real*8 ddt,u_taus,u_taub,z0,ud,vd,aux
	real*8 z0s,z0b,depth

	if( lmax .gt. nlevdi ) stop 'error stop: nlevdi'

	lmax1 = lmax - 1

	ddt = dt
        u_taus = sqrt( taus )	!u_star
        u_taub = sqrt( taub )
	depth = dl(lmax)
	z0s = 0.03
	z0b = 0.03

c compute level differences
c
c z0 is height of surface -> if highly variable -> pass into routine
c
c dzk == dz        dzr == dZ  (Langland and Liou, MWR, May 1996, 905-918)

	z0 = 0.
	dzk(1) = dl(1) - z0
	dzr(0) = 0.
	do l=2,lmax1
	  dzk(l) = dl(l) - dl(l-1)
	  dzr(l-1) = 0.5 * ( dzk(l) + dzk(l-1) )
	end do
	dzk(l) = dl(l) - dl(l-1)	!last level -> l = lmax
	dzr(l-1) = 0.5 * ( dzk(l) + dzk(l-1) )
	dzr(l) = 0.

	do l=1,lmax1
	  dzl(l) = dzk(l)
	  dzk(l) = 1. / dzk(l)
	  dzr(l) = 1. / dzr(l)
	end do
	dzl(l) = dzk(l)
	dzk(l) = 1. / dzk(l)		!last level -> l = lmax

c compute b(uoyancy) and Brunt-Vaisala frequency n2 (N^2)

	aux = -g / rho0
	do l = 1,lmax
	  b(l) = aux * ( rho(l) - rho0 )
	end do

        nn(1) = 0.
        ss(1) = 0.
	do l = 1,lmax1
	  n2(l) = ( b(l) - b(l+1) ) * dzr(l)
	  n2(l) = max(0.D+0,n2(l))
	end do

c compute horizontal shear m2 (M^2)

	do l = 1,lmax1
	  ud = ( u(l) - u(l+1) ) * dzr(l)
	  vd = ( v(l) - v(l+1) ) * dzr(l)
	  m2(l) = ud*ud + vd*vd
	end do

c prepare data for gotm call

	do l=1,lmax
	  lp = lmax - l + 1
	  hh(lp) = dzl(l)
	end do

	do l=1,lmax1
	  lp = lmax-l
	  nn(lp) = n2(l)
	  ss(lp) = m2(l)
	end do
	nn(0) = nn(1)
	ss(0) = ss(1)
	nn(lmax) = nn(lmax1)
	ss(lmax) = ss(lmax1)

	do l=0,lmax
	  lp = lmax-l
	  num(lp) = rkm(l)
	  nuh(lp) = rkt(l)
	  rkin(lp) = k(l)
	  reps(lp) = eps(l)
	  rlen(lp) = len(l)
	end do

c	call gotmturb(lmax,ddt,hh,NN,SS,num,nuh
c     +			,rkin,reps,rlen,u_taus,u_taub)

	call do_gotm_turb(lmax,ddt,depth,u_taus,u_taub,z0s,z0b
     +			,hh,NN,SS,num,nuh
     +			,rkin,reps,rlen)

	do l=0,lmax
	  lp = lmax-l
	  rkm(lp) = num(l)
	  rkt(lp) = nuh(l)
	  k(lp)   = rkin(l)
	  eps(lp) = reps(l)
	  len(lp) = rlen(l)
	end do

	end

c********************************************************************

	subroutine ma(lmax,dt,rho0,taus,taub,dl,u,v,rho,rkm,rkt)

c computes turbulent diffusion coefficients with Munk-Anderson theory

	implicit none

c arguments
	integer lmax		!number of levels [1...lmax]
	real dt			!time step
	real rho0		!reference density
	real taus		!surface stress at new time level
	real taub		!bottom stress at new time level
	real dl(1:lmax)		!bottom depth of single levels
	real u(1:lmax)		!velocity of single levels in x direction
	real v(1:lmax)		!velocity of single levels in y direction
	real rho(1:lmax)	!density of single levels
	real rkm(0:lmax)	!eddy coefficient for momentum		(out)
	real rkt(0:lmax)	!eddy coefficient for scalar quantities (out)

	integer nlevdi
	parameter(nlevdi=200)
	real a0,k0
	parameter(a0=10.e-4,k0=5.e-4)
	real alpha,beta
	parameter(alpha=-0.5,beta=-1.5)
	real aa,bb
	parameter(aa=10,bb=3.33)
	real g,karman
	parameter(g=9.81,karman=0.4)

	real b(nlevdi)			!buoyancy
	real n2(0:nlevdi)		!Brunt-V. frequency
	real m2(0:nlevdi)		!shear production
	real dzk(nlevdi)		!z difference between k (eps) levels
	real dzr(0:nlevdi)		!z difference between rho (u,v) levels

	integer l,lmax1
	real aux
	real z0,ud,vd,ri

	if( lmax .gt. nlevdi ) stop 'error stop: nlevdi'

	lmax1 = lmax - 1

c compute level differences
c
c z0 is height of surface -> if highly variable -> pass into routine
c
c dzk == dz        dzr == dZ  (Langland and Liou, MWR, May 1996, 905-918)

	z0 = 0.
	dzk(1) = dl(1) - z0
	dzr(0) = 0.
	do l=2,lmax1
	  dzk(l) = dl(l) - dl(l-1)
	  dzr(l-1) = 0.5 * ( dzk(l) + dzk(l-1) )
	end do
	dzk(l) = dl(l) - dl(l-1)	!last level -> l = lmax
	dzr(l-1) = 0.5 * ( dzk(l) + dzk(l-1) )
	dzr(l) = 0.

	do l=1,lmax1
	  dzk(l) = 1. / dzk(l)
	  dzr(l) = 1. / dzr(l)
	end do
	dzk(l) = 1. / dzk(l)		!last level -> l = lmax

c compute b(uoyancy) and Brunt-Vaisala frequency n2 (N^2)

	aux = -g / rho0
	do l = 1,lmax
	  b(l) = aux * ( rho(l) - rho0 )
	end do

	do l = 1,lmax1
	  n2(l) = ( b(l) - b(l+1) ) * dzr(l)
	  n2(l) = max(0.,n2(l))
	end do

c compute horizontal shear m2 (M^2)

	do l = 1,lmax1
	  ud = ( u(l) - u(l+1) ) * dzr(l)
	  vd = ( v(l) - v(l+1) ) * dzr(l)
	  m2(l) = ud*ud + vd*vd
	end do

c compute new turbulent viscosities and diffusivities

	do l=1,lmax1
	  ri = 0.
	  if( m2(l) .gt. 0 ) ri = n2(l)/m2(l)
	  rkm(l) = a0*(1.+aa*ri)**alpha  
	  rkt(l) = k0*(1.+bb*ri)**beta  
	  !write(6,*) '*** ',l,n2(l),m2(l),ri,rkm(l),rkt(l)
	end do

c boundary conditions

	rkm(0) = rkm(1)
	rkt(0) = rkt(1)
	rkm(lmax) = rkm(lmax1)
	rkt(lmax) = rkt(lmax1)

	end

c********************************************************************

	subroutine keps(lmax,dt,rho0,taus,taub,dh,u,v,rho,rkm,rkt,k,eps)

c computes turbulent diffusion coefficients with k-epsilon theorie

	implicit none

c arguments
	integer lmax		!number of levels [1...lmax]
	real dt			!time step
	real rho0		!reference density
	real taus		!surface stress at new time level
	real taub		!bottom stress at new time level
	real dh(1:lmax)		!thickness of single levels
	!real dl(1:lmax)	!bottom depth of single levels
	real u(1:lmax)		!velocity of single levels in x direction
	real v(1:lmax)		!velocity of single levels in y direction
	real rho(1:lmax)	!density of single levels
	real rkm(0:lmax)	!eddy coefficient for momentum		(out)
	real rkt(0:lmax)	!eddy coefficient for scalar quantities (out)
	real k(0:lmax)		!turbulent kinetic energy		(out)
	real eps(0:lmax)	!turbulent dissipation			(out)
c parameters
c	include 'param.h'
	integer nlevdi
	parameter(nlevdi=200)

	real g,karman
	real avumol,avtmol,avsmol
	real kmin,epsmin
	real cmu,sqrtcmu
	real sigmak,sigmae
        real ceps1,ceps2
        real ceps3p,ceps3n,ceps3f
	real cmu0,cmu03,cmu06
	real amin,amax
	real a0,a1,a2,a3,b0,b1,b2
	parameter(g=9.81,karman=0.4)
	parameter(avumol=1.3e-6,avtmol=1.4e-7,avsmol=1.1e-9)
	parameter(kmin=1.e-10,epsmin=1.e-12)
c	parameter(kmin=3.e-6,epsmin=5.e-10)
c---------------------------------------------------------------- Patrick
c	parameter(cmu=0.091,sqrtcmu=0.3)
c	parameter(sigmak=1.0,sigmae=1.3)
c	parameter(ceps1=1.51,ceps2=1.92)
c	parameter(ceps3p=0.2,ceps3n=1.0,ceps3f=ceps1)
c	parameter(cmu0=1.0,cmu03=cmu0**3,cmu06=cmu03**2)
c	parameter(amin=-1.37,amax=20.4)
c	parameter(a0=0.023,a1=1.0,a2=0.714,a3=0.067)
c	parameter(b0=0.125,b1=1.0,b2=0.603)
c---------------------------------------------------------------- Hans
	parameter(cmu=0.091,sqrtcmu=0.3)
	parameter(sigmak=1.0,sigmae=1.08)
	parameter(ceps1=1.44,ceps2=1.92)
	parameter(ceps3p=-0.4,ceps3n=1.0,ceps3f=1.0)
	parameter(cmu0=0.5562,cmu03=cmu0**3,cmu06=cmu03**2)
	parameter(amin=-0.0466,amax=0.56)
	parameter(a0=2.182,a1=1.0,a2=20.40,a3=53.12)
	parameter(b0=0.6985,b1=1.0,b2=17.34)
c----------------------------------------------------------------
c local variables
	logical bdebug
	integer l,lmax1,lm1
	real aux,ud,vd
	real ps,gg,forc
	real a,su,sb,k2,reps
	real ceps3

	real b(nlevdi)			!buoyancy
	real n2(0:nlevdi)		!Brunt-V. frequency
	real m2(0:nlevdi)		!shear production
	real dzk(nlevdi)		!z difference between k (eps) levels
	real dzr(0:nlevdi)		!z difference between rho (u,v) levels
	real rke(nlevdi)		!eddy coefficient for epsilon
	real rkk(nlevdi)		!eddy coefficient for kinetic energy

	real*8 aa(0:nlevdi)		!arrays for tri-diagonal algorithm
	real*8 bb(0:nlevdi)		!...
	real*8 cc(0:nlevdi)		!...
	real*8 rr(0:nlevdi)		!...
	real*8 gam(0:nlevdi)		!...
	real*8 gam1(0:nlevdi)		!...

c u,v,rho are evaluated at the center of the layer
c other variables (k,eps,rkm,rkt) at the top and bottom of the layers
c
c dh(1) contains already surface level variation
c dzk == dz        dzr == dZ  (Langland and Liou, MWR, May 1996, 905-918)

c---------------------------------------------------------------------
c start of code
c---------------------------------------------------------------------

	if( lmax .gt. nlevdi ) stop 'error stop: nlevdi'

	bdebug = .true.
	bdebug = .false.
	lmax1 = lmax - 1
	lm1 = lmax + 1

c---------------------------------------------------------------------
c compute level differences
c---------------------------------------------------------------------

	dzk(1) = dh(1)
	dzr(0) = 0.
	do l=2,lmax
	  dzk(l) = dh(l)
	  dzr(l-1) = 0.5 * ( dzk(l) + dzk(l-1) )
	end do
	dzr(l) = 0.
	if( bdebug ) write(6,*) 'keps 12'

	do l=1,lmax1
	  dzk(l) = 1. / dzk(l)
	  dzr(l) = 1. / dzr(l)
	end do
	dzk(l) = 1. / dzk(l)		!last level -> l = lmax

c---------------------------------------------------------------------
c set eddy coefficient for k and eps diffusion
c rkm(0) and rkm(lmax) may not be set -> avoid using them
c---------------------------------------------------------------------

	rkk(1) = rkm(1) / sigmak
	rke(1) = rkm(1) / sigmae
	do l=2,lmax1
	  aux = rkm(l-1) + rkm(l)
	  rkk(l) = ( 0.5 / sigmak ) * aux
	  rke(l) = ( 0.5 / sigmae ) * aux
	end do
	rkk(lmax) = rkm(lmax1) / sigmak
	rke(lmax) = rkm(lmax1) / sigmae

c---------------------------------------------------------------------
c compute b(uoyancy) and Brunt-Vaisala frequency n2 (N^2)
c---------------------------------------------------------------------

	aux = -g / rho0
	do l = 1,lmax
	  b(l) = aux * ( rho(l) - rho0 )
	end do

	do l = 1,lmax1
	  n2(l) = ( b(l) - b(l+1) ) * dzr(l)
	  n2(l) = max(0.,n2(l))
	end do

c---------------------------------------------------------------------
c compute horizontal shear m2 (M^2)
c---------------------------------------------------------------------

	do l = 1,lmax1
	  ud = ( u(l) - u(l+1) ) * dzr(l)
	  vd = ( v(l) - v(l+1) ) * dzr(l)
	  m2(l) = ud*ud + vd*vd
	end do

c---------------------------------------------------------------------
c turbulent dissipation (eps)
c---------------------------------------------------------------------

	do l = 1,lmax-1

	  ps = rkm(l) * m2(l)
	  gg = - rkt(l) * n2(l)

	  if( n2(l) .gt. 0. ) then
	    ceps3 = ceps3f * ceps3p
	  else
	    ceps3 = ceps3f * ceps3n
	  end if

	  forc = ( dt / k(l) ) * ( 
     +		ceps1 * ps + ceps3 * gg - ceps2 * eps(l)
     +				 )
	  aux = dt*dzr(l)
	  aa(l) = - aux * (rke(l)*dzk(l))
	  bb(l) = 1. + aux * ( rke(l)*dzk(l) + rke(l+1)*dzk(l+1) )
	  cc(l) = - aux * (rke(l+1)*dzk(l+1))
	  rr(l) = eps(l)
	  if( forc .lt. 0.0 ) then
	    bb(l) = bb(l) - forc
	  else
	    rr(l) = rr(l) + forc * eps(l)
	  end if

	end do
	if( bdebug ) write(6,*) 'keps 17'

c	---------------------------------
c	impose boundary conditions
c	---------------------------------

	aa(0) = 0.
	bb(0) = 1.
	cc(0) = 0.
	rr(0) = (taus/rho(1))**1.5 * 0.00002 / karman   !ERROR corrected
	rr(0) = (taus/rho(1))**1.5 / (dzk(1)*karman)   !ERROR corrected
	aa(lmax) = 0.
	bb(lmax) = 1.
	cc(lmax) = 0.
	rr(lmax) = (taub/rho(lmax))**1.5 / (dzk(lmax)*karman)

c	---------------------------------
c	solve system
c	---------------------------------

c	call tridag(aa(1),bb(1),cc(1),rr(1),eps(1),gam(1),lmax-1)
	call tridagd(aa(0),bb(0),cc(0),rr(0),eps(0),gam1(0),gam(0),lm1)

	!eps(0) = eps(1)
	!eps(lmax) = eps(lmax1)

c---------------------------------------------------------------------
c turbulent kinetic energy (k)
c---------------------------------------------------------------------

	do l = 1,lmax1

	  ps = rkm(l) * m2(l)
	  gg = - rkt(l) * n2(l)
	  forc = dt * ( ps + gg - eps(l) )
	  aux = dt*dzr(l)
	  aa(l) = - aux * (rkk(l)*dzk(l))
	  bb(l) = 1. + aux * ( rkk(l)*dzk(l) + rkk(l+1)*dzk(l+1) )
	  cc(l) = - aux * (rkk(l+1)*dzk(l+1))
	  rr(l) = k(l)

	  if( forc .lt. 0. ) then
	    bb(l) = bb(l) - forc / k(l)
	  else
	    rr(l) = rr(l) + forc
	  end if

	end do

c	---------------------------------
c	impose boundary conditions
c	---------------------------------

	aa(0) = 0.
	bb(0) = 1.
	cc(0) = 0.
	rr(0) = taus / ( rho(1)*sqrtcmu )
	aa(lmax) = 0.
	bb(lmax) = 1.
	cc(lmax) = 0.
	rr(lmax) = taub / ( rho(lmax)*sqrtcmu )

c	---------------------------------
c	solve system
c	---------------------------------

c	call tridag(aa,bb,cc,rr,k,gam,lmax+1)
	call tridagd(aa,bb,cc,rr,k,gam1,gam,lmax+1)

c---------------------------------------------------------------------
c limit k and epsilon
c
c we need first condition for wind driven layer, second for emergency
c   when no stratification is present
c---------------------------------------------------------------------

	do l=0,lmax
	  k(l) = max(k(l),kmin)
          if(n2(l).gt.0.) then
	    eps(l) = max(eps(l),0.2*k(l)*sqrt(n2(l)))     !ERROR corrected
          else
	    eps(l) = max(eps(l),epsmin)
          end if
	end do

c---------------------------------------------------------------------
c compute the turbulent mixing coefficients from k, eps and n2
c---------------------------------------------------------------------

	do l = 1,lmax1
	  k2 = k(l) * k(l)
	  reps = 1. / eps(l)
	  a = cmu06 * n2(l) * k2 * reps * reps	 	!alpha_N
	  a = max(a,amin)
	  a = min(a,amax)

	  su = ( cmu0 + a0 * a ) / ( a1 + ( a2 + a3 * a ) * a )
	  sb = b0 / ( b1 + b2 * a )

	  aux = cmu03 * k2 * reps
	  rkm(l) = aux * su
	  rkt(l) = aux * sb
	  !rkm(l) = max(rkm(l),avumol)		!molecular viscosity
	  !rkt(l) = max(rkt(l),avsmol)		!molecular diffusivity
	end do

	if( bdebug ) write(6,*) 'keps 21'
	rkm(0) = 0.
	rkm(lmax) = 0.
	rkt(0) = 0.
	rkt(lmax) = 0.

c---------------------------------------------------------------------
c end of routine
c---------------------------------------------------------------------

	end

c************************************************************************

	subroutine kepsin(lmax,k,eps,len,rkm,rkt
     +			,kold,epsold,lenold,rkmold,rktold)

c initializes arrays for keps routine

	implicit none

	integer lmax
	real k(0:lmax),eps(0:lmax),len(0:lmax)
	real rkm(0:lmax),rkt(0:lmax)
	real kold(0:lmax),epsold(0:lmax),lenold(0:lmax)
	real rkmold(0:lmax),rktold(0:lmax)

	real kmin,epsmin,lenmin
        real avumol,avtmol,avsmol
c	parameter(kmin=1.e-10,epsmin=1.e-12,lenmin=0.01)
	parameter(kmin=3.e-6,epsmin=5.e-10,lenmin=0.01)
        parameter(avumol=1.3e-6,avtmol=1.4e-7,avsmol=1.1e-9)

	integer l

	do l=0,lmax
          k(l) = kmin
          eps(l) = epsmin
          len(l) = lenmin
	  rkm(l) = avumol
	  rkt(l) = avsmol

	  kold(l) = k(l)
	  epsold(l) = eps(l)
	  lenold(l) = len(l)
	  rkmold(l) = rkm(l)
	  rktold(l) = rkt(l)
	end do

	rkm(0) = 0.
	rkt(0) = 0.
	rkm(lmax) = 0.
	rkt(lmax) = 0.

	end

c************************************************************************

        subroutine keps_compile_test

	implicit none

	integer lmax
	real dt,rho0,taus,taub
	real dl(1),u(1),v(1),rho(1)

        real k(0:1),eps(0:1),len(0:1),rkm(0:1),rkt(0:1)
        real kold(0:1),epsold(0:1),lenold(0:1),rkmold(0:1),rktold(0:1)

	lmax = 1

	call kepsin(lmax,k,eps,len,rkm,rkt
     +			,kold,epsold,lenold,rkmold,rktold)
	call keps(lmax,dt,rho0,taus,taub,dl,u,v,rho,rkm,rkt,k,eps)

        end

c********************************************************************
c        call keps_compile_test
c        end
c********************************************************************

