!
! $Id: kepsd.f,v 1.4 2008-12-09 11:38:33 georg Exp $
!
! routines for k-epsilon model
!
! revision log :
!
! 25.02.1999	ggu	finished writing routines
! 10.04.2008	ggu&ccf	new gotm routine call integrated
! 01.06.2011	ggu&cpb	call to do_gotm_turb was wrong (BUG)
! 10.02.2015	ggu	in gotm call to do_gotm_turb() with ddt and not dt
!
! notes :
!
! u,v,rho are evaluated at the center of the layer
! other variables (k,eps,rkm,rkt) at the top and bottom of the layers
!
!********************************************************************
!--------------------------------------------------------------------
        module kepsilon
!--------------------------------------------------------------------
        contains
!--------------------------------------------------------------------

	subroutine gotm(lmax,dt,rho0,taus,taub,dl,u,v,rho,rkm,rkt,k,eps,len)

! computes turbulent diffusion coefficients through gotm

	implicit none

! arguments
	integer lmax		!number of levels [1...lmax]
	double precision dt			!time step
	double precision rho0		!reference density
	double precision taus		!surface stress at new time level
	double precision taub		!bottom stress at new time level
	double precision dl(1:lmax)		!bottom depth of single levels
	double precision u(1:lmax)		!velocity of single levels in x direction
	double precision v(1:lmax)		!velocity of single levels in y direction
	double precision rho(1:lmax)	!density of single levels
	double precision rkm(0:lmax)	!eddy coefficient for momentum		(out)
	double precision rkt(0:lmax)	!eddy coefficient for scalar quantities (out)
	double precision k(0:lmax)		!turbulent kinetic energy		(out)
	double precision eps(0:lmax)	!turbulent dissipation			(out)
	double precision len(0:lmax)	!length scale				(out)

	integer nlevdi
	parameter(nlevdi=200)
	double precision g,karman
	parameter(g=9.81,karman=0.4)

	double precision dzk(nlevdi)		!z difference between k (eps) levels
	double precision dzr(0:nlevdi)		!z difference between rho (u,v) levels
	double precision dzl(0:nlevdi)		!needed for gotm
	double precision b(nlevdi)			!buoyancy
	double precision n2(0:nlevdi)		!Brunt-V. frequency
	double precision m2(0:nlevdi)		!shear production

	double precision hh(0:nlevdi)
	double precision nn(0:nlevdi)
	double precision ss(0:nlevdi)
	double precision num(0:nlevdi)
	double precision nuh(0:nlevdi)
	double precision rkin(0:nlevdi)
	double precision reps(0:nlevdi)
	double precision rlen(0:nlevdi)

	integer l,lmax1,lp
	double precision ddt,u_taus,u_taub,z0,ud,vd,aux
	double precision z0s,z0b,depth

	if( lmax .gt. nlevdi ) stop 'error stop: nlevdi'

	lmax1 = lmax - 1

	ddt = dt
        u_taus = sqrt( taus )	!u_star
        u_taub = sqrt( taub )
	depth = dl(lmax)
	z0s = 0.03
	z0b = 0.03

! compute level differences
!
! z0 is height of surface -> if highly variable -> pass into routine
!
! dzk == dz        dzr == dZ  (Langland and Liou, MWR, May 1996, 905-918)

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

! compute b(uoyancy) and Brunt-Vaisala frequency n2 (N^2)

	aux = -g / rho0
	do l = 1,lmax
	  b(l) = aux * ( rho(l) - rho0 )
	end do

	do l = 1,lmax1
	  n2(l) = ( b(l) - b(l+1) ) * dzr(l)
	  n2(l) = max(0.D+0,n2(l))
	end do

! compute horizontal shear m2 (M^2)

	do l = 1,lmax1
	  ud = ( u(l) - u(l+1) ) * dzr(l)
	  vd = ( v(l) - v(l+1) ) * dzr(l)
	  m2(l) = ud*ud + vd*vd
	end do

! prepare data for gotm call

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

!	call gotmturb(lmax,ddt,hh,NN,SS,num,nuh
!     +			,rkin,reps,rlen,u_taus,u_taub)

	call do_gotm_turb(lmax,ddt,depth,u_taus,u_taub,z0s,z0b,hh,NN,SS,num,nuh,rkin,reps,rlen)

	do l=0,lmax
	  lp = lmax-l
	  rkm(lp) = num(l)
	  rkt(lp) = nuh(l)
	  k(lp)   = rkin(l)
	  eps(lp) = reps(l)
	  len(lp) = rlen(l)
	end do

	end

!********************************************************************

	subroutine ma(lmax,dt,rho0,taus,taub,dl,u,v,rho,rkm,rkt)

! computes turbulent diffusion coefficients with Munk-Anderson theory

	implicit none

! arguments
	integer lmax		!number of levels [1...lmax]
	double precision dt			!time step
	double precision rho0		!reference density
	double precision taus		!surface stress at new time level
	double precision taub		!bottom stress at new time level
	double precision dl(1:lmax)		!bottom depth of single levels
	double precision u(1:lmax)		!velocity of single levels in x direction
	double precision v(1:lmax)		!velocity of single levels in y direction
	double precision rho(1:lmax)	!density of single levels
	double precision rkm(0:lmax)	!eddy coefficient for momentum		(out)
	double precision rkt(0:lmax)	!eddy coefficient for scalar quantities (out)

	integer nlevdi
	parameter(nlevdi=200)
	double precision a0,k0
	parameter(a0=10.e-4,k0=5.e-4)
	double precision alpha,beta
	parameter(alpha=-0.5,beta=-1.5)
	double precision aa,bb
	parameter(aa=10,bb=3.33)
	double precision g,karman
	parameter(g=9.81,karman=0.4)

	double precision b(nlevdi)			!buoyancy
	double precision n2(0:nlevdi)		!Brunt-V. frequency
	double precision m2(0:nlevdi)		!shear production
	double precision dzk(nlevdi)		!z difference between k (eps) levels
	double precision dzr(0:nlevdi)		!z difference between rho (u,v) levels

	integer l,lmax1
	double precision aux
	double precision z0,ud,vd,ri

	if( lmax .gt. nlevdi ) stop 'error stop: nlevdi'

	lmax1 = lmax - 1

! compute level differences
!
! z0 is height of surface -> if highly variable -> pass into routine
!
! dzk == dz        dzr == dZ  (Langland and Liou, MWR, May 1996, 905-918)

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

! compute b(uoyancy) and Brunt-Vaisala frequency n2 (N^2)

	aux = -g / rho0
	do l = 1,lmax
	  b(l) = aux * ( rho(l) - rho0 )
	end do

	do l = 1,lmax1
	  n2(l) = ( b(l) - b(l+1) ) * dzr(l)
	  n2(l) = max(0.,n2(l))
	end do

! compute horizontal shear m2 (M^2)

	do l = 1,lmax1
	  ud = ( u(l) - u(l+1) ) * dzr(l)
	  vd = ( v(l) - v(l+1) ) * dzr(l)
	  m2(l) = ud*ud + vd*vd
	end do

! compute new turbulent viscosities and diffusivities

	do l=1,lmax1
	  ri = 0.
	  if( m2(l) .gt. 0 ) ri = n2(l)/m2(l)
	  rkm(l) = a0*(1.+aa*ri)**alpha  
	  rkt(l) = k0*(1.+bb*ri)**beta  
	  !write(6,*) '*** ',l,n2(l),m2(l),ri,rkm(l),rkt(l)
	end do

! boundary conditions

	rkm(0) = rkm(1)
	rkt(0) = rkt(1)
	rkm(lmax) = rkm(lmax1)
	rkt(lmax) = rkt(lmax1)

	end

!********************************************************************

        subroutine keps(lmax,dt,rho0,taus,taub,dh,u,v,rho,rkm,rkt,k,eps)

! computes turbulent diffusion coefficients with k-epsilon theorie

        use keps_util

        implicit none

! arguments
        integer lmax            !number of levels [1...lmax]
        double precision dt                 !time step
        double precision rho0               !reference density
        double precision taus               !surface stress at new time level
        double precision taub               !bottom stress at new time level
        double precision dh(1:lmax)         !thickness of single levels
	!double precision dl(1:lmax)	!bottom depth of single levels
        double precision u(1:lmax)          !velocity of single levels in x direction
        double precision v(1:lmax)          !velocity of single levels in y direction
        double precision rho(1:lmax)        !density of single levels
        double precision rkm(0:lmax)        !eddy coefficient for momentum		(out)
        double precision rkt(0:lmax)        !eddy coefficient for scalar quantities (out)
        double precision k(0:lmax)          !turbulent kinetic energy		(out)
        double precision eps(0:lmax)        !turbulent dissipation			(out)
! parameters
!	include 'param.h'
        integer nlevdi
        parameter(nlevdi=200)

        double precision g,karman
        double precision avumol,avtmol,avsmol
        double precision kmin,epsmin
        double precision cmu,sqrtcmu
        double precision sigmak,sigmae
        double precision ceps1,ceps2
        double precision ceps3p,ceps3n,ceps3f
        double precision cmu0,cmu03,cmu06
        double precision amin,amax
        double precision a0,a1,a2,a3,b0,b1,b2
        parameter(g=9.81,karman=0.4)
        parameter(avumol=1.3e-6,avtmol=1.4e-7,avsmol=1.1e-9)
        parameter(kmin=1.e-10,epsmin=1.e-12)
!	parameter(kmin=3.e-6,epsmin=5.e-10)
!---------------------------------------------------------------- Patrick
!	parameter(cmu=0.091,sqrtcmu=0.3)
!	parameter(sigmak=1.0,sigmae=1.3)
!	parameter(ceps1=1.51,ceps2=1.92)
!	parameter(ceps3p=0.2,ceps3n=1.0,ceps3f=ceps1)
!	parameter(cmu0=1.0,cmu03=cmu0**3,cmu06=cmu03**2)
!	parameter(amin=-1.37,amax=20.4)
!	parameter(a0=0.023,a1=1.0,a2=0.714,a3=0.067)
!	parameter(b0=0.125,b1=1.0,b2=0.603)
!---------------------------------------------------------------- Hans
        parameter(cmu=0.091,sqrtcmu=0.3)
        parameter(sigmak=1.0,sigmae=1.08)
        parameter(ceps1=1.44,ceps2=1.92)
        parameter(ceps3p=-0.4,ceps3n=1.0,ceps3f=1.0)
        parameter(cmu0=0.5562,cmu03=cmu0**3,cmu06=cmu03**2)
        parameter(amin=-0.0466,amax=0.56)
        parameter(a0=2.182,a1=1.0,a2=20.40,a3=53.12)
        parameter(b0=0.6985,b1=1.0,b2=17.34)
!----------------------------------------------------------------
! local variables
        logical bdebug
        integer l,lmax1,lm1
        double precision aux,ud,vd
        double precision ps,gg,forc
        double precision a,su,sb,k2,reps
        double precision ceps3

        double precision b(nlevdi)                  !buoyancy
        double precision n2(0:nlevdi)               !Brunt-V. frequency
        double precision m2(0:nlevdi)               !shear production
        double precision dzk(nlevdi)                !z difference between k (eps) levels
        double precision dzr(0:nlevdi)              !z difference between rho (u,v) levels
        double precision rke(nlevdi)                !eddy coefficient for epsilon
        double precision rkk(nlevdi)                !eddy coefficient for kinetic energy

        double precision aa(0:nlevdi)             !arrays for tri-diagonal algorithm
        double precision bb(0:nlevdi)             !...
        double precision cc(0:nlevdi)             !...
        double precision rr(0:nlevdi)             !...
        double precision gam(0:nlevdi)            !...
        double precision gam1(0:nlevdi)           !...

! u,v,rho are evaluated at the center of the layer
! other variables (k,eps,rkm,rkt) at the top and bottom of the layers
!
! dh(1) contains already surface level variation
! dzk == dz        dzr == dZ  (Langland and Liou, MWR, May 1996, 905-918)

!---------------------------------------------------------------------
! start of code
!---------------------------------------------------------------------

        if( lmax .gt. nlevdi ) stop 'error stop: nlevdi'

        bdebug = .true.
        bdebug = .false.
        lmax1 = lmax - 1
        lm1 = lmax + 1

!---------------------------------------------------------------------
! compute level differences
!---------------------------------------------------------------------

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
        dzk(l) = 1. / dzk(l)            !last level -> l = lmax

!---------------------------------------------------------------------
! set eddy coefficient for k and eps diffusion
! rkm(0) and rkm(lmax) may not be set -> avoid using them
!---------------------------------------------------------------------

        rkk(1) = rkm(1) / sigmak
        rke(1) = rkm(1) / sigmae
        do l=2,lmax1
          aux = rkm(l-1) + rkm(l)
          rkk(l) = ( 0.5 / sigmak ) * aux
          rke(l) = ( 0.5 / sigmae ) * aux
        end do
        rkk(lmax) = rkm(lmax1) / sigmak
        rke(lmax) = rkm(lmax1) / sigmae

!---------------------------------------------------------------------
! compute b(uoyancy) and Brunt-Vaisala frequency n2 (N^2)
!---------------------------------------------------------------------

        aux = -g / rho0
        do l = 1,lmax
          b(l) = aux * ( rho(l) - rho0 )
        end do

        do l = 1,lmax1
          n2(l) = ( b(l) - b(l+1) ) * dzr(l)
          n2(l) = max(0.,n2(l))
        end do

!---------------------------------------------------------------------
! compute horizontal shear m2 (M^2)
!---------------------------------------------------------------------

        do l = 1,lmax1
          ud = ( u(l) - u(l+1) ) * dzr(l)
          vd = ( v(l) - v(l+1) ) * dzr(l)
          m2(l) = ud*ud + vd*vd
        end do

!---------------------------------------------------------------------
! turbulent dissipation (eps)
!---------------------------------------------------------------------

        do l = 1,lmax-1

          ps = rkm(l) * m2(l)
          gg = - rkt(l) * n2(l)

          if( n2(l) .gt. 0. ) then
            ceps3 = ceps3f * ceps3p
          else
            ceps3 = ceps3f * ceps3n
          end if

          forc = ( dt / k(l) ) * (ceps1 * ps + ceps3 * gg - ceps2 * eps(l))
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

!	---------------------------------
!	impose boundary conditions
!	---------------------------------

        aa(0) = 0.
        bb(0) = 1.
        cc(0) = 0.
        rr(0) = (taus/rho(1))**1.5 * 0.00002 / karman   !ERROR corrected
        rr(0) = (taus/rho(1))**1.5 / (dzk(1)*karman)   !ERROR corrected
        aa(lmax) = 0.
        bb(lmax) = 1.
        cc(lmax) = 0.
        rr(lmax) = (taub/rho(lmax))**1.5 / (dzk(lmax)*karman)

!	---------------------------------
!	solve system
!	---------------------------------

!	call tridag(aa(1),bb(1),cc(1),rr(1),eps(1),gam(1),lmax-1)
        call tridagd(aa(0),bb(0),cc(0),rr(0),eps(0),gam1(0),gam(0),lm1)

	!eps(0) = eps(1)
	!eps(lmax) = eps(lmax1)

!---------------------------------------------------------------------
! turbulent kinetic energy (k)
!---------------------------------------------------------------------

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

!	---------------------------------
!	impose boundary conditions
!	---------------------------------

        aa(0) = 0.
        bb(0) = 1.
        cc(0) = 0.
        rr(0) = taus / ( rho(1)*sqrtcmu )
        aa(lmax) = 0.
        bb(lmax) = 1.
        cc(lmax) = 0.
        rr(lmax) = taub / ( rho(lmax)*sqrtcmu )

!	---------------------------------
!	solve system
!	---------------------------------

!	call tridag(aa,bb,cc,rr,k,gam,lmax+1)
        call tridagd(aa,bb,cc,rr,k,gam1,gam,lmax+1)

!---------------------------------------------------------------------
! limit k and epsilon
!
! we need first condition for wind driven layer, second for emergency
!   when no stratification is present
!---------------------------------------------------------------------

        do l=0,lmax
          k(l) = max(k(l),kmin)
          if(n2(l).gt.0.) then
            eps(l) = max(eps(l),0.2*k(l)*sqrt(n2(l)))     !ERROR corrected
          else
            eps(l) = max(eps(l),epsmin)
          end if
        end do

!---------------------------------------------------------------------
! compute the turbulent mixing coefficients from k, eps and n2
!---------------------------------------------------------------------

        do l = 1,lmax1
          k2 = k(l) * k(l)
          reps = 1. / eps(l)
          a = cmu06 * n2(l) * k2 * reps * reps          !alpha_N
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
        rkm(0)          = 0.d0
        rkm(lmax)       = 0.d0
        rkt(0)          = 0.d0
        rkt(lmax)       = 0.d0

!---------------------------------------------------------------------
! end of routine
!---------------------------------------------------------------------

	end

!************************************************************************

	subroutine kepsin(lmax,k,eps,len,rkm,rkt,kold,epsold,lenold,rkmold,rktold)

! initializes arrays for keps routine

	implicit none

	integer lmax
	double precision k(0:lmax),eps(0:lmax),len(0:lmax)
	double precision rkm(0:lmax),rkt(0:lmax)
	double precision kold(0:lmax),epsold(0:lmax),lenold(0:lmax)
	double precision rkmold(0:lmax),rktold(0:lmax)

	double precision kmin,epsmin,lenmin
        double precision avumol,avtmol,avsmol
!	parameter(kmin=1.e-10,epsmin=1.e-12,lenmin=0.01)
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

!************************************************************************

        subroutine keps_compile_test

	implicit none

	integer lmax
	double precision dt,rho0,taus,taub
	double precision dl(1),u(1),v(1),rho(1)

        double precision k(0:1),eps(0:1),len(0:1),rkm(0:1),rkt(0:1)
        double precision kold(0:1),epsold(0:1),lenold(0:1),rkmold(0:1),rktold(0:1)

	lmax = 1

	call kepsin(lmax,k,eps,len,rkm,rkt,kold,epsold,lenold,rkmold,rktold)
	call keps(lmax,dt,rho0,taus,taub,dl,u,v,rho,rkm,rkt,k,eps)

        end

!********************************************************************
!        call keps_compile_test
!        end
!********************************************************************

!--------------------------------------------------------------------
        end module kepsilon
!--------------------------------------------------------------------
