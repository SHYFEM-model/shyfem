!
! $Id: newbcl0.f,v 1.7 1998/09/07 09:42:41 georg Exp $
!
! equation of state routine
!
! contents :
!
! subroutine rhotst			to test function sigma
! double precision function sigma(s,t,p)		computes density
!
! revision log :
!
! revised 30.08.95	$$AUST - austausch coefficient introduced
! revised 11.10.95	$$BCLBND - boundary condition for barocliic runs
! 14.08.1998	ggu	new routine tsmass
! 19.08.1998	ggu	new routines to write NOS file
! 19.08.1998	ggu	call to barcfi changed
! 26.08.1998	ggu	call to bclini changed
! 26.08.1998	ggu	routines deleted: barcpr, chkuvw, vmima
! 26.08.1998	ggu	routines mimari transferred to subssv
! 26.08.1998    ggu     subroutine tsmass transferred to newchk
! 26.08.1998    ggu     subroutine stmima transferred to newut
! 26.08.1998    ggu     subroutines bclini, barcfi, bclbnd deleted
!
!**********************************************************************
!
!	call rhotst
!	end
!
!**********************************************************************
!
!----------------------------------------------------------------------
        module eq_state
!----------------------------------------------------------------------
        contains
!----------------------------------------------------------------------

	subroutine rhotst
!
! to test function sigma
!
	implicit none
!
	integer itot
	parameter (itot=12)
!
	integer i,is,it,ip
	double precision sd,td,pd
	double precision rho,aux
	double precision rhoref(8)
	double precision s(2),t(2),p(2)
	double precision sv(itot),tv(itot),pv(itot)
!
	data sv /0.,0.,35.,35.,0.,0.,35.,35.  ,35.,35.,35.,35./
	data tv /0.,0.,5.,5.,0.,0.,5.,5.  ,5.,5.,5.,5./
	data pv /100.,101.,100.,101.,1000.,1010.,1000.,1010.,0.,0.1, 100.,100.1/
!
	data s /0.,35./
	data t /5.,25./
	data p /0.,1000./
	data rhoref/ 999.96675,1044.12802, 997.04796,1037.90204 &
     &		   ,1027.67547,1069.48914,1023.34306,1062.53817 /
!
	write(6,1000)
!
	i=0
	do is=1,2
	  do it=1,2
	    do ip=1,2
		i=i+1
		sd=s(is)
		td=t(it)
		pd=p(ip)
		rho = sigma(sd,td,pd)
		write(6,1010) sd,td,pd,rho,real(rhoref(i)-1000.)
	    end do
	  end do
	end do
!
	aux=0.
	do i=1,itot
	  rho = sigma(sv(i),tv(i),pv(i))
	  write(6,1010) sv(i),tv(i),pv(i),rho,rho-aux
	  aux=rho
	  if(mod(i,2).eq.0) aux=0.
	end do
!
	stop
 1000	format('       salinity    temperature       pressure','        density','    ref-density')
 1010	format(5f15.6)
	end
!
!*****************************************************************
!
!---------------------------------------------------------------
!
!         d e n s i t y   o f   s e a - w a t e r
!
!           international unesco equation of state 1980
!                  (double precision version)
!
!            test data:
!
!          s=0.  t=5.  p=0.    rho= 999.96675
!          s=0.  t=5.  p=1000. rho=1044.12802
!          s=0.  t=25. p=0.    rho= 997.04796
!          s=0.  t=25. p=1000. rho=1037.90204
!          s=35. t=5.  p=0.    rho=1027.67547
!          s=35. t=5.  p=1000. rho=1069.48914
!          s=35. t=25. p=0.    rho=1023.34306
!          s=35. t=25. p=1000. rho=1062.53817
!
!          1000 m ~ 1010 dbar
!
!            input units (double precision) :
!
!         p - must be given as depth in 0.1*meters (bars)
!         s -       "       as salinity in promille
!         t -       "       as degree centigrade
!
!            output units (double precision) :
!
!         sigma = rho - 1000.   with rho in kg/m**3
!
!------------------------------------------------------------------
!
	double precision function sigma(s,t,p)
!
! computes density
!
	implicit none
!
	double precision s,t,p
	double precision t1,t2,t3,t4,t5
	double precision s1,s2,s32
	double precision p1,p2
	double precision rho0,sbmk
!
	double precision one,thousd
	parameter( one = 1.0d+0 , thousd = 1.0d+3 )
!
	double precision r00
	double precision r01,r02,r03,r04,r05,r06,r07
	double precision r08,r09,r10,r11,r12,r13,r14
!
	parameter( r00 = +999.842594d+0  )
	parameter( r01 = +6.793952d-2 , r02 = -9.095290d-3 )
	parameter( r03 = +1.001685d-4 , r04 = -1.120083d-6 )
	parameter( r05 = +6.536332d-9 , r06 = +8.24493d-1  )
	parameter( r07 = -4.0899d-3   , r08 = +7.6438d-5   )
	parameter( r09 = -8.2467d-7   , r10 = +5.3875d-9   )
	parameter( r11 = -5.72466d-3  , r12 = +1.0227d-4   )
	parameter( r13 = -1.6546d-6   , r14 = +4.8314d-4   )
!
	double precision k00
	double precision k01,k02,k03,k04,k05,k06
	double precision k07,k08,k09,k10,k11
	double precision k13,k14,k15,k16,k17,k18,k19
	double precision k20,k21,k22,k23,k24,k25,k26
!
	parameter( k00 = +19652.21d+0 )
	parameter( k01 = +148.4206d+0 , k02 = -2.327105d+0 )
	parameter( k03 = +1.360477d-2 , k04 = -5.155288d-5 )
	parameter( k05 = +3.239908d+0 , k06 = +1.43713d-3  )
	parameter( k07 = +1.16092d-4  , k08 = -5.77905d-7  )
	parameter( k09 = +8.50935d-5  , k10 = -6.12293d-6  )
	parameter( k11 = +5.2787d-8 )
	parameter( k13 = +54.6746d+0  , k14 = -0.603459d+0 )
	parameter( k15 = +1.09987d-2  , k16 = -6.1670d-5   )
	parameter( k17 = +7.944d-2    , k18 = +1.6483d-2   )
	parameter( k19 = -5.3009d-4   , k20 = +2.2838d-3   )
	parameter( k21 = -1.0981d-5   , k22 = -1.6078d-6   )
	parameter( k23 = +1.91075d-4  , k24 = -9.9348d-7   )
	parameter( k25 = +2.0816d-8   , k26 = +9.1697d-10  )
!
	t1=t
	t2=t1*t1
	t3=t2*t1
	t4=t3*t1
	t5=t4*t1
	s1=s
	s2=s1*s1
!	s32=sqrt(double precision(s2*s1))
	s32=dsqrt(s2*s1)		!ggu
	p1=p
	p2=p1*p1
!
!--------------------> density at standard atmospheric pressure (p=0)
!
	rho0 =   r00 + r01*t1 + r02*t2 + r03*t3 + r04*t4 + r05*t5       &
     &	     + ( r06 + r07*t1 + r08*t2 + r09*t3 + r10*t4 ) * s1         &
     &       + ( r11 + r12*t1 + r13*t2 ) * s32 + r14*s2
!
!--------------------> secant bulk modulus
!
	sbmk = k00+ k01*t1 + k02*t2 + k03*t3 + k04*t4+ ( k05 + k06*t1 + k07*t2 + k08*t3 ) * p1  &
     &       + ( k09 + k10*t1 + k11*t2 ) * p2+ ( k13 + k14*t1 + k15*t2 + k16*t3 ) * s1          &
     &       + ( k17 + k18*t1 + k19*t2 + k23*p1 ) * s32 + ( k20 + k21*t1 + k22*t2 ) * p1 * s1   &
     &       + ( k24 + k25*t1 + k26*t2 ) * p2 * s1
!
!--------------------> effect of sal. and temp. : sigma
!
	sigma = rho0/(one-p/sbmk) - thousd
!
	return
	end
!
!**********************************************************************

	double precision function sigma_JM95(s,t,p)
!
! computes density using equation of state by Jackett and McDougall 1995 
! coded by Mehmet Ilicak (milicak@itu.edu.tr)
! input t for potential temperature
!       s for salinity (psu)
!       p for pressure (dbar)
!
	implicit none
!
	double precision s,t,p
	double precision one,thousd
	parameter( one = 1.0d+0 , thousd = 1.0d+3 )
	double precision press_standard 
	parameter( press_standard = 0.0d0 )
	double precision epsln 
	parameter( epsln = 1.0d-32 )
!
	double precision a0
	double precision a1,a2,a3,a4,a5,a6,a7
	double precision a8,a9,a10,a11

	parameter( a0 = +9.9984085444849347d+02 )
	parameter( a1 = +7.3471625860981584d+00 )
        parameter( a2 = -5.3211231792841769d-02 )
        parameter( a3 = 3.6492439109814549d-04 )
        parameter( a4 = 2.5880571023991390d+00 ) 
        parameter( a5 = -6.7168282786692355d-03 )
        parameter( a6 = 1.9203202055760151d-03 )
        parameter( a7 = 1.1798263740430364d-02 )
        parameter( a8 = 9.8920219266399117d-08 )
        parameter( a9 = 4.6996642771754730d-06 )
        parameter( a10 = -2.5862187075154352d-08 )
        parameter( a11 = -3.2921414007960662d-12 )

	double precision b0
	double precision b1,b2,b3,b4,b5,b6,b7
	double precision b8,b9,b10,b11,b12

        parameter( b0 = 1.0000000000000000d+00 )
        parameter( b1 = 7.2815210113327091d-03 )
        parameter( b2 = -4.4787265461983921d-05 ) 
        parameter( b3 = 3.3851002965802430d-07 )
        parameter( b4 = 1.3651202389758572d-10 )
        parameter( b5 = 1.7632126669040377d-03 )
        parameter( b6 = -8.8066583251206474d-06 )
        parameter( b7 = -1.8832689434804897d-10 )
        parameter( b8 = 5.7463776745432097d-06 )
        parameter( b9 = 1.4716275472242334d-09)
        parameter( b10 = 6.7103246285651894d-06 )
        parameter( b11 = -2.4461698007024582d-17 )
        parameter( b12 = -9.1534417604289062d-18 )

	double precision t1,t2,p1
	double precision s1,sp5,p1t1
        double precision num,den,den2


        t1 = t
        t2 = t1*t1

        s1 = s
        sp5 = sqrt(s1)

        p1 = p - press_standard
        p1t1 = p1*t1

        num = a0 + t1*(a1 + t1*(a2+a3*t1) ) + s1*(a4 + a5*t1  + a6*s1)  &
     &        + p1*(a7 + a8*t2 + a9*s1 + p1*(a10+a11*t2))

        den = b0 + t1*(b1 + t1*(b2 + t1*(b3 + t1*b4))) + s1*(b5 + t1*(b6 + b7*t2) + sp5*(b8 + b9*t2)) &
     &        + p1*(b10 + p1t1*(b11*t2 + b12*p1))


        den2 = one/(epsln+den) 

!
!--------------------> effect of sal. and temp. : sigma
!
	sigma_JM95 = num*den2 - thousd
!
	return
	end

!**********************************************************************
        double precision function sigma_linear(s,t,p)
!
! computes density linear with temperature 
!
        implicit none

        double precision s,t,p
        double precision rho, rho0
        double precision one,thousd
        parameter( one = 1.0d+0 , thousd = 1.0d+3 )
        parameter( rho0 = 1.000d+3  )

!       Thermal expansion coefficient
        double precision alpha
        parameter( alpha = 2.0d-04 )

!       Reference temperature
        double precision t0
        parameter( t0 = 5.0000000000000000d+00 )

        rho = rho0 * ( one - alpha * ( t - t0 ))    
!
        sigma_linear = rho - thousd
!
        return
        end
!********************************************************

        subroutine rhoset_shell

        use shympi

! sets rho iterating to double precision solution

        implicit none

        logical biter
        integer itermax,iter
        double precision eps,resid,resid_old

        itermax = 10
        eps = 1.e-7
        eps = 1.e-14    !ivb
        eps = 1.e-25    !ivb
        eps = 1.e-30    !ivb


        biter   = .true.
        iter    = 0
        resid   = 0.
        resid_old = 0.

        do while( biter )
          resid_old = resid
          call rhoset(resid)
          iter = iter + 1
          if (iter .eq. 1)      biter = .false.      !ivb do only once                 
          if( resid .lt. eps )  biter = .false.
          if( abs(resid-resid_old) .lt. eps ) biter = .false.
          if( iter .gt. itermax ) biter = .false.
        end do

        if( iter .gt. itermax ) then
          write(6,*) '*** warning: max iterations in rhoset_shell',resid
          call tsrho_check
        end if

        end
!********************************************************

        subroutine rhoset(resid)

! computes rhov and bpresv
!
! 1 bar = 100 kPascal ==> factor 1.e-5
! pres = rho0*g*(zeta-z) + bpresv
! with bpresv = int_{z}^{zeta}(g*rho_prime)dz
! and rho_prime = rho - rho_0 = sigma - sigma_0
!
! in bpresv() is bpresv as defined above
! in rhov()   is rho_prime (=sigma_prime)
!
! brespv() and rhov() are given at node and layer interface

        use layer_thickness
        use ts
        use levels
        use basin, only : nkn,nel,ngr,mbw
        use shympi
        use para
        use sigma_admin
        use mud_admin

        implicit none

        double precision resid
! parameter
        include 'param.h'
! common

        include 'femtime.h'
        include 'pkonst.h'


! local
        logical bdebug,debug,bsigma
        integer k,l,lmax
        integer nresid,nsigma,g_nresid
        double precision sigma0,rho0,pres,hsigma
        double precision depth,hlayer,hh
        double precision rhop,presbt,presbc,dpresc,rhop_tmp
        double precision salt
        double precision dresid,g_dresid
! functions
        integer eos_type

        integer ntot

        !eos_type = 1 ! Unesco 1980
        !eos_type = 2 ! Jackett & McDougall 1995 
        !eos_type = 3 ! Linearized Unesco (T) !ivb ! NOT TESTED!!      
        eos_type = nint(getpar('eostype'))      ! Type of Equation of State 

        rho0    = rowass
        sigma0  = rho0 - 1000.

        debug   =.false.
        bdebug  =.false.

        call get_sigma(nsigma,hsigma)
        bsigma = nsigma .gt. 0

        if(debug) write(6,*) sigma0,rowass,rho0

        nresid  = 0
        dresid  = 0.

        ntot    = nkn   !ivb handle for MPI
        if (shympi_partition_on_elements()) ntot = nkn_inner

        do k=1,nkn
          depth  = 0.d0
          presbc = 0.d0
          lmax = ilhkv(k)

          do l=1,lmax
            bsigma = l .le. nsigma

            hlayer = hdkov(l,k)
            if( .not. bsigma ) hlayer = hldv(l)

            hh          = 0.5 * hlayer
            depth       = depth + hh
            rhop        = rhov(l,k)             !rho^prime

            dpresc = rhop * grav * hh           !differential bc. pres.
            presbc = presbc + dpresc            !baroclinic pres. (mid-layer)
            presbt = rho0 * grav * depth        !barotropic pressure

            pres = 1.e-5 * ( presbt + presbc )  !pressure in bars (BUG)
        
            salt = max(0.,saltv(l,k))
            if (eos_type .eq. 1) then
                pres = 1.e-5 * ( presbt + presbc )      !pressure in bars (BUG)
                rhop = sigma(salt,tempv(l,k),pres) - sigma0
            elseif (eos_type .eq. 2) then
                pres = 1.e-4 * ( presbt + presbc )      !pressure in dbars for JM95
                rhop = sigma_JM95(salt,tempv(l,k),pres) - sigma0
            elseif (eos_type .eq. 3) then               !ivb !!! not tested !!!
                rhop = sigma_linear(salt,tempv(l,k),pres) - sigma0
            endif
            call set_rhomud(k,l,rhop)

            rhop_tmp    = rhov(l,k)

            rhov(l,k)   = rhop
            bpresv(l,k) = presbc

            depth       = depth + hh
            presbc      = presbc + dpresc            !baroclinic pres. (bottom-lay.)

            if (k .gt. ntot ) cycle             ! ivb do not overlap calc in MPI case

            nresid = nresid + 1
            dresid = dresid + (rhop_tmp-rhop)**2
          end do
        end do
        
        g_dresid = shympi_sum(dresid)   !ivb
        g_nresid = shympi_sum(nresid)   !ivb

        resid = g_dresid/g_nresid

        return
        end

!*******************************************************************	

	subroutine tsrho_check

! checks values of t/s/rho

	use ts
	use levels
        use chk_NaN
	use basin, only : nkn,nel,ngr,mbw
        use model3d_util

	implicit none

	include 'param.h'


        double precision smin,smax,tmin,tmax,rmin,rmax
	character*30 text

	text = '*** tsrho_check'

	call stmima(saltv,nkn,nlvdi,ilhkv,smin,smax)
	call stmima(tempv,nkn,nlvdi,ilhkv,tmin,tmax)
	call stmima(rhov,nkn,nlvdi,ilhkv,rmin,rmax)

	write(6,*) 'S   min/max: ',smin,smax
	write(6,*) 'T   min/max: ',tmin,tmax
	write(6,*) 'Rho min/max: ',rmin,rmax

	write(6,*) 'checking for Nans...'
        call check2Dr(nlvdi,nlv,nkn,saltv,-1.d0,+70.d0,text,'saltv')
        call check2Dr(nlvdi,nlv,nkn,tempv,-30.d0,+70.d0,text,'tempv')
        call check2Dr(nlvdi,nlv,nkn,rhov,-2000.d0,+2000.d0,text,'rhov')

	end

!*******************************************************************	
!*******************************************************************	



!----------------------------------------------------------------------
        end module eq_state
!----------------------------------------------------------------------
