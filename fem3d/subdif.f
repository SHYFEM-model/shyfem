
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

c routines for diffusion
c
c contents :
c
c subroutine diffstab(dt,rk,istot,v1v,v2v,gamma)        stability of diffusion
c subroutine diffstab1(dt,rkv,istot,v1v,v2v,gamma)      checks diff.  stability
c subroutine diffweight                                 computes diff. weights
c subroutine difflimit(dt,rkv,istot,gammax)             limits diff. param.
c subroutine diffadjust(mode,rkv)                       adjusts diff. coeff.
c
c revision log :
c
c 14.01.2005	ggu	new file for diffusion routines (this file)
c 23.02.2005	ggu	new routines for smagorinski and green
c 15.03.2005	ggu	austau() from newtra copied here
c 08.11.2005	ggu	fixed wrong debug statement in austau()
c 23.03.2006	ggu	changed time step to real
c 27.01.2009	ggu	diffset() deleted
c 12.02.2010	ggu	diffweight() has new method -> idtype=0,1,2
c 17.02.2010	ggu	bug fix in diffweight()
c 23.03.2010	ggu	changed v6.1.1
c 08.04.2010	ggu	better error reporting in diffweight()
c 16.02.2011	ggu	in diffweight() use double precision
c 01.06.2011	ggu	bug fix in green() -> i instead ii
c 05.12.2013	ggu	changed VERS_6_1_70
c 19.12.2014	ggu	changed VERS_7_0_10
c 23.12.2014	ggu	changed VERS_7_0_11
c 19.01.2015	ggu	changed VERS_7_1_3
c 10.07.2015	ggu	changed VERS_7_1_50
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 18.09.2015	ggu	austau() not used anymore - eliminated
c 23.09.2015	ggu	changed VERS_7_2_4
c 01.04.2016	ggu	changed VERS_7_5_7
c 09.10.2017	ggu	changed VERS_7_5_33
c 05.12.2017	ggu	changed VERS_7_5_39
c 03.04.2018	ggu	changed VERS_7_5_43
c 16.02.2019	ggu	changed VERS_7_5_60
c
c*****************************************************************

        subroutine diffstab(dt,rk,istot,v1v,v2v,gamma)

c checks stability of diffusion (old, not used)

	use evgeom
	use basin

        implicit none

        real dt
        real rk
        integer istot
        real v1v(1),v2v(1)
        real gamma              !stability parameter -> must be < 1.

        integer k,ie,ii
        real alpha,beta,area,b,c,bmin,bmax
        real rkmin,rkmax

        do k=1,nkn
          v1v(k) = 0.
          v2v(k) = 0.
        end do

        alpha = 3. * dt * rk

        do ie=1,nel
          area = 12. * ev(10,ie)
          do ii=1,3
            k = nen3v(ii,ie)
            b = ev(3+ii,ie)
            c = ev(6+ii,ie)
            v1v(k) = v1v(k) + area * ( b*b + c*c )
            v2v(k) = v2v(k) + area
          end do
        end do

        bmin = 1.e+30
        bmax = -bmin

        do k=1,nkn
          beta = v1v(k) / v2v(k)
          bmax = max(bmax,beta)
          bmin = min(bmin,beta)
        end do

        rkmax = 1. / (3.*dt*bmin)
        rkmin = 1. / (3.*dt*bmax)

        gamma = alpha * bmax / istot

        write(6,*) 'diffstab: ',dt,rk,rkmin,rkmax,gamma

        end

c*************************************************************

        subroutine diffstab1(dt,rkv,istot,v1v,v2v,gamma)

c checks stability of diffusion (with variable diffusion coef.)

	use evgeom
	use basin

        implicit none

        real dt
        real rkv(1)
        integer istot
        real v1v(1),v2v(1)
        real gamma              !stability parameter -> must be < 1.

        integer k,ie,ii
        real alpha,beta,area,b,c,bmin,bmax
        real rkmin,rkmax

        do k=1,nkn
          v1v(k) = 0.
          v2v(k) = 0.
        end do

        do ie=1,nel
          alpha = 3. * dt * rkv(ie)
          area = 12. * ev(10,ie)
          do ii=1,3
            k = nen3v(ii,ie)
            b = ev(3+ii,ie)
            c = ev(6+ii,ie)
            v1v(k) = v1v(k) + alpha * area * ( b*b + c*c )
            v2v(k) = v2v(k) + area
          end do
        end do

        bmin = 1.e+30
        bmax = -bmin

        do k=1,nkn
          beta = v1v(k) / v2v(k)
          bmax = max(bmax,beta)
          bmin = min(bmin,beta)
        end do

        gamma = bmax / istot

        write(6,*) 'diffstab1: ',dt,bmin,bmax,gamma

        end

c*************************************************************

        subroutine diffweight

c computes weight for diffusion
c
c weights in main diagonal are positive => weights out of diag are negative

	use mod_diff_aux
	use evgeom
	use basin

        implicit none

	logical bdebug,berror
        integer k,ie,ii,iii,i
	integer ia,ib
        integer nchange,idtype
        double precision w,fact,eps
        double precision b(3),c(3)
	double precision bc_orig(3,3)
	double precision bc_adj(3,3)
	double precision bc(3,3)
	double precision wacu_aux(3)
	double precision wacu,wacu_max

	real getpar

c-----------------------------------------------------------------
c initialization
c-----------------------------------------------------------------

	idtype = 1	!0: original  1: adjust  2: new sym weights
	idtype = nint(getpar('idtype'))	!delete after tests
        nchange = 0
	fact = 2./3.
	eps = 1.e-6
	eps = 1.e-3
	bdebug = .false.
	wacu_max = 0.

        write(6,*) 'diffweight: computing weights'

c-----------------------------------------------------------------
c loop over elements
c-----------------------------------------------------------------

        do ie=1,nel

          do ii=1,3
            b(ii) = ev(3+ii,ie)
            c(ii) = ev(6+ii,ie)
          end do

c	  -----------------------------------------------------------------
c 	  type 0
c	  -----------------------------------------------------------------

          do ii=1,3
            do iii=1,3
              bc(iii,ii) = b(iii)*b(ii) + c(iii)*c(ii)
	      bc_orig(iii,ii) = bc(iii,ii)
            end do
	  end do

c	  -----------------------------------------------------------------
c 	  adjust matrix
c	  -----------------------------------------------------------------

	  if( idtype .eq. 1 ) then

c	    -----------------------------------------------------------------
c 	    type 1
c	    -----------------------------------------------------------------

            do ii=1,3
              do iii=1,3
                w = bc(iii,ii)
                if( w .gt. 0. .and. iii .ne. ii ) then
                      i = 6 - ii - iii
                      bc(iii,ii) = 0.
                      bc(i,ii) = bc(i,ii) - w
                      nchange = nchange + 1
                end if
              end do
              do iii=1,3
		bc_adj(iii,ii) = bc(iii,ii)
	      end do
	    end do

	  else if( idtype .eq. 2 ) then	!out of diag are negative

c	    -----------------------------------------------------------------
c 	    type 2
c	    -----------------------------------------------------------------

            do i=1,3
	      ia = 1+mod(i,3)
	      ib = 1+mod(i+1,3)
	      w = 1. / ev(13+i,ie)**2
	      w = w * fact
              bc(ia,ib) = -w		!bug fix 17.2.2010
              bc(ib,ia) = -w
	      bc(i,i) = 0.
	    end do

	    do ii=1,3
	      w = 0.
              do iii=1,3
	        w = w + bc(ii,iii)
	      end do
	      bc(ii,ii) = -w
	    end do

	    do ii=1,3
              do iii=1,3
		bc_adj(ii,iii) = bc(ii,iii)
	      end do
	    end do

	  end if

c	  -----------------------------------------------------------------
c 	  error handling and copy to wdifhv
c	  -----------------------------------------------------------------

	  berror = .false.
	  do ii=1,3
	    wacu = 0.
            do iii=1,3
	      w = bc(ii,iii)
	      wacu = wacu + w
	      if( ii .eq. iii ) then
	        if( w .le. 0 ) berror = .true.
	      else
	        if( w .gt. 0 ) berror = .true.
	      end if
	      if( abs(w-bc(iii,ii)) .gt. eps ) berror = .true.	!symmetric?
	      wdifhv(ii,iii,ie) = bc(ii,iii)
	    end do
	    if( abs(wacu) .gt. eps ) berror = .true.
	    wacu_max = max(wacu_max,abs(wacu))
	    wacu_aux(ii) = wacu
	  end do

	  if( berror .and. idtype .ne. 0 ) then
	  !if( berror ) then
	    write(6,*) 'diffweight: idtype = ',idtype
	    write(6,*) '   ie (intern) = ',ie
	    write(6,*) '   eps = ',eps
	    write(6,*) '   wacu = ',wacu_aux
	    do ii=1,3
	      write(6,*) (bc(iii,ii),iii=1,3)
	    end do
	    stop 'error stop diffweight: error in matrix'
	  end if

	  !bdebug = ie .eq. 100 .or. ie .eq. 101
	  if( bdebug ) then
	    write(6,*) 'diffusion check: ',ie
	    write(6,*) 'diffusion orig'
	    do ii=1,3
	      write(6,*) (bc_orig(iii,ii),iii=1,3)
	    end do
	    write(6,*) 'diffusion adj'
	    do ii=1,3
	      write(6,*) (bc_adj(iii,ii),iii=1,3)
	    end do
	    write(6,*) 'diffusion check end'
	  end if

        end do

c-----------------------------------------------------------------
c end of loop over elements
c-----------------------------------------------------------------

        write(6,*) 'diffweight: total weights changed = ', nchange
        write(6,*) 'diffweight: type of hor diffus    = ', idtype
        write(6,*) 'diffweight: maximum error         = ', wacu_max

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

        end

c*************************************************************

        subroutine difflimit(dt,rkv,istot,gammax)

c limits diffusion parameter

	use evgeom
	use basin

        implicit none

        real dt
        real rkv(1)
        integer istot
        real gammax             !max for stability parameter, should be < 1

        integer k,ie,ii
        integer nchange
        real gamma,rk
        real alpha,beta,area,b,c,bmin,bmax
        real rkmin,rkmax

        nchange = 0
        rkmin = rkv(1)
        rkmax = rkv(1)

        do ie=1,nel
          beta = 0.
          rk = rkv(ie)
          alpha = 3. * dt * rk
          do ii=1,3
            k = nen3v(ii,ie)
            b = ev(3+ii,ie)
            c = ev(6+ii,ie)
            beta = max( beta , alpha * ( b*b + c*c ) )
          end do
          gamma = beta / istot
          if( gamma .le. gammax ) then
            rkv(ie) = rk
          else
            nchange = nchange + 1
            rkv(ie) = rk * gammax / gamma
            rkmin = min(rkmin,rkv(ie))
            rkmax = max(rkmax,rkv(ie))
          end if
        end do

        write(6,*) 'difflimit: ',dt,rk,rkmin,rkmax,nchange

        end

c*************************************************************

        subroutine diffadjust(mode,rkv)

c adjusts diffusion coefficient

	use mod_depth
	use evgeom
	use basin

        implicit none

        integer mode
        real rkv(1)

        integer k,ie,ii
        real h,aux,fact
        real bmin,bmax,beta
        real rmin,rmax
        real hmin,hmax

        real alpha,area

        write(6,*) 'diffadjust : ',mode

        if( mode .le. 0 ) return

        if( mode .le. 3 ) then

        rmin = 1.
        rmax = 2.
        hmin = 10.
        hmax = 2.

        if( mode .eq. 2 ) rmax = 3.
        if( mode .eq. 3 ) rmax = 4.

        aux = (rmax-rmin)/(hmax-hmin)

        do ie=1,nel
          h = hev(ie)
          if( h .lt. hmin ) then
              fact = rmin + aux * (h-hmin)
              rkv(ie) = fact * rkv(ie)
              !write(6,*) ie,h,rkv(ie)
          end if
        end do

        else if( mode .eq. 4 ) then

          alpha = 0.01
          do ie=1,nel
            area = 12. * ev(10,ie)
            rkv(ie) = alpha * area**(2./3.)
          end do

        end if

        call mima(rkv,nel,bmin,bmax)

        write(6,*) 'diffadjust: rkh adjust min/max: ',bmin,bmax

        end

c*************************************************************************** 

        subroutine set_diffusivity

c sets the horizontal diffusion array
c
c idhtyp gives type of diffusion
c	0	constant
c	1	variable with area ( ah = alpha * dx**(4/3) )
c	2	smagorinsky (variable with area and time)

	use mod_diff_visc_fric
	use evgeom
	use levels
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        character*80 file,title
        integer ie,l,lmax
        real dt
        real alpha,ahmax,area,ah
        real dhlen,dhpar,chpar,thpar,shpar,ahpar
	real ve1v(nel)
	real v1v(nkn)
	real v2v(nkn)

        real getpar

        real parmax
        save parmax
        integer idhtyp
        save idhtyp

        integer icall
        save icall
        data icall /0/

c------------------------------------------------------------------
c set params
c------------------------------------------------------------------

        if( icall .lt. 0 ) return

c------------------------------------------------------------------
c time dependent diffusion
c	call this only during time iteration of simulation
c------------------------------------------------------------------

        if( icall .gt. 0 ) then         !time dependent diffusion
	  call smagorinsky
	  return
        end if

c------------------------------------------------------------------
c first call
c------------------------------------------------------------------

        idhtyp = nint(getpar('idhtyp'))
        dhlen = getpar('dhlen')

c       ------------------------------------------------
c       set up diffusion coefficients
c       ------------------------------------------------

        dhpar = getpar('dhpar')
        chpar = getpar('chpar')
        thpar = getpar('thpar')
        shpar = getpar('shpar')
        ahpar = getpar('ahpar')

        if( dhpar .lt. 0. ) dhpar = 0.
        if( chpar .lt. 0. ) chpar = dhpar
        if( thpar .lt. 0. ) thpar = dhpar
        if( shpar .lt. 0. ) shpar = dhpar
        if( ahpar .lt. 0. ) ahpar = 0.

        call putpar('dhpar',dhpar)
        call putpar('chpar',chpar)
        call putpar('thpar',thpar)
        call putpar('shpar',shpar)
        call putpar('ahpar',ahpar)

        parmax = max(dhpar,chpar,thpar,shpar,ahpar)

c       ------------------------------------------------
c       set up area dependent diffusion coefficient
c       ------------------------------------------------

        alpha = 0.
        if( idhtyp .eq. 1 ) then
          if( dhlen .le. 0. ) goto 99
          alpha = 1. / ( dhlen**(4./3.) )
        end if

c       ------------------------------------------------
c       set up time constant diffusion coefficient
c       ------------------------------------------------

        ahmax = 0.
        do ie=1,nel
          area = 12. * ev(10,ie)
          ah = 1.
          if( alpha .gt. 0. ) ah = alpha * area**(2./3.)
          ahmax = max(ahmax,ah)

          lmax = ilhv(ie)
          do l=1,lmax
            difhv(l,ie) = ah
          end do
        end do

c       ------------------------------------------------
c       finished initializing
c       ------------------------------------------------

        if( ahmax * parmax .gt. 1000. ) then
          write(6,*) 'Horizontal diffusion coefficient too high'
          write(6,*) '  ahmax,parmax,ahmax*parmax '
          write(6,*) ahmax,parmax,ahmax*parmax
          stop 'error stop diff_h_set: Horizontal diffusion'
        end if

        icall = -1
        if( idhtyp .ge. 2 ) then
	  icall = 1
	  call smagorinsky
	  write(6,*) 'initializing smagorinski...'
	end if

        write(6,*) 'horizontal diffusion (1): ',parmax,ahmax,dhlen
        write(6,*) 'horizontal diffusion (2): ',idhtyp,icall
        write(6,*) ' dhpar,chpar,thpar,shpar,ahpar : '
        write(6,*) dhpar,chpar,thpar,shpar,ahpar

c       ------------------------------------------------------------------
c       checks stability
c       ------------------------------------------------------------------

	call get_timestep(dt)

c	still to do difhv,parmax
c        call diffstab1(dt,cdifhv,istot,v1v,v2v,gamma)

c       ------------------------------------------------------------------
c       write file
c       ------------------------------------------------------------------

        do ie=1,nel
          ve1v(ie) = parmax * difhv(1,ie)
        end do

!        file = 'rkdiff'
!        title = 'horizontal diffusion coef'
!        call e2n2d(ve1v,v1v,v2v)
!        call wrnos2d(file,title,v1v)

c------------------------------------------------------------------
c end of routine
c------------------------------------------------------------------

        return
   99   continue
        write(6,*) 'dhlen = ',dhlen
        stop 'error stop diff_h_set: dhlen not allowed'
        end

c*************************************************************************** 
c
c subroutines for computing Smagorinsky Horizontal Diffusion
c               Coefficient AHG, AMG
c
c       Ahg=[(CD*Size/Pi)**2ABS(DD)]
c       DD=SQRT[DT**2+DS**2]
c       DT=Ux-Vy
c       DS=Uy+Vx
c
c Ahg is horizontal diffusivity coefficient        
c It is computed considering both the 
c Size of the spatial discretization,
c the horizontal tension DT and 
c the shearing stress DS of the velocity field
c CD is a parameter (Mellor CD=0.5)
c 
c To convert diffusivity coefficient in viscosity coefficient 
c Amg, introduce a CVD=Amg/Ahg factor = 0.2
c
c*************************************************************************** 

        subroutine smagorinsky

	use mod_diff_visc_fric
	use mod_hydro_print
	use evgeom
	use levels
	use basin

        implicit none

        real b(3),c(3),ux,uy,vx,vy
        real dt,ds,dd,dl,aj

	real area,smag
        
        integer k,l,ie,ii,lmax
        
c the FEM method gives for Ux(ie)=unv(ie)*b and for Uy(ie)=unv(ie)*c
       
        do ie=1,nel
          
	 lmax = ilhv(ie)
         area = 12. * ev(10,ie)
         
c compute the spatial derivates of horizontal velocity        
         
         do l=1,lmax

           ux=0.
           uy=0.
           vx=0.
           vy=0.
                     
           do ii=1,3
             k=nen3v(ii,ie)
             b(ii)=ev(ii+3,ie)
             c(ii)=ev(ii+6,ie)
             ux=ux+(uprv(l,k)*b(ii))
             uy=uy+(uprv(l,k)*c(ii))
             vx=vx+(vprv(l,k)*b(ii))
             vy=vy+(vprv(l,k)*c(ii))
           end do
        
	   smag = 2.*ux*ux + 2.*vy*vy + ( uy + vx ) **2
	   smag = area * sqrt(smag)

	   difhv(l,ie) = smag

         end do
        end do

c old part -> deleted
c
c computing of Dt tension strain and Ds shearing strain
c           dt=ux-vy
c           ds=vx+uy
c           dd=sqrt((dt**2)+(ds**2))
c
c computing the length scale of the ie-element
c
c           aj=ev(10,ie)
c           dl=(sqrt(12*aj))/pi
c
c computing horizontal diffusivity Ahg
c
c           ahg(l,ie)=((CD*dl)**2)*dd     
c
c computing horizontal viscosity Amg
c
c           amg(l,ie)=cvd*ahg(l,ie)
c           write(92,*)dl,ahg(l,ie)
         
	return
        end 

c***********************************************************************

	subroutine green(ie,l,ugreen,vgreen)

c solves green identity for reynolds stresses

c ieltv:  >0 element  0: boundary  -1: open boundary

	use mod_geom
	use mod_hydro
	use evgeom
	use levels
	use basin
	use shympi

	implicit none

	integer ie
	integer l
	real ugreen,vgreen	!contribution to integral from green formula

	integer k,i,ii,iii,ienb,i1,i2
	real dl(3)
	real x(3),y(3)
	real u,v,unb,vnb
	real xm,ym,xmb,ymb
	real dist
	real ugrad,vgrad

	real x1,y1,x2,y2
	real distance
	distance(x1,y1,x2,y2) = sqrt((x1-x2)**2+(y1-y2)**2)

	if( shympi_is_parallel() ) then
	  stop 'error stop green: not ready for mpi'
	end if

c------------------------------------------------------------
c get transports
c------------------------------------------------------------

	u = utlov(l,ie)
	v = vtlov(l,ie)

c------------------------------------------------------------
c compute geometric characteristics for element
c------------------------------------------------------------

c	---------------------------------------------
c	center point
c	---------------------------------------------

	do ii=1,3
	  k = nen3v(ii,ie)
	  x(ii) = xgv(k)
	  y(ii) = ygv(k)
	end do

	call baric(ie,xm,ym)

c	---------------------------------------------
c	length of sides
c	---------------------------------------------

	do ii=1,3
	  i1 = mod(ii,3) + 1
	  i2 = mod(i1,3) + 1
	  dl(ii) = distance(x(i1),y(i1),x(i2),y(i2))
	end do

c------------------------------------------------------------
c compute contribution
c------------------------------------------------------------

	do ii=1,3

c	  ------------------------------------------------------------
c	  find neibor element
c	  ------------------------------------------------------------

	  ienb = ieltv(ii,ie)
	  if( ienb .gt. 0 .and. ilhv(ienb) .lt. l ) ienb = 0

c	  ---------------------------------------------
c	  find velocity in neibor element
c	  ---------------------------------------------

	  if( ienb .lt. 0 ) then	!open boundary -> no friction
	    unb = u
	    vnb = v
	  else if( ienb .eq. 0 ) then	!material boundary -> 0 slip
	    unb = 0.
	    vnb = 0.
	    !uncomment next two lines for full slip condition
	    !unb = u
	    !vnb = v
	  else				!neibor element existing
	    unb = utlov(l,ienb)
	    vnb = vtlov(l,ienb)
	  end if

c	  ---------------------------------------------
c	  find a distance for gradient
c	  ---------------------------------------------

	  if( ienb .gt. 0 ) then

c	    ------------------------------------------------
c	    compute distance to center point of neibor element
c	    ------------------------------------------------

	    call baric(ienb,xmb,ymb)
	    dist = distance(xm,ym,xmb,ymb)

	  else

c	    ------------------------------------------------
c	    compute virtual distance
c	    ------------------------------------------------

	    dist = 0.5 * distance(xm,ym,x(ii),y(ii))

	  end if

c	  ---------------------------------------------
c	  compute gradient
c	  ---------------------------------------------

	  ugrad = (unb - u) / dist
	  vgrad = (vnb - v) / dist

c	  ---------------------------------------------
c	  final contribution
c	  ---------------------------------------------

	  ugreen = ugreen + ugrad * dl(ii)
	  vgreen = vgreen + vgrad * dl(ii)

	end do

c------------------------------------------------------------
c end of routine
c------------------------------------------------------------

	end

c***********************************************************************

