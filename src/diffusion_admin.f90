!
! $Id: subdif.f,v 1.11 2010-02-17 11:57:28 georg Exp $
!
! routines for diffusion
!
! contents :
!
! subroutine diffstab(dt,rk,istot,v1v,v2v,gamma)        stability of diffusion
! subroutine diffstab1(dt,rkv,istot,v1v,v2v,gamma)      checks diff.  stability
! subroutine diffweight                                 computes diff. weights
! subroutine difflimit(dt,rkv,istot,gammax)             limits diff. param.
! subroutine diffadjust(mode,rkv)                       adjusts diff. coeff.
!
! revision log :
!
! 14.01.2005    ggu     new file for diffusion routines (this file)
! 23.02.2005    ggu     new routines for smagorinski and green
! 15.03.2005    ggu     austau() from newtra copied here
! 08.11.2005    ggu     fixed wrong debug statement in austau()
! 23.03.2006    ggu     changed time step to double precision
! 27.01.2009    ggu     diffset() deleted
! 12.02.2010    ggu     diffweight() has new method -> idtype=0,1,2
! 17.02.2010    ggu     bug fix in diffweight()
! 08.04.2010    ggu     better error reporting in diffweight()
! 16.02.2011    ggu     in diffweight() use double precision
! 01.06.2011    ggu     bug fix in green() -> i instead ii
! 18.09.2015    ggu     austau() not used anymore - eliminated
!
!*****************************************************************
!---------------------------------------------------------------------
        module diffusion_admin
!---------------------------------------------------------------------
        contains
!---------------------------------------------------------------------

        subroutine diffstab(dt,rk,istot,v1v,v2v,gamma)

! checks stability of diffusion (old, not used)

	use evgeom
	use basin

        implicit none

        double precision dt
        double precision rk
        integer istot
        double precision v1v(1),v2v(1)
        double precision gamma              !stability parameter -> must be < 1.

	include 'param.h'

        integer k,ie,ii
        double precision alpha,beta,areael,b,c,bmin,bmax
        double precision rkmin,rkmax

        do k=1,nkn
          v1v(k) = 0.
          v2v(k) = 0.
        end do

        alpha = 3. * dt * rk

        do ie=1,nel
          areael = 12. * ev(10,ie)
          do ii=1,3
            k = nen3v(ii,ie)
            b = ev(3+ii,ie)
            c = ev(6+ii,ie)
            v1v(k) = v1v(k) + areael * ( b*b + c*c )
            v2v(k) = v2v(k) + areael
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

!*************************************************************

        subroutine diffstab1(dt,rkv,istot,v1v,v2v,gamma)

! checks stability of diffusion (with variable diffusion coef.)

	use evgeom
	use basin

        implicit none

        double precision dt
        double precision rkv(1)
        integer istot
        double precision v1v(1),v2v(1)
        double precision gamma              !stability parameter -> must be < 1.

	include 'param.h'

        integer k,ie,ii
        double precision alpha,beta,areael,b,c,bmin,bmax
        double precision rkmin,rkmax

        do k=1,nkn
          v1v(k) = 0.
          v2v(k) = 0.
        end do

        do ie=1,nel
          alpha = 3. * dt * rkv(ie)
          areael = 12. * ev(10,ie)
          do ii=1,3
            k = nen3v(ii,ie)
            b = ev(3+ii,ie)
            c = ev(6+ii,ie)
            v1v(k) = v1v(k) + alpha * areael * ( b*b + c*c )
            v2v(k) = v2v(k) + areael
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

!*************************************************************

        subroutine diffweight

! computes weight for diffusion
!
! weights in main diagonal are positive => weights out of diag are negative

	use diff_aux
	use evgeom
	use basin
        use para

        implicit none

	include 'param.h'

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

!-----------------------------------------------------------------
! initialization
!-----------------------------------------------------------------

	idtype = 1	!0: original  1: adjust  2: new sym weights
	idtype = nint(getpar('idtype'))	!delete after tests
        nchange = 0
	fact = 2./3.
	eps = 1.e-6
	eps = 1.e-3
	bdebug = .false.
	wacu_max = 0.

        write(6,*) 'diffweight: computing weights'

!-----------------------------------------------------------------
! loop over elements
!-----------------------------------------------------------------

        do ie=1,nel

          do ii=1,3
            b(ii) = ev(3+ii,ie)
            c(ii) = ev(6+ii,ie)
          end do

!	  -----------------------------------------------------------------
! 	  type 0
!	  -----------------------------------------------------------------

          do ii=1,3
            do iii=1,3
              bc(iii,ii) = b(iii)*b(ii) + c(iii)*c(ii)
	      bc_orig(iii,ii) = bc(iii,ii)
            end do
	  end do

!	  -----------------------------------------------------------------
! 	  adjust matrix
!	  -----------------------------------------------------------------

	  if( idtype .eq. 1 ) then

!	    -----------------------------------------------------------------
! 	    type 1
!	    -----------------------------------------------------------------

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

!	    -----------------------------------------------------------------
! 	    type 2
!	    -----------------------------------------------------------------

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

!	  -----------------------------------------------------------------
! 	  error handling and copy to wdifhv
!	  -----------------------------------------------------------------

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

!-----------------------------------------------------------------
! end of loop over elements
!-----------------------------------------------------------------

        write(6,*) 'diffweight: total weights changed = ', nchange
        write(6,*) 'diffweight: type of hor diffus    = ', idtype
        write(6,*) 'diffweight: maximum error         = ', wacu_max

!-----------------------------------------------------------------
! end of routine
!-----------------------------------------------------------------

        end

!*************************************************************

        subroutine difflimit(dt,rkv,istot,gammax)

! limits diffusion parameter

	use evgeom
	use basin

        implicit none

        double precision dt
        double precision rkv(1)
        integer istot
        double precision gammax             !max for stability parameter, should be < 1

	include 'param.h'

        integer k,ie,ii
        integer nchange
        double precision gamma,rk
        double precision alpha,beta,b,c,bmin,bmax
        double precision rkmin,rkmax

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

!*************************************************************

        subroutine diffadjust(mode,rkv)

! adjusts diffusion coefficient

	use depth
	use evgeom
	use basin
        use vec_util

        implicit none

        integer mode
        double precision rkv(1)

	include 'param.h'

        integer k,ie,ii
        double precision h,aux,fact
        double precision bmin,bmax,beta
        double precision rmin,rmax
        double precision hmin,hmax

        double precision alpha,areael

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
            areael = 12. * ev(10,ie)
            rkv(ie) = alpha * areael**(2./3.)
          end do

        end if

        call mima(rkv,nel,bmin,bmax)

        write(6,*) 'diffadjust: rkh adjust min/max: ',bmin,bmax

        end

!*************************************************************************** 

        subroutine set_diffusivity

! sets the horizontal diffusion array
!
! idhtyp gives type of diffusion
!	0	constant
!	1	variable with area ( ah = alpha * dx**(4/3) )
!	2	smagorinsky (variable with area and time)
!	3	Leith (variable with area and time)

        use diffusion
        use evgeom
        use levels
        use basin, only : nkn,nel,ngr,mbw
        use shympi
        use mpi_io_admin
        use para
        use nos_util
        use transforms
        use time_util

        implicit none

        include 'femtime.h'

        character*80 file,title
        integer ie,l,lmax
        double precision dt
        double precision alpha,ahmax,areael,ah
        double precision dhlen,dhpar,chpar,thpar,shpar,ahpar
        double precision ve1v(nel)
        double precision v1v(nkn)
        double precision v2v(nkn)

        double precision parmax
        save parmax
        integer idhtyp
        save idhtyp

        integer icall
        save icall
        data icall /0/

!------------------------------------------------------------------
! set params
!------------------------------------------------------------------

        if( icall .lt. 0 ) return

!------------------------------------------------------------------
! time dependent diffusion
!	call this only during time iteration of simulation
!------------------------------------------------------------------

        if( icall .gt. 0 ) then         !time dependent diffusion
          if( idhtyp .eq. 2 ) then
            call smagorinsky
            return
          elseif( idhtyp .eq. 3 ) then
            call leith
            return
          endif
        end if

!------------------------------------------------------------------
! first call
!------------------------------------------------------------------

        idhtyp = nint(getpar('idhtyp'))
        dhlen = getpar('dhlen')

!       ------------------------------------------------
!       set up diffusion coefficients
!       ------------------------------------------------

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

!       ------------------------------------------------
!       set up area dependent diffusion coefficient
!       ------------------------------------------------

        alpha = 0.
        if( idhtyp .eq. 1 ) then
          if( dhlen .le. 0. ) goto 99
          alpha = 1. / ( dhlen**(4./3.) )
        end if

!       ------------------------------------------------
!       set up time constant diffusion coefficient
!       ------------------------------------------------

        ahmax = 0.
        do ie=1,nel
          areael = 12. * ev(10,ie)
          ah = 1.
          if( alpha .gt. 0. ) ah = alpha * areael**(2./3.)
          ahmax = max(ahmax,ah)

          lmax = ilhv(ie)
          do l=1,lmax
            difhv(l,ie) = ah
          end do
        end do

!       ------------------------------------------------
!       finished initializing
!       ------------------------------------------------

        if( ahmax * parmax .gt. 5000. ) then
          write(6,*) 'Horizontal diffusion coefficient too high'
          write(6,*) '  ahmax,parmax,ahmax*parmax '
          write(6,*) ahmax,parmax,ahmax*parmax
          stop 'error stop diff_h_set: Horizontal diffusion'
        end if

        icall = -1              !this guarantees that it is not called anymore
        if( idhtyp .eq. 2 ) then
          icall = 1
          call smagorinsky
          write(6,*) 'initializing smagorinski...'
        elseif( idhtyp .eq. 3 ) then
          icall = 1
          call leith
          write(6,*) 'initializing Leith...'
        end if

        write(6,*) 'horizontal diffusion (1): ',parmax,ahmax,dhlen
        write(6,*) 'horizontal diffusion (2): ',idhtyp,icall
        write(6,*) ' dhpar,chpar,thpar,shpar,ahpar : '
        write(6,*) dhpar,chpar,thpar,shpar,ahpar

!       ------------------------------------------------------------------
!       checks stability
!       ------------------------------------------------------------------

        call get_timestep(dt)

!	still to do difhv,parmax
!        call diffstab1(dt,cdifhv,istot,v1v,v2v,gamma)

!       ------------------------------------------------------------------
!       write file
!       ------------------------------------------------------------------

        do ie=1,nel
          ve1v(ie) = parmax * difhv(1,ie)
        end do

        file = 'rkdiff'
        title = 'horizontal diffusion coef'
        call e2n2d(ve1v,v1v,v2v)

        call wrnos2d(file,title,v1v)

!------------------------------------------------------------------
! end of routine
!------------------------------------------------------------------

        return
   99   continue
        write(6,*) 'dhlen = ',dhlen
        stop 'error stop diff_h_set: dhlen not allowed'
        end

!*************************************************************************** 
!
! subroutines for computing Smagorinsky Horizontal Diffusion
!               Coefficient AHG, AMG
!
!       Ahg=[(CD*Size/Pi)²ABS(DD)]
!       DD=SQRT[DT²+DS²]
!       DT=Ux-Vy
!       DS=Uy+Vx
!
! Ahg is horizontal diffusivity coefficient        
! It is computed considering both the 
! Size of the spatial discretization,
! the horizontal tension DT and 
! the shearing stress DS of the velocity field
! CD is a parameter (Mellor CD=0.5)
! 
! To convert diffusivity coefficient in viscosity coefficient 
! Amg, introduce a CVD=Amg/Ahg factor = 0.2
!
!*************************************************************************** 

        subroutine smagorinsky

	use diffusion
	use hydro_print
	use evgeom
	use levels
	use basin

        implicit none

        include 'param.h'
        
	include 'femtime.h'
        
        
        double precision b(3),c(3),ux,uy,vx,vy
        double precision dt,ds,dd,dl,aj

	double precision areael,smag
        
        integer k,l,ie,ii,lmax
        
! the FEM method gives for Ux(ie)=unv(ie)*b and for Uy(ie)=unv(ie)*c
       
        do ie=1,nel
          
	 lmax = ilhv(ie)
         areael = 12. * ev(10,ie)
         
! compute the spatial derivates of horizontal velocity        
         
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
	   smag = areael * sqrt(smag)

	   difhv(l,ie) = smag

         end do
        end do

! old part -> deleted
!
! computing of Dt tension strain and Ds shearing strain
!           dt=ux-vy
!           ds=vx+uy
!           dd=sqrt((dt**2)+(ds**2))
!
! computing the length scale of the ie-element
!
!           aj=ev(10,ie)
!           dl=(sqrt(12*aj))/pi
!
! computing horizontal diffusivity Ahg
!
!           ahg(l,ie)=((CD*dl)**2)*dd     
!
! computing horizontal viscosity Amg
!
!           amg(l,ie)=cvd*ahg(l,ie)
!           write(92,*)dl,ahg(l,ie)
         
	return
        end 

!***********************************************************************
! coded by Mehmet Ilicak (milicak@itu.edu.tr)
! subroutines for computing Leith Horizontal Viscosity
!               Coefficient 
!
!       Ahg=[(CD*Size/Pi)^3 * sqrt(|G(curl)|^2+|G(cont)|^2)]
!       G is the grad operator
!       curl=-Uy+Vx
!       cont=Ux+Vy
!
! Ahg is horizontal diffusivity coefficient        
! It is computed considering both the 
! Size of the spatial discretization,
! the curl of the horizontal velocity and 
! the shearing stress DS of the velocity field

        subroutine leith

	use diffusion
	use hydro_print
	use evgeom
	use levels
	use basin
        use area
        use shympi

        implicit none

        include 'param.h'
        
	include 'femtime.h'
        
        
        double precision b(3),c(3),ux,uy,vx,vy
        double precision dt,ds,dd,dl,aj

        double precision areael, vsclth
        double precision q2D(nlv,nel), c2D(nlv,nel)
        double precision q2Dx, c2Dx
        double precision q2Dy, c2Dy
        
        integer k,l,ie,ii,lmax
        
! the FEM method gives for Ux(ie)=unv(ie)*b and for Uy(ie)=unv(ie)*c
       
        do ie=1,nel
          
	 lmax = ilhv(ie)
         
! compute the spatial derivates of horizontal velocity        
         
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
        
	   q2D(l,ie) = ( -uy + vx )
	   c2D(l,ie) = (  ux + vy )

         end do
        end do

! Interpolate relative vorticity (q2D) and continuity (c2D) into nodes        
        q2Dprv = 0.0
        c2Dprv = 0.0
        do ie=1,nel
	 lmax = ilhv(ie)
         areael = 4. * ev(10,ie)
         do l=1,lmax
           do ii=1,3
              k = nen3v(ii,ie)
              q2Dprv(l,k) = q2Dprv(l,k) + areael * q2D(l,ie)
              c2Dprv(l,k) = c2Dprv(l,k) + areael * c2D(l,ie)
           end do
         end do
        end do
        call shympi_exchange_and_sum_3D_nodes(q2Dprv)
        call shympi_exchange_and_sum_3D_nodes(c2Dprv)
        where ( areakv > 0. )
            q2Dprv = q2Dprv / areakv
            c2Dprv = c2Dprv / areakv
        end where
	call shympi_exchange_3d_node(q2Dprv)
	call shympi_exchange_3d_node(c2Dprv)
	!call shympi_barrier

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do ie=1,nel
          
	 lmax = ilhv(ie)
         areael = 12. * ev(10,ie)

         do l=1,lmax
         
           q2Dx=0.
           q2Dy=0.
           c2Dx=0.
           c2Dy=0.
                     
           do ii=1,3
             k=nen3v(ii,ie)
             b(ii)=ev(ii+3,ie)
             c(ii)=ev(ii+6,ie)
             q2Dx=q2Dx+(q2Dprv(l,k)*b(ii))
             q2Dy=q2Dy+(q2Dprv(l,k)*c(ii))
             c2Dx=c2Dx+(c2Dprv(l,k)*b(ii))
             c2Dy=c2Dy+(c2Dprv(l,k)*c(ii))
           end do
           vsclth = (q2Dx**2)+(q2Dy**2) + (c2Dx**2)+(c2Dy**2)
	   vsclth = (areael**1.5) * sqrt(vsclth)

	   difhv(l,ie) = vsclth

         end do
        end do
         
	return
        end 

!***********************************************************************
	subroutine green(ie,l,ugreen,vgreen)

! solves green identity for reynolds stresses

! ieltv:  >0 element  0: boundary  -1: open boundary

	use geom
	use hydro_admin
	use evgeom
	use levels
	use basin
        use fem_util

	implicit none

	integer ie
	integer l
	double precision ugreen,vgreen	!contribution to integral from green formula

	include 'param.h'

! shympi FIXME : will not work with mpi

	integer k,i,ii,iii,ienb,i1,i2
	double precision dl(3)
	double precision x(3),y(3)
	double precision u,v,unb,vnb
	double precision xm,ym,xmb,ymb
	double precision dist
	double precision ugrad,vgrad

	double precision x1,y1,x2,y2
	double precision distance
	distance(x1,y1,x2,y2) = sqrt((x1-x2)**2+(y1-y2)**2)

!------------------------------------------------------------
! get transports
!------------------------------------------------------------

	u = utlov(l,ie)
	v = vtlov(l,ie)

!------------------------------------------------------------
! compute geometric characteristics for element
!------------------------------------------------------------

!	---------------------------------------------
!	center point
!	---------------------------------------------

	do ii=1,3
	  k = nen3v(ii,ie)
	  x(ii) = xgv(k)
	  y(ii) = ygv(k)
	end do

	call baric(ie,xm,ym)

!	---------------------------------------------
!	length of sides
!	---------------------------------------------

	do ii=1,3
	  i1 = mod(ii,3) + 1
	  i2 = mod(i1,3) + 1
	  dl(ii) = distance(x(i1),y(i1),x(i2),y(i2))
	end do

!------------------------------------------------------------
! compute contribution
!------------------------------------------------------------

	do ii=1,3

!	  ------------------------------------------------------------
!	  find neibor element
!	  ------------------------------------------------------------

	  ienb = ieltv(ii,ie)
	  if( ienb .gt. 0 .and. ilhv(ienb) .lt. l ) ienb = 0

!	  ---------------------------------------------
!	  find velocity in neibor element
!	  ---------------------------------------------

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

!	  ---------------------------------------------
!	  find a distance for gradient
!	  ---------------------------------------------

	  if( ienb .gt. 0 ) then

!	    ------------------------------------------------
!	    compute distance to center point of neibor element
!	    ------------------------------------------------

	    call baric(ienb,xmb,ymb)
	    dist = distance(xm,ym,xmb,ymb)

	  else

!	    ------------------------------------------------
!	    compute virtual distance
!	    ------------------------------------------------

	    dist = 0.5 * distance(xm,ym,x(ii),y(ii))

	  end if

!	  ---------------------------------------------
!	  compute gradient
!	  ---------------------------------------------

	  ugrad = (unb - u) / dist
	  vgrad = (vnb - v) / dist

!	  ---------------------------------------------
!	  final contribution
!	  ---------------------------------------------

	  ugreen = ugreen + ugrad * dl(ii)
	  vgreen = vgreen + vgrad * dl(ii)

	end do

!------------------------------------------------------------
! end of routine
!------------------------------------------------------------

	end

!***********************************************************************

!---------------------------------------------------------------------
        end module diffusion_admin
!---------------------------------------------------------------------
