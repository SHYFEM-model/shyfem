
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

! lagrangian diffusion
!
! revision log :
!
! 01.10.2005	aac	written from scratch
! 07.11.2005	ggu	last version integrated in main tree
! 19.06.2006	aac	bugs in diffusion corrected
! 20.06.2006	aac	diffusion coefficient from str file 
! 29.11.2006	ggu	new version integrated into main model
! 14.09.2007	aac	new simplified version more efficient, ran9->1, ran9 del
! 10.11.2007	ggu	renamed random routine from ran1 to ran9 (conflict)
! 23.03.2010	ggu	changed v6.1.1
! 24.01.2012	ggu	changed VERS_6_1_41
! 24.10.2012	aac	do not kill particle if infinite loop -> leave it
! 25.01.2013	ggu	new formula for horizontal diffusivity
! 19.12.2014	ggu	changed VERS_7_0_10
! 16.11.2015	ggu	changed VERS_7_3_14
! 25.05.2017	ccf	vertical diffusion introduced
! 25.10.2018	ggu	changed VERS_7_5_51
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	changed VERS_7_5_61
!
!************************************************************
!
! z = 0		top of layer
! z = 1		bottom of layer
!
!************************************************************

        subroutine diff_body(i,id,ty,x,y,z,xi,iel,lb,time)

! diffusion of one particle: internal coordinates for vertical
!                            external coordinates for horizontal

        use mod_lagrange
        use basin

        implicit none

        integer, intent(in)	:: i            !particle number
        integer, intent(in)	:: id           !particle id
        integer, intent(in)	:: ty           !particle type
        real, intent(inout)	:: x            !x-coordinate
        real, intent(inout)     :: y            !y-coordinate
        real, intent(inout)     :: z            !relative vertical position
        double precision, intent(inout):: xi(3) !internal coordinates
        integer, intent(inout)	:: iel          !element number
        integer, intent(inout)	:: lb           !layer 
        real, intent(in)	:: time         !time to diffuse

	integer			:: mode 	!type of diffusion scheme
	real 			:: dtime
        integer 		:: n
        real 			:: perc
        double precision 	:: xx,yy,zz

        integer, save 		:: nk = 0
        integer, save 		:: nl = 0
        integer, save 		:: nu = 0

        if( iel <= 0 .or. ty < 0) return     !particle out of domain or beached

!---------------------------------------------------------------
! initialize
!---------------------------------------------------------------

        n = 100         !maximum loop count
        zz = z
	dtime = time

!---------------------------------------------------------------
! vertical diffusion 
!---------------------------------------------------------------

        if( bvdiff ) then
	  mode = 3
          do while ( dtime > 0. .and. n > 0 )
            call lgr_diff_vert(mode,id,iel,lb,zz,dtime)
            if( lb < 1 ) exit
            n = n - 1                   !decrement loop counter
          end do
          z = zz 
        end if

!---------------------------------------------------------------
! horizontal diffusion 
!---------------------------------------------------------------

        if( bhdiff ) then
	  call lag_diff_hor(id,time,lb,iel,x,y)
	  xx = x
	  yy = y
	  call xy2xi(iel,xx,yy,xi)
	end if

!---------------------------------------------------------------
! special treatment after diffusion and finish up
!---------------------------------------------------------------

        if( bvdiff .and.  dtime > 0. ) then     !not finished diffusing
          if( n == 0 ) then
            nk = nk + 1
            perc = (100.*nk)/idbdy
            write(6,1000) 'warning particle dif',id,iel,n,dtime,perc
            !iel = -iel
          else if( iel < 1 ) then
            nl = nl + 1
            perc = (100.*nl)/idbdy
            write(6,1000) 'loosing particle dif',id,iel,n,dtime,perc
          else
            nu = nu + 1
            perc = (100.*nu)/idbdy
            write(6,1000) 'unknown error dif',id,iel,n,dtime,perc
          end if
        end if

        if( .not. bback ) then
          if( iel.gt.0 .and. iarv(iel).eq.artype ) iel = -iel
        end if

!---------------------------------------------------------------
! end of routine
!---------------------------------------------------------------

 1000   format(a,3i10,2f10.2)

        end subroutine diff_body

!**********************************************************************
! Vertical diffusion of body id
! Get vertical diffusion velocity
! Mode: 1 naive random walk
!       2 diffusive random walk (Visser 1997)
!       3 Milstein scheme (Grawe 2010)

	subroutine lgr_diff_vert(mode,id,iel,lb,z,time)

	use mod_lagrange

	implicit none

	integer, intent(in)		:: mode !type of diffusion scheme
	integer, intent(in)		:: id	!id of particle
	integer, intent(in)     	:: iel	!element
	integer, intent(inout)		:: lb	!level
	double precision, intent(inout) :: z	!rel vert pos: 0=top, 1=bottom
	real, intent(inout)		:: time	!travel time of particle

	real              		:: dt	!random walk time step
	double precision 		:: w	!vertical velocity
        double precision                :: wd   !vertical velocity from diffusivity
        double precision                :: hd   !layer thickness
	double precision 		:: dz	!vertical displacement
	double precision 		:: dv	!potential vertical dispacement
	double precision 		:: t,tv,tt
	real				:: rr
	logical 			:: bdebug
	logical 			:: bsurf,bbott
	integer 			:: lmax
	integer 			:: lborig

	blgrdebug = id == 23
	blgrdebug = id == 0
	bdebug = .false.
	bdebug = .true.
	bdebug = blgrdebug

	lborig = lb

        call diff_get_vert(mode,id,iel,lb,z,time,dt,lmax,hd,wd)
	w = wd					!vertical velocity due to diffusion

	!-----------------------------------------------
	! handle particles on surface or on bottom
	!-----------------------------------------------

	bsurf = lb == 1 .and. z == 0.d0		!particle on surface
	bbott = lb == lmax .and. z == 1.d0	!particle on bottom

        if( bsurf .and. w > 0.d0 ) w = -w       !reflection at surface
        if( bbott .and. w < 0.d0 ) w = -w       !reflection at bottom

	if ( w > 0.d0 .and. z == 0.d0 ) then
	  lb = lb - 1
	  z  = 1.d0
	end if
	if ( w < 0.d0 .and. z == 1.d0 ) then
	  lb = lb + 1
	  z  = 0.d0
	end if

	!-----------------------------------------------
	! determine time to leave layer (tv)
	!-----------------------------------------------

	tv = 2.*dt

	if( w /= 0. ) then
	  if( w > 0. ) then
	    dv = hd*z
	    tv = dv / w
	  else
	    dv = hd*(1.-z)
	    tv = dv / (-w)
	  end if
	else
	  dv = 0.
	end if

	!-----------------------------------------------
	! what happens first? (tt is total time available)
	!-----------------------------------------------

	tt = time
	t = dt
	if( tv > 0. ) t = min(t,tv)

	if( bdebug ) then
	  write(6,*) 'times for diffusion'
	  write(6,*) tv,tt
	end if

	!-----------------------------------------------
	! handle vertical diffusion
	!-----------------------------------------------

	if( tv == 0. ) then		!body not moving
	  !
	else if( tv > t ) then		!body remains in layer
	  dz = w*t/hd
	  z = z - dz
	else
	  if ( lmax == 1 ) then		!in case of 1 layer z = random(0,1)
            call random_number(rr)
	    z = rr
	    t = tt
	  else
 	    if( w > 0. ) then
	      if( lb > 1 ) then
	        lb = lb - 1
	        z = 1.
	      else
	        z = 0.
	      end if
	    else
	      if( lb < lmax ) then
	        lb = lb + 1
	        z = 0.
	      else
	        z = 1.
	      end if
	    end if
	  end if
	end if

	if( bdebug ) then
	  write(6,*) 'after vertical diffusion'
	  write(6,*) iel
	  write(6,*) lborig,lb,lmax
	  write(6,*) w,z
	end if

	!-----------------------------------------------
	! compute remaining time and check for error
	!-----------------------------------------------

	time = tt - t

	if( bdebug ) then
	  write(6,*) 'final time: ',time
	  write(6,*) iel,lb,z
	  write(6,*) 'lgr_diff_vert end debugging'
	end if

	if( z < 0. .or. z > 1. ) goto 98

	!-----------------------------------------------
	! end of routine
	!-----------------------------------------------

	return 
   98	continue
	write(6,*) 'id,iel ',id,iel
	write(6,*) 'tv,t ',tv,t
	write(6,*) 'dz,w,hd: ',dz,w,hd
	write(6,*) 'z,lb,lmax ',z,lb,lmax
	stop 'error stop lgr_diff_vert: internal error (2)'
   99	continue
	stop 'error stop lgr_diff_vert: internal error (1)'

	end subroutine lgr_diff_vert

!************************************************************
! Compute vertical velocity due to diffusivity 
! Mode: 1 naive random walk
!       2 diffusive random walk (Visser 1997)
!       3 Milstein scheme (Grawe 2010)

        subroutine diff_get_vert(mode,id,iel,lb,z,time,dt,lmax,hd,wd)

        use levels, only : nlvdi,nlv
        use basin
        use mod_lagrange

        implicit none

	integer, intent(in)		:: mode !type of diffusion scheme
        integer, intent(in)		:: id   !id of particle
        integer, intent(in)		:: iel  !element number
        integer, intent(in)		:: lb   !layer
	double precision, intent(in)    :: z	!rel vert pos: 0=top, 1=bottom
        real, intent(in) 	        :: time	!time step 
	real, intent(out)		:: dt	!random walk time step
        integer, intent(out)		:: lmax !maximum layers in element
        double precision, intent(out)	:: hd   !layer thickness
        double precision, intent(out)  	:: wd  	!vertical velocity due to diffusion

        real            	:: rnd		!random number
        real			:: hl(nlv)      !layer thickness
        real			:: hle(0:nlv)   !layer structures
        real 			:: htot         !total depth without zeta
        real                    :: htotz        !total depth with zeta 
	real			:: zp		!depth at which estimate diffusivity
        real                  	:: wdr(0:nlv)  	!vertical diffusivity [real]
	real			:: wdzp		!vertical diffusivity at zp
	double precision	:: dwd		!first derivative of diffusivity 
        integer 		:: iact         !element closest to xe
        real 			:: exxpp        !extrapolated values
	real			:: rs		!random number per sigma

        integer			:: ii,k,l

	!-----------------------------------------------
	! get layer structures and depths
	!-----------------------------------------------
        lmax = nlv
        call lagr_layer_thickness(iel,lmax,hl,htot,htotz)

	hle(0) = 0.
	do l = 1,lmax
	  hle(l) = hle(l-1) + hl(l)
	end do

        if( lb > lmax ) goto 99

	!-----------------------------------------------
	! get derivative of diffusivity on layer lb
	!-----------------------------------------------
        hd = hl(lb)
	dwd = (wde(lb,iel) - wde(lb-1,iel)) / hd

	!-----------------------------------------------
	! Get random walk time step 
	!-----------------------------------------------
        dt = min(dtvd(iel),time)

	!-----------------------------------------------
	! estimate diffusivity at z (+ 0.5*dwd*dt)
	!-----------------------------------------------
	zp = hle(lb-1) + z*hd
	if ( mode == 2 ) zp = zp + 0.5*dwd*dt
	zp = min(zp,htotz)
	zp = max(zp,0.)

	wdr = wde(0:lmax,iel)
	iact = 0
	wdzp = exxpp(2,lmax,hle,wdr,zp,iact)

	!-----------------------------------------------
	! compute vertical velocity due to diffusion
	! vertical velocity is > 0 if directed upward
	!-----------------------------------------------
        call random_number(rnd)
        rnd = (rnd-0.5)/0.5		!normal distribution [-1, 1]
	wd = rnd*sqrt(6.*wdzp/dt)
	if ( mode == 2 ) wd = wd - dwd

	!-----------------------------------------------
	! Milstein scheme 
	!-----------------------------------------------
	if ( mode == 3 ) then
	  call rgauss(1., rnd)		!gaussian distribution
	  rs = rnd*sqrt(dt)
	  wd = -0.5*dwd*(rs*rs)/dt - 0.5*dwd + rs*sqrt(2.*wdzp)/dt
	end if

        return
   99   continue
        write(6,*) 'id,iel,l,lmax: ',id,iel,lb,lmax
        stop 'error stop diff_get_vert: no such layer'

        end subroutine diff_get_vert


!******************************************************************
! Compute vertical diffusivity in element and random walk time step

        subroutine lag_vdiff_ele

        use basin
        use mod_lagrange
        use mod_diff_visc_fric, only : difv
        use levels

	implicit none

        integer                 :: ie,ii,k,l
	real			:: time		!simulation time step
        real                    :: hl(nlv)      !layer thickness
        real                    :: hle(0:nlv)   !layer structures
        real                    :: htot         !total depth without zeta
        real                    :: htotz        !total depth with zeta 
        real                    :: wdl
        double precision        :: dwd2(nlv)    !second derivative of diffusivity 
        real                    :: dmax
	real			:: dtd		!random walk time step
        real			:: hd		!layer thickness
        integer 		:: lmax 	!maximum layers in element
	real, parameter 	:: dtmin = 2.0	!minimum time step
	real, parameter 	:: difmol = 1.0e-06

	call get_timestep(time)

        do ie = 1,nel

	  !-----------------------------------------------
	  ! get layer structures and depths
	  !-----------------------------------------------
          lmax = nlv
          call lagr_layer_thickness(ie,lmax,hl,htot,htotz)
  
  	  hle(0) = 0.
  	  do l = 1,lmax
	    hle(l) = hle(l-1) + hl(l)
	  end do

	  !-----------------------------------------------
	  ! get diffusivity on elements
	  !-----------------------------------------------
	  do l = 0,lmax
	    wdl = difmol
            do ii=1,3
              k = nen3v(ii,ie)
              wdl = wdl + difv(l,k)
	    end do
            wde(l,ie) = wdl / 3.
          end do

	  !-----------------------------------------------
	  ! Compute random walk time step 1/(100*dwd2)
	  !-----------------------------------------------
	  dtd = time
	  if ( lmax > 1 ) then
 	    do l = 1,lmax-1
              hd = hl(l)
	      dwd2(l) = (wde(l+1,ie) + 2.*wde(l,ie) + wde(l-1,ie))/(hd*hd)
	    end do
            dwd2(lmax) = dwd2(lmax-1)
            if ( lmax > 2 ) dwd2(lmax) = dwd2(lmax-1) + 
     +  		(dwd2(lmax-1)-dwd2(lmax-2))/hl(lmax-1)*hl(lmax)
  
   	    dmax = maxval(abs(dwd2))
  	    if ( dmax /= 0. ) dtd = 1./(100.*dmax)
            dtd = max(dtmin, dtd)
	  end if
	  dtvd(ie) = dtd
	end do

	end subroutine lag_vdiff_ele

!******************************************************************
! Horizontal diffusion of body id

        subroutine lag_diff_hor(id,time,lb,iel,x,y)

	use mod_lagrange
        use mod_diff_visc_fric, only : difhv

        implicit none
        
        integer, intent(in)             :: id   !id of particle
        real, intent(in)		:: time !time to diffuse
	integer, intent(in)		:: lb	!level
        integer, intent(inout)          :: iel  !element
        real, intent(inout)		:: x    !x-coordinate
        real, intent(inout)		:: y    !y-coordinate

        real 				:: k     !diffusion coefficient 
        real				:: dx,dy !random displacement
	real				:: xnew,ynew
        integer				:: ieold,ienew
        logical				:: track_xi_has_layer
        integer, save			:: icall = 0	!initialization parameter
        logical, save 			:: bsphe
        integer       			:: isphe
        double precision 		:: ddx,ddy,xc,yc
        double precision 		:: dlat0,dlon0  !center of projection

	if( icall == 0 ) then
          call get_coords_ev(isphe)
          bsphe = isphe .eq. 1
	  icall = 1
	end if

        if( iel <= 0 ) return     !particle out of domain

        !-----------------------------------------------
        ! Type of horizontal diffusion
        ! Mode: 1 Parameter rwhpar
        !       2 Function of type of diffusion (idhtyp)
        !-----------------------------------------------

	if ( rwhpar >= 0. ) then
            k = rwhpar  
	else
            k = -rwhpar*difhv(lb,iel)
        end if

        !-----------------------------------------------
        ! Initialize values
        !-----------------------------------------------
	dx = 0.
	dy = 0.
        ieold = iel

	ienew = 0

        !-----------------------------------------------
        ! Compute horizontal diffusion
        !-----------------------------------------------
	call lag_get_dxdy(time,k,dx,dy)
        xnew = x + dx
        ynew = y + dy

        if ( bsphe ) then
           ddx = dx 
	   ddy = dy
           dlon0 = x
           dlat0 = y
	   call ev_c2g(ddx,ddy,xc,yc,dlon0,dlat0)
           xnew = xc
           ynew = yc
	end if

        call find_elem_from_old(ieold,xnew,ynew,ienew)

	if (ienew == 0 ) then
	   ienew = ieold
	   call particle_on_side(ienew,x,y,xnew,ynew)
	end if

        !-----------------------------------------------
        ! Return new coordinates and element
        !-----------------------------------------------
        if( track_xi_has_layer(abs(ienew),lb) ) then
          x   = xnew
          y   = ynew
          iel = ienew
	end if
	 
        end subroutine lag_diff_hor

!******************************************************************
! Randow walk with rnd having uniform distribution between -1 and 1

        subroutine lag_get_dxdy(time,k,dx,dy)

	use mod_lagrange

        implicit none
        
        real, intent(in)		:: time !time to diffuse
        real, intent(in)		:: k    !diffusion coefficient
        real, intent(out)		:: dx   !x-displacemnet
        real, intent(out)		:: dy   !y-displacement

	real				:: rnd,cf

        cf = sqrt(6.*k*time)

        call random_number(rnd)
        rnd = (rnd-0.5)/0.5
        dx = cf*rnd

        call random_number(rnd)
        rnd = (rnd-0.5)/0.5
        dy = cf*rnd

	end subroutine lag_get_dxdy

!******************************************************************
!Generate random distribution with 0 mean and standard deviation SIGMA.
!Returns one (two) gaussian random number from the same distribution. 

        subroutine rgauss(sigma, y1)

        implicit none

        real, intent(in)        :: sigma   !standard deviation of distribution
        real, intent(out)       :: y1      !random number
        real                    :: y2      !second random number

        real                    :: rnd
        real*8                  :: x1, x2, w

        x1 = 0.
        x2 = 0.
	w  = 0.

        do while ( (w .ge. 1.0) .or. (w .eq. 0.0) )
          call random_number(rnd)
          x1 = 2.0 * rnd - 1.0
          call random_number(rnd)
          x2 = 2.0 * rnd - 1.0
          w = x1 * x1 + x2 * x2
        end do

	w = sigma*sqrt( (-2.0 * log( w ) ) / w )
	y1 = x1 * w
	y2 = x2 * w

	end subroutine rgauss

c*****************************************************************
c******************************************************************
c	OLD ROUTINES
c******************************************************************
c******************************************************************

        subroutine lag_diff(ie,bdy,x,y)

	use mod_lagrange

        implicit none
        
        include 'param.h'
        
	include 'femtime.h'
        real dt,ttime
        
        real k ! coefficiente di diffusione        
        
        integer ie,nlp ! numero elemento in cui si trova il body
        real*8 dx,dy ! spostamento random
        real x,y ! posizione body
        integer bdy ! numero body
      
        real cy,cx,b ! coefficienti traiettoria  
       
        real xold,yold,xnew,ynew,lt_ex
        integer ieold,ienew
	integer icount

        k=rwhpar  

	dx= 0
	dy= 0

	if( ie .le. 0 ) stop 'error stop lag_diff: internal error'

        xold=x
        yold=y
        ieold=ie
        call get_timestep(dt)
        ttime=dt

	ienew = 0
	icount = 20
	do while ( ienew .eq. 0 )

          call lag_rand(k,ttime,dx,dy) !calcolo spostamento 
          xnew=xold+dx
          ynew=yold+dy
          call find_elem_from_old(ieold,xnew,ynew,ienew)
	  icount = icount - 1
	  !if( icount .eq. 0 ) ienew = -ieold
	  if( icount .eq. 0 ) return	!return without diffusion

	end do

        x=xnew
        y=ynew
        ie=ienew
	 
        end 
c*******************************************************************
       
        subroutine lag_rand(k,ttime,dx,dy)

        implicit none

        include 'param.h'

        real ttime
        real k
        real*8 dx,dy
        real a,b,cf
        real nmb 
        real gasdev
        integer idum
        data idum/387/
        save idum
        cf=sqrt(2*k*ttime) 
	cf=sqrt(6*k*ttime)  !new formula 3d
        
        dx=cf*gasdev(idum)
        dy=cf*gasdev(idum)

        end

c*********************************************************************

        FUNCTION gasdev(idum)

        INTEGER idum
        REAL gasdev
c Returns a normally distributed deviate with zero mean and unit variance, using
C ran1(idum)
c as the source of uniform deviates.
        INTEGER iset
        REAL fac,gset,rsq,v1,v2,ran1
        SAVE iset,gset
        DATA iset/0/

        if (idum.lt.0) iset=0
        if (iset.eq.0) then
1        v1=2.*ran1(idum)-1.
         v2=2.*ran1(idum)-1.
         rsq=v1**2+v2**2
         if(rsq.ge.1..or.rsq.eq.0.)goto 1
         fac=sqrt(-2.*log(rsq)/rsq)
         gset=v1*fac
         gasdev=v2*fac
         iset=1
        else
         gasdev=gset
         iset=0
        endif
        return
        END

c**************************************************************************
