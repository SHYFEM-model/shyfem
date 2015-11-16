c
c $Id: lagrange_dif.f,v 1.6 2009-02-13 17:22:44 georg Exp $
c
c lagrangian diffusion
c
c revision log :
c
c 01.10.2005    aac     written from scratch
c 07.11.2005    ggu     last version integrated in main tree
c 19.06.2006    aac     bugs in diffusion corrected
c 20.06.2006    aac     diffusion coefficient from str file 
c 29.11.2006    ggu     new version integrated into main model
c 10.11.2007    ggu     renamed random routine from ran1 to ran9 (conflict)
c 14 09 2007    aac     new simplified version more efficient, ran9->1, ran9 del
c 24 10 2012    aac     do not kill particle if infinite loop -> leave it
c 25 01 2013    ggu     new formula for horizontal diffusivity
c
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
