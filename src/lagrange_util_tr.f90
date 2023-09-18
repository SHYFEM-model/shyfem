!
! $Id: lagrange_util_tr.f,v 1.1 2009-02-13 17:22:44 georg Exp $
!
! utilities for tracking (should be revised)
!
! revision log :
!
! 05.02.2009    ggu     copied from other files (lagrange.f)
!
!**********************************************************************
!-------------------------------------------------------------------------
        module lagrange_util_tr
!-------------------------------------------------------------------------
        contains
!-------------------------------------------------------------------------

	subroutine retta(ext,cy,cx,b,el)
	
! dati gli estremi calcolo i coefficienti a, b della retta passante
	
	use basin

	implicit none
	
	integer ext(2),i,el
	integer p1,p2
	
	include 'param.h'
	
        double precision x1,y1,x2,y2
        double precision cy,cx,b ! parametri retta cyY=cxX+b
	
        double precision ax,ay,b1,b2,b3

	p1=nen3v(ext(1),el)
	p2=nen3v(ext(2),el)
	
	x1=xgv(p1)
	x2=xgv(p2)
	y1=ygv(p1)
        y2=ygv(p2)	

        
        
	if(y2.eq.y1)then
	 cx=0
	 b=y2
	 cy=1
	 goto 332
	elseif(x2.eq.x1)then
	 cx=1
	 b=-x2
	 cy=0
	 goto 332
	endif	
	ax=x2-x1
	ay=y2-y1
	
	cy=1/ay
	cx=1/ax
        b1=y2/ay
        b2=x2/ax
	b=b1-b2
332	continue

	end

!**********************************************************************

       subroutine dirz(pfin,ratio,ext,cy,cx,b,ny,nx,el,pb,ipb)
	
! per il calcolo della retta traiettoria percorsa del body:
! ho bisogno del punto di arrivo, pfin e un punto di partenza pini
!                         | questo punto lo calcolo considerando: 
!         /\ pfin         | differenza deltax tra le coord. x degli estremi
!        /  \             | della retta y=ax+b, applico il rapporto
!       /    \            | (ratio) dei segmenti a tale delta X e  
!      /      \           | trovo la frazione di x da aggiungere  
!     /  ratio \          | all''estremo opposto a cui si riferisce il 
!    <------------------> | ratio. Quindi trovo le nuove coordinate del punto
!	        pini      | pini. Da pini e pfin trovo il fascio di
!_________________________| rette in particolare an il coefficiente
!                           angolare an
	use basin

	implicit none
	
	include 'param.h'
	
	integer pb(2),ipb(2),pi
	
		
	integer pfin ! punto estremo del fascio di rette
        double precision ratio(2) ! frazioni sul segmento
	integer ext(2) ! puntatore estremi del segmento
        double precision cy,cx,b ! parametri retta del segmento
	integer p1,p2 ! estremi del segmento
        double precision an ! coefficient angolare fascio di rette
	integer el,i

        double precision x1,x2,y1,y2
        double precision dist,ax,ay,xi,yi
        double precision nx,ny
        double precision dfrq,aax,aay,axx,ayy
        double precision xn,yn,rr,a
	
! coordinate pfin

	pi=nen3v(pfin,el)
	xi=xgv(pi)
	yi=ygv(pi)

! definizione coordinate estremi del segmento
	
	p1=nen3v(ext(1),el)
	p2=nen3v(ext(2),el)
   	
	x1=xgv(p1)	
	x2=xgv(p2)
	y1=ygv(p1)
	y2=ygv(p2)

! calcola distanza tra 2 punti p1,p2

	ax=x2-x1
	ay=y2-y1
	axx=ax**2
	ayy=ay**2

	dist=sqrt(axx+ayy)

! trasformazione frazione percentuale di 2 in distanza da punto 1

	rr = 0.
	do i=1,2
	 if(ext(1).ne.ipb(i))then
	  dfrq=dist*(ratio(i))
	  rr=ratio(i)
	 endif
	enddo 
        
! calcolo coordinate punto nella retta distante dfrq da p1 e
! compreso tra p1 e p2

	xn = 0.
	yn = 0.
        
	if(cy.eq.0)then
	 xn=x1
	 if(y2.gt.y1)then
	  yn=y1+(rr*abs(ay))
	 elseif(y1.ge.y2)then
	  yn=y1-(rr*abs(ay))
	 endif
	 goto 222
	endif
        
	a=cy/cx
	
	if(a.ge.0)then
	 if(x2.gt.x1)then
	  xn=x1+(rr*abs(ax))
          yn=y1+(rr*abs(ay))
	 elseif(x1.ge.x2)then
	  xn=x1-(rr*abs(ax))
          yn=y1-(rr*abs(ay))
	 endif
	elseif(a.lt.0)then
	 if(x2.gt.x1)then
	  xn=x1+(rr*abs(ax))
          yn=y1-(rr*abs(ay))
	 elseif(x1.ge.x2)then
	  xn=x1-(rr*abs(ax))
	  yn=y1+(rr*abs(ay))
	 endif
	end if
              
              
222	continue

! calcolo coefficiente angolare del fascio di rette

 	ny=xn-xi
        nx=yn-yi
        if(ny.eq.0.)then
	 ny=0
	 nx=1
        endif 
	
        end

!**********************************************************************

	subroutine traj(ny,nx,x,y,tb)

	implicit none

! in questa subroutine viene definita la retta y=an+tb passante per il 
! punto x, y e di coefficiente angolare an

        double precision ny,nx,x,y,tb
        double precision anx,any

	anx=nx*x
	any=ny*y
	tb=any-anx
	
	end

!**********************************************************************

	subroutine interc(ny,nx,b,cy,cx,ib,x,y)
	
! calcolo delle coordinate x e y di intercetta tra la retta y=ax+b
! e la retta y=iax+ib

        implicit none

        double precision ny,nx,b,cx,cy,ib,x,y
        double precision a,aa,aaa,ab,aab,b1,b2,cc,ia

	if(ny.eq.0)then
         x=-(b/nx)
         ia=cx/cy
         cc=ib/cy
         y=(ia*x)+cc
         goto 223
        elseif(cy.eq.0)then
	 x=-(ib/cx)
	 ia=nx/ny
         cc=b/ny
         y=(ia*x)+cc
         goto 223
	endif
	a=nx/ny
	ia=cx/cy
	b1=b/ny
	b2=ib/cy	
        aa=a-ia
        ab=b2-b1
        aab=a*ab
        x=ab/aa
        y=(aab/aa)+b1

223     continue
	
	end

!**********************************************************************
	
	subroutine distp(x,y,ext,d,near,far,el)

! calcolo della distanza (dist) del punto x,y dall'estremo piu distante

	use basin

	implicit none

        double precision x,y ! punto da cui determinare la distanza
	integer ext(2) ! puntatore estremi del segmento
        double precision d ! distanza dal punto piu lontano
	integer near,far ! puntatore dell'estremo piu lontano e piu vicino
			 ! rispetto al punto x,y
	include 'param.h'

        integer p1,p2,el
        double precision x1,x2,y1,y2
        double precision ax,ay,axx,ayy
        double precision dxx1,dxx2
        double precision dist
	
! definizione coordinate estremi del segmento

        p1=nen3v(ext(1),el)
        p2=nen3v(ext(2),el)

        x1=xgv(p1)
        x2=xgv(p2)
        y1=ygv(p1)
        y2=ygv(p2)

! calcola distanza punto x1,y1 e punto x2,y2

        ax=x2-x1
        ay=y2-y1
        axx=ax**2
        ayy=ay**2

        dist=sqrt(axx+ayy)

! calcola distanza punto x,y e punto x1,y1

        ax=x-x1
        ay=y-y1
        axx=ax**2
        ayy=ay**2
        dxx1=sqrt(axx+ayy)
	
! calcola distanza punto x,y e punto x2,y2

	ax=x-x2
        ay=y-y2
        axx=ax**2
        ayy=ay**2
        dxx2=sqrt(axx+ayy)

! individuo estremo piu lontano (far) e piu vicino (near)

	if(dxx2.ge.dxx1)then
	 far=ext(2)
 	 near=ext(1)
	 d=dxx2/dist
	elseif(dxx1.gt.dxx2)then
	 far=ext(1)
	 near=ext(2)
	 d=dxx1/dist
	endif 
        if((dxx1.gt.dist).or.(dxx2.gt.dist))then
!        PRINT*,'STOP BODY IN X ',x,' Y ',y,' NON CONTENUTO IN EL ',el
         el=-el
         return
        endif
        end

!**********************************************************************

	subroutine pnt_inside(d,dd,ny,nx,x1,y1,x2,y2,xn,yn)

! calcolo coordinate di un punto sulla retta nyY=nxX+b distante da x1,y1
! di un valore d e compreso tra x1,y1 e x2,y2

	implicit none
	
        double precision d,a,x1,y1,x2,y2,dd
        double precision xn,yn,ny,nx
	
! calcolo coordinate punto nella retta distante d da x1 e
! compreso tra x1, y1 e x2, y2

        if(ny.eq.0)then
	 xn=x2
	 if(y2.gt.y1)then
	  yn=y1+d	 
	 elseif(y2.lt.y1)then
	  yn=y1-d
	 end if
         goto 444 
	end if 
	a=nx/ny
        if(a.ge.0)then
         if(x2.gt.x1)then
          xn=((d/dd)*(abs(x2-x1)))+x1
          yn=((d/dd)*(abs(y2-y1)))+y1
         elseif(x1.ge.x2)then
          xn=x1-((d/dd)*(abs(x2-x1)))
          yn=y1-((d/dd)*(abs(y2-y1)))
         endif
        elseif(a.lt.0)then
         if(x2.gt.x1)then
	  xn=((d/dd)*(abs(x2-x1)))+x1
          yn=y1-((d/dd)*(abs(y2-y1)))
         elseif(x1.ge.x2)then
          xn=x1-((d/dd)*(abs(x2-x1)))
          yn=((d/dd)*(abs(y2-y1)))+y1
         endif
        end if
444     continue	
	end

!**********************************************************************

!-------------------------------------------------------------------------
        end module lagrange_util_tr
!-------------------------------------------------------------------------
