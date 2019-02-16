
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

c utilities for tracking (should be revised)
c
c revision log :
c
c 05.02.2009    ggu     copied from other files (lagrange.f)
c
c**********************************************************************

	subroutine retta(ext,cy,cx,b,el)
	
c dati gli estremi calcolo i coefficienti a, b della retta passante
	
	use basin

	implicit none
	
	integer ext(2),i,el
	integer p1,p2
	
	include 'param.h'
	
        real x1,y1,x2,y2
        real cy,cx,b ! parametri retta cyY=cxX+b
	
        real ax,ay,b1,b2,b3

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

c**********************************************************************

       subroutine dirz(pfin,ratio,ext,cy,cx,b,ny,nx,el,pb,ipb)
	
c per il calcolo della retta traiettoria percorsa del body:
c ho bisogno del punto di arrivo, pfin e un punto di partenza pini
c                         | questo punto lo calcolo considerando: 
c         /\ pfin         | differenza deltax tra le coord. x degli estremi
c        /  \             | della retta y=ax+b, applico il rapporto
c       /    \            | (ratio) dei segmenti a tale delta X e  
c      /      \           | trovo la frazione di x da aggiungere  
c     /  ratio \          | all''estremo opposto a cui si riferisce il 
c    <------------------> | ratio. Quindi trovo le nuove coordinate del punto
c	        pini      | pini. Da pini e pfin trovo il fascio di
c_________________________| rette in particolare an il coefficiente
c                           angolare an
	use basin

	implicit none
	
	include 'param.h'
	
	integer pb(2),ipb(2),pi
	
		
	integer pfin ! punto estremo del fascio di rette
        real ratio(2) ! frazioni sul segmento
	integer ext(2) ! puntatore estremi del segmento
        real cy,cx,b ! parametri retta del segmento
	integer p1,p2 ! estremi del segmento
        real an ! coefficient angolare fascio di rette
	integer el,i

        real x1,x2,y1,y2
        real dist,ax,ay,xi,yi
        real nx,ny
        real dfrq,aax,aay,axx,ayy
        real xn,yn,rr,a
	
c coordinate pfin

	pi=nen3v(pfin,el)
	xi=xgv(pi)
	yi=ygv(pi)

c definizione coordinate estremi del segmento
	
	p1=nen3v(ext(1),el)
	p2=nen3v(ext(2),el)
   	
	x1=xgv(p1)	
	x2=xgv(p2)
	y1=ygv(p1)
	y2=ygv(p2)

c calcola distanza tra 2 punti p1,p2

	ax=x2-x1
	ay=y2-y1
	axx=ax**2
	ayy=ay**2

	dist=sqrt(axx+ayy)

c trasformazione frazione percentuale di 2 in distanza da punto 1

	rr = 0.
	do i=1,2
	 if(ext(1).ne.ipb(i))then
	  dfrq=dist*(ratio(i))
	  rr=ratio(i)
	 endif
	enddo 
        
c calcolo coordinate punto nella retta distante dfrq da p1 e
c compreso tra p1 e p2

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

c calcolo coefficiente angolare del fascio di rette

 	ny=xn-xi
        nx=yn-yi
        if(ny.eq.0.)then
	 ny=0
	 nx=1
        endif 
	
        end

c**********************************************************************

	subroutine traj(ny,nx,x,y,tb)

	implicit none

c in questa subroutine viene definita la retta y=an+tb passante per il 
c punto x, y e di coefficiente angolare an

        real ny,nx,x,y,tb
        real anx,any

	anx=nx*x
	any=ny*y
	tb=any-anx
	
	end

c**********************************************************************

	subroutine interc(ny,nx,b,cy,cx,ib,x,y)
	
c calcolo delle coordinate x e y di intercetta tra la retta y=ax+b
c e la retta y=iax+ib

        implicit none

        real ny,nx,b,cx,cy,ib,x,y
        real a,aa,aaa,ab,aab,b1,b2,cc,ia

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

c**********************************************************************
	
	subroutine distp(x,y,ext,d,near,far,el)

c calcolo della distanza (dist) del punto x,y dall'estremo piu distante

	use basin

	implicit none

        real x,y ! punto da cui determinare la distanza
	integer ext(2) ! puntatore estremi del segmento
        real d ! distanza dal punto piu lontano
	integer near,far ! puntatore dell'estremo piu lontano e piu vicino
			 ! rispetto al punto x,y
	include 'param.h'

        integer p1,p2,el
        real x1,x2,y1,y2
        real ax,ay,axx,ayy
        real dxx1,dxx2
        real dist
	
c definizione coordinate estremi del segmento

        p1=nen3v(ext(1),el)
        p2=nen3v(ext(2),el)

        x1=xgv(p1)
        x2=xgv(p2)
        y1=ygv(p1)
        y2=ygv(p2)

c calcola distanza punto x1,y1 e punto x2,y2

        ax=x2-x1
        ay=y2-y1
        axx=ax**2
        ayy=ay**2

        dist=sqrt(axx+ayy)

c calcola distanza punto x,y e punto x1,y1

        ax=x-x1
        ay=y-y1
        axx=ax**2
        ayy=ay**2
        dxx1=sqrt(axx+ayy)
	
c calcola distanza punto x,y e punto x2,y2

	ax=x-x2
        ay=y-y2
        axx=ax**2
        ayy=ay**2
        dxx2=sqrt(axx+ayy)

c individuo estremo piu lontano (far) e piu vicino (near)

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
c        PRINT*,'STOP BODY IN X ',x,' Y ',y,' NON CONTENUTO IN EL ',el
         el=-el
         return
        endif
        end

c**********************************************************************

	subroutine pnt_inside(d,dd,ny,nx,x1,y1,x2,y2,xn,yn)

c calcolo coordinate di un punto sulla retta nyY=nxX+b distante da x1,y1
c di un valore d e compreso tra x1,y1 e x2,y2

	implicit none
	
        real d,a,x1,y1,x2,y2,dd
        real xn,yn,ny,nx
	
c calcolo coordinate punto nella retta distante d da x1 e
c compreso tra x1, y1 e x2, y2

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

c**********************************************************************

