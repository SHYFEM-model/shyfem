c
c $Id: lagrange_in.f,v 1.1 2004/08/26 15:12:55 georg Exp $
c
c subroutines per inizializzare posizione flottanti in lagrange.f
c subroutines per scrittura output lagrange.f
c
c contents : subroutine posit           initial t,x,y,el bodies
c            subroutine set_output         write in output file
c
c******************************************************************
	subroutine posit

c iin questa sbroutine viene letto il fiile 
c di input della posizione e del tempo dei flottanti 
c dei flottanti. Quindi viene trovato per ogni flottante
c l elemento di appartenenza 


	implicit none

	include 'param.h'
	include 'lagrange.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        real xgv(1),ygv(1) ! coordinate nodi
        integer nen3v(3,1) ! puntatore sui nodi di ogni elemento
        common /xgv/xgv, /ygv/ygv
        common /nen3v/nen3v

	real x,y,xd(3),yd(3)
	integer p,ie,k,n,in

	integer inline	


c lettura file input 

	open(14,file='posizione.dat',status='unknown',form='formatted')

	
	n=0
987     continue
	read(14,*,end=990)x,y
565     continue        
        do ie=1,nel
         in=0
	 do k=1,3
          p=nen3v(k,ie)
	  xd(k)=xgv(p)
	  yd(k)=ygv(p)
          
          if(x.eq.xd(k) .and. y.eq.yd(k))then !se il body corrisponde
           x=x-10.                            !ad un nodo di el
           y=y-10.                            !vario di 0.1 m le coord. 
           goto 565      
          endif
         end do

         in=inline(3,xd,yd,x,y)
         
         if(in.eq.-999)then !se body sta su lato elemento
          x=x-10.            !vario di 0.1 m le coord. 
          y=y-10.           
          goto 565
         endif
                 
         if(in.ne.0)then
	  n=n+1

	  xst(n)=x
	  yst(n)=y
	  est(n)=ie

	  x_body(n)=x
	  y_body(n)=y
	  ie_body(n)=ie
	  write(78,*)x,y
c         write(78,*)'1 ',n,' 0 ',x,y
	  go to 988	
	 endif
  	end do
988     continue
	goto 987

990	continue 
	nbdy=n	
        print*,'START COMPUTING TRAJECTORIES'
        
	close(14)

	end	
c-----------------------------------------------------------------------------
	function inline(n,x,y,xbar,ybar)

c per ogni elemento ricavo la somma di tutti gli angoli che crea con i punti
c      appartenenti alla linea


        implicit none
        
	integer inline,i,j,n
        real pi,xbar,ybar,x(n),y(n),alfa,ang,summ,x1,x2,y1,y2,x3,y3 

 	pi = 3.14150

          summ=0
          do i=1,n
           if (i.eq.n)then
            j=1
           else
            j=i+1
           end if
           continue
           x1=x(i)
           x3=x(j)
           y1=y(i)
           y3=y(j)
           x2=xbar
           y2=ybar 
           alfa=ang(x1,y1,x2,y2,x3,y3)
           if(alfa.lt.3.15 .and. alfa.gt.3.13)then !check per lati
            inline=-999
            return
           endif
           summ=summ+alfa
          end do
          continue
           summ=summ
	  inline = summ
          return
          end   

c------------------------------------------------------------------------
        function ang(x1,y1,x2,y2,x3,y3)


c  determinazione angolo tra tre punti 1,bar,2

        implicit none
        real x1,x2,y1,y2,c(10),s(10)
        real a(10),pi,adef,ang,x3,y3  
 
	pi = 3.14150  

        call comp(x1,y1,x2,y2,x3,y3,c,s)
        call alfa1(c,s,a,pi) 
        call def(a,adef,pi) 
        ang=adef
        return
        end

c----------------------------------------------------------------------------
        subroutine comp(x1,y1,x2,y2,x3,y3,c,s)  

c componenti cartesiane c(x),s(y) dei due lati formanti langolo  

        real x1,x2,y1,y2
        real x3,y3, c(10),s(10)
         c(1)=x1-x2
         c(2)=x3-x2
         s(1)=y1-y2
         s(2)=y3-y2
        end
c-------------------------------------------------------------------------------
        subroutine alfa1(c,s,a,pi)

c  determinazione angoli tra vertici ed asse cartesiano di riferimento          
c           utilizzando lungezza vertici md, comp.xvertice c(x)
c                e funz. intrinseca arcsen  

        real a(10),md,pi
	real c(10),s(10)

c a(10) angoli tra vertici ed asse x
c md lunghezza vrertici          

        do ii=1,2
        md=sqrt(c(ii)**2+s(ii)**2)
        sn=abs(s(ii)/md)
        a(ii)=asin(sn)
        if((c(ii).le.0).and.(s(ii).gt.0))then
         a(ii)=pi-a(ii)
        else if((c(ii).le.0).and.(s(ii).le.0))then   
         a(ii)=pi+a(ii)
        else if((c(ii).gt.0).and.(s(ii).le.0))then   
         a(ii)=2*pi-a(ii)
        end if
        continue
        end do
        continue
        end         
c---------------------------------------------------------------------------

        subroutine def(a,adef,pi)

c ricavo angolo tra tre punti come differenza tra due angoli, adef=a(2)-a(1) 

        real adef,pi
	real a(10)
        
        adef=(a(2)-a(1))
        if(adef.gt.pi)then
         adef=-((2*pi)-adef)
        else if(adef.lt.(-pi))then
         adef=((2*pi)+adef)
        end if
        end  
c---------------------------------------------------------------------------  
          subroutine elist(nel,m,num,nty,in) 
            
          implicit none
          integer num(1000000),m,nty(100000),in
          integer iii,ie,nel

          do iii=1,m
           do ie=1,nel
            if (ie.eq.num(iii))then
             nty(in)=ie
	     in=in+1
            end if
            continue
           end do
          end do
          return 
          end
c***************************************************************************




        



