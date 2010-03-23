c
c $Id: lagrange_in.f,v 1.11 2009-02-13 17:22:44 georg Exp $
c
c subroutines per inizializzare posizione flottanti in lagrange.f
c subroutines per scrittura output lagrange.f
c
c contents : subroutine posit           initial t,x,y,el bodies
c            subroutine set_output         write in output file
c            subroutine nbody              t,x,y,el bodies
c
c revision log :
c
c 00.00.2003    aac     written from scratch
c 07.11.2005    ggu     last version integrated in main tree
c 01.02.2006    ggu     use xyvar as distance to vary point if on node or side
c 29.11.2006    ggu     new version integrated into main model
c 22.06.2007    ggu     in posit use xydif to avoid endless loop
c 10.11.2007	ggu	renamed ran0 to ran8 (conflict, not used here)
c 13.06.2008	ggu	use insert_particle, bugfix for multiple lines
c
c old routines ... should be eliminated
c
c******************************************************************

        subroutine set_input

c lettura da file input dati su posizione di rilascio flottanti
c costruzione matrice lunghezza lati di ogni elemento

        implicit none

        include 'param.h'
        include 'lagrange.h'

        call filler ! selezione area e inizializzazione flottanti
        call rprs  !lunghezza lati prisma

        end

c****************************************************************

	subroutine posit(xbb,ybb,m)

c in questa subroutine viene trovato per ogni flottante
c l elemento di appartenenza 

	implicit none

	include 'param.h'
	include 'lagrange.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it
        
        real xgv(1),ygv(1) ! coordinate nodi
        integer nen3v(3,1) ! puntatore sui nodi di ogni elemento
        common /xgv/xgv, /ygv/ygv
        common /nen3v/nen3v

	real x,y,xd(3),yd(3),xbb(10*nbdydim),ybb(10*nbdydim)
	real xyvar,xyvareps,xydif
	integer p,ie,k,n,iin
        integer ifileo,ih,m
	real inline,in	

	xyvar = 0.01		!distance to vary point if on node or side
	xyvareps = 1.e-6	!rel. distance to vary point if on node or side

	xydif = 1. - xyvareps
	if( xydif .eq. 1. ) stop 'error stop posit: xyvareps to small'
	
	n=0
        do ih=1,m
	x=xbb(ih)
        y=ybb(ih)
565     continue        
	!write(6,*)  'new point: ',ih,x,y
        do ie=1,nel
         in=0
	 do k=1,3
          p=nen3v(k,ie)
	  xd(k)=xgv(p)
	  yd(k)=ygv(p)
          
          if(x.eq.xd(k) .and. y.eq.yd(k))then !is on node -> vary
           !x=x-xyvar
           !y=y-xyvar
	   x = x * xydif
	   y = y * xydif
           goto 565      
          endif
         end do

         in=inline(3,xd,yd,x,y)
         
         if(in.eq.-999)then !is on side -> vary
          !x=x-xyvar
          !y=y-xyvar
	  x = x * xydif
	  y = y * xydif
          goto 565
         endif
                 
         iin=int(in)
         
         if( iin .ne. 0 )then
	  n=n+1

	  call insert_particle(ie,x,y)

          write(78,*)'1 ',n,' 0 ',x,y
	  go to 988	
	 endif
  	end do
988     continue
        end do
        
	!nbdy=n-1		!why????????????????/
        print*,'TOTAL AMOUNT OF PARTICELS: ',nbdy,' out of ',m
        
	close(14)

	end	

c-----------------------------------------------------------------------------

	function inline(n,x,y,xbar,ybar)

c per ogni elemento ricavo la somma di tutti gli angoli che crea con i punti
c      appartenenti alla linea


        implicit none
        
	integer i,j,n
        real pi,xbar,ybar,x(n),y(n)
        real inline,alfa,ang,summ,x1,x2,y1,y2,x3,y3 

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
c***************************************************************************
        subroutine nbody(iml,it)

        include 'param.h'
        include 'lagrange.h'
                      
        integer i,iml,mbdy,it
        
        mbdy=nbdy
        nbdy=nbdy+iml
        do i=mbdy+1,nbdy
         x_body(i)=xst(i-mbdy)
         y_body(i)=yst(i-mbdy)
         ie_body(i)=est(i-mbdy)
         tin(i)=it
        end do     
        end

c***************************************************************************


        subroutine rnum(x1,x2,x3)

        real ri,rs,x1,x2,x3,r1,r2
        data ri/100./
        data rs/13./
        save rs,ri

        x1=mod(ri,rs)+0.1
        x2=mod(ri,x1)+0.1
        x3=mod(ri,x2)+0.1

        rm=x3
        r1=x1
        r2=x2
        x1=(1-x1/(rs))/2
        x2=(1-x2/(r1))/2
        x3=(1-x3/(r2))/2
        rs =rm
        end

c***************************************************************************

        function ran8(idum)

        integer idum,IA,IM,IQ,IR,MASK
        real ran0,AM
        parameter (IA=16807,IM=2147483647,AM=1./IM,
     *             IQ=127773,IR=2836,MASK=123459876)
        integer k
       

c Minimal random number generator of Park and Miller. Returns a
c uniform random deviate between 0.0 and 1.0. Set or 
c reset idum to any integer value (except the unlikely value MASK)
c to initialize the sequence; idum must not be altered between
c calls for successive deviates in a sequence.

        idum=ieor(idum,MASK)  ! XORing with MASK allows 
                              ! use of zero and other simple
        k=idum/IQ  !  bit patterns for  idum.
        idum=IA*(idum-k*IQ)-IR*k ! Compute without overflows by
        idum=mod(IA*idum,IM)
        if (idum.lt.0)idum=idum+IM  ! Schrage's method.   
        ran8=AM*idum  !   Convert idum to a floating result.
        idum=ieor(idum,MASK)   !  Unmask before return
        return
        end
                      
c****************************************************************

        subroutine filler

c puts points into an specified area
c with a regular distribution

	implicit none

	include 'param.h'
        include 'lagrange.h' 

        character*80 descrp
        common /descrp/ descrp

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real grav,fcor,dcor,dirn,rowass,roluft
        common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft

        real xgv(1), ygv(1)		!coordinate nodes
        common /xgv/xgv, /ygv/ygv

        integer nen3v(3,1)			!element index
        common /nen3v/nen3v

	real hdepth,xm,xn,ym,yn,enx
        integer nx,ny

        integer iapini,nl,io
        integer ie,ii,i,j,iii,ntel

        integer ntotmax
        parameter (ntotmax=10*nbdydim)

        integer ndim,m,n,num(ntotmax),nline,nty(ntotmax)
        parameter (ndim=500)
        real x(ndim), y(ndim),px(ntotmax),py(ntotmax)
        real pi,xbar,ybar,adef,c(10),s(10),a(10),summ

        integer ifileo
        real ttt,xbb(ntotmax),ybb(ntotmax)
        character*80 line,fluxt,lint
        integer lnlag,totnb

        real getpar

c-----------------------------------------------------------------
c lettura file di area e definizione parametri di flusso 
c---------------------------------------------------------------------------
         
        call getfnm('lagra',line)               !lettura file di area
        lnlag=ifileo(0,line,'form','old')
        
        totnb=getpar('nbdy')                    !lettura parametro di flusso

        nline = 0

        n=1 
	m=0
1       continue
        call rdline(ndim,n,x,y,xn,xm,yn,ym,lnlag)  !punti perimetro area
          if( n .le. 0 ) goto 2

          enx=sqrt(real(totnb)) !numero di elementi per lato
          nx=int(enx)
          ny=nx
        
c         ----------------------------------------------------------------
c         fill rectangle with regular points
c         ----------------------------------------------------------------
          
          call fill(ntotmax,xm,xn,ym,yn,nx,ny,ntel,px,py) !riempimento regolare
	  write(6,*) 'after fill'
          nline = nline + 1
	  write(6,*) 'filling line: ',nline,n
          
c         ----------------------------------------------------------------
c         deletes points that are outside of line
c         ----------------------------------------------------------------

          call eline(n,x,y,ntel,px,py,m,num,xbb,ybb) !riempimento area
	  write(6,*) 'after eline'
          
        goto 1
2       continue

c       ----------------------------------------------------------------
c       assigns element number to points
c       ----------------------------------------------------------------

	write(6,*) 'looking for element numbers'
        call posit(xbb,ybb,m)
	write(6,*) 'finished with element numbers'

        end


c----------------------------------------------------------------------------
c******************************************************************************
c      lettura nodi linee
c******************************************************************************88

c----------------------------------------------------------------------------
        subroutine rdline(ndim,n,x,y,xmin,xmax,ymin,ymax,lnlag)
c----------------------------------------------------------------------------
c  lettura file nodi linee e suddivisione degli stessi per appartenenza ad una
c   stessa linea.
c al main arrivano i valori di coord per ogni nodo
c--------------------------------------------------------------------------
        implicit none

        integer ndim
        real x(500), y(500),xmin,xmax,ymin,ymax
        integer isave,lnlag,n,iflag
        real xsave,ysave,xa,ya
        save xsave,ysave,isave

        data isave / 0 /

        n = 0

c acquisizione prima riga dal file se e il primo passaggio di lettura o
c da xsave ysave nel casio in cui si tratti di una linea successiva alla prima

        if( isave .eq. 0 ) then
          read(lnlag,*,end=2) xa,ya,iflag
          n = 1
          x(n) = xa
          y(n) = ya
          xmin=x(n)
          xmax=x(n)
          ymin=y(n)
          ymax=y(n)
        else
          n = 1
          x(n) = xsave
          y(n) = ysave
        end if
        if(x(n).le.xmin)then
         xmin=x(n)
        endif
        if(x(n).ge.xmax)then
         xmax=x(n)
        endif
        if(y(n).le.ymin)then
         ymin=y(n)
        end if
        if(y(n).ge.ymax)then
         ymax=y(n)
        endif
        isave = 0

10      continue

c lettura seconda riga ed assegnazione valori coord. ad x(n),y(n), fino a che   c  iflag e diverso da 1

          read(lnlag,*,end=2) xa,ya,iflag
          if( iflag .eq. 0 ) then
            n = n + 1
            if( n .gt. ndim ) stop 'eror stop rdline: ndim'

            x(n) = xa
            y(n) = ya
c nel caso iflag=1 salvo tali valori come appartenenti alla prima riga del      c successiva serie di nodi
            if(x(n).le.xmin)then
             xmin=x(n)
            end if
            if(x(n).ge.xmax)then
             xmax=x(n)
            endif
            if(y(n).le.ymin)then
             ymin=y(n)
            end if
            if(y(n).ge.ymax)then
             ymax=y(n)
            endif
          
          else
            isave = 1
            xsave = xa
            ysave = ya
            goto 2
          end if

          goto 10
2       continue

        end

c******************************************************************************

        subroutine fill(ntotmax,xm,xn,ym,yn,nx,ny,in,px,py)

        integer ntotmax
        real xm,xn,ym,yn,nnx,nny
        integer nx,ny
        real x,y,px(1),py(1)
        integer in
        
        nny=(ym-yn)/ny
        nnx=(xm-xn)/nx
        x=xn
        in=0

        do i=1,nx
         y=yn
         do j=1,ny
           in=in+1
           if( in .gt. ntotmax ) stop 'error stop fill: ntotmax'
           px(in)=x
           py(in)=y
           write(95,*)' 1 ',in,' 0 ',x,y    
           y=y+nny
         end do
         x=x+nnx
        end do

        end

c******************************************************************************

          subroutine eline(n,x,y,nel,px,py,m,num,xbb,ybb) 

c computes all elements inside a line

          implicit none
                   
          integer n             ! total number of points in line
          real x(1), y(1)       ! coordinates of points forming line
          integer nel
          real px(1),py(1)      ! coordinates of barycenter of elements
          integer m             ! total number of elements found
          integer num(1)        ! list of elements found
          real xbb(1),ybb(1)    ! barycenters of elements found

          integer ie,in,iin
          real xbar,ybar
          integer inlinef
          
          !m=0		!do not initialize, otherwise only one line is filled

          do ie=1,nel
           xbar=px(ie)
           ybar=py(ie)
           in = inlinef(n,x,y,xbar,ybar)
           if( in .ne .0 )then
            m=m+1
            xbb(m)=xbar
            ybb(m)=ybar   
            num(m)=ie
            write(94,*)' 1 ',m,' 0 ',xbar,ybar,iin
           end if
          end do

          end

c****************************************************************************8

	function inlinef(n,x,y,xbar,ybar)

c checks if point (xbar,ybar) is inside line given by (n,x,y)

        implicit none
        
        integer inlinef
        integer n
        real x(1),y(1)
        real xbar,ybar

	integer i,j
        real alfa,ang,summ
        real x1,x2,y1,y2,x3,y3 

        real pi,pi2
        parameter ( pi = 3.14159 , pi2 = 2*pi )

        summ=0

        do i=1,n
           j = mod(i,n) + 1
           x1=x(i)
           x3=x(j)
           y1=y(i)
           y3=y(j)
           x2=xbar
           y2=ybar 
           alfa=ang(x1,y1,x2,y2,x3,y3)
           summ=summ+alfa
        end do

	inlinef = nint(summ/pi2)

        end   

c***************************************************************************

        subroutine bar(ie,xbar,ybar,nen3v,xgv,ygv)

c  computes barycenter of element

        implicit none
       
        integer ie
        real xbar,ybar
	integer nen3v(3,1)
	real xgv(1),ygv(1)

        real xb,yb
        integer k,ii

        xb=0
        yb=0
        do ii=1,3
         k=nen3v(ii,ie)
         xb=xb+xgv(k)
         yb=yb+ygv(k)
        end do
        xbar=xb/3
        ybar=yb/3 

        end

c***************************************************************************
        
