c
c $Id: lagrange.f,v 1.1 2004/08/26 15:12:55 georg Exp $
c
c subroutines for computing lagrangian trajectories
c
c contents :	main_lagrange
c	 subroutine set_input           initial t,x,y bodies
c
c        subroutine setup_fluxes        initializes fx 
c         subroutine getaz	
c	  do  				loop on nodes
c	   function flxtype
c	   subroutine pntfla
c          subroutine mk_rflux           flux through volume k
c          subroutine mk_tflux           flux through vertexes
c          subroutine setup_fx           set up fx(3,neldim)
c	  end do
c       
c	 subroutine drogue(it)          compute trajectories 
c	 do 				loop on drogues
c	  subroutine track_body    
c	   subroutine track_orig
c	   subroutine track_line
c	 end do
c	  subroutine set_output   	write in output file       
c            
c******************************************************************

	subroutine lagrange


	implicit none

        include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

	integer icl,itstart
	data icl/0/
	save icl,itstart 

	if(icl.eq.0)then
	 call set_input
	 itstart=it
	 icl=1
	endif

	call setup_fluxes
	
 	call drogue(idt)
	
  	call set_output(itstart,it)

	end

c****************************************************************

	subroutine set_input
	 
	implicit none

c lettura da file input dati su posizione di rilascio flottanti	
c costruzione matrice lunghezza lati di ogni elemento	
	
 	call posit !inizializzazione flottanti
	
	call rprs  !lunghezza lati prisma

	end

c****************************************************************
	subroutine rprs
	
	implicit none
	
	include 'param.h'
	include 'lagrange.h'

c costruzione matrice lunghezza lati elementi vert(ii,ie)
	
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	
        real xgv(1), ygv(1)
        common /xgv/xgv, /ygv/ygv

        integer nen3v(3,1)
        common /nen3v/nen3v

	real x1,x2,y1,y2,dx,dy,ddx,ddy,d
	integer ie,k,i1,i2,p1,p2

	do ie=1,nel
	 do k=1,3
          i1=mod(k,3)+1
          i2=mod(i1,3)+1
	  p1=nen3v(i1,ie)
	  p2=nen3v(i2,ie)
	  x1=xgv(p1)
	  x2=xgv(p2)
	  y1=ygv(p1)
	  y2=ygv(p2)
	  dx=x2-x1
	  dy=y2-y1
	  ddx=dx**2
	  ddy=dy**2
	  d=sqrt(ddx+ddy)
	  dvert(k,ie)=d
	 end do
	end do

	end

c******************************************************************
        subroutine set_output(itstart,it)

        implicit none

	include 'param.h'
	include 'lagrange.h'

	integer itstart      ! tempo del rilascio
	real x,y,rest,nxp,nno,nnbdy
	integer i,it,mn,f(nbdydim),ie,no,ls,inxp
	save mn,f,no
	data mn,no/0,0/
	data f/nbdydim*0/
        data ls/0/
        save ls
c scrittura in file output dati su posizione e tempo flottanti
c rilasciati

	do i=1,nbdy
         ie=ie_body(i)
         if(ie.lt.0)goto 323
	 x=x_body(i)
	 y=y_body(i)
         if(it-ls.ge.300)then
         write(77,*)x,y
         endif
	 if(i.eq.1491)then
	  mn=mn+1
	  write(82,*)' 1 ',mn,' 0 ',x,y
	 end if
	 goto 334
323      continue
	 if(f(i).eq.0)then
	  no=no+1
          rest=it-itstart
          write(76,*)xst(i),yst(i),rest
	  f(i)=1
	 endif
334      continue
	end do
        if(it-ls.ge.300)then
	nno=real(no)
	nnbdy=real(nbdy)
        nxp=((nno/nnbdy)*100)
        write(77,*)' 0.',' 0.'
	write(77,*)nxp
	write(77,*)it-itstart
cprint*,nxp
        ls=it
        endif
	
	if(no.eq.nbdy)stop

        end

c*****************************************************************

	subroutine setup_fluxes

	implicit none

	include 'param.h'
	include 'lagrange.h'
	
	integer ndim
	parameter (ndim=20)

c Inizializza il campo di moto da fornire come input alla 
c subroutine drogue per il calcolo delle traiettorie	

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	integer k,n,ipf,ipl,i,ii,ik
	integer itype,flxtype
	
	real az
	
	real rflux(ndim)       !flussi attraverso il volume finito k
	real tflux(ndim)       !flussi attraverso i vertici di k
	
        integer el(3)

        data el/4512,4505,4506/

        do ii=1,nel
         do i=1,3
         fx(i,ii)=0
         end do
        end do


	call getaz(az)
        do k=1,nkn
         itype=flxtype(k)
	 call pntfla(k,ipf,ipl)
	 n = ipl-ipf+1
         if( itype .gt. 1 ) n = n + 1   !boundary
         if( n .gt. ndim ) stop 'error stop flxnod: ndim'
	
    	 call mk_rflux(k,n,itype,az,rflux,ipf,ipl)
  	 call mk_tflux(k,n,itype,rflux,tflux)
  	 call setup_fx(k,n,tflux,ipf,ipl)

        end do

        call setup_vl

        end

c******************************************************************

	subroutine mk_rflux(k,n,itype,az,transp,ipf,ipl)

c computes flux through finite volume k
c internal section is defined by:  kbefor - k - kafter
c transp(i) = trasporto dovuto all'elemento i sul volume
c finito k
	
	implicit none

	integer k,itype,ti
	real az

	integer ndim		!must be at least ngr
	parameter (ndim=20)

        integer el(ndim)
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

        real ev(13,1)
        common /ev/ev
        real uov(1),vov(1),unv(1),vnv(1)
        common /uov/uov, /vov/vov, /unv/unv, /vnv/vnv
        real zenv(3,1), zeov(3,1)
        common /zenv/zenv, /zeov/zeov
        integer lenkv(1)
        common /lenkv/lenkv

	integer i,ip,ie,ii,n
	integer ipf,ipl
	real aj,area,dz,uv,rdt
	real b,c
	real azt,tt

	real transp(ndim)
	integer ithis
	 
	rdt = 1./idt
	azt = 1. - az

c compute transports into finite volume of node k -> transp
c computed transports are divergence corrected

	i = 0
	do ip=ipf,ipl
 	 i = i + 1
	 if( i .gt. ndim ) stop 'internal error flxnov: (3)'
	 ie = lenkv(ip)
	 el(i)=ie
         if( ie .le. 0 ) stop 'internal error flxnov: lenkv'
	 ii = ithis(k,ie)
	 aj = ev(10,ie)
	 area = 4. * aj
	 dz = zenv(ii,ie) - zeov(ii,ie)
	 b = ev(3+ii,ie)
	 c = ev(6+ii,ie)
	 uv = az * ( unv(ie) * b + vnv(ie) * c )
	 uv = uv + azt * ( uov(ie) * b + vov(ie) * c )
	 uv = 12. * aj * uv
	 uv = uv - dz * area * rdt
	 transp(i) = uv
	end do
	end

c******************************************************************

        subroutine mk_tflux(k,n,itype,rflux,tflux)

c computes fluxes over sides (tflux) from fluxes into node (rflux)

        implicit none

        integer k               !node
        integer n               !number of sides (tfluxes)
        integer itype           !type of node (1=int,2=bnd,3=BOO,4=OOB,5=OOO)
        real rflux(n)        !fluxes into node (element)
        real tflux(n)        !fluxes through sides (return value)

	integer ipf,ipl      !puntatori vertici

        integer i
        real rr

	if( itype .eq. 1 ) then         !internal node
                rr = 0.
                do i=1,n-1
                  rr = rr + i * rflux(i)
                end do
                rr = rr / n

                tflux(n) = rr
                do i=n-1,1,-1
                  tflux(i) = tflux(i+1) - rflux(i)
                end do
        else if( itype .eq. 2 ) then    !node on material boundary
                tflux(1) = 0.
                do i=2,n-1
                  tflux(i) = tflux(i-1) + rflux(i-1)
                end do
                tflux(n) = 0.
        else if( itype .eq. 3 ) then    !BOO - boundary on left
                tflux(n) = 0.
                do i=n-1,1,-1
                  tflux(i) = tflux(i+1) - rflux(i)
                end do
        else if( itype .eq. 4 ) then    !OOB - boundary on right
                tflux(1) = 0.
                do i=2,n
                  tflux(i) = tflux(i-1) + rflux(i-1)
                end do
        else if( itype .eq. 5 ) then    !node totaly on open boundary
                rr = 0.
                do i=1,n-1
                  rr = rr + i * rflux(i)
                end do
                rr = rr / n
                tflux(n) = rr
                do i=n-1,1,-1
                  tflux(i) = tflux(i+1) - rflux(i)
                end do

        else
                stop 'error stop make_fluxes: internal error (1)'
        end if
	end

c**********************************************************************
        subroutine setup_fx(k,n,tflux,ipf,ipl)


        implicit none

        include 'param.h'
	include 'lagrange.h'

        integer ithis
        integer ibhnd
        integer inext

        integer lenkv(1)
        common /lenkv/lenkv


        integer n,i,ie,ip,ii,k
        integer ipf,ipl,j,n1,n2

        real tflux(n)
	
        i=0
        do ip=ipf,ipl
         i=i+1
         ie=lenkv(ip)
         n1=ibhnd(k,ie)
         n2=inext(k,ie)
         j=i+1
         if(i.eq.n)then
          j=1
         endif
         fx(n2,ie)=fx(n2,ie)-tflux(j)
         fx(n1,ie)=fx(n1,ie)+tflux(i)      
        enddo
        end
c*********************************************************************
        subroutine setup_vl
	
	implicit none

	include 'param.h'
	include 'lagrange.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        real zenv(3,1),zeov(3,1)
        common /zenv/zenv, /zeov/zeov

        real hev(neldim)
        common /hev/hev

	integer ie,k,i1,i2
	real flx,dz,zi1,zi2,dp,ar,dst
	

	do ie=1,nel
	 dz=hev(ie)
         do k=1,3
	  dst=dvert(k,ie)
          i1=mod(k,3)+1
          i2=mod(i1,3)+1
	  zi1=zenv(i1,ie)
	  zi2=zenv(i2,ie)
	  dp=dz+((zi1+zi2)/2)
	  ar=dp*dst	
	  flx=fx(k,ie)          ! flusso su ie
	  vel_ie(k,ie)=flx/ar   ! velocità su ie
         end do
	end do
	
	end


c**********************************************************************

	subroutine drogue(idt)

	implicit none
	
	include 'param.h'
	include 'lagrange.h'

	integer idt
	integer i ! numeratore particelle rilasciate

c in questa subroutine per ogni body (i)
c viene calcolato il percorso effettuato durante un time
c step.
c Il prg si suddivide 2 diversi moduli, il primo in cui dato 
c un punto all'interno di un elemento si calcola il percorso e il tempo
c impiegato per arrivare (se si arriva in un timestep) al lato del
c elemento che contiene il body.
c A questo punto in un secondo viene calcolato il percorso a partire 
c dalla linea.  

	do i=1,nbdy
         call track_body(i,idt)
	end do
	end	

c####################################################################

        subroutine track_body(i,idt)

c in questa subroutine per ogni body (i) calcolo le coordinate x,y, 
c il numero dell'elemento iel che contiene il body, e il tempo (ttime) su 
c cui calcolare il percorso del body. 2 subroutines

c TRACK_ORIG se il body si trova all'interno di un elemento, all'inizio
c di ogni time step

c TRACK_LINE se il body si trova su una linea, se il body in uno stesso
c time step attraversa piu elementi

	implicit none

	include 'param.h'
	include 'lagrange.h'
	
	integer idt,i
	
	integer ttime ! tempo per calcolo del percorso
	real x ! coordinata x particella iesima
	real y ! coordinata y particella iesima
	integer iel ! numero elemento contenente particella iesima

	ttime=idt
	x=x_body(i)
	y=y_body(i)
cprint*,'x1 ',x,'y1 ',y
	iel=ie_body(i)
cprint*,'iel 1 ',iel
	if(iel.le.0)return
        call track_orig(ttime,i,iel,x,y)
	if(iel.le.0)goto 656  ! body i uscito
	do while (ttime.gt.0)
	 call track_line(ttime,i,iel,x,y)
	 if(iel.le.0)goto 656	! body i uscito 
        end do
c       print*,'x ',x,'y ',y
	x_body(i)=x
	y_body(i)=y
656     continue
	ie_body(i)=iel
cprint*,'iel ', iel
	end

c####################################################################
        subroutine track_orig(it,bdy,ie,xbdy,ybdy)

c in questa subroutine si calcola il percorso del body partendo
c da un punto all'interno dell'elemento ie 
c fino al successivo punto interno allo stesso elemento
c se entro il time step il body esce dall'elemento 
c si passa alla successiva subroutine con ttime=it

        implicit none
        
	include 'param.h' 
	include 'lagrange.h'
 
	integer ie,it,bdy ! bdy e' il numero del body
 
        real deltat ! frazione di time step body si muove in ie
 
	integer ieltv(3,1)
        common /ieltv/ieltv !punta ad elemento opposto al nodo 

	real xgv(1),ygv(1) ! coordinate nodi
        integer nen3v(3,1) ! puntatore sui nodi di ogni elemento
        common /xgv/xgv, /ygv/ygv
        common /nen3v/nen3v

        real ulnv(nlvdim,1),vlnv(nlvdim,1)
        common /ulnv/ulnv, /vlnv/vlnv

	real v_ent ! valore mediato tra velocita' int e out
	real v_int ! modulo della velocita' di entrata
	real v_out ! modulo della velocita' di uscita
	integer l_int !lato di entrata numerazione elemento (var. output)
	integer l_out ! lato di uscita numerazione elemento
	real xbdy,ybdy ! coordinate del body
	real nxbdy,nybdy !nuove coordinate del body
	integer newie ! nuovo elemento che contiene il body

c variabili di servizio
	
	integer pb(3),ip,i,exi(2),ipb(2),exio(2)
	real nwdist,distance ! funzione per il calcolo del nwdist
        real dxx,dyy,ddxx,ddyy,ddd
        real itrx,itry
        real nni
        real xit,yit ! coordinate intercetta retta traj e retta enter
	real fsum,fpb(2),ffpb(2)
	real tcomp ! funzione per il calcolo del deltat
  	real ib,aax,icx,icy ! parametri retta del lato da cui entra il body
	real ocy,ocx,ob ! parametri retta del lato da cui esce il body
        real fay,fax ! parametri fascio rette traiettoria
	real db ! parametri retta della traiettoria sguita dal body 
	real dstbdy ! distanza tra lato entrata e uscita lungo la traiettoria
        real bdx ! distanza (frazione) del body da estremo piu vicno retta entrata
	integer pbdx,gbdx ! numero interno dell'estremo piu vicino e piu lontano dal body
        integer ps,ng,nl,nnm,i1,i2
	data nnm/0/
	data pb/3*0/
	save nnm

c inizializzazione parametri

110 	continue
	nnm=nnm+1
        fsum=0

c==========================================================
c INDIVIDUAZIONE CASI
c==========================================================

        ps=0
        ng=0
        nl=0
        do i=1,3
         nni=vel_ie(i,ie)
         if(nni.lt.0)then
          ng=ng+1
          l_out=i
         elseif(nni.gt.0)then
          ps=ps+1
          l_int=i
         elseif(nni.eq.0)then
          nl=nl+1
         endif          
        end do

        v_ent=sqrt((ulnv(1,ie)**2)+(vlnv(1,ie)**2))
        
c_________________________________________________________
c CASO [A]: 1 flusso >0 (entrante), 2 flussi <0 (uscenti).
c Il lato di entrata è quello il cui flusso >0.
c - individuare fascio di rette traiettoria
c - individuare traiettoria del body
c - individuare lato di uscita
c - individuare nuove coordinate del body

        if((ng.eq.2).and.(ps.eq.1))then
         v_int=vel_ie(l_int,ie)
         
         if(v_int.le.0)then
          PRINT*,'STOP!! ERRORE VELOCITA DI ENTRATA NEGATIVA'
	  PRINT*,'CASO [A] TRACK_ORIG'
	  print*,' ng ', ng,' ps ',ps
          print*,' v_int  ',v_int
	  stop
          return
         endif                               
         
         
c calcolo retta del lato da cui è entrato il body
         
         i1=mod(l_int,3)+1
         i2=mod(i1,3)+1
	 exi(1)=i1
	 exi(2)=i2	 	
                                          
         
         call retta(exi,icy,icx,ib,ie)

         
c determinazione dei lati e dei flussi di uscita probabili pb(1,2),fpb(1,2) 
c pb(mi dice il numero interno del lato)


         ip=1
         do i=1,3
          if(l_int.ne.i)then
	   pb(i)=ip
           ipb(ip)=i
           fpb(ip)= vel_ie(i,ie)
           ip=ip+1
          endif
         enddo
         
         
c si calcola la frazione sul lato di entrata di competenza
c del lato pb(1) e del lato pb(2) per la possibile uscita del body         

        fsum=abs(fpb(1))+abs(fpb(2))
        do i=1,2
         ffpb(i)=abs(fpb(i)/(fsum))
        end do
       
       
c calcolo del fascio di rette traiettorie del body, mi serve sapere:
c (l_int) nodo opposto al lato di entrata coordinate del punto di 
c separazione (sep) tra le 2 frazioni sul lato di entrata (ad estremo 1 sommo ffpd(2))

        call dirz(l_int,ffpb,exi,icy,icx,ib,fay,fax,ie,pb,ipb)

        
c determinazione retta traiettoria
         
        call traj(fay,fax,xbdy,ybdy,db)


c deteminazione coordinate punto intercetto retta traiettoria-
c retta di entrata

        call interc(fay,fax,db,icy,icx,ib,xit,yit)


c calcola la distanza (bdx) del punto (xit,yit) dall'estremo piu vicino
c sulla retta di entrata in termini di frazione della lunghezza totale e
c dice qual'e' il numero interno dell'estremo piu vicino (pbdx) e 
c quello piu lontano (gbdx)

        call distp(xit,yit,exi,bdx,pbdx,gbdx,ie)
        
c individuo il lato da cui esce il body (l_out) mediante confronto
c tra la posizione xit,yit  e le distanze sulla retta di entrata
c individuate dalle frazioni (ffpd) 
        
        if(bdx.le.ffpb(pb(pbdx)))then
         l_out=pbdx
        elseif(bdx.gt.ffpb(pb(pbdx)))then
         l_out=gbdx
        end if

        
c determinazione velocità di entrata v_ent del body v_ent
                
        v_out=vel_ie(l_out,ie)

        if(v_out.gt.0)then
         PRINT*,'STOP!! ERRORE VELOCITÀ DI USCITA POSITIVA'
	 PRINT*,'CASO [A] TRACK_ORIG'	 
	 return
        elseif(v_out.eq.0)then
         PRINT*,'STOP!! ERRORE VELOCITÀ DI USCITA NULLA'
	 PRINT*,'CASO [A] TRACK_ORIG'
         return
        endif

        

c estremi del lato da cui esce il body

        i1=mod(l_out,3)+1
        i2=mod(i1,3)+1

        exio(1)=i1
        exio(2)=i2
        
        call retta(exio,ocy,ocx,ob,ie)

        
c calcolo punto intercetto tra retta traiettoria - retta uscita 

        call interc(fay,fax,db,ocy,ocx,ob,itrx,itry)

c___________________________________________________________
c CASO [B]: 2 flussi >0 (entranti) o 1 flusso >0
c 1 flusso =0, 1 flusso <0 (uscente).
c Il lato di uscita è quello il cui flusso <0.
c - individuare fascio di rette traiettoria
c - individuare traiettoria del body
c - individuare lato di entrata
c - individuare nuove coordinate del body

       elseif(((ng.eq.1).and.(ps.eq.2)).or.
     +((ng.eq.1).and.(nl.eq.1).and.(ps.eq.1)))then
        v_out=vel_ie(l_out,ie)

        
        if(v_out.ge.0)then
         PRINT*,'STOP!! ERRORE VELOCITA DI USCITA POSITIVO'
         PRINT*,'CASO [B] TRACK_ORIG'
         print*,' ng ', ng,' ps ',ps,' nl ',nl
         print*,' v_out  ',v_out
	 stop
         return
        endif

c calcolo retta del lato da cui uscirà il body

	i1=mod(l_out,3)+1
        i2=mod(i1,3)+1

                
        exio(1)=i1
        exio(2)=i2
        
        call retta(exio,ocy,ocx,ob,ie)

        
c determinazione dei lati e dei flussi di entrata probabili pb(1,2),
c fpb(1,2), pb(mi dice il numero interno del lato)
        
        ip=1
        do i=1,3
         if(l_out.ne.i)then
          pb(i)=ip
          ipb(ip)=i
          fpb(ip)=vel_ie(i,ie)
          ip=ip+1
         endif
        end do
        
        
c si calcola la frazione sul lato di uscita, di competenza
c del lato pb(1) e del lato pb(2) per la possibile entrata del body
                             
        fsum=abs(fpb(1))+abs(fpb(2))
        do i=1,2
         ffpb(i)=abs(fpb(i)/fsum)
        end do
     
        
c calcolo del fascio di rette traiettorie del body, mi serve sapere:
c (l_out) nodo opposto al lato di uscita coordinate del punto di
c separazione (sep) tra le 2 frazioni sul lato di entrata (ad estremo 1 sommo ffpd(2))

        call dirz(l_out,ffpb,exio,ocy,ocx,ob,fay,fax,ie,pb,ipb)
        
c determinazione retta traiettoria

        call traj(fay,fax,xbdy,ybdy,db)

c deteminazione coordinate punto intercetto retta traiettoria-
c retta di uscita

        call interc(fay,fax,db,ocy,ocx,ob,itrx,itry)
        
c calcola la distanza (bdx) del punto (itrx,itry) dall'estremo piu vicino
c sulla retta di uscita in termini di frazione della lunghezza totale e
c dice qual'e' il numero interno dell'estremo piu vicino (pbdx) e
c quello piu lontano (gbdx)
        
        call distp(itrx,itry,exio,bdx,pbdx,gbdx,ie)

c individuo il lato da cui entra il body (l_int) mediante confronto
c tra la posizione itrx, itry  e le distanze sulla retta di uscita
c individuate dalle frazioni (ffpd)

        if(bdx.le.ffpb(pb(pbdx)))then
         l_int=pbdx
        elseif(bdx.gt.ffpb(pb(pbdx)))then
         l_int=gbdx
        end if


c determinazione velocità di entrata v_ent del body v_ent
         
        v_int=vel_ie(l_int,ie)
                
        if(v_int.lt.0)then
         PRINT*,'STOP!! ERRORE VELOCITÀ DI ENTRATA NEGATIVA'
         PRINT*,'CASO [B] TRACK_ORIG'
	 print*,' ie ',ie,bdy
	 print*,' '
	 print*,' ng ', ng,' ps ',ps,' nl ',nl
	 print*,' '
	 print*,' exio,ocy,ocx,ob,ie '
	 print*,exio,ocy,ocx,ob,ie
         print*,' '
	 print*,'l_out,ffpb,exio,ocy,ocx,ob,fay,fax,ie,pb,ipb '
	 print*,l_out,ffpb,exio,ocy,ocx,ob,fay,fax,ie,pb,ipb
	 print*,' '
	 print*,'fay,fax,xbdy,ybdy,db'
	 print*, fay,fax,xbdy,ybdy,db
	 print*,' '
	 print*,'fay,fax,db,ocy,ocx,ob,itrx,itry'
	 print*,fay,fax,db,ocy,ocx,ob,itrx,itry
	 print*,' ' 
	 print*,'itrx,itry,exio,bdx,pbdx,gbdx,ie'
	 print*,itrx,itry,exio,bdx,pbdx,gbdx,ie
	 print*,' '
	 print*,' v_int  ',v_int 
         stop
         return
        elseif(v_int.eq.0)then
         PRINT*,'STOP!! ERRORE VELOCITÀ DI ENTRATA NULLA'
	 PRINT*,'CASO [B] TRACK_ORIG'
	 print*,' v_int  ',v_int
	 stop
         return
        endif
c____________________________________________________________________________
c CASO [C]:tutti i flussi sono entranti >0 O NULLI =0
c Errore il campo di moto è convergente, problemi con l'eliminazione delle componenti
c DIVERGENTI FERMARE IL CILO

       elseif((ps.eq.3).or.(nl.eq.3))then
        print*,'STOP!! FLUSSI TUTTI POSITIVI O NULLI IN ',ie
	PRINT*,'CASO [C] TRACK_ORIG'
         it=0
        return
       elseif((ng.eq.3))then
        print*,'STOP!! FLUSSI TUTTI NEGATIVI IN ',ie
        PRINT*,'CASO [D] TRACK_ORIG'
         it=0
        return
       elseif((nl.eq.2).and.(ps.eq.1))then
        print*,'STOP!! UN SOLO FLUSSO POSITIVO E 2 NULLI IN ',ie
        PRINT*,'CASO [E] TRACK_ORIG'
        PRINT*,'[elemento di partenza con due lati C.B]'
         it=0
        return
       elseif((nl.eq.2).and.(ng.eq.1))then
        print*,'STOP!! UN SOLO FLUSSO NEGATIVO E 2 NULLI IN ',ie
        PRINT*,'CASO [F] TRACK_ORIG'
	PRINT*,'[elemento di partenza con due lati C.B]'
         it=0
        return
       elseif((nl.eq.1).and.(ps.eq.2))then
        print*,'STOP!! 2 FLUSSI POSITIVI E 1 NULLO IN ',ie
        PRINT*,'CASO [G] TRACK_ORIG'
         it=0
        return
       elseif((nl.eq.1).and.(ng.eq.2))then
        print*,'STOP!! 2 FLUSSI NEGATIVI E 1 NULLO IN ',ie
        PRINT*,'CASO [H] TRACK_ORIG'
         it=0
        return
       end if         
c_________________________________________________________________________

c Introduco una velocita' aggiuntiva u_adj, v_adj.
c ricalcolo  la retta traiettoria e
c il punto intercetto sulla retta d'uscita,
c e la nuova distanza newdist.

c calcolo coordinate punto raggiunto lungo traiettoria
c in it

c        call lpoint(fay,fax,db,xt,yt)

c calcolo punto xa,ya raggiunto da posizione iniziale
c in it con velocita' aggiuntiva

c        call adjvel(xdby,ydby,it,xa,ya)

c calcolo parametri nuova traiettoria acyY=acxX+ab

c        call retnd(xdby,ydby,xa,ya,xt,yt,acy,acx,ab,xf,yf)

c calcolo intercetta tra retta nuova traiettoria e retta uscita

c        call interc(acy,acx,ab,ocy,ocx,ob,itrx,itry)

c calcolo nuova velocita' data da rapporto distanza segmento
c su nuova traiettoria e intervallo di tempo

c        v_ent=sqrt(((xdby-xf)**2)+((ydby-yf)**2))/it

c calcolo distanza massima percorribile all'interno di elemento
       
       dxx=itrx-xbdy
       dyy=itry-ybdy
       ddxx=(dxx)**2
       ddyy=(dyy)**2
       ddd=ddxx+ddyy
       dstbdy=sqrt(ddd)
       nwdist=it*v_ent 
        
c=====================================================================
c CALCOLO NUOVE COORDINATE DEL BODY (nxbdy,nybdy)
c=====================================================================

c calcolo della distanza percorsa in 1 tstep (nwdist)
c se il body arriva esattamente sul lato a ttime=0 (1 iflogico)?
c se si esce dall'elemento (2 iflogico) allora ho coordinate del body sul lato di uscita
c e in piu un deltat > 0 da utilizzare durnate questo stesso tstep.
c se si rimane nell'elemento (3 iflogico) allora calcolo le nuove coordinate 
c (pnt_inside) e passo al tstep successivo

c Se il body finisce a 0.1 m dalla linea decido si farlo passare al
c elemento ad una distanza di 0.1 m. Questa è la precisione da me 
c imposta
               
        if(abs(nwdist-dstbdy).lt.0.01)then      ! body dista meno di 0.1m dal lato
         PRINT*,'STOP!! BODY SU LATO ELEMENTO'
         PRINT*,'track orig'
         nwdist=dstbdy+0.01
        endif
        if(nwdist.gt.dstbdy)then            !2 body su nuovo el
	 nxbdy=itrx
	 nybdy=itry
	 distance=nwdist-dstbdy
	 deltat=distance/v_ent
         newie=ieltv(l_out,ie)
         if(newie.eq.-1)then			!2.a body e uscito dal dominio
          PRINT*,'STOP!! BODY ',bdy,' USCITO '
          PRINT*,'ELEMENTO USCITA ',ie
          ie=-ie ! flag per skip bdy dal calcolo
	  it=0.
          xbdy=nxbdy ! ultime coordinate del body uscito
          ybdy=nybdy ! ultime coordinate del body uscito
          return
         endif
         do i=1,3
          if(ie.eq.ieltv(i,newie))then          !3 body su vecchio el
           lt_body(bdy)=i ! individuazione lato di entrata del body
          endif
         end do
         xbdy=nxbdy ! nuove coordinate del body
	 ybdy=nybdy ! nuove coordinate del body
         ie=newie ! individuazione nuovo elemento
         it=deltat
  	elseif(nwdist.lt.dstbdy)then
	 call pnt_inside(nwdist,dstbdy,fay,fax,xbdy,ybdy,itrx,itry,nxbdy,nybdy)
	 deltat=0.
	 it=deltat
	 ie=ie ! individuazione vecchio elemento
         xbdy=nxbdy ! nuove coordinate del body
         ybdy=nybdy ! nuove coordinate del body
        endif
        end

c#############################################################################
	subroutine track_line(it,bdy,ie,xbdy,ybdy)

c in questa subroutine si calcola il percorso del body partendo
c da un punto sul lato dell'elemento ie
c questa subroutine e' chiamata ogni qualvolta un body
c nello stesso timestep arriva in un nuovo elemento

        implicit none
        
	include 'param.h' 
	include 'lagrange.h' 

	integer ie,it,bdy ! bdy e' il numero del body
 
        real deltat ! frazione di time step body si muove in ie
 
	integer ieltv(3,1)
        common /ieltv/ieltv !punta ad elemento opposto al nodo 

	real xgv(1),ygv(1) ! coordinate nodi
        integer nen3v(3,1) ! puntatore sui nodi di ogni elemento
        common /xgv/xgv, /ygv/ygv
        common /nen3v/nen3v

        real ulnv(nlvdim,1),vlnv(nlvdim,1)
        common /ulnv/ulnv, /vlnv/vlnv
                
	real v_ent ! valore mediato tra velocita' int e out
	real v_int ! modulo della velocita' di entrata
	real v_out ! modulo della velocita' di uscita
	integer l_int !lato di entrata numerazione elemento (var. input)
	integer l_out ! lato di uscita numerazione elemento
	real xbdy,ybdy ! coordinate del body
	real nxbdy,nybdy !nuove coordinate del body
	integer newie ! nuovo elemento che contiene il body
c variabili di servizio
	
	integer pb(3),ip,i,exi(2),ipb(2),exio(2)
	real nwdist,distance ! funzione per il calcolo del nwdist
        real dxx,dyy,ddxx,ddyy,ddd
        real itrx,itry
	real fsum,fpb(2),ffpb(2)
	real tcomp ! funzione per il calcolo del deltat
  	real icy,icx,ib ! parametri retta del lato da cui entra il body
	real ocy,ocx,ob ! parametri retta del lato da cui esce il body
        real fay,fax ! parametri fascio rette traiettoria
        real db ! parametri retta della traiettoria sguita dal body 
	real dstbdy ! distanza tra lato entrata e uscita lungo la traiettoria
        real bdx ! distanza (frazione) del body da estremo piu vicno retta entrata
	integer pbdx,gbdx ! numero interno dell'estremo piu vicino e piu lontano dal body
        integer ps,ng,nl,i1,i2

c inizializzazione parametri

110 	continue

        fsum=0
        
c==========================================================
c CALCOLO TRAIETTORIE
c==========================================================

c calcolo retta del lato da cui entra il body
c exi(1,2) sono i puntatori del nen3v per gli estremi della retta 
c su cui si trova il body

        l_int=lt_body(bdy)


        i1=mod(l_int,3)+1
        i2=mod(i1,3)+1
        exi(1)=i1
        exi(2)=i2
                                          
	call retta(exi,icy,icx,ib,ie)

c velocità di entrata v_int 

	v_int=vel_ie(l_int,ie)

        if(v_int.le.0)then
         PRINT*,'STOP!! ERRORE FLUSSO DI ENTRATA NEGATIVO'
         print*,'track line'
	 print*,' v_int ',v_int
         return
        endif

c determinazione dei lati e dei flussi di uscita probabili pb(1,2),vpb(1,2)
c pb(mi dice il numero interno del lato)

        ip=1
        do i=1,3
	 if(l_int.ne.i)then
 	  pb(i)=ip
	  ipb(ip)=i
          fpb(ip)=vel_ie(i,ie) 
          ip=ip+1  
         endif 
	end do

        v_ent=sqrt((ulnv(1,ie)**2)+(vlnv(1,ie)**2))
         
         
c______________________________________________________________________
c CASO [A] flussi ai lati opposti al l_int entrambi <0
c si calcola la frazione sul lato di entrata di competenza 
c del lato pb(1) e del lato pb(2) per la possibile uscita del body.

      	if((fpb(1).lt.0).and.(fpb(2).lt.0))then

                
        fsum=abs(fpb(1))+abs(fpb(2))
         do i=1,2
	  ffpb(i)=abs(fpb(i)/fsum)
	 end do

         
c calcola la distanza del body (bdx) dall'estremo piu vicino
c in termini di frazione sulla retta totale  e
c dice qual'e' il numero interno dell'estremo piu vicino (pbdx)
c e quello piu lontano (gbdx)

        call distp(xbdy,ybdy,exi,bdx,pbdx,gbdx,ie)
        
c individuo il lato da cui esce il body mediante confronto
c tra posizione del body rispetto alle frazioni  

	 if(bdx.le.ffpb(pb(pbdx)))then	  
	  l_out=pbdx
	 elseif(bdx.gt.ffpb(pb(pbdx)))then
	  l_out=gbdx
	 end if

         
c estremi del lato da cui esce il body

        i1=mod(l_out,3)+1
        i2=mod(i1,3)+1

        exio(1)=i1
        exio(2)=i2

                                          
c equazione della retta del lato di uscita

        call retta(exio,ocy,ocx,ob,ie)


c calcolo del fascio di traiettorie del body, mi serve sapere:
c (l_int) nodo opposto al lato di entrata
c coordinate del punto di separazione (sep) tra le 2 frazioni sul lato
c di entrata (ad estremo 1 sommo ffpd(2))
         
        call dirz(l_int,ffpb,exi,icy,icx,ib,fay,fax,ie,pb,ipb)
        
c determinazione retta traiettoria

        call traj(fay,fax,xbdy,ybdy,db)

c deteminazione coordinate punto intercetto retta traiettoria-
c retta di uscita
         
        call interc(fay,fax,db,ocy,ocx,ob,itrx,itry)

c_____________________________________________________________________________
c CASO [B] flussi ai lati opposti al l_int uno >0 e altro <0
c in questo caso il lato di uscita e' sicuramente quello dove ho la velocita'
c negativa, non serve quindi il calcolo delle rette di partenza e di arrivo
c gli estremi della retta del lato di uscita sono i due estremi non opposti
c alla retta che ha valore di vel negativo

        else if(((fpb(1).lt.0).and.(fpb(2).ge.0)).or.
     + ((fpb(1).ge.0).and.(fpb(2).lt.0)))then
        

         do i=1,2
          if(fpb(i).lt.0)then
           l_out=ipb(i)
          endif
         enddo

c estremi del lato da cui esce il body
	
        i1=mod(l_out,3)+1
        i2=mod(i1,3)+1

        exio(1)=i1
        exio(2)=i2
                                          
                                          
c equazione della retta del lato di uscita

        call retta(exio,ocy,ocx,ob,ie)
        
        
c individuo i 2 flussi di entrata del body nell'elemento per il calcolo
c delle frazioni di competenza sul lato di uscita.
 
         ip=1
         do i=1,3
          if(l_out.ne.i)then
           pb(i)=ip
           ipb(ip)=i
           fpb(ip)=vel_ie(i,ie)
           ip=ip+1
          endif
         end do

         
c calcolo delle frazioni sul lato di uscita per determinare 
c la retta traiettoria del percorso del body
	
         fsum=abs(fpb(1))+abs(fpb(2))
         do i=1,2
          ffpb(i)=abs(fpb(i)/fsum)
         end do


c calcolo del fascio di traiettorie del body, mi serve sapere:
c l_out nodo opposto al lato di uscita
c coordinate del punto di separazione tra le 2 frazioni sul lato 
c di uscita (come prima ad estremo 1 sommo ffpb(2))

         call dirz(l_out,ffpb,exio,ocy,ocx,ob,fay,fax,ie,pb,ipb)
         
c calcolo retta traiettoria 

         call traj(fay,fax,xbdy,ybdy,db)
         
c deteminazione coordinate punto intercetto retta traiettoria-
c retta di uscita

         call interc(fay,fax,db,ocy,ocx,ob,itrx,itry)

c calcolo velocita' in uscita v_out per il calcolo del modulo della
c velocita di percorso v_ent pari alla media v_int, v_out

         v_out=vel_ie(l_out,ie)
        
         if(v_out.gt.0)then
          PRINT*,'STOP!! ERRORE VELOCITÀ DI USCITA POSITIVA'
          PRINT*,'track line'
	  print*,' v_out ',v_out
          return
         elseif(v_out.eq.0)then
          PRINT*,'STOP!! ERRORE VELOCITÀ DI USCITA NULLA'
          PRINT*,'track line'
	  print*,' v_out ',v_out
          return
         endif

        elseif((fpb(1).ge.0).and.(fpb(2).ge.0))then
c________________________________________________________________________________
c CASO [C] tutti i flussi sono entranti o nulli. Errore il campo di moto è
c convergente, problemi con l'eliminazione delle componenti DIVERGENTI
c FERMARE IL CILO

         print*,'STOP!! FLUSSI TUTTI POSITIVI IN ',ie
	 PRINT*,'CASO [C] TRACK_LINE'
         it=0
         return
        elseif((fpb(1).lt.0).and.(fpb(2).lt.0))then
         print*,'STOP!! FLUSSI TUTTI NEGATIVI IN ',ie
         PRINT*,'CASO [D] TRACK_LINE'
         it=0
         return
        elseif((nl.eq.2).and.(ps.eq.1))then
         print*,'STOP!! UN SOLO FLUSSO POSITIVO E 2 NULLI IN ',ie
         PRINT*,'CASO [E] TRACK_LINE'
	 PRINT*,'[elemento con due lati C.B]'
         it=0
         return
        elseif((nl.eq.2).and.(ng.eq.1))then
         print*,'STOP!! UN SOLO FLUSSO NEGATIVO E 2 NULLI IN ',ie
         PRINT*,'CASO [F] TRACK_LINE'
	 PRINT*,'[elemento con due lati C.B]'
         it=0
         return
        elseif((nl.eq.1).and.(ps.eq.2))then
         print*,'STOP!! 2 FLUSSI POSITIVI E 1 NULLO IN ',ie
         PRINT*,'CASO [G] TRACK_LINE'
         it=0
         return
        elseif((nl.eq.1).and.(ng.eq.2))then
         print*,'STOP!! 2 FLUSSI NEGATIVI E 1 NULLO IN ',ie
         PRINT*,'CASO [H] TRACK_LINE'
         it=0
         return
        endif
c_______________________________________________________________________
c Introduco una velocita' aggiuntiva u_adj, v_adj.
c ricalcolo  la retta traiettoria e
c il punto intercetto sulla retta d'uscita,
c e la nuova distanza newdist.

c calcolo coordinate punto raggiunto lungo traiettoria
c in it

c        call lpoint(fay,fax,db,xt,yt)

c calcolo punto xa,ya raggiunto da posizione iniziale
c in it con velocita' aggiuntiva

c        call adjvel(xdby,ydby,it,xa,ya)

c calcolo parametri nuova traiettoria acyY=acxX+ab

c        call retnd(xdby,ydby,xa,ya,xt,yt,acy,acx,ab,xf,yf)

c calcolo intercetta tra retta nuova traiettoria e retta uscita

c        call interc(acy,acx,ab,ocy,ocx,ob,itrx,itry)

c calcolo nuova velocita' data da rapporto distanza segmento
c su nuova traiettoria e intervallo di tempo

c        v_ent=sqrt(((xdby-xf)**2)+((ydby-yf)**2))/it

c calcolo distanza massima percorribile all'interno di elemento

        dxx=itrx-xbdy
        dyy=itry-ybdy
        ddxx=(dxx)**2
        ddyy=(dyy)**2
        ddd=ddxx+ddyy
        dstbdy=sqrt(ddd)

	nwdist=it*v_ent


c=====================================================================
c CALCOLO NUOVE COORDINATE DEL BODY (nxbdy,nybdy)
c=====================================================================

c il body entra sempre da un lato con una velocita' positiva
c calcolo della distanza percorsa in 1 tstep (nwdist)
c se il body arriva esattamente sul lato a ttime=0 (1 iflogico)?
c se si esce dall'elemento (2 iflogico) allora ho coordinate del body sul lato di uscita
c e in piu un deltat > 0 da utilizzare durnate questo stesso tstep.
c se si rimane nell'elemento (3 iflogico) allora calcolo le nuove coordinate 
c (pnt_inside) e passo al tstep successivo

c Se il body finisce a 0.1 m dalla linea decido si farlo passare al
c elemento ad una distanza di 0.1 m. Questa è la precisione da me
c imposta

        if(abs(nwdist-dstbdy).lt.0.01)then     ! body dista meno di 0.1m dal lato    
         PRINT*,'STOP!! BODY SU LATO ELEMENTO'
         PRINT*,'track line'
         nwdist=dstbdy+0.01
        endif
        if(nwdist.gt.dstbdy)then            !2 body su nuovo el
	 nxbdy=itrx
	 nybdy=itry
	 distance=nwdist-dstbdy
         deltat=distance/v_ent
         newie=ieltv(l_out,ie) !
	 print*,' newie ',newie,' it ',it 
         if(newie.eq.-1)then			!2.a body uscito dal dominio
          PRINT*,'STOP!! BODY ',bdy,' USCITO '
          PRINT*,'ELEMENTO USCITA ',ie
          ie=-ie! flag per skip bdy dal calcolo
          xbdy=nxbdy ! ultime coordinate del body uscito
          ybdy=nybdy ! ultime coordinate del body uscito
	  it=0.
	  return     
         endif
         do i=1,3
          if(ie.eq.ieltv(i,newie))then          !3 body su vecchio el
           lt_body(bdy)=i ! individuazione lato di entrata del body
          endif
         end do
         xbdy=nxbdy ! nuove coordinate del body
	 ybdy=nybdy ! nuove coordinate del body
        ie=newie ! individuazione nuovo elemento
        it=deltat
       elseif(nwdist.lt.dstbdy)then
	 call pnt_inside(nwdist,dstbdy,fay,fax,xbdy,ybdy,itrx,itry,nxbdy,nybdy)
	 deltat=0.
	 it=deltat
	 ie=ie ! individuazione vecchio elemento
         xbdy=nxbdy ! nuove coordinate del body
         ybdy=nybdy ! nuove coordinate del body
       endif
       end

c############################################################################
	
	subroutine retta(ext,cy,cx,b,el)
	
c dati gli estremi calcolo i coefficienti a, b della retta passante
	
	implicit none
	
	integer ext(2),i,el
	integer p1,p2
	
	real xgv(1),ygv(1) ! coordinate nodi
        integer nen3v(3,1) ! puntatore sui nodi di ogni elemento
        common /xgv/xgv, /ygv/ygv
        common /nen3v/nen3v
	
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

c____________________________________________________________________________

       subroutine dirz(pfin,ratio,ext,cy,cx,b,ny,nx,el,pb,ipb)
	
c per il calcolo della retta traiettoria percorsa del body:
c ho bisogno del punto di arrivo, pfin e un punto di partenza pini
c                         | questo punto lo calcolo considerando: 
c         /\ pfin         | differenza deltax tra le coord. x degli estremi
c        /  \             | della retta y=ax+b, applico il rapporto
c       /    \            | (ratio) dei segmenti a tale delta X e  
c      /      \           | trovo la frazione di x da aggiungere  
c     /  ratio \          | all'estremo opposto a cui si riferisce il 
c    <------------------> | ratio. Quindi trovo le nuove coordinate del punto
c	        pini      | pini. Da pini e pfin trovo il fascio di
c_________________________| rette in particolare an il coefficiente
c                           angolare an
	implicit none
	
	real xgv(1),ygv(1) ! coordinate nodi
        integer nen3v(3,1) ! puntatore sui nodi di ogni elemento
        common /xgv/xgv, /ygv/ygv
        common /nen3v/nen3v
	
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

	do i=1,2
	 if(ext(1).ne.ipb(i))then
	  dfrq=dist*(ratio(i))
	  rr=ratio(i)
	 endif
	enddo 
        
c calcolo coordinate punto nella retta distante dfrq da p1 e
c compreso tra p1 e p2

        
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


c______________________________________________________________

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

c______________________________________________________________

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

c__________________________________________________________________
	
	subroutine distp(x,y,ext,d,near,far,el)

c calcolo della distanza (dist) del punto x,y dall'estremo piu distante

	implicit none

	real x,y ! punto da cui determinare la distanza
	integer ext(2) ! puntatore estremi del segmento
	real d ! distanza dal punto piu lontano
	integer near,far ! puntatore dell'estremo piu lontano e piu vicino
			 ! rispetto al punto x,y

	real xgv(1),ygv(1) ! coordinate nodi
        integer nen3v(3,1) ! puntatore sui nodi di ogni elemento
        common /xgv/xgv, /ygv/ygv
        common /nen3v/nen3v

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
         PRINT*,'STOP BODY IN X ',x,' Y ',y,' NON CONTENUTO IN EL ',el
        endif
        end

c_______________________________________________________

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

		











 















