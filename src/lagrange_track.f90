!
! $Id: lagrange_track.f,v 1.1 2009-02-13 17:22:44 georg Exp $
!
! routines for tracking particle
!
! revision log :
!
! 05.02.2009    ggu     copied from other files
! 16.12.2011    ggu     write all messages to lunit
! 23.01.2012    ggu     ltbdy is locally passed (not common)
! 02.03.2012    ggu&fra introduced epsggu to avoid nan in time
! 25.01.2013    ggu     error check to avoid segfault (INTERNAL ERROR)
!
!**********************************************************************
!-------------------------------------------------------------------------
        module lagrange_track
!-------------------------------------------------------------------------
        contains
!-------------------------------------------------------------------------

        subroutine track_orig(time,bdy,ie,xbdy,ybdy,zbdy,lybdy,ltbdy)

! in questa subroutine si calcola il percorso del body partendo
! da un punto all'interno dell'elemento ie 
! fino al successivo punto interno allo stesso elemento
! se entro il time step il body esce dall''elemento 
! si passa alla successiva subroutine con time

	use lagrange_data
	use lagrange_vertical
	use geom
	use hydro_vel
	use levels, only : nlvdi,nlv
	use basin
        use regular
	use lagrange_util_tr

        implicit none
        
	include 'param.h' 
 
	double precision time	!total time to travel
	integer bdy	!number of body
	integer ie	!number of element of body
	double precision xbdy,ybdy	!horizontal position of body
	double precision zbdy	!relative vertical position of body in layer
	integer lybdy	!layer of body
	integer ltbdy	!entering side of next element 

	double precision epsggu
	parameter (epsggu = 1.e-7)
 
        double precision deltat ! frazione di time step body si muove in ie
 



	double precision v_ent ! valore mediato tra velocita' int e out
	double precision v_int ! modulo della velocita' di entrata
	double precision v_out ! modulo della velocita' di uscita
	integer l_int !lato di entrata numerazione elemento (var. output)
	integer l_out ! lato di uscita numerazione elemento
	double precision nxbdy,nybdy !nuove coordinate del body
	integer newie ! nuovo elemento che contiene il body

! vertical treatment 
	double precision z0,z1,zn0,zn1  !	cuccomod
	integer l0,l1,addl,lmax
	double precision ztime,lb,subtime,layd
	double precision w
	double precision hl(nlvdi)
	double precision a_int,a_out
        double precision in_d,ou_d
	integer in_dm,in_dx,ou_dm,ou_dx
	
        integer idum ! seme f(x) random

! variabili di servizio
	
	integer pb(3),ip,i,exi(2),ipb(2),exio(2)
	double precision nwdist,distance ! funzione per il calcolo del nwdist
        double precision dxx,dyy,ddxx,ddyy,ddd
        double precision itrx,itry
        double precision drd(3)
        double precision nni
        double precision xit,yit ! coordinate intercetta retta traj e retta enter
	double precision fsum,fpb(2),ffpb(2)
	double precision tcomp ! funzione per il calcolo del deltat
  	double precision ib,aax,icx,icy ! parametri retta del lato da cui entra il body
	double precision ocy,ocx,ob ! parametri retta del lato da cui esce il body
        double precision fay,fax ! parametri fascio rette traiettoria
	double precision db ! parametri retta della traiettoria sguita dal body 
	double precision dstbdy ! distanza tra lato entrata e uscita lungo la traiettoria
        double precision bdx ! distanza (frazione) del body da estremo piu vicno retta entrata
	double precision gf !precision along the side
        integer pbdx,gbdx ! numero interno dell'estremo piu vicino e piu lontano dal body
        integer ps,ng,nl,nnm,i1,i2
	data nnm/0/
	data pb/3*0/
        data idum/76342/
        save nnm
        save idum

! inizializzazione parametri

	l_int = 0
	l_out = 0

	nnm=nnm+1
        fsum=0
	lb = lybdy

!==========================================================
! INDIVIDUAZIONE CASI
!==========================================================
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
        
!_________________________________________________________
! CASO [A]: 1 flusso >0 (entrante), 2 flussi <0 (uscenti).
! Il lato di entrata è quello il cui flusso >0.
! - individuare fascio di rette traiettoria
! - individuare traiettoria del body
! - individuare lato di uscita
! - individuare nuove coordinate del body

        if((ng.eq.2).and.(ps.eq.1))then
         v_int=vel_ie(l_int,ie)
         
         if(v_int.le.0)then
!          write(lunit,*) 'STOP!! ERRORE VELOCITA DI ENTRATA NEGATIVA'
!	  write(lunit,*) 'CASO [A] TRACK_ORIG'
	  !stop
          ie = -ie
	  time=0
	  return
         endif                               
         
         
! calcolo retta del lato da cui è entrato il body
         
         i1=mod(l_int,3)+1
         i2=mod(i1,3)+1
	 exi(1)=i1
	 exi(2)=i2	 	
                                          
         
         call retta(exi,icy,icx,ib,ie)

         
! determinazione dei lati e dei flussi di uscita probabili pb(1,2),fpb(1,2) 
! pb(mi dice il numero interno del lato)


         ip=1
         do i=1,3
	  pb(i) = 0
          if(l_int.ne.i)then
	   pb(i)=ip
           ipb(ip)=i
           fpb(ip)= vel_ie(i,ie)
           ip=ip+1
          endif
         enddo
         
         
! si calcola la frazione sul lato di entrata di competenza
! del lato pb(1) e del lato pb(2) per la possibile uscita del body         

        fsum=abs(fpb(1))+abs(fpb(2))
        do i=1,2
         ffpb(i)=abs(fpb(i)/(fsum))
        end do
       
       
! calcolo del fascio di rette traiettorie del body, mi serve sapere:
! (l_int) nodo opposto al lato di entrata coordinate del punto di 
! separazione (sep) tra le 2 frazioni sul lato di entrata (ad estremo 1 sommo ffpd(2))

        call dirz(l_int,ffpb,exi,icy,icx,ib,fay,fax,ie,pb,ipb)

        
! determinazione retta traiettoria
         
        call traj(fay,fax,xbdy,ybdy,db)


! deteminazione coordinate punto intercetto retta traiettoria-
! retta di entrata

        call interc(fay,fax,db,icy,icx,ib,xit,yit)


! calcola la distanza (bdx) del punto (xit,yit) dall''estremo piu vicino
! sulla retta di entrata in termini di frazione della lunghezza totale e
! dice qual''e'' il numero interno dell''estremo piu vicino (pbdx) e 
! quello piu lontano (gbdx)

        call distp(xit,yit,exi,bdx,pbdx,gbdx,ie)
	
	in_d=bdx
	in_dm=pbdx
	in_dx=gbdx
       
        if(ie.le.0)then
         call find_elem_from_old(-ie,xbdy,ybdy,ie)
         time=0
         return
        end if 
! individuo il lato da cui esce il body (l_out) mediante confronto
! tra la posizione xit,yit  e le distanze sulla retta di entrata
! individuate dalle frazioni (ffpd) 
        
	if(pbdx.lt.1.or.pbdx.gt.3.or.pb(pbdx).le.0) then
!         write(lunit,*) 'STOP!! ERRORE INTERNO (1)'
!	 write(lunit,*) 'CASO [A] TRACK_ORIG'	 
	 time=0
         ie = -ie
	 return
	end if

        if(bdx.le.ffpb(pb(pbdx)))then
         l_out=pbdx
        elseif(bdx.gt.ffpb(pb(pbdx)))then
         l_out=gbdx
        end if

        
! determinazione velocità di entrata v_ent del body v_ent
                
        v_out=vel_ie(l_out,ie)

        if(v_out.gt.0)then
!         write(lunit,*) 'STOP!! ERRORE VELOCITÀ DI USCITA POSITIVA'
!	 write(lunit,*) 'CASO [A] TRACK_ORIG'	 
	 !stop
	 time=0
         ie = -ie
	 return
        elseif(v_out.eq.0)then
!         write(lunit,*) 'STOP!! ERRORE VELOCITÀ DI USCITA NULLA'
!	 write(lunit,*) 'CASO [A] TRACK_ORIG'
	 !stop
	 time=0
         ie = -ie
	 return
        endif

        

! estremi del lato da cui esce il body

        i1=mod(l_out,3)+1
        i2=mod(i1,3)+1

        exio(1)=i1
        exio(2)=i2
        
        call retta(exio,ocy,ocx,ob,ie)

        
! calcolo punto intercetto tra retta traiettoria - retta uscita 

        call interc(fay,fax,db,ocy,ocx,ob,itrx,itry)
	
	call distp(itrx,itrx,exio,bdx,pbdx,gbdx,ie)

	ou_d =bdx
        ou_dm=pbdx
        ou_dx=gbdx
	

!___________________________________________________________
! CASO [B]: 2 flussi >0 (entranti) o 1 flusso >0
! 1 flusso =0, 1 flusso <0 (uscente).
! Il lato di uscita è quello il cui flusso <0.
! - individuare fascio di rette traiettoria
! - individuare traiettoria del body
! - individuare lato di entrata
! - individuare nuove coordinate del body

       elseif(((ng.eq.1).and.(ps.eq.2)).or.((ng.eq.1).and.(nl.eq.1).and.(ps.eq.1)))then
        v_out=vel_ie(l_out,ie)

        
        if(v_out.ge.0)then
!         write(lunit,*) 'STOP!! ERRORE VELOCITA DI USCITA POSITIVA'
!         write(lunit,*) 'CASO [B] TRACK_ORIG'
	 !stop
	 time=0
         ie = -ie
	 return
        endif

! calcolo retta del lato da cui uscirà il body

	i1=mod(l_out,3)+1
        i2=mod(i1,3)+1

                
        exio(1)=i1
        exio(2)=i2
        
        call retta(exio,ocy,ocx,ob,ie)

        
! determinazione dei lati e dei flussi di entrata probabili pb(1,2),
! fpb(1,2), pb(mi dice il numero interno del lato)
        
        ip=1
        do i=1,3
	 pb(i) = 0
         if(l_out.ne.i)then
          pb(i)=ip
          ipb(ip)=i
          fpb(ip)=vel_ie(i,ie)
          ip=ip+1
         endif
        end do
        
        
! si calcola la frazione sul lato di uscita, di competenza
! del lato pb(1) e del lato pb(2) per la possibile entrata del body
                             
        fsum=abs(fpb(1))+abs(fpb(2))
        do i=1,2
         ffpb(i)=abs(fpb(i)/fsum)
        end do
     
        
! calcolo del fascio di rette traiettorie del body, mi serve sapere:
! (l_out) nodo opposto al lato di uscita coordinate del punto di
! separazione (sep) tra le 2 frazioni sul lato di entrata (ad estremo 1 sommo ffpd(2))

        call dirz(l_out,ffpb,exio,ocy,ocx,ob,fay,fax,ie,pb,ipb)
        
! determinazione retta traiettoria

        call traj(fay,fax,xbdy,ybdy,db)

! deteminazione coordinate punto intercetto retta traiettoria-
! retta di uscita

        call interc(fay,fax,db,ocy,ocx,ob,itrx,itry)
        
! calcola la distanza (bdx) del punto (itrx,itry) dall''estremo piu vicino
! sulla retta di uscita in termini di frazione della lunghezza totale e
! dice qual''e'' il numero interno dell''estremo piu vicino (pbdx) e
! quello piu lontano (gbdx)
        
        call distp(itrx,itry,exio,bdx,pbdx,gbdx,ie)

	ou_d=bdx
        ou_dm=pbdx
        ou_dx=gbdx


        if(ie.le.0)then
         call find_elem_from_old(-ie,xbdy,ybdy,ie)
         time=0
         return
        end if
                                                           
! individuo il lato da cui entra il body (l_int) mediante confronto
! tra la posizione itrx, itry  e le distanze sulla retta di uscita
! individuate dalle frazioni (ffpd)

	if(pbdx.lt.1.or.pbdx.gt.3.or.pb(pbdx).le.0) then
!         write(lunit,*) 'STOP!! ERRORE INTERNO (2)'
!	 write(lunit,*) 'CASO [B] TRACK_ORIG'	 
	 time=0
         ie = -ie
	 return
	end if

        if(bdx.le.ffpb(pb(pbdx)))then
         l_int=pbdx
        elseif(bdx.gt.ffpb(pb(pbdx)))then
         l_int=gbdx
        end if

! estremi del lato da cui entra il body

        i1=mod(l_int,3)+1
        i2=mod(i1,3)+1

        exi(1)=i1
        exi(2)=i2

        call retta(exi,ocy,ocx,ob,ie)


! calcolo punto intercetto tra retta traiettoria - retta entrata

        call interc(fay,fax,db,ocy,ocx,ob,xit,yit)

	call distp(xit,yit,exi,bdx,pbdx,gbdx,ie)

	in_d=bdx
        in_dm=pbdx
        in_dx=gbdx

! determinazione velocità di entrata v_ent del body v_ent
         
        v_int=vel_ie(l_int,ie)
                
        if(v_int.lt.0)then
!         write(lunit,*) 'STOP!! ERRORE VELOCITÀ DI ENTRATA NEGATIVA'
!         write(lunit,*) 'CASO [B] TRACK_ORIG'
	 !stop
	 time=0
         ie = -ie
	 return
        elseif(v_int.eq.0)then
!         write(lunit,*) 'STOP!! ERRORE VELOCITÀ DI ENTRATA NULLA'
!	 write(lunit,*) 'CASO [B] TRACK_ORIG'
	 !stop
	 time=0
         ie = -ie
	 return
        endif
!____________________________________________________________________________
! CASO [C]:tutti i flussi sono entranti >0 O NULLI =0
! Errore il campo di moto è convergente, problemi con l''eliminazione delle componenti
! DIVERGENTI FERMARE IL CILO

       elseif((ps.eq.3).or.(nl.eq.3))then
!       write(lunit,*) 'STOP!! FLUSSI TUTTI POSITIVI O NULLI IN ',ie
!       write(lunit,*) 'CASO [C] TRACK_ORIG'
         time=0
!        ie = -ie
        return
       elseif((ng.eq.3))then
!       write(lunit,*) 'STOP!! FLUSSI TUTTI NEGATIVI IN ',ie
!       write(lunit,*) 'CASO [D] TRACK_ORIG'
         time=0
!        ie = -ie
        return
       elseif((nl.eq.2).and.(ps.eq.1))then
!       write(lunit,*) 'STOP!! UN SOLO FLUSSO POSITIVO E 2 NULLI IN ',ie
!       write(lunit,*) 'CASO [E] TRACK_ORIG'
!       write(lunit,*) '[elemento di partenza con due lati C.B]'
         time=0
!        ie = -ie
        return
       elseif((nl.eq.2).and.(ng.eq.1))then
!       write(lunit,*) 'STOP!! UN SOLO FLUSSO NEGATIVO E 2 NULLI IN ',ie
!       write(lunit,*) 'CASO [F] TRACK_ORIG'
!       write(lunit,*) '[elemento di partenza con due lati C.B]'
         time=0
!        ie = -ie
        return
       elseif((nl.eq.1).and.(ps.eq.2))then
!       write(lunit,*) 'STOP!! 2 FLUSSI POSITIVI E 1 NULLO IN ',ie
!       write(lunit,*) 'CASO [G] TRACK_ORIG'
         time=0
!        ie = -ie
        return
       elseif((nl.eq.1).and.(ng.eq.2))then
!       write(lunit,*) 'STOP!! 2 FLUSSI NEGATIVI E 1 NULLO IN ',ie
!       write(lunit,*) 'CASO [H] TRACK_ORIG'
         time=0
!        ie = -ie
        return
       end if         
!_________________________________________________________________________


! calcolo distanza massima percorribile all''interno di elemento
       
       dxx=itrx-xbdy
       dyy=itry-ybdy
       ddxx=(dxx)**2
       ddyy=(dyy)**2
       ddd=ddxx+ddyy
       dstbdy=sqrt(ddd)
       nwdist=time*v_ent 
        
!=====================================================================
! CALCOLO NUOVE COORDINATE DEL BODY (nxbdy,nybdy)
!=====================================================================

! calcolo della distanza percorsa in 1 tstep (nwdist)
! se il body arriva esattamente sul lato a ttime=0 (1 iflogico)?
! se si esce dall''elemento (2 iflogico) allora ho coordinate del body sul lato di uscita
! e in piu un deltat > 0 da utilizzare durnate questo stesso tstep.
! se si rimane nell''elemento (3 iflogico) allora calcolo le nuove coordinate 
! (pnt_inside) e passo al tstep successivo

! Se il body finisce a 0.1 m dalla linea decido si farlo passare al
! elemento ad una distanza di 0.1 m. Questa è la precisione da me 
! imposta
        gf=dvert(1,ie)/100000 
        if(abs(nwdist-dstbdy).lt.gf)then      ! body dista meno di 0.1m dal lato
!          write(lunit,*) 'WARNING!! BODY SU LATO ELEMENTO'
!          write(lunit,*) 'track orig',bdy
          nwdist=dstbdy+2*gf
        endif

	 zn0=zbdy
         l0=lybdy
	 call getalfa(l_int,in_d,in_dm,in_dx,a_int)
	 call getalfa(l_out,ou_d,ou_dm,ou_dx,a_out)
         call getzvel(ie,zn0,l0,l_int,l_out,a_int,a_out,w)
	 if( blgrsurf ) w = 0.
	 lmax = nlvdi
         call lagr_layer_thickness(ie,lmax,hl)
         layd=hl(l0)

! check for new body position 

        if(nwdist.gt.dstbdy)then            !2 body su nuovo el
	 nxbdy=itrx
	 nybdy=itry
	 distance=nwdist-dstbdy
	 if( v_ent .lt. epsggu ) then
	   time = 0.
           xbdy=nxbdy ! nuove coordinate del body
           ybdy=nybdy ! nuove coordinate del body
	   return
	 end if
	   
	 deltat=distance/v_ent ! time within the element 

! VERTICAL TREATMENT

         call vertpos(zn0,deltat,layd,w,zn1,ztime,addl) !compute the new vert. pos. and time

          if(ztime.gt.deltat)then   ! change element keeping the same layer
	   l1=l0
	   zbdy=zn1
           lybdy=l1
	     
	  elseif(ztime.lt.deltat)then ! change layer keeping the same element 
	   l1=l0+addl
           if(l1.le.1)stop 'surface'
	   if(l1.ge.lb)stop 'bottom'
	   zbdy=zn1
	   lybdy=l1
	   nwdist=ztime*v_ent
	   call pnt_inside(nwdist,dstbdy,fay,fax,xbdy,ybdy,itrx,itry,nxbdy,nybdy)
           deltat=time-ztime
           time=deltat
           ie=ie ! individuazione vecchio elemento
           xbdy=nxbdy ! nuove coordinate del body
           ybdy=nybdy ! nuove coordinate del body
           return
	  end if
	 
         newie=ieltv(l_out,ie)

	 do i=1,3
          if(ie.eq.ieltv(i,newie))then
           ltbdy=i ! individuazione lato di entrata del body
          endif
         end do
	 xbdy=nxbdy ! nuove coordinate del body
         ybdy=nybdy ! nuove coordinate del body
         ie=newie ! individuazione nuovo elemento
         time=deltat

         if(newie.eq.-1)then			!2.a body e uscito dal dominio
	  if( .not. bback ) then
!            write(lunit,*) 'STOP!! BODY ',bdy,' USCITO '
!            write(lunit,*) 'ELEMENTO USCITA ',ie
	  end if
          ie=-ie ! flag per skip bdy dal calcolo
	  time=0.
          xbdy=nxbdy ! ultime coordinate del body uscito
          ybdy=nybdy ! ultime coordinate del body uscito
          return
         endif


  	elseif(nwdist.lt.dstbdy)then            !3 body su vecchio el
  	  
	 deltat = time
	 call vertpos(zn0,deltat,layd,w,zn1,ztime,addl) !compute the new vert. pos. and time

          if(ztime.gt.deltat)then   ! keeping the same layer
           l1=l0
           zbdy=zn1
           lybdy=l1
	   call pnt_inside(nwdist,dstbdy,fay,fax,xbdy,ybdy,itrx,itry,nxbdy,nybdy)
           deltat=0.
           time=deltat
           ie=ie ! individuazione vecchio elemento
           xbdy=nxbdy ! nuove coordinate del body
           ybdy=nybdy ! nuove coordinate del body
          elseif(ztime.lt.deltat)then ! change layer keeping the same element 
           l1=l0+addl
           if(l1.le.1)stop 'surface'
           if(l1.ge.lb)stop 'bottom'
           zbdy=zn1
           lybdy=l1
           nwdist=ztime*v_ent
           call pnt_inside(nwdist,dstbdy,fay,fax,xbdy,ybdy,itrx,itry,nxbdy,nybdy)
           deltat=time-ztime
           time=deltat
           ie=ie ! individuazione vecchio elemento
           xbdy=nxbdy ! nuove coordinate del body
           ybdy=nybdy ! nuove coordinate del body
           return
          end if

        endif
        end

!**********************************************************************

	subroutine track_line(time,bdy,ie,xbdy,ybdy,zbdy,lybdy,ltbdy)

! in questa subroutine si calcola il percorso del body partendo
! da un punto sul lato dell''elemento ie
! questa subroutine e'' chiamata ogni qualvolta un body
! nello stesso timestep arriva in un nuovo elemento

	use lagrange_data
	use lagrange_vertical
	use lagrange_util_tr
	use geom
	use hydro_vel
	use levels
	use basin
        use regular

        implicit none
        
	include 'param.h' 

	double precision time	!total time to travel
	integer bdy	!number of body
	integer ie	!number of element of body
	double precision xbdy,ybdy	!horizontal position of body
	double precision zbdy	!relative vertical position of body in layer
	integer lybdy	!layer of body
	integer ltbdy	!entering side of body [1-3]
 
	double precision epsggu
	parameter (epsggu = 1.e-7)

        double precision deltat ! frazione di time step body si muove in ie
 




	double precision v_ent ! valore mediato tra velocita' int e out
	double precision v_int ! modulo della velocita' di entrata
	double precision v_out ! modulo della velocita' di uscita
	integer l_int !lato di entrata numerazione elemento (var. input)
	integer l_out ! lato di uscita numerazione elemento
	double precision nxbdy,nybdy !nuove coordinate del body
	integer newie ! nuovo elemento che contiene il body

! vertical treatment 
        double precision z0,z1,zn0,zn1  !      cuccomod
        integer l0,l1,addl,lmax
        double precision ztime,lb,subtime,layd
        double precision w
        double precision hl(nlvdi)
        double precision a_int,a_out
	double precision in_d,ou_d
        integer in_dm,in_dx,ou_dm,ou_dx

        integer idum ! seme f(x) random
                        

! variabili di servizio
	
	integer pb(3),ip,i,exi(2),ipb(2),exio(2)
	double precision nwdist,distance ! funzione per il calcolo del nwdist
        double precision dxx,dyy,ddxx,ddyy,ddd
        double precision itrx,itry
        double precision drd(3)
        double precision fsum,fpb(2),ffpb(2)
	double precision tcomp ! funzione per il calcolo del deltat
  	double precision icy,icx,ib ! parametri retta del lato da cui entra il body
	double precision ocy,ocx,ob ! parametri retta del lato da cui esce il body
        double precision fay,fax ! parametri fascio rette traiettoria
        double precision db ! parametri retta della traiettoria sguita dal body 
	double precision dstbdy ! distanza tra lato entrata e uscita lungo la traiettoria
        double precision bdx ! distanza (frazione) del body da estremo piu vicno retta entrata
	integer pbdx,gbdx ! numero interno dell'estremo piu vicino e piu lontano dal body
        integer ps,ng,nl,i1,i2
        double precision gf !precision along the side
! inizializzazione parametri
        data idum/947437/
        save idum
        
        fsum=0
	nl = 0
	ng = 0
	ps = 0
	lb = lybdy

	lmax = ilhv(ie)

!==========================================================
! CALCOLO TRAIETTORIE
!==========================================================

! calcolo retta del lato da cui entra il body
! exi(1,2) sono i puntatori del nen3v per gli estremi della retta 
! su cui si trova il body

        
        l_int=ltbdy


        i1=mod(l_int,3)+1
        i2=mod(i1,3)+1
        exi(1)=i1
        exi(2)=i2
                                          
	call retta(exi,icy,icx,ib,ie)

! velocità di entrata v_int 

	v_int=vel_ie(l_int,ie)

        if(v_int.le.0)then
!         write(lunit,*) 'STOP!! ERRORE FLUSSO DI ENTRATA NEGATIVO'
!         write(lunit,*) 'track line'
	 !stop
	 time=0
         ie = -ie
	 return
        endif

! determinazione dei lati e dei flussi di uscita probabili pb(1,2),vpb(1,2)
! pb(mi dice il numero interno del lato)
        
        ip=1
        do i=1,3
 	 pb(i)=0
         if(l_int.ne.i)then
 	  pb(i)=ip
	  ipb(ip)=i
          fpb(ip)=vel_ie(i,ie) 
          ip=ip+1  
         endif 
	end do

        v_ent=sqrt((ulnv(1,ie)**2)+(vlnv(1,ie)**2))
        

	call distp(xbdy,ybdy,exi,bdx,pbdx,gbdx,ie)

        in_d=bdx
        in_dm=pbdx
        in_dx=gbdx
         
!______________________________________________________________________
! CASO [A] flussi ai lati opposti al l_int entrambi <0
! si calcola la frazione sul lato di entrata di competenza 
! del lato pb(1) e del lato pb(2) per la possibile uscita del body.

      	if((fpb(1).lt.0).and.(fpb(2).lt.0))then

                
        fsum=abs(fpb(1))+abs(fpb(2))
         do i=1,2
	  ffpb(i)=abs(fpb(i)/fsum)
	 end do

         
! calcola la distanza del body (bdx) dall''estremo piu vicino
! in termini di frazione sulla retta totale  e
! dice qual''e'' il numero interno dell''estremo piu vicino (pbdx)
! e quello piu lontano (gbdx)

        call distp(xbdy,ybdy,exi,bdx,pbdx,gbdx,ie)

	in_d=bdx
        in_dm=pbdx
        in_dx=gbdx


        if(ie.le.0)then
         call find_elem_from_old(-ie,xbdy,ybdy,ie)
         time=0
         return
        end if
                                                   
        
! individuo il lato da cui esce il body mediante confronto
! tra posizione del body rispetto alle frazioni  

	if(pbdx.lt.1.or.pbdx.gt.3.or.pb(pbdx).le.0) then
!         write(lunit,*) 'STOP!! ERRORE INTERNO (3)'
!	 write(lunit,*) 'CASO [A] TRACK_LINE'	 
	 time=0
         ie = -ie
	 return
	end if

	 if(bdx.le.ffpb(pb(pbdx)))then	  
	  l_out=pbdx
	 elseif(bdx.gt.ffpb(pb(pbdx)))then
	  l_out=gbdx
	 end if

         
! estremi del lato da cui esce il body

        i1=mod(l_out,3)+1
        i2=mod(i1,3)+1

        exio(1)=i1
        exio(2)=i2

                                          
! equazione della retta del lato di uscita

        call retta(exio,ocy,ocx,ob,ie)


! calcolo del fascio di traiettorie del body, mi serve sapere:
! (l_int) nodo opposto al lato di entrata
! coordinate del punto di separazione (sep) tra le 2 frazioni sul lato
! di entrata (ad estremo 1 sommo ffpd(2))
         
        call dirz(l_int,ffpb,exi,icy,icx,ib,fay,fax,ie,pb,ipb)
        
! determinazione retta traiettoria

        call traj(fay,fax,xbdy,ybdy,db)

! deteminazione coordinate punto intercetto retta traiettoria-
! retta di uscita
         
        call interc(fay,fax,db,ocy,ocx,ob,itrx,itry)

	call distp(itrx,itry,exio,bdx,pbdx,gbdx,ie)

	ou_d=bdx
        ou_dm=pbdx
        ou_dx=gbdx	

!_____________________________________________________________________________
! CASO [B] flussi ai lati opposti al l_int uno >0 e altro <0
! in questo caso il lato di uscita e'' sicuramente quello dove ho la velocita''
! negativa, non serve quindi il calcolo delle rette di partenza e di arrivo
! gli estremi della retta del lato di uscita sono i due estremi non opposti
! alla retta che ha valore di vel negativo

        else if(((fpb(1).lt.0).and.(fpb(2).ge.0)).or.((fpb(1).ge.0).and.(fpb(2).lt.0)))then
        

         do i=1,2
          if(fpb(i).lt.0)then
           l_out=ipb(i)
          endif
         enddo

! estremi del lato da cui esce il body
	
        i1=mod(l_out,3)+1
        i2=mod(i1,3)+1

        exio(1)=i1
        exio(2)=i2
                                          
                                          
! equazione della retta del lato di uscita

        call retta(exio,ocy,ocx,ob,ie)
        
        
! individuo i 2 flussi di entrata del body nell''elemento per il calcolo
! delle frazioni di competenza sul lato di uscita.
 
         ip=1
         do i=1,3
 	  pb(i)=0
          if(l_out.ne.i)then
           pb(i)=ip
           ipb(ip)=i
           fpb(ip)=vel_ie(i,ie)
           ip=ip+1
          endif
         end do

         
! calcolo delle frazioni sul lato di uscita per determinare 
! la retta traiettoria del percorso del body
	
         fsum=abs(fpb(1))+abs(fpb(2))
         do i=1,2
          ffpb(i)=abs(fpb(i)/fsum)
         end do


! calcolo del fascio di traiettorie del body, mi serve sapere:
! l_out nodo opposto al lato di uscita
! coordinate del punto di separazione tra le 2 frazioni sul lato 
! di uscita (come prima ad estremo 1 sommo ffpb(2))

         call dirz(l_out,ffpb,exio,ocy,ocx,ob,fay,fax,ie,pb,ipb)
         
! calcolo retta traiettoria 

         call traj(fay,fax,xbdy,ybdy,db)
         
! deteminazione coordinate punto intercetto retta traiettoria-
! retta di uscita

         call interc(fay,fax,db,ocy,ocx,ob,itrx,itry)

	call distp(itrx,itry,exio,bdx,pbdx,gbdx,ie)

        ou_d=bdx
        ou_dm=pbdx
        ou_dx=gbdx     

! calcolo velocita'' in uscita v_out per il calcolo del modulo della
! velocita di percorso v_ent pari alla media v_int, v_out

         v_out=vel_ie(l_out,ie)
        
         if(v_out.gt.0)then
!          write(lunit,*) 'STOP!! ERRORE VELOCITÀ DI USCITA POSITIVA'
!          write(lunit,*) 'track line'
	  !stop
	  time=0
          ie = -ie
	  return
         elseif(v_out.eq.0)then
!          write(lunit,*) 'STOP!! ERRORE VELOCITÀ DI USCITA NULLA'
!          write(lunit,*) 'track line'
	  !stop
	  time=0
          ie = -ie
	  return
         endif

        elseif((fpb(1).ge.0).and.(fpb(2).ge.0))then
!________________________________________________________________________________
! CASO [C] tutti i flussi sono entranti o nulli. Errore il campo di moto è
! convergente, problemi con l''eliminazione delle componenti DIVERGENTI
! FERMARE IL CILO

!        write(lunit,*) 'STOP!! FLUSSI TUTTI POSITIVI IN ',ie
!        write(lunit,*) 'CASO [C] TRACK_LINE'
         time=0
!        ie = -ie
         return
        elseif((fpb(1).lt.0).and.(fpb(2).lt.0))then
!        write(lunit,*) 'STOP!! FLUSSI TUTTI NEGATIVI IN ',ie
!        write(lunit,*) 'CASO [D] TRACK_LINE'
         time=0
!        ie = -ie
         return
        elseif((nl.eq.2).and.(ps.eq.1))then
!        write(lunit,*) 'STOP!! UN SOLO FLUSSO POSITIVO E 2 NULLI IN ',ie
!        write(lunit,*) 'CASO [E] TRACK_LINE'
!        write(lunit,*) '[elemento con due lati C.B]'
         time=0
!        ie = -ie
         return
        elseif((nl.eq.2).and.(ng.eq.1))then
!        write(lunit,*) 'STOP!! UN SOLO FLUSSO NEGATIVO E 2 NULLI IN ',ie
!        write(lunit,*) 'CASO [F] TRACK_LINE'
!        write(lunit,*) '[elemento con due lati C.B]'
         time=0
!        ie = -ie
         return
        elseif((nl.eq.1).and.(ps.eq.2))then
!        write(lunit,*) 'STOP!! 2 FLUSSI POSITIVI E 1 NULLO IN ',ie
!        write(lunit,*) 'CASO [G] TRACK_LINE'
         time=0
!        ie = -ie
         return
        elseif((nl.eq.1).and.(ng.eq.2))then
!        write(lunit,*) 'STOP!! 2 FLUSSI NEGATIVI E 1 NULLO IN ',ie
!        write(lunit,*) 'CASO [H] TRACK_LINE'
         time=0
!        ie = -ie
         return
        endif
!_______________________________________________________________________

! calcolo distanza massima percorribile all''interno di elemento

        dxx=itrx-xbdy
        dyy=itry-ybdy
        ddxx=(dxx)**2
        ddyy=(dyy)**2
        ddd=ddxx+ddyy
        dstbdy=sqrt(ddd)

	nwdist=time*v_ent

	zn0=zbdy
        l0=lybdy
	call getalfa(l_int,in_d,in_dm,in_dx,a_int)
	call getalfa(l_out,ou_d,ou_dm,ou_dx,a_out)
        call getzvel(ie,zn0,l0,l_int,l_out,a_int,a_out,w)
	if( blgrsurf ) w = 0.
	lmax = nlvdi
        call lagr_layer_thickness(ie,lmax,hl)
        layd=hl(l0)

!=====================================================================
! CALCOLO NUOVE COORDINATE DEL BODY (nxbdy,nybdy)
!=====================================================================

! il body entra sempre da un lato con una velocita'' positiva
! calcolo della distanza percorsa in 1 tstep (nwdist)
! se il body arriva esattamente sul lato a ttime=0 (1 iflogico)?
! se si esce dall''elemento (2 iflogico) allora ho coordinate del body sul lato di uscita
! e in piu un deltat > 0 da utilizzare durnate questo stesso tstep.
! se si rimane nell''elemento (3 iflogico) allora calcolo le nuove coordinate 
! (pnt_inside) e passo al tstep successivo

! Se il body finisce a 0.1 m dalla linea decido si farlo passare al
! elemento ad una distanza di 0.1 m. Questa è la precisione da me
! imposta

        gf=dvert(1,ie)/100000    
        if(abs(nwdist-dstbdy).lt.gf)then     ! body dista meno di 0.1m dal lato    
!          write(lunit,*) 'STOP!! BODY SU LATO ELEMENTO'
!          write(lunit,*) 'track line',bdy
          nwdist=dstbdy+2*gf
        endif

        if(nwdist.gt.dstbdy)then            !2 body su nuovo el
	 nxbdy=itrx
	 nybdy=itry
	 distance=nwdist-dstbdy
	 if( v_ent .lt. epsggu ) then
	   time = 0.
           xbdy=nxbdy ! nuove coordinate del body
           ybdy=nybdy ! nuove coordinate del body
	   return
	 end if
         deltat=distance/v_ent

	call vertpos(zn0,deltat,layd,w,zn1,ztime,addl) !compute the new vert. pos. and time

          if(ztime.gt.deltat)then   ! change element keeping the same layer
           l1=l0
           zbdy=zn1
           lybdy=l1
        
          elseif(ztime.lt.deltat)then ! change layer keeping the same element 
           l1=l0+addl
           if(l1.lt.1)stop 'error stop track_line: surface'
           if(l1.gt.lmax)stop 'error stop track_line: bottom'
           zbdy=zn1
           lybdy=l1
           nwdist=ztime*v_ent
           call pnt_inside(nwdist,dstbdy,fay,fax,xbdy,ybdy,itrx,itry,nxbdy,nybdy)
           deltat=time-ztime
           time=deltat
           ie=ie ! individuazione vecchio elemento
           xbdy=nxbdy ! nuove coordinate del body
           ybdy=nybdy ! nuove coordinate del body
           return
          end if


         newie=ieltv(l_out,ie) !
	 do i=1,3
          if(ie.eq.ieltv(i,newie))then          !3 body su vecchio el
           ltbdy=i ! individuazione lato di entrata del body
          endif
         end do
         xbdy=nxbdy ! nuove coordinate del body
         ybdy=nybdy ! nuove coordinate del body
        ie=newie ! individuazione nuovo elemento
        time=deltat



         if(newie.eq.-1)then			!2.a body uscito dal dominio
	  if( .not. bback ) then
!            write(lunit,*) 'STOP!! BODY ',bdy,' USCITO '
!            write(lunit,*) 'ELEMENTO USCITA ',ie
	  end if
          ie=-ie! flag per skip bdy dal calcolo
          xbdy=nxbdy ! ultime coordinate del body uscito
          ybdy=nybdy ! ultime coordinate del body uscito
	  time=0.
	  return     
         endif

       elseif(nwdist.lt.dstbdy)then   !3 body su vecchio el

	 deltat = time
         call vertpos(zn0,deltat,layd,w,zn1,ztime,addl) !compute the new vert. pos. and time

          if(ztime.gt.deltat)then   ! keeping the same layer
           l1=l0
           zbdy=zn1
           lybdy=l1
           call pnt_inside(nwdist,dstbdy,fay,fax,xbdy,ybdy,itrx,itry,nxbdy,nybdy)
           deltat=0.
           time=deltat
           ie=ie ! individuazione vecchio elemento
           xbdy=nxbdy ! nuove coordinate del body
           ybdy=nybdy ! nuove coordinate del body
          elseif(ztime.lt.deltat)then ! change layer keeping the same element 
           l1=l0+addl
           if(l1.le.1)stop 'surface'
           if(l1.ge.lb)stop 'bottom'
           zbdy=zn1
           lybdy=l1
           nwdist=ztime*v_ent
           call pnt_inside(nwdist,dstbdy,fay,fax,xbdy,ybdy,itrx,itry,nxbdy,nybdy)
           deltat=time-ztime
           time=deltat
           ie=ie ! individuazione vecchio elemento
           xbdy=nxbdy ! nuove coordinate del body
           ybdy=nybdy ! nuove coordinate del body
           return
          end if

       endif
       end

!**********************************************************************

!-------------------------------------------------------------------------
        end module lagrange_track
!-------------------------------------------------------------------------
