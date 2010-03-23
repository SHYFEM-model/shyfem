c
c $Id: reactor.f,v 1.8 2003/03/25 14:08:55 georg Exp $
c
c biological reactor from OGS (Cosimo)
c
c revision log :
c
c 13.08.1998	ggu&cs	written from scratch (source Cosimo)
c 27.08.1998	ggu	light included
c 28.08.1998	ggu	new boundary condition handeling
c 07.09.1998	cs&gc	changed to micromols
c 08.09.1998	ggu	some cleaning
c 24.11.1998	ggu	use switch ibio to see if reactor has to be run
c 22.01.1999	ggu	save some not yet saved variables in biocos
c
c****************************************************
c
c	program test
c	call ltest
c	call rtest
c	end
c
c****************************************************

	subroutine ltest

	real luxsup

	do i=0,24*365
	  h = i
	  h = mod(h,24.)
	  id = 1 + i/24
	  rl = luxsup(id,h)
	  write(6,*) i,h,id,rl
	end do

	end

c****************************************************

	subroutine rtest

	implicit none

	include 'reactor.h'

	integer ndim
	parameter(ndim=9)
	integer nlvdim,nkndim
	parameter(nlvdim=1,nkndim=1)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	real hdknv(nlvdim,nkndim)			!thickness of layer
c	common /hdknv/hdknv
	integer ilhkv(nkndim)
c	common /ilhkv/ilhkv

	real e(ndim)
	real einit(ndim),ebound(ndim)
	real temp,salt
	real light(nlvdim,nkndim)
	integer level

	integer it,idt,itend
	integer i
	real dt,t
	real lux
	real oxysat

c------------------------------------------------------------------
c                                             poc  pop 
c                     po4   phy1 bac   zoo    detc detp  doc   dop  oxy
c	data einit  / 0.5,  0.1, 0.05, 0.005, 0.,  0.0,  0.01, 0.1, 8.  /!run1
c	data ebound / 0.8,  0.0, 0.00, 0.000, 0.,  0.0,  0.05, 0.5, 8.  /!run1
c	data einit  / 0.05, 0.8, 0.2,  0.04,  0.,  0.00, 0.8,  3.2, 250./!run3
c	data ebound / 4.6,  0.0, 0.00, 0.00,  0.,  4.8,  0.0,  0.4, 250./!run3
	data einit  / 0.05, 0.8, 0.2,  0.04,  0.,  0.03, 0.8,  0.0, 250./!run4
	data ebound / 4.6,  0.0, 0.00, 0.00,  0.,  4.8,  0.0,  0.4, 250./!run4
c                     po4   phy1 bac   zoo    detc detp  doc   dop  oxy
c                                             poc  pop 
c------------------------------------------------------------------

	idt = 300
	itend = 2139030
	itend = 86400
	itend = 5184000
	dt = idt
	level = 1

	nkn = 1
	hdknv(1,1) = 2.
	ilhkv(1) = 1

	temp = 9.
	salt = 35.

	do i=1,ndim
	  e(i) = einit(i)
	end do

	e(9) = oxysat(temp,salt)

	it = 0.
	t = 0.
	write(6,'(12f12.4)') t,e,temp,salt
	do it=idt,itend,idt
	  t = it/86400.
	  call setlux(it,light,e,hdknv,nlvdim,nkndim,ndim)
	  lux = light(1,1)
          call reactor(it,dt,level,temp,salt,e,lux)
	  write(6,'(12f12.4)') t,e,temp,salt
	end do

	end

c****************************************************

	function oxysat(temp,salt)

c computes oxygen level at saturation

	implicit none

	real oxysat
	real temp,salt

        oxysat = (14.6244-0.367134*temp+4.4972E-3*temp**2
     :     -0.0966*salt+0.00205*salt*temp+2.739E-4*salt**2)
     :    /32.

	end

c****************************************************

	function luxsup(id,h)

c computes light at surface at day id and hour ih

	implicit none

	real luxsup	!light intensity [lux]
	integer id	!day of year [1-365]
	real h		!hour of day [0-24]

	real pi,rad
	parameter( pi = 3.14159 , rad = 2. * pi / 365. )

	real cosinus,luxnoon,period
	real hmin,hmax

	cosinus = cos( rad * ( id + 10 ) )
	luxnoon = 50000. - 14909. * cosinus
	period = 24. * ( 0.5 - 0.125 * cosinus )

	hmin = 12. - 0.5 * period
	hmax = 12. + 0.5 * period

	if( h .le. hmin ) then
	  luxsup = 0.
	else if( h .ge. hmax ) then
	  luxsup = 0.
	else
	  luxsup = luxnoon * cos( pi * (h-12.) / period )
	end if

c	write(6,*) id,h
c	write(6,*) cosinus,luxnoon,period
c	write(6,*) hmin,hmax,luxsup

	end

c****************************************************

	subroutine setlux(it,light,e,hdkn,nlvdim,nkndim,ndim)

c sets light intensity in [lux]

	implicit none

c	include 'compar.inc'
	include 'reactor.h'

	integer it		!time in seconds
	integer nlvdim,nkndim,ndim
	real light(nlvdim,nkndim)
	real e(nlvdim,nkndim,ndim)
	real hdkn(nlvdim,nkndim)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	integer ilhkv(1)
	common /ilhkv/ilhkv

	integer id,ih
	integer k,l,lmax
	real light0,l0
	real fito,pom,h,aux

	real luxsup

	id = 1 + mod(it/86400,365)
	ih = 1 + mod(it/3600,24)
	h = ( it - (id-1) * 86400 ) / 3600.

	if( h .lt. 0. .or. h .gt. 24. ) then
	  write(6,*) 'error in computing hour: '
	  write(6,*) it,id,ih,h
	  stop 'error stop setlux'
	end if

	light0 = luxsup(id,h)

	do k=1,nkn
	  l0 = light0
	  lmax = ilhkv(k)
	  do l=1,lmax
	    fito = e(l,k,2)
	    pom = e(l,k,5)
	    h = hdkn(l,k)
	    aux = kestf * fito + kestd * pom + kestw * h
c
c	sbagliato -> tutto dev'essere moltiplicato con h	!FIXME
c
c meglio:	aux = kestf * fito + kestd * pom + kestw
c		eaux = exp(-aux*h)
c		lmed = (1.-eaux)*l0/(aux*h)	!average available light
c		l0 = l0*eaux			!light at bottom of layer
c
	    aux = 0.5 * aux
	    light(l,k) = l0 * exp ( -aux )	!light halfway down
	    l0 = light(l,k) * exp ( -aux )	!light at bottom of layer
	  end do
	end do
	 
	end

c****************************************************

      subroutine reactor(itime,dts,level,temp,salt,e,lux)
c
c itime		time in seconds
c dts		time step (seconds)
c level		level (1 -> surface, >1 -> below surface)
c temp		temperature
c salt		salinity
c e		state variables
c lux		light in lux

	integer itime
	real dts
	integer level
	real temp,salt
	real e(1)
	real lux

      real       fsedc,fsedp,cp1,ct1,ct2,lass1,lass2,arrt,
     :           kfl1,kfl2,light1,kdefo2
      real       tsf1,tsf2,alfa1,alfa2,kmff1,kmff2,krff1,
     :           kmzz,kgrz1,kgrz2,kescz,knitn,kdecd,kdecn,kdecp
      real	 a2rpcd


	include 'reactor.h'

        real x(nvt),dx(nvt), xold(nvt)

	integer ihdim
	parameter ( ihdim = 366*24 )

	integer icall
	save icall
	data icall /0/

c initialization

	dt = dts / 3600.

C...................................................................
C..............SINKING..............................................
C...................................................................
C         calcola la sedimentazione del detrito,
C         used e' la velocita' di sedimentazione del detrito
C
cJ           fsedc=e(detC)*USED*dt!sedimentazione del detrito carbonioso
cJ           fsedp=e(detp)*USED*dt!sedimentazione del detrito carbonioso
c
C           e(detC)=e(detC)-fsedc
C           e(detp)=e(detp)-fsedp
C           e(sedc)=e(sedc)+fsedc
c           e(sedp)=e(sedp)+fsedp
c...................................................................
c inizializzazione
c ..............................................................

         do l=1,nvd
	  xold(l)=e(l) ! valori correnti delle nvd variabili che diffondono
         end do

c ..................................................................
C     COMPUTE RATE CONSTANT VARIATION
C     WITH TEMPERATURE

        CP1=MIN(TMAX1,temp) !questo serve per azzerare la crescita a valori di T>Tmax

c	write(6,*) cp1,tmax1,temp,TOTT1,COEFF1

        CT1=EXP(COEFF1*(CP1-TOTT1))
        CT2=(((TMAX1-CP1)/(TMAX1-TOTT1))**(COEFF1*(TMAX1-TOTT1)))
        lass1=ct1*ct2         !lassiter per fito1
        CP1=MIN(TMAX2,temp)
        CT1=EXP(COEFF2*(CP1-TOTT2))
        CT2=(((TMAX2-CP1)/(TMAX2-TOTT2))**(COEFF2*(TMAX2-TOTT2)))
        lass2=ct1*ct2          !lassiter per  bacteria
        arrt=(1.07**(temp-20))  !arrhenius per le cinetiche normali

C     COMPUTE LIGHT ATTENUATION
C     AND SELF-SHADING (la luce viene da lux(istep)=luce/50000
C
c	itsec = mod(itime,86400)
c	ihour = 1 +  itime / 3600
c        fluce=lux(ihour)*50000

	fluce = lux
	light1=fluce/Iopt1
	light1=light1*exp(1-light1)
C
C      COMPUTE OXYGEN DEFICIT
C
	if( level .eq. 1 ) then
          KDEFO2=(14.6244-0.367134*temp+4.4972E-3*temp**2-
     :    0.0966*salt+0.00205*salt*temp+2.739E-4*salt**2)-xold(oxy)
	else
	  KDEFO2 = 0.
	end if
C
C
C     PRELIMINARY COMPUTATIONS FOR
C     BIOCHEMICAL REACTIONS
C
c    pieces per fito 1
c
        TSF1=mumax1*LIGHT1*lass1*xold(po4)/(ksp1+xold(po4))
        ALFA1=TSF1*xold(phy1)
        KMFF1=kmf1*arrt*xold(phy1)! va a pom
        KRFF1=krf1*arrt*xold(phy1)! va a dom
c
c     pieces per  bacteri
C
      TSF2=mumax2*lass2*(xold(po4)+xold(dop))/
     :  (ksp2+xold(po4)+xold(dop))*xold(doc)/(xold(doc)+ksdoc2)
        ALFA2=TSF2*xold(bac)
        KMFF2=kmf2*xold(bac)*arrt !include anche resp.
c
c     pieces per zooplancton
c
        KMZZ=kmz*arrt*xold(zoo)   ! va pom
        KESCZ=xold(zoo)*arrt*kexz !escrezione a dom
c
c   grazing tipo II
c
c        Kgrz1=kgr1*xold(zoo)*xold(phy1)/(kfz1+xold(phy1)
c     :   +xold(bac)*alfa)! grazing su fito1
c        Kgrz2=kgr1*xold(zoo)*xold(bac)*alfa/(kfz1+xold(phy1)+
c     :   xold(bac)*alfa)! grazing su bacteri
c
c    grazing tipo III 
c
        Kgrz1=kgr1*xold(zoo)*xold(phy1)*xold(phy1)
     :  /(kfz1+xold(phy1)*xold(phy1)+xold(bac)*xold(bac)*alfa*alfa)! grazing su fito1
c	Kgrz1 = 0.	!ggu
c
        Kgrz2=kgr1*xold(zoo)*xold(bac)*xold(bac)*alfa*alfa
     :  /(kfz1+xold(phy1)*xold(phy1)+xold(bac)*xold(bac)*alfa*alfa)! grazing su bact. 
c
c     pieces other 
c
        kdecd=kdcd*arrt ! pom a dom
        kdecp=kdcp*arrt ! pom a dom
c
c        rpcd=xold(dop)/xold(doc)
c        a2rpcd=alfa2*rpcd
c
      a2rpcd=mumax2*lass2*(xold(po4)+xold(dop))/
     :   (ksp2+xold(po4)+xold(dop))
     :   *xold(dop)/(xold(doc)+ksdoc2)*xold(bac)
c
c     CALCOLA LE DERIVATE
c
c
c     efettua controllo
c
      dx(dop)=-a2rpcd+krff1*rpc1+kescz*rpcz+
     : kdecp*(xold(detp))
     : +(eff2*kgrz2)*(rpc2-rpcz)
c
       if((xold(dop)+dx(dop)).le.0)then
       xrid=xold(dop)/dx(dop)
       alfa2=alfa2*xrid
       a2rpcd=a2rpcd*xrid
       end if
c
       dx(po4)=-alfa1*rpc1-alfa2*rpc2+a2rpcd+kmff2*rpc2
c
       if ((xold(po4)+dx(po4)).le.0) then
        xrid=xold(po4)/dx(po4)
        alfa1=0.
        alfa2=min(alfa2,(a2rpcd/rpc2))
       dx(po4)=-alfa1*rpc1-alfa2*rpc2+a2rpcd+kmff2*rpc2
       end if
c
c     PHOSPHATES (mg/l P ) (16. = rapp. N/P= RNP)
C     DISCIOLTO ORGANICO (DOC E DOP)
      dx(doc)=-alfa2+krff1+kescz+kdecd*(xold(detc))
      dx(dop)=-a2rpcd+krff1*rpc1+kescz*rpcz+
     : kdecp*(xold(detp))
     : +(eff2*kgrz2)*(rpc2-rpcz)
c
c     PHYTOPLANKTON   (carbonio)
c
        dx(phy1)=alfa1-kmff1-krff1-kgrz1
c
c     bacteri   (carbonio)
c
        dx(bac)=alfa2-kmff2-kgrz2
c
c     ZOOPLANKTON  (carbonio)
c
        dx(zoo)=eff1*kgrz1+eff2*kgrz2-kmzz-kescz
c
c
c     DETRITO CARBONIOSO (POM e POP)
c
      dx(detc)=-kdecd*xold(detC)+(1-eff1)*kgrz1+(1-eff2)*kgrz2+
     :      kmzz+kmff1
      dx(detp)=-xold(detp)*kdecp+(kmff1+(1-eff1)*kgrz1)*
     :rpc1+(1-eff2)*kgrz2*rpc2+kmzz*rpcz
c
C
C     OXYGEN
       dx(oxy)=krear*kdefo2+roc1*(alfa1-krff1)-roc2*kmff2
     : -rocd*kdecd*(xold(detC))   
c........................................................................
c....  fine delle derivate. passo ad integrare ..........................
c........................................................................
c,
         do l=1,nvd
	  e(l)=xold(l)+dx(l)*dt! valori correnti delle nvd variabili che diffondono
          if(e(l).lt.0.) then
            e(l)=0.0
          end if
         end do
C
      RETURN
      END
C


c***********************************************

	subroutine biocos(it,idt)

c eco-model cosimo

	implicit none

	include 'param.h'

	integer it	!time in seconds
	integer idt	!time step in seconds

	integer ndim,ntot
	parameter( ndim = 9 , ntot = 3 )

	real e(nlvdim,nkndim,ndim)	!state vector
	real light(nlvdim,nkndim)	!light intensity

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer nlvdi,nlv
        common /level/ nlvdi,nlv

	real saltv(nlvdim,1), tempv(nlvdim,1)
        common /saltv/saltv
        common /tempv/tempv
	integer ilhkv(1)
	common /ilhkv/ilhkv
	real hdknv(nlvdim,nkndim)		!thickness of layer
	common /hdknv/hdknv

        real saux1(nlvdim,1),saux2(nlvdim,1)
        real saux3(nlvdim,1),saux4(nlvdim,1)
        common /saux1/saux1, /saux2/saux2
        common /saux3/saux3, /saux4/saux4
	real zeov(3,1)
	common /zeov/zeov
	real zenv(3,1)
	common /zenv/zenv
	real rsv(1)
	common /rsv/rsv
	real difv(0:nlvdim,1)
	common /difv/difv

	integer k,i,l,lmax
	integer ibio
	real t,s,dt
	real ee(ndim)
	real einit(ndim)
	real ebound(ndim)
	integer icall,iunit
	integer itt,idtt,j
	real lux
	real dtt
	real difmol
	real rkpar,azpar,adpar,aapar
	real rmass0,rmass1,pmass0,pmass1
	integer istot,isact
	real oxysat
	real getpar
	integer iround

	save rkpar,azpar,adpar,aapar,istot
	save difmol

	save einit,ebound,icall
c------------------------------------------------------------------
c                                             poc  pop 
c                     po4   phy1 bac   zoo    detc detp  doc   dop  oxy
c	data einit  / 0.5,  0.1, 0.05, 0.005, 0.,  0.0,  0.01, 0.1, 8. /!cosimo
c	data ebound / 0.8,  0.0, 0.00, 0.000, 0.,  0.0,  0.05, 0.5, 8. /!cosimo
	data einit  / 0.05, 0.8, 0.2,  0.04,  0.,  0.03, 0.8,  0.0, 250./!guido
	data ebound / 4.6,  0.0, 0.00, 0.00,  0.,  4.8,  0.0,  0.4, 250./!guido
c------------------------------------------------------------------
	data icall /0/

c initialization

	if( icall .le. -1 ) return

	if( nlvdi .ne. nlvdim ) stop 'error stop biocos: nlvdim'

	if( icall .eq. 0 ) then
	  ibio = iround(getpar('ibio'))
	  if( ibio .le. 0 ) icall = -1
	  if( icall .le. -1 ) return

	  icall = 1
	  do k=1,nkn		!loop on nodes

	    lmax = ilhkv(k)
	    do l=1,lmax		!loop on levels
	      do i=1,ndim
		e(l,k,i) = einit(i)
	      end do
              t = tempv(l,k)
              s = saltv(l,k)
	      e(l,k,9) = oxysat(t,s)
	    end do
	  end do

          rkpar=getpar('chpar')
          call getaz(azpar)
          adpar=getpar('adpar')
          aapar=getpar('aapar')
          difmol=getpar('difmol')
          istot=iround(getpar('istot'))         !$$istot

	  write(6,*) 'biocos model initialized...'

	end if

c advection and diffusion

c	call tsmass(e(1,1,2),zeov,nlvdim,rmass0)
c	call tsmass(e(1,1,8),zeov,nlvdim,pmass0)

	do i=1,ndim

          do isact=1,istot                !$$istot
            call conz3d(e(1,1,i),saux1,saux2,saux3,saux4,dt,rkpar,difv
     +            ,difmol,azpar,adpar,aapar,istot,isact,nlvdim,nlv)
c	    call biobnd(e(1,1,i),ebound(i),rsv)
          end do

	end do

c	call tsmass(e(1,1,2),zenv,nlvdim,rmass1)
c	call tsmass(e(1,1,8),zenv,nlvdim,pmass1)
c
c	write(6,*) 'biocos: ',rmass0,rmass1,pmass0,pmass1

c normal call

c compute light intensity

	call setlux(it,light,e,hdknv,nlvdim,nkndim,ndim)

	dt = idt
	idtt = idt/ntot
	dtt = idtt

	do k=1,nkn		!loop on nodes

	  lmax = ilhkv(k)
	  do l=1,lmax		!loop on levels
	    t = tempv(l,k)
	    s = saltv(l,k)

	    lux = light(l,k)

	    do i=1,ndim
	      ee(i) = e(l,k,i)
	    end do

	    dtt = dt / ntot
	    do j=1,ntot
	      itt = it - idt + j*idtt
	      call reactor(itt,dtt,l,t,s,ee,lux)
	    end do

	    do i=1,ndim
	      e(l,k,i) = ee(i)
	    end do

	  end do
	end do

c check of nodal values

c	do k=165,171,2
c	  call writee(k,it,k,1,e,tempv,saltv,nlvdim,nkndim,ndim)
c	end do
c	k = 17
c	call writee(k,it,k,1,e,tempv,saltv,nlvdim,nkndim,ndim)

c write of results

	call biobnd0(ndim,e,ebound,dt)
	call bioout(ndim,e)

	end

c*************************************************************

	subroutine bioout(ndim,e)

c outputs results from bio model

	implicit none

        include 'param.h'

	integer ndim			!total number of state vectors
	real e(nlvdim,nkndim,ndim)	!state vector

        character*80 descrp
        common /descrp/ descrp
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it
        integer nlvdi,nlv
        common /level/ nlvdi,nlv
        integer ilhkv(1)
        common /ilhkv/ilhkv
        real hlv(1), hev(1)
        common /hlv/hlv, /hev/hev

        character*80 dir,nam,file
        integer ierr,ivar,i

        integer iround,ifileo
        real getpar

	integer iu
	save iu
        integer idtcon,itmcon,itcon
        save idtcon,itmcon,itcon

        integer icall
        save icall
        data icall /0/

        if(icall.eq.-1) return

c intialize parameters and arrays

        if(icall.eq.0) then
          idtcon=iround(getpar('idtcon'))
          itmcon=iround(getpar('itmcon'))

          if(itmcon.lt.itanf) itmcon=itanf
          if(idtcon.le.0) icall=-1
          if(itmcon+idtcon.gt.itend) icall=-1

          if(icall.eq.-1) return

          itcon=itmcon+idtcon

          call getfnm('datdir',dir)
          call getfnm('runnam',nam)
          call mkname(dir,nam,'.bio',file)

	  iu = ifileo(55,file,'unform','new')

          call wfnos(iu,3,nkn,nel,nlv,ndim,descrp,ierr)
          if(ierr.gt.0) goto 99
          call wsnos(iu,ilhkv,hlv,hev,ierr)
          if(ierr.gt.0) goto 99

          write(6,*) 'bioout : bio-file opened ',it
        end if

c normal call

        icall=icall+1

        if(it.le.itmcon) return

c here do writing

        if(it.lt.itcon) return

	do i=1,ndim
	  ivar = 100 + i
          call wrnos(iu,it,ivar,nlvdi,ilhkv,e(1,1,i),ierr)
          if(ierr.gt.0) goto 99
	end do

        write(6,*) 'bioout : bio-file written ',it

        itcon=itcon+idtcon

        return
99      continue
        write(6,*) 'error opening or writing file .bio : ',ierr
        stop 'error stop bioout'
        end

c*************************************************************

	subroutine writee(iunit,it,k,l,e,t,s,nlvdim,nkndim,ndim)

c formatted write for debug

	implicit none

	integer iunit
	integer it,k,l
	integer nlvdim,nkndim,ndim
	real e(nlvdim,nkndim,ndim)
	real t(nlvdim,nkndim)
	real s(nlvdim,nkndim)

	integer i

	write(iunit,'(i10,11f12.4)') it,
     +			(e(l,k,i),i=1,ndim),
     +			t(l,k),
     +			s(l,k)

	end

c*************************************************************

	subroutine biobnd0(ndim,e,ebound,dt)

c handles open boundary conditions
c
c since this is called AFTER the T/D step, the first time step
c is missed -> nothing is injected in the first time step
c there may be also some inconsistencies, since the water volume
c used is the one of the old time step, but the actual
c water volume injected is determined only afterwards

        implicit none

c parameter
        include 'param.h'
c arguments
	integer ndim
        real e(nlvdim,nkndim,ndim)
        real ebound(ndim)
	real dt
c local
        integer ibc,nbc,ibtyp,levmax,n,j,kn,i
	real rw,vol
c functions
	integer nbnds,itybnd,levbnd,kbnds,nkbnds
	real zvbnds

	nbc = nbnds()

        do ibc=1,nbc

	  ibtyp = itybnd(ibc)
          levmax = levbnd(ibc)

          rw = zvbnds(ibc)
	  vol = rw * dt

	  n = nkbnds(ibc)

          do j=1,n

             kn = kbnds(ibc,j)

	     do i=1,ndim
               call volno0(kn,levmax,nlvdim,e(1,1,i),vol,ebound(i))
	     end do
	  end do

	end do

	end

c*************************************************************

        subroutine biobnd(conz,rbound,rsv)

        implicit none

c parameter
        include 'param.h'
c arguments
        real conz(nlvdim,1)
        real rbound
	real rsv(1)
c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real eps1,eps2,pi,flag,high,hihi
        integer nlvdi,nlv
        integer ilhkv(1)
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /mkonst/ eps1,eps2,pi,flag,high,hihi
        common /level/ nlvdi,nlv
        common /ilhkv/ilhkv
c local
        integer k,l,lmax

        if(nlvdim.ne.nlvdi)stop'error stop : level dimension in biobnd'

        do k=1,nkn
          if( rsv(k) .ne. flag ) then
            lmax=ilhkv(k)
            do l=1,lmax
                conz(l,k) = rbound
            end do
          end if
        end do

        end

c***********************************************

       block data parameter

c  parametri usati nel programma fito2

	include 'reactor.h'

c  parametri per fito1

         DATA
     :   mumax1/0.12/,  !max crescita fito1
     :   tmax1/30/,   !tmax 25.6
     :   tott1/22.5/,   !tott 22.4
     :   coeff1/0.1157/,!coeff lassiter
     :   iopt1/50000/, !iopt steele
     :   ksp1/0.05/,  !semisat. P microM
     :   krf1/0.0025/, !resp.	!ggu changed from 0.007
     :   kmf1/0.0077/  !mort.

c  parametri per fito2

         data
     :   mumax2/0.18/,  !max crescita fito2
     :   tmax2/30/,   !tmax 
     :   tott2/10/,   !tott
     :   coeff2/0.1157/,!coeff lassiter
     :   ksdoc2/3.125/,    !semisat.doc 
     :   ksp2/0.05/,    !semisat. P microM
     :   kmf2/0.01/    !mort.

c   parametri per zoo

         data
     :   kgr1/0.05/,     !max grazing 
     :   kfz1/250/,       !semisat. graz. microM
     :   alfa/3./,       !preferenza graz. (0 solo fito, 5 piu' bac che fito)
     :   eff1/0.7/,      !efficienza conversione
     :   eff2/0.7/,      !efficienza conversione
     :   kmz/0.005/,     !mort.
     :   kexz/0.002/     !escrez.

c   altri parametri

         data
     :   kdcd/0.00477/, !decadimento pom 
     :   kdcp/0.00477/, !decadimento pop 
     :   ksodec/62.5/,    !semisat per OD consumo microM
     :   used/0.00/,   !sedimentazione detrito 0.004
     :   krear/0.04584/, !riareazione
     :   kestf/0.0048/,     !estinzione fito microM 
     :   kestd/0.0048/,     !estinzione detrito microM
     :   kestw/0.1/,     !estinzione acqua
     :   sal/30/,       !salinita'
     :   rpc1/0.00862/,   !rapp P/C in fito1 (uMP/uMC; 1/116)
     :   rpc2/0.0309/,   !rapp P/C in bacteria (uMP/uMC; 1/116)
     :   rpcz/0.00862/,   !rapp P/C in zoo (uMP/uMC; 1/116)
     :   roc1/1/,    !ossigeno consumato da fito1
     :   roc2/1/,    !ossigeno consumato da fito2
     :   rocd/1/     !ossigeno consumato da detrito
         end

         block data variab

c valori del nome delle variabili 
c nome delle variabili e loro indice per il vettore x, xold, dx
c vanno bene anche per e, ma non per efix

	include 'reactor.h'

        data
     :  po4  /1/,
     :  phy1 /2/,
     :  bac /3/,
     :  zoo  /4/,
     :  detc /5/,
     :  detp /6/,
     :  doc /7/,
     :  dop /8/,
     :  oxy  /9/

        end

c***********************************************

