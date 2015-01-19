c
c nosmaniq.f quality map, bod5 map
c donata 4 luglio 2001
c !!!ATTENZIONE!!! decidere quale indice usare: Vismara o TRIX e settare 
c in subroutine   mapq
c set ftrix
c
c****************************************************************

	program nosextr


	include 'param.h'


c--------------------------------------------------
	include 'basin.h'


c--------------------------------------------------

	character*80 title
c	real cv(nkndim)
c	real cv3(nlvdim,nkndim)
	real bod5v(nlvdim,nkndim)
	real quality(nlvdim,nkndim)
	real nh3v(nlvdim,nkndim)
	real noxv(nlvdim,nkndim)
	real opo4v(nlvdim,nkndim)
	real phytov(nlvdim,nkndim)
	real cbodv(nlvdim,nkndim)
	real doxv(nlvdim,nkndim)
	real onv(nlvdim,nkndim)
	real opv(nlvdim,nkndim)
	real zoov(nlvdim,nkndim)
	real salv(nlvdim,nkndim)
	real salvmin(nlvdim,nkndim)
	real salvmax(nlvdim,nkndim)
	real tempv(nlvdim,nkndim)
	real tempvmin(nlvdim,nkndim)
	real tempvmax(nlvdim,nkndim)

	integer ilhkv(nkndim)
	real hlv(nlvdim)
	real hev(neldim)

        real nh3,nox,opo4,phyto,cbod,dox,on,op,zoo
	real sal,temp
	
	logical berror,bbod5
	real par1,par2

c--------------------------------------------------
c--------------------------------------------------
c--------------------------------------------------

c	if bbod5=true   -> bod5 map
c	if bbod5=false  -> quality map

	bbod5 = .false.

	par1 = 0.2
	par2 = 0.1

	nread=0
	rnull=0.

        if(iapini(3,nkndim,neldim,0).eq.0) then
       	stop 'error stop : iapini'
        end if
        
        
        nin=ideffi('datdir','runnam','.bio','unform','old')
        if(nin.le.0) goto 100

        ntin=ideffi('datdir','runnam','.tav','unform','old')
        if(nin.le.0) goto 100

        nsin=ideffi('datdir','runnam','.sav','unform','old')
        if(nin.le.0) goto 100

        nvers=3
	call rfnos(nin,nvers,nkn,nel,nlv,nvar,title,ierr)
        if(ierr.ne.0) goto 100

        write(6,*) 'nvers    : ',nvers
        write(6,*) 'nkn,nel  : ',nkn,nel
        write(6,*) 'nlv,nvar : ',nlv,nvar
        write(6,*) 'title    : ',title
	write(6,*) 'kbod ,knit ',par1,par2

        call dimnos(nin,nkndim,neldim,nlvdim)

	call rsnos(nin,ilhkv,hlv,hev,ierr)
        if(ierr.ne.0) goto 100

	write(6,*) 'Available levels: ',nlv
	write(6,*) (hlv(l),l=1,nlv)

c--------
        nvers=3
        call rfnos(ntin,nvers,nkn,nel,nlv,nvar,title,ierr)
        if(ierr.ne.0) goto 100

        write(6,*) 'nvers    : ',nvers
        write(6,*) 'nkn,nel  : ',nkn,nel
        write(6,*) 'nlv,nvar : ',nlv,nvar
        write(6,*) 'title    : ',title
        write(6,*) 'kbod ,knit ',par1,par2

        call dimnos(ntin,nkndim,neldim,nlvdim)

        call rsnos(ntin,ilhkv,hlv,hev,ierr)
        if(ierr.ne.0) goto 100

        write(6,*) 'Available levels: ',nlv
        write(6,*) (hlv(l),l=1,nlv)
c--------
c--------
        nvers=3
        call rfnos(nsin,nvers,nkn,nel,nlv,nvar,title,ierr)
        if(ierr.ne.0) goto 100

        write(6,*) 'nvers    : ',nvers
        write(6,*) 'nkn,nel  : ',nkn,nel
        write(6,*) 'nlv,nvar : ',nlv,nvar
        write(6,*) 'title    : ',title
        write(6,*) 'kbod ,knit ',par1,par2

        call dimnos(nsin,nkndim,neldim,nlvdim)

        call rsnos(nsin,ilhkv,hlv,hev,ierr)
        if(ierr.ne.0) goto 100

        write(6,*) 'Available levels: ',nlv
        write(6,*) (hlv(l),l=1,nlv)
c--------

	high = 1.e+30

	nout=ifileo(80,'quality.nos','unformatted','new')
	if(nout.le.0) goto 100

        nvers=3
	it = 0
	nvar = 1
	call wfnos(nout,nvers,nkn,nel,nlv,nvar,title,ierr)
	call wsnos(nout,ilhkv,hlv,hev,ierr)


  300   continue
c	prova

        ivar=161
        call rdnos(ntin,it,ivar,nlvdim,ilhkv,tempv,ierr)
c        write(6,*) ': tempv',ivar
        if(ivar.ne.161) write(6,*) 'error in reading var temp : ',ivar
c        ivar=161
        call rdnos(ntin,it,ivar,nlvdim,ilhkv,tempvmin,ierr)
c        write(6,*) ': tempv',ivar
        if(ivar.ne.162) write(6,*) 'error in reading var temp : ',ivar
c        ivar=161
        call rdnos(ntin,it,ivar,nlvdim,ilhkv,tempvmax,ierr)
c        write(6,*) ': tempv',ivar
        if(ivar.ne.163) write(6,*) 'error in reading var temp : ',ivar
        
c        ivar=171
        call rdnos(nsin,it,ivar,nlvdim,ilhkv,salv,ierr)
c        write(6,*) ': salv',ivar
        if(ivar.ne.171) write(6,*) 'error in reading var sal : ',ivar
c        ivar=171
        call rdnos(nsin,it,ivar,nlvdim,ilhkv,salvmin,ierr)
c        write(6,*) ': salv',ivar
        if(ivar.ne.172) write(6,*) 'error in reading var sal : ',ivar
c        ivar=171
        call rdnos(nsin,it,ivar,nlvdim,ilhkv,salvmax,ierr)
c        write(6,*) ': salv',ivar
        if(ivar.ne.173) write(6,*) 'error in reading var sal : ',ivar
 

	ivar=71
	call rdnos(nin,it,ivar,nlvdim,ilhkv,nh3v,ierr)
        if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
        if(ivar.ne.71) write(6,*) 'error in reading var 71 : ',ivar
        if(ierr.ne.0) goto 100
c	write(6,*) 'nh3v: ',ivar

	ivar=72
	call rdnos(nin,it,ivar,nlvdim,ilhkv,noxv,ierr)
c	write(6,*) ': noxv',ivar
        if(ivar.ne.72) write(6,*) 'error in reading var 72 : ',ivar

	ivar=73
	call rdnos(nin,it,ivar,nlvdim,ilhkv,opo4v,ierr)
c	write(6,*) 'opo4v: ',ivar
        if(ivar.ne.73) write(6,*) 'error in reading var 73 : ',ivar

	ivar=74
	call rdnos(nin,it,ivar,nlvdim,ilhkv,phytov,ierr)
c	write(6,*) ': phytov',ivar
        if(ivar.ne.74) write(6,*) 'error in reading var 74 : ',ivar

	ivar=75
	call rdnos(nin,it,ivar,nlvdim,ilhkv,cbodv,ierr)
c	write(6,*) ': cbodv',ivar
        if(ivar.ne.75) write(6,*) 'error in reading var 75 : ',ivar

	ivar=76
	call rdnos(nin,it,ivar,nlvdim,ilhkv,doxv,ierr)
c	write(6,*) ': doxv',ivar
        if(ivar.ne.76) write(6,*) 'error in reading var 76 : ',ivar

	ivar=77
	call rdnos(nin,it,ivar,nlvdim,ilhkv,onv,ierr)
c	write(6,*) ': onv',ivar
        if(ivar.ne.77) write(6,*) 'error in reading var 77 : ',ivar

	ivar=78
	call rdnos(nin,it,ivar,nlvdim,ilhkv,opv,ierr)
c	write(6,*) ': opv',ivar
        if(ivar.ne.78) write(6,*) 'error in reading var 78 : ',ivar

	ivar=79
	call rdnos(nin,it,ivar,nlvdim,ilhkv,zoov,ierr)
c	write(6,*) ': zoov',ivar
        if(ivar.ne.79) write(6,*) 'error in reading var 79 : ',ivar


        nread=nread+1
c       write(6,*) 'time : ',it,ivar

	do k=1,nkn
	  maxlev = ilhkv(k)
	  do l=1,maxlev
	    iaux = iaux + 1
	    nh3=nh3v(l,k)
	    nox=noxv(l,k)
	    opo4=opo4v(l,k)
            phyto=phytov(l,k)
	    cbod =cbodv(l,k) 
	    dox =doxv(l,k) 
	    on =onv(l,k) 
	    op =opv(l,k) 
	    zoo=zoov(l,k)            
	    sal=salv(l,k)
	    temp=tempv(l,k)
c
c	bod5 computation:
c
	         aux1 = 1. - exp( -5. * par1 )
	         aux2 = 1. - exp( -5. * par2 )
	         bod5 = cbod*aux1 + (64./14.) * nh3 * aux2
	         bod5v(l,k) = bod5
c	
c	quality map computation
c
c	    call mapq(dox,cbod,nh3,bod5,quality(l,k))
	call mapq(nh3,nox,opo4,phyto,cbod,dox,on,op,
     $ zoo,sal,temp,bod5,quality(l,k))
	

	  end do
	end do
c
c	write bod5 or quality map
c
	if( bbod5 ) then
	  ivar =80 
	  write(6,*) 'bod5v: ',ivar
	  call wrnos(nout,it,ivar,nlvdim,ilhkv,bod5v,ierr)
	else
	  ivar =81 
c	  write(6,*) 'quality: ',ivar
	  call wrnos(nout,it,ivar,nlvdim,ilhkv,quality,ierr)
	end if

	goto 300

  100	continue


	write(6,*)
	write(6,*) nread,' records read'
	write(6,*)


	write(6,*)
	write(6,*) 'data written to file quality.nos'
	write(6,*)

	end

c***************************************************************

	subroutine mapq(nh3,nox,opo4,phyto,cbod,dox,on,op,zoo,
     $ sal,temp,bod5,quality)


c creates index of water quality

	implicit none

	integer ftrix
	real nh3,nox,opo4,phyto,cbod,dox,on,op,zoo,bod5,quality
	real sal, temp
	real o2sat,o2satp
	real aux1, aux2, aux3
	real trix
	real auxt1,auxt2,auxt3
	real cl,tk,rlncs

c 		aux1 è la variabile punteggio O2 indice vismara
c		aux2 è la variabile punteggio BOD indice vismara
c		aux3 è la variabile punteggio N-NH4 indice vismara



	ftrix=1 	!if ftrix=0 vismara map
			!if ftrix=1 vollenweider, trix map

          CL = SAL/1.80655
          TK = temp + 273.
          RLNCS = - 139.34411 + (1.575701E05/TK) - (6.642308E07/TK**2) +
     1   (1.243800E10/TK**3) - (8.621949E11/TK**4) -
     2   (CL*(3.1929E-02 - (19.428/TK) + (3.8673E03/TK**2)))

            o2sat= EXP (RLNCS)
	    
	if(dox.le.0)then
	dox=0.000005
	end if

		o2satp=(dox*100.)/o2sat
	
	if (ftrix.eq.0.)then

       	    if (o2satp .le. 110. .and. o2satp .gt. 90.) then

	 	aux1=1. 
 	         else if (o2satp .le. 90. .and. o2satp .gt. 70.) then
	        aux1=2. 
       		 else if (o2satp .le.120. .and. o2satp .gt. 110.) then
		aux1=2.           
        	else if (o2satp .le.70. .and. o2satp .gt.50.) then
		 aux1=3. 
        	else if (o2satp .le.130. .and. o2satp .gt.120.) then
	 	 aux1=3.        
        	else if (o2satp .le.50. .and. o2satp .gt. 30.) then
		 aux1=4.        
	        else if (o2satp .le.30. .or. o2satp .ge.130.) then
		 aux1=5.        
	   end if

          if (cbod   .le. 3.) then
		aux2=1.
	  else if (bod5   .le. 6.) then
		aux2=2.
	  else if (bod5   .le. 9.) then
		aux2=3.
          else if (bod5   .le. 15.) then
		aux2=4.
          else if (bod5   .gt. 15.) then
		aux2=5.
	  end if

	     if (nh3 .le. 0.4) then
		 aux3=1.
        	else if (nh3   .le. 1.)	then
		 aux3=2.
	        else if (nh3   .le. 2.)	then
		 aux3=3.
                else if (nh3   .le. 5.) then
		 aux3=4.
                else if (nh3   .gt. 5.)	then
		 aux3=5.
	     end if
c
c	quality vismara
c
	   quality=aux1+aux2+aux3
c
	else if (ftrix.eq.1)then
c
c	l'indice si ricava usando le concentrazioni espresse in mg/m3
c	il firoplancton si converte da carbonio a  clorofilla
c
	auxt1=(phyto/30.)*1000.
	auxt2=(nh3+nox)*1000.
	auxt3=(opo4+op)*1000.
	
	trix=(log10(auxt1*o2satp*auxt2*auxt3+1.5))/1.2
c		trix=(log(((phyto/30.)/1000)*o2satp*((nh3+nox)/1000)*
c     $  ((opo4+op)/1000)+(1.5)))/1.2 	

	quality=trix
	
	 if(quality.lt.0) then
        stop 'error stop : quality'
        end if

	
	end if
c	write(65,*) nh3,nox,opo4,phyto,cbod,dox,on,op,zoo 
	write(65,*)quality
	end

c***************************************************************

