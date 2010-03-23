c
c $Id: sublgr.f,v 1.11 2003/07/31 11:14:10 georg Exp $
c
c lagrangian model routines
c
c contents :
c
c subroutine lgrini(nfldim,nfddim)              float initialization
c subroutine lgrfln                             float tracking 
c subroutine lgrfdn(nfddim)                     float diffusion
c function elnew(ie,xp,yp,xpo,ypo)              finds new element (shell)
c function ielemf(ie,xin,yin)                   finds element, fk method
c function ielema(ie,xin,yin)                   finds element, old method, all
c function ielemo(ie,xin,yin)                   finds element, fk method, all
c subroutine transn(ie,xinou,yinou,tanf,dl)     advection of particle
c subroutine diffus(ie,xinou,yinou,tanf)        diffusion of particle
c subroutine relran(nfddim)                     release random particles
c real function port(k1,k2,h1,h2)               compute discharge through line
c subroutine ranxy(rlenal,kranf,krend,irand,x,y,iee)
c                               gets x,y for a random particle at boundary
c real function pllvel(ie,u,v)
c                               get velocity for element and make parallel
c subroutine intval(ie,xp,yp,u,v,uh,vh)
c                               interpolates u/v onto a point xp,yp in el. ie
c subroutine lgrerr(mode,ip,ie,xp,yp)
c                               error output routine
c real function plline(ie1,ie2,xp,yp,dt)
c                               advect on line
c subroutine relra1(nfddim)
c                               release random particles at one node
c subroutine relra2(nfddim)
c                               release random particles at one node
c
c revision log :
c
c 25.09.1991	ggu	$$new25.09.91 - do not stop advection -> go on
c 19.11.1991	ggu	$$new19.11.91 - take dry areas into account
c 13.12.1991	ggu	$$new13.12.91 - may be bug fix ?
c 25.03.1998	ggu	$$new25.03.98 - initialize iefdv
c 25.03.1998	ggu	integrated changes from technital
c
c notes :
c
c ner is not defined
c
c eps   epsilon for fk method if particle is in element
c eps1  epsilon to compute portion of path gone
c eps2  epsilon to put particle into domain from open boundary
c
c common block /ffloat/ must be saved somewhere		!FIXME
c
c************************************************************
c
	subroutine lgrini(nfldim,nfddim)
c
c float routine initialization
c
c dlfl          0  : no float tracking
c               >0 : float tracking on
c dlfd          0  : no float diffusion
c               >0 : float diffusion on
c
	implicit none
c
c arguments
	integer nfldim,nfddim
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer itanf,itend,idt,nits,niter,it
	integer nfloat,nfdoat
	real tact,dthy,dlfl,dlfd,ddif
	integer itaflt,itafdf
	integer ieflv(1),iefdv(1)
	real xpflv(1),ypflv(1),xpfdv(1),ypfdv(1)
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /femtim/ itanf,itend,idt,nits,niter,it
	common /ffloat/ nfloat,nfdoat
	common /tflfd/ tact,dthy,dlfl,dlfd,itaflt,itafdf,ddif
	common /ieflv/ieflv,/xpflv/xpflv,/ypflv/ypflv
	common /iefdv/iefdv,/xpfdv/xpfdv,/ypfdv/ypfdv
c local
	character*80 name
	integer i,nb,dddif,ier,idum,nverfd
c function
	real getpar
	integer ielemo,ifileo,ckfdf,rffdf,rdfdf
c
	dlfl=getpar('dlfl')
	dlfd=getpar('dlfd')
	itaflt=getpar('itaflt')
	itafdf=getpar('itafdf')
	ddif=getpar('ausfd')
	dddif=getpar('aushfd')
c
	dthy=idt
c
c test variables
c
	call iarini
	call farini
	call farset(3,float(idt))
c
c float tracking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	if(dlfl.gt.0) then
c
c read in supplemental initial distribution
c
	   call getfnm('flinit',name)
	   if(name.ne.' ') then
		nb=ifileo(55,name,'form','old')
		if(nb.le.0) then
			write(6,*) 'Cannot open file flinit :'
			write(6,*) name
			stop 'error stop : lgrini'
		end if
		read(nb,*,end=33) 
     +                  (xpflv(i),ypflv(i),i=nfloat+1,nfldim)
   33           close(nb)
		write(6,*) i-nfloat-1,' supplemental particles read for fl'
		write(6,*)
		nfloat=i-1
	   end if
c
	   write(6,*) nfloat,' total number of particles read for fl'
	   write(6,*)
c
c find element numbers
c
		do i=1,nfloat
		  ieflv(i)=ielemo(0,xpflv(i),ypflv(i))
		end do
c
c               write(6,3040) (ieflv(i),xpflv(i),ypflv(i),i=1,nfloat)
c 3040          format(i10,2f12.2)
c
c particle on land
c
		do i=1,nfloat
		  if(ieflv(i).eq.0) then
			call lgrerr(1,i,0,xpflv(i),ypflv(i))
			ieflv(i)=0
		  end if
		end do
	end if
c
c float diffusion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	if(dlfd.gt.0) then
c
c read in initial distribution
c
	   call getfnm('fdinit',name)
	   if(name.ne.' ') then
		nb=ifileo(55,name,'unform','old')
		if(nb.le.0) then
			write(6,*) 'Cannot open file fdinit'
			write(6,*) name
			stop 'error stop : lgrini'
		end if
                if(ckfdf(nb).gt.0) then
                  nverfd=0
                  nfdoat=nfddim
                  ier=rffdf(nb,6,nverfd,idum,name)
                  ier=rdfdf(nb,6,idum,nfdoat,xpfdv,ypfdv,iefdv)
                  if(ier.ne.0) then
                        write(6,*) 'Error reading fdinit'
                        stop 'error stop : lgrini'
                  else
                        write(6,*) 'title of fdf init file read :'
                        write(6,*) name
                  end if
                else
                  rewind(nb)
                  read(nb) nfdoat
                  if(nfdoat.gt.nfddim) then
                        write(6,*) 'lgrini : too much data'
                        write(6,*) 'Can read only ',nfddim,' particles'
                        nfdoat=nfddim
                  end if
                  read(nb) (xpfdv(i),ypfdv(i),i=1,nfdoat)
		  do i=1,nfdoat		!$$new25.03.98
		    iefdv(i) = 0
		  end do
                end if
                close(nb)
	   end if
c
	   write(6,*) nfdoat,' particles read for fd'
	   write(6,*)
c
c austausch coefficient
c
	   if(dddif.eq.0.) then
c
c               old formula, do not use anymore
c               U' = ausfd/24 ==> AH = (ausfd**2 * dt) /(6 * 24**2)
c               (remember that ddif = 2*U')
c
		ddif=ddif/12.
	   else
c
c               new formula, time step dependent
c               U' = sqrt(6*aushfd/dt) ==> AH = aushfd
c               (remember that ddif = 2*U')
c
		ddif=2.*sqrt(6.*dddif/float(idt))
	   end if
c
	   write(6,*)
	   write(6,*) 'Band width for random velocity : ',0.5*ddif
	   write(6,*)
c
c find element numbers
c
		i=1
		do while(i.le.nfdoat)
		   iefdv(i)=ielemo(iefdv(i),xpfdv(i),ypfdv(i))
		   if(iefdv(i).gt.0) then
			i=i+1
		   else
			xpfdv(i)=xpfdv(nfdoat)
			ypfdv(i)=ypfdv(nfdoat)
			nfdoat=nfdoat-1
		   end if
		end do
	end if
c
	return
	end
c
c***************************************************************
c
	subroutine lgrfln
c
c float tracking program
c
c ie            >0  in system
c               =0  lost at first step
c               =-1 lost through open boundary
c
	implicit none
c
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer itanf,itend,idt,nits,niter,it
	integer nfloat,nfdoat
	real tact,dthy,dlfl,dlfd,ddif
	integer itaflt,itafdf
	integer ieflv(1)
	real xpflv(1),ypflv(1)
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /femtim/ itanf,itend,idt,nits,niter,it
	common /ffloat/ nfloat,nfdoat
	common /tflfd/ tact,dthy,dlfl,dlfd,itaflt,itafdf,ddif
	common /ieflv/ieflv,/xpflv/xpflv,/ypflv/ypflv
c local
        logical berror
        integer ifl,ie,i
        real dl,xn,yn,tanf
c functions
        integer ieext,iround
        integer iarmax,iarget
        real farmax,farget,getpar
c
	if(dlfl.le.0) return
	if(it.lt.itaflt) return
c
	call delta(-4)
c
c test variables
c
c       call iarini
c       call farini
c       call farset(3,float(idt))
c
        berror=iround(getpar('lgrerr')).ge.2    !write err code every time step
c
c cccccccccccccc
c
	dl=dlfl
	tact=it
	tanf=it-idt
c
	do ifl=1,nfloat
	   ie=ieflv(ifl)
	   xn=xpflv(ifl)
	   yn=ypflv(ifl)
c
	   if(ie.gt.0) then             !particle in system
	     call transn(ie,xn,yn,tanf,dl)
	     if(ie.lt.0) then
		call lgrerr(2,ifl,ie,xn,yn)
	     else if(ie.eq.0) then
		call lgrerr(3,ifl,ie,xn,yn)
	     end if
	   end if
c
	   ieflv(ifl)=ie
	   xpflv(ifl)=xn
	   ypflv(ifl)=yn
	end do
c
        if(berror) call lgrerr(0,0,0,0.,0.)
c
	write(ner,*) 'it : ',it,'   particles : ',nfloat
c
	if(niter.eq.nits) then
	  write(ner,*)
	  write(ner,'(1x,a,7i10)') 'ntst :',(iarget(i),i=1,7)
	  if(farget(1).gt.0.) call farset(2,farget(2)/farget(1))
	  write(ner,'(1x,a,e14.6,2f10.3)') 'dt   :',(farget(i),i=1,3)
	  if(.not.berror) call lgrerr(0,0,0,0.,0.)
	  write(ner,*)
	  call iarini
	  call farini
	  call farset(3,float(idt))
	end if
c
	call delta(4)
c
	return
	end
c
c*******************************************************
c
	subroutine lgrfdn(nfddim)
c
c float diffusion program
c
	implicit none
c
c arguments
	integer nfddim
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer itanf,itend,idt,nits,niter,it
	integer nfloat,nfdoat
	real tact,dthy,dlfl,dlfd,ddif
	integer itaflt,itafdf
	integer iefdv(1)
	real xpfdv(1),ypfdv(1)
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /femtim/ itanf,itend,idt,nits,niter,it
	common /ffloat/ nfloat,nfdoat
	common /tflfd/ tact,dthy,dlfl,dlfd,itaflt,itafdf,ddif
	common /iefdv/iefdv,/xpfdv/xpfdv,/ypfdv/ypfdv
c local
        logical berror
        integer i,ifdac,ifd,ie
        real tanf,dl,xn,yn
c function
        integer iround
        integer iarmax,iarget
        real farmax,farget,getpar
c
	if(dlfd.le.0) return
	if(it.lt.itafdf) return
c
	call delta(-5)
c
c test variables
c
c       call iarini
c       call farini
c       call farset(3,float(idt))
c
	berror=iround(getpar('lgrerr')).ge.2    !write err code every time step
c
c cccccccccccccc
c
	dl=dlfd
	tact=it
	tanf=it-idt
	ifdac=0
c
	do ifd=1,nfdoat
	   ie=iefdv(ifd)
	   xn=xpfdv(ifd)
	   yn=ypfdv(ifd)
c
	   call transn(ie,xn,yn,tanf,dl)
	   call diffus(ie,xn,yn,tanf)
c
	   if(ie.gt.0) then
		ifdac=ifdac+1
		iefdv(ifdac)=ie
		xpfdv(ifdac)=xn
		ypfdv(ifdac)=yn
	   end if
	end do
c
	nfdoat=ifdac
c
	call relran(nfddim)
	call relra1(nfddim)
	call relra2(nfddim)
c
        if(berror) call lgrerr(0,0,0,0.,0.)
c
	write(ner,*) 'it : ',it,'   particles : ',nfdoat
c
	if(niter.eq.nits) then
	  write(ner,*)
	  write(ner,'(1x,a,7i10)') 'ntst :',(iarget(i),i=1,7)
	  if(farget(1).gt.0.) call farset(2,farget(2)/farget(1))
	  write(ner,'(1x,a,e14.6,2f10.3)') 'dt   :',(farget(i),i=1,3)
	  if(.not.berror) call lgrerr(0,0,0,0.,0.)
	  write(ner,*)
	  call iarini
	  call farini
	  call farset(3,float(idt))
	end if
c
	if(nfdoat.gt.nfddim) then
		write(6,*) 'Too many particles in fd'
		write(6,*) 'nfddim,nfdoat :',nfddim,nfdoat
		write(6,*) 'Increase dimension of nfddim'
		stop 'error stop : lgrfdn'
	end if
c
	call delta(5)
c
	return
	end
c
c***********************************************************
c
	function elnew(ie,xp,yp,xpo,ypo)
c
c find new element
c
c uses fk method
c follows the path of the particle from one element to its
c neighbour stopping only if the system boundary is reached
c if the particle has not gone out through the open boundary
c the coordinates given back should be always in the element ie
c if no element cant be found at some stage for the computed
c coordinates the original coordinates and element are given back
c
c ie            element number of particle (entry and return)
c xp,yp         new coordinates of particle
c xpo,ypo       old coordinates of particle (are in element ie for sure)
c elnew         [0...1] : ok   -1 : error
c
	implicit none
c
	real elnew
c arguments
	integer ie
	real xp,yp,xpo,ypo
c parameter
	integer itbas
	parameter (itbas=20)
	real eps,eps1
	parameter (eps1=1.e-7,eps=5.e-3)
c common
	real ev(13,1)
	integer ieltv(3,1)
	common /ev/ev
	common /ieltv/ieltv
	integer iwegv(1)
	common /iwegv/iwegv
c local
	integer iee,ifk,ifj,k,k1,k2,ien,i
	real fkf(3)
	real fk,fkold,frel,frel1,frel2
c functions
	integer ieext,ielemf
c statement function
	real a,b,c
	a(k,ie)=ev(k,ie)
	b(k,ie)=ev(k+3,ie)
	c(k,ie)=ev(k+6,ie)
c
	iee=ie
c
			call iarinc(5)
			call delta(-(itbas+5))
	ifk=0
	ifj=0
	do k=1,3
	   fk=a(k,iee)+xp*b(k,iee)+yp*c(k,iee)
	   if(fk.lt.-eps) then
		ifj=ifj+k
		ifk=ifk+1
		fkf(k)=-fk
	   end if
	end do
c
	if(ifk.eq.0) then
c               ok
		elnew=1.
			call delta(itbas+5)
		return
	else if(ifk.eq.1) then
		fkold=a(ifj,iee)+xpo*b(ifj,iee)+ypo*c(ifj,iee)
		if(fkold.lt.eps1) then
			frel=0.
		else
			frel=fkold/(fkf(ifj)+fkold)
		end if
	else if(ifk.eq.2) then
		k=6-ifj
		k1=mod(k,3)+1
		fkold=a(k1,iee)+xpo*b(k1,iee)+ypo*c(k1,iee)
		if(fkold.lt.eps1) then
			frel1=0.
		else
			frel1=fkold/(fkf(k1)+fkold)
		end if
		k2=mod(k1,3)+1
		fkold=a(k2,iee)+xpo*b(k2,iee)+ypo*c(k2,iee)
		if(fkold.lt.eps1) then
			frel2=0.
		else
			frel2=fkold/(fkf(k2)+fkold)
		end if
		if(frel1.lt.frel2) then
			ifj=k1
			frel=frel1
		else
			ifj=k2
			frel=frel2
		end if
	else
                write(6,*) 'value not possible for ifk :',ifk
                write(6,*) iee,ieext(iee),ifj,eps
                write(6,*) xp,yp
                write(6,*) (ieltv(i,iee),i=1,3)
                write(6,*) (a(i,iee),i=1,3)
                write(6,*) (b(i,iee),i=1,3)
                write(6,*) (c(i,iee),i=1,3)
                write(6,*) '  ieltv'
                write(6,*) (ev(i,iee),i=1,13)
                write(6,*) '  gggg'
                write(6,*) (fkf(i),i=1,3)
                write(6,*) '  gggg'
                stop 'error stop : elnew'
	end if
c
	ien=ieltv(ifj,iee)
c
c	if(ien.gt.0) then
        if(ien.gt.0.and.iwegv(ien).eq.0) then   !????? $$new19.11.91
c          o.k.
	   iee=ien
	   elnew=frel
	else if(ien.lt.0) then
c          open boundary
	   ie=ien
	   elnew=-1.
		call delta(itbas+5)
	   return
	else
c          out of system, take old element
c          ...reduce also distance so that part. stops befor boundary
	   frel=frel*0.9
	   elnew=frel
	end if
c
c control if really in element
c
	xp=xpo+(xp-xpo)*frel
	yp=ypo+(yp-ypo)*frel
	iee=ielemf(iee,xp,yp)
c
	if(iee.eq.0) then
		call lgrerr(4,0,iee,xp,yp)
		call lgrerr(8,0,ie,xpo,ypo)
		xp=xpo
		yp=ypo
		iee=ie
		elnew=-1.
	end if
c
	ie=iee
			call delta(itbas+5)
	return
	end
c
c*******************************************************
c
	function ielemf(ie,xin,yin)
c
c finds element to coordinates of particle
c
c uses fk method
c starts from given element and looks to its neighbours
c if it cant find element in this way it starts an
c overall search calling ielemo
c
c ie            >0 : element number that has to be checked first
c               ...if particle is not in element ie ,all other
c               ...elements are checked
c               =0 : check all elements
c xin,yin       coordinates of particle
c ielemf        element number the particle is lying in
c               ...0 if point is in no element
c
	implicit none
c
	integer ielemf
c arguments
	integer ie
	real xin,yin
c parameter
	integer itbas
	parameter (itbas=20)
	real eps
	parameter (eps=5.e-3)
c common
	real ev(13,1)
	integer ieltv(3,1)
	common /ev/ev,/ieltv/ieltv
c local
	integer iee,k
c function
	integer ielemo
c
			call iarinc(6)
			call delta(-(itbas+6))
	iee=ie
c
	do while(iee.gt.0)
c
	do k=1,3
	   if(ev(k,iee)+xin*ev(k+3,iee)+yin*ev(k+6,iee).lt.-eps) goto 10
	end do
c
c       in element
	ielemf=iee
			call delta(itbas+6)
	return
c
c       out of element
   10   iee=ieltv(k,iee)
c
	end do
c
			call iarinc(7)
			call delta(-(itbas+7))
	ielemf=ielemo(ie,xin,yin)
			call delta(itbas+7)
c
			call delta(itbas+6)
	return
	end
c
c*******************************************************
c
	function ielema(ie,xin,yin)
c
c finds element to coordinates of particle
c
c scalar product method (old method)
c looks in all elements starting with given elements
c and proceeding with closest element number
c (e.g. : start=47 ==> 47,46,48,45,49...)
c
c not used currently
c
c ie            >0 : element number that has to be checked first
c               ...if particle is not in element ie ,all other
c               ...elements are checked
c               =0 : check all elements
c xin,yin       coordinates of particle
c ielema        element number the particle is lying in
c               ...0 if point is in no element
c
	implicit none
c
	integer ielema
c arguments
	integer ie
	real xin,yin
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real xgv(1),ygv(1)
	integer nen3v(3,1)
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /xgv/xgv, /ygv/ygv
	common /nen3v/nen3v
c local
	integer idif,ifac,ieo,iee,k1,k2,kn1,kn2
	real xp,yp
	real x1,x2,y1,y2,scal
c
	idif=0
	ifac=1
	ieo=0
	iee=ie
	if(iee.le.0.or.iee.gt.nel) iee=1+nel/2
	xp=xin
	yp=yin
c
	do while(iee.ne.ieo)
c
	do k1=1,3
	   k2=mod(k1,3)+1
	   kn1=nen3v(k1,iee)
	   kn2=nen3v(k2,iee)
	   x1=xgv(kn1)
	   x2=xgv(kn2)
	   y1=ygv(kn1)
	   y2=ygv(kn2)
	   scal=(y1-y2)*(xp-x1)+(x2-x1)*(yp-y1)
	   if(scal.lt.0.) goto 10               !wrong side
	end do
c
	ielema=iee      !element found
	return
c
   10   continue
c
	ieo=iee
	idif=idif+1
	ifac=-ifac
	iee=iee+ifac*idif
	if(iee.le.0) then
		iee=iee+nel
	else if(iee.gt.nel) then
		iee=iee-nel
	end if
c
	end do
c
	ielema=0        !no element found
c
	return
	end
c
c*******************************************************
c
	function ielemo(ie,xin,yin)
c
c finds element to coordinates of particle
c
c uses fk method
c looks in all elements starting with given elements
c and proceeding with closest element number
c (e.g. : start=47 ==> 47,46,48,45,49...)
c
c ie            >0 : element number that has to be checked first
c               ...if particle is not in element ie ,all other
c               ...elements are checked
c               =0 : check all elements
c xin,yin       coordinates of particle
c ielemo        element number the particle is lying in
c               ...0 if point is in no element
c
	implicit none
c
	integer ielemo
c arguments
	integer ie
	real xin,yin
c parameter
	real eps
	parameter (eps=5.e-3)
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real ev(13,1)
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /ev/ev
c local
	integer idif,ifac,ieo,iee,k
	real xp,yp
c
	idif=0
	ifac=1
	ieo=0
	iee=ie
	if(iee.le.0.or.iee.gt.nel) iee=1+nel/2
	xp=xin
	yp=yin
c
	do while(iee.ne.ieo)
c
	do k=1,3
	   if(ev(k,iee)+xp*ev(k+3,iee)+yp*ev(k+6,iee).lt.-eps) goto 10
	end do
c
c       in element
	ielemo=iee
	return
c
c       out of element
   10   continue
c
	ieo=iee
	idif=idif+1
	ifac=-ifac
	iee=iee+ifac*idif
	if(iee.le.0) then
		iee=iee+nel
	else if(iee.gt.nel) then
		iee=iee-nel
	end if
c
	end do
c
	ielemo=0        !no element found
c
	return
	end
c
c*************************************************************
c
	subroutine transn(ie,xinou,yinou,tanf,dl)
c
c advection of particle
c
c advection of particle with coordinates (xinou,yinou)
c ...velocities of element and time step dt computed in routine
c ...returns new coordinates and element in xinou,yinou,ie
c ...if (hydrodynamical absolute) time and element do not
c ...change, element parameters can be reused
c
c ie            element number of particle (entry and return)
c xinou,yinou   coordinates of particle
c tanf          advection is started at that time
c dl            maximum way the particle is advected
c
	implicit none
c
c arguments
	integer ie
	real xinou,yinou,tanf,dl
c parameter
	real eps
	integer itbas
	parameter (eps=5.e-3,itbas=20)
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real tact,dthy,dlfl,dlfd,ddif
	integer itaflt,itafdf
	integer nen3v(3,1)
	integer iwegv(1),ieltv(3,1)
	real xgv(1),ygv(1),ev(13,1)
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /tflfd/ tact,dthy,dlfl,dlfd,itaflt,itafdf,ddif
	common /nen3v/nen3v
	common /xgv/xgv,/ygv/ygv,/ev/ev
	common /iwegv/iwegv,/ieltv/ieltv
c local
	logical buv
	integer iee,ieo,ieold,k,iero
	integer ntrans
	real xp,yp,xpo,ypo
	real told,dtold,dt,trest,rel,fk
	real u(3),v(3),uh,vh,uvm
c       real uh,vh,uvm
c function
c	integer intvel
	real pllvel,elnew,plline
	real farget
c
	save ieold,told,dtold
	data dtold /0./
c
	if(ie.le.0.or.iwegv(ie).gt.0) return
c
		ntrans=0
		call delta(-(itbas+1))
		call iarinc(1)
	iee=ie
	ieo=ie
	iero=0
	xp=xinou
	yp=yinou
	trest=tact-tanf
	dt=dtold
	buv=.true.
c
c loop on rest time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
    1   continue
c
c set values if new element %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	if(tact.ne.told.or.iee.ne.ieold) then
c
			call iarinc(4)
			call delta(-(itbas+4))
		uvm=pllvel(iee,u,v)
			call delta(itbas+4)
c
		if(dthy*uvm.gt.dl) then
			dt=dl/uvm
		else
			dt=dthy
			if(dt.lt.trest) dt=trest
		end if
c               write(6,*) 'uvm,dt',uvm,dt
c
		buv=.true.
		dtold=dt
		told=tact
		ieold=iee
	end if
c
c       statistics on time step
c
	call farinc(1)
	call faradd(2,dt)
	if(dt.lt.farget(3)) call farset(3,dt)
c
c       interpolate velocity only if buv=.true., otherwise already done
c
	if(buv) then
c               k=intvel(iee,xp,yp,uh,vh)
		uh=0.
		vh=0.
		do k=1,3
		   fk=ev(k,iee)+xp*ev(k+3,iee)+yp*ev(k+6,iee)
		   uh=uh+u(k)*fk
		   vh=vh+v(k)*fk
		end do
		buv=.false.
c               write(6,*) 'new vel ',iee,uh,vh
	end if
c
c advect particle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	if(trest.lt.dt) dt=trest
c
	xpo=xp
	ypo=yp
c
	xp=xp+uh*dt
	yp=yp+vh*dt
c
c check if particle still in element and interpolate velocity %%%%%%%%%%
c
	uh=0.
	vh=0.
	do k=1,3
	   fk=ev(k,iee)+xp*ev(k+3,iee)+yp*ev(k+6,iee)
	   if(fk.lt.-eps) goto 7
	   uh=uh+u(k)*fk
	   vh=vh+v(k)*fk
	end do
c
c xp,yp is in old element, loop if trest>0 
c
	trest=trest-dt
	if(trest.gt.0.) goto 1
	goto 5
c
c xp,yp out of element, find new element and trest
c
    7   continue
	buv=.true.
	ieo=iee
	rel=elnew(iee,xp,yp,xpo,ypo)
	trest=trest-dt*rel
c
c control if particle is moving
c
	if(rel.eq.0.) then
		if(iero.eq.iee) then
c                       iero,ieo are elements                   
			iee=iero
			iero=0
			rel=plline(iee,ieo,xp,yp,dt)
                        if(rel.ge.0.) trest=trest-dt*rel   !??? $$new13.12.91
			trest=trest-dt*rel
			call lgrerr(6,0,iee,xp,yp)
			if(rel.eq.0.) call lgrerr(9,0,iee,xp,yp)
c                       write(ner,*) 'rel : ',rel
		else
			iero=ieo
		end if
	else
		iero=0
	end if
c
	ntrans=ntrans+1
	if(ntrans.gt.100) then
		call lgrerr(5,0,iee,xp,yp)
		call intval(ie,xp,yp,u,v,uh,vh)
c               write(ner,*) 'u/v(xp,yp) : ',uh,vh
c               write(ner,'(1x,a,6f10.3)') 'u/v : ',(u(k),v(k),k=1,3)
		rel=-1.
	end if
c
c loop only if no error in elnew (rel >= 0)
c
	if(rel.ge.0.) goto 1    
c
c end of loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
    5   continue
	xinou=xp
	yinou=yp
	ie=iee
			call delta(itbas+1)
c
	return
	end
c
c*************************************************************
c
	subroutine diffus(ie,xinou,yinou,tanf)
c
c diffusion of particle with the quasi particle method
c
c ie            element number of particle (entry and return)
c xinou,yinou   coordinates of particle
c tanf          diffusion is started at that time
c
	implicit none
c
c parameter
	integer itbas
	parameter (itbas=20)
c arguments
	integer ie
	real xinou,yinou,tanf
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer itanf,itend,idt,nits,niter,it
	real tact,dthy,dlfl,dlfd,ddif
	integer itaflt,itafdf
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /femtim/ itanf,itend,idt,nits,niter,it
	common /tflfd/ tact,dthy,dlfl,dlfd,itaflt,itafdf,ddif
c local
c	integer iee,ieo,iseed
	integer iee,iseed
	real xp,yp,xpo,ypo
	real trest,ppp,rel
c	real trest,dt,ppp,rel
c functions
	real ran1,elnew
c
	save iseed
	data iseed /-352677/
c
	if(ddif.le.0.) return
c
        if( iseed.le.0 ) trest = ran1(iseed)

		call iarinc(2)
		call delta(-(itbas+2))
	iee=ie
c	ieo=ie
	xp=xinou
	yp=yinou
	trest=tact-tanf
c	dt=idt
c
c formula for austausch coefficient : p = a / 12
c ...( a = 12 m**2/s  ~  p = 1 m/s )
c ...so for a=0.1 ==> ppp=8.33e-3 ; a=10. ==> ppp=0.833
c ...       a=15. ==> ppp=1.25    ; a=20. ==> ppp=1.67
c
c formula from Kurt Duwe : p = sqrt(6*a/dt)/dt (why dt in formula ?)
c
c       aust=0.1
c       ppp=sqrt(6.*aust/dt)/dt
c
	ppp=ddif
c
c loop on rest time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	do while(trest.gt.0.)
c
c diffuse particle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	xpo=xp
	ypo=yp
c
	xp=xp+trest*ppp*(ran1(iseed)-0.5)
	yp=yp+trest*ppp*(ran1(iseed)-0.5)
c
	rel=elnew(iee,xp,yp,xpo,ypo)
c
	if(rel.lt.0.) rel=1.
	trest=trest*(1.-rel)
c
	end do
c
c end of loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	xinou=xp
	yinou=yp
	ie=iee
		call delta(itbas+2)
c
	return
	end
c
c********************************************************
c
	subroutine relran(nfddim)
c
c release random particles
c
	implicit none
c
c arguments
	integer nfddim
c parameter
	integer itbas
	parameter (itbas=20)
	integer ibndim
	parameter (ibndim=100)
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer itanf,itend,idt,nits,niter,it
	real tact,dthy,dlfl,dlfd,ddif
	integer itaflt,itafdf
	integer nfloat,nfdoat
	integer iefdv(1),irv(1)
	real xpfdv(1),ypfdv(1),bnd(ibndim,1)
	real rhv(1),rlv(1)
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /femtim/ itanf,itend,idt,nits,niter,it
	common /tflfd/ tact,dthy,dlfl,dlfd,itaflt,itafdf,ddif
	common /ffloat/ nfloat,nfdoat
	common /iefdv/iefdv,/irv/irv
	common /xpfdv/xpfdv,/ypfdv/ypfdv
	common /bnd/bnd
	common /rhv/rhv,/rlv/rlv
c local
	integer irand,ifdac,ibc,kranf,krend,iee,i
	real water,rlenal,portal,x,y
	real tddt,tnow,tend
c functions
	integer ielemf,iround
	real port,getpar
c
	save irand
	data irand /34279537/
c
c still to do :
c       counter for water that has not been used in one time step
c       ...to be added to next time step
c       introduce mode to switch between different release techniques
c       ...e.g. on point or in triangle per water or per time...
c       read parameters through nls routines
c
	water=getpar('watfdf')
	if(water.le.0.) return
c
			call delta(-(itbas+3))
			call iarinc(3)
c       water=20000.
c
	ifdac=nfdoat
c
	do 34 ibc=1,nbc
c
	   kranf=iround(bnd(3,ibc))
	   krend=iround(bnd(4,ibc))
c
c          compute total length of boundary line (rlenal) & discharge (portal)
c
	   rlenal=0.
	   portal=0.
	   do 70 i=kranf,krend-1
		rlenal=rlenal+rlv(i)
		portal=portal+port(irv(i),irv(i+1),rhv(i),rhv(i+1))
   70      continue
c
c          enough water for at least one particle ?
c
	   if(portal*idt.lt.water) goto 34
c
c          release from tnow to tend only at inflow with time step tddt
c
	   tddt=water/portal
	   tnow=it-idt+tddt
	   tend=it
c
c          loop to release particles
c
   75      continue
		call ranxy(rlenal,kranf,krend,irand,x,y,iee)
		iee=ielemf(iee,x,y)
		if(iee.gt.0) call transn(iee,x,y,tnow,dlfd)
c
		if(iee.gt.0.) then
			ifdac=ifdac+1
			if(ifdac.gt.nfddim) goto 99
			iefdv(ifdac)=iee
			xpfdv(ifdac)=x
			ypfdv(ifdac)=y
		end if  
c
		tnow=tnow+tddt
	   if(tnow.lt.tend) goto 75
c
   34   continue
c
	nfdoat=ifdac
		call delta(itbas+3)
c
	return
   99   continue
	write(6,*) 'Too many particles for dimension : ',nfddim
	stop 'error stop : relran'
	end
c
c********************************************************
c
	function port(k1,k2,h1,h2)
c
c compute discharge through line between k1 & k2
c
c k1,k2         node defining line
c h1,h2         depth at node
c port          discharge through line (return value)
c
	implicit none
c
	real port
c argument
	integer k1,k2
	real h1,h2
c parameter
	real t24
	parameter (t24=1./24.)
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real xv(3,1),xgv(1),ygv(1)
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /xv/xv,/xgv/xgv,/ygv/ygv
c local
	real u1,u2,v1,v2,z1,z2,uv1,uv2
	real x,y
c
	x=xgv(k2)-xgv(k1)
	y=ygv(k1)-ygv(k2)
	u1=xv(1,k1)+xv(1,k1+nkn)
	u2=xv(1,k2)+xv(1,k2+nkn)
	v1=xv(2,k1)+xv(2,k1+nkn)
	v2=xv(2,k2)+xv(2,k2+nkn)
	z1=xv(3,k1)+xv(3,k1+nkn)+2.*h1
	z2=xv(3,k2)+xv(3,k2+nkn)+2.*h2
c
	uv1 = u1*y + v1*x
	uv2 = u2*y + v2*x
c
	port = t24 * (uv1*(2.*z1+z2) + uv2*(z1+2.*z2))
c
	return
	end
c
c********************************************************
c
	subroutine ranxy(rlenal,kranf,krend,irand,x,y,iee)
c
c gets x,y for a random particle at boundary
c
c rlenal        total length of boundary line
c kranf,krend   start and end of bounday line in irv()
c irand         seed for random function
c x,y           coordinate of particle
c iee           probable element of particle
c
c the particle is put eps2 inside of the boundary line to
c       ...avoid loosing it at first step
c
	implicit none
c
c arguments
	integer kranf,krend,irand,iee
	real rlenal,x,y
c parameters
	real eps2
	parameter (eps2=1.e-5)
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer irv(1),ierv(2,1)
	real rlv(1)
	real xgv(1),ygv(1)
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /irv/irv,/ierv/ierv
	common /rlv/rlv
	common /xgv/xgv,/ygv/ygv
c local
	integer k1,k2,k
	real rrll,rel
	real x1,x2,y1,y2,dx,dy
c functions
	real ran1
c
	rrll=ran1(irand)*rlenal
c
c       find sector of boundary line
c
	do 76 k=kranf,krend-1
	  if(rrll.le.rlv(k)) goto 12
	  rrll=rrll-rlv(k)
   76   continue
c
	k=krend-1
	rrll=rlv(k)
c
   12   continue
c
c       (x,y) is in the k'th part rrll awayfrom start
c
	rel=rrll/rlv(k)
	k1=irv(k)
	k2=irv(k+1)
	x1=xgv(k1)
	y1=ygv(k1)
	x2=xgv(k2)
	y2=ygv(k2)
	dx=x2-x1
	dy=y2-y1
c
	x=x1+dx*rel-eps2*dy
	y=y1+dy*rel+eps2*dx
	iee=ierv(1,k)
c
	return
	end
c
c*************************************************************
c
	function pllvel(ie,u,v)
c
c get velocity for element and make parallel
c
c ie            element number of particle (entry)
c u,v           velocity components (returned)
c
	implicit none
c
	real pllvel
c arguments
	integer ie
	real u(3),v(3)
c parameter
	real frauv,frauv2
	parameter (frauv=0.05,frauv2=frauv/2.)
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real tact,dthy,dlfl,dlfd,ddif
	integer itaflt,itafdf
	integer nen3v(3,1)
	integer iwegv(1),ieltv(3,1)
	integer kantv(2,1),irv(1)
	real xv(3,1),xgv(1),ygv(1),ev(13,1)
	real dxv(1),dyv(1)
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /tflfd/ tact,dthy,dlfl,dlfd,itaflt,itafdf,ddif
	common /nen3v/nen3v
	common /iwegv/iwegv,/ieltv/ieltv
	common /kantv/kantv,/irv/irv 
	common /xv/xv,/xgv/xgv,/ygv/ygv,/ev/ev
	common /dxv/dxv,/dyv/dyv
c local
	integer io,iout,iouts,i,kn,iee
	integer k1,k2,k3
	real x0,y0,dx,dy
	real um,vm,uvm,u3,v3,uv3,un,vn,uvn,u1,v1,uv1,u2,v2,uv2
	real proj,proj1,proj2,projn,projnx,projny
	real x(3),y(3)
c	real fkf(3)
c statement functions
c bout  --> if adjacent element is dry or island
	logical is_internal_node
	logical bout
	bout(i,ie) = ieltv(i,ie).eq.0 .or. 
     +                  ieltv(i,ie).gt.0 .and. iwegv(ieltv(i,ie)).gt.0
c
	iee=ie
c	ixe=3*nkn
c
c       get coordinates and velocities
c
	um=0.
	vm=0.
	io=0
	iout=0
	iouts=0
	do i=1,3
	   kn=nen3v(i,iee)
	   if(bout(i,iee)) then
		iout=iout+1
		iouts=iouts+i
	   else if( .not. is_internal_node(kn) ) then
		io=i
	   end if
	   x(i)=xgv(kn)
	   y(i)=ygv(kn)
	   u(i)=(xv(1,kn)+xv(1,kn+nkn))*.5
	   v(i)=(xv(2,kn)+xv(2,kn+nkn))*.5
	   um=um+u(i)
	   vm=vm+v(i)
	end do
	uvm=sqrt(um*um+vm*vm)/3.
c
c       project velocities if boundary element
c
	if(iout.eq.0) then
c               inner element, control if exactly one node on boundary
c
		if(io.eq.0) goto 1
c
c               control if open boundary node
c               (...this part should be faster...)
c
		kn=nen3v(io,iee)
c
		do i=1,nrb
		  if(irv(i).eq.kn) goto 1
		end do
c
c               project velocities if necessary
c               (...if two boundary nodes without any boundary side ?)
c               (...this algorithm works only for one boundary node)
c
		x0=xgv(kn)
		y0=ygv(kn)
		dx=dxv(kn)
		dy=dyv(kn)
c
		i=mod(io,3)+1
		k1=nen3v(i,ie)
		if((xgv(k1)-x0)*(-dy)+(ygv(k1)-y0)*dx.lt.0.) then
		  k1=kantv(1,kn)
		  u3=xgv(k1)-x0
		  v3=ygv(k1)-y0
		  proj=(u(io)*u3+v(io)*v3)/(u3*u3+v3*v3)
		  u(io)=proj*u3
		  v(io)=proj*v3
		  goto 1
		end if
c
		i=mod(i,3)+1
		k1=nen3v(i,ie)
		if((xgv(k1)-x0)*(-dy)+(ygv(k1)-y0)*dx.lt.0.) then
		  k1=kantv(2,kn)
		  u3=xgv(k1)-x0
		  v3=ygv(k1)-y0
		  proj=(u(io)*u3+v(io)*v3)/(u3*u3+v3*v3)
		  u(io)=proj*u3
		  v(io)=proj*v3
		  goto 1
		end if
	else if(iout.eq.3) then
c               3 boundary element
		uvm=0.
		do 4 i=1,3
		  u(i)=0.
		  v(i)=0.
    4           continue
	else if(iout.eq.1) then
c               1 boundary element
c               iouts is external side (k3)
		k1=mod(iouts,3)+1
		k2=mod(k1,3)+1
		u3=x(k2)-x(k1)
		v3=y(k2)-y(k1)
		uv3=1./(u3*u3+v3*v3)
		proj1=(u(k1)*u3+v(k1)*v3)*uv3
		proj2=(u(k2)*u3+v(k2)*v3)*uv3
		projn=((v(k1)+v(k2))*u3-(u(k1)+u(k2))*v3)
     +                          *uv3*frauv2
c               let velocity point into element
		if(projn.gt.0) then
			projnx=-v3*projn
			projny= u3*projn
		else
			projnx= v3*projn
			projny=-u3*projn
		end if
		u(k1)=proj1*u3+projnx
		v(k1)=proj1*v3+projny
		u(k2)=proj2*u3+projnx
		v(k2)=proj2*v3+projny
	else
c               2 boundary element
c               k3 is internal side
		k3=6-iouts
		k1=mod(k3,3)+1
		k2=mod(k1,3)+1
c               from outer point (un,vn) points out off element
		un=y(k1)-y(k2)
		vn=x(k2)-x(k1)
		uvn=un*un+vn*vn
		projn=frauv*(u(k3)*un+v(k3)*vn)/uvn
c               let velocity point into element
		if(projn.gt.0.) then
			u(k3)=-projn*un
			v(k3)=-projn*vn
		else
			u(k3)=projn*un
			v(k3)=projn*vn
		end if
c               inner points
		u1=x(k3)-x(k2)
		v1=y(k3)-y(k2)
		uv1=u1*u1+v1*v1
		proj1=(u(k2)*u1+v(k2)*v1)/uv1
		u(k2)=proj1*u1
		v(k2)=proj1*v1
		u2=x(k1)-x(k3)
		v2=y(k1)-y(k3)
		uv2=u2*u2+v2*v2
		proj2=(u(k1)*u2+v(k1)*v2)/uv2
		u(k1)=proj2*u2
		v(k1)=proj2*v2
	end if
c
    1   pllvel=uvm
c
	return
	end
c
c*************************************************************
c
	subroutine intval(ie,xp,yp,u,v,uh,vh)
c
c interpolates u/v onto a point xp,yp in element ie
c
c arguments
	integer ie
	real xp,yp
	real u(3),v(3)
	real uh,vh
c common
	real ev(13,1)
	common /ev/ev
c local
	integer k
	real fk
c
	uh=0.
	vh=0.
	do k=1,3
	   fk=ev(k,ie)+xp*ev(k+3,ie)+yp*ev(k+6,ie)
	   uh=uh+u(k)*fk
	   vh=vh+v(k)*fk
	end do
c
	return
	end
c
c************************************************************
c
	subroutine lgrerr(mode,ip,ie,xp,yp)
c
c error output routine
c
	implicit none
c
c argument
	integer mode,ip,ie
	real xp,yp
c parameter
	integer nerr
	parameter (nerr=10)
c common
	integer itanf,itend,idt,nits,niter,it
	common /femtim/ itanf,itend,idt,nits,niter,it
c local
	logical bfull
	integer i,ierr(nerr),jtot,jtotpz,lgrmax,iset
c function
	integer ieext,iround
	real getpar
c data
	save jtot,jtotpz,ierr,lgrmax,iset
	data jtot,jtotpz,ierr,iset /0,0,nerr*0,0/
	data bfull /.false./
c
c lgrmax        0       only stat at end
c               1       full error code and stat at end
c               2       stat every time step
c               >2      stat every lgrmax errors and every time step
c
        if(iset.eq.0) then
                iset=1
                lgrmax=iround(getpar('lgrerr'))
                if(lgrmax.eq.1) bfull=.true.
                if(lgrmax.le.2) lgrmax=-1
        end if

	if(mode.eq.1) then
	    if(bfull) then
		write(ner,*) '$ (1) Particle ',ip,' lost at first step'
		write(ner,*) 'ie,xp,yp : ',ieext(ie),xp,yp
	    end if
	else if(mode.eq.2) then
	    if(bfull) then
		write(ner,*) '$ (2) Particle ',ip,' lost at time ',it
     +                          ,'through open boundary'
		write(ner,*) 'ie,xp,yp : ',ieext(ie),xp,yp
	    end if
	else if(mode.eq.3) then
	    if(bfull) then
		write(ner,*) '$ (3) Particle ',ip,' lost at time ',it
		write(ner,*) 'ie,xp,yp : ',ieext(ie),xp,yp
	    end if
	else if(mode.eq.4) then
	    if(bfull) then
		write(ner,*) '$ (4) Cannot find element'
		write(ner,*) 'ie,xp,yp : ',ieext(ie),xp,yp
	    end if
	else if(mode.eq.5) then
	    if(bfull) then
		write(ner,*) '$ (5) Particle is not moving, end advection'
		write(ner,*) 'ie,xp,yp : ',ieext(ie),xp,yp
	    end if
	else if(mode.eq.6) then
	    if(bfull) then
		write(ner,*) '$ (6) Particle advected along line'
		write(ner,*) 'ie,xp,yp : ',ieext(ie),xp,yp
	    end if
	else if(mode.eq.7) then
	    if(bfull) then
		write(ner,*) '$ (7) Adjacent element is the same'
		write(ner,*) 'ie : ',ieext(ie)
	    end if
	else if(mode.eq.8) then
	    if(bfull) then
		write(ner,*) '$ (8) Old coordinates assumed'
		write(ner,*) 'ie,xp,yp : ',ieext(ie),xp,yp
	    end if
	else if(mode.eq.9) then
	    if(bfull) then
		write(ner,*) '$ (9) Particle advected along line, rel=0.'
		write(ner,*) 'ie,xp,yp : ',ieext(ie),xp,yp
	    end if
        else if(mode.eq.10) then
c           if(bfull) then
                write(ner,*) '$ (10) Elements are not adjacent'
                write(ner,*) 'ie1,ie2 : ',ieext(ip),ieext(ie)
c           end if
	else if(mode.ne.0) then
	    if(bfull) then
		write(ner,*) '$ (0) Unknown error code : ',mode
	    end if
	end if
c
        if(mode.ge.1.and.mode.le.nerr) then
                ierr(mode)=ierr(mode)+1
                jtot=jtot+1
                jtotpz=jtotpz+1
        end if
c
        if(mode.eq.0.or.jtotpz.eq.lgrmax) then
              write(ner,*) '$ (0) Summary of error codes : '
     +                  ,jtot,it
              write(ner,'((1x,7i10))') (ierr(i),i=1,nerr)
              if(mode.eq.0) then
                jtot=0
                do i=1,nerr
                  ierr(i)=0
                end do
              end if
              jtotpz=0
        end if
c
	return
	end
c
c*************************************************************
c
	function plline(ie1,ie2,xp,yp,dt)
c
c advect on line
c
c ie1,ie2       element numbers of particle (entry)
c xp,yp         particle coordinates (emtry and return)
c dt            time step (entry)
c plline        rel of dt
c
	implicit none
c
	real plline
c arguments
	integer ie1,ie2
	real xp,yp,dt
c common
	integer nen3v(3,1)
	integer ieltv(3,1)
	real xv(3,1),xgv(1),ygv(1)
	common /nen3v/nen3v
	common /xv/xv,/xgv/xgv,/ygv/ygv
	common /ieltv/ieltv
c local
	integer k1,k2,i
	real dx,dy,dxy,relx,rely,rel,proj,u,v,xpo,ypo
c functions
	integer ieext
	real elnew
c
c new on 2.5.91 to avoid error stop for ie1=ie2
c
	if(ie1.eq.ie2) then
c               stop advection
		call lgrerr(7,i,ie1,xpo,ypo)
		plline=1.
		return
	end if
c
	do i=1,3
	   if(ieltv(i,ie1).eq.ie2) goto 1
	end do
c
c advection should be stopped here, but lets go on anyway  $$new25.09.91
c
        call lgrerr(10,ie1,ie2,xpo,ypo)
        plline=1.
        return

c	write(6,*) 'elements are not adjacent : ',ieext(ie1),ieext(ie2)
c	stop 'error stop : plline'
c
    1   continue
	k1=nen3v(mod(i,3)+1,ie1)
	k2=nen3v(mod(i+1,3)+1,ie1)
	dx=xgv(k2)-xgv(k1)
	dy=ygv(k2)-ygv(k1)
	relx=(xp-xgv(k1))/dx
	rely=(yp-ygv(k1))/dy
	rel=0.5*(relx+rely)
	u=rel*xv(1,k2)+(1.-rel)*xv(1,k1)
	v=rel*xv(2,k2)+(1.-rel)*xv(2,k1)
	dxy=1./(dx*dx+dy*dy)
	proj=(u*dx+v*dy)*dxy
	u=proj*dx
	v=proj*dy
	xpo=xp
	ypo=yp
	xp=xp+u*dt
	yp=yp+v*dt
c
	plline=elnew(ie1,xp,yp,xpo,ypo)
c
	return
	end
c
c********************************************************
c
	subroutine relra1(nfddim)
c
c release random particles at one node
c
	implicit none
c
c arguments
	integer nfddim
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer itanf,itend,idt,nits,niter,it
	real tact,dthy,dlfl,dlfd,ddif
	integer itaflt,itafdf
	integer nfloat,nfdoat
	integer iefdv(1)
	real xpfdv(1),ypfdv(1)
	real xgv(1),ygv(1)
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /femtim/ itanf,itend,idt,nits,niter,it
	common /tflfd/ tact,dthy,dlfl,dlfd,itaflt,itafdf,ddif
	common /ffloat/ nfloat,nfdoat
	common /iefdv/iefdv
	common /xpfdv/xpfdv, /ypfdv/ypfdv
	common /xgv/xgv, /ygv/ygv
c local
	logical bmanu
	integer nnode,node,iee,iee1,ifdac,ittanf,ittend
	real x,y,ppsec,x1,y1
	real tanf,trest,ddtt
c	real f(10)
c functions
	integer ielemf,iround,ipext,ieext
	real getpar
c data
	save nnode,ddtt,x1,y1,iee1,trest
	save ittanf,ittend
	data bmanu /.false./    !give x1,y1... in subroutine
c       data bmanu /.true./     !give x1,y1... in subroutine
	data nnode,iee1,trest /0,0,0./
c       data x1,y1,ittanf,ittend /28000.,6000.,0,172800/        !marghera
c	data x1,y1,ittanf,ittend /22250.,5000.,0,172800/        !petroli
c       data x1,y1,ittanf,ittend /40000.,625000.,0,172800/      !adria,ven.
c       data x1,y1,ittanf,ittend /100000.,660000.,0,172800/     !adria,trs.
c
c       get internal node number %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	ppsec=getpar('ppsfdf')
	if(ppsec.le.0.) return
c
	if(bmanu.and.iee1.le.0) then
		iee1=ielemf(0,x1,y1)
		if(iee1.le.0) then
			write(6,*) 'Cannot find element'
			stop 'relra1 : error stop'
		end if
		ddtt=1./ppsec
		nnode=-1
		write(6,*) 'release for fdf : ',ieext(iee1)
     +                          ,ddtt,x1,y1
		write(6,*)
	else if(nnode.eq.0) then
		node=iround(getpar('nodfdf'))
		nnode=node
		nnode = ipext(nnode)
		if(nnode.le.0) then
			write(6,*) 'Cannot find node ',node
			stop 'relra1 : error stop'
		end if
		ddtt=1./ppsec
		x1=xgv(nnode)
		y1=ygv(nnode)
		iee1=ielemf(0,x1,y1)
		if(iee1.le.0) then
			write(6,*) 'Cannot find element'
			stop 'relra1 : error stop'
		end if
                ittanf=iround(getpar('itsfdf'))
                ittend=iround(getpar('itefdf'))
                if(ittanf.eq.0.and.ittend.eq.0) then
                        ittanf=itanf
                        ittend=itend
                end if
		write(6,*) 'relra1 : '
		write(6,*) 'node for fdf : ',ipext(nnode),ieext(iee1)
     +                          ,ddtt,x1,y1
		write(6,*) ittanf,ittend
		write(6,*)
	end if
c
c
c	if(bmanu.and.(it.lt.ittanf.or.it.gt.ittend)) return
	if(it.lt.ittanf.or.it.gt.ittend) return
c
	ifdac=nfdoat
c
c       release %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	trest=trest+idt
	tanf=it-trest
	do while(trest-ddtt.ge.0.)
	   iee=iee1
	   x=x1
	   y=y1
	   call transn(iee,x,y,tanf,dlfd)
	   if(iee.gt.0.) then
		ifdac=ifdac+1
		if(ifdac.gt.nfddim) goto 99
		iefdv(ifdac)=iee
		xpfdv(ifdac)=x
		ypfdv(ifdac)=y
	   end if       
	   trest=trest-ddtt
	   tanf=tanf+ddtt
	end do
c
	nfdoat=ifdac
c
	return
   99   continue
	write(6,*) 'Too many particles for dimension : ',nfddim
	stop 'error stop : relra1'
	end
c
c********************************************************
c
	subroutine relra2(nfddim)
c
c release random particles at one node
c
	implicit none
c
c arguments
	integer nfddim
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer itanf,itend,idt,nits,niter,it
	real tact,dthy,dlfl,dlfd,ddif
	integer itaflt,itafdf
	integer nfloat,nfdoat
	integer iefdv(1)
	real xpfdv(1),ypfdv(1)
	real xgv(1),ygv(1)
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /femtim/ itanf,itend,idt,nits,niter,it
	common /tflfd/ tact,dthy,dlfl,dlfd,itaflt,itafdf,ddif
	common /ffloat/ nfloat,nfdoat
	common /iefdv/iefdv
	common /xpfdv/xpfdv, /ypfdv/ypfdv
	common /xgv/xgv, /ygv/ygv
c local
	logical bmanu
	integer nnode,node,iee,iee1,ifdac,ittanf,ittend
	real x,y,ppsec,x1,y1
	real tanf,trest,ddtt
c	real f(10)
c functions
	integer ielemf,iround,ipext,ieext
	real getpar
c data
	save nnode,ddtt,x1,y1,iee1,trest
	save ittanf,ittend
	data bmanu /.false./    !give x1,y1... in subroutine
c       data bmanu /.true./     !give x1,y1... in subroutine
	data nnode,iee1,trest /0,0,0./
c       data x1,y1,ittanf,ittend /28000.,6000.,0,172800/        !marghera
c	data x1,y1,ittanf,ittend /22250.,5000.,0,172800/        !petroli
c       data x1,y1,ittanf,ittend /40000.,625000.,0,172800/      !adria,ven.
c       data x1,y1,ittanf,ittend /100000.,660000.,0,172800/     !adria,trs.
c
c       get internal node number %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	ppsec=getpar('ppsfd2')
	if(ppsec.le.0.) return
c
	if(bmanu.and.iee1.le.0) then
		iee1=ielemf(0,x1,y1)
		if(iee1.le.0) then
			write(6,*) 'Cannot find element'
			stop 'relra2 : error stop'
		end if
		ddtt=1./ppsec
		nnode=-1
		write(6,*) 'release for fdf : ',ieext(iee1)
     +                          ,ddtt,x1,y1
		write(6,*)
	else if(nnode.eq.0) then
		node=iround(getpar('nodfd2'))
		nnode=node
		nnode = ipext(nnode)
		if(nnode.le.0) then
			write(6,*) 'Cannot find node ',node
			stop 'relra2 : error stop'
		end if
		ddtt=1./ppsec
		x1=xgv(nnode)
		y1=ygv(nnode)
		iee1=ielemf(0,x1,y1)
		if(iee1.le.0) then
			write(6,*) 'Cannot find element'
			stop 'relra2 : error stop'
		end if
                ittanf=iround(getpar('itsfd2'))
                ittend=iround(getpar('itefd2'))
                if(ittanf.eq.0.and.ittend.eq.0) then
                        ittanf=itanf
                        ittend=itend
                end if
		write(6,*) 'relra2 : '
		write(6,*) 'node for fdf : ',ipext(nnode),ieext(iee1)
     +                          ,ddtt,x1,y1
		write(6,*) ittanf,ittend
		write(6,*)
	end if
c
c
c	if(bmanu.and.(it.lt.ittanf.or.it.gt.ittend)) return
	if(it.lt.ittanf.or.it.gt.ittend) return
c
	ifdac=nfdoat
c
c       release %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	trest=trest+idt
	tanf=it-trest
	do while(trest-ddtt.ge.0.)
	   iee=iee1
	   x=x1
	   y=y1
	   call transn(iee,x,y,tanf,dlfd)
	   if(iee.gt.0.) then
		ifdac=ifdac+1
		if(ifdac.gt.nfddim) goto 99
		iefdv(ifdac)=iee
		xpfdv(ifdac)=x
		ypfdv(ifdac)=y
	   end if       
	   trest=trest-ddtt
	   tanf=tanf+ddtt
	end do
c
	nfdoat=ifdac
c
	return
   99   continue
	write(6,*) 'Too many particles for dimension : ',nfddim
	stop 'error stop : relra2'
	end

c**************************************************************************

c       integer function ckfdf(nb)
c
c       integer nb
c
c       ckfdf = -1
c
c       return
c       end
c
c       integer function rffdf(nb,nvers,i,idum,name)
c       rffdf=0
c       return
c       end
c
c       integer function rdfdf(nb,nvers,idum,nfdoat,xpfdv,ypfdv,iefdv)
c       rdfdf=0
c       return
c       end
