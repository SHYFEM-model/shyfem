c
c $Id: subcon.f,v 1.26 2003/03/25 14:08:55 georg Exp $
c
c routines for concentration
c
c contents :
c
c subroutine con2sh				shell for conz
c subroutine conz2d(cn,coe,caux,dt,chpar)	computes concentration
c subroutine conzbc(cn,coe,caux,cb,flag)	updates boundary conditions 
c subroutine conzin(iconz,conref,cn,coe)	sets up conz for start
c
c subroutine cmima(nkn,nel,cnv,coev,ckmin,ckmax,cemin,cemax) min/max of conz
c subroutine massconc2d(kvol,ce,res)            total mass of conc
c
c revision log :
c
c 16.06.1998	ggu	new names for some parameters
c			rkpar -> chpar , cref -> conref
c 19.06.1998	ggu	bug fix (ifileo nor declared)
c 20.06.1998	ggu	general cleanup
c 24.06.1998	ggu	heat flux enabled
c 26.01.1999	ggu	new routines (3D) to write files NOS
c 07.06.1999	ggu	bug fix -> bconz... not saved
c 28.10.1999	ggu	names changed
c 16.11.2001	ggu	some old routines deleted
c 11.10.2002	ggu	massconc renamed into massconc2d
c
c*********************************************************************
c
	subroutine con2sh
c
c shell for conz
c
c use negative values for iconz to run special tests
c
c  1	normal run
c  0	no transport/diffusion
c -1	diffusion test
c -2	advection test
c -3	flux test
c -4	dilution test
c
c written 09.01.94 by ggu  (from scratch)
c revised 20.01.94 by ggu  $$iclin - iclin not used to compute volume
c revised 20.01.94 by ggu  $$zov0 - zeov in conz adjusted
c revised 31.01.94 by ggu  $$conzbc - implementation of bc for conz
c revised 03.02.94 by ggu  $$zov00 - zeov in sp159 adjusted
c revised 04.02.94 by ggu  $$azpar - azpar used to compute transport
c revised 07.02.94 by ggu  $$istot - istot for fractional time step
c revised 03.12.97 by ggu  $$TS - temp/salt included
c
	implicit none
c
c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it
        real eps1,eps2,pi,flag,high,higi
        common /mkonst/ eps1,eps2,pi,flag,high,higi

        integer nen3v(3,1)
        common /nen3v/nen3v
	real zenv(3,1),zeov(3,1)
	common /zenv/zenv, /zeov/zeov
	real v1v(1),v2v(1)
	common /v1v/v1v, /v2v/v2v
	real coev(3,1), cnv(1), rcv(1)		!$$conzbc
	common /coev/coev, /cnv/cnv, /rcv/rcv	!$$conzbc
	real toev(3,1), tnv(1), rtv(1)		!$$TS
	common /toev/toev, /tnv/tnv, /rtv/rtv	!$$TS
	real soev(3,1), snv(1), rsv(1)		!$$TS
	common /soev/soev, /snv/snv, /rsv/rsv	!$$TS
c local
	integer icall
c	integer ie,ii,k
	real ckmax,ckmin,cemax,cemin
c	double precision conref,chpar,dt
	real conref,temref,salref
	real chpar,dt,azpar
        real gamma
	logical bdebug
	logical bconz,btemp,bsalt
	integer iconz,itemp,isalt
	integer istot,isact
	integer iu
	integer iuc,ius,iut
	integer itmcon,idtcon
c function
	real getpar
	integer iround
c save & data
	save icall
	save bdebug
	save chpar,dt,azpar
	save conref,temref,salref
	save iconz,itemp,isalt
	save bconz,btemp,bsalt
	save istot
	save iuc,ius,iut
	save itmcon,idtcon

	data icall /0/

	if(icall.eq.-1) return

c set initial values

	if(icall.eq.0) then
	  iconz=iround(getpar('iconz'))
	  itemp=iround(getpar('itemp'))
	  isalt=iround(getpar('isalt'))
	  bconz = iconz .ne. 0
	  btemp = itemp .ne. 0
	  bsalt = isalt .ne. 0

	  if( .not. ( bconz .or. btemp .or. bsalt ) ) icall = -1
	  if( icall .eq. -1 ) return

	  conref=getpar('conref')
	  temref=getpar('temref')
	  salref=getpar('salref')

	  chpar=getpar('chpar')
	  call getaz(azpar)
	  istot=iround(getpar('istot'))		!$$istot
	  dt=idt

	  bdebug = iround(getpar('levdbg')) .ge. 3

	  call conzin(iconz,conref,cnv,coev)
	  call conzin(itemp,temref,tnv,toev)
	  call conzin(isalt,salref,snv,soev)

	  if( btemp ) call qfluxr(it)		!initialize heat flux module

          iu = 55
          itmcon = iround(getpar('itmcon'))
          idtcon = iround(getpar('idtcon'))

	  iuc = iu
          if( bconz ) call confop(iuc,itmcon,idtcon,1,'con')
	  ius = iu
          if( bsalt ) call confop(ius,itmcon,idtcon,1,'sal')
	  iut = iu
          if( btemp ) call confop(iut,itmcon,idtcon,1,'tem')

          call diffstab(dt,chpar,istot,v1v,v2v,gamma)
	end if

	icall=icall+1

	if( bconz ) then
	  do isact=1,istot		!$$istot
	    call conz2d(cnv,coev,v1v,dt,chpar,azpar,istot,isact)	!$$azpar
            call conzbc(cnv,coev,v1v,rcv,flag,azpar) !boundary conditions
	  end do
	end if

	if( bsalt ) then
	  do isact=1,istot		!$$istot
	    call conz2d(snv,soev,v1v,dt,chpar,azpar,istot,isact)	!$$azpar
            call conzbc(snv,soev,v1v,rsv,flag,azpar) !boundary conditions
	  end do
	end if

	if( btemp ) then
	  do isact=1,istot		!$$istot
	    call conz2d(tnv,toev,v1v,dt,chpar,azpar,istot,isact)	!$$azpar
            call conzbc(tnv,toev,v1v,rtv,flag,azpar) !boundary conditions
	  end do
	  call qfluxr(it)		!heat fluxes -> read data
	  call qflux2d(it,idt)		!heat fluxes -> compute temperature
	end if

	if( bdebug ) then
	  if( bconz ) then
	    call cmima(nkn,nel,cnv,coev,ckmin,ckmax,cemin,cemax)
	    write(6,*) 'cemax : ',it,cemin,cemax
	  end if
	  if( bsalt ) then
	    call cmima(nkn,nel,snv,soev,ckmin,ckmax,cemin,cemax)
	    write(6,*) 'semax : ',it,cemin,cemax
	  end if
	  if( btemp ) then
	    call cmima(nkn,nel,tnv,toev,ckmin,ckmax,cemin,cemax)
	    write(6,*) 'temax : ',it,cemin,cemax
	  end if
	end if

	if( bconz ) then
          call confil(iuc,itmcon,idtcon,10,1,cnv)
c	  call conzrs(cnv,coev,v1v)	!FIXME
	end if
	if( bsalt ) then
          call confil(ius,itmcon,idtcon,11,1,snv)
	end if
	if( btemp ) then
          call confil(iut,itmcon,idtcon,12,1,tnv)
	end if

	end

c**************************************************************

        subroutine conz2d(cn,coe,caux,ddt,chpar,azpar,istot,isact)

c computes concentration
c
c cn    new concentration
c coe   old concentration
c caux  aux vector
c ddt    time step
c chpar dispersion parameter
c azpar time weighting parameter
c istot	total inter time steps
c isact	actual inter time step
c
c written 09.01.94 by ggu  (from scratch)
c revised 19.01.94 by ggu  $$flux - flux conserving property
c revised 20.01.94 by ggu  $$iclin - iclin not used to compute volume
c revised 20.01.94 by ggu  $$lumpc - evaluate conz. nodewise
c revised 03.02.94 by ggu  $$itot0 - exception for itot=0 or 3
c revised 04.02.94 by ggu  $$fact3 - factor 3 missing in transport
c revised 04.02.94 by ggu  $$azpar - azpar used to compute transport
c revised 04.02.94 by ggu  $$condry - comute conz also in dry areas
c revised 07.02.94 by ggu  $$istot - istot for fractional time step
c
c solution of purely diffusional part : 
c
c dC/dt = a*laplace(C)    with    c(x,0+)=delta(x)
c
c C(x,t) =  (4*pi*a*t)**(-n/2) * exp( -|x|**2/(4*a*t) )
c
c for n-dimensions and 
c
c C(x,t) =  1/sqrt(4*pi*a*t) * exp( -x**2/(4*a*t) )
c
c for 1 dimension
c
c the solution is normalized, i.e.  int(C(x,t)dx) = 1 over the whole area
c
	implicit none
c
c arguments
        real cn(1),coe(3,1),caux(1)
        real ddt,chpar,azpar			!$$azpar
	integer istot,isact
c parameters
	real drittl
	parameter (drittl=1./3.)
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer itanf,itend,idt,nits,niter,it		!for debug
        common /femtim/ itanf,itend,idt,nits,niter,it	!...
	integer nen3v(3,1)
	real hm3v(3,1)
	real ev(13,1)
	integer iwegv(1)
	real uov(1),vov(1)
	real unv(1),vnv(1)
	real zenv(3,1), zeov(3,1)
	common /nen3v/nen3v
	common /hm3v/hm3v
	common /ev/ev
	common /iwegv/iwegv
	common /uov/uov, /vov/vov
	common /unv/unv, /vnv/vnv
	common /zenv/zenv, /zeov/zeov
c local
	logical bdebug
	integer k,ie,ii
	integer itot,isum	!$$flux
	integer kspez,iespez
	real us,vs
	real zom,znm,hm,cm,cbm,ccm,hom
	real chzm,aux,dt
	real az,azt
	real aj,rk3
	real b(3),c(3),f(3)
	real zo(3),zn(3)
	real rso,rsn,rsot,rsnt,rstot
	double precision ctot,vtot
	double precision cttot
	double precision fctot
	integer ib
c
c	integer ipint,ieint
c
	bdebug=.false.

c	kspez=ipint(2940)
	kspez=-1
c	iespez=ieint(4336)
	iespez=-1

	az=azpar		!$$azpar
	azt=1.-az
	rk3=chpar*3.
c
	rstot=istot		!$$istot
	rso=(isact-1)/rstot
	rsn=(isact)/rstot
	rsot=1.-rso
	rsnt=1.-rsn
c
	dt=ddt/rstot
c	
c    1	continue	!for kerror	!FIXME
c
	ctot=0.d0
	vtot=0.d0
c
	cttot=0.d0
	fctot=0.d0
c
	do k=1,nkn
	  cn(k)=0.
	  caux(k)=0.
	end do
c
        do ie=1,nel
c
	us=az*unv(ie)+azt*uov(ie)		!$$azpar
	vs=az*vnv(ie)+azt*vov(ie)
c
	aj=ev(10,ie)    !area of triangle / 12
c
	if(ie.eq.iespez) then
		write(6,*) '* ',ie,iwegv(ie),aj
		write(6,*) '* ',us,vs
	end if
c
	zom=0.
	znm=0.
	hm=0.
	cm=0.
	cbm=0.
	ccm=0.
	chzm=0.
	itot=0
	isum=0
	ib=0
	do ii=1,3
	  if(nen3v(ii,ie).eq.kspez) ib=ii
	  b(ii)=ev(ii+3,ie)
	  c(ii)=ev(ii+6,ie)
c	  f(ii)=uso*b(ii)+vso*c(ii)	!$$flux
	  f(ii)=us*b(ii)+vs*c(ii)	!$$azpar
	  if(f(ii).lt.0.) then	!flux out of node
	    itot=itot+1
	    isum=isum+ii
	  end if
c	  zo(ii)=zeov(ii,ie)
c	  zn(ii)=zenv(ii,ie)
	  zo(ii)=rso*zenv(ii,ie)+rsot*zeov(ii,ie)		!$$istot
	  zn(ii)=rsn*zenv(ii,ie)+rsnt*zeov(ii,ie)
	  zom=zom+zo(ii)
	  znm=znm+zn(ii)
	  hm=hm+hm3v(ii,ie)
	  cm=cm+coe(ii,ie)
	  cbm=cbm+b(ii)*coe(ii,ie)
	  ccm=ccm+c(ii)*coe(ii,ie)
	  chzm=chzm+(hm3v(ii,ie)+zo(ii))*coe(ii,ie)
	end do
c
c itot=1 -> flux out of one node
c	compute flux with concentration of this node
c itit=2 -> flux into one node
c	for flux use conz. of the other two nodes and
c	minus the sum of these nodes for the flux of this node
c
	if(itot.eq.1) then	!$$flux
	  f(1)=f(1)*coe(isum,ie)
	  f(2)=f(2)*coe(isum,ie)
	  f(3)=f(3)*coe(isum,ie)
	else if(itot.eq.2) then
	  isum=6-isum
	  f(1)=f(1)*coe(1,ie)
	  f(2)=f(2)*coe(2,ie)
	  f(3)=f(3)*coe(3,ie)
	  f(isum) = 0.
	  f(isum) = -(f(1)+f(2)+f(3))
	else			!exception	$$itot0
	  f(1)=0.
	  f(2)=0.
	  f(3)=0.
	end if
	if(ie.eq.iespez) then
		write(6,*) '* ',f(1),f(2),f(3)
		write(6,*) '* ',coe(1,ie),coe(2,ie),coe(3,ie)
	end if
c
	hom=(hm+zom)*drittl	!$$iclin
	hm=hm*drittl		!$$lumpc
c
c	!look out for cttot
c	!if diffusion this is nearly useless (not taken into account)
	if(ib.gt.0) then
	  fctot=fctot+dt*aj*12.d0*f(ib)
	  cttot=cttot+4.d0*aj*(zo(ib)+hm)*coe(ib,ie)
	end if
c
	ctot=ctot+(aj*chzm)
	vtot=vtot+(aj*hom)
c
	if(iwegv(ie).eq.0) then		!wet element
	 do ii=1,3
	  k=nen3v(ii,ie)
	  cn(k) = cn(k) + aj * ( (hm+zo(ii))*coe(ii,ie)
     +				+ dt *  ( 3.*f(ii)		!$$fact3
     +					- b(ii)*rk3*hom*cbm
     +					- c(ii)*rk3*hom*ccm
     +					)
     +			       )  
	  caux(k) = caux(k) + aj * ( hm + zn(ii) )
	 end do
	else				!dry element
	 do ii=1,3						!$$condry
	  aux =                ( (hm+zo(ii))*coe(ii,ie)
     +				+ dt *  ( 3.*f(ii)		!$$fact3
     +					- b(ii)*rk3*hom*cbm
     +					- c(ii)*rk3*hom*ccm
     +					)
     +			       )  
	  coe(ii,ie) = aux / ( hm + zn(ii) )
	 end do
	end if
c
	end do
c
c compute concentration for each node
c
	do k=1,nkn
	  if(caux(k).gt.0.) then
	    cn(k)=cn(k)/caux(k)
	  end if
	end do
c
c put back concentration to element
c
	do ie=1,nel
	  if(iwegv(ie).eq.0) then
	    do ii=1,3
	      coe(ii,ie)=cn(nen3v(ii,ie))
	    end do
	  else	!interpolation into dry nodes
	    aj=ev(10,ie)
	    do ii=1,3
	      k=nen3v(ii,ie)
	      if(caux(k).le.0.) then
		caux(k)=caux(k)-aj
		cn(k)=cn(k)+aj*coe(ii,ie)
	      end if
	    end do
	  end if
	end do
c
	do k=1,nkn
	  if(caux(k).lt.0.) then
	    cn(k) = -cn(k)/caux(k)
	  end if
	end do
c
	do ie=1,nel
	  do ii=1,3
	    if(nen3v(ii,ie).eq.kspez) then
	      aj=ev(10,ie)
	      zn(ii)=rsn*zenv(ii,ie)+rsnt*zeov(ii,ie)
	      cttot=cttot-4.d0*aj*(zn(ii)+hm3v(ii,ie))*coe(ii,ie)
	    end if
	  end do
	end do
c
	if(bdebug) then
	  write(6,*) '************************************'
	  write(6,*) 'conz : ',4.*ctot,12.*vtot
	  write(6,*) fctot,cttot,fctot+cttot
	  !if diffusion last line is nearly useless (not taken into account)
	end if
c
	return
	end
c
c**************************************************************
c
        subroutine conzbc(cn,coe,caux,cb,flag,azpar)
c
c updates boundary conditions and coe for concentration
c
c cn    new concentration
c coe   old concentration
c caux  aux vector
c cb	boundary values
c flag	flag for not-valid boundary values
c
c written 31.01.94 by ggu  (from scratch)
c revised 04.02.94 by ggu  $$azpar - azpar used to compute transport
c
	implicit none
c
c arguments
        real cn(1),coe(3,1),caux(1),cb(1)
        real flag,azpar
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nen3v(3,1)
	real ev(13,1)
	integer iwegv(1),ieltv(3,1)
	real uov(1),vov(1)
	real unv(1),vnv(1)
	common /nen3v/nen3v
	common /ev/ev
	common /iwegv/iwegv, /ieltv/ieltv
	common /uov/uov, /vov/vov
	common /unv/unv, /vnv/vnv
c local
	integer k,ie,ii,i
	integer ibstot,ibntot,isum
	real aj,ff,cc
	real az,azt
	real us,vs
c
	az=azpar		!$$azpar
	azt=1.-az
c
	do k=1,nkn
	  if(cb(k).ne.flag) cn(k)=0.
	  caux(k)=0.
	end do
c
        do ie=1,nel
c
	if(iwegv(ie).gt.0) goto 1
c
	ibstot=0	!sides with open boundary
	ibntot=0	!nodes at open boundary
	isum=0		!to find out which node is which
	do ii=1,3
	  if(ieltv(ii,ie).eq.-1) ibstot=ibstot+1 !open bc
	  if(cb(nen3v(ii,ie)).ne.flag) then
		ibntot=ibntot+1
		isum=isum+ii
	  end if
	end do
c
	if(ibntot.eq.0) goto 1
c
c	maybe useless checks, but better be sure...
c	...we cannot handle two open boundaries in one element
c	...neither can we two open boundary nodes without open
c	...boundary side
c
	if(ibntot.eq.3.or.ibstot.gt.1) goto 99
	if(ibstot.eq.1.and.ibntot.ne.2) goto 99
	if(ibstot.eq.0.and.ibntot.ne.1) goto 99
c
c	now we are sure that we have 
c		- 1 boundary side with two boundary nodes
c		- 0 boundary side with one boundary node
c
	aj=ev(10,ie)
	us=az*unv(ie)+azt*uov(ie)		!$$azpar
	vs=az*vnv(ie)+azt*vov(ie)
c
	if(ibstot.eq.0) then	!no boundary side
	  i=isum		!i is node on boundary
	  k=nen3v(i,ie)		!k is boundary node
	  ff=us*ev(i+3,ie)+vs*ev(i+6,ie)
	  if(ff.gt.0) then	!flux into node -> average 2 internal nodes
	    cc=0.5*(coe(1,ie)+coe(2,ie)+coe(3,ie)-coe(i,ie))
	  else			!flux out of node -> impose bc
	    cc=cb(k)
	  end if
	  cn(k)=cn(k)+cc*aj	!sum to boundary values
	  caux(k)=caux(k)+aj
	else			!one boundary side (two nodes)
	  i=6-isum		!now i is node not (!) on boundary
	  k=nen3v(i,ie)		!k is internal node
	  ff=us*ev(i+3,ie)+vs*ev(i+6,ie)
	  if(ff.gt.0) then	!flux into node -> ingoing flow
	    cc=flag		!flag for later
	  else			!outgoing flow
	    cc=coe(i,ie)
	  end if
	  do ii=1,3		!sum to boundary values
	    if(i.ne.ii) then
		k=nen3v(ii,ie)
		if(cc.eq.flag) then	!value was flagged -> take bc
		  cn(k)=cn(k)+cb(k)*aj
		else
		  cn(k)=cn(k)+cc*aj
		end if
		caux(k)=caux(k)+aj
	    end if
	  end do
	end if
c	write(80,*) ie,ibntot,ibstot,isum,cc
c
    1	continue
	end do
c
	do k=1,nkn
	  if(caux(k).gt.0.) then
	    cn(k)=cn(k)/caux(k)
	  end if
	end do
c
c just copy everything again
c
	do ie=1,nel
	  if(iwegv(ie).eq.0) then
	    do ii=1,3
	      coe(ii,ie)=cn(nen3v(ii,ie))
	    end do
	  end if
	end do
c
	return
   99	continue
	write(6,*) ie,ibstot,ibntot,isum
	stop 'error stop conzbc : invalid boudary element'
	end

c*****************************************************************

	subroutine conzin(iconz,conref,cn,coe)

c sets up conz for start and testing

	implicit none

	integer iconz
	real conref
	real cn(1),coe(3,1)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer nen3v(3,1)
        common /nen3v/nen3v

	integer k,ie,ii
	integer kk1,kk2

	  if(iconz.gt.0) then
	    kk1=1	!initialize whole domain
	    kk2=nkn
          else if(iconz.eq.-1) then
            kk1=109     !one-dim point source
            kk2=117     !...
          else if(iconz.eq.-2) then
            kk1=109     !one-dim blob
            kk2=135     !...
            kk2=117     !...
          else if(iconz.eq.-3) then
            kk1=0               !flux test
            kk2=0               !...
          else if(iconz.eq.-4) then
            kk1=113	!poin source
            kk2=113
	  else
	    kk1=0
	    kk2=0
          end if

          do k=1,nkn
            cn(k)=0.
          end do
          do k=kk1,kk2
            cn(k)=conref
          end do
          do ie=1,nel
            do ii=1,3
              coe(ii,ie)=cn(nen3v(ii,ie))
            end do
          end do

	end

c***********************************************************

	subroutine cmima(nkn,nel,cnv,coev,ckmin,ckmax,cemin,cemax)

c min/max of conz

	implicit none

	integer nkn,nel
	real cnv(1), coev(3,1)
	real ckmin,ckmax,cemin,cemax

	integer k,ie,ii

	ckmax=cnv(1)
	ckmin=ckmax
	do k=1,nkn
	  if(cnv(k).gt.ckmax) ckmax=cnv(k)
	  if(cnv(k).lt.ckmin) ckmin=cnv(k)
	end do

	cemax=coev(1,1)
	cemin=cemax
	do ie=1,nel
	  do ii=1,3
	    if(coev(ii,ie).gt.cemax) cemax=coev(ii,ie)
	    if(coev(ii,ie).lt.cemin) cemin=coev(ii,ie)
	  end do
	end do

	end

c**************************************************************

        subroutine massconc2d(kvol,ce,res)
 
c computes total mass of conc (element values) (2d version)
 
        implicit none
 
c arguments
	integer kvol
        real ce(3,1)
        real res
c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real zeov(3,1),zenv(3,1)
        common /zeov/zeov, /zenv/zenv
	integer nen3v(3,1)
	common /nen3v/nen3v
c local
        integer ie,ii,l,mode
        integer k
        integer nlev,n
        integer iespec
	integer ibase
        real area,h
        double precision sum,mass
        real areaele,depele
 
        iespec = -1
 
        mode = 0        ! no levels
        mass = 0.
 
        do ie=1,nel
          area = areaele(ie) / 3.
          h = depele(ie,mode)
          sum = 0.
          do ii=1,3
	    k = nen3v(ii,ie)
	    if( kvol .le. 0 .or. kvol .eq. k ) then
	!if(kvol.gt.0) then
	!  write(6,*) 'massconc2d element... ',ie,ii,k,h + zenv(ii,ie)
	!end if
	      sum = sum + ce(ii,ie) * ( h + zenv(ii,ie) )
	    end if
          end do
          mass = mass + sum * area
        end do
 
        if( iespec .gt. 0 ) write(88,*) 'tot mass: ',mass
        res = mass
 
        end
 
c*****************************************************************
