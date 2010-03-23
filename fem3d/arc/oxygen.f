c
c $Id: oxygen.f,v 1.5 2003/03/25 14:08:55 georg Exp $
c
c oxygen module for tecnomare
c
c revision log :
c
c 20.01.1999	ggu	adapted from reactor Cosimo
c 26.01.1999	ggu	new write to file with 3D NOS routines
c 17.02.1999	ggu	tested and changed for special areas
c 18.02.1999	ggu	array knoxy -> read it and treat it in check routine
c 24.02.1999	ggu	different parameters, subroutine outoxy
c 06.04.1999	ggu	do not print anything if ioxy = 0
c
c notes :
c
c------------------------------------------------------
c
c files changed for integration of oxygen module ->
c
c	hp.f		oxytar
c	subcst.f	inoxy, ckoxy
c	subnhs.f	rdoxy,proxy,tsoxy
c
c------------------------------------------------------
c
c solve this ->
c
c	*initialization
c	*read
c	*oxy -> oxye
c	 dry?
c	 saturation of own model
c
c------------------------------------------------------
c
c input file format ->
c
c oxydef, cldef, oxyb, ioxy		in oxypar
c
c k,itype,cloro				in oxyarr
c
c------------------------------------------------------
c
c****************************************************

        function oxysat(temp,salt)

c computes oxygen level at saturation
c
c in mg/L , to have mmol/L divide by 32

        implicit none

        real oxysat
        real temp,salt

        oxysat = (14.6244-0.367134*temp+4.4972E-3*temp**2
     :     -0.0966*salt+0.00205*salt*temp+2.739E-4*salt**2)

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

c***********************************************

	subroutine oxytar(it,idt)

c eco-model cosimo

	implicit none

	include 'param.h'

	integer it	!time in seconds
	integer idt	!time step in seconds

	real oxy(nkndim), oxye(3,neldim), rov(nkndim)
	real clalph(nkndim)
	integer iotype(nkndim)
	common /clalph/clalph, /iotype/iotype

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real eps1,eps2,pi,flag,high,higi
        common /mkonst/ eps1,eps2,pi,flag,high,higi

	real soev(3,1), snv(1)
        common /soev/soev, /snv/snv
	real toev(3,1), tnv(1)
        common /toev/toev, /tnv/tnv

	real hkv(1)
	common /hkv/hkv
        real up0v(1), vp0v(1)
        common /up0v/up0v, /vp0v/vp0v
	real v1v(1)
	common /v1v/v1v
	integer nen3v(3,1)
	common /nen3v/nen3v

	integer ioxy
	integer ie,ii,k
	real t,s,dt,u,v
	real umod
	integer iunit
	real secnds,cs,ri,rid
	real oxydef
	integer isact

	real oxysat
	real getpar
	integer iround

	integer istot
	save istot
	real rkpar,azpar,oxyb
	save rkpar,azpar,oxyb
	integer iu,itmcon,idtcon
	save iu,itmcon,idtcon

	integer icall
	save icall
	data icall /0/

c initialization

	if( icall .le. -1 ) return

	if( icall .eq. 0 ) then
	  ioxy = iround(getpar('ioxy'))
	  if( ioxy .le. 0 ) icall = -1
	  if( icall .le. -1 ) return

	  icall = 1

c	  initialize oxygen array

	  oxydef = getpar('oxydef')
	  if( oxydef .lt. 0. ) then	!use saturation value
	    do k=1,nkn
	      oxy(k) = oxysat(tnv(k),snv(k))
	    end do
	  else				!use given value
	    do k=1,nkn
	      oxy(k) = oxydef
	    end do
	  end if

	  call ktoe(oxy,oxye)

	  do k=1,nkn
	    rov(k) = flag
	  end do

c	  other parameters

          oxyb=getpar('oxyb')
          rkpar=getpar('chpar')
          call getaz(azpar)
          istot=iround(getpar('istot'))         !$$istot

c	  initialize file

          iu = 55
          itmcon = iround(getpar('itmcon'))
          idtcon = iround(getpar('idtcon'))
          call confop(iu,itmcon,idtcon,1,'oxy')

	  write(6,*) 'oxygen model initialized...'

	end if

c normal call

	call oxybnd(rov,oxyb)

c advection and diffusion

          do isact=1,istot                !$$istot
            call conz(oxy,oxye,v1v,dt,rkpar,azpar,istot,isact)  !$$azpar
	    call conzbc(oxy,oxye,v1v,rov,flag,azpar)
          end do

c reactions

	dt = idt
	secnds = mod(it,86400)
	cs = 2. * 3.14159 / 86400.
	ri = -cos( cs * secnds )
	rid = cs * sin( cs * secnds )	!derivative

c	k = 100
c	write(6,*) 'rid...: ',ri,rid	!ggu-write
c	write(6,*) 'oxysat: ',tnv(k),snv(k),oxysat(tnv(k),snv(k))

	do k=1,nkn		!loop on nodes
	    u = up0v(k)
	    v = vp0v(k)
	    umod = sqrt( u*u + v*v )
	    call oxyrea(it,dt,tnv(k),snv(k),oxy(k),ri,rid
     +				,iotype(k),hkv(k)
     +				,clalph(k),umod)
	end do

	call ktoe(oxy,oxye)	!put back into elements

c write of results

	call confil(iu,itmcon,idtcon,15,1,oxy)

	end

c*************************************************************

	subroutine oxybnd(rov,value)

c sets boundary conditions for oxygen module

	implicit none

	real rov(1)
	real value

	integer ibc,nbc,ibtyp
	integer nbnds,itybnd

	nbc = nbnds()

	do ibc=1,nbc
	  ibtyp = itybnd(ibc)
	  if( ibtyp .eq. 1 .or. ibtyp .eq. 2 ) then
	    call setbnd(ibc,value,rov)
	  end if
	end do

	end

c*************************************************************

	subroutine oxyrea(it,dt,t,s,oxy,ri,rid,iotype,hdepth
     +				,clalph,umod)

c reactor for oxygen

	implicit none

c arguments
	integer it
	real dt
	real t,s
	real oxy
	real ri,rid
	integer iotype
	real hdepth
	real clalph
	real umod
c parameters ----------------------------------------- change here - GGU
	real fact
	parameter(fact=1./86400.)
	real ogenx,ogenn
	real oalgx,oalgn
	real ooxydf,ooxydm
	real osodf,osodm
	real orespf,orespm
	real k0,epsp
	real krear
	real kts0,qsed,ts0
	real codaux
	real ktd0,qdet,td0
	real torbid
	parameter(ogenx=10.,ogenn=5.)
	parameter(oalgx=0.,oalgn=0.)
	parameter(orespf=4.*fact,orespm=4.*fact)
	parameter(ooxydf=2.*fact,ooxydm=0.)
c	parameter(orespf=1.*fact,orespm=1.*fact)
c	parameter(ooxydf=0.5*fact,ooxydm=0.)
	parameter(osodf=0.,osodm=0.)
	parameter(k0=1.,epsp=0.)
c	parameter(krear=0.5*fact)
	parameter(krear=3.0*fact)
	parameter(kts0=3.*fact,qsed=1.08,ts0=20.)
	parameter(codaux=0.)
	parameter(ktd0=3.*fact,qdet=1.08,td0=20.)
	parameter(torbid=1.)
c -----------------------------------------------------
c local
	real ogen
	real oresp,ooxyd,osod
	real alghe,ke,aux
	real oalg,otot
	real osat,orear
	real ores,osed,odet
	real ksed,kdet
	real d0
c functions
	real oxysat
c internal functions
	real diffwt,tw
	diffwt(tw) = 1.0e-4 * ( 1.2e-05 + 4.58e-07 * tw )

c general formulation

	if( iotype .eq. 0 ) then		!general
	  ogen = 0.5 * (ogenx-ogenn) * rid
	  oxy = oxy + dt * ogen
	  return
	else if( iotype .eq. 1 ) then		!fish
	  oresp = orespf
	  ooxyd = ooxydf
	  osod = osodf
	else if( iotype .eq. 2 ) then		!mussels
	  oresp = orespm
	  ooxyd = ooxydm
	  osod = osodm
	else
	  stop 'error stop oxyrea : internal error iotype (1)'
	end if

c alghe

	alghe = clalph
	ke = k0 + epsp * clalph
	aux = ke * hdepth
	oalg = 0.5 * (oalgx-oalgn) * rid
	oalg = oalg * alghe * (1. - exp(-aux)) / aux
	otot = otot + oalg

c re-areazione

	d0 = diffwt(t)
	osat = oxysat(t,s)
	aux = sqrt(d0*umod/hdepth) / hdepth
	aux = aux + krear
	if( aux * dt .gt. 1. ) write(6,*) 'oxygen forcing too high...'
	orear = aux * (osat-oxy)
	otot = otot + orear

c respirazione

	ores = oresp + ooxyd * 0.5 * ( 1. + ri )
	otot = otot - ores
	
c sediments

	ksed = kts0 * qsed**(t-ts0)
	osed = osod + 2.67 * ksed * codaux
	otot = otot - osed
	
c detritus

	kdet = ktd0 * qdet**(t-td0)
	odet = kdet * torbid
	otot = otot - odet

c advance time step

	oxy = oxy + dt * otot

	call outoxy(iotype,oxy,oalg,orear,ores,osed,odet)
c	write(6,*) 'special: ',oxy,otot,iotype	!ggu-write

	end

c*************************************************************

	subroutine outoxy(is,oxy,oalg,orear,ores,osed,odet)

c writes special output

	implicit none

	integer is
	real oxy
	real oalg,orear,ores,osed,odet
	real o1,o2,o3,o4,o5
	real otot

	otot = abs(oalg) + abs(orear) + abs(ores) + abs(osed) + abs(odet)

	o1 = oalg * 100. / otot
	o2 = orear * 100. / otot
	o3 = ores * 100. / otot
	o4 = osed * 100. / otot
	o5 = odet * 100. / otot

	write(6,'(a,i5,f6.2,e12.4,5f6.1)') 'outoxy: ',is
     +			,oxy,otot,o1,o2,o3,o4,o5

	end

c*************************************************************

	subroutine ktoe(rk,re)

c node values to element values

	implicit none

	real rk(1)
	real re(3,1)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nen3v(3,1)
	common /nen3v/nen3v

	integer ie,ii,k

        do ie=1,nel
          do ii=1,3
            k = nen3v(ii,ie)
            re(ii,ie) = rk(k)
          end do
        end do

	end
	
c***********************************************
c*********************************************** user interface
c***********************************************

	subroutine inoxy

	implicit none

	include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	real clalph(nkndim)
	integer iotype(nkndim)
	integer knoxy(nkndim)
	common /clalph/clalph, /iotype/iotype, /knoxy/knoxy

	integer k
	real oxydef,cldef,oxyb

	oxydef = -1.		!this means we use saturation value later
	cldef = 0.
	oxyb = 0.

	do k=1,nkn
	  knoxy(k) = 0		!only used to preserve node num. until check
	  iotype(k) = 0
	  clalph(k) = -1.	!to recognize from given values
	end do

	call sctpar('oxypar')
	call addpar('oxydef',oxydef)
	call addpar('cldef',cldef)
	call addpar('oxyb',oxyb)
	call addpar('ioxy',0.)

	end

c***********************************************

	subroutine rdoxy

c reads section dealing with node description
c
c $oxyarr
c node	type	[clorophyll]
c ...
c $end
c
c where type means
c
c	0	general
c	1	fish
c	2	mussels

	implicit none

	include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	real clalph(nkndim)
	integer iotype(nkndim)
	integer knoxy(nkndim)
	common /clalph/clalph, /iotype/iotype, /knoxy/knoxy
	save /clalph/, /iotype/, /knoxy/

        character*80 line
        integer ianz,iar,i,j,k
        real f(20)

        integer nrdlin,iscan,iround,ipint

	k = 0

        do while( nrdlin(line) .ne. 0 )
                ianz = iscan(line,1,f)
                if( ianz .gt. 0 ) then
                        if( ianz .gt. 3 .or. ianz .lt. 2 ) goto 86
			k = k + 1
                        knoxy(k) = iround(f(1))
                        if(knoxy(k).le.0) goto 88
			iotype(k) = iround(f(2))
			if( ianz .eq. 3 ) clalph(k) = f(3)
                else if( ianz .lt. 0 ) then
                        goto 98
                end if
        end do

	return
   86   continue
	write(6,*) 'Only 2 or 3 numbers possible on line : '
        write(6,*) line
        stop 'error stop : rdoxy'
   88   continue
        write(6,*) 'No node like this : ',knoxy(k)
        write(6,*) line
        stop 'error stop : rdoxy'
   98   continue
        write(6,*) 'Read error in line :'
        write(6,*) line
        stop 'error stop : rdoxy'
	end

c***********************************************

	subroutine ckoxy

	implicit none

	include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	real clalph(nkndim)
	integer iotype(nkndim)
	integer knoxy(nkndim)
	common /clalph/clalph, /iotype/iotype, /knoxy/knoxy
	real v1v(1), v2v(1), v3v(1)
	common /v1v/v1v, /v2v/v2v, /v3v/v3v

	integer k,iot
	integer kmax,kext,kint
	real cldef

	integer ipext,ipint
	real getpar

	cldef = getpar('cldef')

c copy data read to aux arrays

	k = 0
	do while( knoxy(k+1) .gt. 0 )
	   k = k + 1
	   v1v(k) = knoxy(k)
	   v2v(k) = iotype(k)
	   v3v(k) = clalph(k)
	   knoxy(k) = 0
	   iotype(k) = 0
	   clalph(k) = -1.
	end do
	kmax = k

c now populate arrays at definite place

	do k=1,kmax
	   kext = nint(v1v(k))
	   kint = ipint(kext)
           if(kint.le.0) goto 88
	   iotype(kint) = nint(v2v(k))
	   clalph(kint) = v3v(k)
	end do

c some default values

	do k=1,nkn
	  if( clalph(k) .lt. 0. ) then
	    clalph(k) = cldef
	  end if
	  if( iotype(k) .ne. 0 ) then
	    iot = iotype(k)
	    if( iot .lt. 0 .or. iot .gt. 2 ) then
	      write(6,*) 'Not allowed oxygen type ',iot
     +				,' of node ',ipext(k)
	      stop 'error stop ckoxy : iotype'
	    end if
	  end if
	end do

	return
   88   continue
        write(6,*) 'No node like this : ',kext
        stop 'error stop : ckoxy'
	end

c***********************************************

	subroutine proxy

	implicit none

	include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	real clalph(nkndim)
	integer iotype(nkndim)
	integer knoxy(nkndim)
	common /clalph/clalph, /iotype/iotype, /knoxy/knoxy

	integer k
	integer ioxy
	real oxydef,cldef,oxyb

	integer ipext
	real getpar

	ioxy = nint(getpar('ioxy'))
	oxydef = getpar('oxydef')
	cldef = getpar('cldef')
	oxyb = getpar('oxyb')

	if( ioxy .eq. 0 ) return

	write(6,*) 'Oxygen model : '
	write(6,*) 'oxydef,cldef,oxyb : ',oxydef,cldef,oxyb

	do k=1,nkn
	  if( iotype(k) .ne. 0 ) then
	    write(6,*) 'node,type : ',ipext(k),iotype(k)
	  end if
	end do

	end

c***********************************************

	subroutine tsoxy

	implicit none

	call proxy

	end

c***********************************************

