c
c $Id: subqfxt.f,v 1.15 2009-11-18 16:50:37 georg Exp $
c
c heat flux module (file administration)
c
c contents :
c
c subroutine qflux_init
c		initializes routines
c subroutine qflux3d(it,dt,nkn,nlvdi,temp,dq)
c		computes new temperature (forced by heat flux) - 3d version
c
c revision log :
c
c 01.02.2002    ggu     new as administrative routines
c 11.04.2003    ggu     albedo introduced in heat2t
c 16.08.2004    ggu     some comments on usage, heat2t copied to subqfxu4.f
c 22.02.2005    ggu     subroutine qflux2d deleted
c 23.03.2006    ggu     changed time step to real
c 20.08.2007    ggu     prepared for evaporation - salinity feedback
c 17.04.2008    ggu     save evaporation in evapv for later use
c 17.11.2008    ggu     new call to qfcheck_file() before opening
c 10.03.2009    ggu     new qflux_read(), new call to meteo_get_values()
c 11.03.2009    ggu     new routine qflux_compute()
c 27.08.2009    ggu     new call to heatgill, new routine heatareg
c 11.11.2009    ggu     handle abosrption of short rad in more than one layer
c 04.03.2011    ggu     new routine heatgotm
c 23.03.2011    ggu     new routine check_heat() to check for Nan, new iheat
c 25.03.2011    ggu     new parameters iheat,hdecay,botabs
c 29.04.2014    ccf     read qsens, qlat, long from file
c 16.06.2014	ccf	new routine heatcoare, which also update wind stress
c 20.06.2014	ccf	new routine for computing sea surface skin temperature
c
c notes :
c
c qflux_init is called in subn11.f
c qflux_read is called in subn11.f
c qflux3d is called in newbcl.f
c qfget is called in qflux3d (here)
c
c to use heat flux module insert in section "name" of STR file:
c
c $name
c     qflux = 'qflux.dat'
c $end
c
c if such a file is given the heat flux module is used and works
c   without any other intervention (ibarcl>0)
c
c other notes in subqfxf.f
c
c*****************************************************************************

	subroutine qflux_compute(yes)

c returnes flag -> compute heat flux or not

	implicit none

	integer yes

	character*80 file

        call getfnm('qflux',file)

	if( file .eq. ' ' ) then
	  yes = 0
	else
	  yes = 1
	end if

	end

c*****************************************************************************

	subroutine qflux_init

c initializes routines

	implicit none

	character*80 file

        call getfnm('qflux',file)

	if( file .ne. ' ' ) then
	  call qfcheck_file(file)
	  call qfinit(file)
	end if

	end

c*****************************************************************************

	subroutine qflux_read(it)

c reads new meteo data

	implicit none

	integer it

	call qfmake(it)

	end

c*****************************************************************************

	subroutine qflux3d(it,dt,nkn,nlvddi,temp,dq)

c computes new temperature (forced by heat flux) - 3d version

	use mod_meteo
	use mod_ts
	use levels

	implicit none

        include 'param.h'
        include 'subqfxm.h'

	integer it
	real dt
	integer nkn
	integer nlvddi
	real temp(nlvddi,1)
	double precision dq	!total energy introduced [(W/m**2)*dt*area = J]


c local
	logical bdebug
	integer k
	integer l,lmax,kspec
	integer mode
        integer yes
	integer iheat
	real tm,tnew,hm
	real salt,tfreeze
	real albedo
	real hdecay,adecay,qsbottom,botabs
	real rtot
        real qs,ta,tb,uw,cc,ur,p,e,r,q
	real qsens,qlat,qlong,evap,qrad
	real ev,eeff
	real area

	double precision ddq

	real, save, allocatable :: dtw(:)	!Warm layer temp. diff
	real, save, allocatable :: tws(:)	!Skin temperature (deg C)

	real hb			!depth of modelled T in first layer [m]
	real usw		!surface friction velocity [m/s]
	real rad		!Net shortwave flux 
	real cd			!wind drag coefficient

	integer itdrag
c functions
	real depnode,areanode,getpar
	integer ifemopa
	logical bwind
	save bwind
c save
	integer n93,icall
	save n93,icall
	data n93,icall / 0 , 0 /

	call qflux_compute(yes)
	if( yes .le. 0 ) return

c---------------------------------------------------------
c iheat		1=areg  2=pom  3=gill  4=dejak  5=gotm  6=COARE3.0
c               7=read flux from file
c hdecay	depth of e-folding decay of radiation
c		0. ->	everything is absorbed in first layer
c botabs	1. ->	bottom absorbs remaining radiation
c		0. ->	everything is absorbed in last layer
c---------------------------------------------------------

	!iheat = 1
	!hdecay = 0.
	!botabs = 0.
	iheat = nint(getpar('iheat'))
	hdecay = getpar('hdecay')
	botabs = getpar('botabs')

	if( icall .eq. 0 ) then
	  write(6,*) 'qflux3d: iheat,hdecay,botabs: '
     +				,iheat,hdecay,botabs  

	  allocate(dtw(nkn))
	  allocate(tws(nkn))
	  dtw = 0.
	  tws(:) = temp(1,:)

          itdrag = nint(getpar('itdrag'))
	  bwind = itdrag .eq. 4
	  if ( bwind .and. iheat .ne. 6 ) then
            write(6,*) 'Erroneous value for itdrag = ',itdrag
            write(6,*) 'Use itdrag = 4 only with iheat = 6'
            stop 'error stop qflux3d: itdrag'
	  end if
	end if
	icall = 1

c---------------------------------------------------------
c set other parameters
c---------------------------------------------------------

	mode = +1	!use new time step for depth

	adecay = 0.
	if( hdecay .gt. 0. ) adecay = 1. / hdecay

c---------------------------------------------------------
c loop over nodes
c---------------------------------------------------------

	bdebug = .false.
        ddq = 0.

	do k=1,nkn

	  tm = temp(1,k)
	  salt = saltv(1,k)
	  tfreeze = -0.0575*salt
	  area = areanode(1,k)
	  lmax = ilhkv(k)
	  if( hdecay .le. 0. ) lmax = 1

          call meteo_get_heat_values(k,qs,ta,ur,tb,uw,cc,p)
	  call make_albedo(tm,albedo)
	  rad = qs * (1. - albedo)

	  if( iheat .eq. 1 ) then
	    call heatareg (ta,p,uw,ur,cc,tm,qsens,qlat,qlong,evap)
	  else if( iheat .eq. 2 ) then
	    call heatpom  (ta,p,uw,ur,cc,tm,qsens,qlat,qlong,evap)
	  else if( iheat .eq. 3 ) then
	    call heatgill (ta,p,uw,ur,cc,tm,qsens,qlat,qlong,evap)
	  else if( iheat .eq. 4 ) then
	    call heatlucia(ta,p,uw,tb,cc,tm,qsens,qlat,qlong,evap)
	  else if( iheat .eq. 5 ) then
	    call heatgotm (ta,p,uw,ur,cc,tm,qsens,qlat,qlong,evap)
	  else if( iheat .eq. 6 ) then
	    call get_pe_values(k,r,ev,eeff)
	    call heatcoare(ta,p,uw,ur,cc,tws(k),r,rad,qsens,qlat,
     +                     qlong,evap,cd)
	    if ( bwind ) windcd(k) = cd
	  else if( iheat .eq. 7 ) then
	    qsens = ta
	    qlat  = ur
	    qlong = -cc		!change sign of long wave radiation given by ISAC
	    evap  = qlat / (2.5008e6 - 2.3e3 * tm)	!pom, gill, gotm
	  else
	    write(6,*) 'iheat = ',iheat
	    stop 'error stop qflux3d: value for iheat not allowed'
	  end if

	  bdebug = k .eq. 0 .and. it .ge. 378989200
	  if( bdebug ) then
	    write(6,*) iheat,k
	    write(6,*) ta,p,uw
	    write(6,*) tb,ur,cc,tm
	    write(6,*) qsens,qlat,qlong,evap
	  end if

	  call check_heat(k,tm,qsens,qlat,qlong,evap)

          qrad = - ( qlong + qlat + qsens )
	  rtot = qs + qrad

	  do l=1,lmax
	    hm = depnode(l,k,mode)
	    tm = temp(l,k)
	    if( hdecay .le. 0. ) then
		qsbottom = 0.
	    else
		qsbottom = qs * exp( -adecay*hm )
		if( l .eq. lmax ) qsbottom = botabs * qsbottom
	    end if
            call heat2t(dt,hm,qs-qsbottom,qrad,albedo,tm,tnew)
	    !call check_heat2(k,l,qs,qsbottom,qrad,albedo,tm,tnew)
	    tnew = max(tnew,tfreeze)
	    temp(l,k) = tnew
	    albedo = 0.
	    qrad = 0.
	    qs = qsbottom
	  end do

c         ---------------------------------------------------------
c         Compute sea surface skin temperature
c         ---------------------------------------------------------

	  tm   = temp(1,k)
	  hb   = depnode(1,k,mode) * 0.5
          usw  = max(1.e-5, sqrt(sqrt(tauxnv(k)**2 + tauynv(k)**2)))
          qrad =  -(qlong + qlat + qsens)
	  call tw_skin(rad,qrad,tm,hb,usw,dt,dtw(k),tws(k))

c	  ---------------------------------------------------------
c	  evap is in [kg/(m**2 s)] -> convert it to [m/s]
c	  evap is normally positive -> we are loosing mass
c	  ---------------------------------------------------------

	  evap = evap / rhow			!evaporation in m/s
	  evapv(k) = evap			!positive if loosing mass

	  !if( k .eq. 100 ) write(6,*) 'qflux3d: ',qrad,qs,dt,hm,tm,tnew

          ddq = ddq + rtot * dt * area

	end do

	!call check2Dr(nlvdi,nlvdi,nkn,temp,-10.,+50.,'qflux','tempv')

	dq = ddq

c---------------------------------------------------------
c special output
c---------------------------------------------------------

	k = min(nkn,1000)
	k = -1
	if( k .gt. 0 ) then
	  if( n93 .eq. 0 ) then
	    n93 = ifemopa('opening file 93','.93','form','new')
	  end if
	  write(n93,*) 'qflux3d: ',it,temp(1,k)
	end if

c---------------------------------------------------------
c end of routine
c---------------------------------------------------------

	end

c*****************************************************************************

	subroutine check_heat(k,tm,qsens,qlat,qlong,evap)

	implicit none

	integer k
	real tm,qsens,qlat,qlong,evap

	logical is_r_nan

	if( is_r_nan(tm) ) goto 98
	if( is_r_nan(qsens) ) goto 98
	if( is_r_nan(qlat) ) goto 98
	if( is_r_nan(qlong) ) goto 98
	if( is_r_nan(evap) ) goto 98

	return
   98	continue
	write(6,*) k,tm,qsens,qlat,qlong,evap
	stop 'error stop check_heat: Nan found'
	end

c*****************************************************************************

	subroutine check_heat2(k,l,qs,qsbottom,qrad,albedo,tm,tnew)

	implicit none

	integer k,l
	real qs,qsbottom,qrad,albedo,tm,tnew

	logical is_r_nan

	if( is_r_nan(qs) ) goto 98
	if( is_r_nan(qsbottom) ) goto 98
	if( is_r_nan(qrad) ) goto 98
	if( is_r_nan(albedo) ) goto 98
	if( is_r_nan(tm) ) goto 98
	if( is_r_nan(tnew) ) goto 98

	if( abs(tm) .gt. 100. ) goto 97
	if( abs(tnew) .gt. 100. ) goto 97

	return
   97	continue
	write(6,*) k,l,tm,tnew
	stop 'error stop check_heat2: t out of range'
   98	continue
	write(6,*) k,l,qs,qsbottom,qrad,albedo,tm,tnew
	stop 'error stop check_heat2: Nan found'
	end

c*****************************************************************************

