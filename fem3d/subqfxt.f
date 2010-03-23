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

	subroutine qflux3d(it,dt,nkn,nlvdi,temp,dq)

c computes new temperature (forced by heat flux) - 3d version

	implicit none

	integer it
	real dt
	integer nkn
	integer nlvdi
	real temp(nlvdi,1)
	double precision dq	!total energy introduced [(W/m**2)*dt*area = J]

	integer ilhkv(1)
	common /ilhkv/ilhkv
	real evapv(1)
	common /evapv/evapv
c local
	integer k
	integer l,lmax,kspec
	integer mode,level
        integer yes
	real tm,tnew,hm
	real albedo
	real hdecay,adecay,qsbottom,botabs
	real rtot
	real row
        real qs,ta,tb,uw,cc,ur,p,e,r,q
	real qsens,qlat,qlong,evap,qrad
	real area
	double precision ddq
c functions
	real depnode,areanode
	integer ifemopa
c save
	integer n93
	save n93
	data n93 / 0 /

	call qflux_compute(yes)
	if( yes .le. 0 ) return

	!call qfget(qs,ta,tb,uw,cc,ur,p,e,r,q)

	row = 1026.

	mode = +1	!use new time step for depth
	level = 1	!heat transfer only to first layer

c---------------------------------------------------------
c botabs	1. ->	bottom absorbs remaining radiation
c		0. ->	everything is absorbed in last layer
c hdecay	depth of e-folding decay of radiation
c		0. ->	everything is absorbed in first layer
c---------------------------------------------------------

	hdecay = 0.
	botabs = 0.

	adecay = 0.
	if( hdecay .gt. 0. ) adecay = 1. / hdecay

c---------------------------------------------------------
c loop over nodes
c---------------------------------------------------------

        ddq = 0.

	do k=1,nkn

	  hm = depnode(level,k,mode)
	  tm = temp(level,k)
	  area = areanode(level,k)
	  lmax = ilhkv(k)
	  if( hdecay .le. 0. ) lmax = 1

	  call meteo_get_values(k,qs,ta,ur,tb,uw,cc,p)

	   call heatareg (ta,p,uw,ur,cc,tm,qsens,qlat,qlong,evap)
	  !call heatpom  (ta,p,uw,ur,cc,tm,qsens,qlat,qlong,evap)
	  !call heatgill (ta,p,uw,ur,cc,tm,qsens,qlat,qlong,evap)
	  !call heatlucia(ta,p,uw,tb,cc,tm,qsens,qlat,qlong,evap)

          qrad = - ( qlong + qlat + qsens )
	  rtot = qs + qrad
	  call make_albedo(tm,albedo)

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
	    temp(l,k) = tnew
	    albedo = 0.
	    qrad = 0.
	    qs = qsbottom
	  end do

c	  -----------------------
c	  evap is in [kg/(m**2 s)] -> convert it to [m/s]
c	  evap is normally positive -> we are loosing mass
c	  -----------------------

	  evap = evap / row			!evaporation in m/s
	  evapv(k) = evap			!positive if loosing mass

	  !if( k .eq. 100 ) write(6,*) 'qflux3d: ',qrad,qs,dt,hm,tm,tnew

          ddq = ddq + rtot * dt * area

	end do

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

