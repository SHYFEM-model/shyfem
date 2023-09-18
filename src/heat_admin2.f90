!
! $Id: subqfxt.f,v 1.15 2009-11-18 16:50:37 georg Exp $
!
! heat flux module (file administration)
!
! contents :
!
! subroutine qflux_init
!		initializes routines
! subroutine qflux3d(it,dt,nkn,nlvdi,temp,dq)
!		computes new temperature (forced by heat flux) - 3d version
!
! revision log :
!
! 01.02.2002    ggu     new as administrative routines
! 11.04.2003    ggu     albedo introduced in heat2t
! 16.08.2004    ggu     some comments on usage, heat2t copied to subqfxu4.f
! 22.02.2005    ggu     subroutine qflux2d deleted
! 23.03.2006    ggu     changed time step to real
! 20.08.2007    ggu     prepared for evaporation - salinity feedback
! 17.04.2008    ggu     save evaporation in evapv for later use
! 17.11.2008    ggu     new call to qfcheck_file() before opening
! 10.03.2009    ggu     new qflux_read(), new call to meteo_get_values()
! 11.03.2009    ggu     new routine qflux_compute()
! 27.08.2009    ggu     new call to heatgill, new routine heatareg
! 11.11.2009    ggu     handle abosrption of short rad in more than one layer
! 04.03.2011    ggu     new routine heatgotm
! 23.03.2011    ggu     new routine check_heat() to check for Nan, new iheat
! 25.03.2011    ggu     new parameters iheat,hdecay,botabs
! 29.04.2014    ccf     read qsens, qlat, long from file
! 16.06.2014	ccf	new routine heatcoare, which also update wind stress
! 20.06.2014	ccf	new routine for computing sea surface skin temperature
! 18.09.2015	ccf	do not compute heat fluxes in dry nodes
! 18.09.2015	ccf	checks only for levdbg > 2
! 26.10.2015    ggu     critical omp sections introduced (eliminated data race)
! 07.04.2016    ggu     compute total evaporation
! 04.05.2016    ccf     do not pass albedo into heat2t
! 04.05.2016    ggu     include effect of ice cover
!
! notes :
!
! qflux_init is called in subn11.f
! qflux_read is called in subn11.f
! qflux3d is called in newbcl.f
! qfget is called in qflux3d (here)
!
! to use heat flux module insert in section "name" of STR file:
!
! $name
!     qflux = 'qflux.dat'
! $end
!
! if such a file is given the heat flux module is used and works
!   without any other intervention (ibarcl>0)
!
! other notes in subqfxf.f
!
!*****************************************************************************
!----------------------------------------------------------------------------
        module heat_admin2
!----------------------------------------------------------------------------
        contains
!----------------------------------------------------------------------------

        subroutine compute_heat_flux

! computes heat flux through bulk formulas

        use ts
        use levels, only : nlvdi,nlv
        use basin, only : nkn,nel,ngr,mbw
        use time_util

        implicit none

        include 'param.h'

        include 'femtime.h'

        double precision dq
        double precision dt

        call get_timestep(dt)

        call qflux3d(it,dt,nkn,nlvdi,tempv,dq)  !compute heat flux

        end

!*******************************************************************

	subroutine qflux_compute(yes)

! returnes flag -> compute heat flux or not

        use para

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

!*****************************************************************************

	subroutine qflux_init

! initializes routines

        use para
        use heat_admin

	implicit none

	character*80 file

        call getfnm('qflux',file)

	if( file .ne. ' ' ) then
	  call qfcheck_file(file)
	  call qfinit(file)
	end if

	end

!*****************************************************************************

	subroutine qflux_read(it)

! reads new meteo data

        use heat_admin

	implicit none

	integer it

	call qfmake(it)

	end

!*****************************************************************************

        subroutine qflux3d(it,dt,nkn,nlvddi,temp,dq)

! computes new temperature (forced by heat flux) - 3d version

        use meteo
        use ts
        use levels
        use basin, only : xgv,ygv  
        use hydro_print 
        use para
        use dts
        use elems_dealing
        use meteo_forcing,    only: meteo_get_heat_values,meteo_get_heat_extra,get_pe_values
        use topological
        use heat_temp
        use heat_mfsbulk
        use heat_coare
        use heat_gotm
        use heat_areg,  only: heatpom,heatareg
        use heat_gill
        use heat_default
        use meteo_admin
        use defnames

        implicit none

        include 'subqfxm.h'

        integer it
        double precision dt
        integer nkn
        integer nlvddi
        double precision temp(nlvddi,nkn) 
        double precision dq	!total energy introduced [(W/m**2)*dt*area = J]

! local
        integer levdbg
        logical bdebug
        logical baverevap
        logical bwind
        integer k
        integer l,lmax,kspec
        integer mode
        integer yes
        integer iheat
        integer isolp  
        double precision tm,tnew,hm
        double precision salt,tfreeze
        double precision albedo
        double precision hdecay,adecay,qsbottom,botabs
        double precision qtot
        double precision qs,ta,tb,uw,cc,ur,p,e,r,q
        double precision ddlon,ddlat  
        double precision dp,uuw,vvw  
        double precision qsens,qlat,qlong,evap,qrad
        double precision qswa  
        double precision cice,aice,fice
        double precision ev,eeff
        double precision area
        double precision evaver
        double precision uub,vvb        

        double precision ddq

        double precision, save, allocatable :: dtw(:)	!Warm layer temp. diff
        double precision, save, allocatable :: tws(:)	!Skin temperature (deg C)

        double precision hb			!depth of modelled T in first layer [m]
        double precision usw		!surface friction velocity [m/s]
        double precision qss		!Net shortwave flux 
        double precision cd			!wind drag coefficient

        integer iy,im,id,ih,imn,isec   
        integer days  

        integer itdrag
! save
        integer n93,icall
        save n93,icall
        data n93,icall / 0 , 0 /
        save bdebug,bwind

        call qflux_compute(yes)
        if( yes .le. 0 ) return

!---------------------------------------------------------
! iheat		1=areg  2=pom  3=gill  4=dejak  5=gotm  6=COARE3.0
!               7=read flux from file
! hdecay	depth of e-folding decay of radiation
!		0. ->	everything is absorbed in first layer
! botabs	1. ->	bottom absorbs remaining radiation
!		0. ->	everything is absorbed in last layer
!---------------------------------------------------------

	baverevap = .false.
	aice = 0.	!ice cover for heat: 1: use  0: do not use

	iheat = nint(getpar('iheat'))
        isolp = nint(getpar('isolp'))   
	hdecay = getpar('hdecay')
	botabs = getpar('botabs')

        if( icall .eq. 0 ) then
!$OMP CRITICAL
	  write(6,*) 'qflux3d: iheat,hdecay,botabs: ',iheat,hdecay,botabs  
!$OMP END CRITICAL

          allocate(dtw(nkn))
          allocate(tws(nkn))
          dtw = 0.
          tws(:) = temp(1,:)

          itdrag = nint(getpar('itdrag'))
          bwind = itdrag .eq. 4
          if ( bwind ) then   
            if ( iheat .eq. 6 .or. iheat .eq. 8 ) then   
              write(6,*) 'itdrag = ',itdrag
              write(6,*) 'iheat = ',iheat
            else
              write(6,*) 'Erroneous value for itdrag = ',itdrag
              write(6,*) 'Use itdrag = 4 only with iheat = 6 or 8'
              stop 'error stop qflux3d: itdrag'
            end if
          end if

          levdbg = nint(getpar('levdbg'))
          bdebug = levdbg .ge. 2 

        end if
        icall = 1

!---------------------------------------------------------
! set other parameters
!---------------------------------------------------------

        mode = +1       !use new time step for depth

        adecay = 0.
        if( hdecay .gt. 0. ) adecay = 1. / hdecay

!---------------------------------------------------------
! set date parameters for iheat=8   
!---------------------------------------------------------

        call dts2dt(it,iy,im,id,ih,imn,isec)
        days=id-1+jdmon(iy,im-1)

!---------------------------------------------------------
! loop over nodes
!---------------------------------------------------------

        ddq = 0.

!        if(bmpi) then
!          if(my_id.eq.0) then
!            call rebuild_3d_nodes(temp,outTemp)
!            open(unit=1500, file="temp.txt", action='write')
!            write(1500,*),outTemp(1,:)
!          else
!            call rebuild_3d_nodes(temp)
!          end if
!        else
!
!        end if
!

        do k=1,nkn

          if (is_dry_node(k)) then	!do not compute if node is dry
            dtw(k)   = 0.
            tws(k)   = temp(1,k)
            evapv(k) = 0.
            cycle
          end if

          tm = temp(1,k)
          salt = saltv(1,k)
          call get_ice(k,cice)
          fice = 1. - aice*cice		!fice = 0 if fully ice covered
          tfreeze = -0.0575*salt
          area = areanode(1,k)
          lmax = ilhkv(k)
          if( isolp .eq. 0 .and. hdecay .le. 0. ) lmax = 1   

          call meteo_get_heat_values(k,qs,ta,ur,tb,uw,cc,p)
          call make_albedo(tm,albedo)
          qss = fice * qs * (1. - albedo)

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
            call heatcoare(ta,p,uw,ur,cc,tws(k),r,qss,qsens,qlat,qlong,evap,cd)
            if ( bwind ) windcd(k) = cd
          else if( iheat .eq. 7 ) then
            qsens = ta
            qlat  = ur
            qlong = -cc   !change sign of long wave radiation given by ISAC
            evap  = qlat / (2.5008e6 - 2.3e3 * tm)	!pom, gill, gotm
          else if( iheat .eq. 8 ) then
            ddlon = xgv(k)    
            ddlat = ygv(k)   
            uub = uprv(1,k)  
            vvb = vprv(1,k) 
            call meteo_get_heat_extra(k,dp,uuw,vvw)
            call heatmfsbulk(days,im,ih,ddlon,ddlat,ta,p,uuw,vvw,dp,    &
     &                   cc,tm,uub,vvb,qsens,qlat,qlong,evap,qswa,cd)   
            qss= fice * qswa  !albedo (monthly) already in qshort1 -> qswa  
            if ( bwind ) windcd(k) = cd 
          else
            write(6,*) 'iheat = ',iheat
            stop 'error stop qflux3d: value for iheat not allowed'
          end if

          if (bdebug) call check_heat(k,tm,qsens,qlat,qlong,evap)

          qrad = - fice * ( qlong + qlat + qsens )
          qtot = qss + qrad

          do l=1,lmax
            hm = depnode(l,k,mode)
            tm = temp(l,k)
            if (isolp == 0) then  
              if( hdecay .le. 0. ) then
                qsbottom = 0.
              else
                qsbottom = qss * exp( -adecay*hm )
                if( l .eq. lmax ) qsbottom = botabs * qsbottom
              end if
            else if (isolp == 1) then
              !Jerlov classification (1968) type I - clear water
              qsbottom = qss * (0.58 * exp(-2.8571 * hm ) +  0.42 * exp(-0.0435 * hm ))
              if( l .eq. lmax ) qsbottom = botabs * qsbottom
            else
              write(6,*) 'Erroneous value for isolp = ',isolp
              write(6,*) 'Use isolp = 0 (hdecay) or 1 (Jerlov t-I)'
              stop 'error stop qflux3d: isolp'
            end if
            call heat2t(dt,hm,qss-qsbottom,qrad,tm,tnew)
            if (bdebug) call check_heat2(k,l,qss,qsbottom,qrad,albedo,tm,tnew)
            tnew = max(tnew,tfreeze)
            temp(l,k) = tnew
            albedo = 0.
            qrad = 0.
            qss = qsbottom
          end do

!         ---------------------------------------------------------
!         Compute sea surface skin temperature
!         ---------------------------------------------------------

	  tm   = temp(1,k)
	  hb   = depnode(1,k,mode) * 0.5
          usw  = max(1.e-5, sqrt(sqrt(tauxnv(k)**2 + tauynv(k)**2)))
          qrad =  -(qlong + qlat + qsens)
	  call tw_skin(qss,qrad,tm,hb,usw,dt,dtw(k),tws(k))

!	  ---------------------------------------------------------
!	  evap is in [kg/(m**2 s)] -> convert it to [m/s]
!	  evap is normally positive -> we are loosing mass
!	  ---------------------------------------------------------

	  evap = evap / rhow			!evaporation in m/s
	  evapv(k) = evap			!positive if loosing mass

          ddq = ddq + qtot * dt * area

	end do

	dq = ddq

!---------------------------------------------------------
! compute total evaporation
!---------------------------------------------------------

	if( baverevap ) then
	  !call aver_nodal(evapv,evaver)	!in evaver is average of evaporation m/s
	  write(678,*) it,evaver
	end if

!---------------------------------------------------------
! special output
!---------------------------------------------------------

	k = min(nkn,1000)
	k = -1
	if( k .gt. 0 ) then
	  if( n93 .eq. 0 ) then
	    n93 = ifemopa('opening file 93','.93','form','new')
	  end if
	  write(n93,*) 'qflux3d: ',it,temp(1,k)
	end if

!---------------------------------------------------------
! end of routine
!---------------------------------------------------------

	end

!*****************************************************************************

	subroutine check_heat(k,tm,qsens,qlat,qlong,evap)

        use chk_NaN

	implicit none

	integer k
	double precision tm,qsens,qlat,qlong,evap

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

!*****************************************************************************

	subroutine check_heat2(k,l,qs,qsbottom,qrad,albedo,tm,tnew)

        use chk_NaN

	implicit none

	integer k,l
	double precision qs,qsbottom,qrad,albedo,tm,tnew

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

!*****************************************************************************

!----------------------------------------------------------------------------
        end module heat_admin2
!----------------------------------------------------------------------------
