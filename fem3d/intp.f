c
c $Id: intp.f,v 1.9 2008-11-03 10:42:25 georg Exp $
c
c closure routines
c
c revision log :
c
c 02.09.2004    ggu     new classification routine
c 04.09.2004    ggu     franco added (zfranco)
c 11.10.2008	ggu	franco passed in
c 16.11.2012	ggu	if franco == 99 -> use PS as forecast (init_ps)
c
c**********************************************************************

c	program intp
c	call intp_test
c	end

c**********************************************************************

	subroutine intp_test

	implicit none

	integer ndim
	parameter (ndim=200)

        real rnodes
        parameter( rnodes = 1852. / 3600. )     ! n. miles / hour

	integer it,j,n,iday
	integer iddt,itanf,itend
	integer itcl,iclose,iclass,ih
	integer itotal,ittot
	real wmax,rmax
	real t,value,dz,zlagoon,zlevel
	real zdate,zmax
	real zextra,zsalv,zrise,zfranco
	real psv,dsv,rainv,windv
	real psav(1),dsav(1)
	real vals(ndim)
	character*4 class

	real ps(ndim)
	real ds(ndim)
	real rain(ndim)
	real wind(ndim)

	save itcl,iclose,iclass

	integer icm
	icm(t) = nint(100.*t)

	call exffils('ps.dat',ndim,ps)
	call exffils('ds.dat',ndim,ds)

        zfranco = 0.05    !allow for extra security
        zfranco = 0.10    !allow for extra security
        zfranco = 0.
	zextra = 0.10
	zextra = 0.00
	zrise = 0.00
	zrise = 0.30
	zrise = 0.50
	zsalv = 0.40 + zextra
	zsalv = 1.00 + zextra
	iclose = 0
	itcl = 0

	iddt = 3600
	iddt = 300
	itanf = 3600
	itend = 150000
	itend = 6000000
	itend = 86400 * 366 - 3600
	zlagoon = 0.
	itotal = 0
	ittot = 0

	do it=itanf,itend,iddt
	  t = it
	  iday = 1 + it/86400
	  call exfintp(ps,t,psav)
	  call exfintp(ds,t,dsav)
	  psv = psav(1)
	  dsv = dsav(1)
	  psv = psv + zrise
	  dsv = dsv + zrise

	  call get_prev(ndim,it,zfranco,n,vals)
	  call get_wind_rain(it,windv,rainv,wmax,rmax,dz)

          call class12new(it,n,vals,zsalv,zextra,zrise,wmax,rmax,psv
     +                  ,zmax,zdate,ih,iclass,class)

	  if( iclose .eq. 1 ) then	!already closed
	    dz = dz*iddt/(1000.*3600.)	!convert to [m/timestep]
	    zlagoon = zlagoon + dz
	    !call raise_zeta(dz)
	  end if
	  !call set_zdate(0,zdate)

          zlevel = 0
          if( ih .gt. 0 ) zlevel = zrise+vals(ih)
	  write(66,2300) it,iclass,class,iclose,ih
     +		,icm(zmax),icm(zlevel),icm(psv),icm(dsv),icm(zlagoon)
	  if( iclass .gt. 0 ) then
	    !write(66,2200) (nint(100.*vals(j)),j=1,20)
	  end if

	  if( iclose .eq. 0 .and. psv .ge. zdate ) then
	    !write(66,2500) 'closing inlets: ',class,it,psv,zdate
	    write(6,2500)  'closing inlets: ',class,it,psv,zdate
	    itcl = it
	    iclose = 1
	    zlagoon = psv
	  else if( iclose .eq. 1 .and. dsv .lt. zlagoon ) then
!     +				.and. psv .lt. zdate ) then
	    !write(66,2500) 'opening inlets: ',class,it,zlagoon,zdate,it-itcl
	    write(6,2500) 'opening inlets: '
     +				,class,it,zlagoon,zdate,it-itcl
	    iclose = 0
	    zlagoon = 0.
	    itotal = itotal + 1
	    ittot = ittot + it-itcl
	  end if

	  !write(6,2000) it,iday,psv,dsv,rainv,windv,(vals(j),j=1,3)
	end do

	write(6,*) itotal,ittot,100.*ittot/(365*86400.)

 2200	    format(20i4)
 2100	  format(i10,i2,a5,i2,i3,4f6.2)
 2300	  format(i10,i2,a5,i2,i3,3x,10i4)
 2000	format(i10,i4,7f9.2)
 2500	format(a,a4,i10,2f9.2,i10)
	end

c********************************************************

	subroutine get_prev(nvdim,it,zfranco,if,vals)

c gets previsioni

	implicit none

	integer nvdim		!dimension of vals
	integer it		!time for forecast
        real zfranco		!allow for extra security
	integer if		!filling of vals on return
	real vals(1)		!values on return [m]

	integer ndim,nfore
	parameter(ndim=30000,nfore=150)

	integer n
	integer itime(ndim)
	integer ifill(ndim)
	real pp(nfore,ndim)
	real paux(nfore)

	integer itanf,itaux,i,j
	logical bps,bdebug
	real zfranc

	save n,itime,ifill,pp
	data n / 0 /

	bdebug = .true.
	zfranc = zfranco
	bps = zfranc .gt. .90		!use PS

	if( n .eq. 0 ) then
	  if( bps ) then
	    call init_ps(ndim,nfore,n,itime,ifill,pp,paux)
	  else
	    call init_prev(ndim,nfore,n,itime,ifill,pp,paux)
	  end if
	  if( bdebug ) then
	    write(6,*) 'prev initialized... ', bps
	    do i=1,20
	      if = ifill(i)
	      write(6,*) itime(i),if
	      write(6,*) (pp(j,i),j=1,if)
	    end do
	  end if
	end if
	
	if( bps ) zfranc = 0.

	itanf = itime(1)
	itaux = 3600*(it/3600)
	i = 1 + (it-itanf)/3600
	if( i .lt. 1 .or. i .gt. n ) then
	  write(6,*) i,n
	  write(6,*) it,itanf,itaux
	  stop 'error stop get_prev: (1)'
	end if
	if( itime(i) .ne. itaux ) then
	  write(6,*) it,itanf,itaux,itime(i),i
	  stop 'error stop get_prev: (2)'
	end if

	if = ifill(i)
	if( if .gt. nvdim) then
	  write(6,*) i,n
	  write(6,*) it,itanf,itaux
	  write(6,*) if,nvdim
	  stop 'error stop get_prev: (3)'
	end if

	do j=1,if
	  vals(j) = zfranc + pp(j,i)
	end do

	end

c********************************************************

	subroutine init_ps(ndim,nfore,n,itime,ifill,pp,paux)

c initializes previsioni

	implicit none

	integer ndim,nfore
	integer n
	integer itime(ndim)
	integer ifill(ndim)
	real pp(nfore,ndim)
	real paux(nfore)

	integer iunit
	integer j,ip,i
	integer it,i1,i2,i3,i4,i5,if
	real ps
	integer ifileo
	character*80 file

        iunit = 55
	file = 'psref.dat'
        iunit = ifileo(iunit,file,'form','old')
	if( iunit .le. 0 ) stop 'error stop init_ps'

	if = 10
	n = 0
    1	continue
	  read(iunit,*,end=2) it,ps

	  n = n + 1
	  if( if .gt. nfore ) stop 'error stop init_prev: nfore'
	  if( n .gt. ndim ) stop 'error stop init_prev: ndim'

	  itime(n) = it
	  ifill(n) = if
	  do j=1,if
	    pp(j,n) = ps
	  end do

	  goto 1
    2	continue

	write(6,*) 'init_prev: records read ',n

	close(iunit)

	do i=1,n
	  do j=1,if
	    ip = min(n,i+j)
	    pp(j,i) = pp(j,ip)
	  end do
	end do

	end

c********************************************************

	subroutine init_prev(ndim,nfore,n,itime,ifill,pp,paux)

c initializes previsioni

	implicit none

	integer ndim,nfore
	integer n
	integer itime(ndim)
	integer ifill(ndim)
	real pp(nfore,ndim)
	real paux(nfore)

	integer iunit
	integer j
	integer it,i1,i2,i3,i4,i5,if
	integer ifileo
	character*80 file

        iunit = 55
	file = 'prev.dat'
        iunit = ifileo(iunit,file,'form','old')
	if( iunit .le. 0 ) stop 'error stop init_prev'

c it	fem time (secs)
c i1	day
c i2	month
c i3	year
c i4	julian day
c i5	hour
c if	fill

	n = 0
    1	continue
	  read(iunit,*,end=2) it,i1,i2,i3,i4,i5,if,(paux(j),j=1,if)

	  n = n + 1
	  if( if .gt. nfore ) stop 'error stop init_prev: nfore'
	  if( n .gt. ndim ) stop 'error stop init_prev: ndim'

	  itime(n) = it
	  ifill(n) = if
	  do j=1,if
	    pp(j,n) = 0.01 * paux(j)
	  end do

	  goto 1
    2	continue

	write(6,*) 'init_prev: records read ',n

	close(iunit)

	end

c********************************************************

	subroutine class12(n,vals,zsalv,zextra,zrise,wmax,rmax,psv
     +			,zmax,zdate,ih,iclass,class)

	implicit none

	integer n		!total number of values in vals
	real vals(1)		!values of forecast
	real zsalv		!level of salvaguardia
	real zextra		!extra level for salvaguardia
	real zrise		!water level rise
	real wmax		!maximum wind speed
	real rmax		!maximum rain
	real psv		!actual value at PS
	real zmax		!maximum level in 24 hours (out)
	real zdate		!closing level (out)
	integer ih		!hours of maximum from now (out)
	integer iclass		!class (out)
	character*(*) class	!class (out)

	integer nmax,i
	real vmax,v
	real zd

	nmax = min(n,24)
	vmax = zrise + vals(1)
	ih = 0

	do i=1,nmax
	  v = zrise + vals(i)
	  vmax = max(vmax,v)
	  if( ih .eq. 0 ) then
		if( v .ge. zsalv ) ih = i
	  end if
	end do

	if( psv .ge. zsalv ) ih = 1	!actual value is higher

	iclass = 0
	class = '0'
        zdate = 9.99

	if( vmax - zrise .ge. 1.50 ) then
	    iclass = 2
	    class = '2'
	    zd = 0.55
            if( ih .le. 6 ) zdate = zd + zextra
            if( ih .le. 12 ) zdate = zd + zextra
	else
	    iclass = 1
            if( rmax .ge. 1. ) then   !case B
                if( wmax .ge. 15. ) then
                  zd = 0.65
                  class = '1BV'
                else
                  zd = 0.80
                  class = '1B '
                end if
            else                      !case A
                if( wmax .ge. 15. ) then
                  zd = 0.70
                  class = '1AV'
                else
                  zd = 0.90
                  class = '1A '
                end if
            end if
            if( ih .le. 4 ) zdate = zd + zextra
            if( ih .le. 12 ) zdate = zd + zextra
	end if

	zmax = vmax

	end

c********************************************************

	subroutine class12new(it,n,vals,zsalv,zextra,zrise,wmax,rmax,psv
     +			,zmax,zdate,ih,iclass,class)

	implicit none

        integer it              !time [sec]
	integer n		!total number of values in vals
	real vals(1)		!values of forecast
	real zsalv		!level of salvaguardia
	real zextra		!extra level for salvaguardia
	real zrise		!water level rise
	real wmax		!maximum wind speed
	real rmax		!maximum rain
	real psv		!actual value at PS
	real zmax		!maximum level in 4 or 6 hours (out)
	real zdate		!closing level (out)
	integer ih		!hours of maximum from now (out)
	integer iclass		!class (out)
	character*(*) class	!class (out)

	integer nmax,nmax4,nmax6,i
        integer ih4,ih6
        integer level
        integer ith
	real vmax,v
        real vmax4,vmax6
        real vmax4p,vmax6p
        real vmax4c,vmax6c
	real zd

        integer iclass_mem
        real zdate_mem
        character*4 class_mem
        save iclass_mem,zdate_mem,class_mem
        data iclass_mem,zdate_mem,class_mem /0,999.,'0'/

        integer levels(4)
        integer iclassv(0:5)
        real zdv(0:5)
        character*3 classv(0:5)
        save levels,iclassv,zdv,classv
        data levels  / 0,0,0,0 /
        data zdv     / 9.00 , 0.90 , 0.80 , 0.70 , 0.65 , 0.55 /
        data iclassv /    0 ,    1 ,    1 ,    1 ,    1 ,    2 /
        data  classv / '0  ', '1A ', '1B ', '1AV', '1BV', '2  '/

        ith = it / 3600
        if( mod(it,3600) .eq. 0 ) then
          do i=1,3
            levels(i) = levels(i+1)
          end do
          levels(4) = 0
        end if

	nmax6 = min(n,6)
	nmax4 = min(n,4)
	vmax4 = zrise + vals(1)
	vmax6 = zrise + vals(1)
	ih4 = 0
	ih6 = 0

	do i=1,nmax4
	  v = zrise + vals(i)
	  vmax4 = max(vmax4,v)
	  if( ih4 .eq. 0 ) then
		if( v .ge. zsalv ) ih4 = i
	  end if
	end do

	do i=1,nmax6
	  v = zrise + vals(i)
	  vmax6 = max(vmax6,v)
	  if( ih6 .eq. 0 ) then
		if( v .ge. zsalv ) ih6 = i
	  end if
	end do

        vmax6p = vals(nmax6) + zrise
        vmax4p = vals(nmax4) + zrise

        vmax6c = vmax6
        vmax4c = vmax4

        vmax6c = vmax6p
        vmax4c = vmax4p
        ih6 = nmax6
        ih4 = nmax4

        !if( psv .ge. zsalv ) then
        !  vmax4 = max(vmax4,psv)
        !  vmax6 = max(vmax6,psv)
        !end if

	if( vmax6c - zrise .ge. 1.50 ) then
            level = 5
        else if( vmax4c .ge. zsalv ) then
            if( rmax .ge. 1. ) then   !case B
                if( wmax .ge. 15. ) then
                  level = 4
                else
                  level = 2
                end if
            else                      !case A
                if( wmax .ge. 15. ) then
                  level = 3
                else
                  level = 1
                end if
            end if
        else
            level = 0
	end if

        do i=1,4
          if( level .gt. levels(i) ) levels(i) = level
        end do

        level = levels(1)

        zd = zdv(level)
        zdate = zd + zextra
        iclass = iclassv(level)
        class = classv(level)

        if( iclass .eq. 0 ) then
            ih = 0
            vmax = vmax4c
        else if( iclass .eq. 1 ) then
            ih = ih4
            vmax = vmax4c
        else if( iclass .eq. 2 ) then
            ih = ih6
            vmax = vmax6c
        else
            stop 'error stop: internal error'
        end if

	zmax = vmax

	end

c********************************************************

	subroutine get_max(it,ih,wmax,rmax)

c gets max of wind and rain for ih hours ahead

	implicit none

	integer it,ih
	real wmax,rmax

	integer ndim
	parameter (ndim=30000)

	integer iwind,irain,nwind,nrain
	integer itmax,i,n,nvar
	integer itwind(ndim)
	integer itrain(ndim)
	real wind(ndim)
	real rain(ndim)

	save iwind,irain
	save itwind,itrain,wind,rain
	save nwind,nrain
	data nwind,nrain /0,0/

	if( nwind .eq. 0 ) then
	  nvar = 1
	  call read_ts('wind.dat',ndim,nvar,nwind,itwind,wind)
	  call read_ts('rain.dat',ndim,nvar,nrain,itrain,rain)
	  iwind = 1
	  irain = 1
	end if

	itmax = it + 3600*ih

c wind

	do while( iwind .le. nwind .and. itwind(iwind) .lt. it )
	  iwind = iwind + 1
	end do
	i = iwind
	wmax = wind(i)
	if( itwind(i) .gt. it .and. i .gt. 1 ) wmax = max(wmax,wind(i-1))
	do while( i .le. nwind .and. itwind(i) .le. itmax )
	  wmax = max(wmax,wind(i))
	  i = i + 1
	end do

c rain

	do while( irain .le. nrain .and. itrain(irain) .lt. it )
	  irain = irain + 1
	end do
	i = irain
	rmax = rain(i)
	if( itrain(i) .gt. it .and. i .gt. 1 ) rmax = max(rmax,rain(i-1))
	do while( i .le. nrain .and. itrain(i) .le. itmax )
	  rmax = max(rmax,rain(i))
	  i = i + 1
	end do
	!write(6,*) 'rain ',it,itmax,itrain(irain),irain,i,rmax

	end

c********************************************************

        subroutine leakage(wind,dz)

        implicit none

        real wind       ![m/s]
        real dz         ![mm/h]

        real w1,w2,dz1,dz2,dw,dzz
        parameter( w1=5. , dz1=2.7 , w2=25. , dz2=21. )
        parameter( dw=w2-w1 , dzz=dz2-dz1 )

        if( wind .le. w1 ) then
          dz = dz1
        else if( wind .ge. w2 ) then
          dz = dz2
        else
          dz = dz1 + dzz * ((wind-w1)/dw)**5
        end if

        end

c********************************************************

	subroutine raise_zeta(dz)

c raises water level of lagoon

	use mod_hydro
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	real dz		!water level rise per time step [m]

	include 'param.h'

        integer ie,ii,k

        do k=1,nkn
          znv(k) = znv(k) + dz
        end do

        do ie=1,nel
          do ii=1,3
            zenv(ii,ie) = zenv(ii,ie) + dz
          end do
        end do

	end

c********************************************************

	subroutine get_wind_rain(it,w,r,wmax,rmax,dz)

c wind in [m/s]
c rain in [mm/h]
c leakage dz in [mm/h]

	implicit none

	integer it
	real w,r
	real wmax,rmax
	real dz

        real rnodes
        parameter( rnodes = 1852. / 3600. )     ! n. miles / hour

	integer ndim
	parameter (ndim=200)

	real rain(ndim)
	real wind(ndim)
	save rain,wind

	real t

	integer icall
	save icall
	data icall / 0 /

	if( icall .eq. 0 ) then
	  icall = 1
	  call exffils('rain.dat',ndim,rain)
	  call exffils('wind.dat',ndim,wind)
	end if

	t = it
	call exfintp(rain,t,r)
	call exfintp(wind,t,w)

	call get_max(it,4,wmax,rmax)

	w = w * rnodes			!convert from nodes to m/s
	wmax = wmax * rnodes
	r = r / 24.			!convert from mm/day to mm/hour
	rmax = rmax / 24.

	dz = 0.
        call leakage(w,dz)
	dz = dz + 2. * r		!also drainage basin [mm/h]

	end

c********************************************************

