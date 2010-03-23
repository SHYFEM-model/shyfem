c
c $Id: subqfxa.f,v 1.6 2003/06/20 15:42:22 georg Exp $
c
c heat flux module (administration)
c
c contents :
c
c subroutine qfluxr(it)
c			initializes heat flux module and reads from file
c subroutine qflux2d(it,idt)
c			computes new temperature (forced by heat flux)
c subroutine qfnext(it,qs,ta,tb,uw,cc,ier)
c			reads next record of flux file
c
c revision log :
c
c 01.06.1998	ggu&lz	written from scratch (nearly)
c 24.06.1998	ggu&lz	subroutines from lucia integrated
c 24.05.2000	ggu	no details in element treatment, 3D algorithm
c 28.11.2001	ggu	less debug output
c
c comments :
c
c qs	net solar radiation (reflection already subtracted)
c ta	air temperature
c tb	wet bulb temperature
c uw	wind speed
c cc	cloud cover (0 clear sky, 1 totally covered)
c
c***********************************************************

	subroutine qfluxr(it)

c initializes heat flux module and reads from file

	implicit none

	integer it

	integer iunit,itfold,itfnew
	real qsold,taold,tbold,uwold,ccold
	real qsnew,tanew,tbnew,uwnew,ccnew
	real qs,ta,tb,uw,cc
	common /qflxi/ iunit,itfold,itfnew
	common /qflxro/ qsold,taold,tbold,uwold,ccold
	common /qflxrn/ qsnew,tanew,tbnew,uwnew,ccnew
	common /qflxra/ qs,ta,tb,uw,cc
	save /qflxi/, /qflxro/, /qflxrn/, /qflxra/

	character*80 file
	logical bdebug,bwrite
	integer itact,ier
	integer ifileo
	real r

	integer icall, itperiod, itmany
	save    icall, itperiod, itmany
	data    icall, itperiod, itmany / 0 , 31536000 , 0 /

c itperiod - periodic radiation conditions
c
c the file must contain values for time 0 and itperiod
c these records should be equal
c
c year : 31536000      day : 86400

	if( icall .lt. 0 ) return

	bdebug = .true.
	bdebug = .false.
	bwrite = .true.

c initialize subroutine

	if( icall .eq. 0 ) then
	  call getfnm('qflux',file)

	  if( file .eq. " " ) then
	    icall = -1
	    iunit = -1
	    return
	  end if

	  iunit = ifileo(0,file,'form','old')
	  if( iunit .le. 0 ) then
	    write(6,'(a,a)') 'filename: ',file(1:65)
	    stop 'error stop qfluxin: Cannot open file'
	  end if

	  write(6,*) 'heat flux file opened : '
	  write(6,*) file

	  call qfnext(itfold,qsold,taold,tbold,uwold,ccold,ier)
	  if( ier .ne. 0 ) goto 99

	  call qfnext(itfnew,qsnew,tanew,tbnew,uwnew,ccnew,ier)
	  if( ier .ne. 0 ) goto 99

	end if

c normal call of subroutine -> read from file

	icall = icall + 1

c handle periodic conditions

	itact = it - itmany * itperiod

	if( itact .gt. itperiod ) then
	  rewind(iunit)
	  itmany = itmany + 1
	  itact  = itact  - itperiod
	  itfold = itfold - itperiod
	  itfnew = itfnew - itperiod
	  write(6,*) 'period in qflux: ',it,itact,itfold,itfnew
	end if

c iterate until new record has time greater than actual time

	do while( itact .gt. itfnew )

	  itfold = itfnew
	  qsold = qsnew
	  taold = tanew
	  tbold = tbnew
	  uwold = uwnew
	  ccold = ccnew

	  call qfnext(itfnew,qsnew,tanew,tbnew,uwnew,ccnew,ier)
	  if( ier .ne. 0 ) goto 99
	  if( bdebug .or. bwrite ) then
	    write(6,'(a,i10,5f12.2)') 'qfnext: '
     +			,itfnew,qsnew,tanew,tbnew,uwnew,ccnew
	  end if

	end do

c do interpolation

	r = (itact-itfold) / float(itfnew-itfold)

	qs = r * qsnew + (1.-r) * qsold
	ta = r * tanew + (1.-r) * taold
	tb = r * tbnew + (1.-r) * tbold
	uw = r * uwnew + (1.-r) * uwold
	cc = r * ccnew + (1.-r) * ccold

	if( bdebug ) then
	  write(6,'(a,i10,5f12.2)') 'qfluxr: ',itact,qs,ta,tb,uw,cc
	end if

	return
   99	continue
	write(6,*) it,itact,ier
	stop 'error stop qfluxr: error reading file'
	end

c***********************************************************

	subroutine qflux2d(it,idt)

c computes new temperature (forced by heat flux) - 2d version

	implicit none

	integer it,idt

c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real toev(1)
        common /toev/toev
        real tnv(1)
        common /tnv/tnv
c common flux module
	integer iunit,itfold,itfnew
	common /qflxi/ iunit,itfold,itfnew
	real qs,ta,tb,uw,cc
	common /qflxra/ qs,ta,tb,uw,cc
c local
	integer ie,k,ii
	real tm,tnew,hm,aj,vol
	real rtot
c functions
	real depele,scalele
	
	if( iunit .le. 0 ) return

c loop over elements

	do ie=1,nel

	  hm = depele(ie,+1)		!average depth in element
	  tm = scalele(ie,toev)		!average temperature in element
	  call subtem(idt,hm,qs,uw,cc,ta,tb,tm,tnew,rtot) !new temperature
	  call addsvele(ie,toev,tnew-tm)	!temp to element vertices

	end do

c distribute new temperature to nodes

	call ele2node(toev,tnv)

	k = 3
	k = -1
        if( k .gt. 0 ) write(93,*) 'qflux2d: ',it,tnv(k)

c	write(6,*) 'qfluxc: ',tnv(15),tnv(197)
c	write(93,*) it/3600,tnv(15),tnv(197)

	end

c***********************************************************

	subroutine qflux3d(it,idt,nkn,levmax,temp,dq)

c computes new temperature (forced by heat flux) - 3d version

	implicit none

	integer it,idt
	integer nkn
	integer levmax
	real temp(levmax,1)
	double precision dq	!total energy introduced [(W/m**2)*dt*area = J]

c common flux module
	integer iunit,itfold,itfnew
	common /qflxi/ iunit,itfold,itfnew
	real qs,ta,tb,uw,cc
	common /qflxra/ qs,ta,tb,uw,cc
	integer ilhkv(1)
	common /ilhkv/ilhkv
c local
	integer k
	integer l,lmax,kspec
	integer mode,level
	real tm,tnew,hm
	real rtot
	real cw,row
	double precision heatconold,heatconnew,ddq
c functions
	real depnode,areanode
	
	if( iunit .le. 0 ) return

	mode = +1	!use new time step for depth
	level = 1	!heat transfer only to first layer

	cw  = 3991.
	row = 1026.

c---------------------------------------------------------
c loop over nodes
c---------------------------------------------------------

        heatconold = 0.
        heatconnew = 0.
        ddq = 0.

	do k=1,nkn

	  hm = depnode(level,k,mode)
	  tm = temp(level,k)

	  heatconold = heatconold + tm * areanode(level,k) * cw * row

	  call subtem(idt,hm,qs,uw,cc,ta,tb,tm,tnew,rtot) !new temperature
	  temp(level,k) = tnew

          ddq = ddq + rtot * idt * areanode(level,k)
          heatconnew = heatconnew + tnew * areanode(level,k) * cw * row

	end do

	dq = ddq

c---------------------------------------------------------
c special output
c---------------------------------------------------------

	kspec = 0
	if( kspec .gt. 0 ) then
	write(6,*) 'ATTENTION: in qflux3d homogenizing temperatures...'
	do k=1,nkn
	  lmax = ilhkv(k)
	  do l=1,lmax
	    temp(l,k) = temp(l,kspec)
	  end do
	end do
	end if

	k = 3
	k = -1
	k = 1000
        if( k .gt. 0 ) write(93,*) 'qflux3d_old: ',it,temp(1,k)

c	write(6,*) 'qfluxc: ',tnv(15),tnv(197)
c	write(93,*) it/3600,tnv(15),tnv(197)

c---------------------------------------------------------
c end of routine
c---------------------------------------------------------

	end

c***********************************************************

	subroutine qfnext(it,qs,ta,tb,uw,cc,ier)

c reads next record of flux file
c
c it	time [s]
c qs	solar radiation [W/m**2]
c ta	air temperature [C]
c tb	wet bulb temperature [C]
c ur	relative humidity [%], ([0-100])
c uw	wind speed [m/s]
c cc	cloud cover [0-1], 0 no clouds
c
c format is one of the two following:
c
c it qs ta tb uw cc
c it qs ta ur uw cc

	implicit none

	integer it
	real qs,ta,tb,uw,cc
	real ur
	integer ier

	integer iunit,itfold,itfnew
	common /qflxi/ iunit,itfold,itfnew

	logical bdebug
	integer ios

	integer itold
	save itold
	data itold / 0 /

	bdebug = .false.

c-------------------------------------------------------
c read next record
c-------------------------------------------------------

c	read(iunit,*,iostat=ios) it,qs,ta,tb,uw,cc
	read(iunit,*,iostat=ios) it,qs,ta,ur,uw,cc

c-------------------------------------------------------
c convert relative humidity to wet bulb temperature
c-------------------------------------------------------

	call rh2twb(ta,ur,tb) 

c-------------------------------------------------------
c debug and error handling
c-------------------------------------------------------

	if( bdebug ) then
	  write(6,'(a,i10,5f12.2)') 'qfnext: ',it,qs,ta,tb,uw,cc
	end if

	if( ios .ne. 0 ) then
	  if( ios .gt. 0 ) then
	    write(6,*) 'error reading line close to it = ',itold
	  else if( ios .lt. 0 ) then
	    write(6,*) 'end of file found close to it = ',itold
	  end if
	  stop 'error stop qfnext: no more data or read error'
	end if

	itold = it
	ier = ios

c-------------------------------------------------------
c end of routine
c-------------------------------------------------------

	end

c*****************************************************************************
c*****************************************************************************
c*****************************************************************************

