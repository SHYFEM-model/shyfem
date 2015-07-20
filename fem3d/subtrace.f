c
c $Id: subcus.f,v 1.58 2010-03-08 17:46:45 georg Exp $
c
c trace routines
c
c revision log :
c
c 14.03.2005    ggu     subroutine traccia
c 30.06.2005    ggu     traccia changed, new routines jamal, sedimt
c 01.12.2005    ggu     more changes in traccia
c 16.10.2009    ggu     some changes and documentation for traccia
c 01.12.2010	ggu	routines adjusted for Marco (Romanian coast -> uv)
c 15.12.2010	ggu	routines copied from subcus.f
c
c notes :
c
c still to be done:
c
c	call from main, not from subcus
c	general parameters -> traccia.dat, etc...
c	configure from STR
c
c******************************************************************
c
c check following routines:
c
c	get_next_traccia  write_traccia  interpolate_traccia
c
c check following logical variables:
c
c usedts	time is in YYMMDD etc.. and must be converted to it
c			-> change in get_next_traccia, write_traccia
c bnearpoint	if true find closest point and use this value
c bintmiss	if true use last good value
c		use (xinit,yinit) if no good value has been found yet
c		in this case also (xinit,yinit) have to be set
c bnodry	treat dry areas as not found elements
c
c if one of bnearpoint or bintmiss is .true. a value will be always found
c if both are .false. points outside of the domain return -999.
c
c fantina:	bnearpoint = .true.   bintmiss = .true.  usedts=.true.
c rachel:	bnearpoint = .false.  bintmiss = .true.  usedts=.false.
c marco:	bnearpoint = .false.  bintmiss = .true.  usedts=.false.
c
c for carote di fantina, use usedts=.false.
c
c if there are traccie completely out of domain,
c	use bnearpoint = .false.  bintmiss = .true.
c
c*****************************************************************

	subroutine traccia

	use mod_geom_dynamic
	use mod_hydro

	implicit none

	include 'param.h'

	include 'femtime.h'


	integer itnew,itold,ierr
	real zp,zold,znew
	integer ifileo
	integer year,date,time
	character*80 line

	integer iuin,iuout,itp,np,iep
	double precision xdp,ydp
	real xp,yp
	real uv(4)
	save iuin,iuout,itp,np,iep
	save xp,yp
	save xdp,ydp
	save line

	integer icall
	save icall
	data icall /0/

	if( icall .eq. -1 ) return

c-------------------------------------------------------------------
c initialize
c-------------------------------------------------------------------

	if( icall .eq. 0 ) then
	  iuin = ifileo(55,'traccia.dat','form','old')
	  if( iuin .le. 0 ) goto 99
	  iuout = ifileo(55,'traccia_out.dat','form','new')
	  if( iuout .le. 0 ) goto 99

	  line = ' '
	  write(6,*) 'initializing traccia...'

	  call get_next_traccia(iuin,itp,np,xdp,ydp,line,uv,ierr)
	  xp = xdp
	  yp = ydp
	  if( ierr .ne. 0 ) goto 98

	  iep = 0
	  icall = 1
	end if

c-------------------------------------------------------------------
c prepare time
c-------------------------------------------------------------------

	itnew = it
	itold = it - idt

	if( itp .lt. itold ) goto 89
	ierr = 0

c-------------------------------------------------------------------
c loop on traccie read
c-------------------------------------------------------------------

	do while( ierr .eq. 0 .and. itold .le. itp .and. itp .le. itnew )

	  call interpolate_traccia(iep,itp,itold,itnew,xp,yp,zp,uv)

	  write(6,*) 'traccia: ',itold,itnew,itp,xp,yp,zp,iep
	  write(99,*) 'traccia: ',itold,itnew,itp,xp,yp,zp,iep
     +					,iwegv(abs(iep)),np
	  call write_traccia(iuout,itp,np,xdp,ydp,zp,line,uv)
	  
	  call get_next_traccia(iuin,itp,np,xdp,ydp,line,uv,ierr)
	  xp = xdp
	  yp = ydp

	end do

	if( ierr .gt. 0 ) goto 97
	if( ierr .lt. 0 ) icall = -1

c-------------------------------------------------------------------
c end of routine
c-------------------------------------------------------------------

	return
   99	continue
	stop 'error stop traccia: cannot open file'
   98	continue
	stop 'error stop traccia: first record of traccia'
   97	continue
	write(6,*) itp,np,xp,yp,ierr
	stop 'error stop traccia: read error in traccia'
   89	continue
	write(6,*) itp,itold,itnew
	stop 'error stop traccia: simulation start too late'
        end

c*****************************************************************

	subroutine interpolate_traccia(iep,itp,itold,itnew,xp,yp,zp,uv)

	use mod_geom_dynamic
	use mod_hydro_baro
	use mod_hydro

	implicit none

	integer iep
	integer itp,itold,itnew
	real xp,yp
	real zp			!interpolated value
	real uv(4)		!interpolated value

	real xinit,yinit	!lido for rachel
	parameter (xinit=38889.,yinit=32745.2)

	include 'param.h'

	logical bintmiss
	logical bnearpoint
	logical bnodry
	real zold,znew
	integer iepold,k
	real xold,yold
	save xold,yold,iepold

	integer get_nearest_point

	integer icall
	save icall
	data icall /0/

c bnearpoint	if true find closest point and use this value
c bintmiss	if true use last good value
c		use (xinit,yinit) if no good value has been found yet
c		in this case also (xinit,yinit) have to be set
c bnodry	treat dry areas as not found elements
c
c if one of bnearpoint or bintmiss is .true. a value will be always found
c if both are .false. points outside of the domain return -999.
c
c fantina:	bnearpoint = .true.   bintmiss = .true.
c rachel:	bnearpoint = .false.  bintmiss = .true.

	bintmiss = .true.
	bnearpoint = .false.
	bnodry = .true.

	if( icall .eq. 0 ) then		!get first element
	  xold = xinit
	  yold = yinit
	  icall = 1
	  call find_elem_from_old(0,xold,yold,iepold)
	end if

	call find_elem_from_old(iep,xp,yp,iep)

	if( bnodry .and. iep .gt. 0 .and. iwegv(iep) .gt. 0 ) iep = -iep

	if( iep .gt. 0 ) then		!element found -> interpolate
	    call intp_z(iep,xp,yp,itp,itold,itnew,zp)
	    call intp_uv(iep,xp,yp,itp,itold,itnew,uv(3))
	else if( bnearpoint ) then	!interpolate using nearest point
	    k = get_nearest_point(xp,yp)
	    zold = zov(k)
	    znew = znv(k)
	    zp = zold + (znew-zold)*(itp-itold)/float(itnew-itold)
	    uv(3) = 0.
	    uv(4) = 0.
	    write(6,*) 'using nearest point: ',k,zp
	else if( bintmiss ) then	!interpolate using old point
	    iep = iepold
	    if( iep .le. 0 ) goto 99
	    xp = xold
	    yp = yold
	    call intp_z(iep,xp,yp,itp,itold,itnew,zp)
	    call intp_uv(iep,xp,yp,itp,itold,itnew,uv(3))
	else
	    zp = -999.
	    uv(3) = 0.
	    uv(4) = 0.
	end if

	iepold = iep
	xold = xp
	yold = yp

	return
   99	continue
	write(6,*) iepold
	write(6,*) xold,yold,xp,yp
	write(6,*) itp,itold,itnew
	stop 'error stop interpolate_traccia: no element found'
	end

c*****************************************************************

	subroutine intp_z(iep,xp,yp,itp,itold,itnew,zp)

c interpolates water level

	use mod_hydro

	implicit none

	integer iep
	real xp,yp
	integer itp,itold,itnew
	real zp

	include 'param.h'

	real zold,znew

	call femintp(iep,zeov(1,iep),xp,yp,zold)
	call femintp(iep,zenv(1,iep),xp,yp,znew)

	zp = zold + (znew-zold)*(itp-itold)/float(itnew-itold)

	end

c*****************************************************************

	subroutine intp_uv(iep,xp,yp,itp,itold,itnew,uv)

c interpolates current velocity

	use mod_hydro_baro
	use mod_hydro_vel
	use levels

	implicit none

	include 'param.h'

	integer iep
	real xp,yp
	integer itp,itold,itnew
	real uv(2)


	integer level,lmax
	real uold,unew,vold,vnew
	real ho,hn

	real depele

	level = 1		!0 for barotropic

	lmax = ilhv(iep)

	if( level .le. 0 ) then
	  ho = depele(iep,-1)
	  hn = depele(iep,+1)
	  uold = uov(iep)/ho
	  unew = unv(iep)/hn
	  vold = vov(iep)/ho
	  vnew = vnv(iep)/hn
	else if( level .gt. lmax ) then
	  uold = 0.
	  unew = 0.
	  vold = 0.
	  vnew = 0.
	else
	  uold = ulov(level,iep)
	  unew = ulnv(level,iep)
	  vold = vlov(level,iep)
	  vnew = vlnv(level,iep)
	end if

	uv(1) = uold + (unew-uold)*(itp-itold)/float(itnew-itold)
	uv(2) = vold + (vnew-vold)*(itp-itold)/float(itnew-itold)

	end

c*****************************************************************

	subroutine get_next_traccia(iuin,itp,np,xdp,ydp,line,uv,ierr)

c 2004 07 18 10 0 41 621 37002.0159999998 39108.9670000002
c 23964600  37778.7 39962.7      60 SA14 2005-10-05::08:50:00

	implicit none

	integer iuin,itp,np
	double precision xdp,ydp
	real uv(2)
	
	character*80 line
	integer ierr

	integer year,month,day,hour,min,sec
	integer date,time
	logical dts_initialized

	logical usedts
	!parameter (usedts=.true.)
	parameter (usedts=.false.)

	integer icall
	save icall
	data icall / 0 /

	if( usedts ) then
	  read(iuin,*,iostat=ierr) year,month,day,hour,min,sec,np,xdp,ydp
	else
	  !read(iuin,*,iostat=ierr) itp,xdp,ydp,np,line
	  !write(6,*) 'new traccia read: ',itp,line
	  read(iuin,*,iostat=ierr) itp,xdp,ydp,np,uv(1),uv(2)	!bajo
	  write(6,*) 'new traccia read: ',itp,np
	  line = ' '
	end if

	if( ierr .ne. 0 ) return
	if( .not. usedts ) return

c	------------------------------------------
c	the next part is only executed if we have to convert
c	the time (given as year month...) to fem_time
c	if we alread read fem_time nothing more is to do
c	------------------------------------------

	if( icall .eq. 0 ) then
	  date = 10000*year + 101
	  time = 0
	  call dtsini(date,time)
	  write(6,*) 'get_next_traccia: initialized dts ',date,time
	  icall = 1
	end if

	call dts2it(itp,year,month,day,hour,min,sec)
	write(6,*) 'new traccia read: ',itp,np

	end

c*****************************************************************

	subroutine write_traccia(iuout,itp,np,xdp,ydp,zp,line,uv)

c writes raw traccia

	implicit none

	integer iuout,itp,np
	double precision xdp,ydp
	real zp
	real uv(4)
	character*80 line

	integer year,month,day,hour,min,sec
	integer ilen

	logical usedts
	parameter (usedts=.true.)

	if( usedts ) then
	  call dts2dt(itp,year,month,day,hour,min,sec)
	  !write(iuout,*) year,month,day,hour,min,sec,np,xdp,ydp,zp
	  write(iuout,1000) year,month,day,hour,min,sec,np,xdp,ydp,zp
	else
	  call get_string_length(line,ilen)
	  !write(iuout,*) itp,xdp,ydp,np,zp,' ',line(1:ilen)
	  !write(iuout,2000) itp,xdp,ydp,np,zp,' ',line(1:ilen)
	  write(iuout,2001) itp,xdp,ydp,np,uv		!bajo
	end if

	return
 1000	format(i5,5i3,i10,3f14.6)
 2000	format(i10,2f14.6,i10,f14.6,a,a)
 2001	format(i10,2f14.6,i10,4f14.6)
	end

c*****************************************************************

	function get_nearest_point(xp,yp)

	use basin

	implicit none

	integer get_nearest_point
	real xp,yp


	include 'param.h'

	integer knear,k
	real dist,d
	real dx,dy

	dist=1.e+30
	knear = 0

	do k=1,nkn
	  dx = xgv(k) - xp
	  dy = ygv(k) - yp
	  d = sqrt( dx*dx + dy*dy )
	  if( d .lt. dist ) then
	    dist = d
	    knear = k
	  end if
	end do

	get_nearest_point = knear
	
	end

c*****************************************************************

	subroutine get_string_length(line,ilen)

	implicit none

	character*80 line
	integer ilen

	integer imax,i

	imax = len(line)
	do i=imax,1,-1
	  if( line(i:i) .ne. ' ' ) goto 1
	end do
    1	continue
	ilen = i
	if( ilen .le. 0 ) ilen = 1

	end
	
c*****************************************************************

	subroutine write_traccia0(iuout,itp,np,xdp,ydp,zp,tz)

c old routine -> do not use anymore

	implicit none

	integer iuout,itp,np
	integer time,ihour
	double precision xdp,ydp
	real zp
	integer tz

	double precision x0,y0
	parameter(x0=2330000. - 50000., y0=5000000.)
	integer itype
	parameter(itype=1)	!0=raw  1=fantina  2=rachel

	integer year,month,day,hour,min,sec
	double precision x,y,z

	call dts2dt(itp,year,month,day,hour,min,sec)
	ihour = hour - tz	!GMT to MET

	if( itype .eq. 0 ) then
	  x = xdp
	  y = ydp
	  z = zp
	  write(iuout,*) year,month,day,hour,min,sec,np,xdp,ydp,itp,zp
	else if( itype .eq. 1 ) then
	  x = xdp + x0
	  y = ydp + y0
	  z = 0.23 - zp
	  time = ihour*10000+min*100+sec
	  write(iuout,'(i8,i8,f14.3,2x,f14.3,2x,f10.3)') time,np,x,y,z
	else if( itype .eq. 2 ) then
c	  (rachel) 2319029.6 5032260.39 38835 27 may 2004
	  x = xdp + x0
	  y = ydp + y0
	  z = zp - 0.23
	  time = ihour*3600+min*60+sec
	  write(iuout,*) x,y,time,day,month,year,z
	else 
	  stop 'error stop write_traccia: invalid itype'
	end if

	write(89,'(6f13.3)') xdp,ydp,x0,y0,x,y

	end

c*****************************************************************

