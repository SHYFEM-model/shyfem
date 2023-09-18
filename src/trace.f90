!
! $Id: subcus.f,v 1.58 2010-03-08 17:46:45 georg Exp $
!
! trace routines
!
! revision log :
!
! 14.03.2005    ggu     subroutine traccia
! 30.06.2005    ggu     traccia changed, new routines jamal, sedimt
! 01.12.2005    ggu     more changes in traccia
! 16.10.2009    ggu     some changes and documentation for traccia
! 01.12.2010	ggu	routines adjusted for Marco (Romanian coast -> uv)
! 15.12.2010	ggu	routines copied from subcus.f
!
! notes :
!
! still to be done:
!
!	call from main, not from subcus
!	general parameters -> traccia.dat, etc...
!	configure from STR
!
!******************************************************************
!
! check following routines:
!
!	get_next_traccia  write_traccia  interpolate_traccia
!
! check following logical variables:
!
! usedts	time is in YYMMDD etc.. and must be converted to it
!			-> change in get_next_traccia, write_traccia
! bnearpoint	if true find closest point and use this value
! bintmiss	if true use last good value
!		use (xinit,yinit) if no good value has been found yet
!		in this case also (xinit,yinit) have to be set
! bnodry	treat dry areas as not found elements
!
! if one of bnearpoint or bintmiss is .true. a value will be always found
! if both are .false. points outside of the domain return -999.
!
! fantina:	bnearpoint = .true.   bintmiss = .true.  usedts=.true.
! rachel:	bnearpoint = .false.  bintmiss = .true.  usedts=.false.
! marco:	bnearpoint = .false.  bintmiss = .true.  usedts=.false.
!
! for carote di fantina, use usedts=.false.
!
! if there are traccie completely out of domain,
!	use bnearpoint = .false.  bintmiss = .true.
!
!*****************************************************************
!------------------------------------------------------------------
        module trace
!------------------------------------------------------------------
        contains
!------------------------------------------------------------------

	subroutine traccia

	use geom_dynamic
	use hydro_admin
        use fil

	implicit none

	include 'param.h'

	include 'femtime.h'


	integer itnew,itold,ierr
	double precision zp,zold,znew
	integer year,date,time
	character*80 line

	integer iuin,iuout,itp,np,iep
	double precision xdp,ydp
	double precision xp,yp
	double precision uv(4)
	save iuin,iuout,itp,np,iep
	save xp,yp
	save xdp,ydp
	save line

	integer icall
	save icall
	data icall /0/

	if( icall .eq. -1 ) return

!-------------------------------------------------------------------
! initialize
!-------------------------------------------------------------------

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

!-------------------------------------------------------------------
! prepare time
!-------------------------------------------------------------------

	itnew = it
	itold = it - idt

	if( itp .lt. itold ) goto 89
	ierr = 0

!-------------------------------------------------------------------
! loop on traccie read
!-------------------------------------------------------------------

	do while( ierr .eq. 0 .and. itold .le. itp .and. itp .le. itnew )

	  call interpolate_traccia(iep,itp,itold,itnew,xp,yp,zp,uv)

	  write(6,*) 'traccia: ',itold,itnew,itp,xp,yp,zp,iep
	  write(99,*) 'traccia: ',itold,itnew,itp,xp,yp,zp,iep,iwegv(abs(iep)),np
	  call write_traccia(iuout,itp,np,xdp,ydp,zp,line,uv)
	  
	  call get_next_traccia(iuin,itp,np,xdp,ydp,line,uv,ierr)
	  xp = xdp
	  yp = ydp

	end do

	if( ierr .gt. 0 ) goto 97
	if( ierr .lt. 0 ) icall = -1

!-------------------------------------------------------------------
! end of routine
!-------------------------------------------------------------------

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

!*****************************************************************

	subroutine interpolate_traccia(iep,itp,itold,itnew,xp,yp,zp,uv)

	use geom_dynamic
	use hydro_baro
	use hydro_admin
        use regular

	implicit none

	integer iep
	integer itp,itold,itnew
	double precision xp,yp
	double precision zp			!interpolated value
	double precision uv(4)		!interpolated value

	double precision xinit,yinit	!lido for rachel
	parameter (xinit=38889.,yinit=32745.2)

	include 'param.h'

	logical bintmiss
	logical bnearpoint
	logical bnodry
	double precision zold,znew
	integer iepold,k
	double precision xold,yold
	save xold,yold,iepold

	integer icall
	save icall
	data icall /0/

! bnearpoint	if true find closest point and use this value
! bintmiss	if true use last good value
!		use (xinit,yinit) if no good value has been found yet
!		in this case also (xinit,yinit) have to be set
! bnodry	treat dry areas as not found elements
!
! if one of bnearpoint or bintmiss is .true. a value will be always found
! if both are .false. points outside of the domain return -999.
!
! fantina:	bnearpoint = .true.   bintmiss = .true.
! rachel:	bnearpoint = .false.  bintmiss = .true.

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

!*****************************************************************

	subroutine intp_z(iep,xp,yp,itp,itold,itnew,zp)

! interpolates water level

	use hydro_admin
        use regular

	implicit none

	integer iep
	double precision xp,yp
	integer itp,itold,itnew
	double precision zp

	include 'param.h'

	double precision zold,znew

	call femintp(iep,zeov(1,iep),xp,yp,zold)
	call femintp(iep,zenv(1,iep),xp,yp,znew)

	zp = zold + (znew-zold)*(itp-itold)/float(itnew-itold)

	end

!*****************************************************************

	subroutine intp_uv(iep,xp,yp,itp,itold,itnew,uv)

! interpolates current velocity

	use hydro_baro
	use hydro_vel
	use levels
        use elems_dealing

	implicit none

	include 'param.h'

	integer iep
	double precision xp,yp
	integer itp,itold,itnew
	double precision uv(2)


	integer level,lmax
	double precision uold,unew,vold,vnew
	double precision ho,hn

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

!*****************************************************************

	subroutine get_next_traccia(iuin,itp,np,xdp,ydp,line,uv,ierr)

! 2004 07 18 10 0 41 621 37002.0159999998 39108.9670000002
! 23964600  37778.7 39962.7      60 SA14 2005-10-05::08:50:00

        use dts

	implicit none

	integer iuin,itp,np
	double precision xdp,ydp
	double precision uv(2)
	
	character*80 line
	integer ierr

	integer year,month,day,hour,min,sec
	integer date,time

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

!	------------------------------------------
!	the next part is only executed if we have to convert
!	the time (given as year month...) to fem_time
!	if we alread read fem_time nothing more is to do
!	------------------------------------------

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

!*****************************************************************

	subroutine write_traccia(iuout,itp,np,xdp,ydp,zp,line,uv)

! writes raw traccia

        use dts

	implicit none

	integer iuout,itp,np
	double precision xdp,ydp
	double precision zp
	double precision uv(4)
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

!*****************************************************************

	function get_nearest_point(xp,yp)

	use basin

	implicit none

	integer get_nearest_point
	double precision xp,yp


	include 'param.h'

	integer knear,k
	double precision dist,d
	double precision dx,dy

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

!*****************************************************************

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
	
!*****************************************************************

	subroutine write_traccia0(iuout,itp,np,xdp,ydp,zp,tz)

! old routine -> do not use anymore

        use dts
	implicit none

	integer iuout,itp,np
	integer time,ihour
	double precision xdp,ydp
	double precision zp
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
!	  (rachel) 2319029.6 5032260.39 38835 27 may 2004
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

!*****************************************************************

!------------------------------------------------------------------
        end module trace
!------------------------------------------------------------------
