c
c $Id: supint.f,v 1.9 2009-09-14 08:31:18 georg Exp $
c
c interactive routines for plotsim
c
c revision log :
c
c 12.02.1999  ggu     adapted to auto mode
c 29.01.2002  ggu     new routine getisec()
c 17.03.2004  ggu     new routine okvar()
c 02.03.2005  ggu     new routines set_flag and get_flag
c 17.09.2008  ggu     comments for level = -1
c 06.12.2008  ggu     in extlev set not-existing values to flag
c 14.09.2009  ggu     new way to determine if section plot in getisec()
c
c**********************************************************
c**********************************************************
c**********************************************************
c**********************************************************

	subroutine inilev

c initializes actual level
c
c -1	bottom
c  0	integrated
c >0	level

	implicit none

	integer level3
	common /level3/level3

	real getpar

	integer icall
	save icall
	data icall /0/

	if( icall .gt. 0 ) return
	icall = 1

	level3 = nint(getpar('level'))

	end

c**********************************************************

	subroutine asklev

c asks for actual level

	implicit none

	integer level3
	common /level3/level3

	integer iauto

	integer ideflt
	real getpar

	call inilev

	iauto = nint(getpar('iauto'))

	if( iauto .eq. 0 ) then
	  level3 = ideflt(level3,'Enter level : ')
	else
	  write(6,*) 'Level used : ',level3
	end if

	end

c**********************************************************

	subroutine setlev( level )

c set actual level

	implicit none

	integer level

	integer level3
	common /level3/level3

	call inilev

	level3 = level

	end

c**********************************************************

	function getlev()

c get actual level

	implicit none

	integer getlev
	integer level3
	common /level3/level3

	call inilev

	getlev = level3

	end

c**********************************************************
c**********************************************************
c**********************************************************

        function getisec()

c is it a vertical section

        implicit none

        integer getisec
	real getpar
	character*6 vsect

	call getfnm('vsect',vsect)
	!getisec = nint(getpar('isect'))
	getisec = 0
	if( vsect .ne. ' ' ) getisec = 1

        end

c**********************************************************
c**********************************************************
c**********************************************************

	subroutine inivar

c initializes actual variable

	implicit none

	integer ivar3
	common /ivar3/ivar3
	save /ivar3/

	real getpar

	integer icall
	save icall
	data icall /0/

	if( icall .gt. 0 ) return
	icall = 1

	ivar3 = nint(getpar('ivar'))

c	ivar3 = 0	! 0 -> nothing

	end

c**********************************************************

	subroutine askvar

c asks for actual variable

	implicit none

	integer ivar3
	common /ivar3/ivar3

	integer iauto
	integer ideflt
	real getpar

	call inivar

	iauto = nint(getpar('iauto'))

	if( iauto .eq. 0 ) then
	  ivar3 = ideflt(ivar3,'Enter variable : ')
	else
	  write(6,*) 'Variable used : ',ivar3
	  write(6,*)
	end if

	end

c**********************************************************

	subroutine setvar( ivar )

c set actual variable

	implicit none

	integer ivar

	integer ivar3
	common /ivar3/ivar3

	call inivar

	ivar3 = ivar

	end

c**********************************************************

	function getvar()

c get actual variable

	implicit none

	integer getvar
	integer ivar3
	common /ivar3/ivar3

	call inivar

	getvar = ivar3

	end

c**********************************************************

	function okvar(ivar)

c shall we plot this variable ?

	implicit none

	logical okvar
        integer ivar

	integer ivar3
	common /ivar3/ivar3

	call inivar

	okvar = ivar3 .eq. ivar .or. ivar3 .eq. 0

	end

c**********************************************************
c**********************************************************
c**********************************************************

	subroutine extlev(level,nlvdim,nkn,p3,p2)

c extract level from 3d array

	implicit none

        integer ilhkv(1)
        common /ilhkv/ilhkv

	integer level		!level to extract
	integer nlvdim		!vertical dimension of p3
	integer nkn		!number of nodes
	real p3(nlvdim,nkn)	!3d array
	real p2(nkn)		!2d array filled on return

	integer k
        real flag

	if( level .gt. nlvdim ) then
	  write(6,*) 'level, nlvdim : ',level,nlvdim
	  stop 'error stop extlev: level'
	end if

        call get_flag(flag)

	if( level .le. 0 ) then
	  call intlev(nlvdim,nkn,p3,p2)		!integrate
	else
	  do k=1,nkn
	    p2(k) = flag
            if( level .le. ilhkv(k) ) p2(k) = p3(level,k)
	  end do
	end if

	end

c**********************************************************

	subroutine intlev(nlvdim,nkn,p3,p2)

c extract level from 3d array

	implicit none

	integer nlvdim		!vertical dimension of p3
	integer nkn		!number of nodes
	real p3(nlvdim,nkn)	!3d array
	real p2(nkn)		!2d array filled on return

        integer ilhkv(1)
        common /ilhkv/ilhkv

	integer k,l,lmax
	real value

	do k=1,nkn
	  lmax = ilhkv(k)
	  if( lmax .eq. 1 ) then	!2d
	    p2(k) = p3(1,k)
	  else				!primitive method of averaging
	    if( lmax .gt. nlvdim ) goto 99
	    if( lmax .le. 0 ) goto 99
	    value = 0.
	    do l=1,lmax
	      value = value + p3(l,k)
	    end do
	    p2(k) = value / lmax
	  end if
	end do

	return
   99	continue
	write(6,*) 'lmax,nlvdim : ',lmax,nlvdim
	stop 'error stop intlev : error in lmax'
	end

c**********************************************************
c**********************************************************
c**********************************************************

	subroutine set_flag(flag)

	real flag

	real flagco
	common /flagco/flagco
	save /flagco/

	flagco = flag

	end

	subroutine get_flag(flag)

	real flag

	real flagco
	common /flagco/flagco
	save /flagco/

	flag = flagco

	end

c**********************************************************
