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
c 18.08.2011  ggu     make vsect bigger
c 31.08.2011  ggu     new plotting eos
c 23.02.2012  ccf     allow plotting also for last layer
c 13.06.2013  ggu     scans varnam to decide what to plot
c 05.09.2013  ggu     handle variable choice better
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
	character*80 vsect

	call getfnm('vsect',vsect)
	!getisec = nint(getpar('isect'))
	getisec = 0
	if( vsect .ne. ' ' ) getisec = 1

        end

c**********************************************************
c**********************************************************
c**********************************************************

	subroutine inivar

c initializes actual variable to be plotted (internal routine)

	implicit none

	integer ivar3
	common /ivar3/ivar3
	save /ivar3/

	integer ivar,ivar_name
	character*80 name
	real getpar

	integer icall
	save icall
	data icall /0/

	if( icall .gt. 0 ) return
	icall = 1

	ivar = nint(getpar('ivar'))
	ivar_name = 0
	name = ' '
	call getfnm('varnam',name)
	if( name .ne. ' ' ) then
	  call string2ivar(name,ivar_name)
	end if

	write(6,*) 'var name: ',name(1:50)
	write(6,*) 'var id  : ',ivar_name

	if( ivar .gt. 0 .and. ivar_name .gt. 0 ) then
	  if( ivar .ne. ivar_name ) then
	    write(6,*) 'You cannot give both ivar and varnam'
	    write(6,*) 'ivar = ',ivar
	    write(6,*) 'varnam = ',name(1:30)
	    write(6,*) 'ivar_name = ',ivar_name
	    stop 'error stop inivar: non compatible variables required'
	  end if
	end if

	if( ivar .gt. 0 ) then
	  ivar3 = ivar
	else if( ivar_name .gt. 0 ) then
	  ivar3 = ivar_name
	else
	  ivar3 = 0
	end if

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
	  ivar3 = ideflt(ivar3,'Enter variable id: ')
	end if

	write(6,*) 'askvar: Variable used = ',ivar3
	write(6,*)

	end

c**********************************************************

	subroutine checkvar(ivar)

c checks what variable has to be plotted
c returns in ivar the variable to be plotted

	implicit none

	integer ivar

	integer ivar3
	integer getvar

	call inivar

	if( ivar .gt. 0 ) then	!ivar given - must be equal to ivar3
	  ivar3 = getvar()
	  if( ivar3 .gt. 0 .and. ivar3 .ne. ivar ) goto 99
          call setvar(ivar)
        else
          call askvar
          ivar3 = getvar()
          ivar = getvar()
        end if

	if( ivar .le. 0 ) then
	  write(6,*) 'Do not know what to plot: ivar = ',ivar
	  stop 'error stop checkvar: no ivar given'
	end if

	return
   99	continue
	write(6,*) 'ivar3 = ',ivar3
	write(6,*) 'ivar  = ',ivar
	stop 'error stop checkvar: different values of ivar3 and ivar'
	end

c**********************************************************

	subroutine setvar(ivar)

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

	okvar = ivar3 .eq. ivar .or. ivar3 .le. 0

	end

c**********************************************************

	function read_var()

c reads number or string from STDIN - converts string to number
c is not used actually

	implicit none

	integer read_var

	integer ivar,ios
	character*60 line

	write(6,*) 'Enter variable : '
	read(5,'(a)') line

	read(line,'(i10)',iostat=ios) ivar

	if( ios .ne. 0 ) then
	  call string2ivar(line,ivar)
	end if

	read_var = ivar

	end

c**********************************************************
c**********************************************************
c**********************************************************

	subroutine extnlev(level,nlvddi,nkn,p3,p2)

c extract level from 3d array (nodes)

	use levels

	implicit none

	integer level		!level to extract
	integer nlvddi		!vertical dimension of p3
	integer nkn		!number of nodes
	real p3(nlvddi,nkn)	!3d array
	real p2(nkn)		!2d array filled on return

	include 'param.h'

	call extlev(level,nlvddi,nkn,ilhkv,p3,p2)

	end

c**********************************************************

	subroutine extelev(level,nlvddi,nel,p3,p2)

c extract level from 3d array (elements)

	use levels

	implicit none

	integer level		!level to extract
	integer nlvddi		!vertical dimension of p3
	integer nel		!number of elements
	real p3(nlvddi,nel)	!3d array
	real p2(nel)		!2d array filled on return

	include 'param.h'

	call extlev(level,nlvddi,nel,ilhv,p3,p2)

	end

c**********************************************************

	subroutine extlev(level,nlvddi,n,ilv,p3,p2)

c extract level from 3d array

	implicit none

	integer level		!level to extract
	integer nlvddi		!vertical dimension of p3
	integer n		!number values
        integer ilv(n)
	real p3(nlvddi,n)	!3d array
	real p2(n)		!2d array filled on return

	integer i,lmax,lact
        real flag

	if( level .gt. nlvddi ) then
	  write(6,*) 'level, nlvddi : ',level,nlvddi
	  stop 'error stop extlev: level'
	end if

        call get_flag(flag)

        if( level .lt. -1 ) then
	  stop 'error stop extlev: internal error (1)'
	end if

	if( level .eq. 0 ) then
	  call intlev(nlvddi,n,ilv,p3,p2)		!integrate
	else
	  do i=1,n
	    lmax = ilv(i)
	    lact = level
	    if( lact .eq. -1 ) lact = lmax
	    p2(i) = flag
            if( lact .le. ilv(i) ) p2(i) = p3(lact,i)
	  end do
	end if

	end

c**********************************************************

	subroutine intlev(nlvddi,n,ilv,p3,p2)

c integrate over water column

	implicit none

	integer nlvddi		!vertical dimension of p3
	integer n		!number of nodes
        integer ilv(n)
	real p3(nlvddi,n)	!3d array
	real p2(n)		!2d array filled on return

	integer i,l,lmax
	real value

	do i=1,n
	  lmax = ilv(i)
	  if( lmax .eq. 1 ) then	!2d
	    p2(i) = p3(1,i)
	  else				!primitive method of averaging
	    if( lmax .gt. nlvddi ) goto 99
	    if( lmax .le. 0 ) goto 99
	    value = 0.
	    do l=1,lmax
	      value = value + p3(l,i)
	    end do
	    p2(i) = value / lmax
	  end if
	end do

	return
   99	continue
	write(6,*) 'lmax,nlvddi : ',lmax,nlvddi
	stop 'error stop intlev : error in lmax'
	end

c**********************************************************
c**********************************************************
c**********************************************************

	subroutine set_flag(flag)

	implicit none

	real flag

	real flagco
	common /flagco/flagco
	save /flagco/

	flagco = flag

	end

	subroutine get_flag(flag)

	implicit none

	real flag

	real flagco
	common /flagco/flagco
	save /flagco/

	flag = flagco

	end

c**********************************************************

