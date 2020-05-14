
!--------------------------------------------------------------------------
!
!    Copyright (C) 1999,2002,2004-2005,2008-2011,2013-2016  Georg Umgiesser
!    Copyright (C) 2018-2019  Georg Umgiesser
!    Copyright (C) 2012  Christian Ferrarin
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

c interactive routines for plotsim
c
c revision log :
c
c 12.02.1999	ggu	adapted to auto mode
c 29.01.2002	ggu	new routine getisec()
c 17.03.2004	ggu	new routine okvar()
c 02.03.2005	ggu	new routines set_flag and get_flag
c 17.09.2008	ggu	comments for level = -1
c 06.12.2008	ggu	in extlev set not-existing values to flag
c 14.09.2009	ggu	new way to determine if section plot in getisec()
c 23.03.2010	ggu	changed v6.1.1
c 18.08.2011	ggu	make vsect bigger
c 31.08.2011	ggu	new plotting eos
c 01.09.2011	ggu	changed VERS_6_1_32
c 23.02.2012	ccf	allow plotting also for last layer
c 13.06.2013	ggu	scans varnam to decide what to plot
c 19.06.2013	ggu	changed VERS_6_1_66
c 05.09.2013	ggu	handle variable choice better
c 12.09.2013	ggu	changed VERS_6_1_67
c 28.01.2014	ggu	changed VERS_6_1_71
c 21.10.2014	ggu	changed VERS_7_0_3
c 05.12.2014	ggu	changed VERS_7_0_8
c 19.01.2015	ggu	changed VERS_7_1_2
c 19.01.2015	ggu	changed VERS_7_1_3
c 05.05.2015	ggu	changed VERS_7_1_10
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 18.12.2015	ggu	changed VERS_7_3_17
c 25.05.2016	ggu	changed VERS_7_5_10
c 18.12.2018	ggu	changed VERS_7_5_52
c 21.05.2019	ggu	changed VERS_7_5_62
c
c**********************************************************
c**********************************************************
c**********************************************************
c**********************************************************

!==================================================================
        module mod_plot
!==================================================================

c initializes actual level
c
c -1	bottom
c  0	integrated
c >0	level

	implicit none

	integer, save :: level3 = 0
	integer, save :: ivsect = 0
	integer, save :: ivar3 = 0
	real, save :: flagco = -999.

	integer, save :: icall = 0

!==================================================================
        end module mod_plot
!==================================================================

	subroutine init_plot

	use mod_plot
	use para

	implicit none

	integer ivar
	character*80 name,vsect
	real getpar

	if( icall > 0 ) return

	if( para_has_name('level') ) then
	  level3 = nint(getpar('level'))
	end if

	if( para_has_name('vsect') ) then
	  call getfnm('vsect',vsect)
	  ivsect = 0
	  if( vsect .ne. ' ' ) ivsect = 1
	end if

	ivar = 0
	name = ' '
	if( para_has_name('ivar') ) then
	  ivar = nint(getpar('ivar'))
	end if
	if( para_has_name('varnam') ) then
	  call getfnm('varnam',name)
	end if
	if( ivar > 0 .and. name /= ' ' ) then
	  write(6,*) 'You cannot give both ivar and varnam'
	  stop 'error stop init_plot: non compatible variables'
	end if
	if( name .ne. ' ' ) call string2ivar(name,ivar)
	ivar3 = ivar

	icall = 1

	end

c**********************************************************

	subroutine setlev( level )

c set actual level

	use mod_plot

	implicit none

	integer level

	level3 = level

	end

c**********************************************************

	function getlev()

c get actual level

	use mod_plot

	implicit none

	integer getlev

	getlev = level3

	end

c**********************************************************
c**********************************************************
c**********************************************************

        function getisec()

c is it a vertical section

	use mod_plot

        implicit none

        integer getisec

	getisec = ivsect

        end

c**********************************************************
c**********************************************************
c**********************************************************

	subroutine setvar(ivar)

c set actual variable

	use mod_plot

	implicit none

	integer ivar

	ivar3 = ivar

	end

c**********************************************************

	function getvar()

c get actual variable

	use mod_plot

	implicit none

	integer getvar

	getvar = ivar3

	end

c**********************************************************

	function okvar(ivar)

c shall we plot this variable ?

	use mod_plot

	implicit none

	logical okvar
        integer ivar

	okvar = ivar3 .eq. ivar .or. ivar3 .le. 0

	end

c**********************************************************

       subroutine checkvar(ivar)

c checks what variable has to be plotted
c returns in ivar the variable to be plotted

       implicit none

       integer ivar

       integer ivar3
       integer getvar

       if( ivar .gt. 0 ) then  !ivar given - must be equal to ivar3
         ivar3 = getvar()
         if( ivar3 .gt. 0 .and. ivar3 .ne. ivar ) goto 99
         call setvar(ivar)
       else
         ivar = getvar()
       end if

       if( ivar .le. 0 ) then
         write(6,*) 'Do not know what to plot: ivar = ',ivar
         stop 'error stop checkvar: no ivar given'
       end if

       return
   99  continue
       write(6,*) 'ivar3 = ',ivar3
       write(6,*) 'ivar  = ',ivar
       stop 'error stop checkvar: different values of ivar3 and ivar'
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

	lmax = 0

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

	use mod_plot

	implicit none

	real flag

	flagco = flag

	end

	subroutine get_flag(flag)

	use mod_plot

	implicit none

	real flag

	flag = flagco

	end

c**********************************************************

