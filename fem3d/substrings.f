!
! variable ids for consecutive variables:
!
!	1-199	single variables
!	200	eutro
!	230	waves
!	250	mercury
!	300	conz
!	400	aquabc
!	500	toxi
!	600	bfm
!	
!****************************************************************
!****************************************************************
!****************************************************************

        subroutine string2ivar(string,iv)

! interprets string to associate a variable number iv
!
! see below for possible string names
!
! the special name ivar# can be used to directtly give the variable number #

        implicit none

        character*(*) string
        integer iv

	integer is,isb,i
	integer ie3,ie4,ie5,ie6,ie8,ie11
	integer ichafs
	character*80 s

        iv = -1

	s = string
	do i=1,len(string)
	  if( s(i:i) == '_' ) s(i:i) = ' '	!convert '_' to ' '
	end do
	!write(6,*) 'checking: ',trim(s)

	is = ichafs(s)
	if( is .le. 0 ) is = 1
	isb = is - 1
	ie3 = isb + 3
	ie4 = isb + 4
	ie5 = isb + 5
	ie6 = isb + 6
	ie8 = isb + 8
	ie11 = isb + 11

        if( s(is:ie4) .eq. 'mass' ) then
          iv = 0
        else if( s(is:ie5) .eq. 'level' ) then
          iv = 1
        else if( s(is:ie11) .eq. 'water level' ) then
          iv = 1
        else if( s(is:ie4) .eq. 'zeta' ) then
          iv = 1
        else if( s(is:ie3) .eq. 'vel' ) then
          iv = 2
        else if( s(is:ie5) .eq. 'trans' ) then
          iv = 3
        else if( s(is:ie4) .eq. 'bath' ) then
          iv = 5
        else if( s(is:ie5) .eq. 'depth' ) then
          iv = 5
        else if( s(is:ie3) .eq. 'cur' ) then
          iv = 6
        else if( s(is:ie5) .eq. 'speed' ) then
          iv = 6
        else if( s(is:ie3) .eq. 'dir' ) then
          iv = 7
        else if( s(is:ie4) .eq. 'conc' ) then
          iv = 10
        else if( s(is:ie4) .eq. 'conz' ) then
          iv = 10
        else if( s(is:ie3) .eq. 'sal' ) then
          iv = 11
        else if( s(is:ie4) .eq. 'temp' ) then
          iv = 12
        else if( s(is:ie3) .eq. 'rho' ) then
          iv = 13
        else if( s(is:ie4) .eq. 'dens' ) then
          iv = 13
        else if( s(is:ie4) .eq. 'oxyg' ) then
          iv = 15
        else if( s(is:ie3) .eq. 'rms' ) then
          iv = 18
        else if( s(is:ie4) .eq. 'pres' ) then
          iv = 20
        else if( s(is:ie4) .eq. 'wind' ) then
          iv = 21
        else if( s(is:ie4) .eq. 'sola' ) then
          iv = 22
        else if( s(is:ie3) .eq. 'air' ) then
          iv = 23
        else if( s(is:ie4) .eq. 'humi' ) then
          iv = 24
        else if( s(is:ie4) .eq. 'clou' ) then
          iv = 25
        else if( s(is:ie4) .eq. 'rain' ) then
          iv = 26
        else if( s(is:ie4) .eq. 'evap' ) then
          iv = 27
        else if( s(is:ie3) .eq. 'lgr' ) then
          iv = 80
        else if( s(is:ie3) .eq. 'ice' ) then
          iv = 85
!        else if( s(is:ie3) .eq. 'age' ) then
!          iv = 98
        else if( s(is:ie3) .eq. 'wrt' ) then
          iv = 99
        else if( s(is:ie5) .eq. 'renew' ) then
          iv = 99
        else if( s(is:ie4) .eq. 'resi' ) then
          iv = 99
        else if( s(is:ie6) .eq. 'wave h' ) then
          iv = 231
        else if( s(is:ie8) .eq. 'wave per' ) then
          iv = 232
        else if( s(is:ie6) .eq. 'wave d' ) then
          iv = 233
        else if( s(is:ie6) .eq. 'wave o' ) then
          iv = 234
        else if( s(is:ie8) .eq. 'wave pea' ) then
          iv = 235
        else if( s(is:ie4) .eq. 'ivar' ) then
	  read(s(ie4+1:),'(i5)') iv
        else if( s(is:ie3) .eq. 'var' ) then
	  read(s(ie3+1:),'(i5)') iv
        else if( s(is:ie3) .eq. 'nos' ) then
          !generic - no id
        else if( s(is:ie3) .eq. 'fem' ) then
          !generic - no id
        else if( s(is:ie4) .eq. 'elem' ) then
          !generic - no id
        else if( s .eq. ' ' ) then
          write(6,*) '*** string2ivar: no string given'
        else
          write(6,*) '*** string2ivar: cannot find string description: '
          write(6,*) s
          !write(6,*) is,isb,ie3,ie4,ie5
          !if( string(1:3) .eq. 'fem' ) stop 'error.....'
        end if

	!write(6,*) 'string2ivar: ',string(is:ie4),'   ',iv

        end

!****************************************************************

	subroutine string_direction(string,dir)

c finds direction if vector

	implicit none

	character(*) string,dir

	integer l

	l = len_trim(string)

	if( string(l-1:l) == ' x' ) then
	  dir = 'x'
	else if( string(l-1:l) == ' y' ) then
	  dir = 'y'
	else
	  dir = ' '
	end if

	end

!****************************************************************

        subroutine ivar2string(iv,string)

        implicit none

        integer iv
        character*(*) string

        string = ' '

        if( iv .eq. 0 ) then
          string = 'mass field'
        else if( iv .eq. 1 ) then
          string = 'water level'
        else if( iv .eq. 2 ) then
          string = 'velocity'
        else if( iv .eq. 3 ) then
          string = 'transport'
        else if( iv .eq. 5 ) then
          string = 'bathymetry'
        else if( iv .eq. 10 ) then
          string = 'generic tracer'
        else if( iv .eq. 30 ) then
          string = 'generic tracer'
        else if( iv .eq. 11 ) then
          string = 'salinity'
        else if( iv .eq. 12 ) then
          string = 'temperature'
        else if( iv .eq. 13 ) then
          string = 'density'
        else if( iv .eq. 20 ) then
          string = 'atmospheric pressure'
        else if( iv .eq. 26 ) then
          string = 'rain'
        else if( iv .eq. 85 ) then
          string = 'ice cover'
!        else if( iv .eq. 98 ) then
!          string = 'age'
        else if( iv .eq. 99 ) then
          string = 'renewal time'
        else if( iv .eq. 231 ) then
          string = 'wave height (significant)'
        else if( iv .eq. 232 ) then
          string = 'wave period (mean)'
        else if( iv .eq. 233 ) then
          string = 'wave direction'
        else if( iv .eq. 234 ) then
          string = 'wave orbital velocity'
        else if( iv .eq. 235 ) then
          string = 'wave peak period'
        else if( iv .eq. 335 ) then
          string = 'time over threshold'
        else if( iv > 30 .and. iv < 50 ) then
          string = 'concentration (multi)'
        else
          !string = '*** cannot find description'
          !write(6,*) '*** cannot find description for variable: '
          !write(6,*) iv
	  !stop 'error stop ivar2string: no description'
        end if

        end

!****************************************************************
!****************************************************************
!****************************************************************

	subroutine get_vars_from_string(nvar,strings,ivars)

c gets var numbers from string description

	implicit none

	integer nvar
	character*(*) strings(nvar)
	integer ivars(nvar)

	integer i

	do i=1,nvar
          call string2ivar(strings(i),ivars(i))
	end do

	end

!****************************************************************

