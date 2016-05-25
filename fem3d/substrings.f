
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

	integer is,isb
	integer ie3,ie4,ie5
	integer ichafs

        iv = -1

	is = ichafs(string)
	if( is .le. 0 ) is = 1
	isb = is - 1
	ie3 = isb + 3
	ie4 = isb + 4
	ie5 = isb + 5

        if( string(is:ie4) .eq. 'mass' ) then
          iv = 0
        else if( string(is:ie5) .eq. 'level' ) then
          iv = 1
        else if( string(is:ie4) .eq. 'zeta' ) then
          iv = 1
        else if( string(is:ie3) .eq. 'vel' ) then
          iv = 2
        else if( string(is:ie5) .eq. 'trans' ) then
          iv = 3
        else if( string(is:ie4) .eq. 'bath' ) then
          iv = 5
        else if( string(is:ie5) .eq. 'depth' ) then
          iv = 5
        else if( string(is:ie3) .eq. 'cur' ) then
          iv = 6
        else if( string(is:ie5) .eq. 'speed' ) then
          iv = 6
        else if( string(is:ie3) .eq. 'dir' ) then
          iv = 7
        else if( string(is:ie4) .eq. 'conc' ) then
          iv = 10
        else if( string(is:ie4) .eq. 'conz' ) then
          iv = 10
        else if( string(is:ie3) .eq. 'sal' ) then
          iv = 11
        else if( string(is:ie4) .eq. 'temp' ) then
          iv = 12
        else if( string(is:ie4) .eq. 'oxyg' ) then
          iv = 15
        else if( string(is:ie3) .eq. 'rms' ) then
          iv = 18
        else if( string(is:ie4) .eq. 'pres' ) then
          iv = 20
        else if( string(is:ie4) .eq. 'wind' ) then
          iv = 21
        else if( string(is:ie4) .eq. 'sola' ) then
          iv = 22
        else if( string(is:ie3) .eq. 'air' ) then
          iv = 23
        else if( string(is:ie4) .eq. 'humi' ) then
          iv = 24
        else if( string(is:ie4) .eq. 'clou' ) then
          iv = 25
        else if( string(is:ie4) .eq. 'rain' ) then
          iv = 26
        else if( string(is:ie4) .eq. 'evap' ) then
          iv = 27
        else if( string(is:ie3) .eq. 'lgr' ) then
          iv = 80
        else if( string(is:ie3) .eq. 'ice' ) then
          iv = 85
!        else if( string(is:ie3) .eq. 'age' ) then
!          iv = 98
        else if( string(is:ie3) .eq. 'wrt' ) then
          iv = 99
        else if( string(is:ie5) .eq. 'renew' ) then
          iv = 99
        else if( string(is:ie4) .eq. 'resi' ) then
          iv = 99
        else if( string(is:ie4) .eq. 'ivar' ) then
	  read(string(ie4+1:),'(i5)') iv
        else if( string(is:ie3) .eq. 'var' ) then
	  read(string(ie3+1:),'(i5)') iv
        else if( string(is:ie3) .eq. 'nos' ) then
          !generic - no id
        else if( string(is:ie3) .eq. 'fem' ) then
          !generic - no id
        else if( string(is:ie4) .eq. 'elem' ) then
          !generic - no id
        else if( string .eq. ' ' ) then
          write(6,*) '*** string2ivar: no string given'
        else
          write(6,*) '*** string2ivar: cannot find string description: '
          write(6,*) string
          !write(6,*) is,isb,ie3,ie4,ie5
          !if( string(1:3) .eq. 'fem' ) stop 'error.....'
        end if

	!write(6,*) 'string2ivar: ',string(is:ie4),'   ',iv

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
        else if( iv .eq. 335 ) then
          string = 'time over threshold'
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

