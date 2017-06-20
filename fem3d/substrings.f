!
! variable ids for consecutive variables:
!
!	1-199	single variables
!	230	waves
!	250	mercury
!	300	conz
!	400	aquabc
!	500	toxi
!	600	bfm
!	700	eutro
!	800	sediments
!	
!****************************************************************
!****************************************************************
!****************************************************************

!================================================================
	module shyfem_strings
!================================================================

	implicit none

	type, private :: entry

	  character*80 :: name
	  integer :: ivar
	  integer :: irange

	end type entry

	logical, save, private :: bread = .true.	!still to populate strings
        integer, save, private :: idlast = 0
        integer, save, private :: ndim = 0
	type(entry), save, private, allocatable :: pentry(:)

!================================================================
	contains
!================================================================

        subroutine strings_init_alloc

        type(entry), allocatable :: paux(:)

        if( ndim == 0 ) then
          ndim = 10
          allocate(pentry(ndim))
        else
          allocate(paux(2*ndim))
          paux(1:ndim) = pentry(1:ndim)
          call move_alloc(paux,pentry)
          ndim = ndim*2
        end if

        end subroutine strings_init_alloc

!******************************************************************

        subroutine strings_init_new_id(id)

        integer id

        idlast = idlast + 1
        if( idlast > ndim ) then
          call strings_init_alloc
        end if
        id = idlast

        call strings_init_id(id)

        end subroutine strings_init_new_id

!******************************************************************

        subroutine strings_init_id(id)

        integer id

        if( id > ndim ) then
          stop 'error stop strings_init_id: ndim'
        end if

        pentry(id)%name = ' '
        pentry(id)%ivar = 0
        pentry(id)%irange = 0

        end subroutine strings_init_id

!******************************************************************

	function strings_get_id(name)

	integer strings_get_id
	character*(*) name

	integer id,i
	character(len=len(name)) :: string
	logical compare_svars

	string = adjustl(name)
	do i=1,len(string)
	  if( string(i:i) == '_' ) string(i:i) = ' '
	end do

	do id=1,idlast
	  if( compare_svars(pentry(id)%name,string) ) exit
	end do
	if( id > idlast ) id = 0

	strings_get_id = id

	end function strings_get_id

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine strings_get_name(ivar,name,isub)

	integer ivar
	character*(*) name
	integer isub		!sub-number of multi variables

	integer id,ivmin,ivmax
	logical bdebug

	if( bread ) then
	  call populate_strings
	  bread = .false.
	end if

	bdebug = ivar == -1
	bdebug = .false.
	name = ' '
	isub = 0

	do id=1,idlast
	  ivmin = pentry(id)%ivar
	  ivmax = ivmin + pentry(id)%irange
	  if( bdebug ) write(6,*) id,ivar,ivmin,ivmax
	  if( ivmin == ivar ) exit
	  if( ivmin < ivar .and. ivar < ivmax ) exit
	end do
	if( id > idlast ) return

	name = pentry(id)%name
	isub = ivar - ivmin

	end subroutine strings_get_name

!******************************************************************

	subroutine strings_get_full_name(name,fullname)

	character*(*) name
	character*(*) fullname

	integer ivar,isub

	fullname = ' '
	call strings_get_ivar(name,ivar)
	if( ivar == -1 ) return
	call strings_get_name(ivar,fullname,isub)

	end subroutine strings_get_full_name

!******************************************************************

	subroutine strings_get_ivar(name,ivar)

	character*(*) name
	integer ivar

	integer id

	if( bread ) then
	  call populate_strings
	  bread = .false.
	end if

	ivar = -1
	if( name == ' ' ) return

	id = strings_get_id(name)
	if( id == 0 ) return

	ivar = pentry(id)%ivar

	end subroutine strings_get_ivar

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine strings_add_new(name,ivar,irange)

	character*(*) name
	integer ivar
	integer, optional :: irange

	integer id,irange_local

	id = strings_get_id(name)
	if( id /= 0 ) then
	  write(6,*) id,ivar,'  ',name
	  stop 'error stop strings_add_new: name already present'
	end if

	call strings_init_new_id(id)

	irange_local = 0
	if( present(irange) ) irange_local = irange
	
	pentry(id)%name = name
	pentry(id)%ivar = ivar
	pentry(id)%irange = irange_local

	end subroutine strings_add_new

!================================================================
	end module shyfem_strings
!================================================================

!****************************************************************
!****************************************************************
!****************************************************************

        subroutine string2ivar(string,iv)

        implicit none

        character*(*) string
	integer iv

        call string2ivar_n(string,iv)

	if( iv < 0 ) then
          if( string .eq. ' ' ) then
            write(6,*) '*** string2ivar: no string given'
          else
	    write(6,*) '*** string2ivar: cannot find string description:'
            write(6,*) trim(string),'   (',trim(string),')'
          end if
	end if

	end

!****************************************************************

        subroutine string2ivar_n(string,iv)

	use shyfem_strings

        implicit none

        character*(*) string
	integer iv

	integer ivar

	call string2ivar_intern(string,iv)	!old call - delete
	call strings_get_ivar(string,ivar)	!new call

	if( iv /= ivar .and. iv > 0 ) then
	  write(6,*) 'string: ',trim(string)
	  write(6,*) 'iv,ivar: ',iv,ivar
	  write(6,*) 'error stop string2ivar_n: internal error (1)'
	  !stop 'error stop string2ivar_n: internal error (1)'
	end if

	end

!****************************************************************

        subroutine string2ivar_intern(string,iv)

! old routine - do not use anymore
!
! interprets string to associate a variable number iv
!
! see below for possible string names
!
! the special name ivar# can be used to directtly give the variable number #

        implicit none

        character*(*) string
        integer iv

	integer is,isb,i
	integer ie3,ie4,ie5,ie6,ie8,ie11,ie16
	integer ichafs
	character*80 s

        iv = -1

	s = string
	do i=1,len(string)
	  if( s(i:i) == '_' ) s(i:i) = ' '	!convert '_' to ' '
	end do

	is = ichafs(s)
	if( is .le. 0 ) is = 1
	isb = is - 1
	ie3 = isb + 3
	ie4 = isb + 4
	ie5 = isb + 5
	ie6 = isb + 6
	ie8 = isb + 8
	ie11 = isb + 11
	ie16 = isb + 16

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
        else if( s(is:ie5) .eq. 'index' ) then
          iv = 75
        else if( s(is:ie4) .eq. 'type' ) then
          iv = 76
        else if( s(is:ie3) .eq. 'lgr' ) then
          iv = 80
        else if( s(is:ie3) .eq. 'ice' ) then
          iv = 85
        else if( s(is:ie16) .eq. 'time over thresh' ) then
          iv = 97
        else if( s(is:ie3) .eq. 'age' ) then
          iv = 98
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
        else if( s(is:ie8) .eq. 'bottom stress' ) then
          iv = 238
        else if( s(is:ie4) .eq. 'sedi' ) then
          iv = 800
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

        subroutine ivar2string(iv,string,isub)

	use shyfem_strings

        implicit none

        integer iv
        character*(*) string
	integer isub

	character(len=len(string)) :: s1

	isub = 0
	call ivar2string_intern(iv,s1)			!old call - delete
	call strings_get_name(iv,string,isub)		!new call

	if( s1 /= string .and. s1 /= ' ') then
	  write(6,*) 'ivar = ',iv
	  write(6,*) 's1 = ',trim(s1)
	  write(6,*) 's2 = ',trim(string)
	  write(6,*) 'error stop ivar2string: internal error (1)'
	  !stop 'error stop ivar2string: internal error (1)'
	end if

	end

!****************************************************************

        subroutine ivar2string_intern(iv,string)

! old routine - do not use anymore

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
        else if( iv .eq. 75 ) then
          string = 'index'
        else if( iv .eq. 76 ) then
          string = 'type'
        else if( iv .eq. 85 ) then
          string = 'ice cover'
        else if( iv .eq. 97 ) then
          string = 'time over threshold'
        else if( iv .eq. 98 ) then
          string = 'age'
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
        else if( iv .eq. 238 ) then
          string = 'bottom stress'
        else if( iv > 30 .and. iv < 50 ) then
          string = 'concentration (multi)'
        else if( iv > 700 .and. iv < 720 ) then
          string = 'weutro (pelagic)'
        else if( iv > 720 .and. iv < 730 ) then
          string = 'weutro (sediment)'
        else if( iv > 730 .and. iv < 740 ) then
          string = 'weutro (shell fish)'
        else if( iv >= 800 .and. iv < 900 ) then
          string = 'sediments'
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
          call string2ivar_n(strings(i),ivars(i))
	end do

	end

!****************************************************************

	function compare_svars(s1,s2)

! compares two strings if they indicate the same variable

	implicit none

	logical compare_svars
	character*(*) s1,s2

	integer l1,l2,l

	compare_svars = .false.

	l1 = len_trim(s1)
	l2 = len_trim(s2)
	l = min(l1,l2)

	if( l == 0 ) return

	compare_svars = ( s1(1:l) == s2(1:l) )

	end

!****************************************************************

	function has_direction(name)

! check if directional variable

	use shyfem_strings

	implicit none

	logical has_direction
	character*(*) name

	integer iv

	call strings_get_ivar(name,iv)

	has_direction = ( iv == 2 .or. iv == 3 .or. iv == 21 )

	end

!****************************************************************
!****************************************************************
!****************************************************************

	subroutine populate_strings

! populates string information
!
! put prefered name first if more than one names apply

	use shyfem_strings

	implicit none

	call strings_add_new('mass field',0)
	call strings_add_new('water level',1)
	call strings_add_new('level',1)
	call strings_add_new('zeta',1)
	call strings_add_new('velocity',2)
	call strings_add_new('transport',3)
	call strings_add_new('bathymetry',5)
	call strings_add_new('depth',5)
	call strings_add_new('current speed',6)
	call strings_add_new('speed',6)
	call strings_add_new('current direction',7)
	call strings_add_new('direction',7)
	call strings_add_new('generic tracer',10)
	call strings_add_new('tracer',10)
	call strings_add_new('salinity',11)
	call strings_add_new('temperature',12)
	call strings_add_new('density',13)
	call strings_add_new('rho',13)
	call strings_add_new('oxygen',15)
	call strings_add_new('rms velocity',18)
	call strings_add_new('rms speed',18)

	call strings_add_new('atmospheric pressure',20)
	call strings_add_new('air pressure',20)
	call strings_add_new('pressure',20)
	call strings_add_new('wind velocity',21)
	call strings_add_new('solar radiation',22)
	call strings_add_new('air temperature',23)
	call strings_add_new('humidity (relative)',24)
	call strings_add_new('cloud cover',25)
	call strings_add_new('rain',26)
	call strings_add_new('evaporation',27)

	call strings_add_new('bottom stress',60)

	call strings_add_new('index',75)
	call strings_add_new('type',76)
	call strings_add_new('lgr',80)
	call strings_add_new('ice cover',85)
	call strings_add_new('time over threshold',97)
	call strings_add_new('age',98)
	call strings_add_new('renewal time',99)
	call strings_add_new('residence time',99)
	call strings_add_new('wrt',99)

	call strings_add_new('wave height (significant)',231)
	call strings_add_new('wave period (mean)',232)
	call strings_add_new('wave direction',233)
	call strings_add_new('wave orbital velocity',234)
	call strings_add_new('wave peak period',235)

	call strings_add_new('concentration (multi)',300,100)	!new numbering
	call strings_add_new('concentration (multi old)',30,20)

	call strings_add_new('weutro (pelagic)',700,20)
	call strings_add_new('weutro (sediment)',720,10)
	call strings_add_new('weutro (shell fish)',730,10)

	call strings_add_new('sediments',800,100)

	call strings_add_new('var',-9)		!special treatment
	call strings_add_new('ivar',-9)

	end

!****************************************************************
!****************************************************************
!****************************************************************

