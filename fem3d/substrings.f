!
! handle strings for parameters
!
! revision log :
!
! 31.08.2017    ggu     deleted old versions of subroutines
!
! notes :
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
	ivmin = 0

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

	call strings_get_ivar(string,iv)	!new call

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

	isub = 0
	call strings_get_name(iv,string,isub)		!new call

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

