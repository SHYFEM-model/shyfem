!
! handle strings for parameters
!
! revision log :
!
! 31.08.2017    ggu     deleted old versions of subroutines
! 07.10.2017    ggu     short name introduced, new generic routines
!
! contents :
!
! name		any abbreviation of variable name
! full		full name
! short		short name
! ivar		variable identification
! isub		sub-variable number
! irange	variable range
!
!	function strings_get_id_by_name(name)
!	function strings_get_id_by_ivar(ivar)
!	subroutine strings_get_full(ivar,full,isub)
!	subroutine strings_get_short(ivar,short,isub)
!
!	subroutine strings_get_full_name(name,full)
!	subroutine strings_get_full_name(ivar,full)
!	subroutine strings_get_full_name(ivar,full,isub)
!
!	subroutine strings_get_short_name(name,short)
!	subroutine strings_get_short_name(ivar,short)
!	subroutine strings_get_short_name(ivar,short,isub)
!
!	subroutine strings_get_ivar(name,ivar)
!
!       subroutine strings_add_new(name,ivar,irange)
!       subroutine strings_set_short(ivar,short)
!
!-------------------------------------------------------
!
!       subroutine string2ivar(string,ivar)
!       subroutine ivar2string(ivar,string,isub)
!	subroutine ivar2filename(ivar,filename)
!	subroutine ivar2femstring(ivar,femstring)
!
!       function compare_svars(s1,s2)
!       subroutine string_direction(string,dir)
!       function has_direction(name)
!
!-------------------------------------------------------
!
!       subroutine populate_strings
!
!-------------------------------------------------------
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

	  character*80 :: full
	  character*10 :: short
	  integer :: ivar
	  integer :: irange

	end type entry

        integer, save, private :: idlast = 0
        integer, save, private :: ndim = 0
	type(entry), save, private, allocatable :: pentry(:)

        INTERFACE strings_get_id
        MODULE PROCEDURE         strings_get_id_by_name
     +                          ,strings_get_id_by_ivar
        END INTERFACE

        INTERFACE strings_get_full_name
        MODULE PROCEDURE         strings_get_full_name_by_name
     +                          ,strings_get_full_name_by_ivar
     +                          ,strings_get_full_name_by_ivar_isub
        END INTERFACE

        INTERFACE strings_get_short_name
        MODULE PROCEDURE         strings_get_short_name_by_name
     +                          ,strings_get_short_name_by_ivar
     +                          ,strings_get_short_name_by_ivar_isub
        END INTERFACE

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

        pentry(id)%full = ' '
        pentry(id)%short = ' '
        pentry(id)%ivar = 0
        pentry(id)%irange = 0

        end subroutine strings_init_id

!******************************************************************
!******************************************************************
!******************************************************************

	function strings_get_id_by_name(name)

	integer strings_get_id_by_name
	character*(*) name

	integer id,i
	character(len=len(name)) :: string
	logical compare_svars

	call populate_strings

	string = adjustl(name)
	do i=1,len(string)
	  if( string(i:i) == '_' ) string(i:i) = ' '
	end do

	do id=1,idlast
	  if( compare_svars(pentry(id)%full,string) ) exit
	end do
	if( id > idlast ) id = 0

	strings_get_id_by_name = id

	end function strings_get_id_by_name

!******************************************************************

	function strings_get_id_by_ivar(ivar)

! given ivar finds id and sub-range

	integer strings_get_id_by_ivar
	integer ivar

	integer id,ivmin,ivmax
	logical bdebug

	call populate_strings

	bdebug = ivar == -1
	bdebug = .false.

	do id=1,idlast
	  ivmin = pentry(id)%ivar
	  ivmax = ivmin + pentry(id)%irange
	  if( bdebug ) write(6,*) id,ivar,ivmin,ivmax
	  if( ivmin == ivar ) exit
	  if( ivmin < ivar .and. ivar < ivmax ) exit
	end do
	if( id > idlast ) id = 0

	strings_get_id_by_ivar = id

	end function strings_get_id_by_ivar

!******************************************************************

	subroutine strings_get_full(ivar,full,isub)

! given ivar finds name and sub-range

	integer ivar
	character*(*) full
	integer isub		!sub-number of multi variables

	integer id

	call populate_strings

	full = ' '
	isub = 0

	id = strings_get_id(ivar)
	if( id <= 0 ) return

	full = pentry(id)%full
	isub = ivar - pentry(id)%ivar

	end subroutine strings_get_full

!******************************************************************

	subroutine strings_get_short(ivar,short,isub)

! given ivar finds short and sub-range

	integer ivar
	character*(*) short
	integer isub		!sub-number of multi variables

	integer id

	call populate_strings

	short = ' '
	isub = 0

	id = strings_get_id(ivar)
	if( id <= 0 ) return

	short = pentry(id)%short
	isub = ivar - pentry(id)%ivar

	end subroutine strings_get_short

!******************************************************************

	subroutine strings_get_full_name_by_name(name,full)

	character*(*) name
	character*(*) full

	integer ivar,isub

	full = ' '
	call strings_get_ivar(name,ivar)
	if( ivar == -1 ) return
	call strings_get_full(ivar,full,isub)

	end subroutine strings_get_full_name_by_name

!******************************************************************

	subroutine strings_get_full_name_by_ivar(ivar,full)

	integer ivar
	character*(*) full

	integer isub

	full = ' '
	call strings_get_full(ivar,full,isub)

	end subroutine strings_get_full_name_by_ivar

!******************************************************************

	subroutine strings_get_full_name_by_ivar_isub(ivar,full,isub)

	integer ivar
	character*(*) full
	integer isub

	full = ' '
	call strings_get_full(ivar,full,isub)

	end subroutine strings_get_full_name_by_ivar_isub

!******************************************************************

	subroutine strings_get_short_name_by_name(name,short)

	character*(*) name
	character*(*) short

	integer ivar,isub

	short = ' '
	call strings_get_ivar(name,ivar)
	if( ivar == -1 ) return
	call strings_get_short(ivar,short,isub)

	end subroutine strings_get_short_name_by_name

!******************************************************************

	subroutine strings_get_short_name_by_ivar(ivar,short)

	integer ivar
	character*(*) short

	integer isub

	short = ' '
	call strings_get_short(ivar,short,isub)

	end subroutine strings_get_short_name_by_ivar

!******************************************************************

	subroutine strings_get_short_name_by_ivar_isub(ivar,short,isub)

	integer ivar
	character*(*) short
	integer isub

	short = ' '
	call strings_get_short(ivar,short,isub)

	end subroutine strings_get_short_name_by_ivar_isub

!******************************************************************

	subroutine strings_get_ivar(name,ivar)

	character*(*) name
	integer ivar

	integer id

	call populate_strings

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
	
	pentry(id)%full = name
	pentry(id)%ivar = ivar
	pentry(id)%irange = irange_local

	end subroutine strings_add_new

!******************************************************************

	subroutine strings_set_short(ivar,short)

	integer ivar
	character*(*) short

	integer id

	id = strings_get_id(ivar)
	if( id <= 0 ) return

	pentry(id)%short = short

	end subroutine strings_set_short

!================================================================
	end module shyfem_strings
!================================================================

!****************************************************************
!****************************************************************
!****************************************************************

        subroutine string2ivar(string,ivar)

	use shyfem_strings

        implicit none

        character*(*) string
	integer ivar

	call strings_get_ivar(string,ivar)	!new call

	if( ivar < 0 ) then
          if( string .eq. ' ' ) then
            write(6,*) '*** string2ivar: no string given'
          else
	    write(6,*) '*** string2ivar: cannot find string description:'
            write(6,*) trim(string),'   (',trim(string),')'
          end if
	end if

	end

!****************************************************************

        subroutine ivar2string(ivar,string,isub)

	use shyfem_strings

        implicit none

        integer ivar
        character*(*) string
	integer isub

	isub = 0
	call strings_get_full(ivar,string,isub)		!new call

	end

!****************************************************************

	subroutine ivar2filename(ivar,filename)

! use this to create unique filename

	use shyfem_strings

	implicit none

	integer ivar
	character*(*) filename

	integer isub,i
	character*80 string

	call strings_get_short_name(ivar,string,isub)

	filename = string
	if( isub > 0 ) then
	  write(string,'(i4)') isub
	  do i=1,4
	    if( string(i:i) == ' ' ) string(i:i) = '0'
	  end do
	  string(1:1) = '_'
	  filename = trim(filename) // string(1:4)
	end if

	end

!****************************************************************

	subroutine ivar2femstring(ivar,femstring)

! use this for writing into FEM file

	use shyfem_strings

	implicit none

	integer ivar
	character*(*) femstring

	integer isub,i
	character*80 string

	call strings_get_full_name(ivar,string,isub)

	femstring = string
	if( isub > 0 ) then
	  write(string,'(i4)') isub
	  femstring = femstring // string(1:4)
	end if

	end

!****************************************************************
!****************************************************************
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

	function has_direction(name)

! check if directional variable

	use shyfem_strings

	implicit none

	logical has_direction
	character*(*) name

	integer ivar

	call strings_get_ivar(name,ivar)

	has_direction = ( ivar == 2 .or. ivar == 3 .or. ivar == 21 )

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

	logical, save :: bread = .false.

	if( bread ) return
	bread = .true.

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
	call strings_add_new('salt',11)
	call strings_add_new('temperature',12)
	call strings_add_new('density',13)
	call strings_add_new('rho',13)
	call strings_add_new('oxygen',15)
	call strings_add_new('discharge',16)
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

	call strings_add_new('waves (general)',230)
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

!---------------------------------------------------------------------

	call strings_set_short(0,'mass')
	call strings_set_short(1,'zeta')
	call strings_set_short(2,'vel')
	call strings_set_short(3,'transp')
	call strings_set_short(5,'bathy')
	call strings_set_short(6,'speed')
	call strings_set_short(7,'dir')
	call strings_set_short(10,'conz')
	call strings_set_short(11,'salt')
	call strings_set_short(12,'temp')
	call strings_set_short(13,'rho')
	call strings_set_short(15,'oxy')
	call strings_set_short(16,'disch')
	call strings_set_short(18,'rms')

	call strings_set_short(20,'airp')
	call strings_set_short(21,'wind')
	call strings_set_short(22,'srad')
	call strings_set_short(23,'airt')
	call strings_set_short(24,'rhum')
	call strings_set_short(25,'cc')
	call strings_set_short(26,'rain')
	call strings_set_short(27,'evap')

	call strings_set_short(60,'bstress')
	call strings_set_short(75,'index')
	call strings_set_short(76,'type')
	call strings_set_short(80,'lgr')
	call strings_set_short(85,'ice')
	call strings_set_short(97,'timeot')
	call strings_set_short(98,'age')
	call strings_set_short(99,'wrt')

	call strings_set_short(230,'waves')
	call strings_set_short(231,'wheight')
	call strings_set_short(232,'wper')
	call strings_set_short(233,'wdir')
	call strings_set_short(234,'worb')
	call strings_set_short(235,'wpeak')

	call strings_set_short(300,'conz')
	call strings_set_short(30,'conz')

	call strings_set_short(700,'biop')
	call strings_set_short(720,'bios')
	call strings_set_short(730,'biosf')

	call strings_set_short(800,'sedi')

	end

!****************************************************************
!****************************************************************
!****************************************************************

