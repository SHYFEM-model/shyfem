
!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2018  Georg Umgiesser
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

! handle strings for parameters
!
! revision log :
!
! 25.05.2016	ggu	changed VERS_7_5_10
! 30.05.2016	ggu	changed VERS_7_5_11
! 10.06.2016	ggu	changed VERS_7_5_13
! 14.06.2016	ggu	changed VERS_7_5_14
! 17.06.2016	ggu	changed VERS_7_5_15
! 27.06.2016	ggu	changed VERS_7_5_16
! 09.09.2016	ggu	changed VERS_7_5_17
! 30.09.2016	ggu	changed VERS_7_5_18
! 05.10.2016	ggu	changed VERS_7_5_19
! 13.02.2017	ggu	changed VERS_7_5_23
! 09.05.2017	ggu	changed VERS_7_5_26
! 25.05.2017	ggu	changed VERS_7_5_28
! 13.06.2017	ggu	changed VERS_7_5_29
! 11.07.2017	ggu	changed VERS_7_5_30
! 31.08.2017	ggu	deleted old versions of subroutines
! 02.09.2017	ggu	changed VERS_7_5_31
! 07.10.2017	ggu	short name introduced, new generic routines
! 04.11.2017	ggu	changed VERS_7_5_34
! 14.11.2017	ggu	changed VERS_7_5_36
! 05.12.2017	ggu	changed VERS_7_5_39
! 24.01.2018	ggu	changed VERS_7_5_41
! 22.02.2018	ggu	changed VERS_7_5_42
! 03.04.2018	ggu	changed VERS_7_5_43
! 06.07.2018	ggu	changed VERS_7_5_48
! 25.10.2018	ggu	changed VERS_7_5_51
! 14.02.2019	ggu	changed VERS_7_5_56
! 16.02.2019	ggu	populate_strings declared as recursive
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
!	300	conc
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

        logical, save, private :: bpopulate = .true.    !must still populate?
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

	subroutine strings_check_consistency

	integer id

	do id=1,idlast
	  if( pentry(id)%ivar < 0 ) cycle
	  if( pentry(id)%full==' '.neqv.pentry(id)%short==' ' ) then
	    write(6,*) 'strings entries not consistent'
	    write(6,*) 'either full or short not given'
	    write(6,*) 'full  = ',trim(pentry(id)%full)
	    write(6,*) 'short = ',trim(pentry(id)%short)
	    stop 'error stop strings_check_consistency'
	  end if
	end do

	end subroutine strings_check_consistency

!******************************************************************

	subroutine strings_add_new(name,ivar,irange)

	character*(*) name
	integer ivar
	integer, optional :: irange

	integer id,irange_local

	id = strings_get_id(name)
	if( id /= 0 ) then
	  write(6,*) ivar,'  ',name
	  write(6,*) 'old id: ',id
	  write(6,*) pentry(id)%full
	  write(6,*) pentry(id)%short
	  write(6,*) pentry(id)%ivar
	  write(6,*) pentry(id)%irange
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

	do id=1,idlast
	  if( ivar == pentry(id)%ivar ) then
	    pentry(id)%short = short
	  end if
	end do

	!id = strings_get_id(ivar)
	!if( id <= 0 ) return
	!pentry(id)%short = short

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

	return

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
!****************************************************************
!****************************************************************

	subroutine ivar2filename(ivar,filename)

! use this to create unique filename

	use shyfem_strings

	implicit none

	integer ivar
	character*(*) filename

	integer isub,i
	character*80 string
	character*1, parameter :: dir(2) = (/'x','y'/)
	integer, save :: idir = 0

	logical has_direction

	call strings_get_short_name(ivar,string,isub)

	filename = string
	if( filename == ' ' ) then	!could not find ivar
          write(filename,'(i4)') ivar
          filename = adjustl(filename)
	else if( isub > 0 ) then
	  write(string,'(i4)') isub
	  do i=1,4
	    if( string(i:i) == ' ' ) string(i:i) = '0'
	  end do
	  string(1:1) = '_'
	  filename = trim(filename) // string(1:4)
	else if( has_direction(string) ) then
	  idir = mod(idir,2) + 1
	  filename = trim(filename) // '-' // dir(idir)
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

	function string_is_this_short(short,string)

! checks if string reduces to short

	use shyfem_strings

	implicit none

	logical string_is_this_short
	character*(*) short,string

	integer isvar,ivar

	call strings_get_ivar(short,isvar)
	call strings_get_ivar(string,ivar)

	string_is_this_short = ( isvar == ivar )

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

	!write(6,*) trim(s1),'  ',trim(s2),l1,l2,l
	if( l == 0 ) return

	compare_svars = ( s1(1:l) == s2(1:l) )

	end

!****************************************************************

	subroutine string_direction_and_unit(string,dir,unit)

c finds direction if vector

	implicit none

	character(*) string,dir,unit

	integer l,lu

	l = len_trim(string)

	unit = ' '
	if( string(l:l) == ']' ) then		!has unit
	  lu = index(string,'[',back=.true.)
	  if( lu == 0 ) then
	    write(6,*) 'Cannot parse unit: ',trim(string)
	    stop 'error stop string_direction_and_unit'
	  end if
	  unit = string(lu+1:l-1)
	  l = len_trim(string(1:lu-1))	!pop trailing spaces and unit
	end if

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

	has_direction = (
     +			  ivar == 2 
     +		     .or. ivar == 3 
     +		     .or. ivar == 21
     +		     .or. ivar == 42
     +			)

	end

!****************************************************************

	subroutine list_strings

	use shyfem_strings

	implicit none

	integer ivar,isub,iv
	character*40 full
	character*10 short

	call populate_strings

	do ivar=1,1000
	  if( ivar == 30 ) cycle		!old name - do not use
	  call strings_get_full(ivar,full,isub)
	  call strings_get_short(ivar,short,isub)
	  if( isub > 0 ) cycle
	  if( short /= ' ' ) then
	    write(6,*) ivar,'  ',short,'  ',full
	  end if
	end do

	end

!****************************************************************
!****************************************************************
!****************************************************************

	recursive subroutine populate_strings

! populates string information
!
! put prefered name first if more than one names apply

	use shyfem_strings

	implicit none

	logical, save :: bread = .false.

	if( bread ) return
	bread = .true.

!---------------------------------------------------------------------

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
	call strings_add_new('airpressure',20)
	call strings_add_new('pressure',20)
	call strings_add_new('wind velocity',21)
	call strings_add_new('solar radiation',22)
	call strings_add_new('sradiation',22)
	call strings_add_new('air temperature',23)
	call strings_add_new('airtemperature',23)
	call strings_add_new('tair',23)
	call strings_add_new('humidity (relative)',24)
	call strings_add_new('rhumidity',24)
	call strings_add_new('cloud cover',25)
	call strings_add_new('cc',25)
	call strings_add_new('rain',26)
	call strings_add_new('evaporation',27)
	call strings_add_new('wind speed',28)
	call strings_add_new('windspeed',28)
	call strings_add_new('wind direction',29)
	call strings_add_new('winddir',29)

	call strings_add_new('wet bulb temperature',40)
	call strings_add_new('wetbulbt',40)
	call strings_add_new('dew point temperature',41)
	call strings_add_new('dewpointt',41)
	call strings_add_new('wind stress',42)
	call strings_add_new('wstress',42)
	call strings_add_new('mixing ratio',43)
	call strings_add_new('mixrat',43)
	call strings_add_new('humidity (specific)',44)
	call strings_add_new('shumidity',44)

	call strings_add_new('bottom stress',60)
	call strings_add_new('bstress',60)
	!call strings_add_new('velocity in x-direction',61)
	call strings_add_new('index',75)
	call strings_add_new('type',76)
	call strings_add_new('distance',77)
	call strings_add_new('lgr',80)
	call strings_add_new('lagrangian (general)',80)
	call strings_add_new('lagage',81)
	call strings_add_new('lagrangian age',81)
	call strings_add_new('lagdep',82)
	call strings_add_new('lagrangian depth',82)
	call strings_add_new('lagtyp',83)
	call strings_add_new('lagrangian type',83)
	call strings_add_new('lagcus',84)
	call strings_add_new('lagrangian custom',84)
	call strings_add_new('ice cover',85)
	call strings_add_new('time step',95)
	call strings_add_new('timestep',95)
	call strings_add_new('time over threshold',97)
	call strings_add_new('timeot',97)
	call strings_add_new('age',98)
	call strings_add_new('renewal time',99)
	call strings_add_new('residence time',99)
	call strings_add_new('wrt',99)

	call strings_add_new('waves (general)',230)
	call strings_add_new('wave height (significant)',231)
	call strings_add_new('wheight (significant)',231)
	call strings_add_new('wave period (mean)',232)
	call strings_add_new('wperiod (mean)',232)
	call strings_add_new('wave direction',233)
	call strings_add_new('wdirection',233)
	call strings_add_new('wave orbital velocity',234)
	call strings_add_new('worbital velocity',234)
	call strings_add_new('wave peak period',235)
	call strings_add_new('wpeak period',235)

	call strings_add_new('concentration (multi)',300,100)	!new numbering
	call strings_add_new('concentration (multi old)',30,20)

	call strings_add_new('weutro (pelagic)',700,20)
	call strings_add_new('weutrop',700,20)
	call strings_add_new('weutro (sediment)',720,10)
	call strings_add_new('weutrosd',720,10)
	call strings_add_new('weutro (shell fish)',730,10)
	call strings_add_new('weutrosf',730,10)

	call strings_add_new('suspended sediment concentration',800,50)
	call strings_add_new('ssc',800,50)
	call strings_add_new('erosion-deposition',891)
	call strings_add_new('sederodep',891)
	call strings_add_new('grainsize (average)',892)
	call strings_add_new('bottom shear stress',893)
	call strings_add_new('sbstress',893)
	call strings_add_new('mud fraction',894)
	call strings_add_new('mudfrac',894)
	call strings_add_new('bedload transport',895)

	call strings_add_new('suspended sediments',850)
	call strings_add_new('sssc',850)
	call strings_add_new('bed sediments [kg]',851)
	call strings_add_new('bsedkg',851)
	call strings_add_new('bed sediments [kg/m**2]',852)
	call strings_add_new('bsedka',852)
	call strings_add_new('bed sediments [m]',853)
	call strings_add_new('bsedm',853)

	call strings_add_new('var',-9)		!special treatment
	call strings_add_new('ivar',-9)

!---------------------------------------------------------------------
! here short description
!---------------------------------------------------------------------

	call strings_set_short(0,'mass')
	call strings_set_short(1,'zeta')
	call strings_set_short(2,'vel')
	call strings_set_short(3,'transp')
	call strings_set_short(5,'bathy')
	call strings_set_short(6,'speed')
	call strings_set_short(7,'dir')
	call strings_set_short(10,'tracer')
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
	call strings_set_short(23,'tair')
	call strings_set_short(24,'rhum')
	call strings_set_short(25,'cc')
	call strings_set_short(26,'rain')
	call strings_set_short(27,'evap')
	call strings_set_short(28,'windspeed')
	call strings_set_short(29,'winddir')

	call strings_set_short(40,'wetbulbt')
	call strings_set_short(41,'dewpointt')
	call strings_set_short(42,'wstress')
	call strings_set_short(43,'mixrate')
	call strings_set_short(44,'shum')

	call strings_set_short(60,'bstress')
	call strings_set_short(75,'index')
	call strings_set_short(76,'type')
	call strings_set_short(77,'distance')
	call strings_set_short(80,'lgr')
	call strings_set_short(81,'lagage')
	call strings_set_short(82,'lagdep')
	call strings_set_short(83,'lagtyp')
	call strings_set_short(84,'lagcus')
	call strings_set_short(85,'ice')
	call strings_set_short(95,'timestep')
	call strings_set_short(97,'timeot')
	call strings_set_short(98,'age')
	call strings_set_short(99,'wrt')

	call strings_set_short(230,'waves')
	call strings_set_short(231,'wheight')
	call strings_set_short(232,'wper')
	call strings_set_short(233,'wdir')
	call strings_set_short(234,'worb')
	call strings_set_short(235,'wpeak')

	call strings_set_short(300,'conc')
	call strings_set_short(30,'conc')

	call strings_set_short(700,'weutrop')
	call strings_set_short(720,'weutrosd')
	call strings_set_short(730,'weutrosf')

	call strings_set_short(800,'ssc')
	call strings_set_short(891,'sederodep')
	call strings_set_short(892,'grainsize')
	call strings_set_short(893,'sbstress')
	call strings_set_short(894,'mudfrac')
	call strings_set_short(895,'bedload')

	call strings_set_short(850,'sssc')
	call strings_set_short(851,'bsedkg')
	call strings_set_short(852,'bsedka')
	call strings_set_short(853,'bsedm')

!---------------------------------------------------------------------

	call strings_check_consistency

!---------------------------------------------------------------------

	end

!****************************************************************
!****************************************************************
!****************************************************************

	subroutine test_strings

	use shyfem_strings

	implicit none

	integer ivar,isub,iv
	character*40 full
	character*10 short

	call populate_strings

	do ivar=1,1000
	  if( ivar == 30 ) cycle		!old name - do not use
	  call strings_get_full(ivar,full,isub)
	  call strings_get_short(ivar,short,isub)
	  if( isub > 0 ) cycle
	  if( short /= ' ' ) then
	    write(6,*) ivar,isub,short,'  ',full
	  end if
	  if( short == ' ' .neqv. full == ' ' ) then
	    write(6,*) 'not equivalent: ',short,full
	    stop 'error stop'
	  end if
	  if( full == ' ' ) cycle
	  call strings_get_ivar(full,iv)
	  if( iv /= ivar ) then
	    write(6,*) 'inconsistency full: ',ivar,iv,full
	  end if
	  call strings_get_ivar(short,iv)
	  if( iv /= ivar ) then
	    write(6,*) 'inconsistency short: ',ivar,iv,short
	  end if
	end do

	end

!****************************************************************

	subroutine test_meteo_strings

	use shyfem_strings

	implicit none

	logical bshort
	integer ivar,isub,iv,i
	character*40 full
	character*10 short,unit
	character*3 dir
	character*80 name
	character*80 ss(17)

	logical string_is_this_short

	call populate_strings

        ss(1) = 'wind stress - x [N/m**2]'
        ss(2) = 'wind stress - y [N/m**2]'
        ss(3) = 'wind velocity - x [m/s]'
        ss(4) = 'wind velocity - y [m/s]'
        ss(5) = 'wind speed [m/s]'
        ss(6) = 'wind speed [knots]'
        ss(7) = 'wind direction [deg]'

	ss(8) = 'pressure (atmospheric) [Pa]'
        ss(9) = 'pressure (atmospheric) [mbar]'

        ss(10) = 'rain [mm/day]'
        ss(11) = 'ice cover [0-1]'

        ss(12) = 'solar radiation [W/m**2]'
        ss(13) = 'air temperature [C]'
        ss(14) = 'humidity (relative) [%]'
        ss(15) = 'cloud cover [0-1]'
        ss(16) = 'wet bulb temperature [C]'
        ss(17) = 'dew point temperature [C]'

	do i=1,17
	  name=ss(i)
	  call strings_get_short_name(name,short)
	  if( short == ' ' ) short = '   ****   '
	  call string_direction_and_unit(name,dir,unit)
	  bshort = string_is_this_short(short,name)
	  write(6,*) short,'  ',dir,unit,' ',bshort,'   ',trim(name)
	end do

	end

!****************************************************************

!	program test_strings_main
!	call test_strings
!	call test_meteo_strings
!	end

!****************************************************************

