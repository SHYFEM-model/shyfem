
!--------------------------------------------------------------------------
!
!    Copyright (C) 2016-2020  Georg Umgiesser
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
! 10.09.2019	ggu	condense all strings
! 24.10.2019	ggu	new string grainsizep
! 28.01.2020	ggu	new string vorticity
! 05.04.2020	ggu	review of strings
! 11.04.2020	ggu	review of directional strings
! 17.04.2020	ggu	new routine strings_meteo_convention()
! 18.05.2020	ggu	do not attach direction to file name
! 22.04.2021	ggu	new routines checking and populating strings
!
! contents :
!
! name		any abbreviation of variable name
! full		full name (canonical)
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
!	subroutine strings_get_canonical(name,full)
!
!       subroutine strings_add_new(name,ivar,irange)
!       subroutine strings_set_short(ivar,short)
!
!       subroutine strings_check_consistency
!       subroutine strings_info(id)
!
!-------------------------------------------------------
!
!       subroutine string2ivar(string,ivar)
!       subroutine ivar2string(ivar,string,isub)
!	subroutine ivar2filename(ivar,filename)
!	subroutine ivar2femstring(ivar,femstring)
!
!       function string_is_this_short(short,string)
!
!       subroutine compress_string(string)
!       function compare_svars(s1,s2)
!
!       subroutine string_direction_and_unit(string,dir,unit)
!       function has_direction(name)
!
!       subroutine get_direction_ivars(ivar,ivars,ivard)
!       subroutine strings_pop_direction(name)
!       subroutine list_strings
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

	  character*80 :: search
	  character*80 :: full
	  character*10 :: short
	  integer :: ivar
	  integer :: irange

	end type entry

        logical, save, private :: bpopulated = .false.    !is populated?
        integer, save, private :: idlast = 0
        integer, save, private :: ndim = 0
        integer, save, private :: ivarmax = 0
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

        pentry(id)%search = ' '
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

	integer id,i,ids
	character(len=len(name)) :: string
        logical bdebug
	logical compare_svars
        character*80 name_aux

	call check_populate

	strings_get_id_by_name = 0
        if( name == ' ' ) return

	ids = 0
	string = name
	call compress_string(string)

	do id=1,idlast
	  if( compare_svars(pentry(id)%search,string) ) exit
	  if( compare_svars(pentry(id)%short,string) ) ids = id
	end do
	if( id > idlast ) id = ids

        name_aux = string
        bdebug = ( name_aux(1:6) == 'transp' )
        bdebug = .false.
        if( bdebug ) then
          write(6,*) '-----------------------------'
          write(6,*) 'debug: ',trim(name),' ',trim(string)
          write(6,*) 'id/ids: ',id,ids
          if( id > 0 )  write(6,*) 'id: ',trim(pentry(id)%search)
          if( ids > 0 ) write(6,*) 'ids: ',trim(pentry(ids)%short)
          write(6,*) '-----------------------------'
        end if

	strings_get_id_by_name = id

	end function strings_get_id_by_name

!******************************************************************

	function strings_get_id_by_ivar(ivar)

! given ivar finds id and sub-range

	integer strings_get_id_by_ivar
	integer ivar

	integer id,ivmin,ivmax
	logical bdebug

	call check_populate

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

	call check_populate

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

	call check_populate

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

	ivar = -1

	id = strings_get_id(name)
	if( id == 0 ) return

	ivar = pentry(id)%ivar

	end subroutine strings_get_ivar

!******************************************************************

	subroutine strings_get_canonical(name,full)

! returns canonical name

        implicit none

	character*(*) name
	character*(*) full

	call strings_get_full_name(name,full)

	end subroutine strings_get_canonical

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine strings_add_new(name,ivar,irange)

	character*(*) name
	integer ivar
	integer, optional :: irange

	integer id,irange_local
	character*80 string

	string = name
	call compress_string(string)
	id = strings_get_id(string)
	if( id /= 0 ) then
	  write(6,*) ivar,'  ',name,'  ',string
	  call strings_info(id)
	  stop 'error stop strings_add_new: name already present'
	end if

	call strings_init_new_id(id)

	irange_local = 0
	if( present(irange) ) irange_local = irange
	
	pentry(id)%search = string
	pentry(id)%full = name
	pentry(id)%ivar = ivar
	pentry(id)%irange = irange_local

        ivarmax = max(ivarmax,ivar)

	end subroutine strings_add_new

!******************************************************************

	subroutine strings_set_short(ivar,short)

	integer ivar
	character*(*) short

	integer id

	do id=1,idlast
	  if( ivar == pentry(id)%ivar ) then
            if( pentry(id)%short /= ' ' ) then
              write(6,*) 'ivar,short: ',ivar,short,' ',pentry(id)%short
              stop 'error stop strings_set_short: already set'
            end if
	    pentry(id)%short = short
	  end if
	end do

	end subroutine strings_set_short

!****************************************************************

        subroutine strings_check_consistency(bfull)

        logical, optional :: bfull

        logical bfullcheck
        integer id,iv,ivar
        character*80 name,short,full

        bfullcheck = .false.
        if( present(bfull) ) bfullcheck = bfull

        do id=1,idlast
          if( pentry(id)%ivar < 0 ) goto 99
          if( pentry(id)%full==' ' ) goto 99
          if( pentry(id)%short==' ' ) goto 99
          if( pentry(id)%ivar < 0 ) cycle
        end do

        do id=1,idlast
          ivar = pentry(id)%ivar
          name = pentry(id)%full
	  call strings_get_ivar(name,iv)
          if( ivar /= iv ) goto 98
          name = pentry(id)%short
	  call strings_get_ivar(name,iv)
          if( ivar /= iv ) goto 98
          if( .not. bfullcheck ) cycle
          if( pentry(id)%full == pentry(id)%short ) then
            full = pentry(id)%full
            short = pentry(id)%short
            write(6,*) 'warning: full==short: ',ivar
     +                          ,trim(short),'  ',trim(full)
          end if
        end do

        return
   98   continue
        call strings_info(id)
        write(6,*) 'ivar,iv: ',ivar,iv
        stop 'error stop strings_check_consistency: ivar/=iv'
   99   continue
        call strings_info(id)
        write(6,*) 'some parameter has not been set'
        stop 'error stop strings_check_consistency: consistency'
        end subroutine strings_check_consistency

!******************************************************************

	subroutine strings_info(id)

	integer id

	if( id == 0 ) return

	write(6,*) 'id:     ',id
	write(6,*) 'search: ',trim(pentry(id)%search)
	write(6,*) 'full:   ',trim(pentry(id)%full)
	write(6,*) 'short:  ',trim(pentry(id)%short)
	write(6,*) 'ivar:   ',pentry(id)%ivar
	write(6,*) 'irange: ',pentry(id)%irange

	end subroutine strings_info

!****************************************************************

	subroutine strings_list(bfull)

        logical, optional :: bfull
        !logical bfull

	integer ivar,isub,iv,id
	character*28 full,cname
	character*10 short

        integer, allocatable :: icount(:)

        call check_populate

        allocate(icount(0:ivarmax))
        icount = 0

        do id=1,idlast
          ivar = pentry(id)%ivar
          icount(ivar) = icount(ivar) + 1
        end do

        write(6,'(4a)') ' ivar   short name   long name'
	do ivar=0,ivarmax
	  !if( ivar == 30 ) cycle		!old name - do not use
	  call strings_get_full(ivar,full,isub)
	  call strings_get_short(ivar,short,isub)
	  call strings_get_canonical(full,cname)
	  if( isub > 0 ) cycle
	  if( short /= ' ' ) then
	    !write(6,1000) ivar,'  ',short,'  ',full,'  ',cname
	    write(6,1000) ivar,'   ',short,'   ',full
 1000       format(i5,6a)
	  end if
	end do

        if( .not. present(bfull) ) return
        if( .not. bfull ) return

	write(6,*) 'variables with more than one name:'
        write(6,*) 'ndim,idlast,ivarmax: ',ndim,idlast,ivarmax

        do ivar=0,ivarmax
          if( icount(ivar) > 1 ) then
	    call strings_get_full(ivar,full,isub)
            write(6,'(i5,i3,2x,a)') ivar,icount(ivar),trim(full)
            do id=1,idlast
              iv = pentry(id)%ivar
              if( iv == ivar ) then
                write(6,*) id,pentry(id)%short,trim(pentry(id)%full)
              end if
            end do
          end if
        end do

	end subroutine strings_list

!****************************************************************

        subroutine srings_get_ivarmax(ivar_max)

        integer ivar_max

        ivar_max = ivarmax

        end subroutine srings_get_ivarmax

!****************************************************************

        function srings_exists_id(id)

        logical srings_exists_id
        integer id

        srings_exists_id = ( id >= 1 .and. id <= idlast )

        end function srings_exists_id

!****************************************************************

        subroutine check_populate

	if( .not. bpopulated ) then
	  write(6,*) 'strings have not been set up'
	  write(6,*) 'call populate_strings() needed'
	  stop 'error stop check_populate: not populated'
	end if

        end subroutine check_populate

!****************************************************************

        subroutine set_populate

	bpopulated = .true.

        end subroutine set_populate

!****************************************************************

        function is_populated()

	logical is_populated

	is_populated = bpopulated

        end function is_populated

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
	  !filename = trim(filename) // '-' // dir(idir)
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

	string_is_this_short = ( isvar /= -1 .and. isvar == ivar )

	end

!****************************************************************
!****************************************************************
!****************************************************************

	subroutine compress_string(string)

! deletes spaces and _

	implicit none

	character*(*) string

	integer lmax,l,ll
	character*1 c
	character*80 s,dir,unit

        call pop_direction_and_unit(string,dir,unit)

	s = ' '
	ll = 0
	lmax = len_trim(string)

	do l=1,lmax
	  c = string(l:l)
	  if( c == ' ' .or. c == '_' ) cycle
	  ll = ll + 1
	  s(ll:ll) = c
	end do

	string = s

	end

!****************************************************************

	function compare_svars(s1,s2)

! compares two strings if they indicate the same variable

	implicit none

	logical compare_svars
	character*(*) s1,s2

	integer l1,l2,l
        logical, parameter :: bmin = .false.

	compare_svars = .false.

	l1 = len_trim(s1)
	l2 = len_trim(s2)
	l = min(l1,l2)

	!write(6,*) trim(s1),'  ',trim(s2),l1,l2,l

        if( bmin ) then
	  if( l == 0 ) return
	  compare_svars = ( s1(1:l) == s2(1:l) )
        else
	  compare_svars = ( s1 == s2 )
        end if

	end

!****************************************************************

	subroutine string_direction_and_unit(string,dir,unit)

c finds direction and unit - does not change string

	character*(*) string,dir,unit
        character*80 aux

        aux = string

	call pop_direction_and_unit(aux,dir,unit)

        end

!****************************************************************

	subroutine pop_direction_and_unit(string,dir,unit)

c pops direction and unit and returns cleaned string

	implicit none

	character*(*) string,dir,unit

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

	dir = ' '
	if( string(l-1:l) == ' x' .or. string(l-1:l) == '-x' ) then
	  dir = 'x'
	else if( string(l-1:l) == ' y' .or. string(l-1:l) == '-y' ) then
	  dir = 'y'
	end if

        if( dir /= ' ' ) then
	  l = len_trim(string(1:l-1))	!pop trailing spaces and dir
          if( string(l:l) == '-' ) then
	    l = len_trim(string(1:l-1))	!pop trailing spaces and "-"
          else if( string(l-2:l) == ' in' ) then
	    l = len_trim(string(1:l-2))	!pop trailing spaces and "in"
          end if
        end if

        if( unit /= ' ' .or. dir /= ' ' ) string(l+1:) = ' '

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
     +			  ivar == 2 		!velocity
     +		     .or. ivar == 3 		!transport
     +		     .or. ivar == 21		!wind
     +		     .or. ivar == 42		!wind stress
     +			)

	end

!****************************************************************

        subroutine strings_meteo_convention(ivar,bmeteo)

! checks if we have to use meteo convention for this variable

	integer ivar
	logical bmeteo

	if( ivar == 21 .or. ivar == 28 .or. ivar == 29 ) then
	  bmeteo = .true.
	else
	  bmeteo = .false.
	end if

	end

!****************************************************************

        subroutine get_direction_ivars(ivar,ivars,ivard)

! returns ivars for speed and dir if variable is directional

	use shyfem_strings

	implicit none

        integer ivar,ivars,ivard

        if( ivar == 2 ) then			!velocity
          ivars = 6
          ivard = 7
        else if( ivar == 3 ) then		!transport
          ivars = 8
          ivard = 9
        else if( ivar == 21 ) then		!wind
          ivars = 28
          ivard = 29
        else if( ivar == 42 ) then		!wind stress
          ivars = 45
          ivard = 46
	else
	  ivars = 0
	  ivard = 0
        end if

        end

!****************************************************************

        subroutine strings_pop_direction(name)

! deletes direction (-x,-y) from name

        implicit none

        character*(*) name

        integer l

        if( name == ' ' ) return
        
        l = len_trim(name)
        if( l < 2 ) return

        if( name(l-1:l) == '-x' .or. name(l-1:l) == '-y' ) then
          name(l-1:l) = '  '
        end if

        end

!****************************************************************

        subroutine list_strings

	use shyfem_strings

	implicit none

        call strings_list(.false.)

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

	if( is_populated() ) return
        call set_populate

!---------------------------------------------------------------------

	call strings_add_new('mass field',0)
	call strings_add_new('water level',1)
	call strings_add_new('level',1)
	call strings_add_new('velocity',2)
	call strings_add_new('transport',3)
	call strings_add_new('bathymetry',5)
	call strings_add_new('depth',5)
	call strings_add_new('current speed',6)
	call strings_add_new('speed',6)
	call strings_add_new('current direction',7)
	call strings_add_new('direction',7)
	call strings_add_new('transport speed',8)
	call strings_add_new('transport direction',9)
	call strings_add_new('generic tracer',10)
	call strings_add_new('salinity',11)
	call strings_add_new('temperature',12)
	call strings_add_new('density',13)
	call strings_add_new('oxygen',15)
	call strings_add_new('discharge',16)
	call strings_add_new('rms velocity',18)
	call strings_add_new('rms speed',18)
	call strings_add_new('current vorticity',19)

	call strings_add_new('atmospheric pressure',20)
	call strings_add_new('air pressure',20)
	call strings_add_new('pressure',20)
	call strings_add_new('pressure (atmospheric)',20)
	call strings_add_new('wind velocity',21)
	call strings_add_new('solar radiation',22)
	call strings_add_new('sradiation',22)
	call strings_add_new('air temperature',23)
	call strings_add_new('tair',23)
	call strings_add_new('humidity (relative)',24)
	call strings_add_new('rhumidity',24)
	call strings_add_new('cloud cover',25)
	call strings_add_new('precipitation',26)
	call strings_add_new('evaporation',27)
	call strings_add_new('wind speed',28)
	call strings_add_new('wind direction',29)

	call strings_add_new('wet bulb temperature',40)
	call strings_add_new('dew point temperature',41)
	call strings_add_new('wind stress',42)
	call strings_add_new('mixing ratio',43)
	call strings_add_new('mixrat',43)
	call strings_add_new('humidity (specific)',44)
	call strings_add_new('shumidity',44)
	call strings_add_new('wind stress modulus',45)
	call strings_add_new('wind stress direction',46)

	call strings_add_new('bottom stress',60)
	call strings_add_new('general index',75)
	call strings_add_new('general type',76)
	call strings_add_new('general distance',77)
	call strings_add_new('lagrangian (general)',80)
	call strings_add_new('lagrangian age',81)
	call strings_add_new('lagrangian depth',82)
	call strings_add_new('lagrangian type',83)
	call strings_add_new('lagrangian custom',84)
	call strings_add_new('ice cover',85)
	call strings_add_new('ice thickness',86)
	call strings_add_new('relaxation time',94)
	call strings_add_new('time step',95)
	call strings_add_new('time over threshold',97)
	call strings_add_new('water age',98)
	call strings_add_new('renewal time',99)
	call strings_add_new('residence time',99)

	call strings_add_new('waves (general)',230)
	call strings_add_new('wave height (significant)',231)
	call strings_add_new('wave period (mean)',232)
	call strings_add_new('wave direction',233)
	call strings_add_new('wave orbital velocity',234)
	call strings_add_new('wave peak period',235)

	call strings_add_new('concentration',300,100)	!new numbering
	!call strings_add_new('concentration (multi old)',30,20)

	call strings_add_new('weutro (pelagic)',700,20)
	call strings_add_new('weutro (sediment)',720,10)
	call strings_add_new('weutro (shell fish)',730,10)

	call strings_add_new('suspended sediment concentration',800,50)
	call strings_add_new('erosion-deposition',891)
	call strings_add_new('grainsize (average)',892)
	call strings_add_new('bottom shear stress',893)
	call strings_add_new('mud fraction',894)
	call strings_add_new('bedload transport',895)
	call strings_add_new('grainsize (percentage)',896)

	call strings_add_new('suspended sediments',850)
	call strings_add_new('bed sediments mass',851)
	call strings_add_new('bed sediments volume',852)
	call strings_add_new('bed sediments height',853)

	!call strings_add_new('var',-9)		!special treatment
	!call strings_add_new('ivar',-9)

!---------------------------------------------------------------------
! here short description
!---------------------------------------------------------------------

	call strings_set_short(0,'mass')
	call strings_set_short(1,'zeta')
	call strings_set_short(2,'vel')
	call strings_set_short(3,'transp')
	call strings_set_short(5,'bathy')
	call strings_set_short(6,'curspeed')
	call strings_set_short(7,'curdir')
	call strings_set_short(8,'transspeed')
	call strings_set_short(9,'transdir')
	call strings_set_short(10,'tracer')
	call strings_set_short(11,'salt')
	call strings_set_short(12,'temp')
	call strings_set_short(13,'rho')
	call strings_set_short(15,'oxy')
	call strings_set_short(16,'disch')
	call strings_set_short(18,'rms')
	call strings_set_short(19,'vorticity')

	call strings_set_short(20,'airp')
	call strings_set_short(21,'wind')
	call strings_set_short(22,'srad')
	call strings_set_short(23,'airt')
	!call strings_set_short(23,'tair')
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
	call strings_set_short(45,'wstressmod')
	call strings_set_short(46,'wstressdir')

	call strings_set_short(60,'bstress')
	call strings_set_short(75,'index')
	call strings_set_short(76,'type')
	call strings_set_short(77,'distance')
	call strings_set_short(80,'lgr')
	call strings_set_short(81,'lagage')
	call strings_set_short(82,'lagdep')
	call strings_set_short(83,'lagtyp')
	call strings_set_short(84,'lagcus')
	call strings_set_short(85,'icecover')
	call strings_set_short(86,'icethick')
	call strings_set_short(94,'relaxtime')
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
	!call strings_set_short(30,'conc')

	call strings_set_short(700,'weutrop')
	call strings_set_short(720,'weutrosd')
	call strings_set_short(730,'weutrosf')

	call strings_set_short(800,'ssc')
	call strings_set_short(891,'sederodep')
	call strings_set_short(892,'grainsize')
	call strings_set_short(893,'sbstress')
	call strings_set_short(894,'mudfrac')
	call strings_set_short(895,'bedload')
	call strings_set_short(896,'grainsizep')

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
        integer ivar_max
	character*40 full
	character*10 short

	call populate_strings

        call srings_get_ivarmax(ivar_max)

	do ivar=0,ivar_max
	  if( ivar == 30 ) cycle		!old name - do not use
	  call strings_get_full(ivar,full,isub)
	  call strings_get_short(ivar,short,isub)
	  if( short == ' ' .and. full == ' ' ) cycle
	  if( isub > 0 ) cycle
	  if( short /= ' ' ) then
	    write(6,*) ivar,isub,short,'  ',full
	  end if
	  if( short == ' ' .or. full == ' ' ) then
	    write(6,*) 'not equivalent: ',short,full
	    stop 'error stop'
	  end if
	  call strings_get_ivar(full,iv)
	  if( iv /= ivar ) then
	    write(6,*) 'inconsistency full: ',ivar,iv,full
	  end if
	  call strings_get_ivar(short,iv)
	  if( iv /= ivar ) then
	    write(6,*) 'inconsistency short: ',ivar,iv,short
	  end if
          if( full == short ) then
            write(6,*) 'warning: short = full: ',short,full
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

	integer, parameter :: ndim = 22
	character*80 ss(ndim)

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

        ss(18) = 'rain'
        ss(19) = 'shumidity'
        ss(20) = 'cloud cover'

        ss(21) = 'wind stress in x [N/m**2]'
        ss(22) = 'wind stress in y [N/m**2]'

	do i=1,ndim
	  name=ss(i)
	  call strings_get_short_name(name,short)
	  if( short == ' ' ) short = '   ****   '
	  call string_direction_and_unit(name,dir,unit)
	  bshort = string_is_this_short(short,name)
	  write(6,*) short,'  ',dir,unit,' ',bshort,'   ',trim(name)
	end do

	end

!****************************************************************

	subroutine test_strings_all
	use shyfem_strings
	call populate_strings
	call test_strings
        call strings_check_consistency(.true.)
        call strings_list(.true.)
	call test_meteo_strings
	end

!****************************************************************

!	program test_strings_main
!	call test_strings_all
!	end

!****************************************************************

