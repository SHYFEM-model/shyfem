
!--------------------------------------------------------------------------
!
!    Copyright (C) 1997-1998,2000,2003,2005,2008-2019  Georg Umgiesser
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

c namelist read routines
c
c contents :
c
c revision log :
c
c 01.06.1997	ggu	restructured (localizing nrd functions)
c 17.06.1997	ggu	nlsh, nlsa out of file
c 07.05.1998	ggu	nrdveci, nrdvecr return -1 on error
c 06.08.1998	ggu	nrdtst changed due to compiler warning of AIX
c 01.02.2000	ggu	new routine setsec
c 06.11.2000	ggu	bug found in reading sections ($NRDLIN)
c 05.08.2003	ggu	accept string enclosed in '..' and ".."
c 11.03.2005	ggu	work in double precision in nrdnum()
c 07.11.2005	ggu	better debugging output in nrdpar
c 28.04.2008	ggu	all routines changed to double precision
c 03.09.2008	ggu	nrdvecr slightly changed return value
c 02.12.2008	ggu	bug in nrdnum: kexp was double precision
c 09.03.2009	ggu	bug in nrdsec: use local name to manipolate string
c 26.08.2009	ggu	allow '_' for names (USE_)
c 23.03.2010	ggu	changed v6.1.1
c 09.04.2010	ggu	changed v6.1.3
c 28.09.2010	ggu	changed VERS_6_1_11
c 09.12.2011	ggu	changed VERS_6_1_38
c 30.03.2012	ggu	changed VERS_6_1_51
c 26.06.2012	ggu	changed VERS_6_1_55
c 17.12.2012	ggu	changed VERS_6_1_61a
c 27.02.2013	ggu	handle extra information on section line
c 03.05.2013	ggu	changed VERS_6_1_63
c 20.01.2014	ggu	new routine nrdtable()
c 28.01.2014	ggu	changed VERS_6_1_71
c 26.11.2014	ggu	changed VERS_7_0_7
c 08.01.2015	ggu	new version for nrdvec*()
c 23.01.2015	ggu	changed VERS_7_1_4
c 05.02.2015	ggu	program completely rewritten (modules introduced)
c 08.02.2015	ggu	accept also '!' and '#' for end comment on line
c 26.02.2015	ggu	changed VERS_7_1_5
c 12.05.2015	ggu	new char table
c 21.05.2015	ggu	changed VERS_7_1_11
c 05.06.2015	ggu	changed VERS_7_1_12
c 17.07.2015	ggu	changed VERS_7_1_52
c 01.02.2016	ggu	bug in nls_insert_variable() -> new char variable
c 19.02.2016	ggu	changed VERS_7_5_2
c 25.05.2016	ggu	changed VERS_7_5_10
c 31.03.2017	ggu	changed VERS_7_5_24
c 13.04.2017	ggu	changed VERS_7_5_25
c 26.10.2017	ggu	new isctable read (multiple numbers + description)
c 04.11.2017	ggu	changed VERS_7_5_34
c 14.11.2017	ggu	changed VERS_7_5_36
c 03.04.2018	ggu	changed VERS_7_5_43
c 13.05.2018	ggu	bug fix in nls_copy_isctable()
c 06.07.2018	ggu	changed VERS_7_5_48
c 18.12.2018	ggu	changed VERS_7_5_52
c 16.02.2019	ggu	changed VERS_7_5_60
c 12.03.2019	ggu	before namelist read clean arrays
c
c notes :
c
c structure of calls ----------------------------------
c
c******************************************************************
c******************************************************************
c******************************************************************

!==================================================================
        module nls
!==================================================================

        implicit none

	integer, save, private :: unit = 0
	integer, save, private :: ioff = 0
	integer, save, private :: nline = 0
	integer, save, private :: length = 80
	character*80, save, private :: line = ' '

	integer, save, private :: nlsdim = 0
	double precision, save, private, allocatable :: nls_val(:)
	integer, save, private, allocatable :: nls_table(:,:)
	character*80, save, private, allocatable :: nls_string(:)

	integer, save, private :: snum = 0
	logical, save, private :: sread = .false.
	character*80, save, private :: sname = ' '

	character*80, save, private :: old_name = ' '	!last name read

!==================================================================
        contains
!==================================================================

	subroutine nls_set_section(name,num)

c memorizes section name

	character*(*) name	!section name
	integer num		!number of section

	sname = name
	snum = num
	sread = .false.

	end subroutine nls_set_section

c******************************************

	subroutine nls_get_section(name,num)

c gets section name

	character*(*) name	!section name
	integer num		!number of section

	name = sname
	num = snum 
	ioff = 0

	end subroutine nls_get_section

c******************************************

	function nls_handle_section(name)

c checks if can handle section name

	logical nls_handle_section	!true if can handle section
	character*(*) name	!section name

	if( name .eq. sname ) then
	  sread = .true.
	  nls_handle_section = .true.
	else
	  nls_handle_section = .false.
	end if

	end function nls_handle_section

c******************************************

	function nls_has_read_section()

c actual section has been read ?

	logical nls_has_read_section	!true if section has been read

	nls_has_read_section = sread

	end function nls_has_read_section

c******************************************************************
c******************************************************************
c******************************************************************

	subroutine nls_alloc

	double precision, allocatable :: daux(:)
	integer, allocatable :: iaux(:,:)
	character*80, allocatable :: saux(:)

	if( nlsdim == 0 ) then
	  nlsdim = 10
          allocate(nls_val(nlsdim))
          allocate(nls_table(2,nlsdim))
          allocate(nls_string(nlsdim))
          return
        else
          nlsdim = nlsdim*2
          allocate(daux(nlsdim))
          allocate(iaux(2,nlsdim))
          allocate(saux(nlsdim))
          daux(1:nlsdim/2) = nls_val(1:nlsdim/2)
          iaux(:,1:nlsdim/2) = nls_table(:,1:nlsdim/2)
          saux(1:nlsdim/2) = nls_string(1:nlsdim/2)
          call move_alloc(daux,nls_val)
          call move_alloc(iaux,nls_table)
          call move_alloc(saux,nls_string)
        end if

	end subroutine nls_alloc

c******************************************************************

	subroutine nls_test_alloc(n)

	integer n

	if( n > nlsdim ) call nls_alloc

	end subroutine nls_test_alloc

c******************************************************************

	subroutine nls_init_section

!	empty for now

	end subroutine nls_init_section

c******************************************************************

	subroutine nls_finish_section

	if( nlsdim > 0 ) then
	  nlsdim = 0
          deallocate(nls_val)
          deallocate(nls_table)
          deallocate(nls_string)
	end if

	end subroutine nls_finish_section

c******************************************************************

	subroutine nls_return_line(aline,iline)

	character*(*) aline
	integer iline

	aline = line
	iline = nline

	end subroutine nls_return_line

c******************************************************************
c******************************************************************
c******************************************************************

	function nls_is_nls_file(file)

	logical nls_is_nls_file
	character*(*) file

	integer nsect,iunit
	integer ifileo
	integer num
	character*80 section		!section name
	character*80 extra		!extra information

	nls_is_nls_file = .false.
	nsect = 0

	iunit = ifileo(0,file,'form','old')
	if( iunit <= 0 ) return

	call nls_init(iunit)

	do while( nls_next_section(section,num,extra) )
	  call nls_skip_section
	  nsect = nsect + 1
	end do

	close(iunit)

	nls_is_nls_file = ( nsect > 0 )

	end function nls_is_nls_file

c******************************************************************
c******************************************************************
c******************************************************************

	subroutine nls_init(iunit)

c initializes unit number for name list read

	integer iunit

	unit = iunit
	nline = 0
	ioff = 0

	end subroutine nls_init

c******************************************************************

	function nls_next_line(bcomma)

c reads next line with some information on it
c
c .true. if new line found, else .false.

	logical nls_next_line
	logical bcomma

	integer ios
	!character*80 linaux,name
	character*80 name
	character*128 linaux,lintab

	integer ichafs

	nls_next_line = .false.

	do
	  read(unit,'(a)',iostat=ios) linaux
	  if( ios .gt. 0 ) goto 98			!read error
	  if( ios .lt. 0 ) return			!end of file

	  nline = nline + 1
	  ioff = 1
	  call tablnc(linaux,lintab)
	  call nls_check_line_length(lintab)
	  line = lintab(1:80)
	  !write(6,*) trim(line)
	  if( nls_skip_whitespace_on_line(bcomma) ) exit  !something on line
	end do

	nls_next_line = .true.

	!write(6,*) trim(line)
	!write(6,*) 'nls_next_line: found',ioff,line(ioff:ioff)
	return
   98	continue
	write(6,*) 'error reading from unit: ',unit
	call filna(unit,name)
	write(6,*) 'file name: ',trim(name)
	write(6,*) 'iostat error = ',ios
	write(6,*) 'line number  = ',nline
	stop 'error stop nls_next_line: read error'
	end function nls_next_line

c******************************************************************

	function nls_next_data_line(lineout,bcomma)

c reads next data line and returns it

	logical nls_next_data_line
	character*(*) lineout
	logical bcomma

	character*10 sect

	sect = ' '
	nls_next_data_line = .true.

	if( nls_skip_whitespace(bcomma) ) then
	  if( nls_is_section(sect) ) then
	    if( sect .ne. 'end' ) then
	      write(6,*) 'new section found: ',sect
	      stop 'error stop nls_next_data_line: new section'
	    end if
	    nls_next_data_line = .false.
	  end if
	else
	  nls_next_data_line = .false.
	end if

	lineout = line
	if( sect /= 'end' ) ioff = len_trim(line) + 1

	end function nls_next_data_line

c******************************************************************

	subroutine nls_check_line_length(line)

	implicit none

	character*(*) line

	if( len_trim(line) > 80 ) then
	  write(6,*) 'line longer than 80 chars: ',len_trim(line)
	  write(6,*) trim(line)
	  stop 'error stop nls_check_line_length: line too long'
	end if

	end

c******************************************************************

	function nls_skip_whitespace_on_line(bcomma)

c true if something found on line

	logical nls_skip_whitespace_on_line
	logical bcomma				!also skip commas

	integer l
	character*1 c,ccomma

	ccomma = ' '
	if( bcomma ) ccomma = ','

	nls_skip_whitespace_on_line = .false.
	if( ioff < 1 ) return			!never called

	nls_skip_whitespace_on_line = .true.

	l = len_trim(line)
	do while( ioff <= l )
	  c = line(ioff:ioff)
	  if( c .eq. '!' ) exit					!comment
	  if( c .eq. '#' ) exit					!comment
	  if( c .ne. ' ' .and. c .ne. ccomma ) return
	  ioff = ioff + 1
	end do
	
	nls_skip_whitespace_on_line = .false.

	end function nls_skip_whitespace_on_line

c******************************************************************

	function nls_skip_whitespace(bcomma)

	logical nls_skip_whitespace
	logical bcomma				!also skip commas

	nls_skip_whitespace = .true.

	if( nls_skip_whitespace_on_line(bcomma) ) return
	
	nls_skip_whitespace = nls_next_line(bcomma)

	end function nls_skip_whitespace

c******************************************************************

	subroutine nls_show_line_position(line_opt,ioff_opt)

	character*(*), optional :: line_opt
	integer, optional :: ioff_opt

	character*80 aux
	integer i

	aux = line
	if( present(line_opt) ) aux = line_opt
	write(6,'(a)') trim(aux)

	i = ioff
	if( present(ioff_opt) ) i = ioff_opt
	aux = ' '
	aux(i:i) = '^'
	write(6,'(a)') trim(aux)

	end subroutine nls_show_line_position

c******************************************************************
c******************************************************************
c******************************************************************

	function nls_next_section(section,num,extra)

c finds next section
c
c	finds next section and returns name and number of section
c	if section is not numbered, num = 0
c	if section is found return .true., else .false.

	logical nls_next_section
	integer num
	character*(*) section		!section name
	character*(*) extra		!extra information

c read until '&' or '$' as first non white space char of line found

	section = ' '
	nls_next_section = .false.

	do while( .not. nls_next_section )
	  if( .not. nls_next_line(.false.) ) exit
	  nls_next_section = nls_is_section(section,num,extra)
	end do

	if( nls_next_section ) then
	  if( section == 'end' ) return
	  sname = section
	  snum = num
	  if( .not. nls_next_line(.false.) ) then  !find first item in section
	    write(6,*) 'section found but no content: ',section
	    stop 'error stop nls_next_section: no content'
	  end if
	end if

	end function nls_next_section

c******************************************************************

	function nls_is_section(section,num_opt,extra_opt)

c checks if we are on a new section definition
c
c	returns name and number of section
c	.true. if section is found, else .false.
c	if section is not numbered, num = 0

	logical nls_is_section
	character*(*) section			!section name
	integer, optional :: num_opt		!number of section
	character*(*) , optional ::extra_opt	!extra information

	character*80 name,extra
	character*1 c
	real f(1)
	integer i,ios,ianz,num
	integer istart,iend

	integer nrdvar,iscanf
	logical is_digit

	nls_is_section = .false.
	num = 0
	name = ' '
	section = ' '
	extra = ' '

	c = line(ioff:ioff)
	if( c .ne. '&' .and. c .ne. '$' ) return

c start of section found -> find name and number

	ioff=ioff+1
	istart = ioff
	!write(6,*) ioff,line(ioff:ioff)
	if( .not. nls_read_name(name) ) goto 99
	call to_lower(name)
	iend = len(trim(name))

c name found -> look if there is a number at end of name

	i = iend
	do while( i .gt. 0 .and. is_digit(name(i:i)) )	!number
	  i = i - 1
	end do
	i = i + 1

c if there is a number, strip it and put it in num

	if( i .le. iend ) then			!number to read
	  ianz = iscanf(name(i:iend),f,1)  	!$IREAD
	  if( ianz .ne. 1 ) goto 97
	  num = nint(f(1))
	  name(i:iend) = ' '
	end if

c ok, new section found

	nls_is_section = .true.
	section = name

c now look for extra information

	if( nls_read_name(name) ) extra = name

c copy optional arguments

	if( present(num_opt) ) num_opt = num
	if( present(extra_opt) ) extra_opt = extra

c end of routine

	return
   97	continue
	write(6,*) 'error reading section number in following line'
	write(6,*) line
	stop 'error stop nls_is_section: section number'
   98	continue
	write(6,*) 'error reading from unit ',unit
	stop 'error stop nls_is_section: unit'
   99	continue
	write(6,*) 'error in following line'
	write(6,*) line
	stop 'error stop nls_is_section: read error'
	end function nls_is_section

c******************************************************************

	subroutine nls_skip_section

c skips over data in section

	integer num
	character*10 section,extra,old_section

	old_section = sname(1:10)

	if( nls_next_section(section,num,extra) ) then
	  if( section .eq. 'end' ) return
	  write(6,*) 'new section found: ',section
	  write(6,*) 'old section not ended: ',old_section
	  stop 'error stop nls_skip_section: section not ended'
	else
	  write(6,*) 'EOF found - section not ended: ',sname
	  stop 'error stop nls_skip_section: no end of section'
	end if

	end subroutine nls_skip_section

c******************************************************************

	function nls_next_item(name,value,text)

c returns next item in current section
c
c > 0	type of item found
c	 1 : number variable with name
c	 2 : number variable without name
c	 3 : character variable with name
c	 4 : character variable without name
c = 0	end of section

	integer nls_next_item
	character*(*) name	!name of item read
	double precision value	!value of item if numeric
	character*(*) text	!text of item if string

	integer itype
	character*1 c
	character*10 section

	logical is_letter

	nls_next_item = 0
	name = ' '
	value = 0.
	text = ' '

	if( .not. nls_skip_whitespace(.true.) ) return

	c = line(ioff:ioff)

	if( nls_is_section(section) ) then
	  if( section == 'end' ) return
	  write(6,*) 'new start of section found: ',section
	  write(6,*) 'while still in old section: ',sname
	  stop 'error stop nls_next_item: no end of section'
	else if( is_letter(c) .or. c == '_' ) then
	  call nls_read_assignment(name)
	else if( old_name == ' ' ) then
	  write(6,*) 'No parameter name found in line: '
	  write(6,*) trim(line)
	  stop 'error stop nls_next_item: no parameter name'
	end if

	c = line(ioff:ioff)

	!write(6,*) 'after assignment: ',ioff,c

	if( nls_is_section(section) ) then
	  write(6,*) 'section after assignement found'
	  write(6,*) 'looking for value of parameter ',name
	  stop 'error stop nls_next_item: section instead value found'
	else if( is_letter(c) .or. c == '_' ) then
	  write(6,*) 'parameter name after assignement found'
	  write(6,*) 'looking for value of parameter ',name
	  stop 'error stop nls_next_item: no value found'
	else if( c == '"' .or. c == "'" ) then
	  call nls_read_text(name,text)
	  nls_next_item = 3
	else
	  call nls_read_number(name,value)
	  !write(6,*) 'reading number: ',value
	  nls_next_item = 1
	end if

	if( name == ' ' ) then
	  nls_next_item = nls_next_item + 1
	  name = old_name
	end if

	old_name = name

	!write(6,*) 'end of next_item: ',nls_next_item,name,value

	end function nls_next_item

c******************************************************************

	function nls_read_name(name)

	character*(*) name

	logical nls_read_name
	integer i,itype
	character*1 c

	logical bdummy,is_alpha

	name = ' '
	nls_read_name = .false.

	i = ioff
	if( i < 1 .or. i > length ) return
	do while( i < length )
	  i = i + 1
	  c = line(i:i)
	  if( .not. is_alpha(c) ) exit
	end do

	if( i - ioff < 1 ) return
	nls_read_name = .true.
	name = line(ioff:i-1)
	ioff = i

	bdummy = nls_skip_whitespace_on_line(.false.)

	end function nls_read_name

c******************************************************************

	subroutine nls_read_assignment(name)

	character*(*) name

	name = ' '
	if( .not. nls_read_name(name) ) goto 98

	if( .not. nls_skip_whitespace(.false.) ) goto 99
	if( line(ioff:ioff) .ne. '=' ) goto 99
	ioff = ioff + 1
	if( .not. nls_skip_whitespace(.false.) ) goto 99
	if( line(ioff:ioff) .eq. ',' ) goto 99

	return
   98	continue
	write(6,*) 'cannot read parameter name in line ',nline
	write(6,*) trim(line)
	stop 'error stop nls_read_assignment: no name'
   99	continue
	write(6,*) 'cannot find assignment for parameter: ',trim(name)
	write(6,*) trim(line)
	stop 'error stop nls_read_assignment: no assignement'
	end subroutine nls_read_assignment

c******************************************************************

	subroutine nls_read_text(name,text)

c reads text (must be delimited by with ' ' or " ") 

	character*(*) name
	character*(*) text	!text of item

	integer istos

	text = ' '

	if( istos(line,text,ioff) <= 0 ) then
	  write(6,*) 'Cannot find text for parameter ',name
	  write(6,*) trim(line)
	  write(6,*) 'maybe line is too long for processing (max 80)'
	  write(6,*) 'or the final apostrophe is missing'
	  stop 'error stop nls_read_text: error reading value'
	end if

	end subroutine nls_read_text

c******************************************************************

	subroutine nls_read_number(name,value)

c reads number

	character*(*) name
	double precision value

	integer istod

	value = 0.

	if( istod(line,value,ioff) <= 0 ) then
	  write(6,*) 'Cannot read value for parameter ',name
	  write(6,*) line
	  stop 'error stop nls_read_number: error reading value'
	end if

	end subroutine nls_read_number

c******************************************************************
c******************************************************************
c******************************************************************

	function nls_insert_variable(sect,name,value,text)

c reads and inserts values automatically
c
c does not handle vectors (yet)

	use para

	integer nls_insert_variable
	character*(*) sect
	character*(*) name,text
	double precision value

	integer itspar,iscpar
        character*80 section

	nls_insert_variable = nls_next_item(name,value,text)

	if( nls_insert_variable .le. 0 ) return

	if( name .ne. ' ' ) then
          if( itspar(name) .eq. 0 ) goto 93
          if( iscpar(name,sect) .eq. 0 ) goto 94
	end if

	if( nls_insert_variable == 1 ) then
          call dputpar(name,value)
	else if( nls_insert_variable == 2 ) then
	  if( .not. para_is_array_value(name) ) goto 95
	  call para_put_array_value(name,value)
	else if( nls_insert_variable == 3 ) then
	  call putfnm(name,text)
	else if( nls_insert_variable > 4 ) then
	  write(6,*) 'iwhat = ',nls_insert_variable
	  stop 'error stop nls_insert_variable: internal error (1)'
	end if

	return
   93   continue
        write(6,*) 'no parameter with this name:'
        write(6,*) 'name: ',trim(name)
        write(6,*) 'section: ',trim(sect)
        call check_parameter_values('nrdpar')
        !call parinfo(6)
        stop 'error stop nls_insert_variable: no parameter'
   94   continue
        write(6,*) 'parameter is in wrong section:'
        write(6,*) 'parameter type:  ',nls_insert_variable
        write(6,*) 'parameter name:    ',trim(name)
        write(6,*) 'section: ',trim(sect)
        write(6,*) 'text:    ',trim(text)
        call get_sect_of(name,section)
        write(6,*) 'section for this name found: ',trim(section)
        if( nls_insert_variable .ge. 3 ) call prifnm(6)
        if( nls_insert_variable .le. 2 ) call pripar(6)
        stop 'error stop nls_insert_variable: wrong section'
   95   continue
        write(6,*) 'parameter is not an array'
        write(6,*) 'parameter name:    ',trim(name)
        write(6,*) 'this valiable is defined as a scalar variable'
        write(6,*) 'but is given array values'
        stop 'error stop nls_insert_variable: variable is not an array'
	end function nls_insert_variable

c******************************************************************

	subroutine nls_read_namelist(sect)

c reads parameter section and inserts values automatically
c
c does not handle vectors (yet)

	use para

	character*(*) sect

	character*80 name,text
	double precision value
	integer iwhat

	call nls_init_section
	call para_clean_section(sect)	!this cleans arrays read before...

	do
	  iwhat = nls_insert_variable(sect,name,value,text)
	  !write(6,*) iwhat,sect,name
	  if( iwhat .le. 0 ) exit
	  !if( iwhat == 2 .or. iwhat == 4 ) goto 99
	  if( iwhat == 4 ) goto 99
	end do

	!write(6,*) 'end of namelist: ',line

	call nls_finish_section

	return
   99	continue
	write(6,*) 'cannot insert string array for:'
	write(6,*) name
	write(6,*) sect
	stop 'error stop nls_namelist_read: no array'
	end subroutine nls_read_namelist

c******************************************************************

	function nls_read_vector()

c reads number section, stores numbers in internal array
c
c returns total number of values read
c returns -1 in case of dimension or read error

	integer nls_read_vector		!total number of elements read

	integer n
	double precision value
	character*10 sect

	nls_read_vector = -1
	n = 0

	do
	  call nls_read_number(sname,value)	!stops if read error
	  n=n+1
	  if(n.gt.nlsdim) call nls_alloc
	  nls_val(n)=value
	  if( .not. nls_skip_whitespace(.true.) ) goto 98
	  if( nls_is_section(sect) ) then
	    if( sect /= 'end' ) goto 97
	    exit
	  end if
	end do

	nls_read_vector = n

	return
   97	continue
	write(6,*) 'no end of section found: ',sname
	write(6,*) line
	return
   98	continue
	write(6,*) 'read error in following line'
	write(6,*) line
	return
	end function nls_read_vector

c******************************************************************

	subroutine nls_copy_int_vect(n,ivect)

c copies values read from internal storage to vector ivect

	integer n
	integer ivect(n)

	integer i

	do i=1,n
	  ivect(i) = nint(nls_val(i))
	end do

	end subroutine nls_copy_int_vect

c******************************************************************

	subroutine nls_copy_real_vect(n,rvect)

c copies values read from internal storage to vector rvect

	integer n
	real rvect(n)

	integer i

	do i=1,n
	  rvect(i) = nls_val(i)
	end do

	end subroutine nls_copy_real_vect

c******************************************************************

	subroutine nls_copy_char_vect(n,cvect)

c copies values read from internal storage to vector cvect

	integer n
	character*(*) cvect(n)

	integer i

	do i=1,n
	  cvect(i) = nls_string(i)
	end do

	end subroutine nls_copy_char_vect

c******************************************************************

	function nls_read_table()

c reads table in section
c
c every line is read as character
c
c	string1
c	string2
c	etc..
c
c returns total number of lines read
c returns -1 in case of dimension error
c returns -2 in case of read error

	integer nls_read_table	!total number of lines read

	integer n
	logical bend
	character*80 text
	character*10 sect

	bend = .false.
	nls_read_table = -1
	n = 0

	do while( nls_next_data_line(text,.true.) )
	  !write(6,*) '*** ',trim(text)
	  if( nls_is_section(sect) ) exit
	  n=n+1
	  if(n.gt.nlsdim) call nls_alloc
	  nls_string(n)=text
	end do

	if( nls_is_section(sect) ) then
	  bend = ( sect == 'end' )
	  if( .not. bend ) goto 97
	end if

	nls_read_table = n

c--------------------------------------------------------
c end of routine
c--------------------------------------------------------

	return
   97	continue
	write(6,*) 'no end of section found: ',sname
	write(6,*) 'section line read: ',sect
	write(6,*) line
	stop 'error stop nls_read_table'
	end function nls_read_table

c******************************************************************

	function nls_read_ictable()

c reads table in section
c
c table must have following structure (empty lines are allowed)
c
c	ivalue1 'string1'
c	ivalue2 'string2'
c	etc..
c
c an alternative is also possible, where only numbers are given
c
c	ivalue1 ivalue2 ...
c	etc..
c
c in this case the variable string gets a default value
c
c returns total number of values read
c returns -1 in case of dimension error
c returns -2 in case of read error


	integer nls_read_ictable	!total number of elements read

	integer n
	integer itable
	double precision value
	character*80 text
	character*10 sect
	character*1 c

c--------------------------------------------------------
c find out what table it is - 1: only numbers  2:number with text
c--------------------------------------------------------

	call nls_read_number(sname,value)
	if( .not. nls_skip_whitespace_on_line(.true.) ) then
	  itable = 1		!only numbers
	else
	  c = line(ioff:ioff)
	  if( c == '"' .or. c == "'" ) then
	    itable = 2
	  else
	    itable = 1
	  end if
	end if

	ioff = 1
	if( .not. nls_skip_whitespace_on_line(.true.) ) goto 95

c--------------------------------------------------------
c now read the table
c--------------------------------------------------------

	nls_read_ictable = -1
	n = 0

	do
	  call nls_read_number(sname,value)
	  if( itable == 2 ) then
	    if( .not. nls_skip_whitespace_on_line(.true.) ) goto 95
	    call nls_read_text(sname,text)
	    if( nls_skip_whitespace_on_line(.true.) ) goto 95
	  else
	    text = ' '			!default text (empty)
	  end if
	  n=n+1
	  if(n.gt.nlsdim) call nls_alloc
	  nls_val(n)=value
	  nls_string(n)=text

	  if( .not. nls_skip_whitespace(.true.) ) goto 95
	  if( nls_is_section(sect) ) then
	    if( sect /= 'end' ) goto 97
	    exit
	  end if
	end do

	nls_read_ictable = n

c--------------------------------------------------------
c end of routine
c--------------------------------------------------------

	return
   95	continue
	write(6,*) 'in nls_read_ictable reading section: ',sname
	write(6,*) 'error reading table in following line'
	write(6,*) line
	return
   97	continue
	write(6,*) 'no end of section found: ',sname
	write(6,*) line
	return
	end function nls_read_ictable

c******************************************************************

	subroutine nls_copy_ictable(n,ivect,cvect)

c copies values read from internal storage to vector rvect

	integer n
	integer ivect(n)
	character*(*) cvect(n)

	call nls_copy_int_vect(n,ivect)
	call nls_copy_char_vect(n,cvect)

	end subroutine nls_copy_ictable

c******************************************************************

	subroutine nls_read_isctable(n,ns)

c reads general table with more numerical values per section
c
c table must have following structure (empty lines are allowed)
c
c	iv11 iv12 ... 'string1'
c	iv21 iv22 iv23 ... 'string2'
c	etc..
c
c an alternative is also possible, where only numbers are given
c
c	iv11 iv12 0 iv21 iv22 iv23 0 0 iv31
c	etc..
c
c here the zeros separate one section from the next
c in this case string gets a default value
c
c for both formats, line breaks can be at any point
c

	integer n,ns

	logical bsect,bzero,bstring,bdebug
	logical binsertzero
	integer is,itype
	integer n1,n2,nn,i
	double precision value
	character*80 name,text,sect

c--------------------------------------------------------
c initialize before looping
c--------------------------------------------------------

	binsertzero = .true.	!insert zero to divide sections

	is = 0
	n = 0
	bzero = .false.
	bstring = .false.
	bsect = .false.
	bdebug = .true.
	bdebug = .false.

c--------------------------------------------------------
c loop on input
c--------------------------------------------------------

	do
	  itype = nls_next_item(name,value,text)
	  if( itype <= 0 ) exit
	  if( itype == 2 ) then		!value
	    if( value > 0 ) then
	      n = n + 1
	      call nls_test_alloc(n)
	      if( .not. bsect ) then	!create new section
		if( binsertzero ) then
	          nls_val(n) = 0
	          n = n + 1
	          call nls_test_alloc(n)
		end if
	        is = is + 1
	        call nls_test_alloc(is)
	        bsect = .true.
	        nls_table(1,is) = n
	        nls_table(2,is) = 0
	        nls_string(is) = ' '
	      end if
	      nls_val(n) = value
	    else
	      bzero = .true.
	      if( bsect ) then		!zero - close section
	        bsect = .false.
	        nls_table(2,is) = n
	      end if
	    end if
	  else if( itype == 4 ) then	!string - close section
	    bstring = .true.
	    if( .not. bsect ) goto 92
	    bsect = .false.
	    nls_table(2,is) = n
	    nls_string(is) = text
	  else
	    goto 91
	  end if
	end do

c--------------------------------------------------------
c clean up stuff
c--------------------------------------------------------

        if( nls_is_section(sect) ) then
          if( sect /= 'end' ) goto 97
	end if

	if( bzero .and. bstring ) goto 93
	if( bsect ) nls_table(2,is) = n	!still have to close section
	if( binsertzero ) then
          n = n + 1
          call nls_test_alloc(n)
          nls_val(n) = 0
	end if
	ns = is

c--------------------------------------------------------
c debug
c--------------------------------------------------------

	if( .not. bdebug ) return

	write(6,*) 'section table:'
	do is=1,ns
	  n1 = nls_table(1,is)
	  n2 = nls_table(2,is)
	  nn = n2 - n1 + 1
	  write(6,*) is,nn,n1,n2,'  ',trim(nls_string(is))
	end do
	write(6,*) 'section nodes:'
	do is=1,ns
	  n1 = nls_table(1,is)
	  n2 = nls_table(2,is)
	  nn = n2 - n1 + 1
	  write(6,*) is,(nint(nls_val(i)),i=n1,n2)
	end do
	write(6,*) 'all nodes: ',n
	write(6,*) (nint(nls_val(i)),i=1,n)

c--------------------------------------------------------
c end of routine
c--------------------------------------------------------

	return
   97	continue
	write(6,*) 'no end of section found: ',sect
	write(6,*) line
	stop 'error stop nls_read_isctable: no end section'
	return
   94	continue
	write(6,*) 'strings and zeros used together'
	write(6,*) 'line: ',trim(line)
	stop 'error stop nls_read_isctable: zeros and strings'
   93	continue
	write(6,*) 'strings and zeros used together'
	write(6,*) 'cannot use both together'
	stop 'error stop nls_read_isctable: zeros and strings'
   92	continue
	if( bzero .and. bstring ) goto 94	!make err mes clearer
	write(6,*) 'no section open for string'
	write(6,*) 'node numbers must preceed string'
	write(6,*) 'line: ',trim(line)
	stop 'error stop nls_read_isctable: string out of section'
   91	continue
	write(6,*) 'section cannot contain names'
	write(6,*) 'line: ',trim(line)
	stop 'error stop nls_read_isctable: names on line'
	end subroutine nls_read_isctable

c******************************************************************

	subroutine nls_copy_isctable(n,ns,ivect,itable,cvect)

c copies values read from internal storage to vector rvect

	integer n,ns
	integer ivect(n)
	integer itable(2,ns)
	character*(*) cvect(ns)

	call nls_copy_int_vect(n,ivect)
	call nls_copy_char_vect(ns,cvect)
	itable(:,1:ns) = nls_table(:,1:ns)

	end subroutine nls_copy_isctable

c******************************************************************

!==================================================================
        end module nls
!==================================================================

c******************************************************************
c******************************************************************
c******************************************************************


c******************************************************************
c******************************************************************
c******************************************************************
c compatibility routines
c******************************************************************
c******************************************************************
c******************************************************************

	function nrdveci(ivect,ndim)

c reads integer vector in section (compatibility)

	use nls

	implicit none

	integer nrdveci		!total number of elements read
	integer ndim		!dimension of vector
	integer ivect(ndim)	!vector

	integer n

	n = nls_read_vector()
	if( n > ndim) n = -n			!flag dimension error
	if( n > 0 ) call nls_copy_int_vect(n,ivect)
	nrdveci = n

	end

c******************************************************************

	function nrdvecr(rvect,ndim)

c reads real vector in section (compatibility)

	use nls

	implicit none

	integer nrdvecr		!total number of elements read
	integer ndim		!dimension of vector
	real rvect(ndim)	!vector

	integer n

	n = nls_read_vector()
	if( n > ndim) n = -n			!flag dimension error
	if( n > 0 ) call nls_copy_real_vect(n,rvect)
	nrdvecr = n

	end

c******************************************************************

	function nrdictable(ivect,cvect,ndim)

c reads table in section

	use nls

	implicit none

	integer nrdictable		!total number of elements read
	integer ndim			!dimension of vector
	integer ivect(ndim)		!vector
	character*80 cvect(ndim)	!vector

	integer n

	n = nls_read_ictable()
	if( n > ndim) n = -n			!flag dimension error
	if( n > 0 ) call nls_copy_ictable(n,ivect,cvect)
	nrdictable = n

	end

c******************************************************************

	subroutine setsec(name,num)
	use nls
	implicit none
	character*(*) name
	integer num
	call nls_set_section(name,num)
	end

	function handlesec(name)
	use nls
	implicit none
	logical handlesec
	character*(*) name
	handlesec = nls_handle_section(name)
	end

	function hasreadsec()
	use nls
	implicit none
	logical hasreadsec
	hasreadsec = nls_has_read_section()
	end

	subroutine nrdini(iunit)
	use nls
	implicit none
	integer iunit
	call nls_init(iunit)
	end

	function nrdsec(section,num,extra)
	use nls
	implicit none
	integer nrdsec
	character*(*) section
	integer num
	character*(*) extra
	nrdsec = 0
	if( nls_next_section(section,num,extra) ) nrdsec = 1
	end

	subroutine nrdskp
	use nls
	implicit none
	call nls_skip_section
	end

	function nrdlin(line)
	use nls
	implicit none
	integer nrdlin
	character*(*) line
	nrdlin = 0
	if( nls_next_data_line(line,.false.) ) nrdlin = 1
	end

c***********************

        subroutine nrdins(sect)
	use nls
	implicit none
	character*(*) sect
	call nls_read_namelist(sect)
	end

	function nrdpar(sect,name,value,text)
	use nls
	implicit none
	integer nrdpar
	character*(*) sect,name,text
	double precision value
	nrdpar = nls_insert_variable(sect,name,value,text)
	end

	function nrdnxt(name,value,text)
	use nls
	implicit none
	integer nrdnxt
	character*(*) name,text
	double precision value
	nrdnxt = nls_next_item(name,value,text)
	end

c******************************************************************
c******************************************************************
c******************************************************************

	subroutine nrdtst

c subroutine to test nrd... routines
c
c to use write main as follows :
c
c--------------------------
c	call nrdtst
c	end
c--------------------------

	implicit none


	end

c*******************************************************

	subroutine petras_read

	integer iunit,iw
	integer nrdnxt
	character*80 name,text
	double precision dvalue

	iunit = 5
	call nrdini(iunit)

	do
	  iw = nrdnxt(name,dvalue,text)
	  if( iw .le. 0 ) exit
	  write(6,*) 'name = ',trim(name),' value = ', dvalue
	end do

	end

c*******************************************************

c	program nls_test_main
c	call petras_read
c	end

c*******************************************************

