c
c $Id: subnls.f,v 1.26 2010-02-16 16:21:37 georg Exp $
c
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
c 27.02.2013	ggu	handle extra information on section line
c 20.01.2014	ggu	new routine nrdtable()
c 08.01.2015	ggu	new version for nrdvec*()
c 05.02.2015	ggu	program completely rewritten (modules introduced)
c 08.02.2015	ggu	accept also '!' and '#' for end comment on line
c 12.05.2015	ggu	new char table
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
	character*80, allocatable :: saux(:)

	if( nlsdim == 0 ) then
	  nlsdim = 10
          allocate(nls_val(nlsdim))
          allocate(nls_string(nlsdim))
          return
        else
          nlsdim = nlsdim*2
          allocate(daux(nlsdim))
          allocate(saux(nlsdim))
          daux(1:nlsdim/2) = nls_val(1:nlsdim/2)
          saux(1:nlsdim/2) = nls_string(1:nlsdim/2)
          call move_alloc(daux,nls_val)
          call move_alloc(saux,nls_string)
        end if

	end subroutine nls_alloc

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
	character*80 linaux

	integer ichafs

	nls_next_line = .false.

	do
	  read(unit,'(a)',iostat=ios) linaux
	  if( ios .gt. 0 ) goto 98			!read error
	  if( ios .lt. 0 ) return			!end of file

	  nline = nline + 1
	  ioff = 1
	  call tablnc(linaux,line)
	  !write(6,*) trim(line)
	  if( nls_skip_whitespace_on_line(bcomma) ) exit  !something on line
	end do

	nls_next_line = .true.

	!write(6,*) 'nls_next_line: found',ioff,line(ioff:ioff)
	return
   98	continue
	write(6,*) 'error reading from unit ',unit
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
c	if section is found, nrdsec = 1, else nrdsec = 0

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
	real f(5)
	integer i,ios,ianz,num
	integer istart,iend

	integer nrdvar,itypch,iscan

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
	call uplow(name,'low')
	iend = len(trim(name))

c name found -> look if there is a number at end of name

	i = iend
	do while( i .gt. 0 .and. itypch(name(i:i)) .eq. 1 )	!number
	  i = i - 1
	end do
	i = i + 1

c if there is a number, strip it and put it in num

	if( i .le. iend ) then			!number to read
	  ianz = iscan(name(i:iend),1,f)  	!$IREAD
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

	integer itypch

	nls_next_item = 0
	name = ' '
	value = 0.
	text = ' '

	if( .not. nls_skip_whitespace(.true.) ) return

	c = line(ioff:ioff)
	itype = itypch(c)

	if( nls_is_section(section) ) then
	  if( section == 'end' ) return
	  write(6,*) 'new start of section found: ',section
	  write(6,*) 'while still in old section: ',sname
	  stop 'error stop nls_next_item: no end of section'
	else if( itype == 2 .or. c == '_' ) then
	  call nls_read_assignment(name)
	else if( old_name == ' ' ) then
	  write(6,*) 'No parameter name found in line: '
	  write(6,*) trim(line)
	  stop 'error stop nls_next_item: no parameter name'
	end if

	c = line(ioff:ioff)
	itype = itypch(c)

	!write(6,*) 'after assignment: ',ioff,c

	if( nls_is_section(section) ) then
	  write(6,*) 'section after assignement found'
	  write(6,*) 'looking for value of parameter ',name
	  stop 'error stop nls_next_item: section instead value found'
	else if( itype == 2 .or. c == '_' ) then
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

	logical bdummy
	integer itypch

	name = ' '
	nls_read_name = .false.

	i = ioff
	if( i < 1 .or. i > length ) return
	do while( i < length )
	  i = i + 1
	  c = line(i:i)
	  itype = itypch(c)
	  if( itype /= 1 .and. itype /= 2 .and. c /= '_' ) exit
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

c reads text (must start with ' or ")

	character*(*) name
	character*(*) text	!text of item

	integer istos

	text = ' '

	if( istos(line,text,ioff) <= 0 ) then
	  write(6,*) 'Cannot find text for parameter ',name
	  write(6,*) line(ioff:)
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

	integer nls_insert_variable
	character*(*) sect
	character*(*) name,text
	double precision value

	integer itspar,iscpar

	nls_insert_variable = nls_next_item(name,value,text)

	if( nls_insert_variable .le. 0 ) return

	if( name .ne. ' ' ) then
          if( itspar(name) .eq. 0 ) goto 93
          if( iscpar(name,sect) .eq. 0 ) goto 94
	end if

	if( nls_insert_variable == 1 ) then
          call dputpar(name,value)
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
        write(6,*) 'parameter name:    ',name(1:len_trim(name))
        write(6,*) 'section: ',sect(1:len_trim(sect))
        write(6,*) 'text:    ',text(1:len_trim(text))
        call get_sect_of(name,sect)
        write(6,*) 'section found: ',sect(1:len_trim(sect))
        if( nls_insert_variable .ge. 3 ) call prifnm(6)
        if( nls_insert_variable .le. 2 ) call pripar(6)
        stop 'error stop nls_insert_variable: wrong section'
	end function nls_insert_variable

c******************************************************************

	subroutine nls_read_namelist(sect)

c reads parameter section and inserts values automatically
c
c does not handle vectors (yet)

	character*(*) sect

	character*80 name,text
	double precision value
	integer iwhat

	do
	  iwhat = nls_insert_variable(sect,name,value,text)
	  !write(6,*) iwhat,sect,name
	  if( iwhat .le. 0 ) exit
	  if( iwhat == 2 .or. iwhat == 4 ) goto 99
	end do

	!write(6,*) 'end of namelist: ',line

	return
   99	continue
	write(6,*) 'cannot insert array for:'
	write(6,*) name
	write(6,*) sect
	stop 'error stop nls_namelist_read: no array read yet'
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
	  call nls_read_number(sname,value)
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
	character*80 text
	character*10 sect

	nls_read_table = -1
	n = 0

	do while( nls_next_data_line(text,.true.) )
	  if( nls_is_section(sect) ) then
	    if( sect /= 'end' ) goto 97
	    exit
	  end if

	  n=n+1
	  if(n.gt.nlsdim) call nls_alloc
	  nls_string(n)=text
	end do

	nls_read_table = n

c--------------------------------------------------------
c end of routine
c--------------------------------------------------------

	return
   95	continue
	write(6,*) 'in nls_read_table reading section: ',sname
	write(6,*) 'error reading table in following line'
	write(6,*) line
	return
   97	continue
	write(6,*) 'no end of section found: ',sname
	write(6,*) line
	return
	end function nls_read_table

c******************************************************************

	function nls_read_ictable()

c reads table in section
c
c table must have following structure (empty lines are allowed)
c
c	ivalue1 'char1'
c	ivalue2 'char2'
c	etc..
c
c an alternative is also possible, where only numbers are given
c
c	ivalue1 ivalue2 ...
c	etc..
c
c in this case char1 gets a default value
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

