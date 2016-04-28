!
! handle command line options
!
! revision log :
!
! 24.11.2014    ggu     written from scratch
! 27.11.2014    ggu     is now fully functional
! 14.01.2015    ggu     added new routines clo_add_info/extra
! 26.01.2015    ggu     make argument to clo_parse_options() optional
! 26.01.2015    ggu     new routine clo_check_files()
! 08.02.2015    ggu     bug fix in clo_get_option, new routine clo_info
! 28.05.2015    ggu     new routine clo_add_sep()
!
! notes :
!
! handles automatically options -h|-help and -v|-version
! creates automatically clo_usage and clo_fullusage
! 
! in clo_add/set/get_option numbers can be integer, real or double precision
! in clo_add_option give extra information for number or string (see example)
!
! usage :
!
!        call clo_init('feminf','fem-file','1.2')
!
!	 call clo_add_info('elaborates and rewrites a fem file')
!        call clo_add_option('write',.false.,'write min/max of values')
!        call clo_add_option('out',.false.,'create output file')
!	 call clo_add_sep('other options')
!        call clo_add_option('tmin time',-1
!     +                          ,'only process starting from time')
!        call clo_add_option('tmax time',-1
!     +                          ,'only process up to time')
!	 call clo_add_extra('time is YYYY-MM-DD[::hh:mm:ss]')
!
!        call clo_parse_options(1)       !expecting 1 file
!
!        call clo_get_option('write',bwrite)
!        call clo_get_option('out',bout)
!        call clo_get_option('tmin',tmin)
!        call clo_get_option('tmax',tmax)
!
!        nfile = clo_number_of_files()
!        if( nfile > 0 ) call clo_get_file(1,infile)
!
!**************************************************************

!==================================================================
	module clo
!==================================================================

	implicit none

	type, private :: entry

	  character*80 :: name		! name of option
	  integer :: itype		! type (1: number  2: flag  3: string)
	  double precision :: value	! value if number
	  logical :: flag		! flag if flag
	  character*80 :: string	! string if string
	  character*80 :: text		! description for clo_fullusage
	  character*80 :: textra	! if number or string extra info

	end type entry

	integer, save, private :: idlast = 0
	integer, save, private :: ndim = 0
	type(entry), save, private, allocatable :: pentry(:)

	integer, save, private :: last_option = 0
	integer, save, private :: last_file = 0

	integer, save, private :: ielast = 0
	character*80, save, private :: info = ' '
	character*80, save, private, dimension(10) :: extra = ' '

	character*80, save, private :: routine_name = ' '
	character*80, save, private :: files_name = ' '
	character*80, save, private :: version_name = ' '

	INTERFACE clo_add_option
        MODULE PROCEDURE clo_add_option_d, clo_add_option_r, 
     +			 clo_add_option_i,
     +			 clo_add_option_f, clo_add_option_s
	END INTERFACE

	INTERFACE clo_set_option
        MODULE PROCEDURE clo_set_option_d, clo_set_option_r, 
     +			 clo_set_option_i,
     +			 clo_set_option_f, clo_set_option_s
	END INTERFACE

	INTERFACE clo_get_option
        MODULE PROCEDURE clo_get_option_d, clo_get_option_r, 
     +			 clo_get_option_i,
     +			 clo_get_option_f, clo_get_option_s
	END INTERFACE

!==================================================================
	contains
!==================================================================

	function clo_get_id(name)

	integer clo_get_id
	character*(*) name

	integer id

	clo_get_id = 0

	do id=1,idlast
	  if( pentry(id)%name == name ) then
	    clo_get_id = id
	    return
	  end if
	end do

	end function clo_get_id

!******************************************************************
!******************************************************************
!******************************************************************

        subroutine clo_init_alloc

        type(entry), allocatable :: paux(:)

        if( ndim == 0 ) then
          ndim = 10
          allocate(pentry(0:ndim))
          return
        else
          ndim = ndim*2
          allocate(paux(0:ndim))
          paux(0:ndim/2) = pentry(0:ndim/2)
          call move_alloc(paux,pentry)
        end if

        end subroutine clo_init_alloc

!******************************************************************

	subroutine clo_init_new_id(id)

	integer id

	idlast = idlast + 1
	if( idlast > ndim ) then
          call clo_init_alloc
	end if
	id = idlast

	call clo_init_id(id)
	
	end subroutine clo_init_new_id

!******************************************************************

	subroutine clo_init_id(id)

	integer id

	if( id > ndim ) then
	  stop 'error stop clo_init_id: ndim'
	end if

	pentry(id)%name = ' '
	pentry(id)%itype = 0
	pentry(id)%value = 0.
	pentry(id)%flag = .false.
	pentry(id)%string = ' '
	pentry(id)%text = ' '
	pentry(id)%textra = ' '
	
	end subroutine clo_init_id

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine clo_info(name)

	character*(*), optional :: name
	integer id

	if( present(name) ) then
	  id = clo_get_id(name)
	  if( id > 0 ) then
	    call clo_info_id(id)
	  else
	    write(6,*) 'clo_info: no option with this name: ',trim(name)
	  end if
	else
	  do id=1,idlast
	    call clo_info_id(id)
	  end do
	end if

	end subroutine clo_info

!******************************************************************

	subroutine clo_info_id(id)

	integer id

	write(6,*) 'clo info on id = ',id
	write(6,*) 'name   = ',pentry(id)%name
	write(6,*) 'itype  = ',pentry(id)%itype
	write(6,*) 'value  = ',pentry(id)%value
	write(6,*) 'flag   = ',pentry(id)%flag
	write(6,*) 'string = ',trim(pentry(id)%string)
	write(6,*) 'text   = ',trim(pentry(id)%text)

	end subroutine clo_info_id

!******************************************************************
!******************************************************************
!******************************************************************

	function clo_has_option(name)

	logical clo_has_option
	character*(*) name

	clo_has_option = clo_get_id(name) > 0

	end function clo_has_option

!******************************************************************

	subroutine clo_error(name,text)

	character*(*) name
	character*(*) text

	write(6,*) 'option ',trim(name),' error: ',trim(text)

	stop 'error stop clo_error'
	end subroutine clo_error

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine clo_set_option_d(name,value)

	character*(*) name
	double precision value

	call clo_set_option_n(name,value)

	end subroutine clo_set_option_d

!******************************************************************

	subroutine clo_set_option_r(name,value)

	character*(*) name
	real value
	double precision dvalue

	dvalue = value
	call clo_set_option_n(name,dvalue)

	end subroutine clo_set_option_r

!******************************************************************

	subroutine clo_set_option_i(name,value)

	character*(*) name
	integer value
	double precision dvalue

	dvalue = value
	call clo_set_option_n(name,dvalue)

	end subroutine clo_set_option_i

!******************************************************************

	subroutine clo_set_option_n(name,value)

	character*(*) name
	double precision value

	integer id

	id = clo_get_id(name)
	if( id == 0 ) call clo_error(name,'option not existing')
	if( pentry(id)%itype /= 1 ) call clo_error(name,'wrong type')

	pentry(id)%value = value

	end subroutine clo_set_option_n

!******************************************************************

	subroutine clo_set_option_f(name,flag)

	character*(*) name
	logical flag

	integer id

	id = clo_get_id(name)
	if( id == 0 ) call clo_error(name,'option not existing')
	if( pentry(id)%itype /= 2 ) call clo_error(name,'wrong type')

	pentry(id)%flag = flag

	end subroutine clo_set_option_f

!******************************************************************

	subroutine clo_set_option_s(name,string)

	character*(*) name
	character*(*) string

	integer id

	id = clo_get_id(name)
	if( id == 0 ) call clo_error(name,'option not existing')
	if( pentry(id)%itype /= 3 ) call clo_error(name,'wrong type')

	pentry(id)%string = string

	end subroutine clo_set_option_s

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine clo_add_option_d(name,value,text)

	character*(*) name
	double precision value
	character*(*) text
	optional text

	call clo_add_option_n(name,value,text)

	end subroutine clo_add_option_d

!******************************************************************

	subroutine clo_add_option_r(name,value,text)

	character*(*) name
	real value
	character*(*) text
	optional text

	double precision dvalue

	dvalue = value
	call clo_add_option_n(name,dvalue,text)

	end subroutine clo_add_option_r

!******************************************************************

	subroutine clo_add_option_i(name,value,text)

	character*(*) name
	integer value
	character*(*) text
	optional text

	double precision dvalue

	dvalue = value
	call clo_add_option_n(name,dvalue,text)

	end subroutine clo_add_option_i

!******************************************************************

	subroutine clo_add_option_n(name,value,text)

	character*(*) name
	double precision value
	character*(*) text
	optional text

	integer id
	character*80 name1,name2

	call clo_split_name(name,name1,name2)
	if( name2 == ' ' ) name2 = 'val'

	id = clo_get_id(name1)
	if( id /= 0 ) call clo_error(name1,'option already existing')
	call clo_init_new_id(id)

	pentry(id)%name = name1
	pentry(id)%itype = 1
	pentry(id)%value = value
	pentry(id)%textra = name2
	if( present(text) ) pentry(id)%text = text

	end subroutine clo_add_option_n

!******************************************************************

	subroutine clo_add_option_f(name,flag,text)

	character*(*) name
	logical flag
	character*(*) text
	optional text

	integer id
	character*80 name1,name2

	call clo_split_name(name,name1,name2)
	name2 = ' '

	id = clo_get_id(name1)
	if( id /= 0 ) call clo_error(name1,'option already existing')
	call clo_init_new_id(id)

	pentry(id)%name = name
	pentry(id)%itype = 2
	pentry(id)%flag = flag
	pentry(id)%textra = name2
	if( present(text) ) pentry(id)%text = text

	end subroutine clo_add_option_f

!******************************************************************

	subroutine clo_add_option_s(name,string,text)

	character*(*) name
	character*(*) string
	character*(*) text
	optional text

	integer id
	character*80 name1,name2

	call clo_split_name(name,name1,name2)
	if( name2 == ' ' ) name2 = 'string'

	id = clo_get_id(name1)
	if( id /= 0 ) call clo_error(name1,'option already existing')
	call clo_init_new_id(id)

	pentry(id)%name = name1
	pentry(id)%itype = 3
	pentry(id)%string = string
	pentry(id)%textra = name2
	if( present(text) ) pentry(id)%text = text

	end subroutine clo_add_option_s

!******************************************************************

	subroutine clo_add_sep(text)

	character*(*) text

	integer id

	call clo_init_new_id(id)

	pentry(id)%name = ' '
	pentry(id)%itype = 4
	pentry(id)%text = text

	end subroutine clo_add_sep

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine clo_get_option_d(name,value)

	character*(*) name
	double precision value

	call clo_get_option_n(name,value)

	end subroutine clo_get_option_d

!******************************************************************

	subroutine clo_get_option_r(name,value)

	character*(*) name
	real value

	double precision dvalue

	call clo_get_option_n(name,dvalue)
	value = dvalue

	end subroutine clo_get_option_r

!******************************************************************

	subroutine clo_get_option_i(name,value)

	character*(*) name
	integer value

	double precision dvalue

	call clo_get_option_n(name,dvalue)
	value = dvalue

	end subroutine clo_get_option_i

!******************************************************************

	subroutine clo_get_option_n(name,value)

	character*(*) name
	double precision value

	integer id

	id = clo_get_id(name)
	if( id == 0 ) call clo_error(name,'option not existing')
	if( pentry(id)%itype /= 1 ) call clo_error(name,'wrong type')

	value = pentry(id)%value

	end subroutine clo_get_option_n

!******************************************************************

	subroutine clo_get_option_f(name,flag)

	character*(*) name
	logical flag

	integer id

	id = clo_get_id(name)
	if( id == 0 ) call clo_error(name,'option not existing')
	if( pentry(id)%itype /= 2 ) call clo_error(name,'wrong type')

	flag = pentry(id)%flag

	end subroutine clo_get_option_f

!******************************************************************

	subroutine clo_get_option_s(name,string)

	character*(*) name
	character*(*) string

	integer id

	id = clo_get_id(name)
	if( id == 0 ) call clo_error(name,'option not existing')
	if( pentry(id)%itype /= 3 ) call clo_error(name,'wrong type')

	string = pentry(id)%string

	end subroutine clo_get_option_s

!******************************************************************
!******************************************************************
!******************************************************************

	function clo_is_number(name)

	logical clo_is_number
	character*(*) name

	integer id

	clo_is_number = .false.
	id = clo_get_id(name)
	if( id == 0 ) return
	if( pentry(id)%itype /= 1 ) return
	clo_is_number = .true.

	end function clo_is_number

!******************************************************************

	function clo_is_flag(name)

	logical clo_is_flag
	character*(*) name

	integer id

	clo_is_flag = .false.
	id = clo_get_id(name)
	if( id == 0 ) return
	if( pentry(id)%itype /= 2 ) return
	clo_is_flag = .true.

	end function clo_is_flag

!******************************************************************

	function clo_is_string(name)

	logical clo_is_string
	character*(*) name

	integer id

	clo_is_string = .false.
	id = clo_get_id(name)
	if( id == 0 ) return
	if( pentry(id)%itype /= 3 ) return
	clo_is_string = .true.

	end function clo_is_string

!******************************************************************
!******************************************************************
!******************************************************************

	function clo_has_files()

	logical clo_has_files

	clo_has_files = clo_number_of_files() > 0

	end function clo_has_files

!******************************************************************

	function clo_number_of_files()

	integer clo_number_of_files

	clo_number_of_files = last_file - last_option

	end function clo_number_of_files

!******************************************************************

	subroutine clo_get_file(i,file)

	integer i
	character*(*) file

	file = ' '
	if( i < 1 ) return
	if( last_option + i > last_file ) return

	call get_command_argument(last_option+i,file)

	end subroutine clo_get_file

!**************************************************************

	subroutine clo_check_files(nexpect)

	integer nexpect

	if( nexpect > last_file - last_option ) then
	  call clo_usage
	end if

	end subroutine clo_check_files

!**************************************************************
!**************************************************************
!**************************************************************

	subroutine clo_init(routine,files,version)

	character*(*) routine
	character*(*) files
	character*(*) version

	routine_name = routine
	files_name = files
	version_name = version

	end subroutine clo_init

!**************************************************************

	subroutine clo_parse_options(opt_nexpect)

	integer, optional :: opt_nexpect  !number of expected files (at least)

	integer nexpect
	integer nc,i,n
	character*80 option
	character*80 string
	double precision value

	nexpect = 0
	if( present(opt_nexpect) ) nexpect = opt_nexpect

	nc = command_argument_count()
	last_file = nc
	last_option = 0

	i = 0

	do
	  i = i + 1
	  if( i > nc ) exit
	  option = ' '
	  call get_command_argument(i,option)
	  n = len_trim(option)
	  if( option(1:1) /= '-' ) exit

	  option = option(2:)
	  if( option == 'h' .or. option == 'help' ) then
	    call clo_fullusage
	  else if( option == 'v' .or. option == 'version' ) then
	    call clo_version
	  else if( .not. clo_has_option(option) ) then
	    write(6,*) '*** no such option: ',option
	    call clo_usage
	  else if( clo_is_flag(option) ) then
	    call clo_set_option(option,.true.)
	  else
	    i = i + 1
	    if( i > nc ) then
	      write(6,*) '*** no value for option ',option
	      call clo_usage
	    end if
	    call get_command_argument(i,string)
	    if( clo_is_number(option) ) then
	      call clo_s2d(string,value)
	      call clo_set_option(option,value)
	    else if( clo_is_string(option) ) then
	      call clo_set_option(option,string)
	    else
	      stop 'error stop clo_parse_options: internal error (1)'
	    end if
	    
	  end if
	  last_option = i
	end do

	if( nc-i+1 < nexpect ) call clo_usage

	end subroutine clo_parse_options

!**************************************************************

	subroutine clo_s2d(string,value)

	character*(*) string
	double precision value

	read(string,'(g20.0)',err=99) value

	return
   99	continue
	write(6,*) '*** cannot convert to number: ',string
	stop 'error stop clo_s2d: cannot convert'
	end subroutine clo_s2d

!**************************************************************

	subroutine clo_split_name(name,name1,name2)

	character*(*) name,name1,name2

	integer i

	do i=1,len(name)
	  if( name(i:i) == ' ' ) exit
	end do

	name1 = name(1:i-1)
	name2 = name(i+1:)

	end subroutine clo_split_name

!**************************************************************
!**************************************************************
!**************************************************************

	subroutine clo_add_info(string)

	character*(*) string

	info = string

	end subroutine clo_add_info

!**************************************************************

	subroutine clo_add_extra(string)

	character*(*) string

	ielast = ielast + 1
	if( ielast > 10 ) stop 'error stop clo_add_extra: ielast'

	extra(ielast) = string

	end subroutine clo_add_extra

!**************************************************************
!**************************************************************
!**************************************************************

	subroutine clo_version

	integer nr,nv

	nr = len_trim(routine_name)
	nv = len_trim(version_name)

	write(6,*) routine_name(1:nr),' version ',version_name(1:nv)

	stop
	end subroutine clo_version

!**************************************************************

	subroutine clo_usage

	integer nr,nf

	nr = len_trim(routine_name)
	nf = len_trim(files_name)

	write(6,*) 'Usage: ',routine_name(1:nr)
     +			,' [-h|-help] [-options] ',files_name(1:nf)

	stop
	end subroutine clo_usage

!**************************************************************

	subroutine clo_fullusage

	integer nr,nf,nn,nt,ne,nl,np
	integer itype,id,ie
	integer length,l
	character*80 name
	character*80 text
	character*80 textra

	nr = len_trim(routine_name)
	nf = len_trim(files_name)

	length = 0
	do id=1,idlast
	  l = clo_get_length(pentry(id)%name,pentry(id)%textra)
	  length = max(length,l)
	end do
	length = max(length,clo_get_length('h|-help',' '))
	length = max(length,clo_get_length('v|-version',' '))
	length = length + 5

	write(6,*) 'Usage: ',routine_name(1:nr)
     +			,' [-h|-help] [-options] ',files_name(1:nf)
	if( info /= ' ' ) write(6,*) ' ',info(1:len_trim(info))
	write(6,*) ' options:'
	call clo_write_line(length,'h|-help',' ','this help screen')
	call clo_write_line(length,'v|-version',' ','version of routine')

	do id=1,idlast
	  name = pentry(id)%name
	  textra = pentry(id)%textra
	  text = pentry(id)%text
	  if( name == ' ' ) then
	    write(6,*) ' ',trim(text)
	  else
	    call clo_write_line(length,name,textra,text)
	  end if
	end do

	do ie=1,ielast
	  text = extra(ie)
	  write(6,*) '  ',text(1:len_trim(text))
	end do

	stop
	end subroutine clo_fullusage

!**************************************************************

	subroutine clo_write_line(length,text1,text2,text3)

	integer length
	character*(*) text1,text2,text3

	integer n1,n2,n3,nl,np
	character*80 line,empty

	n1 = len_trim(text1)
	n2 = len_trim(text2)
	n3 = len_trim(text3)

	line = '  -'//text1(1:n1)//' '//text2(1:n2)
	empty = ' '
	nl = len_trim(line)
	np = length - nl

	write(6,*) line(1:nl),empty(1:np),text3(1:n3)

	end subroutine clo_write_line

!**************************************************************

	function clo_get_length(text1,text2)

	integer clo_get_length
	character*(*) text1,text2

	clo_get_length = 1
	clo_get_length = clo_get_length + len_trim(text1)
	clo_get_length = clo_get_length + len_trim(text2)

	end function clo_get_length

!==================================================================
	end module clo
!==================================================================

	subroutine clo_test
	end

!*****************************************************
!	program clo_main
!	call clo_test
!	end
!*****************************************************

