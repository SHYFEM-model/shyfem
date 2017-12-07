c
c $Id: subpar.f,v 1.13 2010-02-16 16:21:37 georg Exp $
c
c parameter management routines
c
c contents :
c
c function manpar(name,value,id,what)   internal routine to access par data
c function getpar(name)                 gets parameter name in value
c function itspar(name)                 tests if name is available
c function iscpar(name,sect)            tests if name is in section
c subroutine get_sect_of(name,sect)	returns section of parameter name
c subroutine sctpar(sect)               sets default name of section
c subroutine putpar(name,value)         puts parameter value in name
c subroutine addpar(name,value)         adds parameter value to name
c function infpar(type)                 info about parameters
c subroutine lstpar(name,value,id)      lists entry id in name and value
c function intpar(name)                 tests if name is integer
c subroutine chapar                     changes parameters interactively
c subroutine pripar                     prints parameter values
c
c revision log :
c
c revised on  01.03.88  by ggu	written
c revised on  30.08.88  by ggu  (default values in rdpar..a/h)
c revised on  08.11.88  by ggu  (addpar,setpar,ittpar,ipadim...)
c revised on  05.12.88  by ggu  (rdpara/h substituted by nlsa/h)
c revised on  30.09.89  by ggu  (putpar,pripar,intpar)
c revised on  26.05.90  by ggu  (newly structured -> manpar)
c revised on  04.02.91  by ggu  (included iar..., far...)
c revised on  15.05.97  by ggu  (nnamdi set to 200)
c revised on  12.06.97  by ggu  (iar..., far... moved to subiar.f)
c revised on  12.06.97  by ggu  (section introduced)
c 18.03.1998    ggu     introduced undocumented feature -> c
c 18.03.1998    ggu     save secpar (bug uncovered by g77)
c 07.11.2005    ggu     helper routine get_sect_of()
c 11.09.2006    ggu     routine chapar removed
c 16.04.2008    ggu     bugfix in pripar (character*79 -> *80)
c 28.04.2008    ggu     all routines changed to double precision
c 28.04.2008    ggu     three new routines: dgetpar, dputpar, daddpar
c 28.07.2010    ggu     new routines (par and fnm together) -> subst. old ones
c 25.06.2012    ggu     debugged
c 13.02.2015    ggu     limit of string raised from 6 to 10
c 13.04.2017    ggu     new array feature
c
c**************************************************************
c**************************************************************
c**************************************************************

!==================================================================
	module para
!==================================================================

	implicit none

	integer, parameter :: charlen_para = 10	!max length of name/section

	integer, parameter, private  :: para_invalid = -9876543

	logical, parameter, private  :: bstrict = .false.

	integer, parameter, private  :: type_deprecated = -1
	integer, parameter, private  :: type_value = 1
	integer, parameter, private  :: type_array_value = 2
	integer, parameter, private  :: type_string = 3

	type, private :: entry

	  character*(charlen_para) :: name
	  character*(charlen_para) :: section
	  integer :: isize
	  integer :: itype
	  double precision :: value
	  double precision, allocatable :: array(:)
	  character*80 :: string
	  !character(len=:), allocatable :: string

	end type entry

	integer, save, private :: idlast = 0
	integer, save, private :: ndim = 0
	type(entry), save, allocatable :: pentry(:)

	character*(charlen_para), save, private :: def_section = ' '

        INTERFACE para_get_array_value
        MODULE PROCEDURE	 para_get_array_value_d
     +				,para_get_array_value_r
     +				,para_get_array_value_i
        END INTERFACE

        INTERFACE para_add_array_value
        MODULE PROCEDURE	 para_add_array_value_d
     +				,para_add_array_value_r
     +				,para_add_array_value_i
        END INTERFACE

!==================================================================
	contains
!==================================================================

	function para_get_id(name,section)

	integer para_get_id
	character*(*) name,section

	integer id,ifound
	integer ln,ls,lc
	character(len=len(name)) :: namarg
	character(len=len(section)) :: s

	ln = len(trim(name))
	ls = len(trim(section))
	lc = charlen_para

	if( ln > lc .or. ls > lc ) then
	  write(6,*) ln,trim(name)
	  write(6,*) ls,trim(section)
	  write(6,*) 'maximum length allowed: ',lc
	  !ls = 0
	  !ifound = ln/ls
	  !write(6,*) ln,ls,ifound
	  stop 'error stop para_get_id: name or section too long'
	end if

        namarg=trim(name)
        call to_lower(namarg)
        s=trim(section)
        call to_lower(s)

	ifound = 0
	para_get_id = 0

	do id=1,idlast
	  if( pentry(id)%name == namarg ) then
	    ifound = ifound + 1
	    if( s == ' ' .or. pentry(id)%section == s ) then
	      para_get_id = id
	      return
	    end if
	  end if
	end do

	if( ifound > 0 ) then
	  write(6,*) name,ifound,' names found in other section'
	end if

	end function para_get_id

!******************************************************************

	function para_get_id_with_error(name,section,routine)

	integer para_get_id_with_error

	character*(*) name,section,routine

	integer id,itype
	character*80 string

	id = para_get_id(name,section)

	if( id .eq. 0 ) then
	  write(6,*) 'Parameter not found: ',trim(name)
	  write(6,*) 'Calling routine: ',routine
	  stop 'error stop para_get_id_with_error: no such name'
	else
	  itype = pentry(id)%itype
	  if( itype == type_deprecated ) then
	    string = pentry(id)%string
	    if( string == ' ' ) then
	      write(6,*) 'Parameter is not supported anymore: ',trim(name)
	    else
	      write(6,*) 'Parameter is deprecated: ',trim(name)
	      write(6,*) 'A new parameter might have substituted it.'
	      write(6,*) 'Please have a look at this parameter: '
     +				,trim(string)
	    end if
	  stop 'error stop para_get_id_with_error: parameter deprecated'
	  end if
	end if

	para_get_id_with_error = id

	end function para_get_id_with_error

!******************************************************************
!******************************************************************
!******************************************************************

        subroutine para_init_alloc

        type(entry), allocatable :: paux(:)

        if( ndim == 0 ) then
          ndim = 10
          allocate(pentry(ndim))
          return
        else
          ndim = ndim*2
          allocate(paux(ndim))
          paux(1:ndim/2) = pentry(1:ndim/2)
          call move_alloc(paux,pentry)
        end if

        end subroutine para_init_alloc

!******************************************************************

	subroutine para_init_new_id(id)

	integer id

	idlast = idlast + 1
	if( idlast > ndim ) then
          call para_init_alloc
	end if
	id = idlast

	call para_init_id(id)
	
	end subroutine para_init_new_id

!******************************************************************

	subroutine para_init_id(id)

	integer id

	if( id > ndim ) then
	  stop 'error stop para_init_id: ndim'
	end if

	pentry(id)%name = ' '
	pentry(id)%section = def_section
	pentry(id)%isize = 0
	pentry(id)%itype = 0
	pentry(id)%value = 0.
	!pentry(id)%array = NULL
	pentry(id)%string = ' '
	
	end subroutine para_init_id

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine para_set_default_section(section)

	character*(*) section
	
	def_section = section

	end subroutine para_set_default_section

!******************************************************************

	function para_has_name(name)

	logical para_has_name
	character*(*) name

	para_has_name = para_get_id(name,' ') > 0

	end function para_has_name

!******************************************************************

	function para_is_name_in_section(name,section)

	logical para_is_name_in_section
	character*(*) name
	character*(*) section

	para_is_name_in_section = para_get_id(name,section) > 0

	end function para_is_name_in_section

!******************************************************************

	subroutine para_find_section_to_name(name,section)

	character*(*) name
	character*(*) section

	integer id

	id = para_get_id(name,' ')
	!call para_info_id(id)
	if( id == 0 ) then
	  section = ' '
	else
	  section = pentry(id)%section
	end if

	end subroutine para_find_section_to_name

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine para_delete_section(section)

	character*(*) section

	integer id

	do id=idlast,1,-1	!we run backwards - section is probably last
	  if( pentry(id)%section == section ) then
	    call para_delete_id(id)
	  end if 
	end do

	end subroutine para_delete_section

!******************************************************************

	subroutine para_delete_name(name)

	character*(*) name

	integer id

	id = para_get_id(name,' ')
	call para_delete_id(id)

	end subroutine para_delete_name

!******************************************************************

	subroutine para_delete_id(id)

	integer id

	if( id < 1 .or. id > idlast ) return

	if( allocated(pentry(id)%array) ) deallocate(pentry(id)%array)

	if( id /= idlast ) then
	  pentry(id) = pentry(idlast)
	end if

	call para_init_id(idlast)
	idlast = idlast - 1

	end subroutine para_delete_id

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine para_get_fill(idfill)

	integer idfill

	idfill = idlast

	end subroutine para_get_fill

!******************************************************************

	subroutine para_get_info(id,name,section,itype,value,string)

	integer id
	character*(*) name,section
	integer itype
	double precision value
	character*(*) string

	name = pentry(id)%name
	section = pentry(id)%section
	itype = pentry(id)%itype
	value = pentry(id)%value
	string = pentry(id)%string

	end subroutine para_get_info

!******************************************************************

	subroutine para_info(id_opt)

	integer, optional :: id_opt

	integer id

	if( present(id_opt) ) then
	  id = id_opt
	  call para_info_id(id)
	else
	  do id=1,idlast
	    call para_info_id(id)
	  end do
	end if

	end subroutine para_info

!******************************************************************

	subroutine para_info_name(name)

	character*(*) name

	integer id

	id = para_get_id(name,' ')
	if( id <= 0 ) then
	  write(6,*) 'para_info_name: no such parameter: ',trim(name)
	end if
	call para_info_id(id)

	end subroutine para_info_name

!******************************************************************

	subroutine para_info_id(id)

	integer id

	integer itype,isize
	double precision value
	character(len=len(pentry(id)%name))    :: name
	character(len=len(pentry(id)%section)) :: section
	character(len=len(pentry(id)%string))  :: string

	name = pentry(id)%name
	section = pentry(id)%section
	isize = pentry(id)%isize
	itype = pentry(id)%itype
	value = pentry(id)%value
	string = pentry(id)%string

	if( isize == 0 ) then
	  write(6,1000) id,itype,isize,trim(section),trim(name)
     +			,value,trim(string)
 1000     format(2i4,i5,2(1x,a),g14.6,1x,a)
	else
	  write(6,1100) id,itype,trim(section),trim(name),isize
 1100     format(2i4,2(1x,a),1x,i5)
	  write(6,1200) pentry(id)%array(1:isize)
 1200     format((5g14.6))
	end if

	end subroutine para_info_id

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine para_check_type(id,it_expect,routine)

	integer id
	integer it_expect
	character*(*) routine

	if( pentry(id)%itype /= it_expect ) then
	  write(6,*) 'Wrong type of parameter: ',pentry(id)%name
	  write(6,*) 'itype = ',pentry(id)%itype
	  write(6,*) 'looking for itype = ',it_expect
	  write(6,*) 'Calling routine: ',routine
	  stop 'error stop para_check_type: itype'
	end if

	end subroutine para_check_type

!******************************************************************

	subroutine para_deprecate(name,subst)

! not yet ready

	character*(*) name,subst

	integer id

	call para_add_string(name,subst)
	id = para_get_id_with_error(name,' ','para_deprecate')
	pentry(id)%itype = type_deprecated

	end subroutine para_deprecate

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine para_get_value(name,value)

	character*(*) name
	double precision value

	integer id

	id = para_get_id_with_error(name,' ','para_get_value')

	if( bstrict ) call para_check_type(id,1,'para_get_value')

	value = pentry(id)%value

	end subroutine para_get_value

!******************************************************************

	subroutine para_add_value(name,value)

	character*(*) name
	double precision value

	integer id

	id = para_get_id(name,' ')

	if( id > 0 ) then
	  write(6,*) 'name is already existing: ',name
	  call para_info(id)
	  stop 'error stop para_add_value: name existing'
	end if

	if( id .eq. 0 ) call para_init_new_id(id)

	pentry(id)%name  = name
	pentry(id)%value = value
	pentry(id)%itype = 1

	end subroutine para_add_value

!******************************************************************

	subroutine para_put_value(name,value)

	character*(*) name
	double precision value

	integer id

	id = para_get_id_with_error(name,' ','para_put_value')

	if( bstrict ) call para_check_type(id,1,'para_put_value')

	if( pentry(id)%itype == type_array_value ) then
	  call para_put_array_value(name,value)
	else
	  pentry(id)%value = value
	end if

	end subroutine para_put_value

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine para_add_scalar(name,value)

	character*(*) name
	double precision value

	integer id

	id = para_get_id(name,' ')

	if( id > 0 ) then
	  write(6,*) 'name is already existing: ',name
	  call para_info(id)
	  stop 'error stop para_add_scalar: name existing'
	end if

	if( id .eq. 0 ) call para_init_new_id(id)

	pentry(id)%name  = name
	pentry(id)%value = value
	pentry(id)%itype = -1

	end subroutine para_add_scalar

!******************************************************************
!******************************************************************
!******************************************************************

	function para_is_array_value(name)

	logical para_is_array_value
	character*(*) name

	integer id

	id = para_get_id_with_error(name,' ','para_is_array_value')

	para_is_array_value = ( pentry(id)%itype == type_array_value )

	end function para_is_array_value

!******************************************************************

	subroutine para_add_array_value_i(name,value)

	character*(*) name
	integer value

	double precision dvalue

	dvalue = value
	call para_add_array_value_d(name,dvalue)

	end subroutine para_add_array_value_i

!******************************************************************

	subroutine para_add_array_value_r(name,value)

	character*(*) name
	real value

	double precision dvalue

	dvalue = value
	call para_add_array_value_d(name,dvalue)

	end subroutine para_add_array_value_r

!******************************************************************

	subroutine para_add_array_value_d(name,value)

	character*(*) name
	double precision value

	integer id,isize

	id = para_get_id(name,' ')

	if( id > 0 ) then
	  write(6,*) 'name is already existing: ',name
	  call para_info(id)
	  stop 'error stop para_add_array_value: name existing'
	end if

	if( id .eq. 0 ) call para_init_new_id(id)

	pentry(id)%name  = name
	pentry(id)%value = value	!this is used as default
	pentry(id)%itype = type_array_value

	end subroutine para_add_array_value_d

!******************************************************************

	subroutine para_put_array_value(name,value)

	character*(*) name
	double precision value

	integer id,isize

	id = para_get_id_with_error(name,' ','para_put_array_value')

	if( bstrict ) then
	  call para_check_type(id,type_array_value,'para_put_array_value')
	end if

	pentry(id)%isize = pentry(id)%isize + 1
	isize = pentry(id)%isize
	call para_alloc_array_value(id,isize)
	pentry(id)%array(isize) = value

	end subroutine para_put_array_value

!******************************************************************

	logical function para_has_array(name)

	character*(*) name

	integer id

	id = para_get_id_with_error(name,' ','para_has_array')

	para_has_array = ( pentry(id)%isize > 0 )

	end function para_has_array

!******************************************************************

	subroutine para_get_array_size(name,n)

	character*(*) name
	integer n

	integer id

	id = para_get_id_with_error(name,' ','para_get_array_size')

	n = pentry(id)%isize
	n = max(1,n)

	end subroutine para_get_array_size

!******************************************************************

	subroutine para_get_array_value_i(name,ndim,n,values)

	character*(*) name
	integer ndim,n
	integer values(ndim)

	double precision dvalues(ndim)

	call para_get_array_value_d(name,ndim,n,dvalues)
	values = nint(dvalues)

	end subroutine para_get_array_value_i

!******************************************************************

	subroutine para_get_array_value_r(name,ndim,n,values)

	character*(*) name
	integer ndim,n
	real values(ndim)

	double precision dvalues(ndim)

	call para_get_array_value_d(name,ndim,n,dvalues)
	values = dvalues

	end subroutine para_get_array_value_r

!******************************************************************

	subroutine para_get_array_value_d(name,ndim,n,values)

	character*(*) name
	integer ndim,n
	double precision values(ndim)

	integer id

	id = para_get_id_with_error(name,' ','para_get_array_value')

	if( bstrict ) then
	  call para_check_type(id,type_array_value,'para_get_array_value')
	end if

	n = pentry(id)%isize
	if( ndim == 0 ) then
	  n = max(1,n)
	  return
	end if

	if( ndim < n ) then
	  write(6,*) 'size of array is too big'
	  write(6,*) ndim,n
	  stop 'error stop para_get_array_value: size too big'
	end if

	if( n == 0 ) then
	  n = 1
	  values(1) = pentry(id)%value
	else
	  values = pentry(id)%array(1:n)
	end if

	end subroutine para_get_array_value_d

!******************************************************************

	subroutine para_alloc_array_value(id,isize)

	integer id
	integer isize

	integer ndim
	double precision, allocatable :: paux(:)

	if( .not. allocated(pentry(id)%array) ) then
          ndim = 10
          allocate(pentry(id)%array(ndim))
	end if

	ndim = size(pentry(id)%array)

        if( isize > ndim ) then
          allocate(paux(2*ndim))
          paux(1:ndim) = pentry(id)%array(1:ndim)
          call move_alloc(paux,pentry(id)%array)
        end if

	end subroutine para_alloc_array_value

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine para_get_string(name,string)

	character*(*) name
	character*(*) string

	integer id

	id = para_get_id_with_error(name,' ','para_get_string')

	if( bstrict ) call para_check_type(id,3,'para_get_string')

	string = pentry(id)%string

	end subroutine para_get_string

!******************************************************************

	subroutine para_add_string(name,string)

	character*(*) name
	character*(*) string

	integer id

	id = para_get_id(name,' ')

	if( id > 0 ) then
	  write(6,*) 'name is already existing: ',name
	  call para_info(id)
	  stop 'error stop para_add_string: name existing'
	end if

	if( id .eq. 0 ) call para_init_new_id(id)

	pentry(id)%name   = name
	pentry(id)%string = trim(string)
	pentry(id)%itype = 3

	end subroutine para_add_string

!******************************************************************

	subroutine para_put_string(name,string)

	character*(*) name
	character*(*) string

	integer id

	id = para_get_id_with_error(name,' ','para_put_string')

	if( bstrict ) call para_check_type(id,3,'para_put_string')

	pentry(id)%string = trim(string)

	end subroutine para_put_string

!==================================================================
	end module para
!==================================================================


c**************************************************************
c**************************************************************
c**************************************************************


c***********************************************************
c***********************************************************
c***********************************************************

	function getpar(name)
	use para
	implicit none
	real getpar
	character*(*) name
	double precision value
	call para_get_value(name,value)
	getpar = value
	end

	function dgetpar(name)
	use para
	implicit none
	double precision dgetpar
	character*(*) name
	double precision value
	call para_get_value(name,value)
	dgetpar = value
	end

	subroutine getfnm(name,string)
	use para
	implicit none
	character*(*) name,string
	call para_get_string(name,string)
	end

c**************************************************************

	subroutine putpar(name,value)
	use para
	implicit none
	character*(*) name
	real value
	double precision dvalue
	dvalue = value
	call para_put_value(name,dvalue)
	end

	subroutine dputpar(name,value)
	use para
	implicit none
	character*(*) name
	double precision value
	call para_put_value(name,value)
	end

	subroutine putfnm(name,string)
	use para
	implicit none
	character*(*) name,string
	call para_put_string(name,string)
	end

c**************************************************************

	subroutine addpar(name,value)
	use para
	implicit none
	character*(*) name
	real value
	double precision dvalue
	dvalue = value
	call para_add_value(name,dvalue)
	end

	subroutine daddpar(name,value)
	use para
	implicit none
	character*(*) name
	double precision value
	call para_add_value(name,value)
	end

	subroutine addfnm(name,string)
	use para
	implicit none
	character*(*) name,string
	call para_add_string(name,string)
	end

c**************************************************************
c**************************************************************
c**************************************************************


c**************************************************************
c**************************************************************
c**************************************************************

	subroutine delete_section(section)
	use para
	implicit none
	character*(*) section
	call para_delete_section(section)
	end

c**************************************************************
c**************************************************************
c**************************************************************

	function haspar(name)
	use para
	implicit none
	logical haspar
	character*(*) name
	haspar=para_get_id(name,' ')>0
	end

	function itspar(name)
	use para
	implicit none
	integer itspar
	character*(*) name
	itspar=para_get_id(name,' ')
	end

	function itsfnm(name)
	use para
	implicit none
	integer itsfnm
	character*(*) name
	itsfnm=para_get_id(name,' ')
	end

	function iscpar(name,section)
	use para
	implicit none
	integer iscpar
	character*(*) name,section
	iscpar = para_get_id(name,section)
	end

	function iscfnm(name,section)
	use para
	implicit none
	integer iscfnm
	character*(*) name,section
	iscfnm = para_get_id(name,section)
	end

c***********************************************************

	subroutine sctpar(section)
	use para
	implicit none
	character*(*) section
	call para_set_default_section(section)
	end

	subroutine sctfnm(section)
	use para
	implicit none
	character*(*) section
	call para_set_default_section(section)
	end

c**************************************************************

	subroutine get_sect_of(name,section)
	use para
	implicit none
	character*(*) name,section
	call para_find_section_to_name(name,section)
	end

c**************************************************************
c**************************************************************
c**************************************************************
c**************************************************************
c**************************************************************
c**************************************************************
c**************************************************************
c**************************************************************

	subroutine parinfo(iunit)

c prints info on parameter values

	use para

	implicit none

	integer iunit

	character*(charlen_para) name,section
	character*80 string
	double precision dvalue
	real value
	integer itype,nlen
	integer id,idfill

	integer check_entry_par,ichanm

	call para_get_fill(idfill)

	write(iunit,*) 'parinfo: ',idfill

        do id=1,idfill
	  call para_get_info(id,name,section,itype,dvalue,string)
	  value = dvalue
	  write(iunit,*) id,itype,value,name,'  ',section
	  if( itype .eq. 3 ) then
            nlen=max(1,ichanm(string))
            write(iunit,*) '    ',nlen,string(1:nlen)
	  end if
	end do

	end

c**************************************************************

        subroutine pripar(iunit)

        implicit none

        integer iunit

        !call pripar1(iunit)
        call pripar3(iunit)
        !call pripar4(iunit)

        end

c*****************************************************

        subroutine pripar1(iunit)

c prints parameter values (1 column)

	use para

        implicit none

        integer iunit

	character*(charlen_para) name,section
	character*80 string
        integer idfill,id,itype
        double precision dvalue
        real value

	call para_get_fill(idfill)

        do id=1,idfill
	  call para_get_info(id,name,section,itype,dvalue,string)
	  value = dvalue
	  if(itype.eq.1) then
            write(iunit,*) name,'  ',value
	  end if
        end do

	end

c*****************************************************

        subroutine pripar3(iunit)

c prints parameter values (3 columns)

	use para

        implicit none

        integer iunit

	character*(charlen_para) name,section
        character*(charlen_para) nameii(3)
	character*80 string
        integer idfill,id,itype,ii
        double precision dvalue
        double precision valueii(3)
        real value

	call para_get_fill(idfill)

	ii = 0
        do id=1,idfill
	  call para_get_info(id,name,section,itype,dvalue,string)
	  value = dvalue
	  if(itype.eq.1) then
	    ii = ii + 1
	    nameii(ii) = name
	    valueii(ii) = value
	  end if
	  if( ii == 3 ) then
            write(iunit,1000) (nameii(ii),' = ',valueii(ii),ii=1,3)
	    ii = 0
	  end if
        end do

	return
 1000	format(3(a8,a3,g14.6,2x))
	end

c**************************************************************

	subroutine pripar4(iunit)

c prints parameter values (4 columns)

	use para

	implicit none

	integer iunit

	logical bflag
	character*80 line	!BUGFIX (was 79)

	integer npara,imod,i,itype
	integer ianf,iend
	integer itspar,infpar,intpar
	integer id,idfill
	character*(charlen_para) name,section
	character*80 string
	double precision dvalue
	real value
	real flag,absval

	line=' '
	bflag=.false.
	flag = -999.

	call para_get_fill(idfill)
	npara=idfill

	imod=0
	do id=1,npara
	  call para_get_info(id,name,section,itype,dvalue,string)
	  value = dvalue
	  if(itype.ne.1) goto 1
	  if(bflag.and.value.eq.flag) goto 1
	  imod=imod+1
	  ianf=20*(imod-1)+1
	  iend=20*imod
	  absval=abs(value)
	  if(intpar(name).eq.1) then
            write(line(ianf:iend),2345)
     +                  name,' =',nint(value),'  '

	  else
	    if(absval.lt.1000.and.absval.ge.0.01
     +                        .or.absval.eq.0.) then
	      write(line(ianf:iend),2347) name,' =',value,'  '
	    else
	      write(line(ianf:iend),2346) name,' =',value,'  '
	    end if
	  end if
	  if(imod.eq.4) then
		write(iunit,*) line(1:79)
		line=' '
		imod=0
	  end if
    1     continue
	end do

	if(imod.ne.4) write(iunit,*) line

	return
 2345   format(a6,a2,i10,a2)
 2346   format(a6,a2,e10.2,a2)
 2347   format(a6,a2,f10.3,a2)
	end

c**********************************************************

        subroutine prifnm(iunit)

c prints parameter values

	use para

        implicit none

        integer iunit

        character*80 string
	character*(charlen_para) name,section
        integer nlen,itype
	integer id,idfill
	double precision dvalue

        integer ichanm

	call para_get_fill(idfill)

        do id=1,idfill
	  call para_get_info(id,name,section,itype,dvalue,string)
          if(itype.eq.3) then
            nlen=max(1,ichanm(string))
            write(iunit,2345) id,nlen,name,section,string(1:nlen)
          end if
        end do

        return
 2345   format(1x,2i4,2(1x,a6,1x),3x,a)
        end
 
c**********************************************************

	subroutine chkparam(iunit)

	use para

	implicit none

	integer iunit

        character*80 string
	character*(charlen_para) name,section
        integer nlen,itype
	integer id,idfill
	double precision dvalue
	real value

        integer ichanm

	call para_get_fill(idfill)

        do id=1,idfill
	  call para_get_info(id,name,section,itype,dvalue,string)
	  value = dvalue
          if(itype.eq.1) then
            write(iunit,2345) id,name,section,itype,value
          else if(itype.eq.3) then
            write(iunit,2346) id,name,section,itype,trim(string)
	  else
	    write(6,*) 'error itype...'
            write(iunit,2345) id,name,section,itype,value
	    stop 'error stop chkparam'
          end if
        end do

        return
 2345   format(1x,i4,2(1x,a6,1x),i4,e12.4)
 2346   format(1x,i4,2(1x,a6,1x),i4,a)
	end

c**********************************************************

	subroutine check_parameter_values(text)

	use para

	implicit none

	character*(*) text

	integer iunit,idfill

	iunit = 117

	return

	call para_get_fill(idfill)

	write(iunit,*) '--------------------------------'
	write(iunit,*) 'start of print in check_parameter_values'
	write(iunit,*) 'info on parameters: ',text
	call chkparam(iunit)
	write(iunit,*) '--------------------------------'
	write(iunit,*) '...printing with pripar...'
        call pripar(iunit)
	write(iunit,*) '...printing with prifnm...'
        call prifnm(iunit)
	write(iunit,*) 'end info on parameters: ',text
	write(iunit,*) 'end of print in check_parameter_values'

	end

c**********************************************************

        function intpar(name)

c tests if name is integer
c
c name          name to test
c intpar        1 if name is integer, 0 if not

        implicit none

        integer intpar
        character*(*) name

        integer i
        character*6 let
        character*1 namein
        data let /'ijklmn'/

        namein=name(1:1)
        call to_lower(namein)

        intpar=0
        do i=1,6
          if(namein.eq.let(i:i)) intpar=1
        end do

        end

c*****************************************************
c	program para_main
c	call test_par
c	end
c*****************************************************

