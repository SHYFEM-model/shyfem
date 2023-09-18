!
! $Id: subpar.f,v 1.13 2010-02-16 16:21:37 georg Exp $
!
! parameter management routines
!
! contents :
!
! function manpar(name,value,id,what)   internal routine to access par data
! function getpar(name)                 gets parameter name in value
! function itspar(name)                 tests if name is available
! function iscpar(name,sect)            tests if name is in section
! subroutine get_sect_of(name,sect)	returns section of parameter name
! subroutine sctpar(sect)               sets default name of section
! subroutine putpar(name,value)         puts parameter value in name
! subroutine addpar(name,value)         adds parameter value to name
! function infpar(type)                 info about parameters
! subroutine lstpar(name,value,id)      lists entry id in name and value
! function intpar(name)                 tests if name is integer
! subroutine chapar                     changes parameters interactively
! subroutine pripar                     prints parameter values
! subroutine impini		initializes parameters for semi-implicit time
! function getimp		gets weight for semi-implicit time-step
! subroutine setimp(it,aweigh)	sets parameters for semi-implicit time-step
!
! revision log :
!
! revised on  01.03.88  by ggu	written
! revised on  30.08.88  by ggu  (default values in rdpar..a/h)
! revised on  08.11.88  by ggu  (addpar,setpar,ittpar,ipadim...)
! revised on  05.12.88  by ggu  (rdpara/h substituted by nlsa/h)
! revised on  30.09.89  by ggu  (putpar,pripar,intpar)
! revised on  26.05.90  by ggu  (newly structured -> manpar)
! revised on  04.02.91  by ggu  (included iar..., far...)
! revised on  15.05.97  by ggu  (nnamdi set to 200)
! revised on  12.06.97  by ggu  (iar..., far... moved to subiar.f)
! revised on  12.06.97  by ggu  (section introduced)
! 18.03.1998    ggu     introduced undocumented feature -> c
! 18.03.1998    ggu     save secpar (bug uncovered by g77)
! 07.02.2003    ggu     routine added: getaz; deleted getaza
! 07.11.2005    ggu     helper routine get_sect_of()
! 11.09.2006    ggu     routine chapar removed
! 16.04.2008    ggu     bugfix in pripar (character*79 -> *80)
! 28.04.2008    ggu     all routines changed to double precision
! 28.04.2008    ggu     three new routines: dgetpar, dputpar, daddpar
! 28.07.2010    ggu     new routines (par and fnm together) -> subst. old ones
! 25.06.2012    ggu     debugged
! 13.02.2015    ggu     limit of string raised from 6 to 10
!
!**************************************************************
!**************************************************************
!**************************************************************

!==================================================================
	module para
!==================================================================

! itype =  1	value
! itype = 10	string

	implicit none

	integer, parameter, private :: lc = 10		!max length of name

	type, private :: entry

	  character*10 :: name
	  character*10 :: section
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

	character*10, save, private :: def_section = ' '

!==================================================================
	contains
!==================================================================

	function para_get_id(name,section)

        use utility

	integer para_get_id
	character*(*) name,section

	integer id,ifound
	integer ln,ls
	character(len=len(name)) :: namarg
	character(len=len(section)) :: s

	ln = len(trim(name))
	ls = len(trim(section))

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
        call uplow(namarg,'low')
        s=trim(section)
        call uplow(s,'low')

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

	subroutine para_info_id(id)

	integer id

	integer itype
	double precision value
	character(len=len(pentry(id)%name))    :: name
	character(len=len(pentry(id)%section)) :: section
	character(len=len(pentry(id)%string))  :: string

	name = pentry(id)%name
	section = pentry(id)%section
	itype = pentry(id)%itype
	value = pentry(id)%value
	string = pentry(id)%string

	!write(6,*) id,itype,trim(name),value,trim(string)
	write(6,1000) id,itype,trim(section),trim(name),value,trim(string)
 1000   format(2i4,2(1x,a),g14.6,1x,a)

	end subroutine para_info_id

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine para_get_value(name,value)

	character*(*) name
	double precision value

	integer id
	logical bstrict

	bstrict = .false.

	id = para_get_id(name,' ')

	if( id .eq. 0 ) then
	  write(6,*) 'Parameter not found: ',name
	  stop 'error stop para_get_value: name'
	end if

	if( bstrict .and. pentry(id)%itype /= 1 ) then
	  write(6,*) 'Wrong type of parameter: ',name
	  write(6,*) 'itype = ',pentry(id)%itype,' looking for itype = 1'
	  stop 'error stop para_get_value: itype'
	end if

	value = pentry(id)%value

	end subroutine para_get_value

!******************************************************************

	subroutine para_add_value(name,value)

	character*(*) name
	double precision value

	integer id

	id = para_get_id(name,' ')

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
	logical bstrict

	bstrict = .false.

	id = para_get_id(name,' ')

	if( id .eq. 0 ) then
	  write(6,*) 'Parameter not found: ',name
	  stop 'error stop para_put_value: name'
	end if

	if( bstrict .and. pentry(id)%itype /= 1 ) then
	  write(6,*) 'Wrong type of parameter: ',name
	  write(6,*) 'itype = ',pentry(id)%itype,' looking for itype = 1'
	  stop 'error stop para_put_value: itype'
	end if

	pentry(id)%value = value

	end subroutine para_put_value

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine para_get_string(name,string)

	character*(*) name
	character*(*) string

	integer id
	logical bstrict

	bstrict = .false.

	id = para_get_id(name,' ')

	if( id .eq. 0 ) then
	  write(6,*) 'Parameter not found: ',name
	  stop 'error stop para_get_string: name'
	end if

	if( bstrict .and. pentry(id)%itype /= 3 ) then
	  write(6,*) 'Wrong type of parameter: ',name
	  write(6,*) 'itype = ',pentry(id)%itype,' looking for itype = 3'
	  stop 'error stop para_get_string: itype'
	end if

	string = pentry(id)%string

	end subroutine para_get_string

!******************************************************************

	subroutine para_add_string(name,string)

	character*(*) name
	character*(*) string

	integer id

	id = para_get_id(name,' ')

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
	logical bstrict

	bstrict = .false.

	id = para_get_id(name,' ')

	if( id .eq. 0 ) then
	  write(6,*) 'Parameter not found: ',name
	  stop 'error stop para_put_string: name'
	end if

	if( bstrict .and. pentry(id)%itype /= 3 ) then
	  write(6,*) 'Wrong type of parameter: ',name
	  write(6,*) 'itype = ',pentry(id)%itype,' looking for itype = 3'
	  stop 'error stop para_put_string: itype'
	end if

	pentry(id)%string = trim(string)

	end subroutine para_put_string


!**************************************************************
!**************************************************************
!**************************************************************


!***********************************************************
!***********************************************************
!***********************************************************

	function getpar(name)
	implicit none
	double precision getpar
	character*(*) name
	double precision value
	call para_get_value(name,value)
	getpar = value
	end

	function dgetpar(name)
	implicit none
	double precision dgetpar
	character*(*) name
	double precision value
	call para_get_value(name,value)
	dgetpar = value
	end

	subroutine getfnm(name,string)
	implicit none
	character*(*) name,string
	call para_get_string(name,string)
	end

!**************************************************************

	subroutine putpar(name,value)
	implicit none
	character*(*) name
	double precision value
	double precision dvalue
	dvalue = value
	call para_put_value(name,dvalue)
	end

	subroutine dputpar(name,value)
	implicit none
	character*(*) name
	double precision value
	call para_put_value(name,value)
	end

	subroutine putfnm(name,string)
	implicit none
	character*(*) name,string
	call para_put_string(name,string)
	end

!**************************************************************

	subroutine addpar(name,value)
	implicit none
	character*(*) name
	double precision value
	double precision dvalue
	dvalue = value
	call para_add_value(name,dvalue)
	end

	subroutine daddpar(name,value)
	implicit none
	character*(*) name
	double precision value
	call para_add_value(name,value)
	end

	subroutine addfnm(name,string)
	implicit none
	character*(*) name,string
	call para_add_string(name,string)
	end

!**************************************************************
!**************************************************************
!**************************************************************


!**************************************************************
!**************************************************************
!**************************************************************

	subroutine delete_section(section)
	implicit none
	character*(*) section
	call para_delete_section(section)
	end

!**************************************************************
!**************************************************************
!**************************************************************

	function haspar(name)
	implicit none
	logical haspar
	character*(*) name
	haspar=para_get_id(name,' ')>0
	end

	function itspar(name)
	implicit none
	integer itspar
	character*(*) name
	itspar=para_get_id(name,' ')
	end

	function itsfnm(name)
	implicit none
	integer itsfnm
	character*(*) name
	itsfnm=para_get_id(name,' ')
	end

	function iscpar(name,section)
	implicit none
	integer iscpar
	character*(*) name,section
	iscpar = para_get_id(name,section)
	end

	function iscfnm(name,section)
	implicit none
	integer iscfnm
	character*(*) name,section
	iscfnm = para_get_id(name,section)
	end

!***********************************************************

	subroutine sctpar(section)
	implicit none
	character*(*) section
	call para_set_default_section(section)
	end

	subroutine sctfnm(section)
	implicit none
	character*(*) section
	call para_set_default_section(section)
	end

!**************************************************************

	subroutine get_sect_of(name,section)
	implicit none
	character*(*) name,section
	call para_find_section_to_name(name,section)
	end

!**************************************************************
!**************************************************************
!**************************************************************
!**************************************************************
!**************************************************************
!**************************************************************
!**************************************************************
!**************************************************************

	subroutine parinfo(iunit)

! prints info on parameter values

        use utility

	implicit none

	integer iunit

	character*6 name,section
	character*80 string
	double precision dvalue
	double precision value
	integer itype,nlen
	integer id,idfill

	integer check_entry_par

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

!**************************************************************

        subroutine pripar(iunit)

        implicit none

        integer iunit

        !call pripar1(iunit)
        call pripar3(iunit)
        !call pripar4(iunit)

        end

!*****************************************************

        subroutine pripar1(iunit)

! prints parameter values (1 column)

        implicit none

        integer iunit

        character*6 name,section
	character*80 string
        integer idfill,id,itype
        double precision dvalue
        double precision value

	call para_get_fill(idfill)

        do id=1,idfill
	  call para_get_info(id,name,section,itype,dvalue,string)
	  value = dvalue
	  if(itype.eq.1) then
            write(iunit,*) name,'  ',value
	  end if
        end do

	end

!*****************************************************

        subroutine pripar3(iunit)

! prints parameter values (3 columns)

        implicit none

        integer iunit

        character*6 name,section
        character*8 nameii(3)
	character*80 string
        integer idfill,id,itype,ii
        double precision dvalue
        double precision valueii(3)
        double precision value

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

!**************************************************************

	subroutine pripar4(iunit)

! prints parameter values (4 columns)

	implicit none

	integer iunit

	logical bflag
	character*80 line	!BUGFIX (was 79)

	integer npara,imod,i,itype
	integer ianf,iend
	integer itspar,infpar
	integer id,idfill
	character*6 name,section
	character*80 string
	double precision dvalue
	double precision value
	double precision flag,absval

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
            write(line(ianf:iend),2345) name,' =',nint(value),'  '
	  else
	    if(absval.lt.1000.and.absval.ge.0.01.or.absval.eq.0.) then
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

!**********************************************************

        subroutine prifnm(iunit)

! prints parameter values

        use utility

        implicit none

        integer iunit

        character*80 string
        character*6 name,section
        integer nlen,itype
	integer id,idfill
	double precision dvalue

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
 
!**********************************************************

	subroutine chkparam(iunit)

        use utility

	implicit none

	integer iunit

        character*80 string
        character*6 name,section
        integer nlen,itype
	integer id,idfill
	double precision dvalue
	double precision value

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

!**********************************************************

	subroutine check_parameter_values(text)

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

!**********************************************************

        function intpar(name)

! tests if name is integer
!
! name          name to test
! intpar        1 if name is integer, 0 if not

        use utility

        implicit none

        integer intpar
        character*(*) name

        integer i
        character*6 let
        character*1 namein
        data let /'ijklmn'/

        namein=name(1:1)
        call uplow(namein,'low')

        intpar=0
        do i=1,6
          if(namein.eq.let(i:i)) intpar=1
        end do

        end

!**********************************************************************

	subroutine getaz(azpar)

! returns actual az

	implicit none

	double precision azpar

	include 'femtime.h'

	double precision ampar

	azpar=getpar('azpar')
	ampar=0.			!dummy

	call changeimp(it,azpar,ampar)

	end

!**********************************************************************

	subroutine getazam(azpar,ampar)

! returns actual az,am

	implicit none

	double precision azpar
	double precision ampar

	include 'femtime.h'

	azpar=getpar('azpar')
	ampar=getpar('ampar')

	call changeimp(it,azpar,ampar)

	end

!**********************************************************************

	subroutine changeimp(it,azpar,ampar)

! changes parameters for semi-implicit time-step if necessary

	implicit none

	integer it
	double precision azpar,ampar

	include 'semi.h'

	call impini

	if( it .le. itimpl ) then
	  azpar = weight
	  ampar = weight
	end if

	end

!**********************************************************************

	subroutine impini

! initializes parameters for semi-implicit time-step

	implicit none

	include 'semi.h'

	include 'femtime.h'

	integer icall
	save icall
	data icall / 0 /

	if( icall .gt. 0 ) return

	icall = 1
	weight = 0.5
	itimpl = itanf

	end


!!**********************************************************************

	subroutine setimp(it,aweigh)

! sets parameters for semi-implicit time-step

	implicit none

	integer it
	double precision aweigh

	include 'semi.h'

	call impini

	itimpl = it
	weight = aweigh

	write(6,*) 'implicit parameters changed: ',itimpl,weight

	end

!**********************************************************************

	function getimp()

! gets weight for semi-implicit time-step

	implicit none

	double precision getimp

	include 'semi.h'

	getimp = weight

	end

!**********************************************************************

	function bimpli(it)

! checks if semi-implicit time-step is active

	implicit none

	logical bimpli
	integer it

	include 'semi.h'

	call impini

	if( it .le. itimpl ) then
	   bimpli = .true.
	else
	   bimpli = .false.
	end if

	end

!**********************************************************************
!==================================================================
	end module para
!==================================================================
!*****************************************************
!	program para_main
!	call test_par
!	end
!*****************************************************

