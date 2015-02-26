c Parameter management routines taken from subpar3.f
c changing the module name to para_aqua and name length
c
c 
c
c parameter management routines
c
c contents :
c
! function para_get_id(name,section)
! subroutine para_init_alloc
! subroutine para_init_new_id(id)
! subroutine para_init_id(id)
! subroutine para_set_default_section(section)
! function para_has_name(name)
! function para_is_name_in_section(name,section)
! subroutine para_find_section_to_name(name,section)
! subroutine para_delete_section(section)
! subroutine para_delete_name(name)
! subroutine para_delete_id(id)
! subroutine para_get_fill(idfill)
! subroutine para_get_info(id,name,section,itype,value,string)
! subroutine para_get_value(name,value)
! subroutine para_add_value(name,value)
! subroutine para_put_value(name,value)
! subroutine para_get_string(name,string)
! subroutine para_add_string(name,string)
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
c
c**************************************************************
c**************************************************************
c**************************************************************

!==================================================================
	module para_aqua
!==================================================================

! itype =  1	value
! itype = 10	string

	implicit none

	type, private :: entry

	  character*80 :: name
	  character*6 :: section
	  integer :: isize
	  integer :: itype
	  double precision :: value
	  double precision, allocatable :: array(:)
	  character*80 :: string

	end type entry

	integer, save, private :: idlast = 0
	!integer, parameter, private :: ndim = 300
	!type(entry), save :: pentry(ndim)
	integer, save, private :: ndim = 0
	type(entry), save, allocatable :: pentry(:)

	character*(6), save, private :: def_section = ' '

!==================================================================
	contains
!==================================================================




	function para_get_id(name,section)

	integer para_get_id
	character*(*) name,section

	integer id,ifound
	character*80 namarg
	character*6 s
    
        namarg=name
        !call uplow(namarg,'low')
        s=section
        !call uplow(s,'low')

	ifound = 0
	para_get_id = 0

	do id=1,idlast
	  !write(6,*) trim(namarg),trim(pentry(id)%name)
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
	  call para_info_id(id_opt)
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
	character*80 name,section,string

	name = pentry(id)%name
	section = pentry(id)%section
	itype = pentry(id)%itype
	value = pentry(id)%value
	string = pentry(id)%string

	!write(6,*) id,itype,trim(name),value,trim(string)

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
	
	!print *, id, trim(name)
	
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
	pentry(id)%string = string
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

	pentry(id)%string = string

	end subroutine para_put_string

!==================================================================
	end module para_aqua
!==================================================================


c**************************************************************
c**************************************************************
c**************************************************************


c***********************************************************
c***********************************************************
c***********************************************************

! All these commented routines are located in subpara3.f:

! 	function getpar(name)
! 	use para_aqua
! 	implicit none
! 	real getpar
! 	character*(*) name
! 	double precision value
! 	call para_get_value(name,value)
! 	getpar = value
! 	end
! 
! 	function dgetpar(name)
! 	use para_aqua
! 	implicit none
! 	double precision dgetpar
! 	character*(*) name
! 	double precision value
! 	call para_get_value(name,value)
! 	dgetpar = value
! 	end
! 
! 	subroutine getfnm(name,string)
! 	use para_aqua
! 	implicit none
! 	character*(*) name,string
! 	call para_get_string(name,string)
! 	end
! 
! c**************************************************************
! 
! 	subroutine putpar(name,value)
! 	use para_aqua
! 	implicit none
! 	character*(*) name
! 	real value
! 	double precision dvalue
! 	dvalue = value
! 	call para_put_value(name,dvalue)
! 	end
! 
! 	subroutine dputpar(name,value)
! 	use para_aqua
! 	implicit none
! 	character*(*) name
! 	double precision value
! 	call para_put_value(name,value)
! 	end
! 
! 	subroutine putfnm(name,string)
! 	use para_aqua
! 	implicit none
! 	character*(*) name,string
! 	call para_put_string(name,string)
! 	end
! 
! c**************************************************************
! 
! 	subroutine addpar(name,value)
! 	use para_aqua
! 	implicit none
! 	character*(*) name
! 	real value
! 	double precision dvalue
! 	dvalue = value
! 	call para_add_value(name,dvalue)
! 	end
! 
! 	subroutine daddpar(name,value)
! 	use para_aqua
! 	implicit none
! 	character*(*) name
! 	double precision value
! 	call para_add_value(name,value)
! 	end
! 
! 	subroutine addfnm(name,string)
! 	use para_aqua
! 	implicit none
! 	character*(*) name,string
! 	call para_add_string(name,string)
! 	end
! 
! c**************************************************************
! c**************************************************************
! c**************************************************************
! 
! 
! c**************************************************************
! c**************************************************************
! c**************************************************************
! 
! 	subroutine delete_section(section)
! 	use para_aqua
! 	implicit none
! 	character*(*) section
! 	call para_delete_section(section)
! 	end
! 
! c**************************************************************
! c**************************************************************
! c**************************************************************
! 
! 	function itspar(name)
! 	use para_aqua
! 	implicit none
! 	integer itspar
! 	character*(*) name
! 	itspar=para_get_id(name,' ')
! 	end
! 
! 	function itsfnm(name)
! 	use para_aqua
! 	implicit none
! 	integer itsfnm
! 	character*(*) name
! 	itsfnm=para_get_id(name,' ')
! 	end
! 
! 	function iscpar(name,section)
! 	use para_aqua
! 	implicit none
! 	integer iscpar
! 	character*(*) name,section
! 	iscpar = para_get_id(name,section)
! 	end
! 
! 	function iscfnm(name,section)
! 	use para_aqua
! 	implicit none
! 	integer iscfnm
! 	character*(*) name,section
! 	iscfnm = para_get_id(name,section)
! 	end
! 
! c***********************************************************
! 
! 	subroutine sctpar(section)
! 	use para_aqua
! 	implicit none
! 	character*(*) section
! 	call para_set_default_section(section)
! 	end
! 
! 	subroutine sctfnm(section)
! 	use para_aqua
! 	implicit none
! 	character*(*) section
! 	call para_set_default_section(section)
! 	end
! 
! c**************************************************************
! 
! 	subroutine get_sect_of(name,section)
! 	use para_aqua
! 	implicit none
! 	character*(*) name,section
! 	call para_find_section_to_name(name,section)
! 	end
! 
! c**************************************************************
! c**************************************************************
! c**************************************************************
! c**************************************************************
! c**************************************************************
! c**************************************************************
! c**************************************************************
! c**************************************************************
! 
! 	subroutine parinfo(iunit)
! 
! c prints info on parameter values
! 
! 	use para_aqua
! 
! 	implicit none
! 
! 	integer iunit
! 
! 	character*6 name,section
! 	character*80 string
! 	double precision dvalue
! 	real value
! 	integer itype,nlen
! 	integer id,idfill
! 
! 	integer check_entry_par,ichanm
! 
! 	call para_get_fill(idfill)
! 
! 	write(iunit,*) 'parinfo: ',idfill
! 
!         do id=1,idfill
! 	  call para_get_info(id,name,section,itype,dvalue,string)
! 	  value = dvalue
! 	  write(iunit,*) id,itype,value,name,'  ',section
! 	  if( itype .eq. 3 ) then
!             nlen=max(1,ichanm(string))
!             write(iunit,*) '    ',nlen,string(1:nlen)
! 	  end if
! 	end do
! 
! 	end
! 
! c**************************************************************
! 
!         subroutine pripar(iunit)
! 
!         implicit none
! 
!         integer iunit
! 
!         call pripar1(iunit)
!         !call pripar4(iunit)
! 
!         end
! 
! c*****************************************************
! 
!         subroutine pripar1(iunit)
! 
! c prints parameter values
! 
! 	use para_aqua
! 
!         implicit none
! 
!         integer iunit
! 
!         character*6 name,section
! 	character*80 string
!         integer idfill,id,itype
!         double precision dvalue
!         real value
! 
! 	call para_get_fill(idfill)
! 
!         do id=1,idfill
! 	  call para_get_info(id,name,section,itype,dvalue,string)
! 	  value = dvalue
! 	  if(itype.eq.1) then
!             write(iunit,*) name,'  ',value
! 	  end if
!         end do
! 
! 	end
! 
! c**************************************************************
! 
! 	subroutine pripar4(iunit)
! 
! c prints parameter values
! 
! 	use para_aqua
! 
! 	implicit none
! 
! 	integer iunit
! 
! 	logical bflag
! 	character*80 line	!BUGFIX (was 79)
! 
! 	integer npara,imod,i,itype
! 	integer ianf,iend
! 	integer itspar,infpar,intpar
! 	integer id,idfill
! 	character*6 name,section
! 	character*80 string
! 	double precision dvalue
! 	real value
! 	real flag,absval
! 
! 	line=' '
! 	bflag=.false.
! 
! 	call para_get_fill(idfill)
! 	npara=idfill
! 
! 	imod=0
! 	do id=1,npara
! 	  call para_get_info(id,name,section,itype,dvalue,string)
! 	  value = dvalue
! 	  if(itype.ne.1) goto 1
! 	  if(bflag.and.value.eq.flag) goto 1
! 	  imod=imod+1
! 	  ianf=20*(imod-1)+1
! 	  iend=20*imod
! 	  absval=abs(value)
! 	  if(intpar(name).eq.1) then
!             write(line(ianf:iend),2345)
!      +                  name,' =',nint(value),'  '
! 
! 	  else
! 	    if(absval.lt.1000.and.absval.ge.0.01
!      +                        .or.absval.eq.0.) then
! 	      write(line(ianf:iend),2347) name,' =',value,'  '
! 	    else
! 	      write(line(ianf:iend),2346) name,' =',value,'  '
! 	    end if
! 	  end if
! 	  if(imod.eq.4) then
! 		write(iunit,*) line(1:79)
! 		line=' '
! 		imod=0
! 	  end if
!     1     continue
! 	end do
! 
! 	if(imod.ne.4) write(iunit,*) line
! 
! 	return
!  2345   format(a6,a2,i10,a2)
!  2346   format(a6,a2,e10.2,a2)
!  2347   format(a6,a2,f10.3,a2)
! 	end
! 
! c**********************************************************
! 
!         subroutine prifnm(iunit)
! 
! c prints parameter values
! 
! 	use para_aqua
! 
!         implicit none
! 
!         integer iunit
! 
!         character*80 string
!         character*6 name,section
!         integer nlen,itype
! 	integer id,idfill
! 	double precision dvalue
! 
!         integer ichanm
! 
! 	call para_get_fill(idfill)
! 
!         do id=1,idfill
! 	  call para_get_info(id,name,section,itype,dvalue,string)
!           if(itype.eq.3) then
!             nlen=max(1,ichanm(string))
!             write(iunit,2345) id,nlen,name,section,string(1:nlen)
!           end if
!         end do
! 
!         return
!  2345   format(1x,2i4,2(1x,a6,1x),3x,a)
!         end
!  
! c**********************************************************
! 
! 	subroutine chkparam(iunit)
! 
! 	use para_aqua
! 
! 	implicit none
! 
! 	integer iunit
! 
!         character*80 string
!         character*6 name,section
!         integer nlen,itype
! 	integer id,idfill
! 	double precision dvalue
! 	real value
! 
!         integer ichanm
! 
! 	call para_get_fill(idfill)
! 
!         do id=1,idfill
! 	  call para_get_info(id,name,section,itype,dvalue,string)
! 	  value = dvalue
!           if(itype.eq.1) then
!             write(iunit,2345) id,name,section,itype,value
!           else if(itype.eq.3) then
!             write(iunit,2346) id,name,section,itype,trim(string)
! 	  else
! 	    write(6,*) 'error itype...'
!             write(iunit,2345) id,name,section,itype,value
! 	    stop 'error stop chkparam'
!           end if
!         end do
! 
!         return
!  2345   format(1x,i4,2(1x,a6,1x),i4,e12.4)
!  2346   format(1x,i4,2(1x,a6,1x),i4,a)
! 	end
! 
! c**********************************************************
! 
! 	subroutine check_parameter_values(text)
! 
! 	use para_aqua
! 
! 	implicit none
! 
! 	character*(*) text
! 
! 	integer iunit,idfill
! 
! 	iunit = 117
! 
! 	!return
! 
! 	call para_get_fill(idfill)
! 
! 	write(iunit,*) '--------------------------------'
! 	write(iunit,*) 'info on parameters: ',text
! 	call chkparam(iunit)
! 	write(iunit,*) '--------------------------------'
! 	write(iunit,*) '...printing with pripar...'
!         call pripar(iunit)
! 	write(iunit,*) '...printing with prifnm...'
!         call prifnm(iunit)
! 	write(iunit,*) '...end of printing in check_parameter_values'
! 
! 	end
! 
! c**********************************************************
! 
!         function intpar(name)
! 
! c tests if name is integer
! c
! c name          name to test
! c intpar        1 if name is integer, 0 if not
! 
!         implicit none
! 
!         integer intpar
!         character*(*) name
! 
!         integer i
!         character*6 let
!         character*1 namein
!         data let /'ijklmn'/
! 
!         namein=name(1:1)
!         call uplow(namein,'low')
! 
!         intpar=0
!         do i=1,6
!           if(namein.eq.let(i:i)) intpar=1
!         end do
! 
!         end
! 
! c*****************************************************
! c	program para_main
! c	call test_par
! c	end
! c*****************************************************

