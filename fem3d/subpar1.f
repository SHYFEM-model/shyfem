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
c
c**************************************************************
c**************************************************************
c**************************************************************

	subroutine init_par

c initializes parameter structure

	implicit none

	include 'subpar.h'

	nentry = 0		!total number of values
	inarr = 0		!total number of arrays
	incha = 0		!total number of char strings
	infarr = 0		!filling of arrays
	!infcha = 0		!filling of chars

	actsec = ' '
	defsec = 'para'

	end

c**************************************************************

	blockdata subpar_blockdata

	implicit none

	include 'subpar.h'

	data nentry /0/
	data inarr /0/
	data incha /0/
	data infarr /0/
	!data infcha /0/
	data actsec /' '/
	data defsec /'para'/

	end

c**************************************************************
c**************************************************************
c**************************************************************

	function find_entry_par(name,sect)

c finds entry for name (internal routine)
c
c returns index where name has been found (negative if none found)
c if section is given, looks for it, else takes first good name

	implicit none

	include 'subpar.h'

	integer find_entry_par
	character*(*) name,sect

	logical bsect
	integer i,ientries
	character*6 namvec,actname,actsect

	integer ichanm

	actname = name
	actsect = sect

	bsect = sect .ne. ' '
	find_entry_par = 0

	!if( name .eq. 'itanf' ) then
	!  write(6,*) '*** looking for itanf: ',nentry
	!end if
	!  write(6,*) '*** looking for param: ',nentry,actname,ichanm(name)

	do i=1,nentry
	  namvec=nampar(i)
	  auxsec=secpar(i)
	  if( actname .eq. namvec ) then
	    if( .not. bsect ) then
	      find_entry_par = i
	      return
	    else if( auxsec .eq. sect ) then
	      find_entry_par = i
	      return
	    else
	      find_entry_par = find_entry_par - 1
	    end if
	  end if
	end do

	end

c***********************************************************

	function check_entry_par(name,sect,itype)

c finds and checks entry for name (internal routine)
c
c checks if type is compatible

	implicit none

	include 'subpar.h'

	integer check_entry_par
	character*(*) name,sect
	integer itype

	integer i,ity
	integer find_entry_par

	i = find_entry_par(name,sect)
	if( i .le. 0 ) then
	  write(6,*) 'cannot find parameter: ',name,' in section ',sect
	  stop 'error stop check_entry_par: no such parameter'
	end if

	ity = itypar(i)
	if( ity .ne. itype ) then
	  write(6,*) 'wrong type of parameter of ',name
	  write(6,*) 'expecting ',itype,' but finding ',ity
	  stop 'error stop check_entry_par: wrong type'
	end if

	check_entry_par = i

	end

c***********************************************************

	function add_entry_par(name,sect,itype)

c finds place where to enter new value

	implicit none

	include 'subpar.h'

	integer add_entry_par
	character*(*) name,sect
	integer itype

	integer i
	integer find_entry_par

	i = find_entry_par(name,' ')
	if( i .gt. 0 ) then
	  write(6,*) 'parameter already defined: ',name,' section ',sect
	  stop 'error stop add_entry_par: parameter already defined'
	end if

	if( nentry .ge. nnamdi ) then
	  write(6,*) 'dimension error in nnamdi: ',nnamdi
	  call parinfo(6)
	  stop 'error stop add_entry_par: nnamdi'
	end if

	nentry = nentry + 1
	nampar(nentry) = name
	secpar(nentry) = actsec
	itypar(nentry) = itype

	add_entry_par = nentry

	end

c***********************************************************
c***********************************************************
c***********************************************************

	function getpar(name)
	implicit none
	real getpar
	character*(*) name
	double precision dgetpar
	getpar = dgetpar(name)
	end

c***********************************************************

	function dgetpar(name)

c gets value of parameter name

	implicit none

	include 'subpar.h'

	double precision dgetpar
	character*(*) name

	integer i
	integer check_entry_par

	i = check_entry_par(name,' ',1)

	!if( i .le. 0 ) then
	!  write(6,*) '*** cannot find: ',name,i,nentry
	!end if

	dgetpar = valpar(i)

	end

c***********************************************************

	subroutine getfnm(name,text)

c gets text of parameter name

	implicit none

	include 'subpar.h'

	character*(*) name,text

	integer i
	integer check_entry_par

	i = check_entry_par(name,' ',3)
	call copy_to_text(i,text)

	end

c**************************************************************
c**************************************************************
c**************************************************************

	subroutine putpar(name,value)
	implicit none
	character*(*) name
	real value
	double precision dvalue
	dvalue = value
	call dputpar(name,dvalue)
	end

c**************************************************************

	subroutine dputpar(name,dvalue)

c puts parameter dvalue in name

	implicit none

	include 'subpar.h'

	character*(*) name
	double precision dvalue

	integer i,ls
	integer check_entry_par
	integer ichanm

	i = check_entry_par(name,' ',1)
	!if( i .le. 0 ) then
	!  write(6,*) '*** cannot find: ',name,i,nentry
	!end if
	valpar(i) = dvalue

	!ls = max(6,ichanm(name))
	!ls = min(ls,len(name))
	!write(15,*) 'putpar: ',name(1:ls),'  ',i,dvalue

	end

c**************************************************************

	subroutine putfnm(name,text)

c puts text into parameter name

	implicit none

	include 'subpar.h'

	character*(*) name,text

	integer i
	integer check_entry_par

	i = check_entry_par(name,' ',3)
	call copy_from_text(i,text)

	end

c**************************************************************
c**************************************************************
c**************************************************************

	subroutine addpar(name,value)
	implicit none
	character*(*) name
	real value
	double precision dvalue
	dvalue = value
	call daddpar(name,dvalue)
	end

c**************************************************************

	subroutine daddpar(name,dvalue)

c adds parameter value to name

	implicit none

	include 'subpar.h'

	character*(*) name
	double precision dvalue

	integer i
	integer add_entry_par

	i = add_entry_par(name,' ',1)

	valpar(i) = dvalue

	end

c**************************************************************

	subroutine addfnm(name,text)

c adds text into parameter name

	implicit none

	include 'subpar.h'

	character*(*) name,text

	integer i
	integer add_entry_par

	i = add_entry_par(name,' ',3)
	incha = incha + 1
	if( incha .gt. nichdi ) then
	  stop 'error stop addfnm: dimension error in nichdi'
	end if
	valpar(i) = incha
	!ip_chapar(1,incha) = 0
	!ip_chapar(2,incha) = 0
	call copy_from_text(i,text)

	end

c**************************************************************
c**************************************************************
c**************************************************************

	subroutine copy_to_text(ientry,text)

	implicit none

	include 'subpar.h'

	integer ientry
	character*(*) text

	integer j

	j = nint(valpar(ientry))
	text = ' '
	text = chapar(j)

	end

	subroutine copy_from_text(ientry,text)

	implicit none

	include 'subpar.h'

	integer ientry
	character*(*) text

	integer j

	j = nint(valpar(ientry))
	chapar(j) = ' '
	chapar(j) = text

	end

c**************************************************************

	subroutine copy_to_text0(ientry,text)

c copies character string "ientry" to text

	implicit none

	include 'subpar.h'

	integer ientry
	character*(*) text

	integer ianf,iend,i,j,ip,ls,lt
	character*80 string

	integer ichanm

	ip = nint(valpar(ientry))

	!ianf = ip_chapar(1,ip)
	!iend = ip_chapar(2,ip)

	string = ' '
	j = 0
	do i=ianf,iend
	  j = j + 1
	  string(j:j) = chapar(i)
	end do
	text = string

	ls = ichanm(string)
	lt = len(text)
	if( lt .lt. ls ) then
	  write(6,*) 'not enough space to copy parameter to text string'
	  write(6,*) 'name: ',nampar(i)
	  write(6,*) 'string: ',string
	  stop 'error stop copy_to_string: text too short'
	end if

	end

c**************************************************************

	subroutine copy_from_text0(ientry,text)

c copies text to character string "ientry"

	implicit none

	include 'subpar.h'

	integer ientry
	character*(*) text

	integer ianf,iend,j,ip,ls,le
	character*80 string

	integer ichanm

	ip = nint(valpar(ientry))

	!ianf = ip_chapar(1,ip)
	!iend = ip_chapar(2,ip)

	string = text
	ls = max(1,ichanm(string))
	le = iend-ianf+1

	if( iend .le. 0 .or. ls .gt. le ) then
	  !ianf = infcha+1
	  !iend = infcha+ls
	  !infcha = infcha + ls
	  !if( infcha .gt. nchadi ) then
	  !  write(6,*) 'no space to put new string: ',infcha
	  !  stop 'error stop copy_from_text: no space for string'
	  !end if
	else
	  iend = ianf+ls-1
	end if

	!ip_chapar(1,ip) = ianf
	!ip_chapar(2,ip) = iend

	j = 0
	do ip=ianf,iend
	  j = j + 1
	  !chapar(ip) = text(j:j)
	end do

	end

c**************************************************************

	subroutine delete_section(sect)

c deletes whole section -> works only if section is last

	implicit none

	include 'subpar.h'

	character*(*) sect

	integer i,istart,ity
	integer ip,ianf,iend,imax
	character*6 name

c-----------------------------------------------
c find first entry for section sect -> istart
c-----------------------------------------------

	istart = 0
	do i=1,nentry
	  if( sect .eq. secpar(i) .and. istart .eq. 0 ) istart = i
	  if( sect .ne. secpar(i) .and. istart .gt. 0 ) goto 99
	end do
	if( istart .eq. 0 ) return		!nothing to delete

c-----------------------------------------------
c now delete parameters, arrays and strings
c-----------------------------------------------

	do i=nentry,istart,-1
	  name = nampar(i)
	  ity = itypar(i)
	  ip = nint(valpar(i))
	  if( ity .eq. 2 ) then			!numeric array
	    if( ip .ne. inarr ) goto 98
	    inarr = inarr - 1
	    !ianf = ip_arrpar(1,ip)
	    !iend = ip_arrpar(2,ip)
	    !if( iend .ne. infarr ) goto 97
	    !infarr = ianf - 1
	  else if( ity .eq. 3 ) then		!string
	    if( ip .ne. incha ) goto 96
	    incha = incha - 1
	    !ianf = ip_chapar(1,ip)
	    !iend = ip_chapar(2,ip)
	    !if( iend .ne. infcha ) goto 95
	    !infcha = ianf - 1
	  end if
	  nampar(i) = 'delete'
	  secpar(i) = 'delete'
	  itypar(i) = -1
	  valpar(i) = -999.
	end do

c-----------------------------------------------
c check new filling of arrays and strings
c-----------------------------------------------

	imax = 0
	do ip=1,inarr
	  iend = ip_arrpar(2,ip)
	  imax = max(imax,iend)
	end do
	!write(15,*) 'delete_section (array): ',infarr,imax
	infarr = imax

	!imax = 0
	!do ip=1,incha
	!  !iend = ip_chapar(2,ip)
	!  imax = max(imax,iend)
	!end do
	!write(15,*) 'delete_section (string): ',infcha,imax
	!infcha = imax

c-----------------------------------------------
c adjourn total number of parameters
c-----------------------------------------------

	nentry = istart - 1

	!write(15,*) 'delete_section: ',nentry,inarr,incha,infarr,infcha

c-----------------------------------------------
c end of routine
c-----------------------------------------------

	return
   95	continue
	write(6,*) 'deleting ',name
	write(6,*) i,ip,incha,ianf,iend!,infcha
	stop 'error stop delete_section: string'
   96	continue
	write(6,*) 'deleting ',name
	write(6,*) i,ip,incha
	stop 'error stop delete_section: string'
   97	continue
	write(6,*) 'deleting ',name
	write(6,*) i,ip,inarr,ianf,iend,infarr
	stop 'error stop delete_section: numeric array'
   98	continue
	write(6,*) 'deleting ',name
	write(6,*) i,ip,inarr
	stop 'error stop delete_section: numeric array'
   99	continue
	write(6,*) 'section must be the last one inserted'
	write(6,*) i,nentry,istart,sect,secpar(i)
	stop 'error stop delete_section: not last section'
	end

c**************************************************************
c**************************************************************
c**************************************************************

	function itspar(name)

c tests if name is available (itspar =|= 0) or not (itspar=0)

	implicit none

	include 'subpar.h'

	integer itspar
	character*(*) name

	integer find_entry_par

	itspar = find_entry_par(name,' ')
	if( itspar .lt. 0 ) itspar = 0

	end

c**************************************************************

	function iscpar(name,sect)

c tests if name is in section (iscpar > 0) or not (iscpar=0)

	implicit none

	include 'subpar.h'

	integer iscpar
	character*(*) name,sect

	integer find_entry_par

	iscpar = find_entry_par(name,sect)
	if( iscpar .lt. 0 ) iscpar = 0

	end

c**************************************************************

	subroutine get_sect_of(name,sect)

c returns section of parameter

	implicit none

	include 'subpar.h'

	character*(*) name,sect

	integer i

	integer find_entry_par

	sect = ' '
	i = find_entry_par(name,' ')
	if( i .gt. 0 ) sect = secpar(i)

	end

c***********************************************************

	subroutine sctpar(sect)

c sets actual name of section

	implicit none

	include 'subpar.h'

	character*(*) sect

	actsec = sect

	end

c**************************************************************
c**************************************************************
c**************************************************************

	function itsfnm(name)
	implicit none
	integer itsfnm
	character*(*) name
	integer itspar
	itsfnm = itspar(name)
	end

	function iscfnm(name,sect)
	implicit none
	integer iscfnm
	character*(*) name,sect
	integer iscpar
	iscfnm = iscpar(name,sect)
	end

	subroutine sctfnm(sect)
	implicit none
	character*(*) sect
	call sctpar(sect)
	end

c**************************************************************
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

	implicit none

	include 'subpar.h'

	integer iunit
	character*6 name,sect
	character*80 text
	real value
	integer itype,i,j,nlen

	integer check_entry_par,ichanm

	write(iunit,*) 'parinfo:'
	write(iunit,*) 'nentry: ',nentry,' of possible ',nnamdi
	write(iunit,*) 'incha:  ',incha, ' of possible ',nichdi

	do i=1,nentry
          name = nampar(i)
          sect = secpar(i)
          value = valpar(i)
          itype = itypar(i)
	  write(iunit,*) i,itype,value,name,'  ',sect
	  if( itype .eq. 3 ) then
	    call copy_to_text(i,text)
            nlen=max(1,ichanm(text))
            write(iunit,*) '    ',nlen,text(1:nlen)
	  end if
	end do

	end

c**************************************************************

	subroutine pripar(iunit)

c prints parameter values

	implicit none

	include 'subpar.h'

	integer iunit

	logical bflag
	character*80 line	!BUGFIX (was 79)
	character*6 name

	integer npara,imod,i,itype
	integer ianf,iend
	integer itspar,infpar,intpar
	real getpar
	real value
	real flag,absval

	line=' '
	bflag=.false.
	if(itspar('flag').ne.0) then
		bflag=.true.
		flag=getpar('flag')
	end if
	npara=nentry

	imod=0
	do i=1,npara
	  name = nampar(i)
	  value = valpar(i)
	  itype = itypar(i)
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

        implicit none

        integer iunit

	include 'subpar.h'

        character*80 text
        character*6 name
        integer i,j,nlen,itype

	integer check_entry_par
        integer ichanm

        do i=1,nentry
	  name = nampar(i)
	  itype = itypar(i)
          if(itype.eq.3) then
	    j = check_entry_par(name,' ',3)
	    call copy_to_text(j,text)
            nlen=max(1,ichanm(text))
            write(iunit,2345) i,nlen,name,auxsec,text(1:nlen)
          end if
        end do

        return
 2345   format(1x,2i4,2(1x,a6,1x),3x,a)
        end
 
c**********************************************************

	subroutine chkparam(iunit)

	implicit none

	include 'subpar.h'

	integer iunit

        character*80 text
        character*6 name,sect
        integer i,j,nlen,itype
	integer ip,ianf,iend
	real value

        integer ichanm
	integer check_entry_par

        do i=1,nentry
	  name = nampar(i)
	  sect = secpar(i)
	  value = valpar(i)
	  itype = itypar(i)
          if(itype.eq.1) then
            write(iunit,2345) i,name,sect,itype,value
          else if(itype.eq.3) then
	    j = check_entry_par(name,' ',3)
	    call copy_to_text(j,text)
            nlen=max(1,ichanm(text))
	    ip = nint(value)
	    !ianf = ip_chapar(1,ip)
	    !iend = ip_chapar(2,ip)
            !write(iunit,2346) i,name,sect,itype
!     +				,ip,ianf,iend,nlen,text(1:nlen)
	  else
            write(iunit,2345) i,name,sect,itype,value
          end if
        end do

        return
 2345   format(1x,i4,2(1x,a6,1x),i4,e12.4,i4,1x,a)
 2346   format(1x,i4,2(1x,a6,1x),i4,4i4,1x,a)
	end

c**********************************************************

	subroutine check_parameter_values(text)

	implicit none

	character*(*) text

	include 'subpar.h'

	integer iunit

	iunit = 15

	return

	write(iunit,*) '--------------------------------'
	write(iunit,*) 'info on parameters: ',nentry,'  ',text
	call chkparam(iunit)
	write(iunit,*) '--------------------------------'
	!write(iunit,*) '...printing with pripar...'
        !call pripar(iunit)
	!write(iunit,*) '...printing with prifnm...'
        !call prifnm(iunit)
	!write(iunit,*) '...end of printing in check_parameter_values'

	end

c**********************************************************
c
c	function ichanm(text)
c	integer ichanm
c	character*(*) text
c	ichanm=0
c	end
c
c	subroutine test_par
c	end
c
c	program main_test_par
c	call test_par
c	end
c
c**********************************************************

c************************************************************************

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
        call uplow(namein,'low')

        intpar=0
        do i=1,6
          if(namein.eq.let(i:i)) intpar=1
        end do

        end

c*****************************************************

