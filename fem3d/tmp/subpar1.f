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
c
c**************************************************************
c**************************************************************
c**************************************************************

	subroutine init_par

	implicit none

	include 'subpar.h'

	nentry = 0
	inarr = 0
	incha = 0
	infarr = 0
	infcha = 0

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
	data infcha /0/
	data actsec /' '/
	data defsec /'para'/

	end

c**************************************************************
c**************************************************************
c**************************************************************

	function find_entry_par(name,sect)

c finds entry for name (internal routine)

	implicit none

	include 'subpar.h'

	integer find_entry_par
	character*(*) name,sect

	logical bsect
	integer i,ientries
	character*6 namvec

	bsect = sect .ne. ' '
	find_entry_par = 0

	do i=1,nentry
	  namvec=nampar(i)
	  auxsec=secpar(i)
	  if( name .eq. namvec ) then
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

	nentry = nentry + 1
	if( nentry .gt. nnamdi ) then
	  write(6,*) 'dimension error in nnamdi: ',nnamdi
	  stop 'error stop add_entry_par: nnamdi'
	end if

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

	integer i
	integer check_entry_par

	i = check_entry_par(name,' ',1)

	valpar(i) = dvalue

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
	valpar(i) = incha
	ip_chapar(1,incha) = 0
	ip_chapar(2,incha) = 0
	call copy_from_text(i,text)

	end

c**************************************************************
c**************************************************************
c**************************************************************

	subroutine copy_to_text(ientry,text)

c copies character string "i" to text

	implicit none

	include 'subpar.h'

	integer ientry
	character*(*) text

	integer ianf,iend,i,j,ip,ls,lt
	character*80 string

	integer ichanm

	ip = nint(valpar(ientry))

	ianf = ip_chapar(1,ip)
	iend = ip_chapar(2,ip)

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

	subroutine copy_from_text(ientry,text)

c copies text to character string "i"

	implicit none

	include 'subpar.h'

	integer ientry
	character*(*) text

	integer ianf,iend,j,ip,ls,le
	character*80 string

	integer ichanm

	ip = nint(valpar(ientry))

	ianf = ip_chapar(1,ip)
	iend = ip_chapar(2,ip)

	string = text
	ls = ichanm(string)
	le = iend-ianf+1

	if( iend .le. 0 .or. ls .gt. le ) then
	  ianf = infcha+1
	  iend = infcha+ls
	  infcha = infcha + ls
	  if( infcha .gt. nchadi ) then
	    write(6,*) 'no space to put new string: ',infcha
	    stop 'error stop copy_from_text: no space for string'
	  end if
	else
	  iend = ianf+ls-1
	end if

	ip_chapar(1,ip) = ianf
	ip_chapar(2,ip) = iend

	j = 0
	do ip=ianf,iend
	  j = j + 1
	  chapar(ip) = text(j:j)
	end do

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

	subroutine pripar

c prints parameter values

	implicit none

	include 'subpar.h'

	logical bflag
	character*80 line	!BUGFIX (was 79)
	character*6 name

	integer npara,imod,i
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
	  if(bflag.and.value.eq.flag) goto 1
	  imod=imod+1
	  ianf=20*(imod-1)+1
	  iend=20*imod
	  absval=abs(value)
	  if(absval.lt.1000.and.absval.ge.0.01
     +                        .or.absval.eq.0.) then
	    write(line(ianf:iend),2347) name,' =',value,'  '
	  else
	    write(line(ianf:iend),2346) name,' =',value,'  '
	  end if
	  if(imod.eq.4) then
		write(6,*) line(1:79)
		line=' '
		imod=0
	  end if
    1     continue
	end do

	if(imod.ne.4) write(6,*) line

	return
 2346   format(a6,a2,e10.2,a2)
 2347   format(a6,a2,f10.3,a2)
	end

c**********************************************************

	function ichanm(text)
	integer ichanm
	character*(*) text
	ichanm=0
	end

	subroutine test_par
	end

	program main_test_par
	call test_par
	end

c**********************************************************

