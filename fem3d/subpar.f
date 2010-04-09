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

	function manpar(name,dvalue,id,what)

c internal routine to access par data
c
c name          name of parameter
c value         contents of parameter
c id            number of entry (only for what='l')
c what          what to do ?
c               g : return contents of name in value
c               t : name is defined ?
c               s : put section of name in auxsec
c               p : put value into name
c               a : add value with name
c               f : how many entries ?
c               l : list name and value of entry id
c		c : check (undocumented)
c manpar        entry number (0 if no entry name is found,
c               ...total number of entries for what='f')
c
c comments:
c
c - uses value in parsec for section of name
c
c problems:
c
c - there is no way to change section for an already given parameter

	implicit none

	integer manpar
	character*(*) name,what
	double precision dvalue
	integer id

c-------------------------------------------
	integer nnamdi
	parameter (nnamdi=200)
c-------------------------------------------
	double precision valpar(nnamdi) !array containing parameter values
	character*6 nampar(nnamdi)      !array containing names of parameters
	character*6 secpar(nnamdi)      !array containing section names
c-------------------------------------------
	character*6 parsec, auxsec
	common /parsco/ parsec, auxsec		!default section
c-------------------------------------------

	character*6 namarg,namvec
	character*1 type
	integer i,ientry

	integer ifirst,nentry
	save ifirst,nentry
	save valpar,nampar,secpar
	data ifirst,nentry / 0 , 0 /

	if(ifirst.eq.0) then
		parsec = ' '
		ifirst=1
	end if

	type=what(1:1)
	call uplow(type,'low')

	namarg=name
	call uplow(namarg,'low')

c find entry name %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	ientry=0
	do i=1,nentry
	   namvec=nampar(i)
	   !call uplow(namvec,'low')
	   if(namarg.eq.namvec) ientry=i
	end do

c see what to do %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if(type.eq.'f') then                    !filling of array
		manpar=nentry
		return
	else if(type.eq.'t') then               !test if name exists
		manpar=ientry
		return
	else if(type.eq.'s') then               !put section of name in auxsec
		auxsec = ' '
		if( ientry .gt. 0 ) then
		  auxsec=secpar(ientry)
		end if
		manpar=ientry
		return
	else if(type.eq.'l') then               !list id entry
		dvalue=0.
		if(id.le.0.or.id.gt.nentry) then
			manpar=0
			name=' '
			return
		else
			ientry=id
			manpar=ientry
			name=nampar(ientry)
		end if
	else if(type.eq.'g') then               !get parameter name
		dvalue=0.
		if(ientry.eq.0) then
			manpar=0
			return
		else
			manpar=ientry
		end if
	else if(type.eq.'a') then               !add name
		if(ientry.eq.0) then
			manpar=nentry+1
		else
			manpar=ientry
		end if
	else if(type.eq.'p') then               !put name
		if(ientry.eq.0) then
			manpar=0
			return
		else
			manpar=ientry
		end if
	else if(type.eq.'c') then		!not documented -> check
		if(id.le.0.or.id.gt.nentry) then
			manpar=0
			dvalue=0.
			name=' '
		else
			manpar = id
			dvalue = valpar(id)
			name = nampar(id)
			auxsec = secpar(id)
			write(6,'(i6,f10.4,2(1x,a6))') id,dvalue,name,auxsec
		end if
	else                                    !not recognized
		write(6,*) 'type not recognized : ',what
		stop 'error stop : manpar'
	end if

c read or write text %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if(type.eq.'g'.or.type.eq.'l') then             !get parameter
		dvalue=valpar(ientry)
	else if(type.eq.'a'.or.type.eq.'p') then        !put parameter
		if(ientry.eq.0) then
			nentry=nentry+1
			ientry=nentry
			if(ientry.gt.nnamdi) then
				write(6,*) 'dimension error nnamdi'
				stop 'error stop : manpar'
			end if
			nampar(ientry)=name
			secpar(ientry)=parsec
			valpar(ientry)=dvalue
		else                                    !put in old place
			valpar(ientry)=dvalue
		end if
	end if

	end

c***********************************************************

	function getpar(name)

c gets parameter name in value

	implicit none

	real getpar
	character*(*) name

	integer i,manpar
	double precision dvalue

	i=manpar(name,dvalue,0,'get')
	if(i.eq.0) then
		write(6,*) 'Parameter not found : ',name
		stop 'error stop : getpar'
	end if

	getpar = dvalue

	end

c***********************************************************

	function dgetpar(name)

c gets parameter name in value

	implicit none

	double precision dgetpar
	character*(*) name

	integer i,manpar
	double precision dvalue

	i=manpar(name,dvalue,0,'get')
	if(i.eq.0) then
		write(6,*) 'Parameter not found : ',name
		stop 'error stop : getpar'
	end if

	dgetpar = dvalue

	end

c***********************************************************

	function itspar(name)

c tests if name is available (itspar =|= 0) or not (itspar=0)

	implicit none

	integer itspar
	character*(*) name

	integer manpar
	double precision dvalue

	itspar=manpar(name,dvalue,0,'test')

	end

c***********************************************************

	function iscpar(name,sect)

c tests if name is in section (iscpar =|= 0) or not (iscpar=0)

	implicit none

	integer iscpar
	character*(*) name,sect

c-------------------------------------------
	character*6 parsec, auxsec
	common /parsco/ parsec, auxsec		!default section
c-------------------------------------------

	character*6 asect
	integer manpar
	double precision dvalue

c	next call puts section of name into auxsec

	iscpar=manpar(name,dvalue,0,'sect')

	if( iscpar .gt. 0 ) then
	    asect = sect
	    call uplow(asect,'low')
	    if( asect .ne. auxsec ) iscpar = 0
	end if

	end

c***********************************************************

	subroutine get_sect_of(name,sect)

c returns section of parameter

	implicit none

	character*(*) name,sect

c-------------------------------------------
	character*6 parsec, auxsec
	common /parsco/ parsec, auxsec		!default section
c-------------------------------------------

	integer iscpar,manpar
	double precision dvalue

c	next call puts section of name into auxsec

	iscpar=manpar(name,dvalue,0,'sect')

	if( iscpar .gt. 0 ) then
	  sect = auxsec
	else
	  sect = ' '
	end if

	end

c***********************************************************

	subroutine sctpar(sect)

c sets default name of section

	implicit none

	character*(*) sect

c-------------------------------------------
	character*6 parsec, auxsec
	common /parsco/ parsec, auxsec		!default section
c-------------------------------------------

	character*6 name
	integer idum,manpar
	double precision dvalue

c	next call is dummy call to be sure that parsco is initialized

	idum=manpar(name,dvalue,0,'fill')
	idum = 2 * idum	!to avoid compiler warnings

	auxsec = sect
	call uplow(auxsec,'low')
	parsec = auxsec

	end

c***********************************************************

	subroutine putpar(name,value)

c puts parameter value in name

	implicit none

	character*(*) name
	real value

	double precision dvalue
	integer i,manpar

	dvalue = value
	i=manpar(name,dvalue,0,'put')
	if(i.eq.0) then
		write(6,*) 'Parameter not found : ',name
		stop 'error stop : putpar'
	end if

	end

c***********************************************************

	subroutine dputpar(name,dvalue)

c puts parameter value in name (double precision version)

	implicit none

	character*(*) name
	double precision dvalue

	integer i,manpar

	i=manpar(name,dvalue,0,'put')
	if(i.eq.0) then
		write(6,*) 'Parameter not found : ',name
		stop 'error stop : putpar'
	end if

	end

c**************************************************************

	subroutine addpar(name,value)

c adds parameter value to name

	implicit none

	character*(*) name
	real value

	double precision dvalue
	integer idum,manpar

	dvalue = value
	idum=manpar(name,dvalue,0,'add')
	idum=idum*2

	end

c**************************************************************

	subroutine daddpar(name,dvalue)

c adds parameter value to name

	implicit none

	character*(*) name
	double precision dvalue

	integer idum,manpar

	idum=manpar(name,dvalue,0,'add')
	idum=idum*2

	end

c**************************************************************

	function infpar(type)

c info about parameters :
c       type = 'f' --> how many entries

	implicit none

	integer infpar
	character*(*) type

	character*1 name
	integer manpar
	double precision dvalue

	if(type(1:1).eq.'f') then       !how many entries
		infpar=manpar(name,dvalue,0,'fill')
	else
		write(6,*) 'type not recognized : ',type
		stop 'error stop : infpar'
	end if

	end

c**************************************************************

	subroutine lstpar(name,value,id)

c lists entry id in name and value

	implicit none

	character*(*) name
	real value
	integer id

	integer idum,manpar
	double precision dvalue

	dvalue = value
	idum=manpar(name,dvalue,id,'list')
	idum=idum*2
	value = dvalue

	end

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

	subroutine pripar(iunit)

c prints parameter values

	implicit none

	integer iunit

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
	npara=infpar('fill')

	imod=0
	do i=1,npara
	  call lstpar(name,value,i)
	  if(bflag.and.value.eq.flag) goto 1
	  imod=imod+1
	  ianf=20*(imod-1)+1
	  iend=20*imod
	  if(intpar(name).eq.1) then
		write(line(ianf:iend),2345) 
     +                  name,' =',nint(value),'  '
 2345           format(a6,a2,i10,a2)
	  else
		absval=abs(value)
		if(absval.lt.1000.and.absval.ge.0.01
     +                                  .or.absval.eq.0.) then
			write(line(ianf:iend),2347)
     +                          name,' =',value,'  '
 2347                   format(a6,a2,f10.3,a2)
		else
			write(line(ianf:iend),2346)
     +                          name,' =',value,'  '
 2346                   format(a6,a2,e10.2,a2)
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

	end

c**********************************************************

	subroutine check_parameter_values
	end
