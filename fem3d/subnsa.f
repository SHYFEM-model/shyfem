c
c $Id: subnsa.f,v 1.21 2009-09-14 08:20:58 georg Exp $
c
c ap utility routines
c
c contents :
c
c subroutine nlsa(iunit)		read of parameter file for pp routines
c subroutine rdtita			reads title section of apn files
c
c revision log :
c
c revised on 25.11.88 by ggu (no parameter, merged sp158k and randk)
c revised on 30.11.88 by ggu (no array ig,iamat any more)
c revised on 08.10.90 by ggu (kantv is ordened with direction)
c 27.03.1998	ggu	eliminated /bnd/, /irv/
c 12.02.1999	ggu	reading title from own subroutine (with sim and bas)
c 06.12.2004	ggu	new section legvar
c 11.03.2005	ggu	write section title to stdout
c 11.09.2009	ggu	new section $sect
c 27.02.2013	ggu     pass what parameter into nlsa, handle extra info
c 13.06.2013	ggu     read also varnam to decide what to plot and read
c 22.08.2013	ggu     new string2ivar() and similar changes
c 05.09.2013	ggu     nlsa now wants integer, better handling of what to read
c 06.06.2014	ggu	deleted sp158k() and sp158kk()
c 07.06.2016	ggu	use both varid and varname to decide on section reading
c 02.02.2017	ggu	nlsa simplified
c
c**********************************************
c
	subroutine nlsa(iu,ivar,bverb)
c
c read of parameter file for post processing routines
c
c iunit		unit number of file

	implicit none

	integer iu		!unit where file is open
	integer ivar		!what type of section to read
	logical bverb

	character*80 name,line,section,extra
	character*20 what0,whatin
	logical bdebug,bverbose
	logical bread_str,bread_iv,bread
	integer num
	integer nrdsec,nrdlin,ichanm
	integer iv_in,iv_read,isub
	integer iunit
	character*80 str_read,str_in
	real getpar

	logical compare_svars

	include 'simul.h'

	bdebug = .true.
	bdebug = .false.
	bverbose = bverb

	iv_in = ivar
	call ivar2string(iv_in,str_in,isub)

	if(iu.eq.0) return

	if( iu > 0 ) call nrdini(iu)
	iunit = abs(iu)

c loop over sections %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	do while( nrdsec(section,num,extra) .eq. 1 )

		iv_read = -1
		str_read = extra
		if( extra .ne. ' ' ) then
		  call string2ivar(str_read,iv_read)
		end if

		bread_iv = iv_read .eq. iv_in
		bread_str = compare_svars(str_read,str_in)
		bread = bread_iv .or. bread_str
		if( bread ) call setsec(section,num)

		if( bverbose .or. bdebug ) then
		  if( bdebug ) then
		    write(6,*) 'section: ',trim(section),' ',trim(extra)
		    write(6,*) 'nlsa : ',bread,iv_in,iv_read
		  end if
		  if( bread ) then
		    write(6,*) 'reading ',trim(section),' ',trim(extra)
		  end if
		end if

		if( .not. bread ) then
			call nrdskp
		else if(section.eq.'title') then
			call rdtita
		else if(section.eq.'para') then
			call nrdins(section)
		else if(section.eq.'color') then
			call nrdins(section)
			!call colrd
		else if(section.eq.'arrow') then
			call nrdins(section)
		else if(section.eq.'legvar') then
			call nrdins(section)
		else if(section.eq.'legend') then
			call legrd
		else if(section.eq.'name') then
			call nrdins(section)
		else if(section.eq.'sect') then
			call nrdins(section)
		else
			goto 97
		end if
	end do

c end of read %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	!if( bdebug ) call prifnm(6)

	return
   97	continue
	write(6,*) 'Not recognized section key word : ',section
	stop 'error stop : nlsa'
	end

c******************************************************************

        subroutine rdtita

c reads title section of apn files

        implicit none

	include 'param.h'
	include 'simul.h'

        character*80 line

        integer nrdlin
	logical bdebug

	bdebug = .true.
	bdebug = .false.

c first line -> title

	if( nrdlin(line) .eq. 0 ) goto 65
	descrp=line
	call putfnm('title',line)
	if( bdebug ) write(6,*) 'rdtita: ',line(1:65)

c maybe more lines ?

	if( nrdlin(line) .eq. 0 ) return	!just one line -> return

c ok, this is simulation

        call putfnm('runnam',line)
	if( bdebug ) write(6,*) 'rdtita: ',line(1:65)

c now basin

        if( nrdlin(line) .eq. 0 ) goto 65
        call putfnm('basnam',line)
	if( bdebug ) write(6,*) 'rdtita: ',line(1:65)

c no more lines allowed

        if( nrdlin(line) .gt. 0 ) goto 65

        return
   65   continue
        write(6,*) 'error in section $title'
        stop 'error stop : rdtitl'
        end

c************************************************************************
c************************************************************************
c************************************************************************

