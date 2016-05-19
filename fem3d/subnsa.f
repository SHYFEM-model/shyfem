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
c
c**********************************************
c
	subroutine nlsa(iu,ivar)
c
c read of parameter file for post processing routines
c
c iunit		unit number of file

	implicit none

	integer iu		!unit where file is open
	integer ivar		!what type of section to read

	character*80 name,line,section,extra
	character*20 what0,whatin
	character*6 sect
	logical bdebug,bread,bverbose,bskip
	integer num
	integer nrdsec,nrdlin,ichanm
	integer iv_in,iv_read
	integer iunit
	real getpar

	include 'param.h'
	include 'simul.h'

	bdebug = .true.
	bdebug = .false.
	bverbose = .true.
	bverbose = .false.

	iv_in = ivar

	if(iu.eq.0) then
c		write(6,*) 'error reading parameter file'
c		write(6,*) 'parameters initialized to default values'
		return
	end if

	if( iu > 0 ) call nrdini(iu)
	iunit = abs(iu)

c loop over sections %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	do while( nrdsec(section,num,extra) .eq. 1 )

		sect = trim(section)
		if( bverbose ) then
		  write(6,*) 'new section: ',trim(section),num,trim(extra)
		end if

		iv_read = -1
		if( extra .ne. ' ' ) then
		  call string2ivar(extra,iv_read)
		end if

		bread = iv_read .eq. iv_in
		bskip = .not. bread
		if( bverbose ) then
		  write(6,*) 'nlsa : ',bread,sect,iv_in,iv_read
		end if

		if( bread ) then
                  call setsec(section,num)                !remember section
		  if( bverbose ) then
		    write(6,'(a,a)') ' reading $',section(1:60)
		  end if
		end if

		if( bskip ) then
			call nrdskp
		  	if( bverbose ) then
		    	  write(6,'(a,a)') ' skipping $',section(1:60)
		  	end if
		else if(section.eq.'title') then
			call rdtita
		else if(section.eq.'para') then
			call nrdins(section)
		else if(section.eq.'color') then
			call colrd
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

	if( bdebug ) call prifnm(6)

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

        subroutine string2ivar_0(string,iv)

c interprets string to associate a variable number iv
c
c see below for possible string names
c
c the special name ivar# can be used to directtly give the variable number #

        implicit none

        character*(*) string
        integer iv

	integer is,isb
	integer ie3,ie4,ie5
	integer ichafs

        iv = -1

	is = ichafs(string)
	if( is .le. 0 ) is = 1
	isb = is - 1
	ie3 = isb + 3
	ie4 = isb + 4
	ie5 = isb + 5

        if( string(is:ie4) .eq. 'mass' ) then
          iv = 0
        else if( string(is:ie5) .eq. 'level' ) then
          iv = 1
        else if( string(is:ie4) .eq. 'zeta' ) then
          iv = 1
        else if( string(is:ie3) .eq. 'vel' ) then
          iv = 2
        else if( string(is:ie5) .eq. 'trans' ) then
          iv = 3
        else if( string(is:ie4) .eq. 'bath' ) then
          iv = 5
        else if( string(is:ie5) .eq. 'depth' ) then
          iv = 5
        else if( string(is:ie3) .eq. 'cur' ) then
          iv = 6
        else if( string(is:ie5) .eq. 'speed' ) then
          iv = 6
        else if( string(is:ie3) .eq. 'dir' ) then
          iv = 7
        else if( string(is:ie4) .eq. 'conc' ) then
          iv = 10
        else if( string(is:ie4) .eq. 'conz' ) then
          iv = 10
        else if( string(is:ie3) .eq. 'sal' ) then
          iv = 11
        else if( string(is:ie4) .eq. 'temp' ) then
          iv = 12
        else if( string(is:ie4) .eq. 'oxyg' ) then
          iv = 15
        else if( string(is:ie3) .eq. 'rms' ) then
          iv = 18
        else if( string(is:ie4) .eq. 'pres' ) then
          iv = 20
        else if( string(is:ie4) .eq. 'wind' ) then
          iv = 21
        else if( string(is:ie4) .eq. 'sola' ) then
          iv = 22
        else if( string(is:ie3) .eq. 'air' ) then
          iv = 23
        else if( string(is:ie4) .eq. 'humi' ) then
          iv = 24
        else if( string(is:ie4) .eq. 'clou' ) then
          iv = 25
        else if( string(is:ie4) .eq. 'rain' ) then
          iv = 26
        else if( string(is:ie4) .eq. 'evap' ) then
          iv = 27
        else if( string(is:ie3) .eq. 'lgr' ) then
          iv = 80
        else if( string(is:ie3) .eq. 'ice' ) then
          iv = 85
!        else if( string(is:ie3) .eq. 'age' ) then
!          iv = 98
        else if( string(is:ie3) .eq. 'wrt' ) then
          iv = 99
        else if( string(is:ie5) .eq. 'renew' ) then
          iv = 99
        else if( string(is:ie4) .eq. 'resi' ) then
          iv = 99
        else if( string(is:ie4) .eq. 'ivar' ) then
	  read(string(ie4+1:),'(i5)') iv
        else if( string(is:ie3) .eq. 'var' ) then
	  read(string(ie3+1:),'(i5)') iv
        else if( string(is:ie3) .eq. 'nos' ) then
          !generic - no id
        else if( string(is:ie3) .eq. 'fem' ) then
          !generic - no id
        else if( string(is:ie4) .eq. 'elem' ) then
          !generic - no id
        else if( string .eq. ' ' ) then
          write(6,*) '*** string2ivar: no string given'
        else
          write(6,*) '*** string2ivar: cannot find string description: '
          write(6,*) string
          write(6,*) is,isb,ie3,ie4,ie5
          !if( string(1:3) .eq. 'fem' ) stop 'error.....'
        end if

	!write(6,*) 'string2ivar: ',string(is:ie4),'   ',iv

        end

c******************************************************

        subroutine ivar2string_0(iv,string)

        implicit none

        integer iv
        character*(*) string

	!call shy_ivar2string(iv,string)
        string = ' '

        if( iv .eq. 0 ) then
          string = 'mass field'
        else if( iv .eq. 1 ) then
          string = 'water level'
        else if( iv .eq. 2 ) then
          string = 'velocity'
        else if( iv .eq. 5 ) then
          string = 'bathymetry'
        else if( iv .eq. 10 ) then
          string = 'generic tracer'
        else if( iv .eq. 30 ) then
          string = 'generic tracer'
        else if( iv .eq. 11 ) then
          string = 'salinity'
        else if( iv .eq. 12 ) then
          string = 'temperature'
        else if( iv .eq. 20 ) then
          string = 'atmospheric pressure'
        else if( iv .eq. 26 ) then
          string = 'rain'
        else if( iv .eq. 85 ) then
          string = 'ice cover'
!        else if( iv .eq. 98 ) then
!          string = 'age'
        else if( iv .eq. 99 ) then
          string = 'renewal time'
        else if( iv .eq. 335 ) then
          string = 'time over threshold'
        else
          string = '*** cannot find description'
          write(6,*) '*** cannot find description for variable: '
          write(6,*) iv
	  !stop 'error stop ivar2string: no description'
        end if

        end

c******************************************************

