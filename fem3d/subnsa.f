
!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2018  Georg Umgiesser
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

c ap utility routines
c
c contents :
c
c subroutine nlsa(iunit)		read of parameter file for pp routines
c subroutine rdtita			reads title section of apn files
c
c revision log :
c
c 25.11.1988	ggu	(no parameter, merged sp158k and randk)
c 30.11.1988	ggu	(no array ig,iamat any more)
c 08.10.1990	ggu	(kantv is ordened with direction)
c 27.03.1998	ggu	eliminated /bnd/, /irv/
c 12.02.1999	ggu	reading title from own subroutine (with sim and bas)
c 06.12.2004	ggu	new section legvar
c 11.03.2005	ggu	write section title to stdout
c 11.09.2009	ggu	new section $sect
c 23.03.2010	ggu	changed v6.1.1
c 09.04.2010	ggu	changed v6.1.3
c 30.03.2012	ggu	changed VERS_6_1_51
c 27.02.2013	ggu	pass what parameter into nlsa, handle extra info
c 03.05.2013	ggu	changed VERS_6_1_63
c 13.06.2013	ggu	read also varnam to decide what to plot and read
c 22.08.2013	ggu	new string2ivar() and similar changes
c 05.09.2013	ggu	nlsa now wants integer, better handling of what to read
c 12.09.2013	ggu	changed VERS_6_1_67
c 25.10.2013	ggu	changed VERS_6_1_68
c 12.11.2013	ggu	changed VERS_6_1_69
c 28.01.2014	ggu	changed VERS_6_1_71
c 07.03.2014	ggu	changed VERS_6_1_72
c 05.05.2014	ggu	changed VERS_6_1_74
c 06.06.2014	ggu	deleted sp158k() and sp158kk()
c 18.06.2014	ggu	changed VERS_6_1_77
c 18.07.2014	ggu	changed VERS_7_0_1
c 21.10.2014	ggu	changed VERS_7_0_3
c 19.01.2015	ggu	changed VERS_7_1_3
c 26.02.2015	ggu	changed VERS_7_1_5
c 05.05.2015	ggu	changed VERS_7_1_10
c 05.06.2015	ggu	changed VERS_7_1_12
c 14.09.2015	ggu	changed VERS_7_2_2
c 29.09.2015	ggu	changed VERS_7_2_5
c 16.12.2015	ggu	changed VERS_7_3_16
c 25.05.2016	ggu	changed VERS_7_5_10
c 30.05.2016	ggu	changed VERS_7_5_11
c 07.06.2016	ggu	use both varid and varname to decide on section reading
c 02.02.2017	ggu	nlsa simplified
c 13.02.2017	ggu	changed VERS_7_5_23
c 09.05.2017	ggu	changed VERS_7_5_26
c 11.07.2017	ggu	changed VERS_7_5_30
c 07.12.2017	ggu	changed VERS_7_5_40
c 19.04.2018	ggu	changed VERS_7_5_45
c 16.02.2019	ggu	changed VERS_7_5_60
c 13.03.2019	ggu	changed VERS_7_5_61
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

