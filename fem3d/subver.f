
!--------------------------------------------------------------------------
!
!    Copyright (C) 1998-2020  Georg Umgiesser
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

c version routines and log
c
c contents :
c
c vers2d
c vers3d
c version
c
c revision log :
c
c 22.01.1998	ggu	version 4.21
c 20.03.1998	ggu	version 4.22
c 26.03.1998	ggu	version 4.23
c 04.05.1998	ggu	version 4.30
c 11.05.1998	ggu	version 4.31
c 20.05.1998	ggu	version 4.32
c 20.05.1998	ggu	version 4.33
c 18.06.1998	ggu	version 4.34
c 19.06.1998	ggu	version 4.34a
c 19.06.1998	ggu	version 4.34b ( versio is character )
c 25.06.1998	ggu	version 4.35
c 13.07.1998	ggu	version 4.36
c 22.07.1998	ggu	version 4.40  ( sub555.f restructured )
c 18.08.1998	ggu	version 4.41
c 07.09.1998	ggu	version 4.42
c 25.01.1999	ggu	version 4.43
c 27.01.1999	ggu	version 4.50
c 31.03.1999	ggu	version 4.51
c 13.04.1999	ggu	version 4.52
c 26.05.1999	ggu	version 4.53
c 22.06.1999	ggu	version 4.54
c 09.12.1999	ggu	version 4.55
c 02.03.2000	ggu	version 4.56
c 26.05.2000	ggu	version 4.60
c 15.01.2001	ggu	version 4.61
c 16.11.2001	ggu	version 4.62
c 07.12.2001	ggu	version 4.63
c 09.10.2002	ggu	version 4.64
c 11.10.2002	ggu	version 4.65
c 12.12.2002	ggu	version 4.70
c 09.01.2003	ggu	version 4.71
c 25.03.2003	ggu	version 4.72
c 20.06.2003	ggu	version 4.73
c 30.07.2003	ggu	version 4.74
c 31.07.2003	ggu	version 4.75
c 12.08.2003	ggu	version 4.75a
c 14.08.2003	ggu	version 4.75b
c 14.08.2003	ggu	version 4.75c
c 14.08.2003	ggu	new routines copyright, femver (from subnsh)
c 14.08.2003	ggu	removed routine hvers
c 14.08.2003	ggu	version 4.75d
c 14.08.2003	ggu	version 4.75e (now version has increased)
c 20.08.2003	ggu	version 4.76
c 20.08.2003	ggu	version 4.77
c 01.09.2003	ggu	version 4.77a
c 01.09.2003	ggu	version 4.77b
c 02.09.2003	ggu	version 4.77c
c 03.09.2003	ggu	version 4.77d
c 03.09.2003	ggu	version 4.78
c 04.09.2003	ggu	version 4.78a
c 04.09.2003	ggu	version 4.78b
c 12.09.2003	ggu	version 4.78c
c 31.10.2003	ggu	version 4.80
c 14.11.2003	ggu	version 4.80a
c 10.03.2004	ggu	version 4.81
c 09.08.2004	ggu	version 4.82
c 26.08.2004	ggu	version 4.83
c 03.09.2004	ggu	version 4.84
c 21.09.2004	ggu	version 4.85
c 02.12.2004	ggu	version 4.86
c 06.12.2004	ggu	version 4.86a
c 17.01.2005	ggu	version 4.87
c 26.01.2005	ggu	version 4.88
c 24.02.2005	ggu	version 4.89
c 15.03.2005	ggu	version 4.90
c 03.11.2005	ggu	version 4.91
c 07.11.2005	ggu	version 4.92
c 07.11.2005	ggu	version 4.93
c 01.02.2006	ggu	version 4.94
c 08.02.2006	ggu	version 4.94a
c 09.02.2006	ggu	version 4.94b
c 22.03.2006	ggu	version 4.95
c 09.06.2006	ggu	version 4.96
c 22.09.2006	ggu	version 4.96a
c 28.09.2006	ggu	version 4.96b
c 18.10.2006	ggu	version 4.97
c 18.10.2006	ggu	version 4.98
c 20.11.2006	ggu	version 4.98a
c 29.11.2006	ggu	version 4.98b
c 20.03.2007	ggu	version 4.99
c 08.06.2007	ggu	version 5.00
c 23.08.2007	ggu	version 5.01
c 27.09.2007	ggu	version 5.02
c 08.11.2007	ggu	version 5.03
c 18.01.2008	ggu	version 5.04
c 17.03.2008	ggu	version 5.05
c 31.03.2008	ggu	version 5.05a
c 09.04.2008	ggu	version 5.05b
c 10.04.2008	ggu	version 5.05c
c 11.04.2008	ggu	version 5.06
c 16.04.2008	ggu	version 5.10
c 17.04.2008	ggu	version 5.11
c 18.04.2008	ggu	version 5.11a
c 22.04.2008	ggu	version 5.12
c 23.04.2008	ggu	version 5.13
c 29.04.2008	ggu	version 5.14
c 29.04.2008	ggu	version 5.14a
c 16.07.2008	ggu	version 5.15
c 22.07.2008	ggu	version 5.16
c 03.09.2008	ggu	version 5.16a
c 10.10.2008	ggu	version 5.17
c 03.11.2008	ggu	version 5.17a
c 20.11.2008	ggu	version 5.18
c 09.12.2008	ggu	version 5.19
c 18.12.2008	ggu	version 5.20
c 19.12.2008	ggu	version 5.20a
c 12.01.2009	ggu	version 5.21
c 13.01.2009	ggu	version 5.22
c 26.01.2009	ggu	version 5.23
c 04.02.2009	ggu	version 5.23a
c 13.02.2009	ggu	version 5.23b
c 11.03.2009	ggu	version 5.24
c 24.03.2009	ggu	version 5.25
c 31.03.2009	ggu	version 5.26
c 31.03.2009	ggu	version 5.26a
c 03.04.2009	ggu	version 5.26b
c 06.04.2009	ggu	version 5.27
c 20.04.2009	ggu	version 5.28
c 21.05.2009	ggu	version 5.28a
c 29.05.2009	ggu	version 5.28b
c 19.06.2009	ggu	version 5.29
c 14.09.2009	ggu	version 5.30
c 14.09.2009	ggu	version 5.30a
c 09.10.2009	ggu	version 5.30b
c 18.11.2009	ggu	version 5.31
c 18.01.2010	ggu	version 5.32
c 16.02.2010	ggu	version 5.33
c 17.02.2010	ggu	version 5.33a
c 22.02.2010	ggu	version 5.34
c 26.02.2010	ggu	version 5.35
c 11.03.2010	ggu	version 5.36
c 22.03.2010	ggu	version 5.37
c 26.03.2010	ggu	version 6.1.2
c 22.04.2010	ggu	version 6.1.5
c 26.04.2010	ggu	version 6.1.7
c 03.05.2010	ggu	version 6.1.8
c 22.07.2010	ggu	version 6.1.9
c 26.07.2010	ggu	version 6.1.10
c 28.09.2010	ggu	version 6.1.11
c 29.09.2010	ggu	version 6.1.12
c 08.10.2010	ggu	version 6.1.13
c 15.12.2010	ggu	version 6.1.14
c 16.12.2010	ggu	version 6.1.15
c 27.01.2011	ggu	version 6.1.17
c 17.02.2011	ggu	version 6.1.18
c 18.02.2011	ggu	version 6.1.19
c 01.03.2011	ggu	version 6.1.20
c 23.03.2011	ggu	version 6.1.21
c 14.04.2011	ggu	version 6.1.22
c 31.05.2011	ggu	version 6.1.23
c 31.05.2011	ggu	version 6.1.24
c 07.06.2011	ggu	version 6.1.25
c 08.06.2011	ggu	version 6.1.26
c 14.07.2011	ggu	version 6.1.27
c 15.07.2011	ggu	version 6.1.29
c 26.08.2011	ggu	version 6.1.30
c 26.08.2011	ggu	version 6.1.31
c 01.09.2011	ggu	version 6.1.32
c 18.10.2011	ggu	version 6.1.33
c 24.10.2011	ggu	version 6.1.34
c 04.11.2011	ggu	version 6.1.35
c 10.11.2011	ggu	version 6.1.36
c 22.11.2011	ggu	version 6.1.37
c 09.12.2011	ggu	version 6.1.38
c 12.12.2011	ggu	version 6.1.39
c 14.12.2011	ggu	version 6.1.40
c 24.01.2012	ggu	version 6.1.41
c 25.01.2012	ggu	version 6.1.42
c 27.01.2012	ggu	version 6.1.43
c 14.02.2012	ggu	version 6.1.44	Valentine day's release	'
c 17.02.2012	ggu	version 6.1.45
c 24.02.2012	ggu	version 6.1.46
c 09.03.2012	ggu	version 6.1.47
c 16.03.2012	ggu	version 6.1.48
c 19.03.2012	ggu	version 6.1.49
c 21.03.2012	ggu	version 6.1.50
c 30.03.2012	ggu	version 6.1.51
c 13.04.2012	ggu	version 6.1.52	Aniversary day's release '
c 01.06.2012	ggu	version 6.1.53
c 21.06.2012	ggu	version 6.1.54
c 26.06.2012	ggu	version 6.1.55
c 29.08.2012	ggu	version 6.1.56
c 12.09.2012	ggu	version 6.1.57
c 08.10.2012	ggu	version 6.1.58
c 25.10.2012	ggu	version 6.1.59
c 05.11.2012	ggu	version 6.1.60
c 19.11.2012	ggu	version 6.1.61
c 17.12.2012	ggu	version 6.1.61a
c 25.01.2013	ggu	version 6.1.62
c 03.05.2013	ggu	version 6.1.63
c 10.05.2013	ggu	version 6.1.64
c 13.06.2013	ggu	version 6.1.65
c 19.06.2013	ggu	version 6.1.66
c 12.09.2013	ggu	version 6.1.67
c 25.10.2013	ggu	version 6.1.68
c 11.11.2013	ggu	version 6.1.69	San Martino's release '
c 12.11.2013	ggu	version 6.1.69a
c 05.12.2013	ggu	version 6.1.70
c 28.01.2014	ggu	version 6.1.71
c 07.03.2014	ggu	version 6.1.72
c 27.03.2014	ggu	version 6.1.73
c 05.05.2014	ggu	version 6.1.74
c 15.05.2014	ggu	version 6.1.75  Frau Flierl's birthday release '
c 30.05.2014	ggu	version 6.1.76  Clara's birthday release '
c 18.06.2014	ggu	version 6.1.77  Clara's maturita' release
c 27.06.2014	ggu	version 6.1.78  Maria Huber release
c 07.07.2014	ggu	version 6.1.79  last of 6 series
c 07.07.2014	ggu	version 7.0.0   first of 7 series
c 18.07.2014	ggu	version 7.0.1   last day before holidays
c 13.10.2014	ggu	version 7.0.2
c 21.10.2014	ggu	version 7.0.3
c 30.10.2014	ggu	version 7.0.4
c 05.11.2014	ggu	version 7.0.5
c 07.11.2014	ggu	version 7.0.6
c 26.11.2014	ggu	version 7.0.7	Peter Epplers Geburtstags release
c 05.12.2014	ggu	version 7.0.8	pre Santa Claus release
c 12.12.2014	ggu	version 7.0.9	pre St Lucia release
c 19.12.2014	ggu	version 7.0.10	pre Christmas release
c 23.12.2014	ggu	version 7.0.11	pre Christmas eve release
c 09.01.2015	ggu	version 7.0.12
c 12.01.2015	ggu	version 7.1.0	Wolfgangs release
c 15.01.2015	ggu	version 7.1.1
c 19.01.2015	ggu	version 7.1.2
c 19.01.2015	ggu	version 7.1.3	huge commit (include for common)
c 23.01.2015	ggu	version 7.1.4
c 25.02.2015	ggu	version 7.1.5
c 26.02.2015	ggu	version 7.1.6
c 01.04.2015	ggu	version 7.1.7	first of aprile version
c 23.04.2015	ggu	version 7.1.8	St. George's release '
c 30.04.2015	ggu	version 7.1.9	Harald's birthday '
c 05.05.2015	ggu	version 7.1.10
c 21.05.2015	ggu	version 7.1.11
c 05.06.2015	ggu	version 7.1.12
c 10.07.2015	ggu	version 7.1.50	big release...
c 13.07.2015	ggu	version 7.1.51
c 17.07.2015	ggu	version 7.1.52
c 17.07.2015	ggu	version 7.1.53
c 17.07.2015	ggu	version 7.1.80
c 20.07.2015	ggu	version 7.1.81
c 24.07.2015	ggu	version 7.1.82
c 30.07.2015	ggu	version 7.1.83
c 31.07.2015	ggu	version 7.1.84
c 31.07.2015	ggu	version 7.2.1	holiday release
c 14.09.2015	ggu	version 7.2.2
c 18.09.2015	ggu	version 7.2.3
c 23.09.2015	ggu	version 7.2.4
c 29.09.2015	ggu	version 7.2.5
c 30.09.2015	ggu	version 7.2.6
c 02.10.2015	ggu	version 7.3.1	first of development branch
c 10.10.2015	ggu	version 7.3.2
c 12.10.2015	ggu	version 7.3.3
c 12.10.2015	ggu	version 7.3.4
c 12.10.2015	ggu	version 7.3.4a
c 13.10.2015	ggu	version 7.3.5
c 19.10.2015	ggu	version 7.3.6
c 22.10.2015	ggu	version 7.3.7	Tante Lores birthday release
c 22.10.2015	ggu	version 7.3.8	Tante Lores birthday release again
c 23.10.2015	ggu	version 7.3.9	after Eric
c 26.10.2015	ggu	version 7.3.10
c 05.11.2015	ggu	version 7.3.11
c 05.11.2015	ggu	version 7.3.12
c 09.11.2015	ggu	version 7.3.13	the wall release
c 16.11.2015	ggu	version 7.3.14
c 20.11.2015	ggu	version 7.3.15	pre Madonna della Salute release
c 16.12.2015	ggu	version 7.3.16
c 18.12.2015	ggu	version 7.3.17	Christmas 2015 release
c 08.01.2016	ggu	version 7.3.18
c 08.01.2016	ggu	version 7.4.0	major stable version
c 08.01.2016	ggu	version 7.5.0	new develop version opened
c 22.01.2016	ggu	version 7.5.1
c 19.02.2016	ggu	version 7.5.2
c 19.02.2016	ggu	version 7.5.3
c 22.02.2016	ggu	version 7.5.4
c 11.03.2016	ggu	version 7.5.5
c 22.03.2016	ggu	version 7.5.6	Alessandra's birthday release '
c 01.04.2016	ggu	version 7.5.7	no April joke
c 15.04.2016	ggu	version 7.5.8	after 20th anniversary
c 28.04.2016	ggu	version 7.5.9
c 25.05.2016	ggu	version 7.5.10	for Leslie
c 30.05.2016	ggu	version 7.5.11	Clara's birthday release '
c 06.06.2016	ggu	version 7.5.12	Monika's 70th birthday release '
c 10.06.2016	ggu	version 7.5.13
c 14.06.2016	ggu	version 7.5.14
c 17.06.2016	ggu	version 7.5.15	Italy-Sweden release
c 27.06.2016	ggu	version 7.5.16	Mutti's birthday release '
c 09.09.2016	ggu	version 7.5.17	after holiday release... sniff
c 30.09.2016	ggu	version 7.5.18
c 05.10.2016	ggu	version 7.5.19
c 11.10.2016	ggu	version 7.5.20
c 12.01.2017	ggu	version 7.5.21
c 20.01.2017	ggu	version 7.5.22	God bless America's release '
c 13.02.2017	ggu	version 7.5.23
c 31.03.2017	ggu	version 7.5.24
c 13.04.2017	ggu	version 7.5.25	pre Good Friday release
c 09.05.2017	ggu	version 7.5.26
c 16.05.2017	ggu	version 7.5.27	Beppe's 50th birthday release '
c 25.05.2017	ggu	version 7.5.28
c 13.06.2017	ggu	version 7.5.29	San Antonio's name day release '
c 11.07.2017	ggu	version 7.5.30	pre holiday release
c 02.09.2017	ggu	version 7.5.31	Memel release
c 26.09.2017	ggu	version 7.5.32	Murcia release
c 09.10.2017	ggu	version 7.5.33
c 04.11.2017	ggu	version 7.5.34	Forze armate release
c 04.11.2017	ggu	version 7.5.35	... and some stupid forgotten things
c 14.11.2017	ggu	version 7.5.36
c 17.11.2017	ggu	version 7.5.37
c 17.11.2017	ggu	version 7.5.38	brown paper bag bug...
c 05.12.2017	ggu	version 7.5.39
c 07.12.2017	ggu	version 7.5.40
c 24.01.2018	ggu	version 7.5.41
c 22.02.2018	ggu	version 7.5.42	post Lithuania release
c 03.04.2018	ggu	version 7.5.43	post Easter 2018 release
c 03.04.2018	ggu	version 7.5.44	small bug fix
c 19.04.2018	ggu	version 7.5.45
c 26.04.2018	ggu	version 7.5.46
c 11.05.2018	ggu	version 7.5.47
c 06.07.2018	ggu	version 7.5.48
c 13.07.2018	ggu	version 7.4.1	stable release
c 31.08.2018	ggu	version 7.5.49
c 16.10.2018	ggu	version 7.5.50
c 25.10.2018	ggu	version 7.5.51
c 18.12.2018	ggu	version 7.5.52
c 21.12.2018	ggu	version 7.5.53	Christmas 2018 edition
c 27.12.2018	ggu	version 7.5.54
c 18.01.2019	ggu	version 7.5.55	penta testing
c 14.02.2019	ggu	version 7.5.55	San Valentine's release '
c 16.02.2019	ggu	version 7.5.60	copyrighted release
c 13.03.2019	ggu	version 7.5.61
c 21.05.2019	ggu	version 7.5.62
c 02.07.2019	ggu	version 7.5.63
c 19.07.2019	ggu	version 7.5.64	
c 31.10.2019	ggu	version 7.5.65	Halloween 2019 edition	
c 25.11.2019	ggu	version 7.5.66
c 20.12.2019	ggu	version 7.5.67	Christmas 2019 edition
c 31.01.2020	ggu	version 7.5.68	Brexit edition
c 06.03.2020	ggu	version 7.5.69	Vincenzo edition
c 19.05.2020	ggu	version 7.5.70	Covid edition
c 27.05.2021	ggu	version 7.5.71	3rd wave Covid edition
c 15.03.2022	ggu	version 7.5.72	Ides of March edition
c
c*****************************************************************

!=================================================================
	module shyfem_version
!=================================================================

c DOCS	START	P_version
c
c \newcommand{\VERSION}{7.5.72}
c \newcommand{\version}{7\_5\_72}
c \newcommand{\COMMIT}{2022-03-16}
c
c DOCS	END

        implicit none

	logical, save		:: bshort = .false.

        character*10, parameter :: version = '7.5.72'
        character*10, parameter :: commit  = '2022-03-16'
        character*17, parameter :: text    = 'SHYFEM VERSION = '

        character*40, parameter :: string = text//version//'  '//commit

	character*50, parameter :: acronym =
     +	    	'System of HydrodYnamic Finite Element Modules'
	character*50, parameter :: copyright =
     +		'Copyright (C) The Shyfem Team 1985-2022'

!=================================================================
	end module shyfem_version
!=================================================================

c*****************************************************************
c*****************************************************************
c*****************************************************************

        subroutine get_shyfem_version_and_commit(vers)

c returns version and commit of model

	use shyfem_version

	implicit none

	character*(*) vers

	vers = 'version '//trim(version)//' - commit '//trim(commit)
	vers = trim(version)//' - '//trim(commit)

	end

c*****************************************************************

        subroutine get_shyfem_version(vers)

c returns version of model

	use shyfem_version

	implicit none

	character*(*) vers

	vers = version

	end

c*****************************************************************

        subroutine get_shyfem_commit(comm)

c returns version/commit of model

	use shyfem_version

	implicit none

	character*(*) comm

	comm = commit

	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

        subroutine shyfem_copyright(routine)

c writes copyright and version/dimension

	use shyfem_version

	implicit none

        character*(*) routine

        character*10 vers
        character*10 comm

	call get_shyfem_version(vers)
	call get_shyfem_commit(comm)

	if( bshort ) then

        write(6,*) 'SHYFEM - '//acronym
        write(6,*) copyright
        write(6,*) routine

	else

        write(6,*)
        write(6,*) '----------------------------------------------'
        write(6,*)
        write(6,*) 'SHYFEM'
        write(6,*) acronym
        write(6,*) copyright
        write(6,*)
        write(6,*) 'version: ',vers
        write(6,*) 'commit : ',comm
        write(6,*) 'routine: ',routine
        write(6,*)
        write(6,*) '----------------------------------------------'
        write(6,*)

	end if

        end

c*****************************************************************

	subroutine shyfem_set_short_copyright(bset)

	use shyfem_version

	implicit none

	logical bset

	bshort = bset

	end

c*****************************************************************

	subroutine shyfem_copyright_test

	implicit none

	call shyfem_copyright('test_routine')

	end

c*****************************************************************

c	program shyfem_copyright_main
c	implicit none
c	call shyfem_copyright_test
c	end

c*****************************************************************
