!
! $Id: subver.f,v 1.136 2010-03-22 15:32:47 georg Exp $
!
! version routines and log
!
! contents :
!
! vers2d
! vers3d
! version
!
! revision log :
!
! 22.01.1998	ggu	version 4.21
! 20.03.1998	ggu	version 4.22
! 26.03.1998	ggu	version 4.23
! 04.05.1998	ggu	version 4.30
! 11.05.1998	ggu	version 4.31
! 20.05.1998	ggu	version 4.32
! 20.05.1998	ggu	version 4.33
! 18.06.1998	ggu	version 4.34
! 19.06.1998	ggu	version 4.34a
! 19.06.1998	ggu	version 4.34b ( versio is character )
! 25.06.1998	ggu	version 4.35
! 13.07.1998	ggu	version 4.36
! 22.07.1998	ggu	version 4.40  ( sub555.f restructured )
! 18.08.1998	ggu	version 4.41
! 07.09.1998	ggu	version 4.42
! 25.01.1999	ggu	version 4.43
! 27.01.1999	ggu	version 4.50
! 31.03.1999	ggu	version 4.51
! 13.04.1999	ggu	version 4.52
! 26.05.1999	ggu	version 4.53
! 22.06.1999	ggu	version 4.54
! 09.12.1999	ggu	version 4.55
! 02.03.2000	ggu	version 4.56
! 26.05.2000	ggu	version 4.60
! 15.01.2001	ggu	version 4.61
! 16.11.2001	ggu	version 4.62
! 07.12.2001	ggu	version 4.63
! 09.10.2002	ggu	version 4.64
! 11.10.2002	ggu	version 4.65
! 12.12.2002	ggu	version 4.70
! 09.01.2003	ggu	version 4.71
! 25.03.2003	ggu	version 4.72
! 20.06.2003	ggu	version 4.73
! 30.07.2003	ggu	version 4.74
! 31.07.2003	ggu	version 4.75
! 12.08.2003	ggu	version 4.75a
! 14.08.2003	ggu	version 4.75b
! 14.08.2003	ggu	version 4.75c
! 14.08.2003	ggu	new routines copyright, femver (from subnsh)
! 14.08.2003	ggu	removed routine hvers
! 14.08.2003	ggu	version 4.75d
! 14.08.2003	ggu	version 4.75e (now version has increased)
! 20.08.2003	ggu	version 4.76
! 20.08.2003	ggu	version 4.77
! 01.09.2003	ggu	version 4.77a
! 01.09.2003	ggu	version 4.77b
! 02.09.2003	ggu	version 4.77c
! 03.09.2003	ggu	version 4.77d
! 03.09.2003	ggu	version 4.78
! 04.09.2003	ggu	version 4.78a
! 04.09.2003	ggu	version 4.78b
! 12.09.2003	ggu	version 4.78c
! 31.10.2003	ggu	version 4.80
! 14.11.2003	ggu	version 4.80a
! 10.03.2004	ggu	version 4.81
! 09.08.2004	ggu	version 4.82
! 26.08.2004	ggu	version 4.83
! 03.09.2004	ggu	version 4.84
! 21.09.2004	ggu	version 4.85
! 02.12.2004	ggu	version 4.86
! 06.12.2004	ggu	version 4.86a
! 17.01.2005	ggu	version 4.87
! 26.01.2005	ggu	version 4.88
! 24.02.2005	ggu	version 4.89
! 15.03.2005	ggu	version 4.90
! 03.11.2005	ggu	version 4.91
! 07.11.2005	ggu	version 4.92
! 07.11.2005	ggu	version 4.93
! 01.02.2006	ggu	version 4.94
! 08.02.2006	ggu	version 4.94a
! 09.02.2006	ggu	version 4.94b
! 22.03.2006	ggu	version 4.95
! 09.06.2006	ggu	version 4.96
! 22.09.2006	ggu	version 4.96a
! 28.09.2006	ggu	version 4.96b
! 18.10.2006	ggu	version 4.97
! 18.10.2006	ggu	version 4.98
! 20.11.2006	ggu	version 4.98a
! 29.11.2006	ggu	version 4.98b
! 20.03.2007	ggu	version 4.99
! 08.06.2007	ggu	version 5.00
! 23.08.2007	ggu	version 5.01
! 27.09.2007	ggu	version 5.02
! 08.11.2007	ggu	version 5.03
! 18.01.2008	ggu	version 5.04
! 17.03.2008	ggu	version 5.05
! 31.03.2008	ggu	version 5.05a
! 09.04.2008	ggu	version 5.05b
! 10.04.2008	ggu	version 5.05c
! 11.04.2008	ggu	version 5.06
! 16.04.2008	ggu	version 5.10
! 17.04.2008	ggu	version 5.11
! 18.04.2008	ggu	version 5.11a
! 22.04.2008	ggu	version 5.12
! 23.04.2008	ggu	version 5.13
! 29.04.2008	ggu	version 5.14
! 29.04.2008	ggu	version 5.14a
! 16.07.2008	ggu	version 5.15
! 22.07.2008	ggu	version 5.16
! 03.09.2008	ggu	version 5.16a
! 10.10.2008	ggu	version 5.17
! 03.11.2008	ggu	version 5.17a
! 20.11.2008	ggu	version 5.18
! 09.12.2008	ggu	version 5.19
! 18.12.2008	ggu	version 5.20
! 19.12.2008	ggu	version 5.20a
! 12.01.2009	ggu	version 5.21
! 13.01.2009	ggu	version 5.22
! 26.01.2009	ggu	version 5.23
! 04.02.2009	ggu	version 5.23a
! 13.02.2009	ggu	version 5.23b
! 11.03.2009	ggu	version 5.24
! 24.03.2009	ggu	version 5.25
! 31.03.2009	ggu	version 5.26
! 31.03.2009	ggu	version 5.26a
! 03.04.2009	ggu	version 5.26b
! 06.04.2009	ggu	version 5.27
! 20.04.2009	ggu	version 5.28
! 21.05.2009	ggu	version 5.28a
! 29.05.2009	ggu	version 5.28b
! 19.06.2009	ggu	version 5.29
! 14.09.2009	ggu	version 5.30
! 14.09.2009	ggu	version 5.30a
! 09.10.2009	ggu	version 5.30b
! 18.11.2009	ggu	version 5.31
! 18.01.2010	ggu	version 5.32
! 16.02.2010	ggu	version 5.33
! 17.02.2010	ggu	version 5.33a
! 22.02.2010	ggu	version 5.34
! 26.02.2010	ggu	version 5.35
! 11.03.2010	ggu	version 5.36
! 22.03.2010	ggu	version 5.37
! 26.03.2010	ggu	version 6.1.2
! 22.04.2010	ggu	version 6.1.5
! 26.04.2010	ggu	version 6.1.7
! 03.05.2010	ggu	version 6.1.8
! 22.07.2010	ggu	version 6.1.9
! 26.07.2010	ggu	version 6.1.10
! 28.09.2010	ggu	version 6.1.11
! 29.09.2010	ggu	version 6.1.12
! 08.10.2010	ggu	version 6.1.13
! 15.12.2010	ggu	version 6.1.14
! 16.12.2010	ggu	version 6.1.15
! 27.01.2011	ggu	version 6.1.17
! 17.02.2011	ggu	version 6.1.18
! 18.02.2011	ggu	version 6.1.19
! 01.03.2011	ggu	version 6.1.20
! 23.03.2011	ggu	version 6.1.21
! 14.04.2011	ggu	version 6.1.22
! 31.05.2011	ggu	version 6.1.23
! 31.05.2011	ggu	version 6.1.24
! 07.06.2011	ggu	version 6.1.25
! 08.06.2011	ggu	version 6.1.26
! 14.07.2011	ggu	version 6.1.27
! 15.07.2011	ggu	version 6.1.29
! 26.08.2011	ggu	version 6.1.30
! 26.08.2011	ggu	version 6.1.31
! 01.09.2011	ggu	version 6.1.32
! 18.10.2011	ggu	version 6.1.33
! 24.10.2011	ggu	version 6.1.34
! 04.11.2011	ggu	version 6.1.35
! 10.11.2011	ggu	version 6.1.36
! 22.11.2011	ggu	version 6.1.37
! 09.12.2011	ggu	version 6.1.38
! 12.12.2011	ggu	version 6.1.39
! 14.12.2011	ggu	version 6.1.40
! 24.01.2012	ggu	version 6.1.41
! 25.01.2012	ggu	version 6.1.42
! 27.01.2012	ggu	version 6.1.43
! 14.02.2012	ggu	version 6.1.44	Valentine day's release	'
! 17.02.2012	ggu	version 6.1.45
! 24.02.2012	ggu	version 6.1.46
! 09.03.2012	ggu	version 6.1.47
! 16.03.2012	ggu	version 6.1.48
! 19.03.2012	ggu	version 6.1.49
! 21.03.2012	ggu	version 6.1.50
! 30.03.2012	ggu	version 6.1.51
! 13.04.2012	ggu	version 6.1.52	Aniversary day's release '
! 01.06.2012	ggu	version 6.1.53
! 21.06.2012	ggu	version 6.1.54
! 26.06.2012	ggu	version 6.1.55
! 29.08.2012	ggu	version 6.1.56
! 12.09.2012	ggu	version 6.1.57
! 08.10.2012	ggu	version 6.1.58
! 25.10.2012	ggu	version 6.1.59
! 05.11.2012	ggu	version 6.1.60
! 19.11.2012	ggu	version 6.1.61
! 17.12.2012	ggu	version 6.1.61a
! 25.01.2013	ggu	version 6.1.62
! 03.05.2013	ggu	version 6.1.63
! 10.05.2013	ggu	version 6.1.64
! 13.06.2013	ggu	version 6.1.65
! 19.06.2013	ggu	version 6.1.66
! 12.09.2013	ggu	version 6.1.67
! 25.10.2013	ggu	version 6.1.68
! 11.11.2013	ggu	version 6.1.69	San Martino's release '
! 12.11.2013	ggu	version 6.1.69a
! 05.12.2013	ggu	version 6.1.70
! 28.01.2014	ggu	version 6.1.71
! 07.03.2014	ggu	version 6.1.72
! 27.03.2014	ggu	version 6.1.73
! 05.05.2014	ggu	version 6.1.74
! 15.05.2014	ggu	version 6.1.75  Frau Flierl's birthday release '
! 30.05.2014	ggu	version 6.1.76  Clara's birthday release '
! 18.06.2014	ggu	version 6.1.77  Clara's maturita' release
! 27.06.2014	ggu	version 6.1.78  Maria Huber release
! 07.07.2014	ggu	version 6.1.79  last of 6 series
! 07.07.2014	ggu	version 7.0.0   first of 7 series
! 18.07.2014	ggu	version 7.0.1   last day before holidays
! 13.10.2014	ggu	version 7.0.2
! 21.10.2014	ggu	version 7.0.3
! 30.10.2014	ggu	version 7.0.4
! 05.11.2014	ggu	version 7.0.5
! 07.11.2014	ggu	version 7.0.6
! 26.11.2014	ggu	version 7.0.7	Peter Epplers Geburtstags release
! 05.12.2014	ggu	version 7.0.8	pre Santa Claus release
! 12.12.2014	ggu	version 7.0.9	pre St Lucia release
! 19.12.2014	ggu	version 7.0.10	pre Christmas release
! 23.12.2014	ggu	version 7.0.11	pre Christmas eve release
! 09.01.2015	ggu	version 7.0.12
! 12.01.2015	ggu	version 7.1.0	Wolfgangs release
! 15.01.2015	ggu	version 7.1.1
! 19.01.2015	ggu	version 7.1.2
! 19.01.2015	ggu	version 7.1.3	huge commit (include for common)
! 23.01.2015	ggu	version 7.1.4
! 25.02.2015	ggu	version 7.1.5
! 26.02.2015	ggu	version 7.1.6
! 01.04.2015	ggu	version 7.1.7	first of aprile version
! 23.04.2015	ggu	version 7.1.8	St. George's release '
! 30.04.2015	ggu	version 7.1.9	Harald's birthday '
! 05.05.2015	ggu	version 7.1.10
! 21.05.2015	ggu	version 7.1.11
! 05.06.2015	ggu	version 7.1.12
! 10.07.2015	ggu	version 7.1.50	big release...
! 13.07.2015	ggu	version 7.1.51
! 17.07.2015	ggu	version 7.1.52
! 17.07.2015	ggu	version 7.1.53
! 17.07.2015	ggu	version 7.1.80
! 20.07.2015	ggu	version 7.1.81
! 24.07.2015	ggu	version 7.1.82
! 30.07.2015	ggu	version 7.1.83
! 31.07.2015	ggu	version 7.1.84
! 31.07.2015	ggu	version 7.2.1	holiday release
! 14.09.2015	ggu	version 7.2.2
! 18.09.2015	ggu	version 7.2.3
! 23.09.2015	ggu	version 7.2.4
! 29.09.2015	ggu	version 7.2.5
! 30.09.2015	ggu	version 7.2.6
! 02.10.2015	ggu	version 7.3.1	first of development branch
! 10.10.2015	ggu	version 7.3.2
! 12.10.2015	ggu	version 7.3.3
! 12.10.2015	ggu	version 7.3.4
! 12.10.2015	ggu	version 7.3.4a
! 13.10.2015	ggu	version 7.3.5
! 19.10.2015	ggu	version 7.3.6
! 22.10.2015	ggu	version 7.3.7	Tante Lores birthday release
! 22.10.2015	ggu	version 7.3.8	Tante Lores birthday release again
! 23.10.2015	ggu	version 7.3.9	after Eric
! 26.10.2015	ggu	version 7.3.10
! 05.11.2015	ggu	version 7.3.11
! 05.11.2015	ggu	version 7.3.12
! 09.11.2015	ggu	version 7.3.13	the wall release
! 16.11.2015	ggu	version 7.3.14
! 20.11.2015	ggu	version 7.3.15	pre Madonna della Salute release
! 16.12.2015	ggu	version 7.3.16
! 18.12.2015	ggu	version 7.3.17	Christmas 2015 release
! 08.01.2016	ggu	version 7.3.18
! 08.01.2016	ggu	version 7.4.0	major stable version
! 08.01.2016	ggu	version 7.5.0	new develop version opened
!
!*****************************************************************
!-----------------------------------------------------------------
        module version
!-----------------------------------------------------------------
 
! DOCS	START	P_version
!
! \newcommand{\VERSION}{7.5.0}
! \newcommand{\version}{7\_5\_0}
!
! DOCS	END

        character*10 vers
        parameter (vers='7.5.0')

        character*30 string
        parameter (string='SHYFEM VERSION = '//vers)

        character*10, save :: shyver=vers
        character*30, save :: shystr=string


!-----------------------------------------------------------------
        contains
!-----------------------------------------------------------------

        subroutine get_shyfem_version(version)

! returns version of model

	implicit none

	character*(*) version

	version = shyver

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

        subroutine shyfem_copyright(routine)

! writes copyright and version/dimension

	implicit none

        character*(*) routine

        character*10 vers

	call get_shyfem_version(vers)

        write(6,*)
        write(6,*) ' ----------------------------------------------'
        write(6,*)
        write(6,*) ' SHYFEM - Finite Element Model for coastal seas'
        write(6,*) ' Copyright (c) The Shyfem Team 1985-2015'
        write(6,*)
        write(6,*) ' version: ',vers
        write(6,*) ' routine: ',routine
        write(6,*)
        write(6,*) ' ----------------------------------------------'
        write(6,*)

        end

!*****************************************************************

	subroutine shyfem_copyright_test

	implicit none

	call shyfem_copyright('test_routine')

	end

!*****************************************************************

!	program shyfem_copyright_main
!	implicit none
!	call shyfem_copyright_test
!	end

!*****************************************************************
!-----------------------------------------------------------------
      end module version
!-----------------------------------------------------------------
