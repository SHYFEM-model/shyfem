c
c $Id: subver.f,v 1.136 2010-03-22 15:32:47 georg Exp $
c
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
c
c*****************************************************************

	subroutine version

c DOCS	START	P_version
c
c \newcommand{\VERSION}{5.37}
c
c DOCS	END

	implicit none

        call addfnm('versio','5.37')      !version of model 

	end

c*****************************************************************

        subroutine copyright

c writes copyright and version/dimension

	implicit none

        character*(10) vers
        integer idim

	call femver(vers,idim)

        write(6,*)
        write(6,*) ' ----------------------------------------------'
        write(6,*)
        write(6,*) ' SHYFEM - Finite Element Model for coastal seas'
        write(6,*) ' Copyright (c) Georg Umgiesser 1985-2010'
        write(6,*)
        write(6,'(i3,a,a10)') idim,'D FEM model,  version ',vers
        write(6,*)
        write(6,*) ' ----------------------------------------------'
        write(6,*)

        end

c*****************************************************************
c*****************************************************************
c*****************************************************************

        subroutine femver(vers,idim)

c returns version and dimensionality of model

	implicit none

        character*(*) vers
        integer idim

        real getpar

        call getfnm('versio',vers)
        idim = nint(getpar('dimens'))

	end

c*****************************************************************

	subroutine vers2d

c sets version and dimensionality (2D)

	implicit none

        call addpar('dimens',2.)        !2D/3D 
	call version
	call copyright

	end

c*****************************************************************

	subroutine vers3d

c sets version and dimensionality (3D)

	implicit none

        call addpar('dimens',3.)        !2D/3D 
	call version
	call copyright

	end

c*****************************************************************

