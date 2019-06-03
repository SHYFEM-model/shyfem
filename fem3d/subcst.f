
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

c routines to set up constants
c
c contents :
c
c subroutine cstinit            set up constants and parameters
c subroutine cstcheck           checks parameters read
c subroutine cstsetup           sets up modules
c subroutine cstfile(strfile)   reads files (str and bas)
c
c subroutine cktime		check time parameters
c subroutine ckdate		check date parameter
c subroutine ckcori		set coriolis parameter
c
c revision log :
c
c 05.08.1992	ggu	$$ibtyp3 - implementation of ibtyp=3
c 31.08.1992	ggu	$$impli - implicit time step
c 23.11.1992	ggu	$$ibtyp11 - implementation of ibtyp=11,51
c 27.10.1993	ggu	$$roger - implementation of ibtyp=70 (nsea)
c 07.04.1995	ggu	!$$conzfl - conz compared to iflag (bug)
c 07.04.1995	ggu	!$$baroc - impl. of baroclinic salt/temp (21/22)
c 02.06.1997	ggu	$$EXTINW - extinw changed to ipint
c 13.06.1997	ggu	!$$kranf - check if kranf <= krend
c 29.06.1997	ggu	no cstdim in file
c 29.04.1998	ggu	module for semi-implicit time-step in own routine
c 12.08.1998	ggu	new parameter dlat -> specify latitude for coriolis
c 03.09.1998	ggu	call bocche to adjust depth at Venice inlets
c 06.11.1998	ggu	call to huniqu to set up hkv and hev
c 22.01.1999	ggu	oxygen modules introduced
c 04.01.2000	ggu	cstset -> cstcheck, new cstsetup
c 24.10.2001	ggu	Write on use of Coriolis
c 14.08.2003	ggu	set depth values transferred to newini
c 23.03.2006	ggu	use routine set_timeunit() to set time unit
c 11.02.2009	ggu	in cstfile set fixed name for STR file
c 23.03.2010	ggu	changed v6.1.1
c 09.04.2010	ggu	changed v6.1.3
c 28.09.2010	ggu	changed VERS_6_1_11
c 15.07.2011	ggu	new call to read basin
c 04.11.2011	ggu	changed VERS_6_1_35
c 14.02.2012	ggu	changed VERS_6_1_44
c 07.03.2014	ggu	changed VERS_6_1_72
c 18.07.2014	ggu	changed VERS_7_0_1
c 20.10.2014	ggu	new routine ckdate()
c 10.11.2014	ggu	time and date routines transfered to subtime.f
c 26.11.2014	ggu	changed VERS_7_0_7
c 23.12.2014	ggu	changed VERS_7_0_11
c 05.06.2015	ggu	changed VERS_7_1_12
c 10.07.2015	ggu	changed VERS_7_1_50
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 30.07.2015	ggu	read str-file from command line
c 31.07.2015	ggu	changed VERS_7_1_84
c 09.05.2017	ggu	changed VERS_7_5_26
c 05.10.2017	ggu	cstfile called with file name
c 05.12.2017	ggu	changed VERS_7_5_39
c 19.04.2018	ggu	changed VERS_7_5_45
c 16.10.2018	ggu	changed VERS_7_5_50
c 16.02.2019	ggu	changed VERS_7_5_60
c 13.03.2019	ggu	changed VERS_7_5_61
c 21.05.2019	ggu	changed VERS_7_5_62
c
c************************************************************************

	subroutine cstinit

c set up constants and parameters

	implicit none

	include 'modules.h'

	include 'mkonst.h'
	include 'pkonst.h'

c parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	call nlsinh
	call check_parameter_values('after nlsinh')

c default names %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	call fnminh
	call check_parameter_values('after fnminh')

c /mkonst/ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	eps1=1.0E-5
	eps2=1.0E-6
	pi=3.141592653
	flag=-9988765.0
	high=1.0E+35
        higi=.21474e+10   ! +/- 2147483647

c /pkonst/ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c these parameters are set one more time in cstset (from STR file)

	grav=9.81
	fcor=0.
	dcor=0.
	dirn=0.
	rowass=1025.
	roluft=1.225

c other modules %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c	call inexta
	call inflxa
	call invola
	call inarea
	call inbnds

	call inclos
c	call inoxy	!oxygen
c	call inlgr	!float tracking

	call modules(M_INIT)

c other parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	call putpar('flag',flag)

	end

c*******************************************************************

	subroutine cstcheck

c sets up and checks parameters read

	use mod_bnd
	use mod_bound_geom
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'modules.h'

	include 'param.h'
	include 'mkonst.h'
	include 'pkonst.h'

	integer ibarcl
        real getpar
        integer nkbnd

c parameters read from STR file

	grav=getpar('grav')
	rowass=getpar('rowass')
	roluft=getpar('roluft')

	call setup_date	!set up and check date parameter
	call setup_time	!set up and check and correct time parameters

c	call ckexta	!extra output points
	call ckflxa	!flux sections
	call ckvola	!flux sections
	call ckarea	!chezy values

c re-allocate boundary arrays

        call mod_bound_geom_reinit(nkn,nrb)
	call mod_bnd_reinit(nbc)
	call ckbnds	!boundary conditions

	call ckclos
c	call ckoxy	!oxygen

	ibarcl=nint(getpar('ibarcl'))
	if( ibarcl == 0 ) then
	  call putpar('itemp',0.)
	  call putpar('isalt',0.)
	end if

	call modules(M_CHECK)

	call ckcori	!sets dlat for coriolis parameter

	end

c*******************************************************************

	subroutine cstsetup

c sets up modules

	implicit none

	include 'modules.h'

	call flxini
	call volini

	write(6,*) 'cstsetup: setting up modules'
	call modules(M_SETUP)
	write(6,*) 'cstsetup: finished setting up modules'

	end

c********************************************************************

	subroutine cstfile(strfile)

c reads files (str and bas)

	use basin

	implicit none

	character*80 strfile

	integer nin,nc
	character*80 basnam

	integer idefbas,ifileo

	if( strfile == ' ' ) then
	  write(6,*) 'Usage: shyfem str-file'
	  stop 'error stop cstfile: internal error (1)'
	else
	  nin = ifileo(0,strfile,'form','old')
	  if( nin < 1 ) then
	    write(6,*) 'error opening STR file: ',trim(strfile)
	    stop 'error stop cstfile: no such STR file'
	  end if
	end if

	call nlsh2d(nin)
	if( nin .ne. 5 ) close(nin)

        call getfnm('basnam',basnam)
	write(6,*) 'reading basin: ',trim(basnam)
        nin = idefbas(basnam,'old')
	call basin_read(nin)
	close(nin)

	end

c********************************************************************

	subroutine ckcori

c set coriolis parameter

	implicit none

	include 'pkonst.h'

	real dlat

	real getpar

	call bas_get_geom(dcor,dirn)	! from basin

	dlat=getpar('dlat')

	if( dlat .le. 90. ) dcor = dlat
	call putpar('dlat',dcor)

	end

c********************************************************************

