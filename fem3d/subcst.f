c
c $Id: subcst.f,v 1.19 2009-02-13 17:22:44 georg Exp $
c
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
c revised 05.08.92 by ggu   $$ibtyp3 - implementation of ibtyp=3
c revised 31.08.92 by ggu   $$impli - implicit time step
c revised 23.11.92 by ggu   $$ibtyp11 - implementation of ibtyp=11,51
c revised 27.10.93 by ggu   $$roger - implementation of ibtyp=70 (nsea)
c revised 02.06.97 by ggu   $$EXTINW - extinw changed to ipint
c revised 29.06.97 by ggu   no cstdim in file
c 29.04.1998    ggu     module for semi-implicit time-step in own routine
c 12.08.1998    ggu     new parameter dlat -> specify latitude for coriolis
c 03.09.1998    ggu     call bocche to adjust depth at Venice inlets
c 06.11.1998    ggu     call to huniqu to set up hkv and hev
c 22.01.1999    ggu     oxygen modules introduced
c 04.01.2000    ggu     cstset -> cstcheck, new cstsetup
c 24.10.2001    ggu     Write on use of Coriolis
c 14.08.2003    ggu     set depth values transferred to newini
c 23.03.2006    ggu     use routine set_timeunit() to set time unit
c 11.02.2009    ggu     in cstfile set fixed name for STR file
c 15.07.2011    ggu     new call to read basin
c 20.10.2014    ggu     new routine ckdate()
c 10.11.2014    ggu     time and date routines transfered to subtime.f
c 30.07.2015    ggu     read str-file from command line
c 05.10.2017    ggu     cstfile called with file name
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

	call putpar('zconst',flag)
	call putpar('flag',flag)

	end

c*******************************************************************

	subroutine cstcheck

c sets up and checks parameters read
c
c revised 07.04.95 by ggu !$$conzfl - conz compared to iflag (bug)
c revised 07.04.95 by ggu !$$baroc - impl. of baroclinic salt/temp (21/22)
c revised 13.06.97 by ggu !$$kranf - check if kranf <= krend

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

