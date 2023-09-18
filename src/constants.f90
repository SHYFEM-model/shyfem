!
! $Id: subcst.f,v 1.19 2009-02-13 17:22:44 georg Exp $
!
! routines to set up constants
!
! contents :
!
! subroutine cstinit            set up constants and parameters
! subroutine cstcheck           checks parameters read
! subroutine cstsetup           sets up modules
! subroutine cstfile            reads files (str and bas)
!
! subroutine cktime		check time parameters
! subroutine ckdate		check date parameter
! subroutine ckcori		set coriolis parameter
!
! revision log :
!
! revised 05.08.92 by ggu   $$ibtyp3 - implementation of ibtyp=3
! revised 31.08.92 by ggu   $$impli - implicit time step
! revised 23.11.92 by ggu   $$ibtyp11 - implementation of ibtyp=11,51
! revised 27.10.93 by ggu   $$roger - implementation of ibtyp=70 (nsea)
! revised 02.06.97 by ggu   $$EXTINW - extinw changed to ipint
! revised 29.06.97 by ggu   no cstdim in file
! 29.04.1998    ggu     module for semi-implicit time-step in own routine
! 12.08.1998    ggu     new parameter dlat -> specify latitude for coriolis
! 03.09.1998    ggu     call bocche to adjust depth at Venice inlets
! 06.11.1998    ggu     call to huniqu to set up hkv and hev
! 22.01.1999    ggu     oxygen modules introduced
! 04.01.2000    ggu     cstset -> cstcheck, new cstsetup
! 24.10.2001    ggu     Write on use of Coriolis
! 14.08.2003    ggu     set depth values transferred to newini
! 23.03.2006    ggu     use routine set_timeunit() to set time unit
! 11.02.2009    ggu     in cstfile set fixed name for STR file
! 15.07.2011    ggu     new call to read basin
! 20.10.2014    ggu     new routine ckdate()
! 10.11.2014    ggu     time and date routines transfered to subtime.f
! 30.07.2015    ggu     read str-file from command line
!
!************************************************************************
!------------------------------------------------------------------------
        module constants
!------------------------------------------------------------------------
        contains
!------------------------------------------------------------------------

	subroutine cstinit

! set up constants and parameters

        use para
        use bnd_admin
        use chezy
        use flux
        use volume
        use modls
        use closing
        use def_para

	implicit none

	include 'modules.h'

	include 'mkonst.h'
	include 'pkonst.h'

! parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	call nlsinh
	call check_parameter_values('after nlsinh')

! default names %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	call fnminh
	call check_parameter_values('after fnminh')

! /mkonst/ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        eps1=1.0D-5
        eps2=1.0D-6
        pi=3.141592653d0
        flag=-9988765.0d0
        high=1.0D+35
        higi=.21474d+10   ! +/- 2147483647

! /pkonst/ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! these parameters are set one more time in cstset (from STR file)

        grav=9.81d0
        fcor=0.d0
        dcor=0.d0
        dirn=0.d0
        rowass=1025.0d0
        roluft=1.225d0

! other modules %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!	call inexta
	call inflxa
	call invola
	call inarea
	call inbnds

	call inclos
!	call inoxy	!oxygen
!	call inlgr	!float tracking

	call modules(M_INIT)

! other parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	call putpar('const',flag)
	call putpar('flag',flag)

	end

!*******************************************************************

	subroutine cstcheck

! sets up and checks parameters read
!
! revised 07.04.95 by ggu !$$conzfl - conz compared to iflag (bug)
! revised 07.04.95 by ggu !$$baroc - impl. of baroclinic salt/temp (21/22)
! revised 13.06.97 by ggu !$$kranf - check if kranf <= krend

	use mod_bnd
	use bnd_geom
	use basin, only : nkn,nel,ngr,mbw
	use shympi
        use para
        use bnd_admin
        use chezy
        use flux
        use volume
        use modls
        use closing
        use time_admin

	implicit none

	include 'modules.h'

	include 'param.h'
	include 'mkonst.h'
	include 'pkonst.h'

! parameters read from STR file

	grav=getpar('grav')
	rowass=getpar('rowass')
	roluft=getpar('roluft')

	call setup_date	!set up and check date parameter
	call setup_time	!set up and check and correct time parameters

!	call ckexta	!extra output points
	call ckflxa	!flux sections
	call ckvola	!flux sections
	call ckarea	!chezy values

! re-allocate boundary arrays
	if(bmpi) then
          call mod_bnd_geom_reinit(nkn,nrb,domain%nodes%totalID)
	else
          call mod_bnd_geom_reinit(nkn,nrb)
	end if
	call mod_bnd_reinit(nbc)
	call ckbnds	!boundary conditions

	call ckclos
!	call ckoxy	!oxygen

	call modules(M_CHECK)

	call ckcori	!sets dlat for  coriolis parameter

	end

!*******************************************************************

	subroutine cstsetup

! sets up modules

        use flux
        use volume
        use modls

	implicit none

	include 'modules.h'

	call flxini
	call volini

	write(6,*) 'cstsetup: setting up modules'
	call modules(M_SETUP)
	write(6,*) 'cstsetup: finished setting up modules'

	end

!********************************************************************

	subroutine cstfile

! reads files (str and bas)

	use basin
        use fil
        use defnames
        use para
        use iostr

	implicit none

	integer nin,nc
	character*80 basnam
	character*80 strfile

	strfile = ' '
        nc = command_argument_count()
        if( nc .gt. 0 ) call get_command_argument(1,strfile)
 
	if( strfile == ' ' ) then
	  write(6,*) 'Usage: shyfem str-file'
	  stop
	  nin = 5
	else
	  nin = ifileo(0,strfile,'form','old')
	end if

	call nlsh2d(nin)
	if( nin .ne. 5 ) close(nin)

        call getfnm('basnam',basnam)
        nin = idefbas(basnam,'old')
	call basin_read(nin)
	close(nin)

	end

!********************************************************************

	subroutine ckcori

! set coriolis parameter

        use para
        use basin

	implicit none

	include 'pkonst.h'

	double precision dlat

	call bas_get_geom(dcor,dirn)	! from basin

	dlat=getpar('dlat')

	if( dlat .le. 90. ) dcor = dlat
	call putpar('dlat',dcor)

	end

!********************************************************************

!------------------------------------------------------------------------
        end module constants
!------------------------------------------------------------------------
