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
        module constants_par
!------------------------------------------------------------------------
        contains
!------------------------------------------------------------------------

	subroutine cstinit

! set up constants and parameters

        use def_para

	implicit none

	include 'mkonst.h'
	include 'pkonst.h'

! parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	call nlsinh
	call check_parameter_values('after nlsinh')

! default names %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	call fnminh
	call check_parameter_values('after fnminh')

! /mkonst/ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        eps1=1.0E-5
        eps2=1.0E-6
        pi=3.141592653
        flag=-9988765.0
        high=1.0E+35
        higi=.21474e+10   ! +/- 2147483647

! /pkonst/ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! these parameters are set one more time in cstset (from STR file)

        grav=9.81
        fcor=0.
        dcor=0.
        dirn=0.
        rowass=1025.0
        roluft=1.225

! other parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	call putpar('const',flag)
	call putpar('flag',flag)

	end

!********************************************************************

	subroutine cstfile

! reads files (str and bas)

	use basin
        use fil
        use para
        use defnames
        use iostr_par

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

!------------------------------------------------------------------------
        end module constants_par
!------------------------------------------------------------------------
