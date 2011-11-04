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
c subroutine cstfile            reads files (str and bas)
c
c subroutine cktime		check time parameters
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
c
c************************************************************************

	subroutine cstinit

c set up constants and parameters

	implicit none

	include 'modules.h'

        real eps1,eps2,pi,flag,high,higi
        common /mkonst/ eps1,eps2,pi,flag,high,higi
	real grav,fcor,dcor,dirn,rowass,roluft
	common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft

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
	call inwnds

	call inclos
c	call inoxy	!oxygen
c	call inlgr	!float tracking

	call modules(M_INIT)

c other parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	call putpar('const',flag)
	call putpar('flag',flag)

	end

c*******************************************************************

	subroutine cstcheck

c checks parameters read
c
c revised 07.04.95 by ggu !$$conzfl - conz compared to iflag (bug)
c revised 07.04.95 by ggu !$$baroc - impl. of baroclinic salt/temp (21/22)
c revised 13.06.97 by ggu !$$kranf - check if kranf <= krend

	implicit none

	include 'modules.h'

        real eps1,eps2,pi,flag,high,higi
        common /mkonst/ eps1,eps2,pi,flag,high,higi
	real grav,fcor,dcor,dirn,rowass,roluft
	common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft
	real getpar

c parameters read from STR file

	grav=getpar('grav')
	rowass=getpar('rowass')
	roluft=getpar('roluft')

c	call ckexta	!extra output points
	call ckflxa	!flux sections
	call ckvola	!flux sections
	call ckarea	!chezy values
	call ckbnds	!boundary conditions
	call ckwnds	!wind

	call ckclos
c	call ckoxy	!oxygen

	call modules(M_CHECK)

	call cktime	!check and correct time parameters
	call ckcori	!set coriolis parameter

	end

c*******************************************************************

	subroutine cstsetup

c sets up modules

	implicit none

	include 'modules.h'

	call flxini
	call volini

	call modules(M_SETUP)

	end

c********************************************************************

	subroutine cstfile(nkndim,neldim)

c reads files (str and bas)

	implicit none

	integer nkndim,neldim
	integer nin
	character*80 file,basnam

	integer idefbas

	nin = 5
	call nlsh2d(nin)

        call getfnm('basnam',basnam)
        nin = idefbas(basnam,'old')
	call sp13rr(nin,nkndim,neldim)
	close(nin)

	end

c********************************************************************

	subroutine cktime

c check time parameters

	implicit none

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

	integer itunit
	integer iround
	real getpar

	itanf = iround(getpar('itanf'))
	itend = iround(getpar('itend'))
	idt = iround(getpar('idt'))

	itunit = iround(getpar('itunit'))
	call set_timeunit(itunit)

	if( idt .le. 0 .or. itanf+idt .gt. itend ) then
	   write(6,*) 'Error in compulsory time parameters'
	   write(6,*) 'itanf,itend,idt :',itanf,itend,idt
	   stop 'error stop : cktime'
	end if

	niter=0
	it=itanf
	nits=(itend-itanf)/idt

	end

c********************************************************************

	subroutine ckcori

c set coriolis parameter

	implicit none

        real eps1,eps2,pi,flag,high,higi
        common /mkonst/ eps1,eps2,pi,flag,high,higi
	real grav,fcor,dcor,dirn,rowass,roluft
	common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft

	integer icor
	real dlat

	integer iround
	real getpar

	icor=iround(getpar('icor'))
	dlat=getpar('dlat')

	if(icor.ne.0) then
	   if( dlat .le. 90. ) dcor = dlat
	   fcor = 2.0 * 0.729E-4 * sin(dcor*pi/180.)
	   call putpar('dlat',dcor)
           write(6,*) 'Coriolis: ',icor,dcor,dlat,fcor
	else
	   fcor = 0.
           write(6,*) 'No Coriolis accelaration'
	end if

	end

c********************************************************************

