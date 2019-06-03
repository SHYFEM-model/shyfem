
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

c bnd administration routines
c
c contents :
c
c subroutine inbnds
c subroutine rdbnds(ibc)
c subroutine ckbnds
c subroutine prbnds
c subroutine tsbnds
c
c function zvbnds(ibc)			 returns value of open boundary
c subroutine stybnd(ibc,ibtyp)		 sets type of open boundary
c subroutine infobnd(ibc,ibtype,...)	 returns important info on boundary ibc
c function itybnd(ibc)			 returns type of open boundary
c function nbnds			 returns total number of open b.
c function nkbnds(ibc)			 returns total number of nodes of ibc
c function kbnd(i)			 returns i th node of all b. nodes
c subroutine kanfend(ibc,kranf,krend)    returns index of first and last bnode
c function kbnds(ibc,i)			 returns i th node of boundary ibc
c subroutine irbnds(ibc,ndim,idim,nodes) returns nodes of boundary ibc
c subroutine setget_bnd_par(ibc,ientry,value,bset) sets/gets value at ientry
c subroutine setbnd(ibc,value,barray)	 sets boundary ibc to value in barray
c
c subroutine setbc(value,array,flag)	 sets all open boundaries to value
c subroutine chkibc(ibc,errtext)	 checks if ibc is in bounds
c
c notes :
c
c for scalars: -999. uses ambient value
c
c what to do when adding a new parameter to boundary section:
c	add parameter in rdbnds() with description
c	add parameter to array in mod_bnd.f
c what to do when adding a new file name to boundary section:
c	add file name in rdbnds() with description
c	add file name to array in mod_bnd.f
c	if needed add file name in get_boundary_file()
c
c revision log :
c
c 01.08.1997	ggu	$$1stnode - first boundary node was not registered
c 23.09.1997	ggu	introduced conzn -> name of file for conz values
c 30.09.1997	ggu	write error in prbnds
c 03.12.1997	ggu	introduced tempn,saltn (see conzn)
c 27.03.1998	ggu	new utility routines
c 25.05.1998	ggu	documentation (DOCS)
c 20.06.1998	ggu	documentation for momentum input
c 13.07.1998	ggu	implemented ibtyp = 4
c 23.07.1998	ggu	documentation
c 24.08.1998	ggu	BC for concentration is bnd(20,..)
c 24.08.1998	ggu	BC for maximum input level is bnd(12,..) -> levmax
c 27.08.1998	ggu	accessor function levbnd for levmax
c 22.01.1999	ggu	new subroutine setbnd
c 08.07.1999	ggu	bug with iqual -> ibtyp not respected
c 20.01.2000	ggu	call to rdbnds without dimension -> use getdim
c 07.04.2000	ggu	new subroutine setbc
c 07.05.2001	ggu	introduced new variable zfact
c 25.09.2001	ggu	introduced bio2dn
c 07.08.2003	ggu	check for nrb in rdbnds
c 15.10.2004	ggu	new boundary types and sedin
c 02.03.2005	ggu	new nbdim for 3D boundary values
c 02.03.2005	ggu	some new helper functions
c 07.11.2005	ccf	introduced sed2dn
c 16.02.2006	ggu	introduced tox3dn
c 07.04.2008	aac	introduced bfm1bc bfm2bc bfm3bc OB condition for ERSEM
c 17.04.2008	ggu	deleted infobnd(), levbnd()
c 28.04.2008	ggu	call to nrdpar in double precision
c 29.04.2008	ggu&aac	new boundary arrays for ERSEM
c 30.05.2008	ggu	eliminated numbers for parameters
c 03.06.2008	ggu	new parameters levmin, kref
c 06.06.2008	ggu	completely restructured
c 02.04.2009	ggu	intpol default is 0, some unused routines deleted
c 20.04.2009	ggu	new variable ztilt, ndim substituted with nbvdim
c 23.03.2010	ggu	changed v6.1.1
c 28.09.2010	ggu	changed VERS_6_1_11
c 17.02.2011	ggu	changed VERS_6_1_18
c 23.02.2011	ggu	new parameters tramp and levflx implemented
c 01.03.2011	ggu	changed VERS_6_1_20
c 21.06.2012	ggu&aar	new file names for mud module
c 26.06.2012	ggu	changed VERS_6_1_55
c 29.11.2013	ggu	allow for non continous boundary numbering
c 05.12.2013	ggu	changed VERS_6_1_70
c 28.01.2014	ggu	changed VERS_6_1_71
c 28.03.2014	ggu	new parameter lgrpps
c 05.05.2014	ggu	changed VERS_6_1_74
c 16.06.2014	ggu	new include file bnd.h (with nbndim)
c 25.06.2014	ggu	new routine exists_bnd_name()
c 21.10.2014	ggu	changed VERS_7_0_3
c 29.10.2014	ccf	include vel3dn boundary file
c 03.11.2014	ggu	nbdim deleted
c 23.12.2014	ggu	changed VERS_7_0_11
c 19.01.2015	ggu	changed VERS_7_1_3
c 23.06.2015	ggu	setbc() deleted, nrz,nrq eliminated
c 10.07.2015	ggu	changed VERS_7_1_50
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 30.07.2015	ggu	changed VERS_7_1_83
c 05.11.2015	ggu	changed VERS_7_3_12
c 18.12.2015	ggu	changed VERS_7_3_17
c 14.01.2016	ggu	check of boundaries considers mpi subdomains
c 15.02.2016	ggu	check if boundary is given twice
c 19.02.2016	ggu	changed VERS_7_5_3
c 22.02.2016	ggu	new files bfmbcn integrated
c 01.04.2016	ggu	restructured - arrays transfered to mod_bnd.f
c 15.04.2016	ggu	changed VERS_7_5_8
c 07.06.2016	ggu	changed VERS_7_5_12
c 30.09.2016	ggu	changed VERS_7_5_18
c 31.03.2017	ggu	changed VERS_7_5_24
c 13.04.2017	ggu	use array feature of para for kbound
c 05.12.2017	ggu	changed VERS_7_5_39
c 03.04.2018	ggu	changed VERS_7_5_43
c 19.04.2018	ggu	changed VERS_7_5_45
c 25.10.2018	ggu	changed VERS_7_5_51
c 18.12.2018	ggu	changed VERS_7_5_52
c 16.02.2019	ggu	changed VERS_7_5_60
c 13.03.2019	ggu	changed VERS_7_5_61
c
c************************************************************************

	subroutine inbnds

c initializes boundary parameters

	use mod_bnd

	implicit none

	nbc = 0
	nrb = 0

	end

c************************************************************************

	subroutine rdbnds(ibc)

c reads boundary info from STR file

	use mod_bnd
	use mod_bound_geom
	use nls
	use para

	implicit none

	integer ibc

	character*80 name,text,file
	double precision dvalue
	real value
	integer n
	integer i,kranf,krend,kref
	integer iweich,id,nbnd
	integer nrdpar
	logical ball

	real getpar

!------------------------------------------------------
! add new elements to boundary arrays
!------------------------------------------------------

	n = max(nbc,ibc)	!allow for non contiguous numbering
	call mod_bnd_adjust(n)
	do i=nbc+1,n
	  call init_bnd_par(i)
	  call init_bnd_file(i)
	end do
	nbc = n

!------------------------------------------------------
! check if boundary section is unique
!------------------------------------------------------

        call get_bnd_ipar(ibc,'kranf',kranf)
	if( kranf > 0 ) then
	  write(6,*) 'ibc = ',ibc
	  stop 'error stop rdbnds: boundary defined twice'
	end if

!------------------------------------------------------
! start initializing boundary section
!------------------------------------------------------

	call sctpar('bound')
	call sctfnm('bound')

c DOCS	START	S_bound
c
c These parameters determine the open boundary nodes and the type
c of the boundary: level or flux boundary. At the first the water levels
c are imposed, on the second the fluxes are prescribed.
c
c There may be multiple sections |bound| in one parameter input file,
c describing all open boundary conditions necessary. Every section
c must therefore be supplied with a boundary number. The numbering
c of the open boundaries must
c be increasing. The number of the boundary must be specified
c directly after the keyword |bound|, such as |bound1| or |bound 1|.
c
c |kbound|	Array containing the node numbers that are part of the
c		open boundary. The node numbers must form one contiguous
c		line with the domain (elements) to the left. This
c		corresponds to an anti-clockwise sense. The type of
c		boundary depends on the	value of |ibtyp|. In case this value
c		is 1 or 2 at least two nodes must be given.

	call para_add_array_value('kbound',0.)

c |ibtyp|	Type of open boundary. 
c		\begin{description}
c		\item[0] No boundary values specified
c		\item[1] Level boundary. At this open boundary
c			 the water level is imposed and the prescribed
c			 values are interpreted as water levels in meters.
c			 If no value for |ibtyp| is specified this
c			 is the default.
c		\item[2] Flux boundary. Here the discharge in \dischargeunit
c			 has to be prescribed.
c		\item[3] Internal flux boundary. As with |ibtyp = 2| a
c			 discharge has to be imposed, but the node where
c			 discharge is imposed can be an internal node
c			 and need not be on the outer boundary of
c			 the domain. For every node in |kbound| the
c			 volume rate specified will be added to the
c			 existing water volume. This behavior is different
c			 from the |ibtyp = 2| where the whole boundary
c			 received the discharge specified.
c		\item[4] Momentum input. The node or nodes may be internal.
c			 This feature can be used to describe local 
c			 acceleration of the water column. 
c			 The unit is force / density [\maccelunit].
c			 In other words it is the rate of volume 
c			 [\dischargeunit] times the velocity [m/s] 
c			 to which the water is accelerated.
c		\end{description}
c |iqual|	If the boundary conditions for this open boundary
c		are equal to the ones of boundary |i|, then
c		setting |iqual = i| copies all the values of
c		boundary |i| to the actual boundary. Note that the
c		value of |iqual| must be smaller than the number
c		of the actual boundary, i.e., boundary |i| must have
c		been defined before. (This feature is temporarily
c		not working; please do not use.)

	call addpar('ibtyp',1.)
	call addpar('iqual',0.)

c The next parameters give a possibility to specify the file name
c of the various input files that are to be read by the model.
c Values for the boundary condition can be given at any time step.
c The model interpolates in between given time steps if needed. The
c grade of interpolation can be given by |intpol|.
c
c All files are in ASCII and share a common format.
c The file must contain two columns, the first giving the
c time of simulation in seconds that refers to the value
c given in the second column. The value in the second
c column must be in the unit of the variable that is given.
c The time values must be in increasing order.
c There must be values for the whole simulation,
c i.e., the time value of the first line must be smaller
c or equal than the start of the simulation, and the time
c value of the last line must be greater or equal than the
c end of the simulation.

c |boundn|	File name that contains values for the boundary condition.
c		The value of the variable given in the second column
c		must be in the unit determined by |ibtyp|, i.e.,
c		in meters for a level boundary, in \dischargeunit for
c		a flux boundary and in \maccelunit for a momentum
c		input.
c |zfact|	Factor with which the values from |boundn|
c		are multiplied to form the final value of the
c		boundary condition. E.g., this value can be used to
c		set up a quick sensitivity run by multiplying
c		all discharges by a factor without generating
c		a new file. (Default 1)

	call addfnm('boundn',' ')
	call addpar('zfact',1.)

c |levmin, levmax|	A point discharge normally distributes its discharge
c		over the whole water column. If it is important that in
c		a 3D simulation the water mass discharge is concentrated 
c		only in some levels, the parameters |levmin| and |levmax|
c		can be used. They indicate the lowest and deepest level over
c		which the discharge is distributed. Default values are 0, which
c		indicate that the discharge is distributed over the
c		whole water column. Setting only |levmax| distributes from
c		the surface to this level, and setting only |levmin|
c		distributes from the bottom to this level.

	call addpar('levmin',0.)
	call addpar('levmax',0.)

c |conzn, tempn, saltn|	File names that contain values for the respective
c			boundary condition, i.e., for concentration,
c			temperature and salinity. The format is the same
c			as for file |boundn|. The unit of the values
c			given in the second column must the ones of the
c			variable, i.e., arbitrary unit for concentration,
c			centigrade for temperature and psu (per mille)
c			for salinity.

	call addfnm('conzn',' ')
	call addfnm('tempn',' ')
	call addfnm('saltn',' ')

c |vel3dn|	File name that contains current velocity values for the 
c		boundary condition.  The format is the same as for file 
c		|tempn| but it has two variables: 
c	 	current velocity in x and current velocity in y.
c		Velocity can be nudged or imposed depending on the value 
c		of |tnudge| (mandatory). The unit is [m/s].

	call addfnm('vel3dn',' ')

c |tnudge|	Relaxation time for nudging of boundary velocity.
c		For |tnudge| = 0. velocities are imposed, for
c		|tnudge| > 0. velocities are nudged. The
c		default is -1 which means do nothing. Unit is [s].
c		(Default -1)

        call addpar('tnudge',-1.)

c The next variables specify the name of the boundary value file
c for different modules. Please refer to the documentation of the
c single modules for the units of the variables.

c |bio2dn|	File name that contains values for the ecological
c		module (EUTRO-WASP).
c |sed2dn|	File name that contains values for the sediment
c		transport module.
c		The unit of the values given
c		in the second and following columns (equal to the 
c		number of defined grainsize in parameter |sedgrs|).

c |mud2dn|	File name that contains values for the fluid mud
c		module.
c |lam2dn|	File name that contains values for the fluid mud
c		module (boundary condition for the structural parameter, 
c		to be implemented).
c |dmf2dn|	File name that contains values for the fluid mud
c		module (boundary conditions for the advection of flocsizes,
c		to be implemented).
c |tox3dn|	File name that contains values for the toxicological
c		module.

	call addfnm('bio2dn',' ')
	call addfnm('sed2dn',' ')
	call addfnm('mud2dn',' ')
	call addfnm('lam2dn',' ')
	call addfnm('dmf2dn',' ')
	call addfnm('tox3dn',' ')

cc File name for OB condition in ERSEM MODULE - undocumented
cc ... will be removed sooner or later ...

	call addfnm('bfm1bc',' ')
	call addfnm('bfm2bc',' ')
	call addfnm('bfm3bc',' ')

c |bfmbcn|	File name that contains values for the bfm module.

	call addfnm('bfmbcn',' ')

c |mercn|	File name that contains values for the mercury module.

	call addfnm('mercn',' ')

c |intpol|	Order of interpolation for the boundary values read
c		in files. Use for 1 for stepwise (no) interpolation,
c		2 for linear and 4 for cubic interpolation. 
c		The default is linear interpolation, except for
c		water level boundaries (|ibtyp=1|) where cubic
c		interpolation is used.

	call addpar('intpol',0.)

c The next parameters can be used to impose a sinusoidal water level
c (tide) or flux at the open boundary. These values are used if no
c boundary file |boundn| has been given. The values must be in the unit
c of the intended variable determined by |ibtyp|.

c |ampli|	Amplitude of the sinus function imposed. (Default 0)
c |period|	Period of the sinus function. (Default 43200, 12 hours)
c |phase|	Phase shift of the sinus function imposed. A positive value
c		of one quarter of the period reproduces a cosine
c		function. (Default 0)
c |zref|	Reference level of the sinus function imposed. If only
c		|zref| is specified (|ampli = 0|) a constant value
c		of |zref| is imposed on the open boundary.

	call addpar('ampli',0.)
	call addpar('period',43200.)
	call addpar('phase',0.)
	call addpar('zref',0.)

c With the next parameters a constant value can be imposed for the 
c variables of concentration, temperature and salinity. In this case
c no file with boundary values has to be supplied. The default for all
c values is 0, i.e., if no file with boundary values is supplied and
c no constant is set the value of 0 is imposed on the open boundary.
c A special value of -999 is also allowed. In this case the value
c imposed is the ambient value of the parameter close to the boundary.

c |conz, temp, salt|	Constant boundary values for concentration,
c			temperature and salinity respectively. If these
c			values are set no boundary file has to be supplied.
c			(Default 0)

	call addpar('conz',0.)			!$$conz
	call addpar('temp',0.)			!$$baroc
	call addpar('salt',0.)			!$$baroc

c The next two values are used for constant momentum input. 
c This feature can be used to describe local acceleration of the
c water column. The values give the input of momentum 
c in x and y direction. The unit is force / density (\maccelunit).
c In other words it is the rate of volume (\dischargeunit) times
c the velocity (m/s) to which the water is accelerated.
c
c These values are used if 
c boundary condition |ibtyp = 4| has been chosen and
c no boundary input file has been given.
c If the momentum input is varying then it may be specified with
c the file |boundn|. In this case the file |boundn| must contain
c three columns, the first for the time, and the other two for
c the momentum input in $x,y$ direction.
c
c Please note that this feature is temporarily not available.
c
c |umom, vmom|		Constant values for momentum input. (Default 0)

	call addpar('umom',0.)
	call addpar('vmom',0.)

c The next two values can be used
c to achieve the tilting of the open boundary if only one water level value
c is given. If only |ktilt| is given then the boundary values
c are tilted to be in equilibrium with the Coriolis force. This may avoid
c artificial currents along the boundary. |ktilt| must be a boundary node
c on the boundary.
c
c If |ztilt| is given the tilting of the boundary is explicitly set
c to this value. The tilting of the first node of the boundary is set 
c to $-|ztilt|$
c and the last one to $+|ztilt|$. The total amount of tilting is
c therefore is $2 \cdot |ztilt|$. If |ktilt| is not specified
c then a linear interpolation between the first and the last boundary
c node will be carried out. If also |ktilt| is specified then
c the boundary values are arranged that the water levels are 
c tilted around |ktilt|, e.g., $-|ztilt|$ at the first boundary node,
c 0 at |ktilt|, and $+|ztilt|$ at the last boundary node.
c
c |ktilt|		Node of boundary around which tilting should
c			take place. (Default 0, i.e., no tilting)
c |ztilt|		Explicit value for tilting (unit meters).
c			(Default 0)

	call addpar('ktilt',0.)
	call addpar('ztilt',0.)

c Other parameters:

c |igrad0|		If different from 0 a zero gradient boundary
c			condition will be implemented. This is already the
c			case for scalars under outflowing conditions. However,
c			with |igrad0| different from 0 this conditions
c			will be used also for inflow conditions. (Default 0)

	call addpar('igrad0',0.)	!use 0 gradient for scalars

c |tramp|		Use this value to start smoothly a discharge
c			boundary condition. If set it indicates the
c			time (seconds) that will be used to increase
c			a discharge from 0 to the desired value (Default 0)

	call addpar('tramp',0.)		!start smoothly for discharge

c |levflx|		If discharge is depending on the water level
c			(e.g., lake outflow) then this parameter indicates to
c			use one of the possible outflow curves. Please
c			note that the flow dependence on the water level
c			must be programmed in the routine $|level\_flux()|$.
c			(Default 0)

	call addpar('levflx',0.)	!use level-discharge relationship

c |nad|			On the open boundaries it is sometimes convenient
c			to not compute the non-linear terms in the momentum
c			equation because instabilities may occur. Setting 
c			the parameter |nad| to a value different from 0
c			indicates that in the first |nad| nodes from the
c			boundary the non linear terms are switched off.
c			(Default 0)

	call addpar('nad',-1.)		!no advective terms for this boundary

c |lgrpps|		Indicates the number of particles released at
c			the boundary for the lagrangian module. If positive
c			it is the number of particles per second released
c			along the boundary. If negative its absolute
c			value indicates the particles per volume flux
c			(unit \dischargeunit) released along the boundary.
c			(Default 0)

	call addpar('lgrpps',0.)	!particles per second for lagrange
					!if negative parts per volume flux

c DOCS	END

cc undocumented
	call addpar('kref',0.)		!not working...
c here add dummy variables
	call addpar('zval',0.)
	call addpar('kmanf',0.)
	call addpar('kmend',0.)

!------------------------------------------------------
! start reading loop
!------------------------------------------------------

	kranf = nrb + 1
	call set_bnd_ipar(ibc,'kranf',kranf)	!position of starting node
	call addpar('kranf',float(kranf))

	iweich=1
	do while(iweich.ne.0)
	    iweich=nls_insert_variable('bound',name,dvalue,text)
	    value = dvalue
	    if( iweich .lt. 0 ) goto 92
	    if( iweich .eq. 2 .and. name .ne. 'kbound' ) goto 93
	    if( iweich .eq. 4 ) goto 93
	    if( name .eq. 'kbound' ) then	!$$1stnode
		!nrb=nrb+1
		!call mod_irv_init(nrb)
		!irv(nrb)=nint(value)
	    else if( iweich .eq. 3 ) then	!file name
                call set_bnd_file(ibc,name,text)
	    !else if( iweich .ne. 0 ) then
	    else if( iweich .eq. 1 ) then
                call set_bnd_par(ibc,name,value)
	    end if
	end do

	call para_get_array_size('kbound',n)
	call mod_irv_init(nrb+n)
	call para_get_array_value('kbound',n,n,irv(nrb+1:))
	nrb=nrb+n

	krend = nrb
	call set_bnd_ipar(ibc,'krend',krend)	!position of end node
	call addpar('krend',float(krend))

!------------------------------------------------------
! copy all values and files to private boundary arrays
!------------------------------------------------------

	call copy_bnd_par(ibc)
	call copy_bnd_file(ibc)

!------------------------------------------------------
! write out all values of boundary section (only for debug)
!------------------------------------------------------

	ball = .true.
        !call check_bnd_par_entries(ibc,ball)
        !call check_bnd_file_entries(ibc,ball)

!------------------------------------------------------
! delete public section
!------------------------------------------------------

	call check_parameter_values('before deleting section')
	call delete_section('bound')
	call check_parameter_values('after deleting section')

!------------------------------------------------------
! end of routine
!------------------------------------------------------

	return
   92	continue
	write(6,*) 'Error in name list read : ',iweich
	stop 'error stop : rdbnds'
   93	continue
	write(6,*) 'Variable not allowed in vector context : ',name
	stop 'error stop : rdbnds'
   94   continue
        write(6,*) 'Boundary conditions out of order : ',ibc
        stop 'error stop : rdbnds'
	end

c********************************************************************

	subroutine ckbnds

c checks boundary information read from STR

	use mod_bnd
	use mod_bound_geom
	use shympi

	implicit none

	logical bstop
	integer istop
	integer i,k,ibc
	integer iqual,ibtyp,kranf,krend,ktilt,knode,kref
	integer levmax,levmin
	integer intpol
	integer kmanf,kmend,krtot,kmtot
	real period
	real ztilt
	integer ipint
	character*80 file

	bstop = .false.

	do i=1,nbc

	 ibc = i

         call get_bnd_ipar(ibc,'iqual',iqual)
	 if(iqual.ge.i) then
	   write(6,'(a,i2,a)') 'Section BOUND ',i,' :'
	   write(6,*) '   Value for iqual must be'
	   write(6,*) '   smaller than actual boundary number'
	   write(6,*) '   iqual = ',iqual
	   bstop=.true.
	 else if(iqual.lt.0) then
	   write(6,'(a,i2,a)') 'Section BOUND ',i,' :'
	   write(6,*) '   Value for iqual must be'
	   write(6,*) '   greater than 0'
	   write(6,*) '   iqual = ',iqual
	   bstop=.true.
	 else if( iqual .gt. 0 ) then
           call get_bnd_ipar(iqual,'ibtyp',ibtyp)
           call set_bnd_ipar(ibc,'ibtyp',ibtyp)
           call get_bnd_ipar(iqual,'levmax',levmax)
           call set_bnd_ipar(ibc,'levmax',levmax)
           call get_bnd_ipar(iqual,'levmin',levmin)
           call set_bnd_ipar(ibc,'levmin',levmin)
	 end if

         call get_bnd_ipar(ibc,'ibtyp',ibtyp)
         if(ibtyp.ge.0.and.ibtyp.le.3) then !$$ibtyp3
	 else if(ibtyp.eq.4) then		 !$$ibtyp11
	 else if(ibtyp.eq.11) then		 !$$ibtyp11
	 else if(ibtyp.ge.30.and.ibtyp.le.33) then	 !$$ibtyp11
	 else if(ibtyp.ge.50.and.ibtyp.le.53) then	 !$$ibtyp11
	 else if(ibtyp.eq.70) then	 !$$roger
	 else
	   write(6,'(a,i2,a)') 'Section BOUND ',i,' :'
	   write(6,*) '   Value not allowed for ibtyp'
	   write(6,*) '   ibtyp = ',ibtyp
	   bstop=.true.
	 end if

	 if( ibtyp .gt. 0 ) then
	   call get_boundary_file(ibc,'zeta',file)
           call get_bnd_par(ibc,'period',period)
	   if( period .le. 0. .and. file .eq. ' ' ) then
		write(6,'(a,i2,a)') 'section BOUND ',i,' :'
		write(6,*) '   Period must be > 0'
		write(6,*) '   period = ',period
		bstop=.true.
	   end if
	 end if

         call get_bnd_ipar(ibc,'ktilt',ktilt)
	 if( ktilt .gt. 0 ) then
	   knode=ipint(ktilt)		!$$EXTINW
	   if(knode.le.0) then
		write(6,'(a,i2,a)') 'Section BOUND ',i,' :'
		write(6,*) '   tilt node not found ',ktilt
		bstop=.true.
	   end if
           call set_bnd_ipar(ibc,'ktilt',knode)
	 end if

         call get_bnd_ipar(ibc,'kref',kref)
	 if( kref .gt. 0 ) then
	   knode=ipint(kref)		!$$EXTINW
	   if(knode.le.0) then
		write(6,'(a,i2,a)') 'Section BOUND ',i,' :'
		write(6,*) '   reference node not found ',kref
		bstop=.true.
	   end if
           call set_bnd_ipar(ibc,'kref',knode)
	 end if

         call get_bnd_ipar(ibc,'intpol',intpol)
	 if(intpol.lt.0.or.intpol.gt.4) then
		write(6,'(a,i2,a)') 'Section BOUND ',i,' :'
		write(6,*) '   erroneous value for INTPOL : ',intpol
		bstop=.true.
	 end if

         call get_bnd_ipar(ibc,'kranf',kranf)
         call get_bnd_ipar(ibc,'krend',krend)
	 if( kranf .gt. krend ) then	!$$kranf
	   write(6,'(a,i2,a)') 'section BOUND ',i,' :'
	   write(6,*) '   No nodes given for boundary'
	   bstop=.true.
	 end if

	 kmanf = 0
	 kmend = 0
	 krtot = krend-kranf+1

	 istop = 0
	 do k=kranf,krend
	   if( k == 0 ) cycle
	   knode=ipint(irv(k))		!$$EXTINW
	   if(knode.le.0) then
             if( .not. bmpi ) then
               write(6,'(a,i2,a)') 'Section BOUND ',i,' :'
               write(6,*) '   boundary node not found ',irv(k)
	     end if
           else if( .not. shympi_is_inner_node(knode) ) then
             !knode = 0		!we keep ghost node
	   end if
	   if( knode /= 0 ) then
	     if( kmanf == 0 ) kmanf = k
	     kmend = k
	   else
	     istop = istop + 1
	   end if
	   irv(k)=knode
	 end do

	 kmtot = kmend-kmanf+1
         call set_bnd_ipar(ibc,'kmanf',kmanf)
         call set_bnd_ipar(ibc,'kmend',kmend)
!         write(6,'(a,10i5)') 'boundary: ',my_id,ibc,istop
!     +			,kranf,krend,kmanf,kmend

         if( istop > 0 ) then
           write(6,'(a,10i5)') 'boundary: ',my_id,ibc,istop
     +			,kranf,krend,kmanf,kmend
           if( shympi_is_parallel() ) then
             if( istop == krtot ) then
               write(6,*) 'boundary completely in one domain... ok'
               call set_bnd_ipar(ibc,'ibtyp',0)
             else if( istop == krtot-kmtot ) then
               write(6,*) 'boundary in more than one domain... ok'
             else
	       stop 'error stop ckbnds: internal error'
             end if
           else
	     stop 'error stop ckbnds: no MPI and missing nodes'
           end if
         end if

	end do

	!call shympi_exit(0)
	if( bstop ) stop 'error stop: ckbnds'

	end

c********************************************************************

	subroutine prbnds

	use mod_bnd
	use mod_bound_geom

	implicit none

	integer i,ibc
	integer ibtyp,kranf,krend
	integer intpol
	integer ipext
	logical, parameter :: ball = .false.	!write all file names

	if( nbc .le. 0 ) return

	write(6,*)
	write(6,*) '====== info on open boundaries ========='
	write(6,*)

	do ibc=1,nbc
          call get_bnd_ipar(ibc,'kranf',kranf)
          call get_bnd_ipar(ibc,'krend',krend)

	  write(6,'(a,i9)') 'boundary: ',ibc
	  if( kranf > 0 .and. krend > 0 ) then
	    write(6,'(a,i9)') ' boundary nodes: ',krend-kranf+1
	    write(6,'(8i9)') (ipext(irv(i)),i=kranf,krend)
	  end if

	  call check_bnd_par_entries(ibc,ball)
	  call check_bnd_file_entries(ibc,ball)
	end do

	write(6,*)
	write(6,*) '==== end info on open boundaries ======='
	write(6,*)

	return
	end

c********************************************************************

	subroutine tsbnds

	use mod_bnd
	use mod_bound_geom

	implicit none

	integer j,ibc
	logical, parameter :: ball = .true.	!write all file names

	write(6,*) '/bnd/'
	write(6,*) nbc
	do ibc=1,nbc
	  call check_bnd_par_entries(ibc,ball)
	  call check_bnd_file_entries(ibc,ball)
	end do

	write(6,*) '/irv/'
	write(6,*) nrb
	do j=1,nrb
	  write(6,*) irv(j)
	end do

	end

c********************************************************************
c********************************************************************
c********************************************************************
c      utility routines
c********************************************************************
c********************************************************************
c********************************************************************

	function zvbnds(ibc)

c returns value of open boundary

	implicit none

	real zvbnds
	integer ibc

	real zval

        call chkibc(ibc,'zvbnds:')

        call get_bnd_par(ibc,'zval',zval)
	zvbnds = zval

	end

c********************************************************************

	subroutine setbnds(ibc,zval)

c sets value of open boundary

	implicit none

	integer ibc
	real zval

        call chkibc(ibc,'setbnds:')

        call set_bnd_par(ibc,'zval',zval)

	end

c********************************************************************

	subroutine stybnd(ibc,ibtyp)

c sets type of open boundary

	implicit none

	integer ibc
	integer ibtyp

        call chkibc(ibc,'stybnd:')

        call set_bnd_ipar(ibc,'ibtyp',ibtyp)

	end

c********************************************************************

	function itybnd(ibc)

c returns type of open boundary

	implicit none

	integer itybnd
	integer ibc

	integer ibtyp

        call chkibc(ibc,'itybnd:')

        call get_bnd_ipar(ibc,'ibtyp',ibtyp)
	itybnd = ibtyp

	end

c********************************************************************

	function nbnds()

c returns total number of open boundaries

	use mod_bnd

	implicit none

	integer nbnds


	nbnds = nbc

	end

c********************************************************************

	function nkbnd()

c returns total number of open boundary nodes

	use mod_bnd

	implicit none

	integer nkbnd


	nkbnd = nrb

	end

c********************************************************************

	function nkbnds(ibc)

c returns total number of nodes of boundary ibc

	implicit none

	integer nkbnds
	integer ibc

	integer kranf,krend

        call chkibc(ibc,'nkbnds:')

        call get_bnd_ipar(ibc,'kranf',kranf)
        call get_bnd_ipar(ibc,'krend',krend)

	nkbnds = krend - kranf + 1
	if( kranf == 0 .or. krend == 0 ) nkbnds = 0

	end

c********************************************************************

	function kbnd(i)

c returns i th node of all boundary nodes

	use mod_bnd
	use mod_bound_geom

	implicit none

	integer kbnd
	integer i

        if( i .gt. nrb ) then
            write(6,*) 'i, nrb, nbc : ',i,nrb,nbc
            stop 'error stop kbnd: i out of bounds'
        end if

	kbnd = irv(i)

	end

c********************************************************************

	subroutine kanfend(ibc,kranf,krend)

c returns index of first and last boundary node of boundary ibc

	implicit none

	integer ibc
        integer kranf,krend

        call chkibc(ibc,'kanfend:')

        call get_bnd_ipar(ibc,'kranf',kranf)
        call get_bnd_ipar(ibc,'krend',krend)

	end

c********************************************************************

	subroutine kmanfend(ibc,kmanf,kmend)

c returns index of first and last boundary node of boundary ibc
c mpi version

	implicit none

	integer ibc
        integer kmanf,kmend

        call chkibc(ibc,'kmanfend:')

        call get_bnd_ipar(ibc,'kmanf',kmanf)
        call get_bnd_ipar(ibc,'kmend',kmend)

	end

c********************************************************************

	function kbndind(ibc,i)

c returns global index of i th node of boundary ibc

	implicit none

	integer kbndind
	integer ibc
	integer i

        integer kranf,krend,n

        call chkibc(ibc,'kanfend:')

        call get_bnd_ipar(ibc,'kranf',kranf)
        call get_bnd_ipar(ibc,'krend',krend)

	n = krend - kranf + 1

	if( i .lt. 1 .or. i .gt. n ) then
	    write(6,*) 'i, imax, ibc : ',i,n,ibc
	    stop 'error stop kbndind: i out of bounds'
	end if

	kbndind = kranf + i - 1

	end

c********************************************************************

	function kbnds(ibc,i)

c returns i th node of boundary ibc

	use mod_bound_geom

	implicit none

	integer kbnds
	integer ibc
	integer i

	integer kranf,krend,n

        call chkibc(ibc,'kbnds:')

        call get_bnd_ipar(ibc,'kranf',kranf)
        call get_bnd_ipar(ibc,'krend',krend)

	n = krend - kranf + 1

	if( i .lt. 1 .or. i .gt. n ) then
	    write(6,*) 'i, imax, ibc : ',i,n,ibc
	    stop 'error stop kbnds: i out of bounds'
	end if

	kbnds = irv(kranf+i-1)

	end

c********************************************************************

	subroutine irbnds(ibc,ndim,idim,nodes)

c returns nodes of boundary ibc (maximum ndim)

	use mod_bound_geom

	implicit none

	integer ibc			!number of open boundary  (in)
	integer ndim			!dimension of nodes()     (in)
	integer idim			!total number of nodes    (out)
	integer nodes(ndim)		!boundary nodes           (out)

	integer i,imaxi
	integer kranf,krend

        call chkibc(ibc,'irbnds:')

        call get_bnd_ipar(ibc,'kranf',kranf)
        call get_bnd_ipar(ibc,'krend',krend)

	idim = krend - kranf + 1

	imaxi = min(ndim,idim)

	do i=1,imaxi
	  nodes(i) = irv(i+kranf-1)
	end do

	end

c********************************************************************

        subroutine ksinget(ibc,ampli,period,phase,zref)
        call chkibc(ibc,'ksinget:')
        call get_bnd_par(ibc,'ampli',ampli)
        call get_bnd_par(ibc,'period',period)
        call get_bnd_par(ibc,'phase',phase)
        call get_bnd_par(ibc,'zref',zref)
        end

c********************************************************************

	subroutine setbnd(ibc,value,barray)

c sets boundary ibc to value in barray (apparently not used)

	use mod_bound_geom

	implicit none

	integer ibc
	real value
	real barray(1)

	integer kranf,krend,k,kn

        call chkibc(ibc,'setbnd:')

        call get_bnd_ipar(ibc,'kranf',kranf)
        call get_bnd_ipar(ibc,'krend',krend)

	do k=kranf,krend
	  kn = irv(k)
	  barray(kn) = value
	end do

	end

c********************************************************************

        subroutine chkibc(ibc,errtext)
 
c checks if ibc is in bounds
 
	use mod_bnd

        implicit none
 
        integer ibc
	character*(*) errtext
 
        if( ibc .lt. 1 .or. ibc .gt. nbc ) then
	    write(6,*) errtext
            write(6,*) 'ibc, nbc : ',ibc,nbc
            stop 'error stop bnd: ibc out of bounds'
        end if

	end

c********************************************************************

	subroutine is_closed(ibc,nbc,icl)

	implicit none

	integer ibc,nbc,icl

	integer i
        integer nbnds,itybnd

        icl=0
        nbc = nbnds()

	if( ibc .gt. 0 ) then
	  if( itybnd(ibc) .lt. 0 ) icl = icl + 1
	else
          do i=1,nbc
            if( itybnd(i) .lt. 0 ) icl = icl + 1
          end do
	end if

	end

c********************************************************************

	subroutine get_oscil(ibc,dtime,zvalue)

	implicit none

	integer ibc
	double precision dtime
	real zvalue

        real pi
        parameter( pi=3.141592653 )

	real ampli,period,phase,zref

	call ksinget(ibc,ampli,period,phase,zref)

        zvalue = zref+ampli*sin(2.*pi*(dtime+phase)/period)

	end

c********************************************************************
c********************************************************************
c********************************************************************
c
c routines dealing with parameter values
c
c********************************************************************
c********************************************************************
c********************************************************************

	subroutine init_bnd_par(ibc)

c initializes boundary ibc

	use mod_bnd

	implicit none

	integer ibc

	bnd(:,ibc) = 0.

	end

c********************************************************************

        subroutine setget_bnd_par(ibc,ientry,name,value,bset)

c sets/gets value at entry ientry

	use mod_bnd

        implicit none

        integer ibc
        integer ientry
	character*(*) name
        real value
        logical bset

        call chkibc(ibc,'setget_bnd_par:')

        if( ientry < 1 .or. ientry > nbvdim ) then
          write(6,*) 'no such entry: ',ientry,ibc,name
	  stop 'error stop setget_bnd_par: no such entry'
        end if

        if( bset ) then
          bnd(ientry,ibc) = value
        else
          value = bnd(ientry,ibc)
        end if

        end

c********************************************************************
 
        subroutine set_bnd_par(ibc,name,value)
        character*(*) name
        id = iget_bnd_par_id(name,.true.)
        call setget_bnd_par(ibc,id,name,value,.true.)
        end

        subroutine get_bnd_par(ibc,name,value)
        character*(*) name
        id = iget_bnd_par_id(name,.true.)
        call setget_bnd_par(ibc,id,name,value,.false.)
        end

        subroutine set_bnd_ipar(ibc,name,ivalue)
        character*(*) name
        id = iget_bnd_par_id(name,.true.)
	value = ivalue
        call setget_bnd_par(ibc,id,name,value,.true.)
        end

        subroutine get_bnd_ipar(ibc,name,ivalue)
        character*(*) name
        id = iget_bnd_par_id(name,.true.)
        call setget_bnd_par(ibc,id,name,value,.false.)
	ivalue = nint(value)
        end

c********************************************************************

	subroutine copy_bnd_par(ibc)

c copies parameter values

	use mod_bnd

	implicit none

	integer ibc

	integer id
	real value
	character*6 name

	real getpar

	do id=1,nbvdim
	  call get_bnd_par_name(id,name)
	  value = getpar(name)
          call set_bnd_par(ibc,name,value)
	end do

	end

c********************************************************************

	subroutine get_bnd_npar(nbnd)

	use mod_bnd

        implicit none

	integer nbnd

	nbnd = nbvdim

	end

c********************************************************************

        function exists_bnd_par(name)

c tests if parameter name exists

        implicit none

        logical exists_bnd_par
        character*(*) name

	integer iget_bnd_par_id

	exists_bnd_par = iget_bnd_par_id(name,.false.) > 0

	end

c********************************************************************

        function iget_bnd_par_id(name,berror)

c gets id given a parameter name

	use mod_bnd

        implicit none

        integer iget_bnd_par_id
        character*(*) name
	logical berror		!raises error if name not existing

        integer id
        character*6 bname

        do id=1,nbvdim
	  call get_bnd_par_name(id,bname)
          if( bname .eq. name ) then
            iget_bnd_par_id = id
            return
          end if
        end do

	if( berror ) then
          write(6,*) 'unknown parameter name for boundary: ',name
          stop 'error stop iget_bnd_par_id: name'
	end if

	iget_bnd_par_id = 0

        end

c********************************************************************

	subroutine get_bnd_par_name(id,name)

c gets parameter name given id

	use mod_bnd

	implicit none

	integer id
        character*(*) name

	if( id < 1 .or. id > nbvdim ) then
	  name = ' '
	else
	  name = bnd_par_names(id)
	end if

	end

c********************************************************************

	subroutine check_bnd_par_entries(ibc,ball)

c writes parameter values for given boundary

	use mod_bnd

	implicit none

	integer ibc
	logical ball	!write all parameter names

        integer i
	real value
	character*6 name

	!write(6,*) 'check_bnd_par_entries: ',ibc
        do i=1,nbvdim
	  call get_bnd_par_name(i,name)
          call setget_bnd_par(ibc,i,' ',value,.false.)
	  if( ball .or. value /= 0. ) then
	    write(6,*) name,value
	  end if
	end do

        end

c********************************************************************
c********************************************************************
c********************************************************************
c
c routines dealing with file names
c
c********************************************************************
c********************************************************************
c********************************************************************

	subroutine init_bnd_file(ibc)

c initializes boundary ibc

	use mod_bnd

	implicit none

	integer ibc

	bnd_file(:,ibc) = ' '

	end

c********************************************************************

        subroutine setget_bnd_file(ibc,ientry,name,file,bset)

c sets/gets file at entry ientry

	use mod_bnd

        implicit none

        integer ibc
        integer ientry
	character*(*) name
	character*(*) file
        logical bset

        call chkibc(ibc,'setget_bnd_file:')

        if( ientry < 1 .or. ientry > nbfdim ) then
          write(6,*) 'no such entry: ',ientry,ibc,name
	  stop 'error stop setget_bnd_file: no such entry'
        end if

        if( bset ) then
          bnd_file(ientry,ibc) = file
        else
          file = bnd_file(ientry,ibc)
        end if

        end

c********************************************************************
 
        subroutine set_bnd_file(ibc,name,file)
        character*(*) name,file
        id = iget_bnd_file_id(name,.true.)
        call setget_bnd_file(ibc,id,name,file,.true.)
        end

        subroutine get_bnd_file(ibc,name,file)
        character*(*) name,file
        id = iget_bnd_file_id(name,.true.)
        call setget_bnd_file(ibc,id,name,file,.false.)
        end

c********************************************************************

	subroutine copy_bnd_file(ibc)

c copies file names

	use mod_bnd

	implicit none

	integer ibc

	integer id
	character*80 file
	character*6 name

	do id=1,nbfdim
	  call get_bnd_file_name(id,name)
	  call getfnm(name,file)
          call set_bnd_file(ibc,name,file)
	end do

	end

c********************************************************************

	subroutine get_bnd_nfile(nbnd)

	use mod_bnd

        implicit none

	integer nbnd

	nbnd = nbfdim

	end

c********************************************************************

        function exists_bnd_file(name)

c tests if file name exists

        implicit none

        logical exists_bnd_file
        character*(*) name

	integer iget_bnd_file_id

	exists_bnd_file = iget_bnd_file_id(name,.false.) > 0

	end

c********************************************************************

        function iget_bnd_file_id(name,berror)

c gets id given a file name

	use mod_bnd

        implicit none

        integer iget_bnd_file_id
        character*(*) name
	logical berror		!raises error if name not existing

        integer id
        character*6 bname

        do id=1,nbfdim
	  call get_bnd_file_name(id,bname)
          if( bname .eq. name ) then
            iget_bnd_file_id = id
            return
          end if
        end do

	if( berror ) then
          write(6,*) 'unknown file name for boundary: ',name
          stop 'error stop iget_bnd_file_id: name'
	end if

	iget_bnd_file_id = 0

        end

c********************************************************************

	subroutine get_bnd_file_name(id,name)

c gets file name given id

	use mod_bnd

	implicit none

	integer id
        character*(*) name

	if( id < 1 .or. id > nbfdim ) then
	  name = ' '
	else
	  name = bnd_file_names(id)
	end if

	end

c********************************************************************

	subroutine check_bnd_file_entries(ibc,ball)

c writes file names for given boundary

	use mod_bnd

	implicit none

	integer ibc
	logical ball	!write all file names

        integer i
	character*6 name
	character*80 file

	!write(6,*) 'check_bnd_file_entries: ',ibc
        do i=1,nbfdim
	  call get_bnd_file_name(i,name)
          call setget_bnd_file(ibc,i,' ',file,.false.)
	  if( ball .or. file /= ' ' ) then
	    write(6,*) name,'  ',trim(file)
	  end if
	end do

        end

c********************************************************************
c********************************************************************
c********************************************************************

	subroutine get_boundary_file(ibc,what,file)

	implicit none

	integer ibc
	character*(*) what
	character*(*) file

	character*6 name

	logical exists_bnd_file

        if( what .eq. 'zeta' ) then
	  name = 'boundn'
        else if( what .eq. 'conz' ) then
	  name = 'conzn'
        else if( what .eq. 'temp' ) then
	  name = 'tempn'
        else if( what .eq. 'salt' ) then
	  name = 'saltn'
        else if( what .eq. 'vel' ) then
	  name = 'vel3dn'
        else if( what .eq. 'sedt' ) then
	  name = 'sed2dn'
        else if( what .eq. 'lagvebio' ) then
	  name = 'bio2dn'
        else if( what .eq. 'toxi' ) then
	  name = 'tox3dn'
        else if( what .eq. 'bfm1' ) then
	  name = 'bfm1bc'
        else if( what .eq. 'bfm2' ) then
	  name = 'bfm2bc'
        else if( what .eq. 'bfm3' ) then
	  name = 'bfm3bc'
        else if( what .eq. 'bfm' ) then
	  name = 'bfmbcn'
        else if( what .eq. 'mercury' ) then
	  name = 'mercn'
        else
          if( exists_bnd_file(what) ) then	!use name given in what
	    name = what
	  else
            write(6,*) 'keyword not recognized: ',what
            write(6,*) 'boundary: ',ibc
            stop 'error stop get_boundary_file'
	  end if
        end if

        call get_bnd_file(ibc,name,file)

	end

c********************************************************************
c********************************************************************
c********************************************************************

