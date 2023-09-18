!
! $Id: subbnd.f,v 1.35 2009-05-21 09:24:00 georg Exp $
!
! bnd administration routines
!
! contents :
!
! subroutine inbnds
! subroutine rdbnds(ibc)
! subroutine ckbnds
! subroutine prbnds
! subroutine tsbnds
!
! function zvbnds(ibc)			 returns value of open boundary
! subroutine stybnd(ibc,ibtyp)		 sets type of open boundary
! subroutine infobnd(ibc,ibtype,...)	 returns important info on boundary ibc
! function itybnd(ibc)			 returns type of open boundary
! function nbnds			 returns total number of open b.
! function nkbnds(ibc)			 returns total number of nodes of ibc
! function kbnd(i)			 returns i th node of all b. nodes
! subroutine kanfend(ibc,kranf,krend)    returns index of first and last bnode
! function kbnds(ibc,i)			 returns i th node of boundary ibc
! subroutine irbnds(ibc,ndim,idim,nodes) returns nodes of boundary ibc
! subroutine setget_bnd_par(ibc,ientry,value,bset) sets/gets value at ientry
! subroutine setbnd(ibc,value,barray)	 sets boundary ibc to value in barray
!
! subroutine setbc(value,array,flag)	 sets all open boundaries to value
! subroutine chkibc(ibc,errtext)	 checks if ibc is in bounds
!
! function get_discharge(ibc)		returns discharge through boundary ibc
! subroutine adjust_mass_flux		adjusts mass flux for dry nodes
!
! notes :
!
! for scalars: -999. uses ambient value
!
! nbndim is the maximum dimension used for boundary information 
! nbvdim is the actual filling of the variables
!
! revision log :
!
! 01.08.1997	ggu	$$1stnode - first boundary node was not registered
! 23.09.1997	ggu	introduced conzn -> name of file for conz values
! 30.09.1997	ggu	write error in prbnds
! 03.12.1997	ggu	introduced tempn,saltn (see conzn)
! 27.03.1998	ggu	new utility routines
! 25.05.1998	ggu	documentation (DOCS)
! 20.06.1998	ggu	documentation for momentum input
! 13.07.1998	ggu	implemented ibtyp = 4
! 23.07.1998	ggu	documentation
! 24.08.1998    ggu     BC for concentration is bnd(20,..)
! 24.08.1998    ggu     BC for maximum input level is bnd(12,..) -> levmax
! 27.08.1998    ggu     accessor function levbnd for levmax
! 22.01.1999    ggu     new subroutine setbnd
! 08.07.1999    ggu     bug with iqual -> ibtyp not respected
! 20.01.2000    ggu     call to rdbnds without dimension -> use getdim
! 07.04.2000    ggu     new subroutine setbc
! 07.05.2001    ggu     introduced new variable zfact
! 25.09.2001    ggu     introduced bio2dn
! 07.08.2003    ggu     check for nrb in rdbnds
! 15.10.2004    ggu     new boundary types and sedin
! 02.03.2005    ggu     new nbdim for 3D boundary values
! 02.03.2005    ggu     some new helper functions
! 07.11.2005    ccf     introduced sed2dn
! 16.02.2006    ggu     introduced tox3dn
! 07.04.2008    acc     introduced bfm1bc bfm2bc bfm3bc OB condition for ERSEM
! 17.04.2008    ggu     deleted infobnd(), levbnd()
! 28.04.2008    ggu     call to nrdpar in double precision
! 29.04.2008    ggu&aac new boundary arrays for ERSEM
! 30.05.2008	ggu	eliminated numbers for parameters
! 03.06.2008	ggu	new parameters levmin, kref
! 06.06.2008	ggu	completely restructured
! 27.03.2009    ggu     new routine adjust_mass_flux() for dry nodes
! 02.04.2009	ggu	intpol default is 0, some unused routines deleted
! 20.04.2009	ggu	new variable ztilt, ndim substituted with nbvdim
! 23.02.2011    ggu     new parameters tramp and levflx implemented
! 21.06.2012    ggu&aar new file names for mud module
! 29.11.2013    ggu	allow for non continous boundary numbering
! 28.03.2014    ggu	new parameter lgrpps
! 16.06.2014    ggu	new include file bnd.h (with nbndim)
! 25.06.2014    ggu	new routine exists_bnd_name()
! 29.10.2014    ccf	include vel3dn boundary file
! 03.11.2014    ggu	nbdim deleted
! 23.06.2015    ggu	setbc() deleted, nrz,nrq eliminated
! 14.01.2016    ggu	check of boundaries considers mpi subdomains
!
!************************************************************************
!------------------------------------------------------------------------
        module bnd_admin
!------------------------------------------------------------------------
        contains
!------------------------------------------------------------------------

	subroutine inbnds

! initializes boundary parameters

	use mod_bnd

	implicit none


	nbc = 0
	nrb = 0

	end

!************************************************************************

	subroutine rdbnds(ibc)

! reads boundary info from STR file

	use mod_bnd
	use bnd_geom
        use para
        use nls

	implicit none

	integer ibc

	include 'param.h'

	include 'bound_names.h'
!AR mud

	character*80 name,text
	double precision dvalue
	double precision value
	integer n
	integer nbcdi
	integer i,kranf,krend,kref
	integer iweich,id,nbnd

	integer icall
	save icall
	data icall / 0 /

	nbcdi = nbc_dim

	if( icall .eq. 0 ) then		!initialize bnd array
	  do i=1,nbcdi

	    !call bnd_init(i)

	    boundn(i) = ' '
	    conzn(i) = ' '
	    tempn(i) = ' '
	    saltn(i) = ' '
	    vel3dn(i) = ' '

	    bio2dn(i) = ' '
	    sed2dn(i) = ' '
	    mud2dn(i) = ' '
	    dmf2dn(i) = ' '
	    lam2dn(i) = ' '
	    tox3dn(i) = ' '

	    bfm1bc(i) = ' '
	    bfm2bc(i) = ' '
	    bfm3bc(i) = ' '

	  end do
	  icall = 1
	end if

	n = max(nbc,ibc)	!allow for non contiguous numbering
	if(n.gt.nbcdi) goto 77
	call mod_bnd_adjust(n)
	do i=nbc+1,n
	  call bnd_init(i)
	end do
	nbc = n

	call sctpar('bound')
	call sctfnm('bound')

! DOCS	START	S_bound
!
! These parameters determine the open boundary nodes and the type
! of the boundary: level or flux boundary. At the first the water levels
! are imposed, on the second the fluxes are prescribed.
!
! There may be multiple sections |bound| in one parameter input file,
! describing all open boundary conditions necessary. Every section
! must therefore be supplied with a boundary number. The numbering
! of the open boundaries must
! be increasing. The number of the boundary must be specified
! directly after the keyword |bound|, such as |bound1| or |bound 1|.
!
! |kbound|	Array containing the node numbers that are part of the
!		open boundary. The node numbers must form one contiguous
!		line with the domain (elements) to the left. This
!		corresponds to an anti-clockwise sense. At least
!		two nodes must be given.

	call addpar('kbound',0.d0)

! |ibtyp|	Type of open boundary. 
!		\begin{description}
!		\item[0] No boundary values specified
!		\item[1] Level boundary. At this open boundary
!			 the water level is imposed and the prescribed
!			 values are interpreted as water levels in meters.
!			 If no value for |ibtyp| is specified this
!			 is the default.
!		\item[2] Flux boundary. Here the discharge in \dischargeunit
!			 has to be prescribed.
!		\item[3] Internal flux boundary. As with |ibtyp = 2| a
!			 discharge has to be imposed, but the node where
!			 discharge is imposed can be an internal node
!			 and need not be on the outer boundary of
!			 the domain. For every node in |kbound| the
!			 volume rate specified will be added to the
!			 existing water volume. This behavior is different
!			 from the |ibtyp = 2| where the whole boundary
!			 received the discharge specified.
!		\item[4] Momentum input. The node or nodes may be internal.
!			 This feature can be used to describe local 
!			 acceleration of the water column. 
!			 The unit is force / density [\maccelunit].
!			 In other words it is the rate of volume 
!			 [\dischargeunit] times the velocity [m/s] 
!			 to which the water is accelerated.
!		\end{description}
! |iqual|	If the boundary conditions for this open boundary
!		are equal to the ones of boundary |i|, then
!		setting |iqual = i| copies all the values of
!		boundary |i| to the actual boundary. Note that the
!		value of |iqual| must be smaller than the number
!		of the actual boundary, i.e., boundary |i| must have
!		been defined before.

	call addpar('ibtyp',1.d0)
	call addpar('iqual',0.d0)

! The next parameters give a possibility to specify the file name
! of the various input files that are to be read by the model.
! Values for the boundary condition can be given at any time step.
! The model interpolates in between given time steps if needed. The
! grade of interpolation can be given by |intpol|.
!
! All files are in ASCII and share a common format.
! The file must contain two columns, the first giving the
! time of simulation in seconds that refers to the value
! given in the second column. The value in the second
! column must be in the unit of the variable that is given.
! The time values must be in increasing order.
! There must be values for the whole simulation,
! i.e., the time value of the first line must be smaller
! or equal than the start of the simulation, and the time
! value of the last line must be greater or equal than the
! end of the simulation.

! |boundn|	File name that contains values for the boundary condition.
!		The value of the variable given in the second column
!		must be in the unit determined by |ibtyp|, i.e.,
!		in meters for a level boundary, in \dischargeunit for
!		a flux boundary and in \maccelunit for a momentum
!		input.
! |zfact|	Factor with which the values from |boundn|
!		are multiplied to form the final value of the
!		boundary condition. E.g., this value can be used to
!		set up a quick sensitivity run by multiplying
!		all discharges by a factor without generating
!		a new file. (Default 1)

	call addfnm('boundn',' ')
	call addpar('zfact',1.d0)

! |levmin, levmax|	A point discharge normally distributes its discharge
!		over the whole water column. If it is important that in
!		a 3D simulation the water mass discharge is concentrated 
!		only in some levels, the parameters |levmin| and |levmax|
!		can be used. They indicate the lowest and deepest level over
!		which the discharge is distributed. Default values are 0, which
!		indicate that the discharge is distributed over the
!		whole water column. Setting only |levmax| distributes from
!		the surface to this level, and setting only |levmin|
!		distributes from the bottom to this level.

	call addpar('levmin',0.d0)
	call addpar('levmax',0.d0)

! |conzn, tempn, saltn|	File name that contains values for the respective
!			boundary condition, i.e., for concentration,
!			temperature and salinity. The format is the same
!			as for file |boundn|. The unit of the values
!			given in the second column must the ones of the
!			variable, i.e., arbitrary unit for concentration,
!			centigrade for temperature and psu (per mille)
!			for salinity.

	call addfnm('conzn',' ')
	call addfnm('tempn',' ')
	call addfnm('saltn',' ')

! |vel3dn|	File name that contains current velocity values for the 
!		boundary condition.  The format is the same as for file 
!		|tempn| but it has two variables: 
!	 	current velocity in x and current velocity in y.
!		Velocity can be nudged or imposed depending on the value 
!		of |tnudge| (mandatory). The unit is [m/s].

	call addfnm('vel3dn',' ')

! |tnudge|	Relaxation time for nudging of boundary velocity.
!		For |tnudge| = 0. velocities are imposed, for
!		|tnudge| > 0. velocities are nudged. The
!		default is -1 which means do nothing. Unit is [s].
!		(Default -1)

        call addpar('tnudge',-1.d0)

! The next variables specify the name of the boundary value file
! for different modules. Please refer to the documentation of the
! single modules for the units of the variables.

! |bio2dn|	File name that contains values for the ecological
!		module (EUTRO-WASP).
! |sed2dn|	File name that contains values for the sediment
!		transport module.
!		The unit of the values given
!		in the second and following columns (equal to the 
!		number of defined grainsize in parameter |sedgrs|).

! |mud2dn|	File name that contains values for the fluid mud
!		module.
! |lam2dn|	File name that contains values for the fluid mud
!		module (boundary condition for the structural parameter, 
!		to be implemented).
! |dmf2dn|	File name that contains values for the fluid mud
!		module (boundary conditions for the advection of flocsizes,
!		to be implemented).
! |tox3dn|	File name that contains values for the toxicological
!		module.

	call addfnm('bio2dn',' ')       !HACK
	call addfnm('sed2dn',' ')       !HACK
	call addfnm('mud2dn',' ')       !HACK
	call addfnm('lam2dn',' ')       !HACK
	call addfnm('dmf2dn',' ')       !HACK
	call addfnm('tox3dn',' ')       !HACK

!c File name for OB condition in ERSEM MODULE - undocumented

	call addfnm('bfm1bc',' ')
	call addfnm('bfm2bc',' ')
	call addfnm('bfm3bc',' ')

! |intpol|	Order of interpolation for the boundary values read
!		through files. Use for 1 for stepwise (no) interpolation,
!		2 for linear and 4 for cubic interpolation. 
!		The default is linear interpolation, except for
!		water level boundaries (|ibtyp=1|) where cubic
!		interpolation is used.

	call addpar('intpol',0.d0)

! The next parameters can be used to impose a sinusoidal water level
! (tide) or flux at the open boundary. These values are used if no
! boundary file |boundn| has been given. The values must be in the unit
! of the intended variable determined by |ibtyp|.

! |ampli|	Amplitude of the sinus function imposed. (Default 0)
! |period|	Period of the sinus function. (Default 43200, 12 hours)
! |phase|	Phase shift of the sinus function imposed. A positive value
!		of one quarter of the period reproduces a cosine
!		function. (Default 0)
! |zref|	Reference level of the sinus function imposed. If only
!		|zref| is specified (|ampli = 0|) a constant value
!		of |zref| is imposed on the open boundary.

	call addpar('ampli',0.d0)
	call addpar('period',43200.d0)
	call addpar('phase',0.d0)
	call addpar('zref',0.d0)

! With the next parameters a constant value can be imposed for the 
! variables of concentration, temperature and salinity. In this case
! no file with boundary values has to be supplied. The default for all
! values is 0, i.e., if no file with boundary values is supplied and
! no constant is set the value of 0 is imposed on the open boundary.
! A special value of -999 is also allwoed. In this case the value
! imposed is the ambient value of the parameter close to the boundary.

! |conz, temp, salt|	Constant boundary values for concentration,
!			temperature and salinity respectively. If these
!			values are set no boundary file has to be supplied.
!			(Default 0)

	call addpar('conz',0.d0)			!$$conz
	call addpar('temp',0.d0)			!$$baroc
	call addpar('salt',0.d0)			!$$baroc

! The next two values are used for constant momentum input. 
! This feature can be used to describe local acceleration of the
! water column. The values give the input of momentum 
! in x and y direction. The unit is force / density (\maccelunit).
! In other words it is the rate of volume (\dischargeunit) times
! the velocity (m/s) to which the water is accelerated.
!
! These values are used if 
! boundary condition |ibtyp = 4| has been chosen and
! no boundary input file has been given.
! If the momentum input is varying then it may be specified with
! the file |boundn|. In this case the file |boundn| must contain
! three columns, the first for the time, and the other two for
! the momentum input in $x,y$ direction.
!
! |umom, vmom|		Constant values for momentum input. (Default 0)

	call addpar('umom',0.d0)
	call addpar('vmom',0.d0)

!c undocumented

	call addpar('ktilt',0.d0)
	call addpar('ztilt',0.d0)
	call addpar('kref',0.d0)
	call addpar('igrad0',0.d0)	!use 0 gradient for scalars

	call addpar('tramp',0.d0)		!start smoothly for discharge
	call addpar('levflx',0.d0)	!use level-discharge relationship
	call addpar('nad',-1.d0)		!no advective terms for this boundary
        call addpar('lgrpps',0.d0)	!parts per second for lagrange
					!if negative parts per volume flux

! DOCS	END

! here add dummy variables

	call addpar('zval',0.d0)

!ccccccccccccccccccccccccccccccccccccccccccccccccc
! here we start the reading loop
!ccccccccccccccccccccccccccccccccccccccccccccccccc

	kranf = nrb + 1
	call set_bnd_ipar(ibc,'kranf',kranf)	!position of starting node
	call addpar('kranf',dble(kranf))

	iweich=1
	do while(iweich.ne.0)
	    !iweich=nrdpar('bound',name,value,text)
	    iweich=nrdpar('bound',name,dvalue,text)
	    value = dvalue
	    if( iweich .lt. 0 ) goto 92
	    if( iweich .eq. 2 .and. name .ne. 'kbound' ) goto 93
	    if( iweich .eq. 4 ) goto 93
	    if( name .eq. 'kbound' ) then	!$$1stnode
		nrb=nrb+1
		call mod_irv_init(nrb)
		irv(nrb)=nint(value)
	    else if( iweich .eq. 3 ) then	!file name
	        ! must be handeled later
	    else if( iweich .ne. 0 ) then
                call set_bnd_par(ibc,name,value)
	    end if
	end do

	krend = nrb
	call set_bnd_ipar(ibc,'krend',krend)	!position of end node
	call addpar('krend',dble(krend))

	call get_bnd_nbnd(nbnd)
	do id=1,nbnd
	  call get_bnd_name(id,name)
	  value = getpar(name)
          call set_bnd_par(ibc,name,value)
	end do

	call getfnm('boundn',boundn(ibc))
	call getfnm('conzn',conzn(ibc))
	call getfnm('tempn',tempn(ibc))
	call getfnm('saltn',saltn(ibc))
	call getfnm('vel3dn',vel3dn(ibc))

	call getfnm('bio2dn',bio2dn(ibc))
	call getfnm('sed2dn',sed2dn(ibc))
	call getfnm('mud2dn',mud2dn(ibc))
        call getfnm('dmf2dn',dmf2dn(ibc))
        call getfnm('lam2dn',lam2dn(ibc))
	call getfnm('tox3dn',tox3dn(ibc))

	call getfnm('bfm1bc',bfm1bc(ibc))
	call getfnm('bfm2bc',bfm2bc(ibc))
	call getfnm('bfm3bc',bfm3bc(ibc))

	!call check_bnd_entries(ibc)

	call check_parameter_values('before deleting section')
	call delete_section('bound')
	call check_parameter_values('after deleting section')

	return
   77	continue
	write(6,*) 'Dimension error for nbc_dim'
	write(6,*) 'nbc_dim :',nbcdi
	stop 'error stop : rdbnds'
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

!********************************************************************

	subroutine ckbnds

! checks boundary information read from STR

	use mod_bnd
	use bnd_geom
	use shympi
        use fem_util

	implicit none

	include 'param.h'

	logical bstop
        integer istop
	integer i,k,ibc
	integer iqual,ibtyp,kranf,krend,ktilt,knode,kref
	integer levmax,levmin
	integer intpol
	double precision period
	double precision ztilt
	character*80 file
        integer,dimension(:),allocatable :: temp
        integer,dimension(:),allocatable :: before,after
        integer icount,gid
	logical last
	integer first
	logical bmpielems

	bstop = .false.
        bounds%rankfirst = 0
	bmpielems = shympi_partition_on_elements()

        bounds%nob = nbc
        if (bounds%nob .gt. 0) then
          allocate(bounds%bstart(bounds%nob))
          allocate(bounds%bend(bounds%nob))
          allocate(bounds%tnob(bounds%nob))
          allocate(temp(size(irv)))
          allocate(before(size(irv)))
          allocate(after(size(irv)))
          !bounds%niob = 0
          bounds%bstart(1) = 1
          bounds%nnob = 0
	  before=0
	  after=0
	  bounds%nneighbors=0
        end if

	do i=1,nbc

	 ibc = i

         bounds%iobound=ibc

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

         call get_bnd_ipar(ibc,'kranf',kranf)
         call get_bnd_ipar(ibc,'krend',krend)
	 if( kranf .gt. krend ) then	!$$kranf
	   write(6,'(a,i2,a)') 'section BOUND ',i,' :'
	   write(6,*) '   No nodes given for boundary'
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

         istop = 0
         icount = 0
         bounds%tnob(ibc)=krend-kranf+1
	 first=0
	 last=.false.
	 do k=kranf,krend
	    if( k == 0 ) cycle
	    if(bmpielems) then
	      knode=ipint_mpi(irv(k))
	    else
	      knode=ipint(irv(k))		!$$EXTINW
	    end if
	    if((knode.le.0) .or. (bmpielems.and.knode.gt.domain%nodes%domainID)) then
              !write(6,'(a,i2,a)') 'Section BOUND ',i,' :'
              !write(6,*) '   boundary node not found ',irv(k)
              istop = istop + 1
	    end if
	    irv(k)=knode
            if((knode.gt.0).and.(bmpielems)) then
              if((ibc.eq.1).and.(k.eq.kranf)) then
                bounds%rankfirst = my_id
              end if
	      if(knode .le. domain%nodes%domainID) then
		last=.true.
                icount = icount + 1
                temp(icount+bounds%nnob)=knode
		if(first.gt.0) then
		  before(icount+bounds%nnob)=first
		  first=0
		end if
	      else
		irv(k)=0
		if(last) then
		  after(icount+bounds%nnob)=knode
		  last=.false.
		else
		  first=knode
		end if
	      end if
	    else
	      last=.false.
            end if
	 end do

         if( istop > 0 ) then
           if( shympi_is_parallel() ) then
             if( istop == krend-kranf+1 ) then
              !write(6,*) 'boundary completely in one domain... ok',my_id
               call set_bnd_ipar(ibc,'ibtyp',0)
               bounds%nob = bounds%nob - 1
             end if
           else
	     stop 'error stop ckbnds: no MPI and missing nodes'
           end if
         end if
         if(bounds%iobound .gt. 1) then
           bounds%bstart(bounds%iobound) = bounds%bend(bounds%iobound-1)+1
         end if
         bounds%nnob = bounds%nnob + icount
         bounds%bend(bounds%iobound)=bounds%bstart(bounds%iobound) + icount-1
	end do


	if(bmpielems) then
	  allocate(bounds%bneighbors(nbc))
	  bounds%bneighbors=0
	  icount=0
          allocate(bounds%niob(bounds%nnob))
          allocate(bounds%first(bounds%nnob))
          allocate(bounds%last(bounds%nnob))
          do i=1,bounds%nnob
           bounds%niob(i) = temp(i)
           bounds%first(i) = before(i)
           bounds%last(i) = after(i)
	   if(bounds%first(i).ne.0) bounds%nneighbors=bounds%nneighbors+1
	   if(bounds%last(i).ne.0) bounds%nneighbors=bounds%nneighbors+1
          end do
	  allocate(bounds%neighbors(bounds%nneighbors))
	  do i=1,bounds%nnob
	    if(bounds%first(i).ne.0) then
	     do k=1,nbc
	       if((i.le.(bounds%bend(k)-bounds%bstart(k)+1)).and.(i.ge.bounds%bstart(k))) then
	         bounds%bneighbors(k)=bounds%bneighbors(k)+1
	       end if
	     end do
	     icount=icount+1
	     bounds%neighbors(icount)=bounds%first(i)
	    end if
	    if(bounds%last(i).ne.0) then
	     do k=1,nbc
	       if((i.le.(bounds%bend(k)-bounds%bstart(k)+1)).and.(i.ge.bounds%bstart(k))) then
	         bounds%bneighbors(k)=bounds%bneighbors(k)+1
	       end if
	     end do
	     icount=icount+1
	     bounds%neighbors(icount)=bounds%last(i)
	    end if
	  end do
          if(allocated(temp)) deallocate(temp)
          if(allocated(before)) deallocate(before)
          if(allocated(after)) deallocate(after)
          nrb = bounds%nnob
          do i=1,nbc
            call set_bnd_ipar(i,'kranf',bounds%bstart(i))
            call set_bnd_ipar(i,'krend',bounds%bend(i))
          end do
	end if

	if( bstop ) stop 'error stop: ckbnds'

	end

!********************************************************************

	subroutine prbnds

        use shympi
	use mod_bnd
	use bnd_geom
        use fem_util

	implicit none

	include 'param.h'

	include 'bound_names.h'


	integer i,ibc
	integer ibtyp,kranf,krend
	integer intpol

	if( nbc .le. 0 ) return

	write(6,*)
	write(6,*) '====== info on open boundaries =========',my_id
	write(6,*)
	write(6,*) ' inlet,type,intpol,nnodes ... nodes,files : ',my_id

	do ibc=1,nbc
          call get_bnd_ipar(ibc,'ibtyp',ibtyp)
          call get_bnd_ipar(ibc,'intpol',intpol)
          call get_bnd_ipar(ibc,'kranf',kranf)
          call get_bnd_ipar(ibc,'krend',krend)

	  write(6,'(4i9)') ibc,ibtyp,intpol,krend-kranf+1,my_id
	  if( kranf > 0 .and. krend > 0 ) then
            if(bmpi) then
             write(6,'(8i9)')(ipext(bounds%niob(i)),i=kranf,krend),my_id
            else
             write(6,'(8i9)')(ipext(irv(i)),i=kranf,krend),my_id
            end if
	  end if

	  call print_filename(boundn(ibc))
	  call print_filename(conzn(ibc))
	  call print_filename(tempn(ibc))
	  call print_filename(saltn(ibc))
	  call print_filename(vel3dn(ibc))
	  call print_filename(bio2dn(ibc))
	  call print_filename(sed2dn(ibc))
	  call print_filename(mud2dn(ibc))
          call print_filename(lam2dn(ibc))
          call print_filename(dmf2dn(ibc))
	  call print_filename(tox3dn(ibc))
	  call print_filename(bfm1bc(ibc))
	  call print_filename(bfm2bc(ibc))
	  call print_filename(bfm3bc(ibc))
	end do

	write(6,*)
	write(6,*) '========================================',my_id
	write(6,*)

	return
	end

!********************************************************************

	subroutine print_filename(name)
	character*(*) name
	character*79 local
	local = name
        if( local .ne. ' ' ) then
                write(6,'(a79)') local
        end if
	end 

!********************************************************************

	subroutine tsbnds

	use mod_bnd
	use bnd_geom

	implicit none

	include 'param.h'

	include 'bound_names.h'



	integer j,i

	write(6,*) '/bnd/'
	write(6,*) nbc
	do j=1,nbc
	  write(6,*) 'ibc = ',j
	  write(6,*) boundn(j)
	  write(6,*) conzn(j)
	  write(6,*) tempn(j)
	  write(6,*) saltn(j)
	  write(6,*) vel3dn(j)
	  write(6,*) bio2dn(j)
	  write(6,*) sed2dn(j)
	  write(6,*) mud2dn(j)
          write(6,*) dmf2dn(j)
          write(6,*) lam2dn(j)
	  write(6,*) tox3dn(j)
	  write(6,*) bfm1bc(j),bfm2bc(j),bfm3bc(j)
	  write(6,*) (bnd(i,j),i=1,nbvdim)
	end do

	write(6,*) '/irv/'
	write(6,*) nrb
	do j=1,nrb
	  write(6,*) irv(j)
	end do

	end

!********************************************************************
!********************************************************************
!********************************************************************
!      utility routines
!********************************************************************
!********************************************************************
!********************************************************************

	function zvbnds(ibc)

! returns value of open boundary

	implicit none

	double precision zvbnds
	integer ibc

	double precision zval

        call chkibc(ibc,'zvbnds:')

        call get_bnd_par(ibc,'zval',zval)
	zvbnds = zval

	end

!********************************************************************

	subroutine stybnd(ibc,ibtyp)

! sets type of open boundary

	implicit none

	integer ibc
	integer ibtyp

        call chkibc(ibc,'stybnd:')

        call set_bnd_ipar(ibc,'ibtyp',ibtyp)

	end

!********************************************************************

	function itybnd(ibc)

! returns type of open boundary

	implicit none

	integer itybnd
	integer ibc

	integer ibtyp

        call chkibc(ibc,'itybnd:')

        call get_bnd_ipar(ibc,'ibtyp',ibtyp)
	itybnd = ibtyp

	end

!********************************************************************

	function nbnds()

! returns total number of open boundaries

	use mod_bnd

	implicit none

	integer nbnds


	nbnds = nbc

	end

!********************************************************************

	function nkbnd()

! returns total number of open boundary nodes

	use mod_bnd

	implicit none

	integer nkbnd


	nkbnd = nrb

	end

!********************************************************************

	function nkbnds(ibc)

! returns total number of nodes of boundary ibc

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

!********************************************************************

	function kbnd(i)

! returns i th node of all boundary nodes

	use mod_bnd
	use bnd_geom
        use shympi

	implicit none

	integer kbnd
	integer i

	include 'param.h'

        if( i .gt. nrb ) then
            write(6,*) 'i, nrb, nbc : ',i,nrb,nbc
            stop 'error stop kbnd: i out of bounds'
        end if

        if(bmpi) then
	  kbnd = bounds%niob(i)
        else
	  kbnd = irv(i)
        end if

	end

!********************************************************************

	subroutine kanfend(ibc,kranf,krend)

! returns index of first and last boundary node of boundary ibc

	implicit none

	integer ibc
        integer kranf,krend

        call chkibc(ibc,'kanfend:')

        call get_bnd_ipar(ibc,'kranf',kranf)
        call get_bnd_ipar(ibc,'krend',krend)

	end

!********************************************************************

	function kbndind(ibc,i)

! returns global index of i th node of boundary ibc

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

!********************************************************************

	function kbnds(ibc,i)

! returns i th node of boundary ibc

	use bnd_geom
        use shympi

	implicit none

	integer kbnds
	integer ibc
	integer i

	include 'param.h'

	integer kranf,krend,n

        call chkibc(ibc,'kbnds:')

        call get_bnd_ipar(ibc,'kranf',kranf)
        call get_bnd_ipar(ibc,'krend',krend)

	n = krend - kranf + 1

	if( i .lt. 1 .or. i .gt. n ) then
	    write(6,*) 'i, imax, ibc : ',i,n,ibc
	    stop 'error stop kbnds: i out of bounds'
	end if

        if(bmpi) then
          kbnds = bounds%niob(kranf+i-1)
        else
	  kbnds = irv(kranf+i-1)
        end if

	end

!********************************************************************

	subroutine irbnds(ibc,ndim,idim,nodes)

! returns nodes of boundary ibc (maximum ndim)

	use bnd_geom
        use shympi

	implicit none

	integer ibc			!number of open boundary  (in)
	integer ndim			!dimension of nodes()     (in)
	integer idim			!total number of nodes    (out)
	integer nodes(1)		!boundary nodes           (out)

	include 'param.h'

	integer i,imaxi
	integer kranf,krend

        call chkibc(ibc,'irbnds:')

        call get_bnd_ipar(ibc,'kranf',kranf)
        call get_bnd_ipar(ibc,'krend',krend)

	idim = krend - kranf + 1

	imaxi = min(ndim,idim)

        if(bmpi) then
	  do i=1,imaxi
	    nodes(i) = bounds%niob(i+kranf-1)
	  end do
        else
	  do i=1,imaxi
	    nodes(i) = irv(i+kranf-1)
	  end do
        end if

	end

!********************************************************************

        subroutine ksinget(ibc,ampli,period,phase,zref)

        double precision zref,ampli,period,phase
        call chkibc(ibc,'ksinget:')
        call get_bnd_par(ibc,'ampli',ampli)
        call get_bnd_par(ibc,'period',period)
        call get_bnd_par(ibc,'phase',phase)
        call get_bnd_par(ibc,'zref',zref)
        end

!********************************************************************

	subroutine setbnd(ibc,value,barray)

! sets boundary ibc to value in barray (apparently not used)

	use bnd_geom

	implicit none

	integer ibc
	double precision value
	double precision barray(1)

	include 'param.h'

	integer kranf,krend,k,kn

        call chkibc(ibc,'setbnd:')

        call get_bnd_ipar(ibc,'kranf',kranf)
        call get_bnd_ipar(ibc,'krend',krend)

	do k=kranf,krend
	  kn = irv(k)
	  barray(kn) = value
	end do

	end

!********************************************************************

        subroutine chkibc(ibc,errtext)
 
! checks if ibc is in bounds
 
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

!********************************************************************

	subroutine is_closed(ibc,nbc,icl)

	implicit none

	integer ibc,nbc,icl

	integer i

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

!********************************************************************

	subroutine get_oscil(ibc,rit,zval)

	implicit none

	integer ibc
	double precision rit
	double precision zval

        double precision pi
        parameter( pi=3.141592653 )

	double precision ampli,period,phase,zref

	call ksinget(ibc,ampli,period,phase,zref)

        zval = zref+ampli*sin(2.*pi*(rit+phase)/period)

	end

!********************************************************************

	subroutine get_boundary_file(ibc,what,file)

	implicit none

	integer ibc
	character*(*) what
	character*(*) file

	include 'param.h'
	include 'bound_names.h'

        if( what .eq. 'zeta' ) then
          file = boundn(ibc)
        else if( what .eq. 'conz' ) then
          file = conzn(ibc)
        else if( what .eq. 'temp' ) then
          file = tempn(ibc)
        else if( what .eq. 'salt' ) then
          file = saltn(ibc)
        else if( what .eq. 'vel' ) then
          file = vel3dn(ibc)
        else if( what .eq. 'sedt' ) then
          file = sed2dn(ibc)
        else if( what .eq. 'lagvebio' ) then
          file = bio2dn(ibc)
        else if( what .eq. 'toxi' ) then
          file = tox3dn(ibc)
        else if( what .eq. 'bfm1' ) then
          file = bfm1bc(ibc)
        else if( what .eq. 'bfm2' ) then
          file = bfm2bc(ibc)
        else if( what .eq. 'bfm3' ) then
          file = bfm3bc(ibc)
        else
          write(6,*) 'keyword not recognized: ',what
          stop 'error stop get_boundary_file'
        end if

	end

!********************************************************************

	subroutine bnd_init(ibc)

! initializes boundary ibc

	use mod_bnd

	implicit none

	integer ibc

	include 'param.h'

	integer i

	do i=1,nbvdim
	  bnd(i,ibc) = 0.
	end do

	end

!********************************************************************
!********************************************************************
!********************************************************************

        subroutine setget_bnd_par(ibc,ientry,name,value,bset)

! sets/gets value at entry ientry

	use mod_bnd

        implicit none

        integer ibc
        integer ientry
	character*(*) name
        double precision value
        logical bset

	include 'param.h'

        call chkibc(ibc,'setget_bnd_par:')

        if( ientry .le. 0 ) then
          write(6,*) 'no such entry: ',ientry,ibc,name
	  stop 'error stop setget_bnd_par: no such entry'
        end if

        if( bset ) then
          bnd(ientry,ibc) = value
        else
          value = bnd(ientry,ibc)
        end if

        end

!********************************************************************
 
        subroutine set_bnd_par(ibc,name,value)
        character*(*) name
        double precision value
        id = iget_bnd_id(name,.true.)
        call setget_bnd_par(ibc,id,name,value,.true.)
        end

        subroutine get_bnd_par(ibc,name,value)
        character*(*) name
        double precision value
        id = iget_bnd_id(name,.true.)
        call setget_bnd_par(ibc,id,name,value,.false.)
        end

        subroutine set_bnd_ipar(ibc,name,ivalue)
        character*(*) name
        double precision value
        id = iget_bnd_id(name,.true.)
	value = ivalue
        call setget_bnd_par(ibc,id,name,value,.true.)
        end

        subroutine get_bnd_ipar(ibc,name,ivalue)
        character*(*) name
        double precision value
        id = iget_bnd_id(name,.true.)
        call setget_bnd_par(ibc,id,name,value,.false.)
	ivalue = nint(value)
        end

!********************************************************************
!********************************************************************
!********************************************************************
!
! to insert new parameter increase nbvdim below and insert name in data
!
!********************************************************************
!********************************************************************
!********************************************************************

	subroutine get_bnd_nbnd(nbnd)

	use mod_bnd

        implicit none

	integer nbnd


	nbnd = nbvdim

	end

!********************************************************************

        function exists_bnd_name(name)

        implicit none

        logical exists_bnd_name
        character*(*) name

	exists_bnd_name = iget_bnd_id(name,.false.) > 0

	end

!********************************************************************

        function iget_bnd_id(name,berror)

	use mod_bnd

        implicit none

        integer iget_bnd_id
        character*(*) name
	logical berror		!raises error if name not existing


        integer id
        character*6 bname

        do id=1,nbvdim
	  call get_bnd_name(id,bname)
          if( bname .eq. name ) then
            iget_bnd_id = id
            return
          end if
        end do

	if( berror ) then
          write(6,*) 'unknown name for boundary: ',name
          stop 'error stop iget_bnd_id: name'
	end if

	iget_bnd_id = 0

        end

!********************************************************************

	subroutine get_bnd_name(id,name)

! all variables used in section bound
!
! to add a variable:
!	add to names
!	increase nbvdim (more instances in file)

	use mod_bnd

	implicit none

	integer id
        character*(*) name


        character*6 names(nbvdim)
        save names
        data names      /'iqual','ibtyp','kranf','krend','zval','ampli','period','phase','zref','ktilt' &
     &                  ,'intpol','levmax','igrad0','zfact','conz','temp','salt','levmin','kref'        &
     &			,'ztilt','tramp','levflx','nad','lgrpps','tnudge'/

	if( id .gt. nbvdim ) then
	  name = ' '
	else
	  name = names(id)
	end if

	end

!********************************************************************

	subroutine check_bnd_entries(ibc)

	use mod_bnd

	implicit none

	integer ibc


        integer i
	double precision value

	write(6,*) 'check_bnd_entries: ',ibc
        do i=1,nbvdim
          call setget_bnd_par(ibc,i,' ',value,.false.)
	  write(6,*) i,value
	end do

        end

!********************************************************************
 
	function get_discharge(ibc)

! returns discharge through boundary ibc for points sources
! for z-boundaries 0 is returned

	use bnd_dynamic

	implicit none

	double precision get_discharge
	integer ibc

	integer itype,nk,i,k
	double precision acc

	get_discharge = 0.

	itype = itybnd(ibc)
	if( itype .le. 1 .or. itype .gt. 3 ) return

	acc = 0.
        nk = nkbnds(ibc)
        do i=1,nk
          k = kbnds(ibc,i)
	  acc = acc + rqpsv(k)
        end do

	get_discharge = acc

	end

!**********************************************************************

	subroutine adjust_mass_flux

! adjusts mass flux for dry nodes

	use geom_dynamic
	use bnd_dynamic
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer k,l

        logical iskout
        iskout(k) = inodv(k).eq.-2

	do k=1,nkn
	  if( iskout(k) ) then
	    do l=1,nlv
	      mfluxv(l,k) = 0.
	    end do
	    rqv(k) = 0.
	    rqdsv(k) = 0.
	  end if
	end do

	end

!**********************************************************************

!------------------------------------------------------------------------
        end module bnd_admin
!------------------------------------------------------------------------
