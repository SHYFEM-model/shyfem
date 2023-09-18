!--------------------------------------------------------------------
        module bnd_admin_par
!--------------------------------------------------------------------
        contains
!--------------------------------------------------------------------

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

!--------------------------------------------------------------------
!
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
	real value
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
!
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
!
! |bio2dn|	File name that contains values for the ecological
!		module (EUTRO-WASP).
! |sed2dn|	File name that contains values for the sediment
!		transport module.
!		The unit of the values given
!		in the second and following columns (equal to the 
!		number of defined grainsize in parameter |sedgrs|).
!
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

!         File name for OB condition in ERSEM MODULE - undocumented

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
!
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
!
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

! |umom, vmom|		Constant values for momentum input. (Default 0)

	call addpar('umom',0.d0)
	call addpar('vmom',0.d0)

! undocumented

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

	subroutine mkbnds

! makes  boundary information on different processes 

	use mod_bnd
	use bnd_geom
	!use shypart
        use mpi_communication_struct
        use fem_util

	implicit none

	include 'param.h'

	logical bstop
        integer istop
	integer i,k,ibc
	integer iqual,ibtyp,kranf,krend,ktilt,knode,kref
	integer levmax,levmin
	integer intpol
	real period
	real ztilt
	character*80 file
        integer,dimension(:),allocatable :: temp
        integer icount,gid
        integer,dimension(:),allocatable,save :: totkranf,totkrend
        integer,save :: icall


	bstop = .false.

        bounds%nob = nbc
        if (bounds%nob.gt.0 .and.(.not.allocated(bounds%bstart))) then
          allocate(bounds%bstart(bounds%nob))
          allocate(bounds%bend(bounds%nob))
          allocate(bounds%tnob(bounds%nob))
        end if
        allocate(temp(size(irv)))
        bounds%niob = 0
        bounds%bstart(1) = 1
        bounds%nnob = 0

        if(icall .eq.0) then
          allocate(totkranf(nbc))
          allocate(totkrend(nbc))
	  do ibc=1,nbc
            call get_bnd_ipar(ibc,'kranf',kranf)
            call get_bnd_ipar(ibc,'krend',krend)
            totkranf(ibc) = kranf
            totkrend(ibc) = krend
          end do
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

         !call get_bnd_ipar(ibc,'kranf',kranf)
         !call get_bnd_ipar(ibc,'krend',krend)
         kranf = totkranf(ibc)
         krend = totkrend(ibc)
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
	 do k=kranf,krend
	    if( k == 0 ) cycle
            if(icall .le. 0) then
	      knode=ipint(irv(k))		!$$EXTINW
	      irv(k)=knode
            else
              knode=irv(k)
            end if
	    if(knode.le.0) then
              !write(6,'(a,i2,a)') 'Section BOUND ',i,' :'
              !write(6,*) '   boundary node not found ',irv(k)
              istop = istop + 1
	    end if
            if(knode .gt. 0) then
              icount = icount + 1
              temp(icount+bounds%nnob)=knode
            end if
	 end do

         if( istop > 0 ) then
           !if( bmpi ) then
           !  if( istop == krend-kranf+1 ) then
           !   !write(6,*) 'boundary completely in one domain... ok',my_id
           !    call set_bnd_ipar(ibc,'ibtyp',0)
           !    bounds%nob = bounds%nob - 1
           !  end if
           !else
	   !  stop 'error stop mkbnds: no MPI and missing nodes'
           !end if
         end if
         if(bounds%iobound .gt. 1) then
           bounds%bstart(bounds%iobound) = &
           bounds%bend(bounds%iobound-1)+1
         end if
         bounds%nnob = bounds%nnob + icount
         bounds%bend(bounds%iobound)=bounds%bstart(bounds%iobound)+icount-1
	end do

        if(.not. allocated(bounds%niob)) then
          allocate(bounds%niob(bounds%nnob))
        else
          deallocate(bounds%niob)
          allocate(bounds%niob(bounds%nnob))
        end if
        do i=1,bounds%nnob
          bounds%niob(i) = temp(i)
        end do
        if(allocated(temp)) deallocate(temp)
        nrb = bounds%nnob
        do i=1,nbc
          call set_bnd_ipar(i,'kranf',bounds%bstart(i))
          call set_bnd_ipar(i,'krend',bounds%bend(i))
        end do

        icall = icall+1
	if( bstop ) stop 'error stop: mkbnds'

        !if(bmpi) then
          if(allocated(bounds%last)) then
            deallocate(bounds%last)
            allocate(bounds%last(bounds%nnob))
          else
            allocate(bounds%last(bounds%nnob))
          end if
          bounds%last=0
          do i=1,size(irv)-1
            if((irv(i).ne.0).and.(irv(i+1).eq.0)) then
              do k=1,bounds%nnob
                if(bounds%niob(k).eq.irv(i)) then
                  bounds%last(k)=1
                end if
              end do
            end if
          end do
	!end if


	end

!##########################################################################!

        subroutine get_bnd_ipar(ibc,name,ivalue)
        character*(*) name
        id = iget_bnd_id(name,.true.)
        call setget_bnd_par(ibc,id,name,value,.false.)
	ivalue = nint(value)
        end

!##########################################################################!

        subroutine set_bnd_ipar(ibc,name,ivalue)
        character*(*) name
        id = iget_bnd_id(name,.true.)
	value = ivalue
        call setget_bnd_par(ibc,id,name,value,.true.)
        end

!##########################################################################!

        subroutine get_bnd_par(ibc,name,value)
        character*(*) name
        id = iget_bnd_id(name,.true.)
        call setget_bnd_par(ibc,id,name,value,.false.)
        end

!##########################################################################!

	subroutine get_bnd_nbnd(nbnd)

	use mod_bnd

        implicit none

	integer nbnd


	nbnd = nbvdim

	end

!##########################################################################!
 
        subroutine set_bnd_par(ibc,name,value)
        character*(*) name
        id = iget_bnd_id(name,.true.)
        call setget_bnd_par(ibc,id,name,value,.true.)
        end

!##########################################################################!

        subroutine setget_bnd_par(ibc,ientry,name,value,bset)

! sets/gets value at entry ientry

	use mod_bnd

        implicit none

        integer ibc
        integer ientry
	character*(*) name
        real value
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

	function nbnds()

! returns total number of open boundaries

	use mod_bnd

	implicit none

	integer nbnds


	nbnds = nbc

	end function

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

	end function

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

	end function

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

!************************************************************************

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

        end function

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
        data names      / &
     &                   'iqual','ibtyp','kranf','krend','zval' &
     &                  ,'ampli','period','phase','zref','ktilt' &
     &                  ,'intpol','levmax','igrad0','zfact' &
     &			,'conz','temp','salt','levmin','kref' &
     &			,'ztilt','tramp','levflx','nad','lgrpps' &
     &			,'tnudge' &
     &                  /

	if( id .gt. nbvdim ) then
	  name = ' '
	else
	  name = names(id)
	end if

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


!##########################################################################!
!##########################################################################!
!##########################################################################!

!--------------------------------------------------------------------
        end module bnd_admin_par
!--------------------------------------------------------------------
