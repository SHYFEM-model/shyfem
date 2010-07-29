c
c $Id: subbnd.f,v 1.35 2009-05-21 09:24:00 georg Exp $
c
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
c function kbnd(i)			 returns i'th node of all b. nodes
c subroutine kanfend(ibc,kranf,krend)    returns index of first and last bnode
c function kbnds(ibc,i)			 returns i'th node of boundary ibc
c subroutine irbnds(ibc,ndim,idim,nodes) returns nodes of boundary ibc
c subroutine bndsetget(ibc,ientry,value,bset) sets/gets value at entry ientry
c subroutine setbnd(ibc,value,barray)	 sets boundary ibc to value in barray
c
c subroutine setbc(value,array,flag)	 sets all open boundaries to value
c subroutine chkibc(ibc,errtext)	 checks if ibc is in bounds
c
c notes :
c
c for scalars: -999. uses ambient value
c
c ibndim is the maximum dimension used for boundary information 
c ibndim is also declared in ht
c nbvdim is the actual filling of the variables
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
c 24.08.1998    ggu     BC for concentration is bnd(20,..)
c 24.08.1998    ggu     BC for maximum input level is bnd(12,..) -> levmax
c 27.08.1998    ggu     accessor function levbnd for levmax
c 22.01.1999    ggu     new subroutine setbnd
c 08.07.1999    ggu     bug with iqual -> ibtyp not respected
c 20.01.2000    ggu     call to rdbnds without dimension -> use getdim
c 07.04.2000    ggu     new subroutine setbc
c 07.05.2001    ggu     introduced new variable zfact
c 25.09.2001    ggu     introduced bio2dn
c 07.08.2003    ggu     check for nrbdim in rdbnds
c 15.10.2004    ggu     new boundary types and sedin
c 02.03.2005    ggu     new nbdim for 3D boundary values
c 02.03.2005    ggu     some new helper functions
c 07.11.2005    ccf     introduced sed2dn
c 16.02.2006    ggu     introduced tox3dn
c 07.04.2008    acc     introduced bfm1bc bfm2bc bfm3bc OB condition for ERSEM
c 17.04.2008    ggu     deleted infobnd(), levbnd()
c 28.04.2008    ggu     call to nrdpar in double precision
c 29.04.2008    ggu&aac new boundary arrays for ERSEM
c 30.05.2008	ggu	eliminated numbers for parameters
c 03.06.2008	ggu	new parameters levmin, kref
c 06.06.2008	ggu	completely restructured
c 02.04.2009	ggu	intpol default is 0, some unused routines deleted
c 20.04.2009	ggu	new variable ztilt, ndim substituted with nbvdim
c
c************************************************************************

	subroutine inbnds

c initializes boundary parameters

	implicit none

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	nbc = 0
	nrb = 0
	nrz = 0
	nrq = 0

	end

c************************************************************************

	subroutine rdbnds(ibc)

c reads boundary info from STR file

	implicit none

	integer ibc

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer irv(1)
        common /irv/irv

	character*80 boundn(1)
        common /boundn/ boundn
	character*80 conzn(1)
        common /conzn/ conzn
	character*80 tempn(1)
        common /tempn/ tempn
	character*80 saltn(1)
        common /saltn/ saltn
	character*80 bio2dn(1)
        common /bio2dn/ bio2dn
        character*80 sed2dn(1)
        common /sed2dn/ sed2dn
        character*80 tox3dn(1)
        common /tox3dn/tox3dn
	character*80 bfm1bc(1)
        common /bfm1bc/bfm1bc
        character*80 bfm2bc(1)
        common /bfm2bc/bfm2bc
        character*80 bfm3bc(1)
        common /bfm3bc/bfm3bc

	character*80 name,text
	double precision dvalue
	real value
	integer nbcdim,nrbdim
	integer i,kranf,krend,kref
	integer iweich,id,nbnd
	integer nrdpar

	real getpar

	call getdim('nbcdim',nbcdim)
	call getdim('nrbdim',nrbdim)

	nbc = nbc + 1

	if(nbc.ne.ibc) goto 94
	if(nbc.gt.nbcdim) goto 77

	call bnd_init(nbc)	!sets everything to 0

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
c		corresponds to an anti-clockwise sense. At least
c		two nodes must be given.

	call addpar('kbound',0.)

c |ibtyp|	Type of open boundary. 
c		\begin{description}
c		\item[0] No boundary values specified
c		\item[1] Level boundary. At this open boundary
c			 the water level is imposed and the prescribed
c			 values are interpreted as water levels in meters.
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
c		been defined before.

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

c |conzn, tempn, saltn|	File name that contains values for the respective
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

c The next variables specify the name of the boundary value file
c for different modules. Please refer to the documentation of the
c single modules for the units of the variables.

c |bio2dn|	File name that contains values for the ecological
c		module (EUTRO-WASP).
c |sed2dn|	File name that contains values for the sediment
c		transport module
c |tox3dn|	File name that contains values for the toxicological
c		module.

	call addfnm('bio2dn',' ')       !HACK
	call addfnm('sed2dn',' ')       !HACK
	call addfnm('tox3dn',' ')       !HACK

cc File name for OB condition in ERSEM MODULE - undocumented

	call addfnm('bfm1bc',' ')
	call addfnm('bfm2bc',' ')
	call addfnm('bfm3bc',' ')

c |intpol|	Order of interpolation for the boundary values read
c		through files. Use for 1 for stepwise (no) interpolation,
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
c |umom, vmom|		Constant values for momentum input. (Default 0)

	call addpar('umom',0.)
	call addpar('vmom',0.)

cc undocumented

	call addpar('nbdim',0.)
	call addpar('ktilt',0.)
	call addpar('ztilt',0.)
	call addpar('kref',0.)
	call addpar('igrad0',0.)		!use 0 gradient for scalars

c DOCS	END

c here add dummy variables

	call addpar('zval',0.)

cccccccccccccccccccccccccccccccccccccccccccccccccc
c here we start the reading loop
cccccccccccccccccccccccccccccccccccccccccccccccccc

	kranf = nrb + 1
	call set_bnd_ipar(ibc,'kranf',kranf)	!position of starting node
	call addpar('kranf',float(kranf))

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
		if( nrb .gt. nrbdim ) goto 76
		irv(nrb)=nint(value)
	    else if( iweich .eq. 3 ) then	!file name
	        ! must be handeled later
	    else if( iweich .ne. 0 ) then
                call set_bnd_par(ibc,name,value)
	    end if
	end do

	krend = nrb
	call set_bnd_ipar(ibc,'krend',krend)	!position of end node
	call addpar('krend',float(krend))


	call get_nbnd(nbnd)
	do id=1,nbnd
	  call get_bnd_name(id,name)
	  value = getpar(name)
          call set_bnd_par(ibc,name,value)
	end do

	call getfnm('boundn',boundn(nbc))
	call getfnm('conzn',conzn(nbc))
	call getfnm('tempn',tempn(nbc))
	call getfnm('saltn',saltn(nbc))

	call getfnm('bio2dn',bio2dn(nbc))
	call getfnm('sed2dn',sed2dn(nbc))
	call getfnm('tox3dn',tox3dn(nbc))

	call getfnm('bfm1bc',bfm1bc(nbc))
	call getfnm('bfm2bc',bfm2bc(nbc))
	call getfnm('bfm3bc',bfm3bc(nbc))

	!call check_bnd(ibc)

	call check_parameter_values('before deleting section')
	call delete_section('bound')
	call check_parameter_values('after deleting section')

	return
   76	continue
	write(6,*) 'Dimension error for nrbdim'
	write(6,*) 'nrbdim :',nrbdim
	stop 'error stop : rdbnds'
   77	continue
	write(6,*) 'Dimension error for nbcdim'
	write(6,*) 'nbcdim :',nbcdim
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

c********************************************************************

	subroutine ckbnds

c checks boundary information read from STR

	implicit none

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer irv(1)
	common /irv/irv
	character*80 boundn(1)
        common /boundn/ boundn

	logical bstop
	integer i,k,ibc
	integer iqual,ibtyp,kranf,krend,ktilt,knode,kref
	integer levmax,levmin
	integer intpol
	real period
	real ztilt
	integer ipint

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

         call get_bnd_ipar(ibc,'kranf',kranf)
         call get_bnd_ipar(ibc,'krend',krend)
	 if( kranf .gt. krend ) then	!$$kranf
	   write(6,'(a,i2,a)') 'section BOUND ',i,' :'
	   write(6,*) '   No nodes given for boundary'
	   bstop=.true.
	 else
	   if(ibtyp.eq.1) then
		nrz=nrz+krend-kranf+1
	   else if(ibtyp.eq.2) then
		nrq=nrq+krend-kranf+1
	   end if
	 end if

         call get_bnd_par(ibc,'period',period)
	 if( period .le. 0. .and. boundn(i) .eq. ' ' ) then
		write(6,'(a,i2,a)') 'section BOUND ',i,' :'
		write(6,*) '   Period must be > 0'
		write(6,*) '   period = ',period
		bstop=.true.
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

	 do k=kranf,krend
	    knode=ipint(irv(k))		!$$EXTINW
	    if(knode.le.0) then
              write(6,'(a,i2,a)') 'Section BOUND ',i,' :'
              write(6,*) '   boundary node not found ',irv(k)
              bstop=.true.
	    end if
	    irv(k)=knode
	 end do

	end do

	if( bstop ) stop 'error stop: ckbnds'

	end

c********************************************************************

	subroutine prbnds

	implicit none

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	character*80 boundn(1)
        common /boundn/ boundn
	character*80 conzn(1)
        common /conzn/ conzn
	character*80 tempn(1)
        common /tempn/ tempn
	character*80 saltn(1)
        common /saltn/ saltn
	character*80 bio2dn(1)
        common /bio2dn/ bio2dn
	character*80 sed2dn(1)
        common /sed2dn/ sed2dn
        character*80 tox3dn(1)
        common /tox3dn/tox3dn
	character*80 bfm1bc(1)
        common /bfm1bc/bfm1bc
        character*80 bfm2bc(1)
        common /bfm2bc/bfm2bc
        character*80 bfm3bc(1)
        common /bfm3bc/bfm3bc

	integer irv(1)
	common /irv/irv

	integer i,ibc
	integer ibtyp,kranf,krend
	integer intpol
	integer ipext

	if( nbc .le. 0 ) return

	write(6,*)
	write(6,1003)

	do ibc=1,nbc
          call get_bnd_ipar(ibc,'ibtyp',ibtyp)
          call get_bnd_ipar(ibc,'intpol',intpol)
          call get_bnd_ipar(ibc,'kranf',kranf)
          call get_bnd_ipar(ibc,'krend',krend)
	  write(6,1036) ibc,ibtyp,intpol,(ipext(irv(i)),i=kranf,krend)

	  call print_filename(boundn(ibc))
	  call print_filename(conzn(ibc))
	  call print_filename(tempn(ibc))
	  call print_filename(saltn(ibc))
	  call print_filename(bio2dn(ibc))
	  call print_filename(sed2dn(ibc))
	  call print_filename(tox3dn(ibc))
	  call print_filename(bfm1bc(ibc))
	  call print_filename(bfm2bc(ibc))
	  call print_filename(bfm3bc(ibc))
	end do

	return
 1003   format(' inlet,type,intpol,file,nodes : ')
 1036   format(i7,3x,i5,2x,i5/(10x,10i6))
	end

c********************************************************************

	subroutine print_filename(name)
	character*(*) name
	character*79 local
	local = name
        if( local .ne. ' ' ) then
                write(6,'(a79)') local
        end if
	end 

c********************************************************************

	subroutine tsbnds

	implicit none

	integer ibndim
	parameter (ibndim=100)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	character*80 boundn(1)
        common /boundn/ boundn
	character*80 conzn(1)
        common /conzn/ conzn
	character*80 tempn(1)
        common /tempn/ tempn
	character*80 saltn(1)
        common /saltn/ saltn
	character*80 bio2dn(1)
        common /bio2dn/ bio2dn
	character*80 sed2dn(1)
        common /sed2dn/ sed2dn
        character*80 tox3dn(1)
        common /tox3dn/tox3dn
	character*80 bfm1bc(1)
        common /bfm1bc/bfm1bc
        character*80 bfm2bc(1)
        common /bfm2bc/bfm2bc
        character*80 bfm3bc(1)
        common /bfm3bc/bfm3bc

        real bnd(ibndim,1)
        common /bnd/bnd
	integer irv(1)
	common /irv/irv

	integer j,i

	write(6,*) '/bnd/'
	write(6,*) nbc
	do j=1,nbc
	  write(6,*) 'ibc = ',j
	  write(6,*) boundn(j)
	  write(6,*) conzn(j)
	  write(6,*) tempn(j)
	  write(6,*) saltn(j)
	  write(6,*) bio2dn(j)
	  write(6,*) sed2dn(j)
	  write(6,*) tox3dn(j)
	  write(6,*) bfm1bc(j),bfm2bc(j),bfm3bc(j)
	  write(6,*) (bnd(i,j),i=1,ibndim)
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

	implicit none

	integer nbnds

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	nbnds = nbc

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

	end

c********************************************************************

	function kbnd(i)

c returns i'th node of all boundary nodes

	implicit none

	integer kbnd
	integer i

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer irv(1)
	common /irv/irv

	integer krend

        call get_bnd_ipar(nbc,'krend',krend)  !here we use nbc (last boundary)

        if( i .gt. krend ) then
            write(6,*) 'i, krend, nbc : ',i,krend,nbc
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

	function kbndind(ibc,i)

c returns global index of i'th node of boundary ibc

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

c returns i'th node of boundary ibc

	implicit none

	integer kbnds
	integer ibc
	integer i

	integer irv(1)
	common /irv/irv

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

	implicit none

	integer ibc			!number of open boundary  (in)
	integer ndim			!dimension of nodes()     (in)
	integer idim			!total number of nodes    (out)
	integer nodes(1)		!boundary nodes           (out)

	integer irv(1)
	common /irv/irv

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

	implicit none

	integer ibc
	real value
	real barray(1)

	integer irv(1)
	common /irv/irv

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

	subroutine setbc(value,array,flag)

c sets all open boundaries to value

	implicit none

	real value
	real array(1)
	real flag

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer irv(1)
	common /irv/irv

	integer ibc
	integer ibtyp,kranf,krend
	integer k,kn

	do k=1,nkn
	  array(k) = flag
	end do

	do ibc=1,nbc

          call get_bnd_ipar(ibc,'ibtyp',ibtyp)
          call get_bnd_ipar(ibc,'kranf',kranf)
          call get_bnd_ipar(ibc,'krend',krend)

          do k=kranf,krend

             kn=irv(k)

             if( ibtyp .eq. 1 .or. ibtyp .eq. 2 ) then
               array(kn)=value
	     end if

	  end do
	end do

	end

c********************************************************************

        subroutine chkibc(ibc,errtext)
 
c checks if ibc is in bounds
 
        implicit none
 
        integer ibc
	character*(*) errtext
 
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
 
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

	subroutine get_oscil(ibc,rit,zval)

	implicit none

	integer ibc
	real rit
	real zval

        real pi
        parameter( pi=3.141592653 )

	real ampli,period,phase,zref

	call ksinget(ibc,ampli,period,phase,zref)

        zval = zref+ampli*sin(2.*pi*(rit+phase)/period)

	end

c********************************************************************

	subroutine get_boundary_file(ibc,what,file)

	implicit none

	integer ibc
	character*(*) what
	character*(*) file

        character*80 boundn(1)
        common /boundn/boundn
        character*80 conzn(1)
        common /conzn/conzn
        character*80 tempn(1)
        common /tempn/ tempn
        character*80 saltn(1)
        common /saltn/ saltn

        if( what .eq. 'zeta' ) then
          file = boundn(ibc)
        else if( what .eq. 'conz' ) then
          file = conzn(ibc)
        else if( what .eq. 'temp' ) then
          file = tempn(ibc)
        else if( what .eq. 'salt' ) then
          file = saltn(ibc)
        else
          write(6,*) 'keyword not recognized: ',what
          stop 'error stop get_boundary_file'
        end if

	end

c********************************************************************

	subroutine bnd_init(ibc)

c initializes boundary ibc

	implicit none

	integer ibc

	integer ibndim
	parameter (ibndim=100)

        real bnd(ibndim,1)
	common /bnd/bnd

	integer i

	do i=1,ibndim
	  bnd(i,ibc) = 0.
	end do

	end

c********************************************************************

        subroutine bndsetget(ibc,ientry,value,bset)

c sets/gets value at entry ientry

        implicit none

        integer ibc
        integer ientry
        real value
        logical bset

	integer ibndim
	parameter (ibndim=100)

        real bnd(ibndim,1)
	common /bnd/bnd

        call chkibc(ibc,'bndsetget:')

        if( bset ) then
          bnd(ientry,ibc) = value
        else
          value = bnd(ientry,ibc)
        end if

        end

c********************************************************************
 
        subroutine set_bnd_par(ibc,name,value)
        character*(*) name
        id = iget_bnd_id(name)
        call bndsetget(ibc,id,value,.true.)
        end

        subroutine get_bnd_par(ibc,name,value)
        character*(*) name
        id = iget_bnd_id(name)
        call bndsetget(ibc,id,value,.false.)
        end

        subroutine set_bnd_ipar(ibc,name,ivalue)
        character*(*) name
        id = iget_bnd_id(name)
	value = ivalue
        call bndsetget(ibc,id,value,.true.)
        end

        subroutine get_bnd_ipar(ibc,name,ivalue)
        character*(*) name
        id = iget_bnd_id(name)
        call bndsetget(ibc,id,value,.false.)
	ivalue = nint(value)
        end

c********************************************************************

	subroutine get_nbnd(nbnd)

        implicit none

	integer nbnd

        integer nbvdim
        parameter(nbvdim=21)

	nbnd = nbvdim

	end

c********************************************************************

        function iget_bnd_id(name)

        implicit none

        integer iget_bnd_id
        character*(*) name

        integer nbvdim
        parameter(nbvdim=21)

        integer id
        character*6 bname

	!write(6,*) 'iget_bnd_id: looking for ',name

        do id=1,nbvdim
	  call get_bnd_name(id,bname)
          if( bname .eq. name ) then
            iget_bnd_id = id
            return
          end if
        end do

        write(6,*) 'unknown name for boundary: ',name
        stop 'error stop iget_bnd_id: name'

        end

c********************************************************************

	subroutine get_bnd_name(id,name)

	implicit none

	integer id
        character*(*) name

        integer nbvdim
        parameter(nbvdim=21)

        character*6 names(nbvdim)
        save names
        data names      /
     +                   'iqual','ibtyp','kranf','krend','zval'
     +                  ,'ampli','period','phase','zref','ktilt'
     +                  ,'intpol','levmax','igrad0','zfact','nbdim'
     +			,'conz','temp','salt','levmin','kref'
     +			,'ztilt'
     +                  /

	if( id .gt. nbvdim ) then
	  name = ' '
	else
	  name = names(id)
	end if

	end

c********************************************************************

	subroutine check_bnd(ibc)

	implicit none

	integer ibc

        integer nbvdim
        parameter(nbvdim=21)

        integer i
	real value

	write(6,*) 'check_bnd: ',ibc
        do i=1,nbvdim
          call bndsetget(ibc,i,value,.false.)
	  write(6,*) i,value
	end do

        end

c********************************************************************

