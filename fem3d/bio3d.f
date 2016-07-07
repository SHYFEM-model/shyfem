c
c $Id: bio3d.f,v 1.33 2008-10-10 09:29:54 georg Exp $
c
c bio3d - EUTRO in SHYFEM
c
c contents :
c
c revision log :
c
c 20.06.2003    ggu&dmk new routine for sediments
c 20.08.2003    ggu	new routine bio_av_shell (aver/min/max)
c 03.09.2003    ggu	bug fix for sediments -> not saved (SAVESED)
c 03.09.2003    ggu	new routine check_bio
c 09.09.2003    ggu	call to scal3sh changed -> 3D version
c 19.12.2003    ggu	sediments added, init for taranto
c 14.01.2004    dmk     call tsmass per calcolo conserv massa
c 03.03.2004    ggu	decay function for bacteria decad_bio()
c 04.03.2004    ggu	changes from donata integrated
c 05.03.2004    ggu	initialization from file
c 22.03.2004    ggu	change in call to massconc (bug)
c 30.03.2004    ggu	bug fix -> call to confil with nlvdi
c 20.07.2004    dmk	new routine sed_av_shell (aver/min/max)
c 24.08.2004    ggu	new version from donata (jeanmichel)
c 24.08.2004	ggu     new check_es, changes in check_bio
c 24.08.2004	ggu     ivect(8) in bio_av_shell, 
c 24.08.2004	ggu     sedload deleted -> substituted by setload_new
c 24.08.2004	ggu	loicz moved to sediment routine
c 25.08.2004    ggu	setsedload moved to sedim routines
c 25.08.2004	ggu	light routines deleted
c 25.08.2004	ggu     new call to eutro0d implemented, read lux
c 25.08.2004	ggu     rluxaux introduced (for ntot>0)
c 17.01.2005    ggu	new horizontal diffusion
c 07.11.2005    ggu     sinking velocity wsink introduced in call to scal3sh
c 17.02.2006    ggu     pass what to subroutines to see calling routine
c 23.03.2006    ggu     ntot eliminated
c 23.03.2006    ggu     changed time step to real
c 18.10.2006    ggu     new routine custom_restime
c 18.10.2006    ggu     introduce bresi,breact,bdecay
c 17.03.2008    ggu     new open boundary routines introduced
c 08.04.2008    ggu     treatment of boundaries slightly changed
c 22.04.2008    ggu     advection parallelized
c 23.04.2008    ggu     call to bnds_set_def() changed
c 09.10.2008    ggu     new call to confop
c 08.05.2014    ggu     bug in call to inicfil for es -> must be inic2fil
c 21.10.2014    ggu     converted to new boundary treatment
c 17.05.2015    dmc     Insert benthic feeders 
c 17.06.2016    dmc     light from shyfem get_light (Watt/m2) 
c 17.06.2016    dmc     link to shyfem 7_5_13 
c 23.06.2016    ggu     bug fix: forgot to initialize eload
c
c notes :
c
c cambiamenti fatti da ggu
c
c weutro:
c bio3d:
c
c       controllare bsedim and einit
c	controllare bloicz
c       scommentare setload()
c       dati in setload sono cambiati
c       a subroutine decad_bio has been added 
c               -> (is not used normally, so ignore)
c       ho integrato il conrollo di massa, ma c'e' da controllare

c********************************************************************
c********************************************************************
c
c notes :
c
c State variables used: (Wasp)
c
c nh3		71	1
c no3		72	2
c opo4		73	3
c phyto		74	4
c cbod		75	5
c do		76	6
c on		77	7
c op		78	8
c zoo		79	9
c
c opsed         91      1
c onsed         92      2
c
c shellfarm     93      density of benthic filter feeding      
c shellsize     94      size of each individual
c shelldiag     95      diagnostic variable
c
c State variables used: (Haka)
c
c php		81	1
c zoo		82	2
c det		83	3
c dop		84	4
c dip		85	5
c
c********************************************************************

!====================================================================
        module eutro
!====================================================================

        implicit none

	integer, parameter :: nstate = 9
	integer, parameter :: nsstate = 2
	integer, parameter :: nshstate = 3

	real, save, allocatable :: e(:,:,:)	!state vector
	real, save, allocatable :: eload(:,:,:)	!loadings
	real, save, allocatable :: eseed(:,:,:)	!seed benthic filters
	real, save, allocatable :: es(:,:)	!sediment state vector
	real, save, allocatable :: esh(:,:)	!benthic filters state vector

        integer, save :: ia_out(4)
        double precision, save :: da_out(4)

        integer, save :: iubp,iubs,iubh

	logical, save :: bsedim = .false.
        logical, save :: bshell = .false.

!====================================================================
        end module eutro
!====================================================================

        subroutine ecological_module(it,dt)

c general interface to ecological module

        implicit none

        integer it
        real dt

        call bio3d_eutro(it,dt)

        end

c********************************************************************

	subroutine bio3d_eutro(it,dt)

c eco-model cosimo

	use mod_diff_visc_fric
	use levels
	use basin
	use eutro

	implicit none

	integer it	!time in seconds
	real dt		!time step in seconds

! Leslie: what is eseed? please document, does it has to be 3D?

	include 'mkonst.h'

        character*10 what,whataux
	character*2 whatn

	integer k,i,l,lmax
	integer ibio
        integer id
	integer nintp,ivar
	integer nbc
	real t,s
	real u,v

	real eaux(nstate)
	real esaux(nsstate)
        real eshaux(nshstate)
	real elaux(nstate)

	real, save :: einit(nstate)
	real, save :: esinit(nsstate)
        real, save :: eshinit(nshstate)
        real, save :: elinit(nstate)
        real, save :: ebound(nstate)

	integer, save, allocatable :: idbio(:)

        real tstot(nstate)              !for mass test
        real tsstot(nsstate)

	integer icall,iunit
	integer j
	real rlux,rluxaux,itot,fday
	real dtday
	real area,vol
	real oxysat
	real getpar
	integer iround
	integer ieint,ipint
	logical has_output,next_output

        integer mode
        real ai,lsurf

	logical bcheck,bspec
	logical bresi,breact,bdecay
	integer ie,ii
	integer kspec
	integer itanf,nvar
	double precision dtime0,dtime
	real d
	real cbod,nh3,krear,sod
	real vel
	real windspeed,tempair
	real tday,t0,tsec
	real stp
        real mass
	real wsink
        real shellfarm
        real qrad       !solar radiation Watt/m2

	integer nbnds

	integer iespecial,inspecial
	save iespecial,inspecial
	real rkpar,difmol
	save rkpar,difmol

        save icall

c------------------------------------------------------------------
c	initial and boundary conditions  [mg/l]			??
c	initial loadings                 [g/m**2/sec]		??
c
c ebound is used in case no values are given in STR file
c
c laguna di Venezia
c        data einit /0.05, 0.4, 0.01, 0.05, 2.,   11.,0.2,0.01,0.015/
c 	 data einit /0.0, 0., 0.0, 0.0, 0.,   0.,0.,0.0,0.0/
c 	 data einit /1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0/
c 	 data ebound /10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0/
c 	 data ebound /1.0, 2., 3.0, 4.0, 5.,   6.,7.,8.0,9.0/

c                     nh3 no2 opo4 phyto cbod do  on  op  zoo 
 	 data ebound  /0., 0., 0.,   0.,  0.,  0., 0., 0., 0./
 	 data einit   /0., 0., 0.,   0.,  0.,  0., 0., 0., 0./
 	 data elinit  /0., 0., 0.,   0.,  0.,  0., 0., 0., 0./
	 data esinit  /0.,0./
         data eshinit /0., 0.,0. /

c mare di taranto
c        data einit /0.042,0.355,0.009,0.0342,3.15,7.78,0.2,0.01,0.015/
c 	 data einit /0.0, 0., 0.0, 0.0, 0.,   0.,0.,0.0,0.0/
c	data esinit /0.,0./
c	data esinit /0.01,0.01/

c hakata reactor
c                     phy     zoo     det     dop     dip
c        data einit  / 1.6E-3, 1.0E-4, 1.8E-3, 3.5E-3, 2.7E-3 /
c        data einit  / 0.,     0.,     0.,     0.,     0.     /
c        data ebound / 1.6E-3, 1.0E-4, 1.8E-3, 3.5E-3, 2.7E-3 /
c        data elinit / 0.,     0.,     0.,     0.,     0.     /

c------------------------------------------------------------------
	data icall /0/
c------------------------------------------------------------------

	bresi  = .false.	!computes residence times with custom_restime
	breact = .true.		!use reactor
	bdecay = .false.	!imposes decay through decad_bio
	bcheck = .true.		!checks for out of bound values

        what = 'lagvebio'

	kspec = 0
	!kspec = 4

c-------------------------------------------------------------------
c initialization
c-------------------------------------------------------------------

	if( icall .le. -1 ) return

	if( icall .eq. 0 ) then
	  ibio = iround(getpar('ibio'))
	  if( ibio .le. 0 ) icall = -1
	  if( icall .le. -1 ) return
	  icall = 1

c         --------------------------------------------------
c	  initialize state variables with einit
c         --------------------------------------------------

	  allocate(e(nlvdi,nkndi,nstate))
	  allocate(eload(nlvdi,nkndi,nstate))
	  allocate(es(nkndi,nsstate))
          allocate(eseed(nlvdi,nkndi,nshstate))	!Leslie - do we need 3D here?
          allocate(esh(nkndi,nshstate))

	  do i=1,nstate
	    e(:,:,i) = einit(i)
	    eload(:,:,i) = elinit(i)
	  end do

	  do i=1,nsstate
	    es(:,i) = esinit(i)
          end do

          do i=1,nshstate
            esh(:,i) = eshinit(i)
          end do

c         --------------------------------------------------
c	  initialize state variables from external file
c         --------------------------------------------------

          call inicfil('bio',e,nstate)
          call inic2fil('bios',es,nsstate)

c         --------------------------------------------------
c	  set loadings in the interal areas
c         --------------------------------------------------

          call setseed_new(eseed) !Seeding for benthic filters feeding

          do i=1,nshstate
            esh(:,i) = eseed(1,:,i)	!Leslie - not clear, esh already set
          end do

c         --------------------------------------------------
c	  set boundary conditions for all state variables
c         --------------------------------------------------

          nbc = nbnds()
          allocate(idbio(nbc))
          idbio = 0

	  call get_first_time(itanf)
          dtime0 = itanf
          nintp = 2
	  nvar = nstate
          call bnds_init_new(what,dtime0,nintp,nvar,nkn,nlv
     +                          ,ebound,idbio)

c         --------------------------------------------------
c	  initialize eco model
c         --------------------------------------------------

	  call eutroini

c         --------------------------------------------------
c	  parameters for transport/diffusion resolution
c         --------------------------------------------------

          rkpar=getpar('chpar')
          difmol=getpar('difmol')

c         --------------------------------------------------
c	  initialize output 
c         --------------------------------------------------

          call eutro_init_file_output
          call eutro_write_file_output(dtime0)

	  write(6,*) 'bio3d model initialized...'
	  if( bsedim ) write(6,*) 'sediment module active...'
	  if( bshell ) write(6,*) 'shellfish module active...'

	  call loicz1(0,0.,0.)

	end if

c-------------------------------------------------------------------
c custom computation of residence times
c-------------------------------------------------------------------

	if( bresi ) call custom_restime(it,dt,e)

c-------------------------------------------------------------------
c normal call
c-------------------------------------------------------------------

	wsink = 0.

c	-------------------------------------------------------------------
c	time management
c	-------------------------------------------------------------------

	t0 = 0.
	dtday = dt / 86400.

	tsec = it
	tday = it / 86400. + t0		!time in days, FEM 0 is day t0

c	-------------------------------------------------------------------
c	loop on elements for biological reactor
c	-------------------------------------------------------------------

        mode = +1               !new time level for volume and depth

	if( bcheck ) call check_bio('before eutro',e,es)
c	call check_es(es)

	if( breact ) then	!use reactor ?

	do k=1,nkn		!loop on nodes

          lmax = ilhkv(k)
	  bspec = k .eq. kspec

          call get_light(k,qrad)

          do l=1,lmax
            call dvanode(l,k,mode,d,vol,area)   !gets depth, volume and area
            call getts(l,k,t,s)                 !gets temp and salt
            call getuv(l,k,u,v)                 !gets velocities u/v
            vel = sqrt(u*u+v*v)

            id = 1000*k+l

	    eaux(:) = e(l,k,:)
	    elaux(:) = eload(l,k,:)

	    if( bspec ) write(6,*) 'bio3d 1: ',eaux
	    if( bspec ) write(6,*) 'bio3d 1a: ',elaux

	    call eutro0d(id,tday,dtday,vol,d,vel,t,s,qrad,eaux,elaux)
            !call haka0d(tsec,dt,vol,d,t,ai,eaux,elaux)

	    if( bspec ) write(6,*) 'bio3d 3: ',eaux

	    e(l,k,:) = eaux(:)
          end do

	  l = lmax

          if( bsedim ) then
	    esaux(:) = es(k,:)
	    if( bspec ) write(6,*) 'before wsedim: ',eaux,esaux
	    call wsedim(k,tday,dtday,vol,d,vel,t,eaux,esaux)
	    if( bspec ) write(6,*) 'after wsedim: ',eaux,esaux
	    e(l,k,:) = eaux(:)
	    es(k,:) = esaux(:)
            es(k,:) = esaux(:)
          end if

          if( bshell ) then
            shellfarm=eseed(1,k,1)	!FIXME - Leslie - not clear
            if (shellfarm.gt.0) then
              eshaux(:)=esh(k,:)
              call wshell(k,tday,dtday,vol,d,vel,t,eaux,eshaux)
              esh(k,:) = eshaux(:)
              e(l,k,:) = eaux(:)
              esh(k,:) = eshaux(:)
            end if
          end if

	end do

	end if	!breact

c	-------------------------------------------------------------------
c	simplified decay for bacteria etc.
c	-------------------------------------------------------------------

	if( bdecay ) call decad_bio(e,dt)

c	-------------------------------------------------------------------
c	advection and diffusion
c	-------------------------------------------------------------------

	if( bcheck ) call check_bio('before advection',e,es)

	dtime = it
	call bnds_read_new(what,idbio,dtime)

!$OMP PARALLEL PRIVATE(i)
!$OMP DO SCHEDULE(DYNAMIC)	

	do i=1,nstate

          call scal_adv(what,i
     +                          ,e(1,1,i),idbio
     +                          ,rkpar,wsink
     +                          ,difhv,difv,difmol)

          call tsmass (e(1,1,i),1,nlvdi,tstot(i)) !mass control

	end do

!$OMP END DO NOWAIT	
!$OMP END PARALLEL	

	if( bcheck ) call check_bio('after advection',e,es)

        do i=1,nsstate
          call scalmass(es(1,i),0.1,tsstot(i))   !mass ctrl sed
	end do

        !call pn_tot(it,nstate,nsstate,tstot,tsstot)  !writes to unit 17,18,19

c	-------------------------------------------------------------------
c	write of results (file BIO)
c	-------------------------------------------------------------------

	if( bcheck ) call check_bio('before write',e,es)

        call eutro_write_file_output(dtime)

	call bio_av_shell(e)		!aver/min/max of state vars
	call sed_av_shell(es)		!aver/min/max of sed var

c	call loicz1(0,0.,0.)

	if( bcheck ) call check_bio('at end',e,es)

c	-------------------------------------------------------------------
c	debug output
c	-------------------------------------------------------------------

        !k = 100
        !l = 1
        !call getts(l,k,t,s)
        !call writee(95,it,k,l,e,t,s,nlvdi,nkndi,nstate)
        !call bioprint(it,e,nstate)

c	-------------------------------------------------------------------
c	end of routine
c	-------------------------------------------------------------------

	end

c*************************************************************

	subroutine writee(iunit,it,k,l,e,t,s,nlvddi,nknddi,nstate)

c formatted write for debug

	implicit none

	integer iunit
	integer it,k,l
	integer nlvddi,nknddi,nstate
	real e(nlvddi,nknddi,nstate)
	real t
	real s

	integer i

	write(iunit,'(i10,11f12.4)') it,
     +			(e(l,k,i),i=1,nstate),
     +			t,s

	end

c*************************************************************

	subroutine bio_av_shell(e)

c computes and writes average/min/max of bio variables
c
c id = 260
c
c e(1) average	== 261
c e(1) min	== 262
c e(1) max	== 263
c e(2) average	== 264
c ...

	use basin
	use levels

	implicit none

c parameter

	integer nstate
	parameter( nstate = 9 )

	real e(nlvdi,nkndi,nstate)	!state vector

c local
	integer idtc,itmc,itsmed
	integer id,nvar
c function
	real getpar
c save
	double precision, save, allocatable :: bioacu(:,:,:)
	real, save, allocatable :: biomin(:,:,:)
	real, save, allocatable :: biomax(:,:,:)

	save ivect
	integer ivect(8)

	integer icall
	save icall

	data icall / 0 /

	if( icall .lt. 0 ) return

	if( icall .eq. 0 ) then

          itsmed=nint(getpar('itsmed'))
          if( itsmed .le. 0 ) then
            icall = -1
            return
          end if

	  allocate(bioacu(nlvdi,nkndi,nstate))
	  allocate(biomin(nlvdi,nkndi,nstate))
	  allocate(biomax(nlvdi,nkndi,nstate))

	  idtc=nint(getpar('idtcon'))
	  itmc=nint(getpar('itmcon'))

	  nvar = nstate

	  id = 260
	  call cmed_init('bav',id,nvar,nlvdi,idtc,itmc
     +				,bioacu,biomin,biomax,ivect)

	  icall = 1
	end if

	call cmed_accum(nlvdi,e,bioacu,biomin,biomax,ivect)

	end

c*************************************************************

	subroutine check_bio(title,e,es)

c checks bio vars

	use levels, only : nlvdi,nlv
	use basin

	implicit none

	integer nstate
	parameter( nstate = 9 )
	integer nsstate
	parameter( nsstate = 2 )

	character*(*) title
	real e(nlvdi,nkndi,nstate)	!state vector
	real es(nkndi,nsstate)		!sediment state variables


        character*20 text
	integer i

c	write(6,*) 'check_bio: ',title

        text = '*** bio check e     '
	do i=1,nstate
          write(text(18:19),'(i2)') i
          call check2Dr(nlvdi,nlv,nkn,e(1,1,i),0.,1.e+20,text,title)
	end do

        text = '*** bio check es    '
	do i=1,nsstate
          write(text(18:19),'(i2)') i
          call check1Dr(nkn,es(1,i),0.,1.e+20,text,title)
	end do

	end

c*************************************************************

	subroutine decad_bio(e,dt)

c simulates decay for virus and bacteria

	use levels
	use basin

	implicit none

	integer nstate
	parameter( nstate = 9 )

        real dt
	real e(nlvdi,nkndi,nstate)	!state vector


        integer k,l,i,lmax
        real aux,tau

c----------------------------------------------------
        tau = 2.                !decay time in days
        tau = 0.                !tau = 0 => no decay
c----------------------------------------------------

        if( tau .le. 0. ) return

        tau = tau * 86400.
        aux = exp(-dt/tau)

	do k=1,nkn		!loop on nodes
          lmax = ilhkv(k)
          do l=1,lmax
	    do i=1,nstate
	      e(l,k,i) = aux * e(l,k,i)
            end do
	  end do
        end do

	end

c*************************************************************
c
c DOCS  START   S_biopar_h 
c
c DOCS  COMPULS         Compulsory bio parameters
c
c These parameters are compulsory parameters that define if the
c water quality module is run and what kind of output is written.
c
c |ibio|	Flag if the computation on the temperature is done.
c		The model writes at each time step the state 
c		variable values in the the .bio output file 
c |itsmed|	Flag if the average, minimum, maximum file of variables 
c		bio, salinity, temperature is done.
c		if |itsmed=1| the model writes |.sav|, |.tav| output files
c		of the corresponding variables.
c
cc        call addfnm('ibio',' ')
cc        call addfnm('itsmed',' ')
c
c DOCS  BIONAME		Boundary conditions
c
c Boundary conditions have to be given in a file in every 
c section |$bound|.
c
c |bio2dn|	File name that contains boundary conditions for concentration 
c		of the water quality state variables. 
c		The format is the same as for the file |boundn|. 
c		The unit of the values given in the second 
c		and following column (9 data columns for EUTRO)
c		must the ones of the variable.
c
cc        call addfnm('bio2dn',' ')
c
c DOCS  FILENAME        Initial conditions
c
c Initialization of variables are done by file. The files can be created
c by the progam |laplap|. They have to be given in
c section |$name|.
c
c |bio|		File with concentration values of water quality variable
c		to be used for the initialization.
c |salt, temp|	Files with salinity concentration values [psu] and
c		Temperature values [deg C] for the initialization.
c |conz| 	Files with tracer concentration values [%] 
c		for the initialization.
c
cc        call addfnm('bio',' ')
cc        call addfnm('salt',' ')
cc        call addfnm('temp',' ')
c
cc FIXME
c 
cc        call addfnm('conz',' ')         !file with values in time -> imposed
cc        call addfnm('salt',' ')
cc        call addfnm('temp',' ')
cc        call addfnm('bio',' ')
cc        call addfnm('bios',' ')
cc
cc        call addfnm('conzin',' ')       !not yet implemented    FIXME
cc        call addfnm('saltin',' ')
cc        call addfnm('tempin',' ')
c
c DOCS	END
c
c*************************************************************

	subroutine custom_restime(it,dt,e)

c custom routine to compute residence times with bio variables
c
c reactor must be commented
c einit must be 1.

	use levels
	use basin

	implicit none

	integer nstate
	parameter( nstate = 9 )

	integer it
	real dt
	real e(nlvdi,nkndi,nstate)	!state vector

	integer k,lmax,l,i
	integer istate,itper
	real e0,einit
	real remnant
	real secs_in_day
	real v1v(nkn)

	integer ifemop

	real restim(nstate)
	real perc(nstate)
	double precision mass(nstate)
	double precision emass(nstate)
	real mass0(nstate)
	save emass,mass0

	integer iu
	save iu

	integer icall
	save icall
	data icall / 0 /

c---------------------------------------------------------------------
c parameters to be set
c---------------------------------------------------------------------

c e0		initial concentration
c itper		period for re-initialization

	e0 = 1.
	itper = 30 * 86400
	itper = 86400

	secs_in_day = 86400.

c---------------------------------------------------------------------
c initialize accumulated mass and open file
c---------------------------------------------------------------------

	if ( icall .eq. 0 ) then
	  iu = ifemop('.jam','formatted','new')
	  do i=1,nstate
	    emass(i) = 0.
	  end do
	end if

c---------------------------------------------------------------------
c compute flag for valid nodes where residence time has to be computed
c---------------------------------------------------------------------

	call valid_node(v1v)

	istate = 0

c---------------------------------------------------------------------
c see if we must re-initialize
c---------------------------------------------------------------------

	if( mod(it,itper) .eq. 0 ) then

	istate = it/itper
	istate = mod(istate,nstate) + 1

	write(6,*) 'new initialization of state vector ',istate,itper

	do k=1,nkn		!loop on nodes
	  if( v1v(k) .ne. 0 ) then
	    einit = e0
	  else
	    einit = 0.
	  end if

          lmax = ilhkv(k)
          do l=1,lmax
	        e(l,k,istate) = einit
	  end do
        end do

	end if

c---------------------------------------------------------------------
c compute mass and initialize mass0 if necessary (first call)
c---------------------------------------------------------------------

	call comp_tot_mass(e,v1v,mass)

	if( icall .eq. 0 ) then
	  do i=1,nstate
	    mass0(i) = mass(i)
	  end do
	end if

c---------------------------------------------------------------------
c if concentrations have been re-initialized -> initialize emass and mass0
c---------------------------------------------------------------------

	if( istate .gt. 0 ) then
	  emass(istate) = 0.
	  mass0(istate) = mass(istate)
	end if

	do i=1,nstate
	  remnant = mass(i)/mass0(i)
	  emass(i) = emass(i) + remnant * dt
	  restim(i) = emass(i) / secs_in_day
	  perc(i) = 100.*remnant
	end do
	
c---------------------------------------------------------------------
c write results
c---------------------------------------------------------------------

	write(iu,'(i10,9f7.2,9f7.2)') it,(perc(i),i=1,nstate)
     +				,(restim(i),i=1,nstate)

	icall = icall + 1

c---------------------------------------------------------------------
c end of routine
c---------------------------------------------------------------------

	end

c*************************************************************

	subroutine valid_node(val)

c computes valied nodes (nodes that are inside lagoon)
c
c must be customized

	use basin

	implicit none

	real val(1)

	integer k,ie,ii,ia
	integer iaout

c-------------------------------
	iaout = -1	!area code of elements outside lagoon
c-------------------------------

        do k=1,nkn
          val(k) = 0.
        end do

        do ie=1,nel
          ia = iarv(ie)
          if( ia .ne. iaout ) then		!here put your area
              do ii=1,3
                k = nen3v(ii,ie)
                val(k) = 1.
              end do
          end if
        end do

	end

c*************************************************************

	subroutine comp_tot_mass(e,v1v,mass)

c computes total mass of state variables (only where v1v is not 0)

	use levels
	use basin

	implicit none

	integer nstate
	parameter( nstate = 9 )

	real e(nlvdi,nkndi,nstate)	!state vector
	real v1v(1)
	double precision mass(1)


	integer k,lmax,l,i
	real vol,conz

	real volnode

	do i=1,nstate
	  mass(i) = 0.
	end do

	do k=1,nkn		!loop on nodes
	  if( v1v(k) .ne. 0 ) then
            lmax = ilhkv(k)
            do l=1,lmax
	      vol = volnode(l,k,+1)
	      do i=1,nstate
	        conz = e(l,k,i)
	        mass(i) = mass(i) + vol*conz
	      end do
	    end do
	  end if
        end do

	end

c*************************************************************

c****************************************************************

        subroutine bioprint(it,e,nstate)

	use levels
	use basin

        implicit none

        integer it,nstate
	real e(nlvdi,nkndi,nstate)	!state vector

        integer i,k,n
        logical berror


        integer icall
        save icall
        data icall / 0 /

        integer ndim
        parameter (ndim=5)
        integer nodes(ndim)
        save nodes
        !data nodes / 984, 4860, 4636, 4585 /
        data nodes / 984, 4860, 4636, 4585 , 3452 /

        if( icall .eq. 0 ) then
          icall = 1
          call n2int(ndim,nodes,berror)
          if( berror ) stop 'error stop cprint'
        end if

        write(84,*)
        write(84,*) 'time = ',it
        write(84,*)

        do i=1,ndim
          k = nodes(i)
          write(84,*) i,k,(e(1,k,n),n=1,nstate)
        end do

        end

c****************************************************************

        subroutine eutro_init_file_output

        use basin
        use levels
        use eutro

        implicit none

        integer ishyff,nvar,id
        logical has_output,has_output_d
        real getpar

        ishyff = nint(getpar('ishyff'))

          call init_output('itmcon','idtcon',ia_out)
          if( ishyff == 1 ) ia_out = 0
          if( has_output(ia_out) ) then
            call open_scalar_file(ia_out,nlv,nstate,'bio')
            iubp = ia_out(4)
	    if( bsedim ) then
              call open_scalar_file(ia_out,1,nsstate,'sed')
              iubs = ia_out(4)
	    end if
	    if( bshell ) then
              call open_scalar_file(ia_out,1,nshstate,'she')
              iubh = ia_out(4)
	    end if
          end if

	  nvar = nstate
	  if( bsedim ) nvar = nvar + nsstate
	  if( bshell ) nvar = nvar + nshstate

          call init_output_d('itmcon','idtcon',da_out)
          if( ishyff == 0 ) da_out = 0
          if( has_output_d(da_out) ) then
            call shyfem_init_scalar_file('eutro',nvar,.false.,id)
            da_out(4) = id
          end if

        end

c*************************************************************

        subroutine eutro_write_file_output(dtime)

        use basin
        use levels
        use eutro

        implicit none

        double precision dtime

        integer nvar,id,idc,i
        logical next_output,next_output_d

        if( next_output(ia_out) ) then

          ia_out(4) = iubp
          do i=1,nstate
            idc = 200 + i
            call write_scalar_file(ia_out,idc,nlvdi,e(1,1,i))
          end do

	  if( bsedim ) then
            ia_out(4) = iubs
            do i=1,nsstate
              idc = 220 + i
              call write_scalar_file(ia_out,idc,1,es(1,i))
            end do
	  end if

	  if( bshell ) then
            ia_out(4) = iubh
            do i=1,nshstate
              idc = 230 + i
              call write_scalar_file(ia_out,idc,1,esh(1,i))
            end do
	  end if

        end if

        if( next_output_d(da_out) ) then

          id = nint(da_out(4))
          do i=1,nstate
            idc = 200 + i
            call shy_write_scalar_record(id,dtime,idc,nlvdi
     +                                          ,e(1,1,i))
          end do

	  if( bsedim ) then
            do i=1,nsstate
              idc = 220 + i
              call shy_write_scalar_record(id,dtime,idc,1
     +                                          ,es(1,i))
            end do
	  end if

	  if( bshell ) then
            do i=1,nsstate
              idc = 230 + i
              call shy_write_scalar_record(id,dtime,idc,1
     +                                          ,esh(1,i))
            end do
	  end if

        end if

        end

c*************************************************************

