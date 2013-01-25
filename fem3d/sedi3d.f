! ***********************************
! ------ SUBROUTINE SEDI3D ---------
! ***********************************
!
! This routine manages the sediment transport computation.
! It calls the subrourine SEDTRANS05 which computes
! the sediment transport rate and the suspended sediment
! concentration in case of current or current and wave
! (see SUBWAVES routine). The routine SCAL3SH handles the
! avdection and diffusion of suspended sediment.
!
C 16/11/2004 Christian Ferrarin ISMAR-CNR
C
C revision log :
C
C Mar, 2005     ccf     (sedi3d_f1.f) coming from sedi3d_e3_multi7.f
C                       	new cohesive routine
C Mar, 2005     ccf     (sedi3d_f2.f) add mixing thickness
C Mar, 2005     ccf     (sedi3d_f3.f) merge layer 1 and 2 if bedn(1) < bmix
C Mar, 2005     ccf     (sedi3d_f4.f) update element depth
C Mar, 2005     ccf     (sedi3d_f5.f) get viscosity from new routine
C Mar, 2005     ccf     (sedi3d_f6.f) active layer = bottom roughness heigh
C Apr, 2005     ccf     (sedi3d_f7.f) style changes, create sedt05
C Apr, 2005     ccf     (sedi3d_f8.f) initialize percbd from file
C Apr, 2005     ccf     (sedi3d_f8.f) change rhos per cohesive sed
C Jun, 2005     ccf     (sedi3d_f13.f) adapt to 3D, bottom layer and total depth
C                       	bugfix error in computing bedn(1,3), 
C				dimension in tuek
C Jun, 2005     ccf     (sedi3d_f15.f) change deposition density,
C				add consolidation
C Jun, 2005     ccf     (sedi3d_f16.f) change units, bug fix in TCONC1
C Jun, 2005     ccf     (sedi3d_f17.f) bug fix in upedepth
C Jun, 2005     ccf     (sedi3d_f18.f) bug fix in checkbed, 
C				adapt to sedtrans05-h5.f
C Jul, 2005     ccf     (sedi3d_f20.f) adapt to sedtrans05-h6.f
C Jul, 2005     ccf     (sedi3d_f21.f) bug fix in updepth (as subdry.f)
C Jul, 2005     ccf     (sedi3d_f22.f) bug fix in bedload
C Aug, 2005     ccf     (sedi3d_f23.f) eliminate ripple variables
C Aug, 2005     ccf     (sedi3d_f24.f) bed slope contribution. adjust updepth
C Aug, 2005     ccf     (sedi3d_f25.f) change deposition characteristics
C Sep, 2005     ccf     (sedi3d_f26.f) bug fix in bedman,
C				change boundary for cohesive
C Sep, 2005     ccf     (sedi3d_f27.f) separete erosion and deposition in bedman
C                       	number of layer computed each time
C Nov, 2005     ccf     (sedi3d_f28f) adapt for 3D version
C Nov, 2005     ccf     (sedi3d_f29f) bed slope on threshold, bug fix in updepth
C Nov, 2005     ccf     (sedi3d_f30f) bug fix in bedslope and in getmd
C Nov, 2005     ccf     (sedi3d_f31f) bed slope by gradients
C Nov, 2005     ccf     (sedi3d_f32f) last layer not smaller than 0.1 m
C Nov, 2005     ccf     (sedi3d_f33f) compute vertical mixing coeffcients
C Jan, 2006     ccf     (sedi3d_f34f) bug fix in upedepth and blimit
C Feb, 2006     ccf     (sedi3d_f35f) new suspco routine
C Feb, 2006     ccf     (sedi3d_f36f) adapt to Van Rijn (C0) in sedtrans05-h8.f
C Feb, 2006     ccf     (sedi3d_f37f) correct bugs in checkbed and other things
C May, 2006     ccf     (sedi3d_f38f) correct bugs in line 1119. Introduce KCOES
C May, 2006     ccf     (sedi3d_f39f) bugs nonco. non used edr in cohse.
C                       	no limit percbd. pers(1)>0.1 in line 1120
C May, 2006     ccf     (sedi3d_f40f) limit BEDCHA(1,2) to 0.1. 
C				bugfix in suspco, better conc in cohse
C Jun, 2006     ccf     (sedi3d_f41f) no transport in case of depth< 0.1m
C Jun, 2006     ccf     (sedi3d_f42f) read constants from file
C Jun, 2006     ccf     (sedi3d_f43f) add limcoh
C Jun, 2006     ccf     (sedi3d_f44f) smooth bed elevation change, get_timestep
C Jul, 2006     ccf     (sedi3d_f45f) read angle of repose, 
C				limit shear velocity, write bathymetry
c 11.04.2008    ggu&ccf treatment of boundaries slightly changed
c 16.04.2008    ggu&ccf bugfix calling lin (must pass double precision 0.)
c 22.04.2008    ggu     advection parallelized
c 27.04.2008    ccf     new check, other changes
c 28.04.2008    ggu     call to nrdpar, nrdnls read with double precision
c 29.04.2008    ccf     save temref and salref in common
c 09.10.2008    ggu     new call to confop
c 28.04.2009    ggu     links re-structured
c 25.01.2013    ggu     mixed real/double common blocks spearated
c
!****************************************************************************

      subroutine sedi(it,dt)

      implicit none

      integer it			!time in seconds
      real dt				!time step in seconds

      include 'param.h'
      include 'sed_param.h'

! -------------------------------------------------------------
! fem variables
! -------------------------------------------------------------

      integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
      common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
      real eps1,eps2,pi,flag,high,higi
      common /mkonst/ eps1,eps2,pi,flag,high,higi

      integer ilhkv(nkndim)		!number of element and node level
      common /ilhkv/ilhkv
      real difv(0:nlvdim,1)
      common /difv/difv
      real difhv(nlvdim,1)
      common /difhv/difhv
      real sedkpar,difmol

      real hkv(nkndim)
      common /hkv/hkv
      real v1v(nkndim)
      common /v1v/v1v

      integer nlvdi,nlv			!number of levels
      common /level/ nlvdi,nlv
      integer lmax			!number of element levels
      integer ius1,ius2,itmcon,idtcon 	!output parameter

! -------------------------------------------------------------
! local sediment variables
! -------------------------------------------------------------

      integer icall
      integer isedi			!call parameter
      integer nscls			!number of grainsize classes
      integer nbed			!initial number of bed layers
      integer i,is,k,l,j		!counters
      integer nintp,ivar

      double precision gs(nsdim)	!SEDIMENT GRAIN DIAMETER (M)
      double precision ws(nsdim)	!SETTLING VELOCITY FOR NON-COHESIVE SEDIMENT (M/SEC)
      double precision bdh(nkndim)      !total elevation change [>0depo,<0ero]
      real bh(nkndim)		    	!bottom height variation [m]
      real tcn(nlvdim,nkndim) 		!total sediment concentration [kg/m3]
      real gskm(nkndim)			!AVERAGE SEDIMENT GRAIN DIAMETER ON NODES (M)
      real sbound(nsdim)                !boundary vector [kg/m3]
      real scn(nlvdim,nkndim,nsdim)     !non-cohesive suspended sediment conc (kg/m3)
      real scc(nlvdim,nkndim,nsdim)     !cohesive suspended sediment conc (kg/m3)
      real eps(0:nlvdim,nkndim,nsdim)	!vertical mixing coefficient
      real thick			!initial thickness [m]
      real conref			!initial concentration [kg/m3]
      real dist(nsdim)			!for cohesive boundary condition
      real wsink			!settling velocity
      real gdx(nkndim),gdy(nkndim)	!slope gradients
      real tao(nkndim)                  !wave-current shear stress
      double precision riph(nkndim)     !ripple height [m]
      double precision ripl(nkndim)     !ripple length [m]
      double precision timedr		!Duration of call to Sedtrans (s)
      double precision percbd(nlbdim,nsdim,nkndim)        !fraction of sediment [0,1]
      double precision bedn(nlbdim,3,nkndim)    !bed characteristics in 3 column table
                                                ! (1) depth below sediment surface (m)
                                                ! (2) critical erosion stress (Pa)
                                                ! (3) dry bulk density (kg/m**3)
      common /percbd/percbd
      common /bedn/bedn
      double precision sflx(nsdim,nkndim)	!flux of suspend sediment [m3/m2]
      double precision bflx(nsdim,nkndim)	!flux of bedload sediment [m3/m2]
      double precision sedx(nsdim,nkndim)	!x bedload component [kg/m2s]
      double precision sedy(nsdim,nkndim)       !y bedload component [kg/m2s]

      real sedpa(6)                     !sediment parameter vector
      common /sedpa/sedpa
      real gsc(nsdim)                   !grainsize class
      common /gsc/gsc

      real bnd3_sed(nb3dim,0:nbcdim)	!array containing boundary state
      real bnd3_aux(nb3dim)		!aux array for boundaries
      save bnd3_sed

      character*80 sed2dn(1)
      common /sed2dn/sed2dn

      real const3d(0:nlvdim,nkndim)
      common /const3d/const3d

      real hzoff
      integer isstart

      real salref,temref		!salinity [psu] and temperature [C]
      common /temsal/salref,temref

      character*10 what
      real tsec
      real z0bk(nkndim)                   !bottom roughenss on nodes
      common /z0bk/z0bk

      common /sedaux1/scn,scc,eps,sflx,sedx,sedy,gdx,gdy,tao,gskm
      common /sedaux1d/riph,ripl
      common /sedaux2/tcn,bh
      common /sedaux2d/bflx,bdh

! function
      integer iround
      real getpar

! save and data
      save ius1,ius2,itmcon,idtcon
      save isstart
      save sedkpar,difmol
      save gs,sbound
      save nscls,dist
      save hzoff
      save /sed2dn/
      save /gsc/
      save /sedpa/
      save /percbd/
      save /bedn/
      save /sedaux1/,/sedaux1d/
      save /sedaux2/,/sedaux2d/
      save /temsal/

      save icall
      data icall /0/

! ----------------------------------------------------------
! Documentation
! ----------------------------------------------------------
!
! nsdim		dimension of grain size classes
! nscls		number of different grain size classes
! nbcc		number of cohesive grain size classes (automatic)
! gs(i)		value for single grain size classes
!		only first entry of gs may be cohesive
!		if this is the fact, then nbcc has a default value of 7
!		e.g., 7 cohesive grain size sub-classes will be used
!		but boundary conditions will be specified with one value
!		...for the first class (gs(1))
!		if gs(1) is non cohesive, no cohesive classes will be used
!		nbcc = 0
! dist(i)	distribution of grain size sub-classes in cohesive class
!
! sedx/y	-> maybe invert indices
!
! ----------------------------------------------------------
! Initialization
! ----------------------------------------------------------
! This section is called only the first time step when ICALL = 0

        if( icall .le. -1 ) return

        isedi = nint(sedpa(1))
	if( isedi .le. 0 ) then
	  icall = -1
	  return
	end if

        if( icall .eq. 0 ) then
          itmcon = iround(getpar('itmcon'))
          idtcon = iround(getpar('idtcon'))
          if( it .lt. itmcon ) return
          icall = 1

!         --------------------------------------------------
!	  Initializes state variables and constants
!         --------------------------------------------------

          call readsedconst		!initialize constants

          conref = sedpa(2)		!initial concentration
          sedkpar = sedpa(3)		!diffusion parameter
          nscls = nint(sedpa(4))	!number of grain size classes
          difmol = getpar('difmol')	!molecular vertical diffusivity [m**2/s]
	  hzoff = getpar('hzon')
          nbed = 5			!initial number of bed layer
          thick = 0.0			!initial bed layer thickness [m]
          temref = getpar('temref')
          salref = getpar('salref')

          do is = 1,nsdim
           gs(is) = 0.
           ws(is) = 0.
          end do

          do is = 1,nscls
           gs(is) = gsc(is)
	   sbound(is) = 0.
          end do
	  if (gs(1) .gt. limcoh) nbcc = 0

	  isstart = 1
	  if (nbcc .gt. 0 ) isstart = 2

          do k = 1,nkn
            bh(k) = 0.
            gdx(k) = 0.
            gdy(k) = 0.
	    riph(k) = 0.
	    ripl(k) = 0.
            do l=1,nlvdim
              tcn(l,k) = 0.			!total suspended sediment
              do is=1,nsdim
                scn(l,k,is) = conref		!suspended sediment (sand)
                scc(l,k,is) = conref		!suspended sediment (cohesive)
                sedx(is,k) = 0.			!bedload
                sedy(is,k) = 0.			!bedload
                eps(l,k,is) = difv(l,k)		!vertical diffusion for sand
                eps(0,k,is) = difv(0,k)		!vertical diffusion
              end do
            end do
          end do

!         --------------------------------------------------
!         Initializes bed configuration
!         --------------------------------------------------

          call inibed(gs,nscls,nbed,thick,gskm)

!         --------------------------------------------------
!         Sets boundary conditions for all state variables
!         --------------------------------------------------

          do i = 1,nsdim
            dist(i) = 0.
          end do
          do i = 1,ndistr
            j=(wsclay+1)-meddistr+i-1
            j=max(min(j,nbcc),1)
            dist(j) = distr(i)
          end do

	  nintp = 2
          call bnds_init(what,sed2dn,nintp,nscls,nb3dim,bnd3_sed,sbound)
	  !call bnds_set_def(what,nb3dim,bnd3_sed,bnd3_aux)

!         --------------------------------------------------
!	  Initializes output
!         --------------------------------------------------

          ius1 = 0
          call confop(ius1,itmcon,idtcon,1,4,'sed')

          ius2 = 0
          call confop(ius2,itmcon,idtcon,nlv,1,'sco')

          write(6,*) 'sediment model initialized...'

	endif

! -------------------------------------------------------------------
! Normal call
! -------------------------------------------------------------------

        what = 'sedt'

        TIMEDR = dt
        tsec = it

!       -------------------------------------------------------------
!       Computes bed slope gradient
!       -------------------------------------------------------------

        call tvd_grad_2d(hkv,gdx,gdy,v1v)

!       -------------------------------------------------------------------
!       Computes sediment transport on nodes
!       -------------------------------------------------------------------

        call sed_loop(it,timedr,nscls,gs,ws,hzoff)

!       -------------------------------------------------------------------
!       Computes bedload transport
!       -------------------------------------------------------------------

        if(gs(1).gt.limcoh .or. nscls.gt.1) call bedload(nscls,sedx,
     @		sedy,dt,bflx)

!	-------------------------------------------------------------------
!	Updates bed configuration and compute bed elevation change      
!	-------------------------------------------------------------------

        call bedman(gs,nscls,timedr,bflx,sflx,gdx,gdy,tao,bh,gskm,bdh)

!	-------------------------------------------------------------------
!       Smooths node total depth variation
!	-------------------------------------------------------------------

        call smooth_node(bh,bdh,smooth,gdx,gdy,angrep)

!	-------------------------------------------------------------------
!       Updates total depth and velocities
!	-------------------------------------------------------------------

        call upedepth(bdh)

!       -------------------------------------------------------------
!       Boundary condition for suspended sediment
!       -------------------------------------------------------------

	call scal_bnd(what,tsec,bnd3_sed)

!       -------------------------------------------------------------------
!       Transport and diffusion for each sediment class
!       -------------------------------------------------------------------

!       ------------------------
!       Cohesive sediments
!       ------------------------

	ivar = 1	!this indicates only where to take boundary condition

!$OMP PARALLEL PRIVATE(is,wsink)
!$OMP DO SCHEDULE(DYNAMIC)

        do is = 1,nbcc

	  wsink = wsi(is)	!ggu error -> should be wsi(is)??
          call scal_adv_fact(what,ivar,dist(is)
     +                          ,scc(1,1,is),bnd3_sed
     +                          ,sedkpar,wsink,const3d
     +                          ,difhv,difv,difmol)

        end do

!$OMP END DO NOWAIT
c !$OMP END PARALLEL

!       ------------------------
!       Noncohesive sediments
!       ------------------------

c !$OMP PARALLEL PRIVATE(is,wsink)
!$OMP DO SCHEDULE(DYNAMIC)

        do is = isstart,nscls

          wsink = ws(is)
          call scal_adv(what,is
     +                          ,scn(1,1,is),bnd3_sed
     +                          ,sedkpar,wsink
     +                          ,difhv,eps(1,1,is),difmol)

        end do

!$OMP END DO NOWAIT
!$OMP END PARALLEL

        call totcon(nscls,scn,scc,tcn)			!total susp conc.

!       -------------------------------------------------------------------
!       Writes of results (files SED and SCO)
!       -------------------------------------------------------------------

        call confil(ius1,itmcon,idtcon,80,1,bh)
        call confil(ius1,itmcon,idtcon,81,1,gskm)
        call confil(ius1,itmcon,idtcon,82,1,tao)
        call confil(ius1,itmcon,idtcon,83,1,hkv)

        call confil(ius2,itmcon,idtcon,84,nlvdim,tcn)

! -------------------------------------------------------------------
! End of routine
! -------------------------------------------------------------------

	end

! ********************************************************************
! SUBROUTINE READGS
! This subroutine reads the simulation sediment parameter from the .str file: 
! grainsize classes, the initial sediment distribution (percentage) and the 
! critical erosion threshol value from the .str file

        subroutine readsed

        implicit none

        include 'sed_param.h'

        character*80 name		!name of item read
        character*80 text		!text of item if string
        character*80 line 
        real value			!value of item if numeric
	double precision dvalue		!value of item if numeric
        integer iweich
        integer nrdpar
        real getpar
        integer i

! --- output variables
        real sedpa(6)			!sediment parameter vector
        common /sedpa/sedpa
        integer nrs			!number of grainsize classes
        real gsc(nsdim)			!grainsize class
        common /gsc/gsc
        integer npi			!
        integer nps
        real prin(0:nsdim,nsdim)        !initial percetage [0,1]
        common /prin/prin
        integer nrst			!number of sediment tipe
        real tue(nsdim)			!initial erosion threshold
        common /tue/tue        

c DOCS  START	S_sedpar_h
c 
c DOCS  COMPULS        Compulsory sediment parameters
c
c These parameters are compulsory parameters that define if the
c sediment transport module is run and what kind of output is written.
c
c |sedtr|	Sediment section name.

        call sctpar('sedtr')             !sets default section
        call sctfnm('sedtr')

c |isedi|	Flag if the computation on the sediment is done:
c		\begin{description}
c		\item[0] Do nothing (default)
c		\item[1] Compute sediment transport
c		\end{description}

        call addpar('isedi',0.)

c |sedref|	Sediment initial reference concentration [kg/m3]
c		(Default 0).

        call addpar('sedref',0.)

c |sedhpar|	Sediment diffusion coefficient (Default 0).

        call addpar('sedhpar',0.)

c |sedgrs|	Sediment grainsize class vector [m]. Values has be 
c		ordered from the finest to the more coarse. \\
c		example: sedgrs = 0.0001 0.0002 0.0003 0.0004

        call addpar('sedgrs',0.)

c |percin|	Initial sediment distribution [0,1] for each 
c		grainsize class. The sum of percin must be equal to 1. \\
c		example: percin = 0.25 0.25 0.25 0.25 \\
c		If percin is not selected the model impose equal
c		percentage for each grainsize class (percin = 1/nrs).
c		In case of spatial differentiation of the sediment
c		distribution set a number of percin equal to the number
c		of grainsize classes per the number of area type. \\
c		example:  \\
c			percin = 0.25 0.25 0.25 0.25  \\
c				0.20 0.20 0.30 0.30 \\
c				0.45 0.15 0.15 0.15 \\

        call addpar('percin',0.)

c |tauin|	Initial dry density or TAUCE. In function of the value: \\
c		0-50 : critical erosion stress (Pa) \\
c		\textgreater 50  : dry bulk density of the surface (kg/m**3). \\
c		In case of spatial differentiation set a number of tauin
c		equal to the number of area type. \\
c		example: tauin = 0.9 1.4 2.5 1.1

        call addpar('tauin',0.5)

c |sedp|	File containing the initial sediment distribution. Values
c		are in percentage of each class.

        call addfnm('sedp',' ')

c |sedt|	File containing the initial critical erosion stress (Pa)
c		or dry bulk density (kg/m3).

        call addfnm('sedt',' ')

c |sedcon|	File containing the constants used in sediment model

        call addfnm('sedcon',' ')

c DOCS  FILENAME        Boundary conditions
c
c Boundary conditions have to be given in a file in every
c section |bound|.
c
c |sed2dn|	File name that contains boundary conditions for
c        	concentration of the sediment grainsize variables.
c        	The format is the same as for the file |boundn|.
c        	The unit of the values given in the second
c        	and following columns (nscls data columns for SEDTRANS)
c        	must the ones of the variable.
c
c DOCS  FILENAME        Initial conditions
c
c Initialization of variables are done by file. The files can be created
c by the progam |laplap|. They have to be given in section |name|.
c
c |sed|		File with concentration values of suspended sediment
c     		concentration to be used for the initialization.
c
c DOCS  END

!       --------------------------------------------------
!       Initializes variables
!       --------------------------------------------------

        do i = 1,nsdim
         tue(i) = 0.
        end do

        npi = 1
        nps = 0
        nrs = 0
        nrst = 0

!       --------------------------------------------------
!       Starts the reading loop
!       --------------------------------------------------

        iweich = 1
        do while(iweich.ne.0)
          iweich=nrdpar('sedtr',name,dvalue,text)
	  value = dvalue

!         --------------------------------------------------
!         Reads sediment grainsize vector with dimension nrs
!         --------------------------------------------------
          if( name .eq. 'sedgrs' ) then
              nrs=nrs+1
              if( nrs .gt. nsdim ) goto 35
              gsc(nrs)=value
          end if

!         --------------------------------------------------
!         Reads sediment percentage matrix for each type
!         --------------------------------------------------
          if( name .eq. 'percin' ) then
              nps=nps+1
              if( nps .gt. nsdim ) goto 35
              if( npi .gt. nsdim ) goto 35
              if( nps .gt. nrs ) then
                nps = 1
                npi = npi + 1
              end if
              prin(npi,nps)=value
          end if

!         --------------------------------------------------
!         Reads initial erosion threshold vector
!         --------------------------------------------------
          if( name .eq. 'tauin' ) then
              nrst=nrst+1
              if( nrst .gt. nsdim ) goto 35
              tue(nrst)=value
          end if
        end do

        if ( nrs .eq. 0.) goto 30

!       --------------------------------------------------
!       Case with only one sediment class
!       --------------------------------------------------
        if( nrs .eq. 1) then
          npi = 1
          nps = 1
          prin(npi,nps) = 1
        end if

!       --------------------------------------------------
!       Case with no initial distribution selected
!	impose the same percentage for each class
!       --------------------------------------------------
        if ( nps .eq. 0. ) then
          do i = 1,nrs
            prin(npi,i) = 1./nrs
          end do
          nps = nrs
        end if

!       --------------------------------------------------
!       Stores parameters in variable sedpa()
!       --------------------------------------------------

        sedpa(1) = getpar('isedi')
        sedpa(2) = getpar('sedref')
        sedpa(3) = getpar('sedhpar')
        sedpa(4) = nrs
        sedpa(5) = npi
        sedpa(6) = nrst

        return

  30    continue
        write(6,*) 'No sediment class selected'
        write(6,*) 'in the .str file'
        stop 'error stop : readsed'

  35    continue
        write(6,*) 'Dimension error for nsdim'
        write(6,*) 'nrbdim  :',nsdim
        write(6,*) 'nscls :',nrs
        stop 'error stop : readsed'

        end

! ********************************************************************
! SUBROUTINE READCONST
! This subroutine read constants from file sedcon and save them in the
! common block. The format of the file is: variable = value

        subroutine readsedconst

        implicit none

        include 'sed_param.h'

        integer nrdnls,iw,ioff
        integer ifileo
        integer iunit
        real value
	double precision dvalue
        integer nc		!number of variables
        parameter(nc = 38)
        character*80 file,name,text,line,cname(nc)
        double precision cvalue(nc)
        integer j

	iw = 0
	iunit = 0
	value = 0
        file = ''

!       --------------------------------------------------------
!       Assigns constants to cname and default values to cvalue
!       --------------------------------------------------------

        cname(1) = 'CSULVA'	!coefficient for the solid transmitted stress by Ulva
        cvalue(1) = 159.4d0
        cname(2) = 'TMULVA'	!threshold of motion of Ulva (Pa)
        cvalue(2) = 1.054d-3
        cname(3) = 'TRULVA'	!threshold of full resuspension of Ulva (Pa)
        cvalue(3) = 0.0013d0
        cname(4) = 'E0'		!Minimum erosion rate
        cvalue(4) = 1.95d-5
        cname(5) = 'RKERO'	!Erosion proportionality coefficient
        cvalue(5) = 5.88d0
        cname(6) = 'WSCLAY'	!primary median Ws class (must be in the range 1:NBCONC)
        cvalue(6) = 5
        cname(7) = 'CDISRUPT'	!constant for turbulent floc disruption during erosion
        cvalue(7) = 0.001d0
        cname(8) = 'CLIM1'	!lower limit for flocculation (kg/m**3)
        cvalue(8) = 0.1d0
        cname(9) = 'CLIM2'	!limit between simple and complex flocculation equation (kg/m**3)
        cvalue(9) = 2d0
        cname(10) = 'KFLOC'   	!constant K for flocculation equation
        cvalue(10) = 0.001d0
        cname(11) = 'MFLOC'   	!constant M for flocculation equation
        cvalue(11) = 1d0
        cname(12) = 'RHOCLAY' 	!density of clay mineral
        cvalue(12) = 2600d0
        cname(13) = 'CTAUDEP' 	!scaling factor for TAUCD
        cvalue(13) = 1d0
        cname(14) = 'PRS'     	!Resuspension probability (range 0-1)
        cvalue(14) = 0d0
        cname(15) = 'RHOMUD'  	!Density of the freshly deposited mud
        cvalue(15) = 50d0	!50d0
        cname(16) = 'DPROFA'	!Constants for density profile
        cvalue(16) = 470d0	!FINAL(I) = A - B * EXP(-C*ABOVE(I)) - D*EXP(-E*ABOVE(I))
        cname(17) = 'DPROFB'	!A : final deep density
        cvalue(17) = 150d0	!A-B-D : final surface density
        cname(18) = 'DPROFC'	!define the shape (in conjunction with B and C)
        cvalue(18) = 0.015d0
        cname(19) = 'DPROFD'
        cvalue(19) = 0d0
        cname(20) = 'DPROFE'
        cvalue(20) = 0d0
        cname(21) = 'CONSOA'	!time constant of consolidation
        cvalue(21) = 1d-5
        cname(22) = 'TEROA'	!Constants for erosion threshold from density and overlaying mass
        cvalue(22) = 6d-10
        cname(23) = 'TEROB'
        cvalue(23) = 3d0
        cname(24) = 'TEROC'
        cvalue(24) = 3.47d0
        cname(25) = 'TEROD'
        cvalue(25) = -1.915d0
        cname(26) = 'KCOES'	!Fraction of mud for sediment to be choesive
        cvalue(26) = 0.15
        cname(27) = 'CDRAGRED'	!constant for the drag reduction formula
        cvalue(27) = -0.0893d0
        cname(28) = 'Z0COH'	!BED ROUGHNESS LENGHT FOR COHESIVE SEDIMENTS
        cvalue(28) = 2.0D-4
        cname(29) = 'FCWCOH'	!FRICTION FACTOR FOR COHESIVE SEDIMENTS
        cvalue(29) = 2.2D-3
        cname(30) = 'LIMCOH'	!Limit of cohesive sediment grainsize
        cvalue(30) = 0.000063
        cname(31) = 'SMOOTH'	!smoothing factor for morphodynamic
        cvalue(31) = 1.0
        cname(32) = 'ANGREP'	!angle of repose
        cvalue(32) = 32.0
        cname(33) = 'IOPT'	!SEDIMENT TRANSPORT FORMULA OPTION NUMBER
        cvalue(33) = 5
        cname(34) = 'MORPHO'	!Morphological acceleration factor
        cvalue(34) = 1
        cname(35) = 'RHOSED'    !sediment grain density
        cvalue(35) = 2650
        cname(36) = 'POROS'	!bed porosity [0,1]
        cvalue(36) = 0.4
        cname(37) = 'DOCOMPACT'	!if not zero, call COMPACT routine
        cvalue(37) = 0
        cname(38) = 'NBCC'	!number of Ws classes = elements in CONC and in WSI
        cvalue(38) = 7

!       --------------------------------------------------------
!       Gets sedcon file name
!       --------------------------------------------------------

        call getfnm('sedcon',file)
        write(*,*)

        if( file .ne. ' ' ) then
           iunit = ifileo(0,file,'form','old')
           write(*,*)'Sediment constants initialized from file: ',
     @file(1:40)
           if( iunit .le. 0 ) then
            write(6,'(a,a)') 'filename: ',file(1:65)
            stop 'error stop readcon: Cannot open file'
           end if

!         --------------------------------------------------------
!         Reads first line in file
!         --------------------------------------------------------

          read(iunit,'(a)') line
          ioff = 1

!         --------------------------------------------------------
!         Loop on lines
!         --------------------------------------------------------

    1     continue
            iw = nrdnls(name,dvalue,text,line,ioff)
	    value = dvalue
            if( iw .le. 0 ) then
              read(iunit,'(a)',end=2) line
              ioff = 1
            else
              if( iw .eq. 1 ) text = ' '
              if( iw .eq. 2 ) value = 0
              do j = 1,nc
                if (name .eq. cname(j)) cvalue(j) = value
              end do
            end if
            goto 1
    2     continue

        end if

!       --------------------------------------------------------
!       First intializes constants from SUBROUTINE INICONST in sedtrans
!       --------------------------------------------------------

        NBCC = cvalue(38)

        call iniconst(NBCC)

!       --------------------------------------------------------
!       Assigns new values to constants
!       --------------------------------------------------------
        
        CSULVA = cvalue(1)
        TMULVA = cvalue(2)
        TRULVA = cvalue(3)
        E0 = cvalue(4)
        RKERO = cvalue(5)
        WSCLAY = cvalue(6)
        CDISRUPT = cvalue(7)
        CLIM1 = cvalue(8)
        CLIM2 =  cvalue(9)
        KFLOC = cvalue(10)
        MFLOC = cvalue(11)
        RHOCLAY = cvalue(12)
        CTAUDEP = cvalue(13)
        PRS = cvalue(14)
        RHOMUD = cvalue(15)
        DPROFA = cvalue(16)
        DPROFB = cvalue(17)
        DPROFC = cvalue(18)
        DPROFD = cvalue(19)
        DPROFE = cvalue(20)
        CONSOA = cvalue(21)
        TEROA = cvalue(22)
        TEROB = cvalue(23)
        TEROC = cvalue(24)
        TEROD = cvalue(25)
        KCOES = cvalue(26)
        CDRAGRED = cvalue(27)
        Z0COH = cvalue(28)
        FCWCOH = cvalue(29) 
        LIMCOH = cvalue(30) 
        SMOOTH = cvalue(31) 
        ANGREP = cvalue(32) / (45./atan (1.))	!from deg to rad
        IOPT = cvalue(33)
        MORPHO = cvalue(34)
        RHOSED = cvalue(35)
        POROS = cvalue(36)
	DOCOMPACT = cvalue(37)

!       --------------------------------------------------------
!       Writes constant value on the screen
!       --------------------------------------------------------

        write(6,*)'Constants for sediment transport:'
	write(6,*)
        do j = 1,nc
          write(6,44)cname(j),cvalue(j)
        end do
        write(6,*)

44      format(3x,a10,f12.7)
        end    

! ********************************************************************
! SUBROUTINE INIBED
! This subroutine set the initial bed conformation in function of
! sediment grainsize classes and the initial percentage defined in
! the .str file. The distribution is the same in each bed level (no
! vertical stratification). 

        subroutine inibed(gs,nscls,nbed,thick,gskm)

        implicit none

        include 'param.h'
        include 'sed_param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer iarv(1)                 	!element code number
        common /iarv/iarv
        integer nen3v(3,1)           	        !node number
        common /nen3v/nen3v

        real sedpa(6)	                        !sediment parameter vector
        common /sedpa/sedpa
        real prin(0:nsdim,nsdim)                !initial percetage [0,1]
        common /prin/prin
        real tue(nsdim)             		!initial erosion threshold/density
        common /tue/tue
        double precision gs(nsdim)              !sediment grain diameter (m)
        double precision percbd(nlbdim,nsdim,nkndim)        !fraction of sediment [0,1]
        double precision bedn(nlbdim,3,nkndim)  !bed characteristics in 3 column table
                                                ! (1) depth below sediment surface (m)
                                                ! (2) critical erosion stress (Pa)
                                                ! (3) dry bulk density (kg/m**3)
        common /percbd/percbd
	common /bedn/bedn
        integer nscls				!number of grainsize classes
        integer nbed				!initial number of bed layers
        real thick				!initial thickness [m]
        real gskm(nkndim)			!average grainsize per node [m]
        real gsa
	real ptot
        double precision rhosa			!initial erosion threshold
        real rhossand				!function for density from % sand
        integer itype   	                !element code number
        real pers(nlvdim,nkndim,nsdim)		!percbd initialized from file
        real tuek(nlvdim,nkndim,1)		!tauce initialized from file
        integer npi
        integer k,ib,is,ie,ii,l

	is = 0
        npi = sedpa(5)

!       -------------------------------------------------------------------
!       Initializes sediment distribution in function of bed type
!       -------------------------------------------------------------------

        do ie = 1,nel
          itype = iarv(ie)
          if ( npi.eq.1 ) itype = 1

          do ii = 1,3
            k = nen3v(ii,ie)
            tuek(1,k,1) = tue(itype)
            do ib = 1,nbed
              do is = 1,nscls
                pers(ib,k,is) = prin(itype,is)
              end do
            end do
          end do
        end do

!       -------------------------------------------------------------------
!       Initializes sediment fraction and tuek from external file
!       -------------------------------------------------------------------

        call inicfil('sedp',pers,nscls)
        call inicfil('sedt',tuek,1)

!       -------------------------------------------------------------------
!       Initializes bed thickness if thick .ne. 0
!       -------------------------------------------------------------------

        if (thick.gt.0.) then
          do k = 1,nkn
            bedn(1,1,k) = 0.
            do ib = 2,nbed
             bedn(ib,1,k) = bedn(ib-1,1,k) + thick
            end do
          end do
        end if

!       -------------------------------------------------------------------
!       Initializes bed characteristics and compute initial average grainsize
!       -------------------------------------------------------------------

        do k = 1,nkn
          gsa = 0.
	  ptot = 0.
          do is = 1,nscls
            do ib = 1,nbed
             percbd(ib,is,k) = pers(1,k,is)
            end do
            ptot = ptot + percbd(1,is,k)
            gsa = gsa + gs(is)*percbd(1,is,k)
          end do
          if(ptot.lt.0.99.or.ptot.gt.1.01) go to 130
          gskm(k) = gsa

          rhosa = rhossand(gs(1),percbd(1,1,k),1.2)
          if(tuek(1,k,1).gt.0.) rhosa = tuek(1,k,1)

          call bedini(bedn(1,1,k),nlbdim,rhosa)
        end do

	return

 130    continue
        write(6,*) 'Error in computing the sediment fraction: 3'
        write(6,*) 'total percentage:',ptot, 'node:',k
        do is = 1,nscls
          write(6,*) 'classnumber:',is,'percentage:',percbd(1,is,k)
        enddo
        stop 'error stop : percbd'

        end

! *******************************************************
! FUNCTION RHOSSAND
! returns the dry bulk density from sand fraction (Allersma, 1988)

        function rhossand(gs,pers,consc)

        implicit none

	include 'sed_param.h'

        double precision gs                     !grainsize of first class (m)
        double precision pers                   !fraction of first class
        real consc                              !consolidation coefficient [0-2.4]
        real persand                            !sand fraction [0,1]
        real rhossand                           !density function [kg/m3]

        persand = 1.
        if(gs.lt.limcoh) then 
          persand = 1. - pers
          rhossand = 480.*consc + (1300. - 280.*consc)*(persand**0.8)
        else
          rhossand = RHOSED * (1 - POROS)
        endif

        end

! *******************************************************
! SUBROUTINE ZVEL
! This subroutine computes z, the height of uz (the main current 
! velocity) above seafloor, assuming a log profile of the velocity 
! on the depth, based on the law: U = (USTAR/K)*LN(Z/ZO)

	subroutine zvel(uz,d,gd,z0bk,z)

	implicit none

! --- input variables
	double precision d	!water depth (m)
	double precision uz	!current at height z above the seafloor (m/s)
	double precision gd	!sediment grain diameter (m)
	real z0bk		!bottom roughness length from sedtrans [m]

! --- output variables
	double precision z	!height of uz above seafloor [m]

! --- local variables
	real ustar		!shear velocity 
	real br			!bottom roughness height [m]
	real z0			!bottom roughness length [m]
	real k			!von karman costant
	parameter(k=0.4)

        br = 2.5*gd		!see amos article
        z0 = br/30.
	z0 = max(z0,z0bk)
        ustar = (uz*k*d/z0)/(d/z0*(log(d/z0)-1)+1)
        z = z0*exp(uz*k/ustar)
        if (uz .lt. 0.0001) z=d/2.
	
	end

! *******************************************************
! SUBROUTINE GETMD
! This subroutine compute velocity and direction from u and v

        subroutine getmd(u,v,m,d)

        implicit none

! --- input variables
        real u,v		!x,y components

! --- output variables
        double precision m      !velocity
        double precision d	!direction

! --- local variables
        real rad
        real alfa

        rad=90./asin(1.)

        m = sqrt( u**2 + v**2 )		!get m

        alfa=atan(v/u)*(90./asin(1.))
        if(u.eq.0.) then
          alfa = 0.
          if(v.lt.0.) alfa = 180.
        end if
        if(u.gt.0.and.v.gt.0.or.u.gt.0.and.v.lt.0) alfa=90.-alfa
        if(u.lt.0.and.v.gt.0.or.u.lt.0.and.v.lt.0) alfa=270.-alfa

        d = alfa

        end

! *******************************************************
! SUBROUTINE VMIXCOEF
! This subroutine compute vertical mixing coefficient coefficient
! following Van Rijn 1993

        subroutine vmixcoef(l,dep,h,wsink,ustc,ub,dcw,ht,per,dx,eps)

        implicit none

        include 'param.h'

! --- input variables

        integer l			!number of levels
        real dep(nlvdim)		!thickness of layer [m]
        double precision h		!total water depth [m]
        double precision wsink		!settling velocity [m/s]
        double precision ustc		!current shear velocity [m/s]
        double precision ub 		!wave orbital velocity [s]
        double precision dcw		!height of the wave-current boundary layer
        double precision ht		!significant wave height [m]
        double precision per		!significant wave persion [s]
        double precision dx		!dimensionless particle diameter

! --- output variables

        real eps(0:nlvdim)

! --- local variables

        integer m
        real d				!depth from the bottom [m]
        real bet			!ratio of sediment and fluid mixing [1-1.5]
        real epsc			!current related vertical mixcoef 
        real epscmax			!max current related vertical mixcoef 
        real epsw			!wave related vertical mixcoef 
        real epswb			!wave related vertical mixcoef near the bed
        real epswmax			!wave related vertical mixcoef in the upper layer

!       ------------------------------------------------------------
!       Initializes variables
!       ------------------------------------------------------------

        d = 0.
        epsc = 0.
        epscmax = 0.
        epsw = 0.
        epswb = 0.
        epswmax= 0.

        bet = 1. + 2.* (wsink/ustc)**2
        bet = max(bet,1.)
        bet = min(bet,1.5)
        epscmax = 0.25 * bet * ustc * h
        epswb = 0.004 * dx * dcw * ub
        epswmax = 0.035*h*ht/per
        if(ub.lt.0.1) epswmax = 0.

        do m = l,1,-1
          eps(m) = 0.
          d = d + dep(m)

!         ------------------------------------------------------------
!         Current related mixing
!         ------------------------------------------------------------

          if (d .lt. 0.5*h) then
            epsc = epscmax - epscmax*(1.-(2.*d/h)**2)
          else
            epsc = epscmax
          end if

!         ------------------------------------------------------------
!         Wave related mixing
!         ------------------------------------------------------------

          if (d .le. dcw) then
            epsw = epswb
          elseif(d .gt. dcw .and. d .le. 0.5*h) then
            epsw = epswb + (epswmax-epswb)*((d-dcw)/(0.5*h-dcw))
          else
            epsw = epswmax
          end if 

!         ------------------------------------------------------------
!         Combined current and wave mixing coefficient
!         ------------------------------------------------------------

          eps(m) = sqrt(epsc**2 + epsw**2)

        end do

        eps(0) = eps(1)

        end

! *******************************************************
! SUBROUTINE SED_LOOP

        subroutine sed_loop(it,timedr,nscls,gs,ws,hzoff)

        implicit none

        include 'param.h'
        include 'sed_param.h'

! -------------------------------------------------------------
! local variables
! -------------------------------------------------------------
	integer it
        integer nscls			!number of grainsize classes
        integer i,is,k,l,m		!counters
        double precision dl		!thickness of bottom layer 
        double precision gs(nsdim)	!SEDIMENT GRAIN DIAMETER (M)
        double precision ws(nsdim)	!SETTLING VELOCITY FOR NON-COHESIVE SEDIMENT (M/SEC)
        real scn(nlvdim,nkndim,nsdim)   !non-cohesive suspended sediment conc (kg/m3)
        real scc(nlvdim,nkndim,nsdim)   !cohesive suspended sediment conc (kg/m3)
        double precision riph(nkndim)	!ripple height [m]
        double precision ripl(nkndim)	!ripple length [m]
        real gskm(nkndim)		!AVERAGE SEDIMENT GRAIN DIAMETER ON NODES (M)
        real eps(0:nlvdim,nkndim,nsdim)	!vertical mixing coefficient
        real u,v			!x and y current components
        real gdx(nkndim),gdy(nkndim)	!slope gradients
        real tao(nkndim)		!wave-current shear stress
        double precision sflx(nsdim,nkndim)	!flux of suspend sediment [m3/m2]
        double precision sedx(nsdim,nkndim)	!x bedload component [kg/m2s]
        double precision sedy(nsdim,nkndim)    	!y bedload component [kg/m2s]
        double precision percbd(nlbdim,nsdim,nkndim)        !fraction of sediment [0,1]
        double precision bedn(nlbdim,3,nkndim)	!bed characteristics in 3 column table
                                                ! (1) depth below sediment surface (m)
                                                ! (2) critical erosion stress (Pa)
                                                ! (3) dry bulk density (kg/m**3)
        common /percbd/ percbd
        common /bedn/ bedn
        double precision scnd(nsdim)    	!non-cohesive suspended sediment conc (kg/m3)
        double precision sccd(nbcc)	   	!cohesive suspended sediment conc (kg/m3)
	real hzoff
        real z0bk(nkndim)			!bottom roughenss on nodes

        common /z0bk/z0bk
        common /sedaux1/scn,scc,eps,sflx,sedx,sedy,gdx,gdy,tao,gskm
        common /sedaux1d/riph,ripl
 
        real salref,temref		!salinity [psu] and temperature [C]
	common /temsal/salref,temref

! -------------------------------------------------------------
! fem variables
! -------------------------------------------------------------

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer ilhkv(nkndim)		!number of element and node level
        common /ilhkv/ilhkv
        real hdknv(nlvdim,nkndim)
        common /hdknv/hdknv
        real getpar
        real salt,temp			!salinity [psu] and temperature [C]

! -------------------------------------------------------------
! wave variables
! -------------------------------------------------------------

        real waveh(nkndim)        !wave height [m]
        real wavep(nkndim)        !wave period [s]
        real waved(nkndim)        !wave direction (same as wind direction)
        real waveov(nkndim)	!orbital velocity
        common /waveh/waveh, /wavep/wavep, /waved/waved, /waveov/waveov

! -------------------------------------------------------------
! sedtrans05 variables
! -------------------------------------------------------------

        DOUBLE PRECISION D		!WATER DEPTH (M)
        DOUBLE PRECISION UZ		!AMBIENT CURRENT AT HEIGHT Z ABOVE THE SEAFLOOR (M/S)
        DOUBLE PRECISION Z		!HEIGHT OF UZ ABOVE SEAFLOOR
        DOUBLE PRECISION CDIR		!DIRECTION OF AMBIENT CURRENT (DEGREES TRUE)
        DOUBLE PRECISION HT		!WAVE HEIGHT (M)
        DOUBLE PRECISION PER		!WAVE PERIOD (S)
        DOUBLE PRECISION WDIR		!WAVE PROPAGATION DIRECTION (DEGREES TRUE)
        DOUBLE PRECISION GD		!SEDIMENT GRAIN DIAMETER (M) 
        DOUBLE PRECISION BETA		!BED SLOPE (DEGREE)
        DOUBLE PRECISION TEM		!in degrees Celsius
        DOUBLE PRECISION SAL		!practical salinity scale (seawater=35)
        DOUBLE PRECISION TIMEDR		!Duration of call to Sedtrans (s)
        DOUBLE PRECISION Z0       	!BED ROUGHNESS LENGTH (M)

! cohesive
        DOUBLE PRECISION AULVA		!PERCENTAGE OF AREA COVERED BY THE ALGAE 'ULVA' (%)

!       -------------------------------------------------------------
!       Phisical parameter
!       -------------------------------------------------------------

        BETA  = 0.			! [degree]
        AULVA  = 0.                     ! [%]

!       -------------------------------------------------------------------
!       Starts loop on nodes
!       -------------------------------------------------------------------

        do k = 1,nkn

          D = 0.
          l = ilhkv(k)				!bottom layer
          m = l

!         -------------------------------------------------------------------
!         Initializes input arguments
!         -------------------------------------------------------------------

          HT = waveh(k)                         !wave height
          PER = wavep(k)			!wave period
          WDIR = waved(k)                       !wave direction

          do i = 1,l
           D = D + hdknv(i,k)			!gets total water depth
          end do
          DL = hdknv(l,k)			!gets depth of bottom layer

          if ( DL .lt. 0.20 .and. l .gt. 1) then
            m = l - 1
            DL = DL + hdknv(m,k)
          end if

          call getuv(m,k,u,v)                   !gets velocities u/v
          call getmd(u,v,UZ,CDIR)		!gets UZ, CDIR
          GD = gskm(k)				!gs50
          if (GD .gt.0.1) GD = 0.0025		!rocks
          call zvel(UZ,DL,GD,z0bk(k),Z)		!get Z
          call getts(m,k,temp,salt)             !gets temp and salt
          if (temp.eq.0. .and. salt.eq.0.) then
            temp = temref
            salt = salref
          end if
          TEM = temp
          SAL = salt

!         -------------------------------------------------------------------
!         Calculates sediment transport rate for each sediment class in node k
!         -------------------------------------------------------------------

          do is = 1,nscls
           scnd(is) = scn(l,k,is)
          end do

          do i = 1,nbcc
            sccd(i) = scc(l,k,i)
          end do

        call sedt05(it,k,D,DL,UZ,Z,CDIR,HT,PER,WDIR,GD,riph(k),ripl(k),
     @BETA,TEM,SAL,bedn(1,1,k),percbd(1,1,k),AULVA,TIMEDR,nscls,gs,
     @hzoff,scnd,sccd,sedx(1,k),sedy(1,k),sflx(1,k),ws,gdx(k),
     @gdy(k),l,eps,tao(k),Z0)

	  z0bk(k) = Z0

          do is = 1,nscls
            scn(l,k,is) = scnd(is)
          end do

          do i = 1,nbcc
            scc(l,k,i) = sccd(i)
          end do

!       -------------------------------------------------------------------
!       End of node loop
!       -------------------------------------------------------------------

        end do

        end

! ********************************************************************
! SUBROUTINE SEDT05

      subroutine sedt05(it,k,D,DL,UZ,Z,CDIR,HT,PER,WDIR,GD,RHINP,RLINP,
     @BETA,TEM,SAL,BEDCHA,percbd,AULVA,TIMEDR,nscls,gs,hzoff,scn,CONC,
     @sedx,sedy,sflx,ws,gdx,gdy,l,eps,tao,Z0)

        implicit none

	include 'param.h'
        include 'sed_param.h'

	integer it
! ------------ INPUT VARIABLES -----------------
        DOUBLE PRECISION D        	!WATER DEPTH (M)
        DOUBLE PRECISION DL		!DEPTH OF BOTTOM LAYER (M)
        DOUBLE PRECISION UZ       	!AMBIENT CURRENT AT HEIGHT Z ABOVE THE SEAFLOOR (M/S)
        DOUBLE PRECISION Z        	!HEIGHT OF UZ ABOVE SEAFLOOR
        DOUBLE PRECISION CDIR     	!DIRECTION OF AMBIENT CURRENT (DEGREES TRUE)
        DOUBLE PRECISION HT       	!WAVE HEIGHT (M)
        DOUBLE PRECISION PER      	!WAVE PERIOD (S)
        DOUBLE PRECISION WDIR     	!WAVE PROPAGATION DIRECTION (DEGREES TRUE)
        DOUBLE PRECISION GD       	!SEDIMENT GRAIN DIAMETER (M)
        DOUBLE PRECISION RHINP    	!INPUT RIPPLE HEIGHT (M)
        DOUBLE PRECISION RLINP    	!INPUT RIPPLE LENGTH (M)
        DOUBLE PRECISION BETA     	!BED SLOPE (DEGREE)
        DOUBLE PRECISION TEM      	!in degrees Celsius
        DOUBLE PRECISION SAL      	!practical salinity scale (seawater=35)
        DOUBLE PRECISION TIMEDR   	!Duration of call to Sedtrans (s)
        double precision percbd(nlbdim,nsdim)   !fraction of sediment [0,1]
        DOUBLE PRECISION BEDCHA(nlbdim,3) ! bed characteristics in 3 column table
                                        ! (1) depth below sediment surface (m)
                                        ! (2) critical erosion stress (Pa)
                                        ! (3) dry bulk density (kg/m**3) user input
        integer nscls			!number of grainsize classes
	real hzoff

! cohesive
        DOUBLE PRECISION AULVA    	!PERCENTAGE OF AREA COVERED BY THE ALGAE 'ULVA' (%)

! ------------ OUTPUT VARIABLES -----------------
        double precision sedx(nsdim)
        double precision sedy(nsdim)
        double precision sflx(nsdim)	!flux of suspend sediment [m3/m2]
        double precision ws(nsdim)      !SETTLING VELOCITY FOR NON-COHESIVE SEDIMENT (M/SEC)
        DOUBLE PRECISION CONC(NBCC) 	!COHESIVE SUSPENDED SEDIMENT CONCENTRATION (kg/m3)
        real tao			!wave-current shear stress

! ------------ LOCAL VARIABLES -----------------
        double precision bmix		!thickness of active layer
        DOUBLE PRECISION UB       	!MAXIMUM WAVE INDUCED ORBITAL VELOCITY AT THE BOTTOM (M/S)
        DOUBLE PRECISION FCW      	!BOTTOM (SKIN) FRICTION FACTOR
        DOUBLE PRECISION UA       	!CURRENT SPEED TO BE USED IN BOTTOM STRESS CALC. (M/SEC)
        DOUBLE PRECISION U100     	!CURRENT SPEED AT 1 M. ABOVE SEABED (M/SEC)
        DOUBLE PRECISION USTCWS   	!COMBINED SKIN-FRICTION SHEAR VELOCITY OF GM
        DOUBLE PRECISION USTCW    	!COMBINED TOTAL SHEAR VELOCITY OF GM
        DOUBLE PRECISION Z0       	!BED ROUGHNESS LENGTH (M)
        DOUBLE PRECISION RHOW     	!DENSITY OF FLUID (WATER)  (KG/M**3)
        DOUBLE PRECISION VISC   	!dynamic viscosity of the (sea)water (KG/M*SEC (OR N.S/M**2))
        DOUBLE PRECISION VISK     	!KINEMAMIC VISCOSITY OF THE FLUID (KG/M*SEC (OR N.S/M**2))
        DOUBLE PRECISION RHOS           !DENSITY OF SEDIMENT MINERAL(S)   (KG/M**3)
        DOUBLE PRECISION TCONC1		!total cohesive suspended sediment concentration (kg/m**3)
        DOUBLE PRECISION TCONC2  	!total SSC (sum of CONC) after deposition (kg/m**3)
        DOUBLE PRECISION DMASS		!deposited mass (kg)
        DOUBLE PRECISION TAU0		!effective skin friction shear stress (Pa)
        DOUBLE PRECISION EMASS		!eroded mass (kg/m2)
        DOUBLE PRECISION ZS       	!Height change in bed surface (m)
                                  	!(positive: erosion, negative: deposition)
        DOUBLE PRECISION AB	        !EXCURSION LENGTH OF BOTTOM WAVE ORBIT (M)
        DOUBLE PRECISION WL      	!WAVE LENGTH (M)
        DOUBLE PRECISION Z0C		!APPARENT BED ROUGHNESS LENGTH (M) 
        DOUBLE PRECISION PHIB		!ANGLE BETWEEN WAVE AND CURRENT DIRECTIONS (RADIANS)
        DOUBLE PRECISION PHI100		!ANGLE BETWEEN WAVE AND CURRENT DIRECTIONS AT 1 M. ABOVE SEABED
        DOUBLE PRECISION DELTACW	!HEIGHT OF THE WAVE-CURRENT BOUNDARY LAYER
        DOUBLE PRECISION USTCS		!CURRENT SKIN-FRICTION SHEAR VELOCITY OF GM
        DOUBLE PRECISION USTWS		!WAVE SKIN-FRICTION SHEAR VELOCITY OF GM
        DOUBLE PRECISION USTCWSB	!TRANSPORT-RELATED COMBINED SHEAR VELOCITY
        DOUBLE PRECISION USTC		!TOTAL CURRENT SHEAR VELOCITY OF GM
        DOUBLE PRECISION USTW		!TOTAL WAVE SHEAR VELOCITY OF GM
        DOUBLE PRECISION RPLCOEF	!RIPPLE COEFFICIENT FOR SHEAR VELOCITY CONVERSION
        DOUBLE PRECISION USTBF		!CRITICAL SHEAR VELOCITY FOR RIPPLE BREAKOFF

        double precision pers(nsdim)
        double precision scn(nsdim) 	!non-cohesive suspended sediment conc (kg/m3)
        double precision usb(nsdim)	!CRITICAL SHEAR VELOCITY FOR INITIATION OF BEDLOAD
        double precision uss(nsdim)	!CRITICAL SHEAR VELOCITY FOR INITIATION OF SUSPENDED
        double precision ust(nsdim)	!CRITICAL SHEAR VELOCITY FOR INITIATION OF SHEET FLOW
        double precision dxx(nsdim)	!DIMENSIONLESS GRAIN SIZE
        double precision gs(nsdim)      !SEDIMENT GRAIN DIAMETER (M)
        integer l,is,isd,k		!counters
        real gdx,gdy			!slope gradients
        real alph                       !slope effect
        real eps(0:nlvdim,nkndim,nsdim)	!vertical mixing coefficient
	real dmin
        real hdknv(nlvdim,nkndim)
        common /hdknv/hdknv

!       --------------------------------------------------
!       Initializes variables
!       --------------------------------------------------

        do is = 1,nscls
         usb(is) = 0.
         uss(is) = 0.
         ust(is) = 0
         dxx(is) = 0.
         sedx(is) = 0.
         sedy(is) = 0.
         sflx(is) = 0.
         pers(is) = percbd(1,is)
        end do

        TCONC1 = 0.

        ZS = 0.
        EMASS = 0.
        bmix = 0.
        RHOS = RHOSED
        TCONC2 = 0.
	dmin = 100.
        isd = 1

!       -------------------------------------------------------------------
!       Calculates slope effect
!       -------------------------------------------------------------------

        call bedslope(cdir,gdx,gdy,angrep,alph)

!       -------------------------------------------------------------------
!       Calculates water density and viscosity
!       -------------------------------------------------------------------

        call densvisc(tem,sal,rhow,visc)	!get density and viscosity
        VISK = visc/rhow 		        !get viscosity

!       -------------------------------------------------------------------
!       Calculates wave-induced bottom particle velocity and orbital diameter
!       -------------------------------------------------------------------

        CALL OSCIL(HT,PER,D,UB,AB,WL)

!       -------------------------------------------------------------------
!       Calculates threshold criteria for bedload, suspended load and sheet flow
!       -------------------------------------------------------------------

        do is = 1,nscls
         CALL THRESH(VISK,gs(is),RHOS,RHOW,ws(is),usb(is),uss(is),
     @ust(is),dxx(is))

!        -------------------------------------------------------------------
!        Updates threshold criteria in function of consolidation and mud fraction
!        -------------------------------------------------------------------

         if ( gs(1) .lt. limcoh .and. pers(1).lt. KCOES .and.
     @pers(1).gt.0.1 ) then     !Van Rijn, 1993
           usb(is) = usb(is)*(pers(1)*100.)**0.5
           uss(is) = uss(is)*(pers(1)*100.)**0.5
           ust(is) = ust(is)*(pers(1)*100.)**0.5
         end if

!        -------------------------------------------------------------------
!        Corrects bedload threshold criteria if function of bed slope
!        -------------------------------------------------------------------

         usb(is) = usb(is)*sqrt(alph)
         uss(is) = max(usb(is),uss(is))

         if(abs(gs(is)-GD).lt.dmin) then
           dmin = abs(gs(is)-GD)
           isd = is
         end if

        end do

!       -------------------------------------------------------------------
!       Calculates bottom stresses for d50
!       -------------------------------------------------------------------

        CALL FRICFAC(D,UZ,Z,CDIR,UB,PER,WDIR,AB,GD,VISK,RHOW,RHOS,
     @Z0C,PHIB,PHI100,DELTACW,USTCS,USTWS,USTCWSB,USTC,USTW,RPLCOEF,
     @USTBF,usb(isd),ust(isd),Z0,FCW,UA,U100,USTCWS,USTCW,RHINP,RLINP)

        if(ustcw.eq.0.) ustcw=max(ustc,ustw)
        if(ustcws.eq.0.) ustcws=max(ustcs,ustws)
        tao = rhow * ustcws**2

!       -------------------------------------------------------------------
!       Calculates vertical mixing coefficient
!       -------------------------------------------------------------------

        do is = 1,nscls
         call vmixcoef(l,hdknv(1,k),D,ws(is),ustc,ub,deltacw,ht,per,
     @                 dxx(is),eps(1,k,is))
        end do

!       -------------------------------------------------------------------
!       Returns if Depth < hzoff
!       -------------------------------------------------------------------

        if (D.LT.hzoff) RETURN

!       -------------------------------------------------------------------
!       Checks if surficial layer is smaller than active layer (bmix = Z0)
!       -------------------------------------------------------------------

        bmix = Z0

        if(BEDCHA(2,1).lt.bmix .and. (RHOW*USTCWS**2 .gt. BEDCHA(1,2)
     @ .or. USTCWS. gt. usb(isd))) then
          call bedmerge(nscls,bmix,BEDCHA,percbd)
          do is = 1,nscls
           pers(is) = percbd(1,is)
          end do
        end if

!       -------------------------------------------------------------------
!       Calculates sediment transport rate for each sediment class in node k
!       -------------------------------------------------------------------

        if( gs(1) .lt. limcoh .and. pers(1) .ge. KCOES ) then

!         -------------------------------------------------------------------
!         COHESIVE SEDIMENT
!         -------------------------------------------------------------------

          call cohse(k,l,BEDCHA,TIMEDR,USTCWS,USTCW,DL,RHOW,VISK,
     @AULVA,UB,CDIR,Z0,CONC,nscls,scn,ws,usb,uss,ust,gs,dxx,
     @Z0C,DELTACW,USTCS,USTWS,USTCWSB,USTC,USTW,RPLCOEF,USTBF,
     @pers,eps,RHOS,TAU0,EMASS,sflx,TCONC1)

          scn(1) = TCONC1
        else

!         -------------------------------------------------------------------
!         NON-COHESIVE SEDIMENT
!         -------------------------------------------------------------------

          call nonco(it,k,l,BEDCHA,D,DL,UA,UB,U100,HT,PER,CDIR,
     @RHOW,USTCWS,USTCW,usb,uss,ust,Z0,BETA,VISK,RHOS,FCW,WDIR,TIMEDR,
     @Z0C,DELTACW,USTCS,USTWS,USTCWSB,USTC,USTW,RPLCOEF,USTBF,
     @nscls,pers,gs,dxx,ws,eps,scn,sedx,sedy,sflx,EMASS)

        end if

!       -------------------------------------------------------------------
!       Deposition for cohesive suspended sediments
!       -------------------------------------------------------------------
 
        if(TCONC1.gt.0.) then 
          call depos(DL,TIMEDR,NBCC,CONC,TAU0,RHOW,TCONC1,TCONC2,VISK,
     @DMASS)

          ZS = DMASS*2./(RHOMUD+BEDCHA(1,3))
          sflx(1) = sflx(1) + ZS
          scn(1) = TCONC2
        end if

!       -------------------------------------------------------------------
!       Puts cohesive eroded sediment in suspension
!       -------------------------------------------------------------------

        if(EMASS .gt. 0.) then
         call eros2(BEDCHA(1,2),DL,RHOW,TAU0,NBCC,CONC,EMASS,TCONC2)

         scn(1) = TCONC2
        end if

        end

! ********************************************************************
! SUBROUTINE COHSE
! This subroutine computes the sediment transport for cohesive
! sediments

        subroutine cohse(k,l,BEDCHA,TIMEDR,USTCWS,USTCW,D,RHOW,VISK,
     @AULVA,UB,CDIR,Z0,CONC,nscls,scn,ws,usb,uss,ust,gs,dxx,
     @Z0C,DELTACW,USTCS,USTWS,USTCWSB,USTC,USTW,RPLCOEF,USTBF,
     @pers,eps,RHOS,TAU0,EMASS,sflux,TCONC1)

        implicit none

	include 'param.h'
        include 'sed_param.h'

! ------------ INPUT VARIABLES -----------------
        DOUBLE PRECISION BEDCHA(nlbdim,3) ! bed characteristics in 3 column table
                                  	! (1) depth below sediment surface (m)
                                  	! (2) critical erosion stress (Pa)
                                  	! (3) dry bulk density (kg/m**3)          user input
        DOUBLE PRECISION TIMEDR   	!Duration of call to Sedtrans (s)
        DOUBLE PRECISION USTCWS   	!Combined skin friction stress            FRICFRAC
        DOUBLE PRECISION D        	!WATER DEPTH (M)
        DOUBLE PRECISION RHOW     	!DENSITY OF FLUID (WATER)  (KG/M**3)
        DOUBLE PRECISION VISK     	!KINEMAMIC VISCOSITY OF THE FLUID (KG/M*SEC (OR N.S/M**2))
        DOUBLE PRECISION AULVA    	!PERCENTAGE OF AREA COVERED BY THE ALGAE 'ULVA' (%)
        DOUBLE PRECISION UB       	!MAXIMUM WAVE INDUCED ORBITAL VELOCITY AT THE BOTTOM (M/S)
        DOUBLE PRECISION USTCW    	!COMBINED TOTAL SHEAR VELOCITY OF GM
        DOUBLE PRECISION CDIR     	!DIRECTION OF AMBIENT CURRENT (DEGREES TRUE)
        DOUBLE PRECISION Z0       	!BED ROUGHNESS LENGTH (M)
        DOUBLE PRECISION CONC(NBCC) 	!COHESIVE SUSPENDED SEDIMENT CONCENTRATION (kg/m3)

        integer nscls			!number of grainsize classes
        double precision scn(nsdim) 	!cohesive suspended sediment conc (kg/m3)
        double precision pers(nsdim)	!fraction of each class
        double precision gs(nsdim)      !SEDIMENT GRAIN DIAMETER (M)
        double precision usb(nsdim)     !CRITICAL SHEAR VELOCITY FOR INITIATION OF BEDLOAD
        double precision uss(nsdim)     !CRITICAL SHEAR VELOCITY FOR INITIATION OF SUSPENDED
        double precision ust(nsdim)     !CRITICAL SHEAR VELOCITY FOR INITIATION OF SHEET FLOW
        double precision dxx(nsdim)     !DIMENSIONLESS GRAIN SIZE
        double precision ws(nsdim)      !SETTLING VELOCITY FOR NON-COHESIVE SEDIMENT (M/SEC)
        real eps(0:nlvdim,nkndim,nsdim) !vertical mixing coefficient

        DOUBLE PRECISION Z0C      !APPARENT BED ROUGHNESS LENGTH (M)
        DOUBLE PRECISION DELTACW  !HEIGHT OF THE WAVE-CURRENT BOUNDARY LAYER
        DOUBLE PRECISION USTCS    !CURRENT SKIN-FRICTION SHEAR VELOCITY OF GM
        DOUBLE PRECISION USTWS    !WAVE SKIN-FRICTION SHEAR VELOCITY OF GM
        DOUBLE PRECISION USTCWSB  !TRANSPORT-RELATED COMBINED SHEAR VELOCITY
        DOUBLE PRECISION USTC     !TOTAL CURRENT SHEAR VELOCITY OF GM
        DOUBLE PRECISION USTW     !TOTAL WAVE SHEAR VELOCITY OF GM
        DOUBLE PRECISION RPLCOEF  !RIPPLE COEFFICIENT FOR SHEAR VELOCITY CONVERSION
        DOUBLE PRECISION USTBF    !CRITICAL SHEAR VELOCITY FOR RIPPLE BREAKOFF

! ------------ OUTPUT VARIABLES -----------------
        DOUBLE PRECISION TAU0		!effective skin friction shear stress (Pa)
        DOUBLE PRECISION EMASS		!eroded mass (kg/m2)
        double precision sflux(nsdim)	!flux of suspended sediment [m3/m2]
        double precision sandf(nsdim)	!flux of suspended sediment [m3/m2]

! ------------ LOCAL VARIABLES -----------------
        INTEGER NBED         		!number of rows in BEDCHA that are used
        DOUBLE PRECISION ZS       	!Height change in bed surface (m)
                                  	!(positive: erosion, negative: deposition)
        DOUBLE PRECISION TAOS		!solid transmitted stress due to Ulva (Pa)
        DOUBLE PRECISION TCONC1		!total SSC (sum of CONC) before deposition (kg/m**3)
        DOUBLE PRECISION RHOS		!dry bulk density at depth ZS (kg/m**3)
        DOUBLE PRECISION TAUCE		!critical erosion stress at depth ZS (Pa)
        DOUBLE PRECISION MAXEMASS	!Maximum erodable sediment mass to avoid problem with drag reduction
        DOUBLE PRECISION C0A     	!DEPTH AVERAGED REFERENCE CONCENTRATION AT Z0 (KG/M^3)
        DOUBLE PRECISION QS		!SUSPENDED SEDIMENT TRANSPORT RATE (KG/M/S)
        DOUBLE PRECISION QSDIR		!DIRECTION OF SUSPENDED SEDIMENT TRANSPORT (DEGREE)
        integer l,i,is,k		!counters
	double precision massa
       
!       --------------------------------------------------
!       Initializes variables
!       --------------------------------------------------
 
        emass = 0.

!       -------------------------------------------------------------------
!       Computes the total SSC before deposition and erosion
!       -------------------------------------------------------------------

        TCONC1 = 0.
        DO 210, I=1,NBCC
          TCONC1 = TCONC1 + CONC(I)
210     CONTINUE

!       -------------------------------------------------------------------
!       Computes TAU0, the effective skin friction stress acting on the bed
!       drag reduction and solid transmitted stress due to ulva
!       -------------------------------------------------------------------

        TAU0=RHOW*USTCWS**2

        CALL DRAGRED(TAU0,TCONC1,D,MAXEMASS)
  
        CALL SOLULVA(TAU0,AULVA,TAOS)

!       -------------------------------------------------------------------
!       EROSION: Test if effective bed shear stress (TAU0) is larger or 
!                equal to surface critical erosion stress (BEDCHA(1,2)
!       -------------------------------------------------------------------

        if( TAU0.GT.BEDCHA(1,2) .OR. ( TAU0.EQ.BEDCHA(1,2).AND.
     &    ( TAU0.GE.BEDCHA(2,2) ) ) ) THEN

          DO 10, I=2,nlbdim
           IF (BEDCHA(I,1).EQ.0.0) GOTO 20
10        CONTINUE
20        NBED=I-1

          CALL EROS1(NBED,nlbdim,BEDCHA,TAU0,MAXEMASS,EMASS,ZS,TIMEDR,
     @TAUCE)

          do is = 1,nscls
            sflux(is) = -ZS * pers(is)
            massa = emass * pers(is)
            scn(is) = scn(is) + massa/D
          end do
	
          emass = emass * pers(1)

        end if

!       -------------------------------------------------------------------
!       DEPOSITION: check if non-cohesive sediments deposit
!       -------------------------------------------------------------------        

        do is = 2,nscls

          sandf(is) = 0.
          if( scn(is) .gt. 0. .and. gs(is).gt.0.01) then

          CALL PROFL(D,RHOW,ws(is),UB,CDIR,USTCWS,USTCW,usb(is),
     @uss(is),ust(is),Z0C,DELTACW,USTCS,USTWS,USTCWSB,USTC,USTW,
     @RPLCOEF,USTBF,Z0,gs(is),dxx(is),RHOS,QS,QSDIR,C0A)

            !C0A = C0A * pers(is)
            C0A = MIN(C0A,scn(is))

            call suspco(k,bedcha,timedr,ws(is),Z0,C0A,eps(l,k,is),
     @pers(is),scn(is),sandf(is))

            sflux(is) = sflux(is) + sandf(is)

          end if

        end do

        end 

! ********************************************************************
! SUBROUTINE NONCO
! This subroutine computes the sediment transport for non-cohesive
! sediments

        subroutine nonco(it,k,l,BEDCHA,D,DL,UA,UB,U100,HT,PER,CDIR,RHOW,
     @USTCWS,USTCW,usb,uss,ust,Z0,BETA,VISK,RHOS,FCW,WDIR,TIMEDR,
     @Z0C,DELTACW,USTCS,USTWS,USTCWSB,USTC,USTW,RPLCOEF,USTBF,
     @nscls,pers,gs,dxx,ws,eps,scn,sedx,sedy,sflux,EMASS)

        implicit none

	include 'param.h'
        include 'sed_param.h'

	integer it
! ------------ INPUT VARIABLES -----------------
        DOUBLE PRECISION BEDCHA(nlbdim,3) ! bed characteristics in 3 column table
                                        ! (1) depth below sediment surface (m)
                                        ! (2) critical erosion stress (Pa)
                                        ! (3) dry bulk density (kg/m**3)          user input
        DOUBLE PRECISION D        !WATER DEPTH (M)
        DOUBLE PRECISION DL       !DEPTH OF BOTTOM LAYER (M)
        DOUBLE PRECISION UA       !CURRENT SPEED TO BE USED IN BOTTOM STRESS CALC. (M/SEC)
        DOUBLE PRECISION UB       !MAXIMUM WAVE INDUCED ORBITAL VELOCITY AT THE BOTTOM (M/S)
        DOUBLE PRECISION U100     !CURRENT SPEED AT 1 M. ABOVE SEABED (M/SEC)
        DOUBLE PRECISION HT       !WAVE HEIGHT (M)
        DOUBLE PRECISION PER      !WAVE PERIOD (S)
        DOUBLE PRECISION CDIR     !DIRECTION OF AMBIENT CURRENT (DEGREES TRUE)
        DOUBLE PRECISION RHOW     !DENSITY OF FLUID (WATER)  (KG/M**3)
        DOUBLE PRECISION USTCWS   !COMBINED SKIN-FRICTION SHEAR VELOCITY OF GM
        DOUBLE PRECISION USTCW    !COMBINED TOTAL SHEAR VELOCITY OF GM
        DOUBLE PRECISION Z0       !BED ROUGHNESS LENGTH (M)
        DOUBLE PRECISION BETA     !BED SLOPE (DEGREE)
        DOUBLE PRECISION VISK     !KINEMAMIC VISCOSITY OF THE FLUID (KG/M*SEC (OR N.S/M* *2))
        DOUBLE PRECISION RHOS     !DENSITY OF SEDIMENT MINERAL(S)   (KG/M**3)
        DOUBLE PRECISION FCW      !BOTTOM (SKIN) FRICTION FACTOR
        DOUBLE PRECISION WDIR     !WAVE PROPAGATION DIRECTION (DEGREES TRUE)
        DOUBLE PRECISION TIMEDR   !Duration of call to Sedtrans (s)
        DOUBLE PRECISION Z0C      !APPARENT BED ROUGHNESS LENGTH (M)
        DOUBLE PRECISION DELTACW  !HEIGHT OF THE WAVE-CURRENT BOUNDARY LAYER
        DOUBLE PRECISION USTCS    !CURRENT SKIN-FRICTION SHEAR VELOCITY OF GM
        DOUBLE PRECISION USTWS    !WAVE SKIN-FRICTION SHEAR VELOCITY OF GM
        DOUBLE PRECISION USTCWSB  !TRANSPORT-RELATED COMBINED SHEAR VELOCITY
        DOUBLE PRECISION USTC     !TOTAL CURRENT SHEAR VELOCITY OF GM
        DOUBLE PRECISION USTW     !TOTAL WAVE SHEAR VELOCITY OF GM
        DOUBLE PRECISION RPLCOEF  !RIPPLE COEFFICIENT FOR SHEAR VELOCITY CONVERSION
        DOUBLE PRECISION USTBF    !CRITICAL SHEAR VELOCITY FOR RIPPLE BREAKOFF
        DOUBLE PRECISION TB1      !TIME AT WHICH BEDLOAD TRANSPORT CEASES (SEC) 
        DOUBLE PRECISION TB2      !TIME AT WHICH BEDLOAD TRANSPORT RECOMMENCES (SEC)
        DOUBLE PRECISION TS1      !TIME AT WHICH SUSPENDED LOAD TRANSPORT CEASES (SEC)
        DOUBLE PRECISION PERBED   !PERCENTAGE OF TIME SPENT IN ONLY BEDLOAD TRANSPORT PHASE
        DOUBLE PRECISION PERSUSP  !PERCENTAGE OF TIME SPENT IN SUSPENDED LOAD TRANSPORT PHASE
        DOUBLE PRECISION PHIB     !ANGLE BETWEEN WAVE AND CURRENT DIRECTIONS (RADIANS)
        DOUBLE PRECISION PHI100   !ANGLE BETWEEN WAVE AND CURRENT DIRECTIONS AT 1 M. ABOVE SEABED

        integer nscls			!number of grainsize classes
        double precision pers(nsdim)	!fraction of sediment [0,1]
        double precision gs(nsdim)      !SEDIMENT GRAIN DIAMETER (M)
        double precision usb(nsdim)     !CRITICAL SHEAR VELOCITY FOR INITIATION OF BEDLOAD
        double precision uss(nsdim)     !CRITICAL SHEAR VELOCITY FOR INITIATION OF SUSPENDED
        double precision ust(nsdim)     !CRITICAL SHEAR VELOCITY FOR INITIATION OF SHEET FLOW
        double precision dxx(nsdim)     !DIMENSIONLESS GRAIN SIZE
        double precision ws(nsdim)      !SETTLING VELOCITY FOR NON-COHESIVE SEDIMENT (M/SEC)
        double precision scn(nsdim) 	!non-cohesive suspended sediment conc (kg/m3)
        real eps(0:nlvdim,nkndim,nsdim) !vertical mixing coefficient

! ------------ OUTPUT VARIABLES -----------------
        double precision sedx(nsdim)	!x bedload component [m3/m]
        double precision sedy(nsdim)	!y bedload component [m3/m]
        double precision sflux(nsdim)	!flux of suspend sediment [m3/m2]
        DOUBLE PRECISION EMASS		!eroded mass [kg/m2]

! ------------ LOCAL VARIABLES -----------------
        DOUBLE PRECISION QS       !SUSPENDED SEDIMENT TRANSPORT RATE (KG/M/S)
        DOUBLE PRECISION QSDIR    !DIRECTION OF SUSPENDED SEDIMENT TRANSPORT (DEGREE)
        DOUBLE PRECISION SED      !TIME-AVERAGED NET SEDIMENT TRANSPORT AS VOLUME (M**3/S/M)
        DOUBLE PRECISION SEDM     !TIME-AVERAGED NET SEDIMENT TRANSPORT AS MASS (KG/S/M)
        DOUBLE PRECISION SEDDIR   !DIRECTION OF NET SEDIMENT TRANSPORT (AZIMUTH,DEGREES)
        DOUBLE PRECISION C0A      !DEPTH AVERAGED REFERENCE CONCENTRATION AT Z0 (KG/M^3)
        integer l,is,k		  !counters
	real totsedm

	totsedm = 0.
        do is = 1,nscls

          C0A = 0.

          if (pers(is).gt.0.) then          
	  if (gs(is) .gt. 0.01) goto 120

!          -------------------------------------------------------------------
!          Calculates the duration of the different sediment transport phases
!          -------------------------------------------------------------------
           CALL TIMING(RHOW,UA,UB,PER,U100,usb(is),uss(is),USTCWS,
     @PHIB,USTCS,USTWS,RPLCOEF,TB1,TB2,TS1,PERBED,PERSUSP)

!          -------------------------------------------------------------------
!          Calculates sediment transport rate and direction
!          -------------------------------------------------------------------
           CALL TRANSPO(D,UA,UB,U100,PER,gs(is),BETA,RHOS,RHOW,
     @usb(is),ust(is),PHIB,PHI100,USTCS,USTWS,USTCWSB,RPLCOEF,TB1,
     @TB2,TS1,PERBED,PERSUSP,USTBF,FCW,USTCWS,HT,CDIR,WDIR,IOPT,
     @dxx(is),ws(is),SED,SEDM,SEDDIR)

           if(SEDDIR.le.90.and.SEDDIR.ge.0)   SEDDIR=90.-SEDDIR
           if(SEDDIR.le.360.and.SEDDIR.ge.90) SEDDIR=450.-SEDDIR
           SEDDIR = SEDDIR / (45. / atan (1.))                 !rad

           SEDM = SEDM * pers(is)
	   if(k.eq.6710) totsedm = totsedm + sedm
           SEDM = 2. * SEDM / (BEDCHA(1,3)+BEDCHA(2,3)) !m3/ms

           sedx(is) = SEDM * cos(SEDDIR)
           sedy(is) = SEDM * sin(SEDDIR)

          end if

120       continue

!         -------------------------------------------------------------------
!         Calculates velocity profile, suspended sediment concentration profile
!         -------------------------------------------------------------------

	  if (IOPT.eq.1 .or. IOPT.eq.3 .and. gs(is) .gt. 0.01) goto 130
          CALL PROFL(DL,RHOW,ws(is),UB,CDIR,USTCWS,USTCW,usb(is),
     @uss(is),ust(is),Z0C,DELTACW,USTCS,USTWS,USTCWSB,USTC,USTW,
     @RPLCOEF,USTBF,Z0,gs(is),dxx(is),RHOS,QS,QSDIR,C0A)
!          C0A = C0A * pers(is)
130       continue

          call suspco(k,bedcha,timedr,ws(is),Z0,C0A,eps(l,k,is),
     @pers(is),scn(is),sflux(is))

          if(gs(is).lt.limcoh.and.sflux(is).lt.0.) 
     @EMASS=-sflux(is)*bedcha(1,3)

        end do
	if(k.eq.6710) write(92,*)it,totsedm

        end

! ********************************************************************
! SUBROUTINE SUSPCO
! This subroutine computes the suspended concentration for non-cohesive
! sediment due to the sediment flux from the bottom. If edr < 0
! (erosion) use an explicit scheme, otherwise if edr > 0 (deposition)
! use an implicit scheme.

        subroutine suspco(k,bedcha,ddt,fal,z0,ceq,eps,pers,cn,edr)

        implicit none

        include 'param.h'
        include 'sed_param.h'

        integer ilhkv(nkndim)           !number of node level
        common /ilhkv/ilhkv
        real hdknv(nlvdim,nkndim)
        common /hdknv/hdknv
        real hdkov(nlvdim,nkndim)
        common /hdkov/hdkov

! --- input variables
        double precision bedcha(nlbdim,3) ! bed characteristics in 3 column table
                                        ! (1) depth below sediment surface (m)
                                        ! (2) critical erosion stress (Pa)
                                        ! (3) dry bulk density (kg/m**3)          user input
        double precision ddt            !time step [s]
        double precision fal            !settling velocity [m/s]
        double precision z0		!bed roughness length (m)
        real eps			!vertical mixing coefficient
        double precision ceq            !equilibrium concentration [kg/m3]
        double precision cn             !new susp concetration [kg/m3]
	double precision pers		!percent

! --- output variables
        double precision edr            !deposition - erosion [m3/m2]

! --- local variables
        double precision co             !old susp concetration [kg/m3]
        double precision cne            !susp conc after erosion [kg/m3]
        real hnew                       !new depth [m]
        real hold                       !old depth [m]

        double precision dz,dvz
        integer k,l

!       -------------------------------------------------------------------
!       Initializes valieables
!       -------------------------------------------------------------------

        edr = 0.
        co = cn
        l = ilhkv(k)
        hnew = hdknv(l,k)
        hold = hdkov(l,k)
        dz = hnew/2. - z0
        dvz = eps/dz

!       -------------------------------------------------------------------
!       Source term - erosion - esplicit scheme
!       -------------------------------------------------------------------

        cne = (co*hold + ceq*dvz*ddt) / hnew

!       -------------------------------------------------------------------
!       Sink term - deposition - implicit scheme
!       -------------------------------------------------------------------

        cn = cne*hold / (hnew + (dvz+fal)*ddt)

        if(pers.eq.0.) cn = min(cn,co)

!       -------------------------------------------------------------------
!       Computes erosion - deposition and convert to m3/m2
!       -------------------------------------------------------------------

        if(cn.lt.10e-10)cn = 0.

        edr = (co - cn) * hnew

        if(edr.gt.0.) then
          edr = edr / (RHOSED*(1-POROS))
        elseif(edr.lt.0.) then
          edr = 2. * edr / (bedcha(1,3)+bedcha(2,3))
        end if

!       -------------------------------------------------------------------
!       If other classes present, limit erosion
!       -------------------------------------------------------------------

        if (-edr .gt. bedcha(2,1)*pers .and. pers.lt.0.998) then
          edr = - bedcha(2,1)*pers
          cn = (co*hold - edr*bedcha(1,3)) / hnew
        end if

        end

! ********************************************************************
! SUBROUTINE BEDLOAD
! This subroutine computes the bedload sediment transport using the 
! sediment continuity equation

        subroutine bedload(nscls,sedx,sedy,dt,bflx)

        implicit none

        include 'param.h'
        include 'sed_param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	include 'ev.h'
        integer nen3v(3,1)        		!node number
        common /nen3v/nen3v
        integer ilhkv(nkndim)			!number of element and node level
        common /ilhkv/ilhkv

! --- input variables
        integer nscls				!number grainsize classes
        double precision sedx(nsdim,nkndim)	!bedload transport in x direction [m3/ms]
        double precision sedy(nsdim,nkndim)	!bedload transport in y direction [m3/ms]
	real dt					!time step

! --- output variables
        double precision bflx(nsdim,nkndim)	!flux of bedload sediment [m3/m2]

! --- local variables
        integer k,ie,ii,is,l
        double precision bflux
        double precision sexe(nsdim)           	!element transport in x direction
        double precision seye(nsdim)           	!element transport in x direction
        real b,c           		        !x and y derivated form function [1/m]
        real area,areanode			!area of node
        double precision ss(nsdim),sm(nsdim)

!       -------------------------------------------------------------------
!       Initializes local variables
!       -------------------------------------------------------------------

        do k = 1,nkn
         do is = 1,nscls
           bflx(is,k) = 0.
         end do
        end do

        do ie = 1,nel

          area = 4.*ev(10,ie)                    !area of element/3

          do is = 1,nscls
            sm(is) = 0.
            ss(is) = 0.
          end do

!         -------------------------------------------------------------------
!         Converts node transport value to element value
!         -------------------------------------------------------------------

          do ii=1,3
            k = nen3v(ii,ie)
            do is = 1,nscls
              sm(is) = sm(is) + sedx(is,k)
              ss(is) = ss(is) + sedy(is,k)
            end do
          end do

          do is = 1,nscls
            sexe(is) = sm(is) / 3.
            seye(is) = ss(is) / 3.
          end do

!         -------------------------------------------------------------------
!         Loop over elements vertex
!         -------------------------------------------------------------------

          do ii = 1,3
            k = nen3v(ii,ie)
            b = ev(3+ii,ie)                                   !1/m
            c = ev(6+ii,ie)                                   !1/m

            do is = 1,nscls
              bflux = (dt*area) * (b*sexe(is) + c*seye(is))   !m3/m2
              bflx(is,k) = bflx(is,k) + bflux                 !m3/m2
            end do
          end do

        end do

!       -------------------------------------------------------------------
!       Loop over area node
!       -------------------------------------------------------------------

          do k = 1,nkn
            l = ilhkv(k)                          		!bottom layer
	    area = areanode(l,k)
            do is = 1,nscls
              bflx(is,k) = bflx(is,k) / area
            end do
          end do
         
        end

! ********************************************************************
! SUBROUTINE BEDSLOPE
! Compute the bed slope on nodes in the direction of the transport (current)
! direction and slope effect (alph). When bed slope equal 0 ---> alph = 1

        subroutine bedslope(cdir,gdx,gdy,beta,alph)

        implicit none

! --- input variables
        double precision cdir           !current direction [degree]
        real gdx,gdy			!depth gradient in x and y

! --- output variables
        real alph                       !slope effect

! --- local variables
        real ksl                        !node slope angle [radian]
        real asp	        	!node upslope aspect angle [radian]
        real phi                        !angle between asp and slope
        real dir			!current direction [radians]
        real pi
        real beta                       !angle of repose
	real bmin			!min
        real degtorad

!       -------------------------------------------------------------------
!       Initializes variables
!       -------------------------------------------------------------------

        alph = 1.
        pi = 2.*asin(1.)
        degtorad = (45./atan (1.))
        dir = cdir / degtorad   !from degree to rad
        bmin = 90. / degtorad 

!       -------------------------------------------------------------------
!       Computes slope angle and slope aspect
!       -------------------------------------------------------------------

        ksl = atan(sqrt(gdx**2 + gdy**2))
        asp = atan(gdx/gdy)
        if (gdy.eq.0.) asp = 0.
        if (gdy.gt.0.) asp = pi + asp
        if (asp.lt.0.) asp = 2.*pi - asp

!       -------------------------------------------------------------------
!       Computes bed slope effect (Soulsby,97, eq.80a)
!       -------------------------------------------------------------------

        phi = abs(asp - dir)
        phi = acos(cos(phi))

        if(ksl.ge.beta) then	!avalanche occur
          if (phi.gt.bmin) then
           alph = 1.7
          else
           alph = 0.2
          end if
          return
        end if

        alph = (cos(phi)*sin(ksl) + sqrt((cos(ksl))**2.*(tan(beta))**2.
     @ - (sin(phi))**2.*(sin(ksl))**2.)) / tan(beta)

        alph = max(alph,0.2)

        end

! ********************************************************************
! SUBROUTINE BEDMAN
! This subroutine manages the bed composition in function of what is eroded
! or deposited

        subroutine bedman(gs,nscls,timedr,bflux,sflux,gdx,
     @gdy,tao,bh,gskm,bdh)

        implicit none

        include 'param.h'
        include 'sed_param.h'

        integer nscls				!number of grainsize class
        double precision gs(nsdim)		!grainsize class
        double precision timedr                	!time step [s]
        double precision bflux(nsdim,nkndim)	!bedload sediment contribution [m3/m2]
        double precision sflux(nsdim,nkndim)	!suspended sediment contribution [m3/m2]
        real gdx(nkndim),gdy(nkndim)		!slope gradients
        real bh(nkndim)				!bed elevation
        real gskm(nkndim)			!average sediment grainsize in node k
        real tao(nkndim)			!combined wave-current shear stress

        double precision percbd(nlbdim,nsdim,nkndim)        !fraction of sediment [0,1]
        double precision bedn(nlbdim,3,nkndim)  !bed characteristics in 3 column table
                                                ! (1) depth below sediment surface (m)
                                                ! (2) critical erosion stress (Pa)
                                                ! (3) dry bulk density (kg/m**3)
        common /percbd/ percbd
        common /bedn/ bedn

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	double precision bdh(nkndim)		!total elevation change [>0depo,<0ero]
        double precision dzco			!elevation variation due to compaction [m]
        double precision tau0 			!effective skin friction shear stress (Pa)
        double precision flx 			!bedload + suspended flux [m]
        real gsa
        double precision gsmax
        double precision gsmin
        integer k,is

!       -------------------------------------------------------------------
!       Initializes valiables
!       -------------------------------------------------------------------

        gsmax = 0.
        gsmin = 100.

        do is = 1,nscls
          gsmax = max(gsmax,gs(is))*1.01
          gsmin = min(gsmin,gs(is))*0.99
        end do

!       -------------------------------------------------------------------
!       Loop over nodes --> only bottom layer
!       -------------------------------------------------------------------

        do k = 1,nkn
         dzco = 0.
         bdh(k) = 0.
	 tau0 = tao(k)

!        -------------------------------------------------------------------
!        Rearrange bed characteristics
!        -------------------------------------------------------------------

         do is = 1,nscls

          flx = bflux(is,k) + sflux(is,k)

          call checkbed(k,is,nscls,gs(1),bedn(1,1,k),percbd(1,1,k),
     +			tau0,flx)

          bdh(k) = bdh(k) + flx

         end do

!        -------------------------------------------------------------------
!        Calculates compaction
!        -------------------------------------------------------------------

         if(DOCOMPACT.ne.0)call compact(bedn(1,1,k),nlbdim,timedr,dzco)

!        -------------------------------------------------------------------
!        Computes bottom thikness variation
!        -------------------------------------------------------------------

         bdh(k) = bdh(k) - dzco
         bdh(k) = bdh(k) * MORPHO

         bh(k) = bh(k) + bdh(k) 

!        -------------------------------------------------------------------
!        Updates average sediment grainsize (D50)
!        -------------------------------------------------------------------

         gsa = 0.
         do is = 1,nscls
          gsa = gsa + gs(is)*percbd(1,is,k)
         end do
         gskm(k) = gsa

!       -------------------------------------------------------------------
!       End of node loop
!       -------------------------------------------------------------------

        end do

	do k = 1,nkn
         if(gskm(k).gt.gsmax.or.gskm(k).lt.gsmin) go to 120
	end do

        return

 120    continue
        write(6,*) 'Error in computing the average grainsize'
        write(6,*) 'gsaver:',gskm(k),'node:',k
        do is = 1,nscls
          write(6,*)'classnumber:',is,'grainsize:',gs(is)
        enddo
        stop 'error stop : gsaver'

        end

! ********************************************************************
! SUBROUTINE BEDMERGE
! If upper layer is smaller than the active layer --> create a new layer
! merging layer 1 and layer 2

        subroutine bedmerge(nscls,bmix,BEDCHA,percbd)

        implicit none

        include 'sed_param.h'

        double precision percbd(nlbdim,nsdim)	!fraction of sediment [0,1]
        double precision BEDCHA(nlbdim,3)         !bed characteristics in 3 column table
                                                ! (1) depth below sediment surface (m)
                                                ! (2) critical erosion stress (Pa)
                                                ! (3) dry bulk density (kg/m**3)
        integer nlbd				!number of bed layer in k
        integer nscls                           !number of grainsize class
        double precision bmix                   !thickness of active layer [m]
        real ptot
        integer ib,is
        double precision lin,x1,x2,y1,y2,x	!statement function which interpolate 
						!linearly the value Y of (X,Y), located 
						!on the line between (X1,Y1) and (X2,Y2)

        lin(x1,x2,y1,y2,x) = y1 + (x-x1)/(x2-x1)*(y2-y1)

!       -------------------------------------------------------------------
!       Finds the number of rows that are used. This is indicated by 
!       z=0 in the row after the last used row.
!       -------------------------------------------------------------------

        do 10, ib=2,nlbdim
          if (BEDCHA(ib,1).eq.0.0) goto 20
10      continue
20      nlbd=ib-1

!       -------------------------------------------------------------------
!       If level 3 is lower than bmix merge layer 1 and 2
!       -------------------------------------------------------------------

        if (BEDCHA(3,1).lt.bmix) then
125       nlbd = nlbd - 1

          if ( nlbd .lt. 3 ) then		!create the third layer
            nlbd = nlbd + 1
            do is = 1,nscls
              percbd(nlbd+1,is) = percbd(nlbd,is)
            end do
            BEDCHA(nlbd+1,1) = 2.*BEDCHA(nlbd,1) !- BEDCHA(nlbd-1,1)
            BEDCHA(nlbd+1,2) = BEDCHA(nlbd,2)*1.1
            !BEDCHA(nlbd+1,2) = max(BEDCHA(nlbd+1,2),0.1)
            BEDCHA(nlbd+1,3) = BEDCHA(nlbd,3)*1.1
          end if

          do is = 1,nscls			!update percbd in layer 1
            percbd(1,is) = (percbd(2,is)*(BEDCHA(3,1)-BEDCHA(2,1)) +
     @percbd(1,is)*BEDCHA(2,1))/BEDCHA(3,1)
          end do

          do ib = 2,nlbd			!shift one layer up
           do is = 1,nscls
            percbd(ib,is) = percbd(ib+1,is)
           end do
           BEDCHA(ib,1) = BEDCHA(ib+1,1)
           BEDCHA(ib,2) = BEDCHA(ib+1,2)
           BEDCHA(ib,3) = BEDCHA(ib+1,3)
          end do

          if (BEDCHA(3,1).lt.bmix) goto 125

        end if

!       -------------------------------------------------------------------
!       Updates layer 1 (level 2) by linear interpolation
!       -------------------------------------------------------------------

        do is = 1,nscls
          percbd(1,is) = (percbd(1,is)*BEDCHA(2,1) + percbd(2,is)*
     @(bmix-BEDCHA(2,1)))/bmix
        end do

        BEDCHA(2,3) = lin(BEDCHA(2,1),BEDCHA(3,1),BEDCHA(2,3),
     @BEDCHA(3,3),bmix)
        BEDCHA(2,2) = lin(BEDCHA(2,1),BEDCHA(3,1),BEDCHA(2,2),
     @BEDCHA(3,2),bmix)
        !BEDCHA(2,2) = max(BEDCHA(2,2),0.1)
        BEDCHA(2,1) = bmix

!       -------------------------------------------------------------------
!       Adjusts percbd and percentage check ---> ptot must be = 1 +-0.1
!       -------------------------------------------------------------------

        ptot = 0.
        do is = 1,nscls
         ptot = ptot + percbd(1,is)
        end do
        if(ptot.lt.0.99.or.ptot.gt.1.01) go to 130

        do is = 1,nscls
         percbd(1,is) = percbd(1,is)/ptot
        end do

!       -------------------------------------------------------------------
!       Resets valiables of not used rows
!       -------------------------------------------------------------------

        BEDCHA(1,1) = 0.

        do ib = nlbd+1,nlbdim
         BEDCHA(ib,1) = 0.
         BEDCHA(ib,2) = 0.
         BEDCHA(ib,3) = 0.
        end do

        return

 130    continue
        write(6,*) 'Error in computing the sediment fraction: 1'
        write(6,*) 'total percentage:',ptot
        do is = 1,nscls
          write(6,*) 'classnumber:',is,'percentage:',percbd(1,is)
        enddo
        stop 'error stop : percbd'

        end

! ********************************************************************
! SUBROUTINE CHECKBED
! The bed is divided in several layers, each with its own sediment 
! composition. At the beginning all the layers has the same volume.
! Layer = 1 is the top layer (close to water surface). 

        subroutine checkbed(k,iss,nscls,gs,BEDCHA,percbd,tau0,flux)

        implicit none

        include 'sed_param.h'
        
        double precision percbd(nlbdim,nsdim)	!fraction of sediment [0,1]
        double precision BEDCHA(nlbdim,3)  	!bed characteristics in 3 column table
                                               	! (1) depth below sediment surface (m)
                                               	! (2) critical erosion stress (Pa)
                                               	! (3) dry bulk density (kg/m**3)
        integer nlbd				!number of bed layer in k
        integer nscls				!number of grainsize class
        double precision gs			!grainsize of first class
        double precision tau0 			!effective skin friction shear stress (Pa)
        double precision flux			!flux of sediment [m] (<0 ero, >0 depo)
        double precision bnet			!new net sediment in the layer [m]
        double precision btot			!total sediment [m]
        double precision bpre                   !sediment is present in the layer [m]
        real ptot				!total percentage
        integer k,ib,is,iss

!       -------------------------------------------------------------------
!       Finds the number of rows that are used (NBED). This is indicated by
!       z=0 in the row after the last used row.
!       -------------------------------------------------------------------

        do 10, ib=2,nlbdim
          if (BEDCHA(ib,1).eq.0.0) goto 20
10      continue
20      nlbd=ib-1

!       -------------------------------------------------------------------
!       Initializes valiables
!       -------------------------------------------------------------------

        btot = BEDCHA(2,1) + flux
        bpre = BEDCHA(2,1)*percbd(1,iss)
        bnet = bpre  + flux

!       -------------------------------------------------------------------
!       Computes erosion or deposition
!       -------------------------------------------------------------------

        if (flux.gt.0.) then

          call depbed(iss,nlbd,nscls,BEDCHA,percbd,tau0,flux,bnet,
     +		      btot,gs)

        elseif(flux.lt.0.) then

          call erobed(iss,nlbd,nscls,BEDCHA,percbd,flux,bnet,btot,bpre)

        end if

!       -------------------------------------------------------------------
!       Adjusts percbd and percentage check ---> ptot must be = 1 +-0.1
!       -------------------------------------------------------------------

        if (nscls.eq.1) percbd(1,nscls) = 1.

        ptot = 0.
        do is = 1,nscls
         ptot = ptot + percbd(1,is)
        end do
        if(ptot.lt.0.99.or.ptot.gt.1.01) go to 130

        do is = 1,nscls
         percbd(1,is) = percbd(1,is)/ptot
        end do

!       -------------------------------------------------------------------
!       Resets valiables of not used rows
!       -------------------------------------------------------------------

        BEDCHA(1,1) = 0.

        do ib = nlbd+1,nlbdim
         BEDCHA(ib,1) = 0.
         BEDCHA(ib,2) = 0.
         BEDCHA(ib,3) = 0.
        end do

        return

 130    continue
        write(6,*) 'Error in computing the sediment fraction: 2'
        write(6,*) 'total percentage:',ptot,'node:',k
        do is = 1,nscls
          write(6,*) 'classnumber:',is,'percentage:',percbd(1,is)
        enddo
        stop 'error stop : percbd'

        end

! ********************************************************************

        subroutine depbed(iss,nlbd,nscls,BEDCHA,percbd,tau0,flux,
     +			  bnet,btot,gs)

        implicit none

        include 'sed_param.h'

        double precision percbd(nlbdim,nsdim)   !fraction of sediment [0,1]
        double precision BEDCHA(nlbdim,3)         !bed characteristics in 3 column table
                                                ! (1) depth below sediment surface (m)
                                                ! (2) critical erosion stress (Pa)
                                                ! (3) dry bulk density (kg/m**3)
        integer nlbd                            !number of bed layer in k
        integer nscls                           !number of grainsize class
        double precision gs                     !grainsize of first class
        double precision tau0 			!effective skin friction shear stress (Pa)
        double precision flux                   !flux of sediment [m] (<0 ero, >0 depo)
        double precision bnet                   !new net sediment in the layer [m]
        double precision btot                   !total sediment [m]
        double precision bres                   !sediment that is left [m]
        double precision blimit                 !maximum thikness of layer [m]
        double precision rhosa                  !non consolidated density [kg/m3]
	double precision taunew			!new critical shear stress
	double precision tausurf		!average surface critical shear stress
	double precision rhosurf		!average surface density
        integer ib,is,iss

!       -------------------------------------------------------------------
!       Initializes valiables
!       -------------------------------------------------------------------

        rhosa = RHOSED * (1 - POROS)
        if (gs.lt.limcoh) rhosa = RHOMUD
        blimit = BEDCHA(3,1)
        taunew = teroa*(rhosa**terob)
        tausurf = (BEDCHA(1,2) + BEDCHA(2,2))/2.
        rhosurf = (BEDCHA(1,3) + BEDCHA(2,3))/2.

!       -------------------------------------------------------------------
!       Normal deposition
!       -------------------------------------------------------------------

        do is = 1,nscls
          percbd(1,is) = percbd(1,is)*BEDCHA(2,1)/btot
        end do

        percbd(1,iss) = bnet / btot

        BEDCHA(1,2) = (taunew*flux + tausurf*BEDCHA(2,1)) / btot
        BEDCHA(1,2) = max(BEDCHA(1,2),taunew)
        BEDCHA(1,2) = min(BEDCHA(1,2),BEDCHA(2,2))

        BEDCHA(1,3) = (rhosa*flux + rhosurf*BEDCHA(2,1))/ btot
        BEDCHA(1,3) = min(BEDCHA(1,3),BEDCHA(2,3))

        do ib = 2,nlbd
         BEDCHA(ib,1) = BEDCHA(ib,1) + flux
        end do

        bres = blimit - btot
        if ( bres .lt. 0. ) then

!         -------------------------------------------------------------------
!         Excess deposition, create new layer and upgrade level 2
!         -------------------------------------------------------------------

          nlbd = nlbd + 1			!shift one layer down
          do ib = nlbd,2,-1
           do is = 1,nscls
            percbd(ib,is) = percbd(ib-1,is)
           end do
           BEDCHA(ib,1) = BEDCHA(ib-1,1)
           BEDCHA(ib,2) = BEDCHA(ib-1,2)
           BEDCHA(ib,3) = BEDCHA(ib-1,3)
          end do

          BEDCHA(2,1) = -bres
          BEDCHA(1,2) = min(BEDCHA(1,2),TAU0*1.1d0)
          BEDCHA(1,2) = max(BEDCHA(1,2),taunew)
          BEDCHA(1,3) = rhosa
          nlbd = min(nlbdim-1,nlbd)

        end if

        end

! ********************************************************************
      
        subroutine erobed(iss,nlbd,nscls,BEDCHA,percbd,flux,bnet,btot,
     @bpre)

        implicit none

        include 'sed_param.h'

	double precision d0
	parameter (d0 = 0.d0)

        double precision percbd(nlbdim,nsdim)   !fraction of sediment [0,1]
        double precision BEDCHA(nlbdim,3)         !bed characteristics in 3 column table
                                                ! (1) depth below sediment surface (m)
                                                ! (2) critical erosion stress (Pa)
                                                ! (3) dry bulk density (kg/m**3)
        integer nlbd                            !number of bed layer in k
        integer nscls                           !number of grainsize class
        double precision flux                   !flux of sediment [m] (<0 ero, >0 depo)
        double precision bnet                   !new net sediment in the layer [m]
        double precision btot                   !total sediment [m]
        double precision bpre                   !sediment is present in the layer [m]
        integer ib,is,iss

        double precision lin,x1,x2,y1,y2,x      !statement function which interpolate
                                                !linearly the value Y of (X,Y), located
                                                !on the line between (X1,Y1) and (X2,Y2)

        lin(x1,x2,y1,y2,x) = y1 + (x-x1)/(x2-x1)*(y2-y1)

        if (bnet.le.0.) go to 134

!       -------------------------------------------------------------------
!       Normal erosion
!       -------------------------------------------------------------------

        do is = 1,nscls
          percbd(1,is) = percbd(1,is)*BEDCHA(2,1) / btot
        end do

        percbd(1,iss) = bnet / btot

        BEDCHA(1,2) = lin(d0,BEDCHA(2,1),BEDCHA(1,2),BEDCHA(2,2),-flux)
        BEDCHA(1,3) = lin(d0,BEDCHA(2,1),BEDCHA(1,3),BEDCHA(2,3),-flux)

        do ib = 2,nlbd
          BEDCHA(ib,1) = BEDCHA(ib,1) + flux
        end do

        return

134     continue

!       -------------------------------------------------------------------
!       Excess erosion
!       -------------------------------------------------------------------

        if (percbd(1,iss) .lt. 0.998) then

!         -------------------------------------------------------------------
!         If other classes present, limit erosion
!         -------------------------------------------------------------------

          flux = -bpre
          btot = BEDCHA(2,1) + flux
          do is = 1,nscls
            percbd(1,is) = percbd(1,is)*BEDCHA(2,1) / btot
          end do

          percbd(1,iss) = 0.
          BEDCHA(1,2)=lin(d0,BEDCHA(2,1),BEDCHA(1,2),BEDCHA(2,2),-flux)
          BEDCHA(1,3)=lin(d0,BEDCHA(2,1),BEDCHA(1,3),BEDCHA(2,3),-flux)

          do ib = 2,nlbd
            BEDCHA(ib,1) = BEDCHA(ib,1) + flux
          end do

        else

!         -------------------------------------------------------------------
!         If no other classes present, go to lower layer
!         -------------------------------------------------------------------

          call dellayer(iss,nlbd,nscls,BEDCHA,percbd,btot,bnet)

          if (bnet .eq. 0.) return

          if ( percbd(1,iss) .ne. 0. ) then

!           -------------------------------------------------------------------
!           Sediment is present in the lower layer
!           -------------------------------------------------------------------

            bpre = BEDCHA(2,1) - btot

            do is = 1,nscls
              percbd(1,is) = percbd(1,is)*BEDCHA(2,1) / btot
            end do

            percbd(1,iss) = bnet / btot

            BEDCHA(1,2)=lin(d0,BEDCHA(2,1),BEDCHA(1,2),BEDCHA(2,2),bpre)
            BEDCHA(1,3)=lin(d0,BEDCHA(2,1),BEDCHA(1,3),BEDCHA(2,3),bpre)

            do ib = 2,nlbd
              BEDCHA(ib,1) = BEDCHA(ib,1) - bpre
            end do

          else
!           -------------------------------------------------------------------
!           No is sediment in the lower layer
!           -------------------------------------------------------------------

            flux = -bpre

          end if

        end if

        end

! ********************************************************************
! SUBROUTINE DELLAYER
! Delete layers and update bed characteristics

        subroutine dellayer(iss,nlbd,nscls,BEDCHA,percbd,btot,bnet)

        implicit none

        include 'sed_param.h'

        double precision percbd(nlbdim,nsdim)   !fraction of sediment [0,1]
        double precision BEDCHA(nlbdim,3)         !bed characteristics in 3 column table
                                                ! (1) depth below sediment surface (m)
                                                ! (2) critical erosion stress (Pa)
                                                ! (3) dry bulk density (kg/m**3)
        integer nlbd                            !number of bed layer in k
        integer nscls                           !number of grainsize class
        double precision bnet                   !new net sediment in the layer [m]
        double precision btot                   !total sediment [m]
        double precision thick
        integer ib,is,iss

100     continue

        thick = BEDCHA(2,1)

        nlbd = nlbd - 1

        if ( nlbd .lt. 3 ) then              !create the third layer
          nlbd = nlbd + 1
          do is = 1,nscls
            percbd(nlbd+1,is) = percbd(nlbd,is)
          end do
          BEDCHA(nlbd+1,1) = 2.*BEDCHA(nlbd,1)! - BEDCHA(nlbd-1,1)
          BEDCHA(nlbd+1,2) = BEDCHA(nlbd,2)*1.1
          BEDCHA(nlbd+1,3) = BEDCHA(nlbd,3)*1.1
        end if

        do ib = 1,nlbd                       !shift one layer up
         do is = 1,nscls
          percbd(ib,is) = percbd(ib+1,is)
         end do
         BEDCHA(ib,1) = BEDCHA(ib+1,1) - thick
         BEDCHA(ib,2) = BEDCHA(ib+1,2)
         BEDCHA(ib,3) = BEDCHA(ib+1,3)
        end do

        BEDCHA(1,1) = 0.

        if (bnet .eq. 0.) return

        bnet = percbd(1,iss)*BEDCHA(2,1) + btot
        btot = BEDCHA(2,1) + btot
        if (bnet .lt. 0.) go to 100

        end

! ********************************************************************
! SUBROUTINE TOTCON
! Compute the total node concentration

        subroutine totcon(nscls,scon,scc,tcon)

        implicit none

        include 'param.h'
        include 'sed_param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer ilhkv(nkndim)			!number of element and node level
        common /ilhkv/ilhkv

        integer nscls				!number of grainsize class
        real scon(nlvdim,nkndim,nsdim)  	!suspended concentration
        real scc(nlvdim,nkndim,nsdim)		!cohesive suspended sediment conc (kg/m3)
        real tcon(nlvdim,nkndim)		!total concentration

        real aux
        integer k,l,is,lmax,isstart

	isstart = 1
	if (nbcc .gt. 0 ) isstart = 2

        do k = 1,nkn
         lmax = ilhkv(k)

         do l = 1,lmax
          aux = 0.

          do is = isstart,nscls		!sum non-cohesive
           aux = aux + scon(l,k,is)
          end do
          do is = 1,nbcc		!sum cohesive
           aux = aux + scc(l,k,is)
          end do

          tcon(l,k) = aux
         end do

        end do
            
        end

! ********************************************************************
! SUBROUTINE UPEDEPTH
! Update the element depth in function of erosion/deposition

        subroutine upedepth(bdh)

        implicit none

        include 'param.h'

        double precision bdh(nkndim)          !total elevation change [>0depo,<0ero]

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer nen3v(3,1)
        common /nen3v/nen3v
        real hm3v(3,1)
        common /hm3v/hm3v
        real zenv(3,1)
        common /zenv/zenv
        real hkv(nkndim), hev(neldim)
        common /hkv/hkv, /hev/hev
        real hdknv(nlvdim,nkndim)
        common /hdknv/hdknv
        real hdenv(nlvdim,neldim)
        common /hdenv/hdenv
        real areakv(nlvdim,nkndim)
        common /areakv/areakv
        real v1v(nkndim)
        common /v1v/v1v
	integer iwegv(neldim)
	common /iwegv/iwegv

        real evdep				!element depth variation
        integer ie,ii,k,is,iw

!       -------------------------------------------------------------------
!       Updates element depth
!       -------------------------------------------------------------------

        do ie = 1,nel

          evdep = 0.

          do ii = 1,3
            k = nen3v(ii,ie)
            evdep = evdep + bdh(k)
          end do

          evdep = evdep / 3.
	  if( iwegv(ie) .gt. 0 ) evdep = 0.

          hev(ie) = hev(ie) - evdep

          do ii = 1,3
            hm3v(ii,ie) = hm3v(ii,ie) - evdep
!            zenv(ii,ie) = zenv(ii,ie) - evdep
          end do

        end do

!       ------------------------------------------------------------------
!       Sets up depth vectors
!       ------------------------------------------------------------------

        call makehkv(hkv,v1v)

!        call set_ilhv
!        call set_last_layer
!        call set_ilhkv

        call setarea(nlvdim,areakv)
        call setdepth(nlvdim,hdknv,hdenv,zenv,areakv)
        call setweg(0,iw)
!        call setznv

!       ------------------------------------------------------------------
!       Sets up velocity
!       ------------------------------------------------------------------

!        call setuvd
        call ttov

        end

! ********************************************************************
! Smooths variable averaging over the neibor nodes
! weight is the node`s area. Smooths only if bed elevation changes occour
! or if bed slope is greater than angle of repose

        subroutine smooth_node(kvalue,kdiff,smooth,gdx,gdy,angrep)

        implicit none

        include 'param.h'

        real kvalue(nkndim)             !variable
        double precision kdiff(nkndim)	!istantaneous variable difference
        real smooth			!smoothing factor for morphodynamic [0-1]
        real gdx(nkndim),gdy(nkndim)	!slope gradients
        real angrep			!angle of repose [rad]
        real ksl                        !node slope angle [radian]
        real areanode,area
        real baux(nkndim)
        real saux,aaux,smm

        integer k,kn,n,i,l

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer ilhkv(nkndim)
	common /ilhkv/ilhkv
	include 'links.h'

        do k = 1,nkn
          saux = 0.
          aaux = 0.
          smm = smooth

          ksl = atan(sqrt(gdx(k)**2 + gdy(k)**2))

          if (smm.lt.1. .and. abs(kdiff(k)).gt.1e-5) then
            if(ksl.ge.angrep) smm = 0.5
	    call set_node_links(k,n)
            do i = 1,n
              kn = lnk_nodes(i)           !kn is number of neibor node
              l = ilhkv(kn)
              area = areanode(l,kn)
              saux = saux + kvalue(kn)*area
              aaux = aaux + area
            end do
            baux(k) = kvalue(k)*smm + (saux / aaux)*(1. - smm)
  
          else
            baux(k) = kvalue(k)
          end if
        end do

        do k = 1,nkn
          kdiff(k) = kdiff(k) + (baux(k)-kvalue(k))
          kvalue(k) = baux(k)
        end do

        end

!******************************************************************
