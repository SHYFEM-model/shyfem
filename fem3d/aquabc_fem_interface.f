c
c $Id: aquabc_fem_interface.f
c
c aquabc interface routines
c
c Interface between SHYFEM and AQUABC (Aquatic Biochemical Cycling)
c Adapted by Petras Zemlys and Ali Erturk from an interface between 
C EUTRO and SHYFEM called BIO3D 
c
c contents :
c
c   subroutine aquabc_fem_interface - Interface between SHYFEM and AQUABC (Aquatic Biochemical Cycling)
c   subroutine setload
c   subroutine check_var  - checking strange values. Atention! nstate and nsstate are fixed in the code
C   subroutine check2Dbio
c   function max_int
c   subroutine inicfil_petras_ali - initial state of WC
c   subroutine calc_time     - not used anymore
c   subroutine get_ITOT      - reads light for aquabc. Calculation of light should be implemented also. Fixme
c   subroutine biotser_init  - initialisation of ASCII output for given nodes
C   subroutine biotser_write - writing of ASCII output for given nodes
C   subroutine cur_wmeteo    - interface for wind and air temperarature with SHYFEM (temporary solution, promissed to fix by Georg)
c   subroutine cur_param_read - reading WC model parameters(constants)
c   subroutine cur_param_sed_read - reading BS model parameters(constants)
c   subroutine auquabcini    - does intialisation for WC model
c   routines to read forcing time series for water temp, pH and interpolation of their values  
c
c revision log :
c
c 20.07.2004    The structure has been modified to allow dynamic
c               forcings to read from external files
c    07.2006    EUTRO changed by new eutrofication module AQUABC with 
c               14 state variables
c    07.2009    Water column eutrophication module changed to ALUKAS
c               (Advanced Level nUtrient Kinetics for Aquatic Systems)
c 14.03.2012    ggu     fixed too long lines for gfortran
c 14.03.2012    ggu     added dummies for restart
c 21.10.2014    ggu     converted to new boundary treatment
c
c notes :
c
c                  ccccccccccccccccccccccccc
c                  c    WC state variables c
c                  ccccccccccccccccccccccccc
c        NO                              STATE VARIABLE   NAME IN ALUKAS
c        --    ----------------------------------------   --------------       
c         1                           AMMONIUM NITROGEN             NH4N
c         2                            NITRATE NITROGEN             NO3N
c         3                   ORTHOPHOSPHATE PHOSPHORUS             PO4P
c         4             PHYTOPLANKTON CARBON FOR GREENS          PHYTC_G
c         5       EXT. LABILE DISSOLVED DETRITUS CARBON        EXLADDETC
c         6                            DISSOLVED OXYGEN             DOXY
c         7     EXT. LABILE PARTICULATE DETRITUS CARBON        EXLAPDETC
c         8       EXT. REFRACTORY DISS. DETRITUS CARBON        EXREDDETC
c         9                          ZOOPLANKTON CARBON             ZOOC
c        10                          NONBIOGENIC SILICA              ISI
c        11       EXT. REFRACTORY PART. DETRITUS CARBON        EXREPDETC
c        12            PHYTOPLANKTON CARBON FOR DIATOMS          PHYTC_D
c        13      PHYTOPLANKTON CARBON FOR CYANOBACTERIA          PHYTC_C
c        14                            INORGANIC CARBON           INCARB
c        15           GREENS. DISSOLVED DETRITUS CARBON        GPHYDDETC
c        16         GREENS. PARTICULATE DETRITUS CARBON        GPHYPDETC
c        17            DIATOM DISSOLVED DETRITUS CARBON        DPHYDDETC
c        18          DIATOM PARTICULATE DETRITUS CARBON        DPHYPDETC
c        19         CYANOBACTERIA DISS. DETRITUS CARBON        CPHYDDETC
c        20         CYANOBACTERIA PART. DETRITUS CARBON        CPHYPDETC
c        21           ZOOPLANKTON DISS. DETRITUS CARBON        ZOOPDDETC
c        22           ZOOPLANKTON PART. DETRITUS CARBON        ZOOPPDETC
c        23          BS BASED ORGANIC CARBON                   SEDB_DOC    
c        24          BS BASED ORGANIC Nitrogen                 SEDB_NOC
c        25          BS BASED ORGANIC CARBON                   SEDB_POC
c
c                  ccccccccccccccccccccccccc
c                  c    BS state variables c
c                  ccccccccccccccccccccccccc
c        NO    NAME IN BS model       STATE VARIABLE   
c        --    ----------------------------------------  
c        1     SED_NH4N         BS total amonia (in solutes and solids?) mg/l of sediments 
c        2     SED_NO3N         BS nitrates,                   mg/l in pore water 
c        3     SED_DON          BS dissolved organic nitrogen, mg/l in pore water  
c        4     SED_PON          BS particulate organic nitrogen
c        5     SED_PO4P         BS dissolved ortophosphates phophorus
c        6     SED_DOP          BS dissolved organic phophorus
c        7     SED_POP          BS particulate organic phosphorus
c        8     SED_DOXY         BS dissolved oxygen
c        9     SED_DOC          BS dissolved organic carbon  
c       10     SED_POC          BS particulate organic carbon  
c       11     SED_DSi          BS dissolved silica
c       12     SED_PSi          BS particulate silica
c
cccccccccccccccccccccccccccccccccccccccccccccccc
c Output:
c  *.bio - water column variables     (unit iub)
c  *.bs  - bottom sediments variables (unit iubs)
cccccccccccccccccccccccccccccccccccccccccccccccc
c
c**************************************************************

        subroutine ecological_module(it,dt)

c general interface to ecological module

        implicit none

        integer it
        real dt

        call aquabc_fem_interface(it,dt)

        end

c**************************************************************

	subroutine aquabc_fem_interface(it,dt)

      implicit none

      include 'param.h'       !Geometrical parameters(from SHYFEM)
      include 'aquabc.h'      !Main aquabc parameters
      include 'aquabc_aout.h' ! variables and parameters for output to ASCII files

      integer it	!time in seconds
      real dt	!time step in seconds!495vers real dt remove idt 

	
      integer ls

	integer narr
	parameter( narr = 100 )
	
	integer nnode
	
      real par(800)        !aquabc model WC parameters array
        save par
      real par_sed(100)    !aquabc model BC parameters array
        save par_sed

	real e(nlvdim,nkndim,nstate)	!state vector
c	real eb(nlvdim,nkndim,nstate)	!boundary vector of state vectors
	real eload(nlvdim,nkndim,nstate)!vector of loadings
	

c        real tstot(nstate)             !for mass test
c        real tsstot(nsstate)

       real es(NOSLAY,nkndim,nsstate)		!sediment state variables
       real sed_output(NOSLAY,nkndim,nsoutput)     ! BS output
       save e, es, eload				    !SAVESED

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        real eps1,eps2,flag,high,higi
        common /mkonst/ eps1,eps2,pi,flag,high,higi
        integer nlvdi,nlv
        common /level/ nlvdi,nlv

      character*80 bio2dn(1)
      common /bio2dn/bio2dn

      integer ilhkv(1)   !Number of WC layers for each node
      common /ilhkv/ilhkv
	
	integer ilhkv_sed(nkndim)	
       save ilhkv_sed     ! Needed for compatibility of WC and BS writing in biotser_write


      
        real difv(0:nlvdim,1)
        common /difv/difv
        real difhv(nlvdim,1)
        common /difhv/difhv

        character*10 what,whataux
        character*2 whatn

       integer k,i,l,lmax
       integer ibio
       integer id
       integer nintp, ivar
        
      real t,s 
      real u,v

      real eaux(nstate)  ! to keep WC variables for one node, time step and layer
      real esaux(NOSLAY,nstate) !  to keep BS variables for one node, time step and layer
      real esaux_out(NOSLAY,nsoutput)
      
      real elaux(nstate)

      real einit(nstate)       ! for initial values of WC variables

      real esinit(nsstate)     ! for initial values of BS variables 1-st layer
      real esinit2(nsstate)    ! for initial values of WC variables 2-st layer
      real esinit3(nsstate)    ! for initial values of WC variables 3-st layer
 
       real elinit(nstate)
       real ebound(nstate)
       save einit,esinit,elinit,ebound

       real tstot(nstate)   !for mass test
       real tsstot(nsstate) 
 
       integer icall,iunit
       integer j
	
      real rlux,rluxaux 
	
      real dtday
      real area ,vol,vel

       real pi

      real oxysat
      real getpar
      integer iround
      integer ieint,ipint

       integer mode
       real ai,lsurf

       integer ibsedim !1 - process BS. Duplicates ibsedim (easier to read for constants file) 
       save ibsedim
       logical bsedim  !true -process BS 	
       save bsedim

       logical bcheck
       
       logical bresi,breact,bdecay

	
      integer ie,ii
      integer kspec
	
      real depth
	
      real windspeed,tempair	
		
      real tday,t0,tsec
      real stp
        real mass
        real wsink

      integer iespecial,inspecial
      save iespecial,inspecial
      real rkpar,difmol
      save rkpar,difmol
      integer iub,itmcon,idtcon
      save iub,itmcon,idtcon
      integer iubs,itmcons,idtcons
      save iubs,itmcons,idtcons
	
	double precision dtime0
	integer itanf,nvar
	integer idbio(nbcdim)
	save idbio

	
       save icall
        
      real    PHTAB(1000), pH
      integer PHTIME(1000)
      save    PHTIME, PHTAB
      real    GET_PH
      
      real    TEMPTAB(1000),temp
      INTEGER TEMPTIME(1000)
      save    TEMPTIME, TEMPTAB
      real    GET_TEMP
      
      real a_ITOT(nkndim, 2)!total incident solar radiation, ly/day and light fraction of the day 
      real ITOT	        !total incident solar radiation, ly/day
      real FDAY

	
       integer ipv(nkndim)
       common /ipv/ipv


C     ADDED BY ALI
C     CORPI, 17 August 2004
C     07 2009 Changed by P.Zemlys. All declarations for output moved to aquabc_aout.h 

c Variables used for WC state variables output to ASCII files 
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
       save    NBIOTS     !declared in aquabc_aout
       save    BIOTSFUN  !declared in aquabc_aout
       save    BIOTSNOD  !declared in aquabc_aout

c Variables used for BS state variables output to ASCII files 
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
       save    NBIOTS_sed     !declared in aquabc_aout
       save    BIOTSFUN_sed  !declared in aquabc_aout
       save    BIOTSNOD_sed  !declared in aquabc_aout

c      variables used for intermediate WQ results(processes) output(processes) to ASCII files
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
       integer dg_count	        !do loop parameter
       real    dg  (nlvdim,nkndim,nstate,NDIAGVAR)

       save    NDGTS   !declared in aquabc_aout
       save    DGTSFUN !declared in aquabc_aout
       save    DGTSNOD !declared in aquabc_aout     
      
c      variables used for intermediate bottom sediments results(processes) output(processes)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 	
       real    dg_sed  (NOSLAY,nkndim,nsstate,NDIAGVAR) ! diagnostic data

       save    NDGTS_sed   !declared in aquabc_aout
       save    DGTSFUN_sed !declared in aquabc_aout
       save    DGTSNOD_sed !declared in aquabc_aout
 

      integer ulogbio,ifileo
      save ulogbio
	
C Variables to calculate state variable derivatives with dynamic volumes	
      real vololds(nlvdim,nkndim)
      save vololds
      real volold
      
      integer ivfirst
      data ivfirst /1/
      save ivfirst 
      
      integer ifirst_sed   !indicator of first call for BS, 1 - if first, 0 - if not
      data ifirst_sed /1/
      save ifirst_sed    
       
      real SED_DENSITIES(NOSLAY) ! temporary. fixme!
      real SED_POROSITIES(NOSLAY)
      real SOLID_CONCS(NOSLAY)
      real SOLUTE_FRACTIONS_NH4N(NOSLAY)
      real SOLUTE_FRACTIONS_PO4P(NOSLAY)
      real WATER_DENSITY
         
c------------------------------------------------------------------
c	initial and boundary conditions  [mg/l]			??
c	initial loadings                 [g/m**2/sec]		??
c------------------------------------------------------------------

C  Initialisation of BS state variables.  Sandy sed are assumed.! Fixme!
cccccccccccccccccccccccccccccccccccccccccccccccc
!data esinit  /0.1,0.1,0.1,0.1,0.02,0.1,0.1,0.1,0.1,0.1,1.,0.1/ old initial conditions
c                   1      2      3      4      5      6      7     
      data esinit /0.5000,0.0319,1.000,1000.0,0.1800,0.1400,10.0,
     *      0.0141,105.0,1400.0 ,3.500,10.000/
c            8       9     10      11     12
c                   1      2      3      4      5      6      7  
      data esinit2 /1.2000,0.0432,0.707 ,1500.,0.160,0.1200,130.0,
     *      0.000001,90.0  ,1300.,8.000,130.0/
c            8       9     10      11     12     
c                   1      2      3      4      5      6      7               
      data esinit3 /1.5000,0.0509,0.5000,1500.,0.150,0.4000,165.0,
     *      0.000001,90.0  ,1000.,8.000,165.0/       
c            8       9     10      11     12

c------------------------------------------------------------------
	data icall /0/
c------------------------------------------------------------------
        bresi = .false.
        breact = .true.
        bdecay = .false.
c       bsedim = .false.      !.true. if sediment dynamics is simulated
c       bsedim = .true.      ! Initializiation moved after aquabcinit

c	 bcheck = .false.
        bcheck = .true.       !.true. if state variables are checked    

        what = 'lagvebio'
c----------------------------------------------------------
c----------------------------------------------------------
c----------------------------------------------------------
c Initialization section (executed only first time step)
c----------------------------------------------------------
c----------------------------------------------------------
c----------------------------------------------------------

       if( icall .le. -1 ) return 
       
       
c----------------------------------------------------------
       if( icall .eq. 0 ) then    

	    ibio = iround(getpar('ibio'))
	    if( ibio .le. 0 ) icall = -1
	    if( icall .le. -1 ) return

	    icall = 1

c        ---------------------------------------------
c	  parameters for transport/diffusion resolution
c        ---------------------------------------------


          rkpar=getpar('chpar')
          difmol=getpar('difmol')

          do i=1,nkndim
           ilhkv_sed(i)= NOSLAY
          enddo

	  
c-----------------------------------------------------
c     Open file for bio variables strange values check
c-----------------------------------------------------
         if(bcheck) then
	       ulogbio = ifileo(60,'logbio.txt','f','u')
	      if( ulogbio .le. 0 ) then      
             write(6,*) 'Cannot open/create file logbio.txt'
             stop
            else
             write(6,*) 'Aquabc_fem: File logbio.txt on unit: ',
     +        ulogbio
             write(ulogbio,*)
     +          '      STRANGE VALUES FOR EUTRO VARIABLES'
             write(ulogbio,*)
     +          'Var.No Value  Node Level Time  Checking moment' 
            endif
           endif
c         --------------------------------------------------
c	  initialize state variables with einit
c         --------------------------------------------------

      do k=1,nkn		!loop on nodes
            lmax = ilhkv(k)
            do l=1,lmax
             do i=1,nstate
	          e(l,k,i) = einit(i)
             end do
           end do
      end do

c----------------------------------------------------------
c	  initialize WQ state variables from external file
c----------------------------------------------------------
          call inicfil_petras_ali('bio',e,nstate)

       


c         --------------------------------------------------
c	  set boundary conditions for all WC state variables
c         --------------------------------------------------

          !nintp=2
          !call bnds_init(what,bio2dn,nintp,nstate,nb3dim,bioarr,ebound)
	  !call bnds_print('initialization',narr,bioarr)

	  call get_first_time(itanf)
          dtime0 = itanf
          nintp = 2
	  nvar = nstate
          ebound = 0.
          call bnds_init_new(what,dtime0,nintp,nvar,nkn,nlv
     +                          ,ebound,idbio)

c         --------------------------------------------------
c	  initialize eco model output to ascii files
c         --------------------------------------------------

	        call aquabcini(par,par_sed,PHTIME,PHTAB,TEMPTIME,TEMPTAB,
     *                       NBIOTS, BIOTSNOD, BIOTSFUN,
     *                       NBIOTS_sed, BIOTSNOD_sed, BIOTSFUN_sed,	          
     *                       NDGTS, DGTSNOD, DGTSFUN,
     *                       NDGTS_sed,DGTSNOD_sed,DGTSFUN_sed)
     
c         Writing initial conditions to text file for defined nodes
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          
         idtcon = iround(getpar('idtcon')) !time step for concentration output
         itmcon = iround(getpar('itmcon')) !final time for concentration output

         print *,'WRITING WC STATE VARIABLE INITIAL',
     *              ' VALUES FOR SELECTED NODES'    
     	 call biotser_write(1, 'wc',
     *                  e, noutput,nstate, dg,NDIAGVAR,
     *	                 ilhkv,nlvdim,
     *	                 itmcon,idtcon,
     *                  NBIOTSMX,NBIOTS,BIOTSNOD,BIOTSFUN,
     *                  NDGTSMX, NDGTS, DGTSNOD, DGTSFUN)

     
         ibsedim = par(790)
         
C ADDITIONAL PREPARATION FOR BS ***************************************  
         if (ibsedim.eq.1) then 
         
           bsedim = .true.
               
c         -----------------------------
c	      Temporary initialize state BS variables Fixme!
c         -----------------------------
            
           do k=1,nkn
            do i=1,nsstate
	        
	         es(1,k,i) = esinit(i)
	         es(2,k,i) = esinit2(i)
	         es(3,k,i) = esinit3(i)	              
	        
	        end do
	       enddo 
	       
c   Recalculate initial amonia and phosphates to total per sed volume	       
        if (ifirst_sed .eq. 1) then
        
C       Densities of BS layers(later should come from se. transport or input )
c       They are assigned already in aquabc! Fixme!

C       Porosities of BS layers (later should come from se. transport) Fixme!
c       They are assigned already in aquabc! Fixme!          
           SED_POROSITIES(1) = 0.40
           SED_POROSITIES(2) = 0.30
           SED_POROSITIES(3) = 0.25
                 
           SED_DENSITIES(1) = 1.75
           SED_DENSITIES(2) = 1.75
           SED_DENSITIES(3) = 1.75
           
           WATER_DENSITY = 1.
           
          do k=1,nkn 	      
           do i=1,NOSLAY 
           
             SOLID_CONCS(I) = (SED_DENSITIES(I) - 
     *                     (WATER_DENSITY * SED_POROSITIES(I)))          

             SOLUTE_FRACTIONS_NH4N(i) = 
     *      1.0 / (1.0 + (SOLID_CONCS(I) * par_sed(41)))

             SOLUTE_FRACTIONS_PO4P(I) = 
     *       1.0 / (1.0 + (SOLID_CONCS(I) * par_sed(42)))
     
             sed_output(i,k,13) = es(i,k,1)
             es(i,k,1) = es(i,k,1)*SED_POROSITIES(i)/
     *                                        SOLUTE_FRACTIONS_NH4N(i)
             sed_output(i,k,14) = es(i,k,5) 
             es(i,k,5) = es(i,k,5)*SED_POROSITIES(i)/
     *                                        SOLUTE_FRACTIONS_PO4P(i)
             do j = 1, nsstate
               sed_output(i,k,j) = es(i,k,j) 
             end do
                       
           end do !NOSLAY
           
          end do  !nkn
           
            ifirst_sed =0           
        end if   

           print *,'WRITING BS STATE VARIABLE INITIAL',
     *              ' VALUES FOR SELECTED NODES'    

              
         call biotser_write(1, 'bs',sed_output,nsoutput,nsstate,
     *              dg_sed,NDIAGVAR_sed,
     *              ilhkv_sed,NOSLAY,
     *              itmcon,idtcon,
     *              NBIOTSMX_sed,NBIOTS_sed,BIOTSNOD_sed,BIOTSFUN_sed,
     *              NDGTSMX_sed, NDGTS_sed, DGTSNOD_sed, DGTSFUN_sed)
     
     
C        Calculate solute fractions of NH4 and PO4      
 
 
         else
               bsedim =.false.
         end if
c    end of preparation for BS


c         --------------------------------------------------
c	  initialize output to binary files
c         --------------------------------------------------

	   iub = 0
       itmcon = iround(getpar('itmcon'))
       idtcon = iround(getpar('idtcon'))

       call confop(iub,itmcon,idtcon,nlv,nstate,'bio')

	   write(6,*) 'bio3d model initialized...'
	  
       if( bsedim ) then 
	      iubs = 0
          itmcons = iround(getpar('itmcon'))
          idtcons = iround(getpar('idtcon'))

          call confop(iubs,itmcons,idtcons,NOSLAY,nsstate,'bs')

	      write(6,*) 'bio3d sediment model initialized...'
	   end if
	   

      end if  ! 
      
c-------------------------------------------------------------------
c-------------------------------------------------------------------
c-------------------------------------------------------------------
c end of initialisation
c-------------------------------------------------------------------
c-------------------------------------------------------------------
c-------------------------------------------------------------------

c-------------------------------------------------------------------
c-------------------------------------------------------------------
c-------------------------------------------------------------------
c normal call (every time step)
c-------------------------------------------------------------------
c-------------------------------------------------------------------
c-------------------------------------------------------------------

         what = 'lagvebio'

	  kspec = 4350
	  kspec = -1

         wsink = 0.


C       GET TIME LIHGT AND PH SERIES VALUES FOR EACH TIME STEP
C      -------------------------------------------------------       
        call get_ITOT(a_ITOT, it, int(dt))
        
        pH = GET_PH(it, PHTIME, PHTAB)
        temp = GET_TEMP(it,TEMPTIME, TEMPTAB)        
        
c         --------------------------------------------------
c	  set loadings in the interal areas
c         --------------------------------------------------

C       MODIFIED BY ALI
C       CORPI, 19 July 2004
C       
C       NEW LOCATION OF SETLOAD
c	write(*,*) 'Bug with the INTEL8.1 compiler -->fixme'
c	call setload(eload, it, int(dt)) 
	
c	-------------------------------------------------------------------
c	time management
c	-------------------------------------------------------------------
			      
	t0    = 0.	
       dtday = dt / 86400        
	tsec  = it
	tday  = it / 86400. + t0		!time in days, FEM 0 is day t0	

       mode = +1               !new time level for volume and depth

	if( bcheck ) call check_var('BEFORE aquabc',it,ulogbio,e,es)
	
       call cur_wmeteo(tempair,windspeed)       !FIXED BY PETRAS 18.08.2004
       
c------------------------------------------------------------------------	
c	-------------------------------------------------------------------
c	loop on elements for biological reactor
c	-------------------------------------------------------------------
c------------------------------------------------------------------------       
       
	
!$OMP PARALLEL PRIVATE
!$OMP* (k,i,l,ls,lmax,dg_count,eaux,esaux,esaux_out,
!$OMP* dgar,dgar_sed,elaux,s,t,
!$OMP* vol,volold,depth,vel,u,v,area,nnode,ITOT,FDAY)
!$OMP DO SCHEDULE(STATIC)
 

	do k=1,nkn		!loop on nodes
	
	

           lmax = ilhkv(k)  !maximum number of levels for the node
           
           rlux = 1.


c*****************************************************************
c          get time series values for each node
          		          
	     ITOT = a_ITOT(k, 1)
	     FDAY = a_ITOT(k, 2)



ccccccccccccccccccccccccc
c          Loop on levels 
ccccccccccccccccccccccccc
          
        do l=1,lmax
          
	       call dvanode(l,k,mode,depth,vol,area)   !gets depth, volume and area

c              print *,'depth:',depth

        if ((.not.(depth.lt.0).and..not.(depth.ge.0)).or.depth.eq.0)
     +       then
                 print*,'CURON: Depth is NaN or zero:', depth,
     +          'on level: ', l, 'on node: ',k
                 stop
              endif

        if ((.not.(vol.lt.0).and..not.(vol.ge.0)).or.vol.eq.0)
     +        then
                 print*,'CURON: Volume is NaN or zero:', vol,
     +          'on level: ', l, 'on node: ',k
                 stop
              endif

         if ((.not.(area.lt.0).and..not.(area.ge.0)).or.area.eq.0) 
     +       then
                 print*,'CURON: Area is NaN or zero:', area,
     +          'on level: ', l, 'on node: ',k
                 stop
              endif


	    
	    	if (ivfirst.eq.1) then	        
	        volold =  vol	        
	      else
	        volold = vololds(l,k)
	      endif
	      
            call getts(l,k,t,s)     !gets temp and salt
            
            call getuv(l,k,u,v)     !gets velocities u/v
            vel = sqrt(u*u+v*v)

            id = 1000*k+l 

C Changing WC old values  to new in eaux for the current step calculations
            do i=1,nstate
	        eaux(i) = e(l,k,i)
	        elaux(i) = eload(l,k,i)
            end do
            
C Changing BS old values  to new in esaux for the current step calculations
         if(bsedim) then
           do ls=1,NOSLAY
	        do i=1,nsstate
	         esaux(ls,i) = es(ls,k,i)
	        end do
           end do
         endif 


c            if( k .eq. kspec ) write(6,*) 'bio3d 1: ', eaux
c  call eutro0d(id,tday,dtday,vol,depth,vel,t,s,rlux,eaux,elaux)

              nnode=ipv(k)              
c             t=temp
                            
	        CALL AQUABC(nnode,l,lmax,
     *                  tday,dtday,
     *                  vol,volold,area,
     *	                 depth,vel,t,windspeed,tempair,s, 
     *	                 pH,
     *                  ITOT,FDAY, 
     *                  elaux,
     *                  par,
     *                  eaux,nstate,
     *                  dgar,NDIAGVAR,
     *                  ibsedim,
     *                  par_sed,
     *                  esaux,nsstate,NOSLAY,
     *                  esaux_out, nsoutput,     
     *                  dgar_sed, NDIAGVAR_sed)



     
             vololds(l,k)=vol 
              
            if( k .eq. kspec ) write(6,*) 'bio3d 3: ',eaux
c
C Assignement of calculated WC values for the next step
            do i=1,nstate
	         e(l,k,i) = eaux(i)	        
             do dg_count = 1, NDIAGVAR
                  dg(l,k,i,dg_count) = dgar(i, dg_count)
	         end do	      
            end do          

      end do 

cccccccccccccccccccccccccccccccccc
c            !End loop on levels
cccccccccccccccccccccccccccccccccc
c------------------------------------------------------
	    l = lmax


          if( bsedim ) then     
     
C Assignement of calculated BS values for the next step for each node          
	    
	      do i=1,nsstate
	       do ls=1,NOSLAY
	        es(ls,k,i) = esaux(ls,i)
	           do	dg_count = 1, NDIAGVAR_sed
                  dg_sed(ls,k,i,dg_count) = dgar_sed(ls,i, dg_count)
	           end do	   
	       end do
	      end do
	      
	      
	      do i=1,nsoutput
	       do ls=1,NOSLAY
	        sed_output(ls,k,i) = esaux_out(ls,i)	             
	       end do
	      end do	      
	      
          end if !bsedim
c----------------------------------------------

	end do !k
	
!$OMP END DO NOWAIT	
!$OMP END PARALLEL
 		
c---------------------------------------------------------------------	
c       --------------------------------------------------------------	
c	end of loop on nodes for biological reactor
c       --------------------------------------------------------------
c---------------------------------------------------------------------
      
      ivfirst=0 !for dynamic volumes        
      




c      write(17,*) 'AFTER eutro0d:', 'time=',it
c      write(17,1000) (ipv(k),(e(1,k,i),i=1,9),k=1,nkn)
	 
c	-------------------------------------------------------------------
c	advection and diffusion
c	-------------------------------------------------------------------

	if( bcheck ) call check_var('BEFORE advection',it,ulogbio,e,es)
       
       !call scal_bnd(what,tsec,bioarr)
	
!$OMP PARALLEL PRIVATE(i)
!$OMP DO SCHEDULE(STATIC)	

	do i=1,nstate	  
	
          call scal_adv(what,i
     +                          ,e(1,1,i),idbio
     +                          ,rkpar,wsink
     +                          ,difhv,difv,difmol)

          call tsmass (e(1,1,i),1,nlvdim,tstot(i)) !mass control     
   
       
	end do
	
!$OMP END DO NOWAIT	
!$OMP END PARALLEL	

! Commented temporary for experiment with sediments (check bio3d!). es is 3d now. Is this check necessary?
!      do i=1,nsstate
!          call scalmass(es(1,i),0.1,tsstot(i))   !mass ctrl sed
!      end do
      
!      write(18,*) it,(tstot(i),i=1,nstate)
!      write(17,*) it,(tsstot(i),i=1,nsstate)
!  21  format (9(f18.4,2x))     

c	    -------------------------------------------------------------------
c	    write of results (file binary BIO and text)
c	    -------------------------------------------------------------------
      
c         write(17,*) 'AFTER scal3sh:','time=',it
c         write(17,1000) (ipv(k),(e(1,k,i),i=1,9),k=1,nkn)
1000      format((I4,1x,9(F8.4,1x)))
       
	    do i=1,nstate
		   call confil(iub,itmcon,idtcon,70+i,nlvdim,e(1,1,i))
	    end do

C         ADDED BY ALI AND PETRAS
C         CORPI, 9 August 2004

c            print *,'WRITING WC STATE VARIABLE ',
c     *              ' VALUES FOR SELECTED NODES'        
	     call biotser_write(0, 'wc',e,noutput, nstate, dg,NDIAGVAR,
     *	                ilhkv,nlvdim,
     *	                itmcon,idtcon,
     *                  NBIOTSMX,NBIOTS,BIOTSNOD,BIOTSFUN,
     *                  NDGTSMX,NDGTS, DGTSNOD, DGTSFUN)
     
 
c          sediment ouput

          if( bsedim ) then    

	     do i=1,nsstate
            call confil(iubs,itmcon,idtcon,100+i,NOSLAY,es(1,1,i))
	     end do

c            print *,'WRITING BS STATE VARIABLE ',
c     *              ' VALUES FOR SELECTED NODES'  

          call biotser_write(0, 'bs',sed_output,nsoutput,nsstate,
     *                   dg_sed,NDIAGVAR_sed,
     *                    ilhkv_sed,NOSLAY,
     *                    itmcon,idtcon,
     *            NBIOTSMX_sed,NBIOTS_sed,BIOTSNOD_sed,BIOTSFUN_sed,
     *            NDGTSMX_sed, NDGTS_sed, DGTSNOD_sed, DGTSFUN_sed)
           end if

       
c--------------------------------------------------------------
c--------------------------------------------------------------      
c     Averages, min,max
c--------------------------------------------------------------
c--------------------------------------------------------------      
	call bio_av_shell(e)		!aver/min/max of state vars
c	call sed_av_shell(es)		!aver/min/max of sed var

	if( bcheck ) call check_var('AFTER advection',it,ulogbio,e,es)

c	-------------------------------------------------------------------
c	debug output
c	-------------------------------------------------------------------

c	-------------------------------------------------------------------
c	end of aquabc_fem_interface
c	-------------------------------------------------------------------

	end

c*************************************************************

c*************************************************************

	subroutine setload(eload, it, idt)

c sets up eload which is loading for specified areas
c
c the computed loadings in eload are in [g/(m**3 day)] == [mg/(l day)]
c the specified loadings in areaload are in [kg/day]
c
c variables to be specified:
c
c nimmis        total number of areas for which loading is specified
c nodes         total number of nodes used to identify all areas
c karee         node numbers that specify the areas of loading
c iaree         area numbers for the nodes [1-nimmis]
c areaload      total loadings [kg/day] for areas
c
c the node numbers in karee are external node numbers
c
c
C SUBROUTINE UPDATED BY ALI AND PETRAS TO ALLOW DYNAMIC LOADINGS,
C WHICH WILL BE READ FROM AN EXTERNAL FILE.
C
C CORPI, 20 July 2004 -----> Main updates on subroutine SETLOAD
C
C
C CORPI, 22 July 2004 -----> - SETLOAD corrected to overjump loading
C                              time intervals before simulation start
C                             
C                            - New header lines are added to the loading
C                              file to fill in some usefull information
C
C                            - A second (alternative) file format and
C                              structure has been developed. The new 
C                              structure is a better alternative if
C                              time series with different time intervals 
C                              are to be read for each load.
C  
C CORPI, 23 July 2004 -----> - Error TAKING ONE DAYS LOADING FOR EACH 
C                              TIME STEP has been fixed
     

	implicit none

      include 'param.h'
	integer nstate
	parameter(nstate=9)

C     MODIFIED BY ALI
C     CORPI, 15 July 2004
C     TAKE CARE
C     eload(3,neldim,nstate) -----> eload(nlvdim,nkndim,nstate)	
	real eload(nlvdim,nkndim,nstate)

      integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
      common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	include 'ev.h'
	integer nen3v(3,1)
	common /nen3v/nen3v
	real hm3v(3,1)
	common /hm3v/hm3v

C	ADDED BY PETRAS AND ALI
C     CORPI, 19 July 2004
C     New integer variable added nimmax
C     
C     nimmax : Maximum number of loadings allowed for the compiled
C              executable image. For this exectuteable 50. If more
C              loadings are needed, increase nimmax and recompile.
C
	integer nimmax
	parameter (nimmax=50)

C     Actual mumber of loadings. 1..nimmax
	integer nimmis
	save nimmis

	real volaux(nimmax)
	real areaload(nimmax,nstate)
	save volaux,areaload


C     ADDED BY ALI
C     CORPI, 22 July 2004 
C     New variable pareaload, narealoaf
C
C     nareaload(j, jj) : Next loading for load j, state variable jj

	real nareaload(nimmax,nstate)
      save nareaload

C     ADDED BY ALI
C     CORPI, 19 July 2004
C
C     New variables it, idt, itload, preitl
C     
C	it	   : time in seconds
C	idt	   : time step in seconds
C	itload : time of load interval
C	preitl : time of prevois load interval

	integer it
	integer idt
	integer itload
	integer preitl

	save itload, preitl
      

C     ADDED BY ALI
C     CORPI, 22 July 2004
C
C     New variables it, idt, itload, preitl
C     
C	itloa2(j) : time of load interval for 2nd type loading file load j
C	preit2(j) : time of prevois load interval for second type load j
	                    
	integer itloa2(nimmax) !time of load interval for 2nd type loading file
	integer preit2(nimmax) !time of prevois load interval for second type
	                       
	save itloa2, preit2

    	integer aree(nkndim)


C	ADDED BY PETRAS AND ALI
C     CORPI, 19 July 2004
C     New integer variable added nimmax
C     
C     nimmax : Maximum number of nodes with loadings allowed for 
C              the compiled executable image. For this exectuteable 
C              5000. If moreloadings are needed, increase nimmax 
C              and recompile.
C
      integer nodmax
	parameter(nodmax=5000)

C     Actual mumber of nodes with loadings. 1..nodmax

	integer nodes
	save nodes !ADDED BY PETRAS 12-10-2004

C     ADDED BY ALI
C     CORPI, 22 July 2004
C     New integer variable added lftype
C
C     lftype : Loading file type
C
C              lftype = 1 ----> Loading file keeps all the loading
C                               information
C
C              lftype = 2 ----> Loading file keeps basic loading
C                               information and names of the time series
C                               files for each load
      integer lftype

	integer karee(nodmax)
	integer iaree(nodmax)
	save karee,iaree
	

C	ADDED BY PETRAS AND ALI
C     CORPI, 19 July 2004
C     New integer variable added
C     
C     icall : Is it necessary to read information about 
C                  loading areas from main loads  file 
C             icall = 0 ---> Need to read
C             icall = 1 ---> information is already readed

      integer icall	
	save icall
      data icall /0/


C	ADDED BY ALI
C     CORPI, 19 July 2004
C     New integer variable added
C     
C     ifirst : Is it the the first time for reading loading data
C              ifirst = 0 ---> Reading loading data for the first time
C              ifirst = 1 ---> Reading loading data not for the first time

      integer ifirst	
	save ifirst
      data ifirst /0/


C	ADDED BY ALI
C     CORPI, 22 July 2004
C     New integer variable added
C
C     
C     ifirs2 : Is it necessary to read the next time intervall for loads
C
C              ifirs2(j) = 0 ---> Need to read the next time intervall 
C                                 for load j 
C
C              ifirs2(j) = 1 ---> Do not need to next time intervall
C                                 for load j

      integer ifirs2(nimmax)	
	save ifirs2


C	ADDED BY ALI
C     CORPI, 19 July 2004
C     New integer variable added
C     
C     inext : Is it necessary to read the next time intervall for loads
C             inext = 0 ---> Need to read the next time intervall
C             inext = 1 ---> Do not need to next time intervall

      integer inext	
	save inext
      data inext /0/

C	ADDED BY ALI
C     CORPI, 22 July 2004
C     New integer variable added
C
C     
C     inext2 : Is it necessary to read the next time intervall for loads
C
C              inext2(j) = 0 ---> Need to read the next time intervall 
C                                 for load j 
C
C              inext2(j) = 1 ---> Do not need to next time intervall
C                                 for load j

      integer inext2(nimmax)	
	save inext2


C	ADDED BY ALI
C     CORPI, 19 July 2004
C     New integer variable added
C     
C     ilast : Is last loading time interval read
C             ilast = 0 ---> Last loading time interval not read
C             inext = 1 ---> Last loading time interval read

      integer ilast	
	save ilast
      data ilast /0/

C	ADDED BY ALI
C     CORPI, 22 July 2004
C     New integer variable added
C
C     
C     ilast2 : Is it necessary to read the next time intervall for loads
C
C              ilast2(j) = 0 ---> Last loading time interval not read for load j
C              ilast2(j) = 1 ---> Last loading time interval read for load j

      integer ilast2(nimmax)	
	save ilast2

	logical berror
	integer k,ie,ii,ia,i
	integer itype
	real area
	real litri,kgs
	real fact,rlfact
	real hdep,vol,load

	real getpar


C	ADDED BY ALI
C     CORPI, 02 August 2004
      integer ipv(nkndim)	!external node numbers
      common  /ipv/ipv


C	ADDED BY PETRAS AND ALI
C     CORPI, 19 July 2004
C     New variables added : nb, file, irec
C     
C     nb     : File number for the loadings file
C     header : Header information
C     irec   : Record number
C     ivar   : Variable number 
C     j, jj  : General purposed counter for array indices, ... 
C
      integer nb
	save nb
      
	character*90 header 
      
      integer irec, ivar
	integer j, jj

C     ADDED BY ALI
C     CORPI, 22 July 2004
C     New variables added : ltsferr, ltsfun, ltsfnm 
C        
C     ltsfer : Used for error checkong when opening loading time series file
C     ltsfun : Loading time series file units
C     ltsfnm : Loading time series file names

      integer ltsfer
      data ltsfer /0/

      integer ltsfun(nimmax)
	save ltsfun
	
	character*256 ltsfnm(nimmax)
      save ltsfnm

c loading is kg/day
c
c
c	loading for areas [kg/day]
c 	
	integer ifileo
        integer max_int
        
        character*80 file
	 
C     FIX FILE NAME FOR THIS VERSION
C      file = 'INPUT/loads.dat'


C	ADDED BY PETRAS AND ALI
C     CORPI, 19 July 2004
C     
	if( icall .eq. 0 ) then

C         OPEN THE MAIN LOADINGS FILE

            call getfnm('bioload',file)    
	    nb = ifileo(90,file,'f','old')

	    if( nb .le. 0 ) goto 97
      

C         Initialize the arrays
          do j = 1, nodmax
              iaree(j) = 0
	        karee(j) = 0
	    end do

	    do j = 1, nimmax

	        do jj = 1, nstate
	            areaload (j,jj) = 0.0
                  nareaload(j,jj) = 0.0
	        end do
              
              ltsfun(j) = 0
	        ltsfnm(j) = ''
              ifirs2(j) = 0
			inext2(j) = 0
              ilast2(j) = 0
              itloa2(j) = 0
              preit2(j) = 0

	    end do
          
		  write(6,*) 'READING LOADING FILE...'

	    preitl = 0
	    ivar = 0

C         Read the RECORD 1
          irec = 1


C         Read two header lines
          read(nb, 5040, err=101) header
          read(nb, 5040, err=101) header
	    
		  write(6,*) 'RECORD 1 OF THE LOADING FILE READ'

C         Read the RECORD 2
          irec = 2

C         Read two header lines
          read(nb, 5040, err=101) header
          read(nb, 5040, err=101) header
		  write(6,*) 'RECORD 2 OF THE LOADING FILE READ'

C         Read the RECORD 3
          irec = 3

	    read(nb, 5010, err=90) nimmis, nodes, lftype

	    if(nimmis.gt.nimmax) goto 95
	    if(nodes.gt.nodmax) goto 96
          
	    write(6,*) ''
		  write(6,*) 'RECORD 3 OF THE LOADING FILE READ'
          write(6,*) '---------------------------------'
	    write(6,*) 'TOTAL NUMBER OF LOADINGS        : ', nimmis
	    write(6,*) 'NUMBER OF NODES RECIEVING LOADS : ', nodes
	    write(6,*) 'TYPE OF THE LOADING FILES       : ', lftype
	    write(6,*) ''

C         Read the RECORD 4
          irec = 4

C         Read two header lines
          read(nb, 5040, err=101) header
          read(nb, 5040, err=101) header
		  write(6,*) 'RECORD 4 OF THE LOADING FILE READ'

      
C         Read the RECORD 5
          irec = 5
          
	    do j=1, nodes

	        read(nb, 5010, err=93) iaree(j), karee(j)
              
	        if(iaree(j).lt.1.or.iaree(j).gt.nimmis) then
			    goto 91
	        end if
			
			if(karee(j).lt.1.or.karee(j).gt.max_int(ipv, nkn)) then
			    goto 92
	        end if

          end do  
	        
		  write(6,*) 'RECORD 5 OF THE LOADING FILE READ'

C         Read the RECORD 6
          irec = 6

C         Read two header lines
          read(nb, 5040, err=101) header
          read(nb, 5040, err=101) header
		  write(6,*) 'RECORD 6 OF THE LOADING FILE READ'

	    if(lftype.eq.2) then
          
C             Read the RECORD 7 of LOADING FILE TYPE 2
              irec = 7

		      write(6,*) ''
		      write(6,*) 'READING RECORD 7 OF THE LOADING FILE'
              write(6,*) '===================================='
              write(6,*) ''
          
		      do j = 1, nimmis
C                 Read loading time series file name 
	            read(nb, 5050, err=201) ltsfnm(j)
	            
C				  Open the loading time series file 
                  ltsfun(j) = ifileo((80 + j),ltsfnm(j), 'f','old')

C                 Write loading time series file information on terminal
                  write(6,*) 'FOR LOAD ', j
	            write(6,*) '------------------'

	            if(ltsfun(j).le.0) then
                      write(6,*) 'PROBLEMS ENCOUNTERED WHEN OPENING ',
     +                           'TIME SERIES FILE FOR LOADING ', j
	                write(6,*) '' 
			        ltsfer = 1
	            else
                      write(6,*) 'Unit of the loading time series file:'
     +                           ,ltsfun(j)
                      write(6,*) 'Name of the loading time series file:'
     +                           ,ltsfnm(j)
	                write(6,*) ''
	            end if

	        end do
	            
			  close(nb)

	        if(ltsfer.eq.1) goto 200
		
		end if

c         extern to intern	    
		call n2int(nodes,karee,berror)

	    if( berror) stop 'error stop: loading'
		
		icall = 1

	end if

	if((ifirst.eq.0).and.(lftype.eq.1)) then

    1     continue
C         Read the RECORD 7 of LOADING FILE TYPE 1
          irec = 7

          read(nb, 5030, err=98) itload
		  write(6,*) ''
		  write(6,*) 'RECORD 7 OF THE LOADING FILE READ'
		
		  preitl = itload

	    
C         ADDED BY ALI
C         CORPI, 22 July 2004
C         
C         The following if structure overjumps loading time interval
C         before simulation start
		  if(itload.lt.(it-idt)) then
              
              if(itload.lt.0) then
                  write(6,*) 'FOR THIS SIMULATION NO LOADS WILL BE READ'
                  write(6,*) 'ZERO LOADING ASSUMED'
	            ilast = 1
                  close(nb)
	        else
                  write(6,*) 'SIMULATION START : ', (it - idt), 
     +			' START OF THE LOAD INTERVAL : ', itload

                  write(6,*)'NEXT LOADING TIME INTERVAL WILL BE READ...'

C                 Read the RECORD 8 of LOADING FILE TYPE 1
                  irec = 8
		        
		          do j = 1, nimmis
	                read(nb, 5020, err=94)
     +				(areaload(j,ivar), ivar=1,nstate)
		          end do 
          
		          write(6,*) 'RECORD 8 OF THE LOADING FILE READ'	            
                  
                  goto 1

	        end if
	    
		end if
	    
		ifirst = 1

      end if

	if(((inext.eq.0).and.(ilast.eq.0)).and.lftype.eq.1) then
		
C         Read the RECORD 8 of LOADING FILE TYPE 1
          irec = 8
		        
		  do j = 1, nimmis
	        read(nb, 5020, err=94)(areaload(j,ivar), ivar=1,nstate)
		  end do 
          
		  write(6,*) 'RECORD 8 OF THE LOADING FILE READ'

		  preitl = itload
		  inext = 1

C         Read the RECORD 7 of LOADING FILE TYPE 1 - NEXT TIME INTERVAL
          irec = 7
		 
          read(nb, 5030, err=98) itload
		  write(6,*) 'RECORD 7 OF THE LOADING FILE READ'
          
	    if(itload.le.preitl) then
		
              if(itload.lt.0) then
	            write(6,*) 'NO MORE LOADING TIME INTERVALS'	            
				  ilast = 1
	            close(nb)
	        else
		          goto 100
              end if
	    
		end if
      
	end if

C     CHECK IF TIME FOR THE NEXT LOADING INTERVAL
      if((((it+idt).ge.itload).and.(ilast.eq.0)).and.lftype.eq.1) then
	    write(6,*) 'NEW LOADING INTERVAL STARTING NEXT TIME STEP'
		  inext = 0
	end if


C     READ LOADING DATA FOR LOADIND FILE TYPE 2
      if(lftype.eq.2) then
      
	    do j = 1, nimmis
          
	        if(ifirs2(j).eq.0) then
    
    2             continue
                  
C                 Read the RECORD 8 of LOADING FILE TYPE 2
                  irec = 8

	            read(ltsfun(j), 5060, err=94)
     +			  itloa2(j), (areaload(j,ivar), ivar=1,nstate)

		          write(6,*) 'RECORD 8 OF THE LOADING FILE READ'
		        
				  if(itloa2(j).lt.(it-idt)) then
              
                      if(itloa2(j).lt.0) then
                          write(6,*) 'FOR THIS SIMULATION NO LOADS WILL'
     +					         , ' BE READ FOR LOADING ', j
                          write(6,*) 'ZERO LOADS ASSUMED FOR ', j
	                    
						  do jj=1, nstate
                              areaload(j,jj) = 0.0
                          end do
						
						ilast2(j) = 1
                          close(ltsfun(j))
	                else
                          write(6,*) 'SIMULATION START : ', (it - idt), 
     +			        ' START OF THE LOAD INTERVAL FOR LOADING ', j, 
     +                    ' : ', itloa2(j)

                          write(6,*)'READING THE NEXT LOADING TIME ',
     +					          'INTERVAL'

                          goto 2

	                end if
	    
		        end if

                  ifirs2(j) = 1
				 
              end if
			
			if((inext2(j).eq.0).and.(ilast2(j).eq.0)) then

		        preit2(j) = itloa2(j)
		        
C                 Read the RECORD 8 of LOADING FILE TYPE 2
                  irec = 8

	            read(ltsfun(j), 5060, err=94)
     +			itloa2(j), (nareaload(j,ivar), ivar=1,nstate)
                  
	            inext2(j) = 1
     	    
                  if(itloa2(j).le.preit2(j)) then
		
                      if(itloa2(j).lt.0) then
	                    write(6,*) 'NO MORE LOADING TIME INTERVALS'	            
				        ilast2(j) = 1
	                    
						close(ltsfun(j))
	                else
		                goto 100
                      end if
	    
		        end if
						    
	        end if

C             CHECK IF TIME FOR THE NEXT LOADING INTERVAL
              if(((it+idt).ge.itloa2(j)).and.(ilast2(j).eq.0)) then
	            write(6,*) 'NEW LOADING INTERVAL STARTING NEXT TIME ',
     +			           'STEP FOR LOAD ', j

				  do jj=1, nstate
                       areaload(j,jj) = nareaload(j,jj)
                  end do

		        inext2(j) = 0
              end if

          end do

	end if
		
 5010 FORMAT(3I10)
 5020 FORMAT(9F10.0)
 5030 FORMAT(I10)
 5040 FORMAT(A90)
 5050 FORMAT(A256)
 5060 FORMAT(I10, 9F10.0)

	litri = 1000.	!litri in m*3
	kgs  = 10.e+6	!mg in kg

c	rlfact = getpar('flbio')
	rlfact = 1.


C     MODIFIED BY ALI
C     CORPI, 23 July 2004
C     ERROR ABOUT TAKING ONE DAYS LOADING FOR ECAH TIME STEP FIXED
C
C     fact = (rlfact*kgs/litri) 
C     fact = -------> (rlfact*kgs/litri) * (idt/86400.) 

C     [kg/m**3] -----> [mg/l] ------> kg/day ---------> (kg/s)*idt
	fact = (rlfact*kgs/litri) * (idt/86400.)
     


c     intialize

	do i=1,nimmis
	  volaux(i) = 0.
	end do

	do k=1,nkn
	  aree(k) = 0
	end do

	do i=1,nodes
	  k = karee(i)
	  itype = iaree(i)
	  aree(k) = itype
	end do

c     compute total volume for all areas given -> store in volaux

	do ie=1,nel
	  area = 4. * ev(10,ie)	  
	  do ii=1,3
	    k = nen3v(ii,ie)
	    hdep = hm3v(ii,ie)
	    ia = aree(k)
	    if( ia .gt. nimmis ) stop 'error stop ia'
	    if( ia .gt. 0 ) then
		   volaux(ia) = volaux(ia) + area * hdep
	    end if
	  end do
	end do

c     compute and set loading in eload [g/(m**3 day)] == [mg/(l day]

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    ia = aree(k)
	    vol = volaux(ia)
	    do i=1,nstate
c	 Corrected by Petras to avoid division by zero:     
		  if( ia .le. 0 ) then
		    load = 0.        
	    else
		    load = fact * areaload(ia,i) / vol
	    endif
	      
C           MODIFIED BY ALI
C           CORPI, 15 July 2004
C           TAKE CARE
C           eload(ii,ie,i) -----> eload(1,k,i)
		  eload(1,k,i) = load
	    end do
	  end do
	end do

   
      return   
    
   90	continue
	write(6,*) 'read error in record = ',irec
	write(6,*) '... reading file ', 'loads.dat'
      write(6,*) 'Please check your loading file'
	stop 'error stop setload'

   91 continue
      write(6,*) 'error in record = ', irec, ' at row = ', j, 'node = ',
     +            karee(j)
	write(6,*) 'undefined loading area'
      write(6,*) 'Please check your loading file'
	stop 'error stop setload'

   92 continue
      write(6,*) 'error in record = ', irec,   ' at row = ', j,
     +           ' loading area = ', iaree(j), ' node = ', karee(j)
	write(6,*) 'undefined nodes for loadings'
	write(6,*) 'Please check your loading file'
	stop 'error stop setload'

   93	continue
	write(6,*) 'read error in record = ',irec,' at row = ', j
	write(6,*) '... reading file ','loads.dat'
      write(6,*) 'Please check your loading file'
	stop 'error stop setload'

   94	continue
	write(6,*) 'read error in record = ',irec,' at row = ', j,
     +           ' ivar = ', ivar
	write(6,*) '... reading file ','loads.dat'
      write(6,*) 'Please check your loading file'
	stop 'error stop setload'
    
   95 continue
      write(6,*) 'Array dimension error :'
	write(6,*) 'This executeable image was compiled for ', nimmax,
     + ' loadings but you use ', nimmis , ' loadings. Please deacrease',
     + ' the number of loadings in the loads input file RECORD 3 or ',
     + ' change nimmax parameter in SUBROUTINE SETLOAD to ', nimmis, 
     + ' or greater and recompile.'
	stop 'error stop setload'
	   
   96 continue
      write(6,*) 'Array dimension error :'
	write(6,*) 'This executeable image was compiled for ', nodmax,
     + ' nodes for loadings but you use ', nodes , ' nodes.', 
     + ' Please deacrease the number of nodes in the loads input ', 
     + 'file RECORD 3 or change nodmax parameter in SUBROUTINE SETLOAD', 
     + ' to ', nodes, ' or greater and recompile.'
	stop 'error stop setload'

   97	continue
	write(6,*) 'Cannot open loadings file ','loads.dat'
	stop 'error stop setload'
 
   98 continue
	
	if(preitl.eq.0) then
	    write(6,*) 'read error in record = ',irec,' Please check ',
     +               'RECORD 7 of the first loading interval.'
	else
	    write(6,*) 'read error in record = ',irec,' Please check ',
     +               'RECORD 7 of the loading interval next to the ',
     +               'interval staring at ', preitl, ' secs.'
      end if
	stop 'error stop setload'
   
  100 continue
      write(6,*) 'Time error :'
      write(6,*) 'Next loading time interval starts before the ',
     +           'current time interval. Please check RECORD 7 ',
     +           'after the time interval starting at ', preitl,
     +           ' secs.'
	stop 'error stop setload'
	
  101	continue
	write(6,*) 'read error in record = ',irec
	write(6,*) 'Please check the header information format'
	stop 'error stop setload'

  200	continue
	write(6,*) 'Error when reading loading time series file(s)'
      write(6,*) 'Please check your loading file'   
	stop 'error stop setload'

  201	continue
	write(6,*) 'read error in record = ',irec
	write(6,*) 'Plese check the time series file name format'
	stop 'error stop setload'		
	
	end !setload

c*************************************************************
c*************************************************************

	subroutine check_var(title,it,ulog,e,es)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c checks WC and BS vars for strange values

C Uses:
C   subroutine check2Dbio - cheks for strange values for given ranges
C                           and variable numbers 
C Developped by G.Umgiesser
C Modified for use for biogeochemical variables by P.Zemlys
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  

	implicit none

      include 'param.h'
      include 'aquabc.h'
	
      integer ulog,it
	
      
	
      integer nslayer
	
      

      character*(*) title
      real e(nlvdim,nkndim,nstate)	        !state vector
      real es(NOSLAY,nkndim,nsstate)		!sediment state variables

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer nlvdi,nlv
        common /level/ nlvdi,nlv

        character*16 text
        integer i
        
        nslayer = NOSLAY

	!write(6,*) 'check_var: ',title

        text = 'BIO CHECK'
C Changed by Petras 11 January 2005 to pass number of variable 
       do i=1,nstate	  
        call check2Dbio(it,ulog,nlv,nkn,70+i,e(1,1,i),
     *                     -0.001,25.,text,title)
       end do
	
		do i=1,nsstate	  
        call check2Dbio(it,ulog,nslayer,nkn,100+i,es(1,1,i),
     *                     -0.001,100000.,text,title)
        end do

      end

c***************************************************************

      subroutine check2Dbio(it,ul,nlv,n,nvar,a,vmin,vmax,textgen,text)

c tests array for nan and strange values for EUTRO variables
c Made from check2Dr by P.Zemlys 11 January 2005

        implicit none

        include 'param.h'

        integer ipv(nkndim) !external numbers of nodes
        common /ipv/ipv
        integer nlv,n
        integer nvar         !state variable number
        real a(nlvdim,1)     !array of variable values for levels and nodes
        real vmin,vmax       !minimal and maximal allowed values
        character*(*) textgen,text ! '***BIO CHECK' and text indicating the
                                   ! time of checking

        logical debug,bval
        integer inan,iout,i,l,ul,it
        real val

        logical is_r_nan

        bval = vmin .lt. vmax
        debug = .true.
        inan = 0
        iout = 0
        
        return

        do i=1,n
          do l=1,nlv
            val = a(l,i)
            if( is_r_nan(val) ) then
              inan = inan + 1
              if( debug ) write(ul,'(I3,G11.3,I5,I2,I10,1X,A16)') nvar,
     +                        val,ipv(i),l,it,text
            else if(bval .and. (val .lt. vmin .or. val .gt. vmax)) then
              iout = iout + 1
              if( debug ) write(ul,'(I3,G11.3,I5,I2,I10,1x,A16)') nvar,
     +                        val,ipv(i),l,it,text
            end if
          end do
        end do

        if( inan .gt. 0 .or. iout .gt. 0 ) then
          write(6,'(2A13,A4,A16,A2,2(A5,I4))') 'CHECK2DBIO: ',
     +  textgen," (",text,") ",
     +  ' NaN=',inan,' OUT=',iout
        end if

        end

c*************************************************************


      function max_int(array, index)
      
      integer array(index)
	integer maximum

	integer index
	integer i
	
      maximum = array(1)

	do i=2, index
	 
          if(array(i).ge.maximum) then
	        maximum = array(i)
	    end if

	end do

      max_int = maximum

	end
	
c*********************************************************

       subroutine inicfil_petras_ali(name,var,nvar)
c initializes nodal value variable from file

        implicit none

        include 'param.h'

        character*(*) name          !name of variable
        real var(nlvdim,nkndim,1)   !variable to set
        real varval
        integer nvar
        integer ftype               !type of  1- homogeneous initial cond.
                                    !2-heterogeneous initial conditions. Added by Petras 12-12-2004

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real hlv(1)       !hlv(i)absolute depth of bottom of layer i
        common /hlv/hlv

        character*80 file
        integer nb,irec
        integer nkk,lmax
        integer l,k
        integer ivars,ivar
        real val
        real rlaux(nlvdim) !absolute depth of bottom of layer i
                           !nlvdim - total number of vertical levels

        integer ifileo

c-------------------------------------------------------
c get file name
c-------------------------------------------------------

        call getfnm(name,file)

        if( file .eq. ' ' ) return      !nothing to initialize

c-------------------------------------------------------
c open file
c-------------------------------------------------------

C     MODIFIED BY ALI AND PETRAS
C     CORPI, 16 July 2004
C     This modification is related to change the way how inital value
C     files for the nodes are read. Before this modification it was not
C     exactly clear how to format this file. By this modification
C     the file will be read formatted
C
C     nb = ifileo(55,file,'unform','old')
C      -------> nb = ifileo(55,file,'f','old')
C
        nb = ifileo(55,file,'f','old')
        if( nb .le. 0 ) goto 97

c-------------------------------------------------------
c read first record
c-------------------------------------------------------

      ivar = 0

        irec = 1

C     MODIFIED BY ALI AND PETRAS
C     CORPI, 16 July 2004
C     This modification is related to change the way how inital value
C     files for the nodes are read. Before this modification it was not
C     exactly clear how to format this file. By this modification
C     the file will be read formatted
C
C     read(nb,err=90) nkk,lmax,ivars
C      -----> read(nb, 5000, err=90) nkk,lmax,ivars
C
C Skip 4 lines. Added by Petras 12-12-2004
        read(nb,*)
        read(nb,*)
        read(nb,*)
        read(nb,*)
C Read control information
        read(nb, 5010, err=90) nkk,lmax,ivars,ftype
 5010 FORMAT(4I5)

      if( nkk .ne. nkn .or. lmax .gt. nlvdim ) goto 99
      if( ivars .ne. nvar ) goto 96
      if(ftype .ne. 1 .and. ftype .ne. 2) goto 91

C********************************************************
C FILE TYPE 2 (spatialy heterogeneous initial conditions)
C********************************************************
      if(ftype .eq. 1) goto 1000
c-------------------------------------------------------
c read second record (only if lmax > 0)
c-------------------------------------------------------             !

        if( lmax .gt. 1 ) then          !changed from 0 to 1 (5.3.2004) LMAX
          irec = 2
C       MODIFIED BY ALI AND PETRAS
C       CORPI, 16 July 2004
C       This modification is related to change the way how inital value
C       files for the nodes are read. Before this modification it was not
C       exactly clear how to format this file. By this modification
C       the file will be read formatted
C
C
C       read(nb, err=90) (rlaux(l),l=1,lmax)
C       -----> read(nb, 5020, err=90) (rlaux(l),l=1,lmax)
C
          read(nb, 5020, err=90) (rlaux(l),l=1,lmax)
 5020   FORMAT(8F10.0)

          do l=1,lmax
            if( hlv(l) .ne. rlaux(l) ) goto 98
          end do
        else
          lmax = 1
        end if

c-------------------------------------------------------
c read data records
c-------------------------------------------------------

        irec = 3

        do ivar=1,nvar

C     MODIFIED BY ALI AND PETRAS
C     CORPI, 16 July 2004
C     This modification is related to change the way how inital value
C     files for the nodes are read. Before this modification it was not
C     exactly clear how to format this file. By this modification
C     the file will be read formatted
C
C     read(nb, err=90) ((var(l,k,ivar),l=1,lmax),k=1,nkn)
C --> read(nb, 5030, err=90) ((var(l,k,ivar),l=1,lmax),k=1,nkn)
C

        read(nb, 5030, err=90) ((var(l,k,ivar),l=1,lmax),k=1,nkn)
 5030 FORMAT(F5.0)

      end do
      goto 1001

C******************************************************
C FILE TYPE 1 (spatialy homogeneous initial conditions)
C******************************************************
 1000 continue
      irec=3
      do ivar=1,nvar
       read(nb, *, err=90) varval
       print *,'Initial condition: ','var= ',ivar,'value=', varval
       do l=1,lmax
        do k=1,nkn
         var(l,k,ivar)=varval
        end do
       end do
      end do

c-------------------------------------------------------
c reading done -> close file
c-------------------------------------------------------
 1001 continue
        close(nb)

c-------------------------------------------------------
c initialize the other levels if only surface is given
c-------------------------------------------------------

        do ivar=1,nvar
          do k=1,nkn
            val = var(lmax,k,ivar)
            do l=lmax+1,nlvdim
              var(l,k,ivar) = val
            end do
          end do
        end do

c-------------------------------------------------------
c end of routine
c-------------------------------------------------------

        write(6,*) 'Succesfull initialization for ',name,' from file '
        write(6,*) file

        return
   90   continue
        write(6,*) 'read error in record = ',irec,' ivar = ',ivar
        write(6,*) '... reading file',file
        stop 'error stop inicfil'
   91   continue
        write(6,*) 'bad file type descriptor: value 1 or 2 is allowed'
        write(6,*) '... reading file',file
        stop 'error stop inicfil'
   96   continue
        write(6,*) 'ivars not compatible with nvar: ',ivars,nvar
        stop 'error stop inicfil'
   97   continue
        write(6,*) 'Cannot open file ',file
        stop 'error stop inicfil'
   98   continue
        write(6,*) 'levels are not the same from init file ',file
        write(6,*) (hlv(l),l=1,lmax)
        write(6,*) (rlaux(l),l=1,lmax)
        stop 'error stop inicfil'
   99   continue
        write(6,*) 'parameters are not the same from init file ',file
        write(6,*) 'nkn, lmax from file  : ',nkk,lmax
        write(6,*) 'nkn, lmax from model : ',nkn,nlvdim
        stop 'error stop inicfil'
        end

c*******************************************************************

      subroutine calc_time (startdate,starttime)

c prints stats after last time step

	implicit none
        character  date*8, time*10           !Added by Petras 11.09.2004
        character  startdate*8, starttime*10
         
	include 'femtime.h'

	write(6,1035) it,niter
 1035   format(' program stop at time =',i10,' seconds'/
     +         ' iterations = ',i10)
     
         write(6,*)
         write(6,*) 'START DATE: ', startdate !Added by Petras 11.09.2004
         write(6,*) 'START TIME: ', starttime !Added by Petras 11.09.2004
         write(6,*)  
	   call date_and_time(date,time)
         write(6,*) 'END DATE: ', date !Added by Petras 11.09.2004
         write(6,*) 'END TIME: ', time !Added by Petras 11.09.2004
         end



c********************************************************************


	subroutine get_ITOT(a_ITOT, it, idt)

C     Subroutine to read daily total light from external file
C     Developed by ALI August 2004
C      REVISIONS:
C                Corrected error output messages, 01.09.2004 by Petras

	implicit none

      include 'param.h'

      real a_ITOT(nkndim,2)    !Keeps values for total daily ligth for each node

      integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
      common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

      integer ipv(nkndim)	!external node numbers
      common  /ipv/ipv
      

C     nimmax : Maximum number of regions allowed for the compiled
C              executable image. For this exectuteable 50. If more
C              loadings are needed, increase nimmax and recompile.
C
	integer nimmax
	parameter (nimmax=50)

C     Actual number of regions 1..nimmax
	integer nimmis
	save nimmis

	real light(nimmax)
	real fday(nimmax)
	save light, fday

	real nlight(nimmax)
	real nfday(nimmax)
      save nlight, nfday

	integer it      !it     : time in seconds
	integer idt     !idt    : time step in seconds
	integer itligh  !itligh : time of ligth interval
	integer preitl  !preitl : time of ligth load interval

	save itligh, preitl
      
	integer itloa2(nimmax) !time of light interval for 2nd type loading file
	integer preit2(nimmax) !time of previous load light for second type
	save itloa2, preit2
    	
	integer aree(nkndim)

      integer nodmax
	parameter(nodmax=nkndim)

C     lftype : Total daily light file type
C
C              lftype = 1 ----> File keeps all the light information
C
C              lftype = 2 ----> File keeps basic information and names of files
	integer karee(nodmax)
	integer iaree(nodmax)
	save karee,iaree

        character*80 file
	

C     icall : Is it necessary to read the light series from file
C             icall = 0 ---> Need to call
C             icall = 1 ---> Do not need to call
      integer icall	
	save icall
      data icall /0/


C     ifirst : Is it the the first time for regional light data
C              ifirst = 0 ---> Reading light data for the first time
C              ifirst = 1 ---> Reading light data not for the first time
      integer ifirst	
	save ifirst
      data ifirst /0/


C     ifirs2 : Is it necessary to read the next time intervall for light
C     ifirs2(j) = 0 ---> Need to read the next time intervall for region j 
C     ifirs2(j) = 1 ---> Do not need to next time intervall for region j
      integer ifirs2(nimmax)	
	save ifirs2


C     inext : Is it necessary to read the next time intervall for regions
C     inext = 0 ---> Need to read the next time intervall
C     inext = 1 ---> Do not need to next time intervall

      integer inext	
	save inext
      data inext /0/

C     inext2 : Is it necessary to read the next time intervall for light
C     inext2(j) = 0 ---> Need to read the next time intervall for region j 
C     inext2(j) = 1 ---> Do not need to next time intervall for region j
      integer inext2(nimmax)	
	save inext2


C     ilast : Is last loading time interval read
C     ilast = 0 ---> Last light time interval not read
C     inext = 1 ---> Last light time interval read
      integer ilast	
	save ilast
      data ilast /0/

C     ilast2 : Is it necessary to read the next time intervall for light
C     ilast2(j) = 0 ---> Last light time interval not read for region j
C     ilast2(j) = 1 ---> Last light time interval read for region j
      integer ilast2(nimmax)	
	save ilast2

	logical berror
	integer k,ie,ii,ia,i
	integer itype
	real area

	real getpar

C     nb     : File number for the loadings file
C     header : Header information
C     irec   : Record number
C     ivar   : Variable number(not necessary) 
C     j, jj  : General purposed counter for array indices, ... 
C      integer nb
	save nb
      
	character*90 header 
      
      integer irec, ivar
	integer j, jj


C     ltsfer : Used for error checkong when opening light time series file
C     ltsfun : Light time series file units
C     ltsfnm : Light time series file names
C     filej  : Light time series current file name
      integer ltsfer
      data ltsfer /0/

      integer ltsfun(nimmax)
	save ltsfun
	
	character*256 ltsfnm(nimmax)
      save ltsfnm
        character*256 filej
	integer ifileo
	integer nb
      
	integer lftype
	save lftype


      integer max_int

C     FIX FILE NAME FOR THIS VERSION
      
	if( icall .eq. 0 ) then
        call getfnm('biolight',file)

C     OPEN THE LIGHT CONTROL FILE
	nb = ifileo(90,file,'f','old')

	if( nb .le. 0 ) goto 97
      

C         Initialize the arrays
          do j = 1, nodmax
              iaree(j) = 0
	        karee(j) = 0
	    end do

	    do j = 1, nimmax

              light(j)  = 0.0
              nlight(j) = 0.0
              fday(j)   = 0.0
              nfday(j)  = 0.0

              ltsfun(j) = 0
	        ltsfnm(j) = ''
              ifirs2(j) = 0
			inext2(j) = 0
              ilast2(j) = 0
              itloa2(j) = 0
              preit2(j) = 0

	    end do
          
		write(6,*) 'READING TOTAL DAILY LIGHT CONTROL FILE...'

	    preitl = 0
	    ivar = 0

C         Read the RECORD 1
          irec = 1


C         Read two header lines
          read(nb, 5040, err=101) header
          read(nb, 5040, err=101) header
	    
		write(6,*) 'RECORD 1 OF THE ITOT FILE READ'

C         Read the RECORD 2
          irec = 2

C         Read two header lines
          read(nb, 5040, err=101) header
          read(nb, 5040, err=101) header
		write(6,*) 'RECORD 2 OF THE ITOT CONTROL FILE READ'

C         Read the RECORD 3
          irec = 3

	    read(nb, 5010, err=90) nimmis, lftype

	    if(nimmis.gt.nimmax) goto 95
	    if((lftype.ne.1).and.(lftype.ne.2)) goto 96
          
	    write(6,*) ''
		write(6,*) 'RECORD 3 OF THE ITOT CONTROL FILE READ'
          write(6,*) '------------------------------'
	    write(6,*) 'TOTAL NUMBER OF REGIONS : ', nimmis
	    write(6,*) 'TYPE OF THE ITOT FILE   : ', lftype
	    write(6,*) ''

          if(nimmis.eq.1) then
              write(6,*) 'FOUND ONLY ONE REGION, ALL NODES IN REGION 1'
	        write(6,*) 'RECORD 4 OF THE ITOT CONTROL FILE IS SKIPPED'
	        write(6,*) 'RECORD 5 OF THE ITOT CONTROL FILE IS SKIPPED'

	        do j=1, nodmax
                  iaree(j) = 1
                  karee(j) = ipv(j)
	        end do

          else

C             Read the RECORD 4
              irec = 4

C             Read two header lines
              read(nb, 5040, err=101) header
              read(nb, 5040, err=101) header
		    write(6,*) 'RECORD 4 OF THE ITOT FILE READ'
      
C             Read the RECORD 5
              irec = 5
          
	        do j=1, nkn

	            read(nb, 5010, err=93) iaree(j), karee(j)
	            
	            if(iaree(j).lt.1.or.iaree(j).gt.nimmis) then
				    goto 91
	            end if
	            
				if(karee(j).lt.1.or.karee(j).gt.max_int(ipv, nkn))then
				    goto 92
	            end if

              end do
	        
		    write(6,*) 'RECORD 5 OF THE ITOT FILE READ'

	    end if

C         Read the RECORD 6
          irec = 6

C         Read two header lines
          read(nb, 5040, err=101) header
          read(nb, 5040, err=101) header
		write(6,*) 'RECORD 6 OF THE ITOT FILE READ'

	    if(lftype.eq.2) then
          
C             Read the RECORD 7 of LOADING CONTROL FILE TYPE 2
              irec = 7

		    write(6,*) ''
		    write(6,*) 'READING RECORD 7 OF THE ITOT FILE'
              write(6,*) '================================='
              write(6,*) ''
          
		    do j = 1, nimmis
C                 Read ITOT time series file name 
	            read(nb, 5050, err=201) ltsfnm(j)
		    filej=ltsfnm(j)
	            
C				Open the ITOT time series file 
                  ltsfun(j) = ifileo((70 + j),ltsfnm(j), 'f','old')

C                 Write ITOT time series file information on terminal
                  write(6,*) 'FOR REGION ', j
	            write(6,*) '------------------'

	            if(ltsfun(j).le.0) then
                      write(6,*) 'ERROR ENCOUNTERED WHEN OPENING ',
     +                 'LIGHT TIME SERIES FILE', ltsfnm(j)
	                write(6,*) '' 
			        ltsfer = 1
	            else
                      write(6,*) 'Unit of the ITOT time series file:'
     +                           ,ltsfun(j)
                      write(6,*) 'Name of the ITOT time series file:'
     +                           ,ltsfnm(j)
	                write(6,*) ''
	            end if

	        end do
	            
			close(nb)

	        if(ltsfer.eq.1) goto 200
		
		end if

c         extern to intern	    
		call n2int(nodmax,karee,berror)

	    if( berror) stop 'error stop: ITOT'
		
		icall = 1

	end if

	if((ifirst.eq.0).and.(lftype.eq.1)) then

    1 continue
C         Read the RECORD 7 of LOADING FILE TYPE 1
          irec = 7

          read(nb, 5030, err=98) itligh
		write(6,*) ''
		write(6,*) 'RECORD 7 OF THE ITOT FILE READ'
		
		preitl = itligh

	    
C         The following if structure overjumps ITOT time interval
C         before simulation start
		if(itligh.lt.(it-idt)) then
              
              if(itligh.lt.0) then
                  write(6,*) 'FOR THIS SIMULATION NO LIGHT WILL BE READ'
                  write(6,*) 'ZERO LIGHT ASSUMED'
	            ilast = 1
                  close(nb)
	        else
                  write(6,*) 'SIMULATION START : ', (it - idt), 
     +			' START OF THE ITOT INTERVAL : ', itligh

                  write(6,*)'NEXT LIGHT TIME INTERVAL WILL BE READ...'

C                 Read the RECORD 8 of ITOT FILE TYPE 1
                  irec = 8
		        
		        do j = 1, nimmis
	                read(nb, 5020, err=94)light(j), fday(j)
		        end do 
          
		        write(6,*) 'RECORD 8 OF THE ITOT FILE READ'	            
                  
                  goto 1

	        end if
	    
		end if
	    
		ifirst = 1

      end if

	if(((inext.eq.0).and.(ilast.eq.0)).and.lftype.eq.1) then
		
C         Read the RECORD 8 of LOADING FILE TYPE 1
          irec = 8
		        
		do j = 1, nimmis
	        read(nb, 5020, err=94)light(j), fday(j)
		end do 
          
		write(6,*) 'RECORD 8 OF THE ITOT FILE READ'

		preitl = itligh
		inext = 1

C         Read the RECORD 7 of LOADING FILE TYPE 1 - NEXT TIME INTERVAL
          irec = 7
		 
          read(nb, 5030, err=98) itligh
		write(6,*) 'RECORD 7 OF THE ITOT FILE READ'
          
	    if(itligh.le.preitl) then
		
              if(itligh.lt.0) then
	            write(6,*) 'NO MORE LIGHT TIME INTERVALS'	            
				ilast = 1
	            close(nb)
	        else
		        goto 100
              end if
	    
		end if
      
	end if

C     CHECK IF TIME FOR THE NEXT ITOT INTERVAL
      if((((it+idt).ge.itligh).and.(ilast.eq.0)).and.lftype.eq.1) then
	    write(6,*) 'NEW LIGHT INTERVAL STARTING NEXT TIME STEP'
		inext = 0
	end if


C     READ LIGHT DATA FOR LIGHT FILE TYPE 2
      if(lftype.eq.2) then
      
	    do j = 1, nimmis
          
	        if(ifirs2(j).eq.0) then
    
    2             continue
                  
C                 Read the RECORD 8 of LIGHT FILE TYPE 2
                  irec = 8

	            read(ltsfun(j), 5060, err=94)
     +			itloa2(j), light(j), fday(j)
	            
		        write(6,*) 'RECORD 8 OF THE ITOT FILE READ'
		        
				if(itloa2(j).lt.(it-idt)) then
              
                      if(itloa2(j).lt.0) then
                          write(6,*) 'FOR THIS SIMULATION NO LIGHT DATA'
     +					         , ' WILL BE READ FOR LIGHT ', j
                          write(6,*) 'ZERO LIGHT ASSUMED FOR ', j
	                    
						
                          light(j) = 0.0
                          fday(j)  = 0.0
						
						ilast2(j) = 1
                          close(ltsfun(j))
	                else
                          write(6,*) 'SIMULATION START : ', (it - idt), 
     +			        ' START OF THE LIGHT INTERVAL FOR LOADING ', j, 
     +                    ' : ', itloa2(j)

                          write(6,*)'READING THE NEXT LIGHT TIME ',
     +					          'INTERVAL'

                          goto 2

	                end if
	    
		        end if

                  ifirs2(j) = 1
				 
              end if
			
			if((inext2(j).eq.0).and.(ilast2(j).eq.0)) then

		        preit2(j) = itloa2(j)
		        
C                 Read the RECORD 8 of ITOT FILE TYPE 2
                  irec = 8

	            read(ltsfun(j), 5060, err=94)
     +				itloa2(j),nlight(j),nfday(j)
                  
	            inext2(j) = 1
     	    
                  if(itloa2(j).le.preit2(j)) then
		
                      if(itloa2(j).lt.0) then
	                    write(6,*) 'NO MORE LIGHT TIME INTERVALS'	            
				        ilast2(j) = 1
	                    
						close(ltsfun(j))
	                else
		                goto 100
                      end if
	    
		        end if
						    
	        end if

C             CHECK IF TIME FOR THE NEXT ITOT INTERVAL

              if(((it+idt).ge.itloa2(j)).and.(ilast2(j).eq.0)) then
	            write(6,*) 'NEW LIGHT INTERVAL STARTING NEXT TIME ',
     +			           'STEP FOR LIGHT ', j

				write(6,*) '(it+idt)  : ', it+idt
	            write(6,*) 'itloa2(j) : ', itloa2(j)

                  light(j) = nlight(j)
                  fday (j) = nfday (j)                  

		        inext2(j) = 0
              end if

          end do

	end if
		
 5010 FORMAT(2I10)
 5020 FORMAT(2F10.0)
 5030 FORMAT(I10)
 5040 FORMAT(A90)
 5050 FORMAT(A80)
 5060 FORMAT(I10, 2F10.0)
 
C     intialize

	do k=1,nodmax
	  aree(k) = 0
	end do

	do i=1,nodmax
	  k = karee(i)
	  itype = iaree(i)
	  aree(k) = itype
	end do


	do ie=1,nodmax
	  
	    ia = aree(ie)

	    a_ITOT(ie, 1) = light(ia)
	    a_ITOT(ie, 2) = fday(ia)

	end do

   
      return   
    
   90	continue
	write(6,*) 'read error in record = ',irec
	write(6,*) '... reading file ', 'itot.dat'
      write(6,*) 'Please check your loading file'
	stop 'error stop get_ITOT'

   91 continue
      write(6,*) 'error in record = ', irec, ' at row = ', j, 'node = ',
     +            karee(j)
	write(6,*) 'undefined light area'
      write(6,*) 'Please check your light file'
	stop 'error stop get_ITOT'

   92 continue
      write(6,*) 'error in record = ', irec,   ' at row = ', j,
     +           ' region = ', iaree(j), ' node = ', karee(j)
	write(6,*) 'undefined nodes for ITOT'
	write(6,*) 'Please check your ITOT file'
	stop 'error stop get_ITOT'

   93	continue
	write(6,*) 'read error in record = ',irec,' at row = ', j
	write(6,*) '... reading file ','itot.dat'
      write(6,*) 'Please check your ITOT file'
	stop 'error stop get_ITOT'

   94	continue
	write(6,*) 'read error in record = ',irec,' at row = ', j    
	write(6,*) '... reading light series file ',filej
        write(6,*) 'Please check your file'
	stop 'error stop get_ITOT'
    
   95 continue
      write(6,*) 'Array dimension error :'
	write(6,*) 'This executeable image was compiled for ', nimmax,
     + ' regions but you use ', nimmis , ' regions. Please decrease',
     + ' the number of regions in the ITOT control file RECORD 3 or ',
     + ' change nimmax parameter in SUBROUTINE get_ITOT to ', nimmis, 
     + ' or greater and recompile.'
	stop 'error stop get_ITOT'
	   
   96 continue
      write(6,*) 'Array dimension error :'
	write(6,*) 'This executeable image was compiled for ', nodmax,
     + ' nodes for light but you use ', nodmax , ' nodes.', 
     + ' Please deacrease the number of nodes in the light input ', 
     + 'file RECORD 3 or change nodmax parameter in', 
     + ' SUBROUTINE get_ITOT to ', nodmax, ' or greater and recompile.'
	stop 'error stop get_ITOT'

   97	continue
	write(6,*) 'Cannot open ITOT file ','itot.dat'
	stop 'error stop get_ITOT'
 
   98 continue
	
	if(preitl.eq.0) then
	    write(6,*) 'read error in record = ',irec,' Please check ',
     +               'RECORD 7 of the first loading interval.'
	else
	    write(6,*) 'read error in record = ',irec,' Please check ',
     +               'RECORD 7 of the loading interval next to the ',
     +               'interval staring at ', preitl, ' secs.'
      end if
	stop 'error stop get_ITOT'
   
  100 continue
      write(6,*) 'Time error :'
      write(6,*) 'Next loading time interval starts before the ',
     +           'current time interval. Please check RECORD 7 ',
     +           'after the time interval starting at ', preitl,
     +           ' secs.'
	stop 'error stop get_ITOT'
	
  101	continue
	write(6,*) 'read error in record = ',irec
	write(6,*) 'Please check the header information format'
	stop 'error stop get_ITOT'

  200	continue
	write(6,*) 'Error when reading light time series file(s)'
      write(6,*) 'Please check your loading file'   
	stop 'error stop get_ITOT'

  201	continue
	write(6,*) 'read error in record = ',irec
	write(6,*) 'Plese check the time series file name format'
	stop 'error stop get_ITOT'		
	
	end !get_ITOT

c*********************************************************************	
c*********************************************************************

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC	
C     SUBROUTINE DEVELOPPED BY ALI AND PETRAS
C     CORPI, 9 August 2004
C
C     THIS SUBROUTINE INITIALIZES EUTRO TIME SERIES ASCII OUTPUTS for 	
C     given nodes (stations)

C     REVISIONS:
C  
C     CORPI, 17 August 2004
C     CORPI, 31 August 2004, by Petras: Bug fix for case when output
C                                       is not required    
C                               
C	---------------------
C     Subroutine extended to read diagnostic otuput control information
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	subroutine biotser_init(nb,nstate,
     *                       NBIOTSMX,NBIOTS,
     *                       BIOTSNOD,BIOTSFUN,
     *                       NDGTSMX,NDIAGVAR,
     *                       NDGTS,DGTSNOD,DGTSFUN)

        
      integer nstate  	
	integer NDIAGVAR          !Maximum number of different intermediate variables in output	
	integer NDGTSMX           !Total number of nodes for intermediates of each state variable 	
	integer NDGTS(nstate)     !Actual number of nodes for intermediates of each state variable      
	integer DGTSFUN(NDGTSMX, nstate) !Output file units for EUTRO diagnostics time series(nodes) plots	
	integer DGTSNOD(NDGTSMX, nstate) !Nodes for intermediate variables in output
       
      integer NBIOTSMX            !Total number of time series for EUTRO
      integer NBIOTS              !Actual number of time series for EUTRO		
	integer BIOTSFUN(NBIOTSMX)  !Array keeping output file units for EUTRO time series plots	
	integer BIOTSNOD(NBIOTSMX)  ! Array keeping nodes of EUTRO time series plots
     

	character*90 header
      
	character*256 BTSERNAME(NBIOTSMX)
	character*256 DTSERNAME(NDGTSMX,nstate)

	integer i
	integer itroub
      integer DODIAG
      integer DUMMY(NDGTSMX)
      logical berror

C     INITIALIZE

      itroub = 0
      
	do i= 1, NBIOTSMX
	    BIOTSNOD(i) = 0
	    BIOTSFUN(i) = 0
      end do

C     OPEN THE MAIN BIO TIME SERIES FILE    
	
      if( nb .le. 0 ) goto 97

C      Read RECORD 1 - Read two lines
       read(nb, 5010, err=98) header
       read(nb, 5010, err=98) header
C       Read RECORD 2
       read(nb, 5020, err=99) NBIOTS, DODIAG

	
	if(NBIOTS.EQ.0) then
	    write(6,*) 'NO NODES FOR TIME SERIES OUTPUT TO ASCII FILES'
	    goto 1
	end if

	
	if(NBIOTS.LT.0) then
          write(6,*) 'NUMBER NODES FOR TIME SERIES OUTPUTS', 
     +              'CANNOT BE NEGATIVE'
	    stop 'error stop biotser_init'
	end if

	
	if (NBIOTS.GT.NBIOTSMX) then
          write(6,*) 'This executable image was compiled for ',
     +              NBIOTSMX, ' number of AQUABC time series outputs.'
          write(6,*) 'please change  NBIOTSMX in aquabc.h from'
     +              , NBIOTSMX, ' to ', NBIOTS, ' and recompile.'
	    stop 'error stop biotser_init'      
	end if

       write(6,*) 'AQUABC TIME SERIES OUTPUT TO ASCII ', 
     +            'FILES WILL BE CREATED'
     
      
	write(6,*) 'NUMBER OF  NODES WITH TIME SERIES OUTPUTS : ', NBIOTS	
      
C     Read RECORD 3
      
	do i=1, NBIOTS	    
		read(nb, 5030, err=100, end=200) BIOTSNOD(i), BTSERNAME(i)          
		BIOTSFUN(i) = ifileo(70+i,BTSERNAME(i),'f','u')          
	    if(BIOTSFUN(i).le.0) then
	        write(6,*) 'Cannot open/create the AQUABC time series ',
     +                   ' output file ', i
	        write(6,*) 'UNIT : ', BIOTSFUN(i)
	        itroub = itroub + 1
	    else
	      write(6,*) 'AQUABC time series output file ', i, ' opened.' 
              write(6,*) 'File unit : ', BIOTSFUN(i)
              write(6,*) 'File name : ', BTSERNAME(i)
	    end if
	end do
      
	if(itroub.gt.0) then
	    goto 101
	end if 
      
      call n2int(NBIOTS, BIOTSNOD, berror)

	if(berror) then
	    write(6,*) 'Problem in converting external nodes to internal'
	    stop 'error stop: biotser_init'
	end if
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC	
    1 continue
      if( DODIAG .EQ. 0) then
	    write(6,*) 'NO DIAGNOSTIC TIME SERIES OUTPUT TO ASCII FILES'
	    return
	else
	    write(6,*) 'DIAGNOSTIC TIME SERIES OUTPUT TO ASCII FILES ',
     +                'WILL BE CREATED'
	end if

C     Read RECORD 4 - Read two lines
      read(nb, 5010, err=102) header
      read(nb, 5010, err=102) header

      
      do j=1, nstate
C         Read RECORD 5     
          read(nb, 5010, err=103) header
          write(*,*) 'RECORD 5 Read for state variable ', j

c          write(6,*)'bioser_init:', 'j=',j, 'NDGTS: ',NDGTS

C         Read RECORD 6
          read(nb, 5040, err=104) NDGTS(j) 
          write(*,*) 'RECORD 6 Read for state variable ', j
	    if(NDGTS(j).LT.0) then
           write(6,*) 'NUMBER OF DIAGNOSTIC TIME SERIES OUTPUTS FOR ',
     +             'STATE VARIABLE ', j, ' IS ENTERED LESS THAN ZERO'
           write(6,*) 'NUMBER OF DIAGNOSTIC TIME SERIES OUTPUTS ',
     +		           'CANNOT BE NEGATIVE'
	        stop 'error stop biotser_init'
	    end if

	
	    if(NDGTS(j).GT.NBIOTSMX) then
              write(6,*) 'This executable image was compiled for ',
     +             NDGTSMX, ' number of EUTRO diagnostic time series ',
     +                      'outputs.'
              write(6,*) 'please change parameter NDGTSMX ',
     +                   'in AQUADIAG.H ',
     +            'from', NDGTSMX, ' to ', NDGTS(j), ' and recompile.'
	        stop 'error stop biotser_init'      
	    end if
            
	    write(6,*) 'NUMBER OF NODES FOR DIAGNOSTIC TIME SERIES ',
     +               'OUTPUTS FOR STATE VARIABLE ', j, ' : ', NDGTS(j)
	    if(NDGTS(j).GT.0) then
C             Read RECORD 7
	        do i=1, NDGTS(j)	    
		        read(nb, 5030, err=105) DGTSNOD(i,j), DTSERNAME(i,j)
          		DGTSFUN(i,j)=ifileo(70+(j*i),DTSERNAME(i,j),'f','u')          
	            if(DGTSFUN(i,j).le.0) then	            
				    write(6,*) 'Cannot open/create EUTRO diagnostic',
     +                ' time series file ', i, ' for state variable ', j	            
				    write(6,*) 'UNIT : ', DGTSFUN(i,j)	            
				    itroub = itroub + 1	        
			    else	            
				    write(6,*) 'EUTRO diagnostic time series file ', i 
     +			              ,' for state variable ', j, ' opened.'
                      write(6,*) 'File unit : ', DGTSFUN(i,j)
                      write(6,*) 'File name : ', DTSERNAME(i,j)
	        
			    end if
              
                  DUMMY(i) = DGTSNOD(i,j)

	        end do

	        if(itroub.gt.0) then
	            goto 106
	        end if


C             CONVERT EXTERNAL NODES TO INTERNAL 
              call n2int(NDGTS(j), DUMMY, berror)              
	        if(berror) then
	            write(6,*) 'Problem in converting external nodes to ',
     +		           'internal, when processing state variable ', j	        
			    stop 'error stop: biotser_init'	    
		    end if
	        do i=1, NDGTS(j)
                  DGTSNOD(i,j) = DUMMY(i)
	        end do          
	    end if          
      end do
	close(nb)

 5010 FORMAT(A90)
 5020 FORMAT(2I10)
 5030 FORMAT(I10,A30)
 5040 FORMAT(I10)     
      
      return     

   97 continue
	write(6,*) 'Cannot open AQUABC time series output'
     +	       ,' information file biotser.dat'
	stop 'error stop biotser_init'   
   98 continue
      write(6,*) 'Read error in RECORD 1 of AQUABC time series ',
     +           'ASCII output  information file '
	stop 'error stop biotser_init'
   99 continue
      write(6,*) 'Read error in RECORD 2 of AQUABC time series ',
     +           'ASCII output  information file '
	stop 'error stop biotser_init'
  100 continue
      write(6,*) 'Read error in RECORD 3 of AQUABC time series   ',
     +           'ASCII output  information file : ', i
	stop 'error stop biotser_init'  
  101 continue
      write(6,*) 'One or several AQUABC time series ASCII output',
     +           '  files could not be opened/created.'
	stop 'error stop biotser_init'
  102 continue
      write(6,*) 'Read error in RECORD 4 of AQUABC time series',
     +           ' ASCII output  information file '
	stop 'error stop biotser_init'
  103 continue
      write(6,*) 'Read error in RECORD 5 of AQUABC time series',
     +           ' ASCII output  information file '
	stop 'error stop biotser_init'

  104 continue
      write(6,*) 'Read error in RECORD 6 of AQUABC time series ',
     +           'ASCII output  information file'
	stop 'error stop biotser_init'
  105 continue
      write(6,*) 'Read error in RECORD 7 of AQUABC time series ',
     +           'ASCII output  information file :  on line', i, 
     +           'state variable ', j
	stop 'error stop biotser_init'

  106 continue
      write(6,*) 'One or several EUTRO diagnostic time ',
     +           'series output files could not be opened/created.'
	stop 'error stop biotser_init'
  200 continue
      write(6,*) 'ERROR, LESS TIME SERIES PLOT DATA AVAILABLE ',
     +            'THEN REQUIRED'
      stop 'error stop biotser_init'
	return	
	end !biotser_init

c*******************************************************************************
c*******************************************************************************
c*******************************************************************************
c*******************************************************************************
c*******************************************************************************

	subroutine biotser_write(initial,what,
     *                         e,noutput,nstate,dg,NDIAGVAR,
     *                         ilhkv,nlvd,
     *                         itmcon,idtcon,
     *                         NBIOTSMX,NBIOTS,BIOTSNOD,BIOTSFUN,
     *	                        NDGTSMX, NDGTS, DGTSNOD, DGTSFUN)
     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
c Writes state variable(water column or sediments) values for defined nodes to text files
c
c Control information is read by biotser_ini 
c
c 16 July 2006 Routine is started to be reworked by P.Zemlys  to process 
c              water column and sediment kinetic variables for 3D
c    July 2009 Finished updates to process WC and BS variables for 3D
c          
c Notes: Water column variable names are used  for bottom sediments
c        All these variables are local
c
c        To have a nice output correct format statements 5010,5011,5012,5013 for state variables
c                                                        5...
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      include 'param.h'
	include 'femtime.h'

      real e(nlvd,nkndim,noutput)	!state variables
      real dg  (nlvd,nkndim,nstate,NDIAGVAR) !intermediate variables
	
	
	 integer nstate              !number of state variables
	 integer noutput            !number of variables for the output
	   
	 integer NBIOTSMX            !Total number of time series for EUTRO
        integer NBIOTS              !Actual number of time series for EUTRO   	     
	 integer BIOTSFUN(NBIOTSMX)  !Array keeping output file units for EUTRO time series plots	
	 integer BIOTSNOD(NBIOTSMX)  !Array keeping nodes of EUTRO time series plots
	 
        integer NDIAGVAR	  !Maximum number of intermediate variables displayed in diagnostics file
        integer NDGTSMX             !Total number of time series(nodes) for EUTRO	    
        integer NDGTS(nstate)       !Actual number of time series(nodes) for EUTROC      
	 integer DGTSFUN(NDGTSMX,nstate) !Array keeping output file units for EUTRO diagnostics time series(nodes) plotsC      
	 integer DGTSNOD(NDGTSMX,nstate) !Array keeping nodes of EUTRO diagnostics time series plots
	
	 integer ilhkv(1)            !number of layers for each node
		
      integer i
      integer istate

      integer iend	
      save iend
      data iend /0/

      character*2 what !'wc' - for water column, 'bs' - for bottom sediments

      integer initial
      
      
      if(initial.eq.1) goto 1000 ! writing initial conditions

      if( it .lt. itmcon ) then
       return
      end if

      if( mod(it-itmcon,idtcon) .ne. 0 ) then
       return
      end if
	
1000  continue  ! going here if initial conditions writing only 	
	
	if((NBIOTS.EQ.0).or.(iend.eq.1)) then
	    return
	end if
             
C Writing state variables to files      
	do i=1, NBIOTS
	
             if (ilhkv(BIOTSNOD(i)).le.1) then
              if(what.eq.'wc') write(BIOTSFUN(i),5010) it,     !different format for WC and BS
     *                           (e(1,BIOTSNOD(i),j),j=1,noutput)
              if(what.eq.'bs') write(BIOTSFUN(i),5012) it,
     *                           (e(1,BIOTSNOD(i),j),j=1,noutput)
             else
              do k=1,ilhkv(BIOTSNOD(i))
               if(what.eq.'wc') write(BIOTSFUN(i),5011) it,k,
     *                               (e(k,BIOTSNOD(i),j),j=1,noutput)
               if(what.eq.'bs') write(BIOTSFUN(i),5013) it,k,
     *                               (e(k,BIOTSNOD(i),j),j=1,noutput)
              enddo
             endif

C Write to the standard output        
          write(6,*) 'biotser_write: for node ',i,' written at ',it

          if((it+idtcon).gt.itend) then
	        close(BIOTSFUN(i))
	        iend = 1
          end if
	    
	end do
	
      if (initial.eq.1) return ! Writing initial conditions only


C Writing diagnostics (auxilary variables, rate components) to files      
      
      do istate=1,nstate
       
	 do i=1, NDGTS(istate)
	 
	  if (ilhkv(DGTSNOD(i,istate)).le.1) then             
          if(what.eq.'wc') write(DGTSFUN(i,istate), 5020) it, 
     +        (dg(1,DGTSNOD(i,istate),istate,j),j=1,NDIAGVAR)
          if(what.eq.'bs') write(DGTSFUN(i,istate), 5022) it, 
     +        (dg(1,DGTSNOD(i,istate),istate,j),j=1,NDIAGVAR)

        else 
         do k=1,ilhkv(DGTSNOD(i,istate)) 
          if(what.eq.'wc') write(DGTSFUN(i,istate),5021) it,k,
     *              (dg(k,DGTSNOD(i,istate),istate,j),
     *                                  j=1,NDIAGVAR)
          if(what.eq.'wc') write(DGTSFUN(i,istate),5023) it,k,
     *              (dg(k,DGTSNOD(i,istate),istate,j),
     *                                  j=1,NDIAGVAR)
         enddo
        endif 
  
C Write to the standard output
         write(6,*) 'biotser_write: STATE VAR. ',istate,
     +     ' auxilaries for node ',i,
     +	           ' written at ',it         
         
         if((it+idtcon).gt.itend) then
	        close(DGTSFUN(i,istate))
	        iend = 1
         end if

	 end do
	end do

C state variables
 5010 FORMAT(I15,25F10.4)    ! for water column 1 layer
 5011 FORMAT(I15,I5,25F10.4) ! for water column many layers

 5012 FORMAT(I15,25F10.4)    !for bottom sediments 1 layer
 5013 FORMAT(I15,I5,25F10.4) !for bottom sediments many layers

C diagnostics (processes)
 5020 FORMAT(I15,20E10.3)   ! for water column 1 layer
 5021 FORMAT(I15,I5,20E10.3)! for water column many layers
  
 5022 FORMAT(I15,20E10.3)   !for bottom sediments 1 layer 
 5023 FORMAT(I15,I5,20E10.3)!for bottom sediments many layers
  
       return	
	end !biotser_write
	
c********************************************************************

c********************************************************************
c********************************************************************

c
	subroutine cur_wmeteo(tempair,windspeed)

c sets meteo parameters just fort he case when meteorological
c information is spatialy heterogeneous
c
c tempair is the air temperature in degrees C
c windspeed is the wind speed in m/s

	implicit none


C COPIED BY PETRAS FROM  subfx.h TO GET airtemp AND wind, 18.08.2004
C Wind spead is spatialy homgeneous in this version and might be
C slightly different from used for hydrodynamics.
     
      real qsact,taact,tbact,uwact,ccact
     +			,uract,pact,eact,ract,qact
     
C
c	taact  - air temperature, C
c	uwact  - wind speed, m/s
c   qsact  - solar radiation, W/m2
c   tbact  - wet bulb temperature, C
c   ccact  - cloud cover (0 clear sky, 1 totally covered)
	

	real tempair
	real windspeed


  
	
      tempair = taact      
      windspeed = uwact
C
	end

c********************************************************************
c********************************************************************
c********************************************************************
c********************************************************************

      subroutine cur_param_read(par)
      
C lun lug  3 10:00:57 CEST 2006
c changed from Christian and Michol
c 
c reads parameters from file with biocon flag in control file *.str  
c setting default parameters 
c
C---------------------------------------------------------------
C     DEFAULT VALUES FOR MODEL CONSTANTS ARE INITIALIZED
C---------------------------------------------------------------
C     SUBROUTINE MODIFIED BY ALI AND PETRAS, CORPI
C     TO READ MODEL CONSTANTS FROM FILE 
C     29 July 2004
C
C     2006 modified by Petras Michol and Christian to read as SHYFEM par file
C     ------------
C           
C     INITIALIZATION OF AQUABC WC MODEL CONSTANTS, PARAMETERS AND
C     KINETIC TIME FUNCTIONS WITH DEFAULT VALUES
C



      implicit none

      integer nrdnls,iw,ioff
      integer ifileo
      integer iunit
      double precision value
      integer nc              !number parameters
      parameter(nc = 800)
      character*80 file,name,text,line,bname(nc)
      integer j
	integer ifile
      

  
C     ARRAY FOR MODEL CONSTANTS
      real par(nc)

C     MODEL CONSTANTS DATA TYPE
      integer ctype

C     GENERAL PURPOSE COUNTER FOR ARRAY INDICES, ....
      integer i

C     INITIALIZE ARRAYS BEFORE READING
     
      do i=1, nc
          par(i)   = 0.0
          bname(i) = ' '
      end do

C     INITIALIZE ARRAYS WITH DEFAULT VALUES
     
C WATER COLUMN KINETICS PARAMETERS
C     1 = small, 2 = medium, 3 = large
      bname(1)  = 'WTYPE' 
      par(1) = 3.0
C      WCID(1) = 1

C     kcnit  NH4->NO3 Nitrification rate at 20C, day-1
      bname(11) = 'K1320C'
      par(11) = 0.3
C      WCID(11) = 11

C     NH4->NO3 Nitrification rate temp. const.
      bname(12) = 'K1320T'
      par(12) = 1.08
C      WCID(12) = 12

C     NH4->NO3 Half sat. for nitrific., mgO2/l
      bname(13) = 'knit'
      par(13) = 2.0
C      WCID(13) = 13

C     NO3-> Denitrification rate at 20C, day-1day-1
      bname(21) = 'K140C'
      par(21) = 0.1
C      WCID(21) = 21

C     NO3-> Denitrification rate temp. const.unitless
      bname(22) = 'K140T'
      par(22) = 1.045
C      WCID(22) = 22

C     NO3-> Half sat. for denitrific., mgO2/L
      bname(23) = 'KNO3'
      par(23) = 0.1
C      WCID(23) = 23

C     Greens growth rate constant at 20C, day-1
      bname(41) = 'k1c'
      par(41) = 8.0 
C      WCID(41) = 41

C     Diatoms growth rate constant at 20C, day-1
      bname(241) = 'k1c_2'
      par(241) = 6.0
C      WCID(241) = 241

C     Cyanobacteria growth rate constant at 20C, day-1
      bname(341) = 'k1c_3'
      par(341) = 7.0 
C      WCID(341) = 341

C     PHY   Carbon half saturation, mgC/l
      bname(42) = 'kc'
      par(42) = 0.3
C      WCID(42) = 42

C     PHY_2 Carbon half saturation, mgC/l
      bname(242) = 'kc_2'
      par(242) = 0.3
C      WCID(242) = 242

C     PHY_3 Carbon half saturation, mgC/l
      bname(342) = 'kc_3'
      par(342) = 0.3
C      WCID(342) = 342 

C     1Di Toro 2Smith  3 Steele
      bname(43) = 'LGHTSW'
      par(43) = 2.0
C      WCID(43) = 43

C     PHY   Quantum yield const. mg C/mole photon, LGHTSW2
      bname(44) = 'PHIMX'
      par(44) = 720
C      WCID(44) = 44

C     PHY_2 Quantum yield const. mg C/mole photon, LGHTSW2
      bname(244) = 'PHIMX_2'
      par(244) = 720
C      WCID(244) = 244

C     PHY_3 Quantum yield const. mg C/mole photon, LGHTSW2
      bname(344) = 'PHIMX_3'
      par(344) = 720
C      WCID(344) = 344

C     Chloroph. extinction, ( mcg Chla/l/m, LGHTSW2
      bname(45) = 'XKC'
      par(45) = 0.011
C      WCID(45) = 45

C     PHYT   Carbon to chlorophyl ratio
      bname(46) = 'CCHL'
      par(46) = 30.0
C      WCID(46) = 46

C     PHY_2  Carbon to chlorophyl ratio
      bname(246) = 'CCHL_2'
      par(246) = 30.0
C      WCID(246) = 246

C     PHYT_3 Carbon to chlorophyl ratio
      bname(346) = 'CCHL_3'
      par(346) = 30.0
C      WCID(346) = 346

C     PHY   Nitrogen half saturation, mgN/l
      bname(48) = 'KMNG1'
      par(48) = 0.1
C      WCID(48) = 48

C     PHY_2 Nitrogen half saturation, mgN/l
      bname(248) = 'KMNG1_2'
      par(248) = 0.1
C      WCID(248) = 248

C     PHY_3 Nitrogen half saturation, mgN/l
      bname(348) = 'KMNG1_3'
      par(348) = 0.1
C      WCID(348) = 348

C     PHY Phosphorus half saturation for phyto.
      bname(49) = 'KMPG1'
      par(49) = 0.005
C      WCID(49) = 49

C     PHY_2 Phosphorus half saturation for phyto.
      bname(249) = 'KMPG1_2'
      par(249) = 0.005
C      WCID(249) = 249

C     PHY_3 Phosphorus half saturation for phyto.
      bname(349) = 'KMPG1_3'
      par(349) = 0.004
C      WCID(349) = 349

C     Silica half saturation coefficient for phyto
      bname(51) = 'KMNSI'
      par(51) = 0.1
C      WCID(51) = 51

C     PHY Algal respiration rate
      bname(50) = 'k1rc'
      par(50) = 0.01
C      WCID(50) = 50

C     PHY_2 Algal respiration rate
      bname(250) = 'k1rc_2'
      par(250) = 0.01
C      WCID(250) = 250

C     PHY_3 Algal respiration rate
      bname(350) = 'k1rc_3'
      par(350) = 0.01
C      WCID(350) = 350

C     PHY     Phytoplankton death rate,day-1
      bname(52) = 'k1d'
      par(52) = 0.18
C      WCID(52) = 52

C     PHY_2   Phytoplankton death rate,day-1
      bname(252) = 'k1d_2'
      par(252) = 0.30
C      WCID(252) = 252

C     PHY_3   Phytoplankton death rate,day-1
      bname(352) = 'k1d_3'
      par(352) = 0.22
C      WCID(352) = 352

C     PHYT Nutrient lim. option ( 1-min, 0-product
      bname(54) = 'NUTLIM'
      par(54) = 1.0
C      WCID(54) = 54

C     PHY   Phosporus to Carbon Ratio - - PCRBOP,  mgP/mgC
      bname(57) = 'PCRB'
      par(57) = 0.024
C      WCID(57) = 57

C     PHY_2 Phosporus to Carbon Ratio - - PCRBOP,  mgP/mgC
      bname(257) = 'PCRB_2'
      par(257) = 0.024
C      WCID(257) = 257

C     PHY_3 Phosporus to Carbon Ratio - - PCRBOP,  mgP/mgC
      bname(357) = 'PCRB_3'
      par(357) = 0.024
C      WCID(357) = 357

C     PHY   Nitrogen to Carbon Ratio,  mgP/mgC
      bname(58) = 'NCRB'
      par(58) = 0.176
C      WCID(58) = 58

C     PHY_2 Nitrogen to Carbon Ratio,  mgP/mgC
      bname(258) = 'NCRB_2'
      par(258) = 0.176
C      WCID(258) = 258

C     PHY_3 Nitrogen to Carbon Ratio,  mgP/mgC
      bname(358) = 'NCRB_3'
      par(358) = 0.176
C      WCID(358) = 358

C     PHY   fraction of photorespiration in primary production
      bname(60) = 'fpr'
      par(60) = 0.1
C      WCID(60) = 60

C     PHY_2 fraction of photorespiration in primary production
      bname(260) = 'fpr_2'
      par(260) = 0.1
C      WCID(260) = 260

C     PHY_3 fraction of photorespiration in primary production
      bname(360) = 'fpr_3'
      par(360) = 0.1
C      WCID(360) = 360

C     PHY  excretion rate
      bname(61) = 'kexcr'
      par(61) = 0.01
C      WCID(61) = 61

C     PHY_2 excretion rate
      bname(261) = 'kexcr_2'
      par(261) = 0.02
C      WCID(261) = 261

C     PHY_3 excretion rate
      bname(361) = 'kexcr_3'
      par(361) = 0.02
C      WCID(361) = 361

C     PHY  excretion and respiration temperature coefficient
      bname(62) = 'kert'
      par(62) = 1.045
C      WCID(62) = 62

C     PHY_2 excretion and respiration temperature coefficient
      bname(262) = 'kert_2'
      par(262) = 1.045
C      WCID(262) = 262

C     PHY_3 excretion and respiration temperature coefficient
      bname(362) = 'kert_3'
      par(362) = 1.045
C      WCID(362) = 362

C     PHY  death temperature coefficient
      bname(63) = 'kdt'
      par(63) = 1.05
C      WCID(63) = 63

C     PHY_2 death temperature coefficient
      bname(263) = 'kdt_2'
      par(263) = 1.05
C      WCID(263) = 263

C     PHY_3 death temperature coefficient
      bname(363) = 'kdt_3'
      par(363) = 1.05
C      WCID(363) = 363

C     Ratio between N fixing and nonfixing growth rate
      bname(65) = 'R_fix'
      par(65) = 0.8
C      WCID(65) = 65

C     Effectivity parameter of switching to N fixation
      bname(66) = 'K_fix'
      par(66) = 0.02
C      WCID(66) = 66

C     ZOO grazing Preference coefficient for greens
      bname(67) = 'p_1'
      par(67) = 0.30
C      WCID(67) = 67

C     ZOO Grazing Preference coefficient for diatoms
      bname(267) = 'p_2'
      par(267) = 0.30
C      WCID(267) = 267

C     Grazing Preference coefficient for cyanobacteria
      bname(367) = 'p_3'
      par(367) = 0.0
C      WCID(367) = 367

C     Grazing Preference coefficient for detritus
      bname(467) = 'p_4'
      par(467) = 0.40
C      WCID(467) = 467
      bname(468) = 'p_5'
      par(468)   =   0.15     !468 Grazing Preference coefficient for greens based part. det. carbon
      bname(469) = 'p_6'
      par(469)   =   0.15     !469 Grazing Preference coefficient for diatoms based part. det. carbon
      bname(470) = 'p_7'
      par(470)   =   0.35     !470 Grazing Preference coefficient for cyanobacteria based part. det. carbon
C     Minimum food concentration for zoo
      bname(68) = 'FOOD_min'
      par(68) = 0.0
C      WCID(68) = 68

C     Rearation rate at 20C, day-1,if k2=0 then use kawind, kahydra
      bname(82) = 'K2'
      par(82) = 0
C      WCID(82) = 82

C     PHYT light extinction coef. m-1
      bname(83) = 'ke'
      par(83) = 1.0
C      WCID(83) = 83

C     pHminm Mineralization  rate pH multiplier parameter
      bname(105) = 'pHminm'
      par(105) = 5.0
C      WCID(105) = 105

C     pHmaxm Mineralization  rate pH multiplier parameter
      bname(106) = 'pHmaxm'
      par(106) = 8.8
C      WCID(106) = 106

C     pHminn Nitrification   rate pH multiplier parameter
      bname(107) = 'pHminn'
      par(107) = 7.5
C      WCID(107) = 107

C     pHmaxn Nitrification   rate pH multiplier parameter
      bname(108) = 'pHmaxn'
      par(108) =9.0
C      WCID(108) = 108

C     pHmindn Denitrification rate pH multiplier parameter
      bname(109) = 'pHmindn'
      par(109) = 5.0
C      WCID(109) = 109

C     pHmaxdn Denitrification rate pH multiplier parameter
      bname(110) = 'pHmaxdn'
      par(110) = 9.0
C      WCID(110) = 110

C     PHYT->ZOO Phytoplankton grazing rate by Zoo  0.32
      bname(112) = 'kgrz'
      par(112) = 0.32
C      WCID(112) = 112

C     PHYT->ZOO half saturation constant for grazing
      bname(113) = 'KPHYZ'
      par(113) = 0.1
C      WCID(113) = 113

C     PHYT->ZOO zoo-phyto digestion efficiency
      bname(114) = 'eff'
      par(114) = 0.8
C      WCID(114) = 114

C     ZOO death
      bname(115) = 'kdz'
      par(115) = 0.12
C      WCID(115) = 115

C     ZOO respiration rate
      bname(120) = 'krz'
      par(120) = 0.03
C      WCID(120) = 120

C     ZOO excretion fraction in respiration
      bname(121) = 'kexcrz'
      par(121) = 0.02
C      WCID(121) = 121

C     Silica to carbon ratio ( asc
      bname(211) = 'asc'
      par(211) = 0.25
C      WCID(211) = 211

C     ZOO growth Q10 factor
      bname(130) = 'q10z'
      par(130) = 2.0
C      WCID(130) = 130

C     ZOO growth maximal temperature
      bname(131) = 'tmaxz'
      par(131) = 30.0
C      WCID(131) = 131

C     ZOO growth optimal temperature
      bname(132) = 'toptz'
      par(132) = 18.0
C      WCID(132) = 132

C     ZOO growth maximal acclimation temperature( delay
      bname(133) = 'tamaxz'
      par(133) = 0.0
C      WCID(133) = 133

C     ZOO growth maximal acclimation approaching rate
      bname(134) = 'kamaxz'
      par(134) = 0.0
C      WCID(134) = 134

C     ZOO growth temperature below wich no acclimation
      bname(135) = 'taminz'
      par(135) = 0.0
C      WCID(135) = 135

C     ZOO respiration temperature coefficient
      bname(136) = 'krtz'
      par(136) = 1.045
C      WCID(136) = 136

C     ZOO death temperature coefficient
      bname(137) = 'kdtz'
      par(137) = 1.05
C      WCID(137) = 137

C     PHY growth Q10 factor
      bname(180) = 'q10'
      par(180) = 2.1
C      WCID(180) = 180

C     PHY growth maximal temperature
      bname(181) = 'tmax'
      par(181) = 32.0
C      WCID(181) = 181

C     PHY growth optimal temperature
      bname(182) = 'topt'
      par(182) = 20.0
C      WCID(182) = 182

C     PHY growth maximal acclimation temperature( delay
      bname(183) = 'tamax'
      par(183) = 0.0
C      WCID(183) = 183

C     PHY growth maximal acclimation temperature( delay
      bname(184) = 'kamax'
      par(184) = 0.0
C      WCID(184) = 184

C     PHY growth temperature below wich no acclimation
      bname(185) = 'tamin'
      par(185) = 0.0
C      WCID(185) = 185

C     PHY_2 growth Q10 factor
      bname(280) = 'q10_2'
      par(280) = 2.5
C      WCID(280) = 280

C     PHY_2 growth maximal temperature
      bname(281) = 'tmax_2'
      par(281) = 20.0
C      WCID(281) = 281

C     PHY_2 growth optimal temperature
      bname(282) = 'topt_2'
      par(282) = 12.0
C      WCID(282) = 282

C     PHY_2 growth maximal acclimation temperature( delay
      bname(283) = 'tamax_2'
      par(283) = 0.0
C      WCID(283) = 283

C     PHY_2 growth maximal acclimation approaching rate
      bname(284) = 'kamax_2'
      par(284) = 0.0
C      WCID(284) = 284

C     PHY_2 growth temperature below wich no acclimation
      bname(285) = 'tamin_2'
      par(285) = 0.0
C      WCID(285) = 285

C     PHY_3 growth Q10 factor
      bname(380) = 'q10_3'
      par(380) = 2.5
C      WCID(380) = 380

C     PHY_3 growth maximal temperature
      bname(381) = 'tmax_3'
      par(381) = 32.0
C      WCID(381) = 381

C     PHY_3 growth optimal temperature
      bname(382) = 'topt_3'
      par(382) = 20.0
C      WCID(382) = 382

C     PHY_3 growth maximal acclimation temperature( delay
      bname(383) = 'tamax_3'
      par(383) = 0.0
C      WCID(383) = 383

C     PHY_3 growth maximal acclimation approaching rate
      bname(384) = 'kamax_3'
      par(384) = 0.0
C      WCID(384) = 384

C     PHY_3 growth temperature below wich no acclimation
      bname(385) = 'tamin_3'
      par(385) = 0.0
C      WCID(385) = 385

C     Rate constant for ZOOPPDETC dissolution
      bname(501) = 'C_DISS_ZOOPPDETC'
      par(501) = 0.7
C      WCID(501) = 501

C     Temperature correction for ZOOPPDETC dissolution
      bname(502) = 'T_DISS_ZOOPPDETC'
      par(502) = 1.03
C      WCID(502) = 502

C     Rate constant for GPHYPDETC dissolution
      bname(503) = 'C_DISS_GPHYPDETC'
      par(503) = 0.70
C      WCID(503) = 503

C     Temperature correction for GPHYPDETC dissolution
      bname(504) = 'T_DISS_GPHYPDETC'
      par(504) = 1.03
C      WCID(504) = 504

C     Rate constant for DPHYPDETC dissolution
      bname(505) = 'C_DISS_DPHYPDETC'
      par(505) = 0.70
C      WCID(505) = 505

C     Temperature correction for DPHYPDETC dissolution
      bname(506) = 'T_DISS_DPHYPDETC'
      par(506) = 1.03
C      WCID(506) = 506

C     Rate constant for CPHYPDETC dissolution
      bname(507) = 'C_DISS_CPHYPDETC'
      par(507) = 0.70
C      WCID(507) = 507

C     Temperature correction for CPHYPDETC dissolution
      bname(508) = 'T_DISS_CPHYPDETC'
      par(508) = 1.03
C      WCID(508) = 508

C     Rate constant for EXLAPDETC dissolution
      bname(509) = 'C_DISS_EXLAPDETC'
      par(509) = 0.10
C      WCID(509) = 509

C     Temperature correction for EXLAPDETC dissolution
      bname(510) = 'T_DISS_EXLAPDETC'
      par(510) = 1.03
C      WCID(510) = 510

C     N:C ratio for EXLAPDETC
      bname(511) = 'NC_EXLAPDETC'
      par(511) = 0.09
C      WCID(511) = 511

C     P:C ratio for EXLAPDETC
      bname(512) = 'PC_EXLAPDETC'
      par(512) = 0.03
C      WCID(512) = 512

C     Si:C ratio for EXLAPDETC
      bname(513) = 'SiC_EXLAPDETC'
      par(513) = 0.03
C      WCID(513) = 513

C     Rate constant for EXREPDETC dissolution
      bname(514) = 'C_DISS_EXREPDETC'
      par(514) = 0.01
C      WCID(514) = 514

C     Temperature correction for EXREPDETC dissolution
      bname(515) = 'T_DISS_EXREPDETC'
      par(515) = 1.03
C      WCID(515) = 515

C     N:C ratio for EXREPDETC
      bname(516) = 'NC_EXREPDETC'
      par(516) = 0.09
C      WCID(516) = 516

C     P:C ratio for EXREPDETC
      bname(517) = 'PC_EXREPDETC'
      par(517) = 0.03
C      WCID(517) = 517

C     Si:C ratio for EXREPDETC
      bname(518) = 'SiC_EXREPDETC'
      par(518) = 0.03
C      WCID(518) = 518

C     Rate constant for ZOOPDDETC oxidation
      bname(519) = 'C_OX_ZOOPDDETC'
      par(519) = 0.70
C      WCID(519) = 519

C     Temprature correction factor for ZOOPDDETC oxidation
      bname(520) = 'T_OX_ZOOPDDETC'
      par(520) = 1.03
C      WCID(520) = 520

C     Minimum pH for ZOOPDDETC oxidation
      bname(521) = 'pHminm_OX_ZOOPDDETC'
      par(521) = 5.0
C      WCID(521) = 521

C     Maximum pH for ZOOPDDETC oxidation
      bname(522) = 'pHmaxm_OX_ZOOPDDETC'
      par(522) = 9.0
C      WCID(522) = 522 

C     Half saturation for oxygen for ZOOPDDETC oxidation
      bname(523) = 'k_OX_ZOOPDDETC'
      par(523) = 2.0
C      WCID(523) = 523

C     O:C ratio for ZOOPDDETC
      bname(524) = 'OC_ZOOPDDETC'
      par(524) = 2.6
C      WCID(524) = 524

C     Rate constant for GPHYDDETC oxidation
      bname(525) = 'C_OX_GPHYDDETC'
      par(525) = 0.70
C      WCID(525) = 525

C     Temprature correction factor for GPHYDDETC oxidation
      bname(526) = 'T_OX_GPHYDDETC'
      par(526) = 1.03
C      WCID(526) = 526

C     Minimum pH for GPHYDDETC oxidation
      bname(527) = 'pHminm_OX_GPHYDDETC'
      par(527) = 5.0
C      WCID(527) = 527

C     Maximum pH for GPHYDDETC oxidation
      bname(528) = 'pHmaxm_OX_GPHYDDETC'
      par(528) = 9.0
C      WCID(528) = 528

C     Half saturation for oxygen for GPHYDDETC oxidation
      bname(529) = 'k_OX_GPHYDDETC'
      par(529) = 2.0
C      WCID(529) = 529

C     O:C ratio for GPHYDDETC
      bname(530) = 'OC_GPHYDDETC'
      par(530) = 2.6
C      WCID(530) = 530

C     Rate constant for DPHYDDETC oxidation
      bname(531) = 'C_OX_DPHYDDETC'
      par(531) = 0.70
C      WCID(531) = 531

C     Temprature correction factor for DPHYDDETC oxidation
      bname(532) = 'T_OX_DPHYDDETC'
      par(532) = 1.03
C      WCID(532) = 532

C     Minimum pH for DPHYDDETC oxidation
      bname(533) = 'pHminm_OX_DPHYDDETC'
      par(533) =  5.0
C      WCID(533) = 533

C     Maximum pH for DPHYDDETC oxidation
      bname(534) = 'pHmaxm_OX_DPHYDDETC'
      par(534) = 9.0
C      WCID(534) = 534

C     Half saturation for oxygen for DPHYDDETC oxidation
      bname(535) = 'k_OX_DPHYDDETC'
      par(535) = 2.0
C      WCID(535) = 535

C     O:C ratio for DPHYDDETC
      bname(536) = 'OC_DPHYDDETC'
      par(536) = 2.6
C      WCID(536) = 536

C     Rate constant for CPHYDDETC oxidation  
      bname(537) = 'C_OX_CPHYDDETC'
      par(537) = 0.70
C      WCID(537) = 537

C     Temprature correction factor for CPHYDDETC oxidation
      bname(538) = 'T_OX_CPHYDDETC'
      par(538) = 1.03
C      WCID(538) = 538

C     Minimum pH for CPHYDDETC oxidation
      bname(539) = 'pHminm_OX_CPHYDDETC'
      par(539) = 5.0
C      WCID(539) = 539

C     Maximum pH for CPHYDDETC oxidation
      bname(540) = 'pHmaxm_OX_CPHYDDETC'
      par(540) = 9.0
C      WCID(540) = 540

C     Half saturation for oxygen for CPHYDDETC oxidation
      bname(541) = 'k_OX_CPHYDDETC'
      par(541) = 2.0
C      WCID(541) = 541

C     O:C ratio for CPHYDDETC
      bname(542) = 'OC_CPHYDDETC'
      par(542) = 2.6
C      WCID(542) = 542

C     Rate constant for EXLADDETC oxidation
      bname(543) = 'C_OX_EXLADDETC'
      par(543) = 0.7
C      WCID(543) = 543

C     Temprature correction factor for EXLADDETC oxidation
      bname(544) = 'T_OX_EXLADDETC'
      par(544) = 1.03
C      WCID(544) = 544

C     Minimum pH for EXLADDETC oxidation
      bname(545) = 'pHminm_OX_EXLADDETC'
      par(545) = 5.0
C      WCID(545) = 545

C     Maximum pH for EXLADDETC oxidation
      bname(546) = 'pHmaxm_OX_EXLADDETC'
      par(546) = 9.0
C      WCID(546) = 546

C     Half saturation for oxygen for EXLADDETC oxidation
      bname(547) = 'k_OX_EXLADDETC'
      par(547) = 2.0
C      WCID(547) = 547

C     O:C ratio for EXLADDETC
      bname(548) = 'OC_EXLADDETC'
      par(548) = 2.6
C      WCID(548) = 548

C     Rate constant for EXREDDETC oxidation
      bname(549) = 'C_OX_EXREDDETC'
      par(549) = 0.03
C      WCID(549) = 549

C     Temprature correction factor for EXREDDETC oxidation
      bname(550) = 'T_OX_EXREDDETC'
      par(550) = 1.03
C      WCID(550) = 550

C     Minimum pH for EXREDDETC oxidation
      bname(551) = 'pHminm_OX_EXREDDETC'
      par(551) = 5.0
C      WCID(551) = 551

C     Maximum pH for EXREDDETC oxidation
      bname(552) = 'pHmaxm_OX_EXREDDETC'
      par(552) = 9.0
C      WCID(552) = 552

C     Half saturation for oxygen for EXREDDETC oxidation
      bname(553) = 'k_OX_EXREDDETC'
      par(553) = 2.0
C      WCID(553) = 553

C     O:C ratio for EXREDDETC
      bname(554) = 'OC_EXREDDETC'
      par(554) = 2.6
C      WCID(554) = 554

C    Temporary used for bottom sediment management(0 - no BS, 1-BS)
      bname(790) =   'bsedim'


C Making all values zero to control reading
c         do i=1, nc
c          par(i) = 0.0          
c        end do
        

c       --------------------------------------------------------
c       Get biocon file name
c       --------------------------------------------------------
c
        call getfnm('biocon',file)
        write(*,*)
c
c ----------------- reading constants file--------------------
        if( file .ne. ' ' ) then
           iunit = ifileo(0,file,'form','old')
	   call trimline(file,ifile)
           write(*,*)'Constants initialized from file: ',file(1:ifile)
           if( iunit .le. 0 ) then
            write(6,'(a,a)') 'filename: ',file(1:ifile)
            stop 'PARAM_READ: error stop: Cannot open parameters file'
           end if
c
c         --------------------------------------------------------
c         Read first line in file
c         --------------------------------------------------------

          read(iunit,'(a)') line
          
          write(6,*) line
          
          ioff = 1
c
c         --------------------------------------------------------
c         Loop on lines
c         --------------------------------------------------------

    1     continue
            iw = nrdnls(name,value,text,line,ioff)

c   nrdnls     type of variable read :
c			-1 : error
c			 0 : end of line, nothing read
c			 1 : number variable with name
c			 2 : number variable without name
c			 3 : character variable with name
c			 4 : character variable without name
            
c            write(6,*)'iw= ', iw
            
            if( iw .le. 0 ) then
              read(iunit,'(a)',end=2) line
              
c              write(6,*) line
              
              ioff = 1
            else
c            print *,'name=',name,' value= ', value
       
              if( iw .eq. 1 ) text = ' '
              if( iw .eq. 2 ) value = 0.0D+0
              do j = 1,nc
                call triml(bname(j))
                if (name .eq. bname(j)) par(j) = value                
                 
!                 write(6,*) '%%%%%%%%%%%%%%%%%%%%%%%%%',
!     *                    'j=',j,name, bname(j),value, par(j)               
                
              end do

            end if
            goto 1
    2     continue
c-------end of loop on lines---------------------
 
        end if
c----------- end of reading constants file-------

c       --------------------------------------------------------
c       Writes constant value on the screen
c       --------------------------------------------------------



        write(6,*) 'Maximal number of parameters ',nc
        write(*,*)'Constants for water quality module:'       
        
        do j = 1,nc
          if (bname(j). ne. ' ' )write(*,44)bname(j),par(j)
        end do
        write(*,*)

44      format(3x,a20,f12.7)
        end !param_read
      

c********************************************************************
c********************************************************************
c********************************************************************
c********************************************************************

      subroutine cur_param_sed_read(par)
      
c 
c     reads parameters from file with bioscon flag in control file *.str
c     for BS model No1 
C     default values for model constants are initialized
C     Adapted for sediments from cur_param by Petras Zemlys 20010 06
C           





      implicit none

      integer nrdnls,iw,ioff
      integer ifileo
      integer iunit
      double precision value
      integer nc              !number parameters
      parameter(nc = 100)
      character*80 file,name,text,line,bname(nc)
      integer j
      integer ifile
      

  
C     ARRAY FOR MODEL CONSTANTS
      real par(nc)

C     MODEL CONSTANTS DATA TYPE
      integer ctype

C     GENERAL PURPOSE COUNTER FOR ARRAY INDICES, ....
      integer i

C     INITIALIZE ARRAYS BEFORE READING
     
      do i=1, nc
          par(i)   = 0.0
          bname(i) = ' '
      end do

C     INITIALIZE ARRAYS WITH DEFAULT VALUES
     
C Bottom sediment kinetics parameters
      bname(1 ) = 'K_OXIC_DISS_POC    ' 
        par(1 ) = 0.05   !1: Dissolution rate constant of particulate organic carbon at 20 C (aerobic) - 1/day
      bname(2 ) = 'K_ANOXIC_DISS_POC  ' 
        par(2 ) = 0.01   !2: Dissolution rate constant of particulate organic carbon at 20 C (anoxic) - 1/day
      bname(3 ) = 'THETA_DISS_POC     ' 
        par(3 ) = 1.04   !3: Temperature correction for dissolution of particulate organic carbon
      bname(4 ) = 'KHS_DISS_POC       ' 
        par(4 ) = 1.00   !4: ** Half saturation concentration of POC for dissolution
      bname(5 ) = 'K_OXIC_DISS_PON    ' 
        par(5 ) = 0.05   !5: Dissolution rate constant of particulate organic nitrogen at 20 C (aerobic) - 1/day
      bname(6 ) = 'K_ANOXIC_DISS_PON  ' 
        par(6 ) = 0.01   !6: Dissolution rate constant of particulate organic nitrogen at 20 C (anoxic) - 1/day
      bname(7 ) = 'THETA_DISS_PON     ' 
        par(7 ) = 1.04   !7: Temperature correction for dissolution of particulate organic nitrogen
      bname(8 ) = 'KHS_DISS_PON       ' 
        par(8 ) = 0.50   !8: ** Half saturation concentration of PON for dissolution   
      bname(9 ) = 'K_OXIC_DISS_POP    ' 
        par(9 ) = 0.05   !9: Dissolution rate constant of particulate organic phosphorus at 20 C (aerobic) - 1/day
      bname(10) = 'K_ANOXIC_DISS_POP  ' 
        par(10) = 0.01   !10: Dissolution rate constant of particulate organic phosphorus at 20 C (anoxic) - 1/day
      bname(11) = 'THETA_DISS_POP     ' 
        par(11) = 1.04   !11: Temperature correction for dissolution of particulate organic phosphorus
      bname(12) = 'KHS_DISS_POP       ' 
        par(12) = 0.50   !12: ** Half saturation concentration of POP for dissolution
      bname(13) = 'K_DISS_PSi         ' 
       par(13) = 0.05   !13: Dissolution rate constant of particulate silicon at 20 C (aerobic) - 1/day
      bname(14) = 'K_ANOXIC_DISS_PSi  ' 
        par(14) = 0.01   !14: Dissolution rate constant of particulate silicon at 20 C (anoxic) - 1/day
      bname(15) = 'THETA_DISS_PSi     ' 
        par(15) = 1.04   !15: Temperature correction for dissolution of particulate silicon
      bname(16) = 'KHS_DISS_PSi       ' 
        par(16) = 1.00   !16: ** Half saturation concentration of PSi for dissolution 
      bname(17) = 'K_OXIC_MINER_DOC   ' 
        par(17) = 0.10   !17: ** Mineralization rate constant of dissolved organic carbon at 20 C (aerobic) - 1/day
      bname(18) = 'K_ANOXIC_MINER_DOC ' 
        par(18) = 0.05   !18: ** Mineralization rate constant of dissolved organic carbon at 20 C (anoxic) - 1/day
      bname(19) = 'THETA_MINER_DOC    ' 
        par(19) = 1.04   !19: ** Temperature correction for dissolution of dissolved organic carbon
      bname(20) = 'KHS_MINER_DOC      ' 
        par(20) = 1.00   !20: ** Half saturation concentration of DOC for mineralization
      bname(21) = 'K_OXIC_MINER_DON   ' 
        par(21) = 0.10   !21: ** Mineralization rate constant of dissolved organic nitrogen at 20 C (aerobic) - 1/day
      bname(22) = 'K_ANOXIC_MINER_DON ' 
        par(22) = 0.05   !22: ** Mineralization rate constant of dissolved organic nitrogen at 20 C (anoxic) - 1/day
      bname(23) = 'THETA_MINER_DON    ' 
        par(23) = 1.04   !23: ** Temperature correction for dissolution of dissolved organic nitrogen
      bname(24) = 'KHS_MINER_DON      ' 
        par(24) = 1.00   !24: ** Half saturation concentration of DON for mineralization
      bname(25) = 'K_OXIC_MINER_DOP   ' 
        par(25) = 0.10   !25: ** Mineralization rate constant of dissolved organic phosphorus at 20 C (aerobic) - 1/day
      bname(26) = 'K_ANOXIC_MINER_DOP ' 
        par(26) = 0.05   !26: ** Mineralization rate constant of dissolved organic phosphorus at 20 C (anoxic) - 1/day
      bname(27) = 'THETA_MINER_DOP    ' 
        par(27) = 1.04   !27: ** Temperature correction for dissolution of dissolved organic phosphorus
      bname(28) = 'KHS_MINER_DOP      ' 
        par(28) = 1.00   !28: ** Half saturation concentration of DOP for mineralization
      bname(29) = 'O_TO_C             ' 
        par(29) = 0.2    !29: Oxygen to carbon ratio
      bname(30) = 'K_NITR             ' 
        par(30) = 0.15   !30: Nitrification rate constant at 20 C - 1/day
      bname(31) = 'THETA_NITR         ' 
        par(31) = 1.04   !31: Temperature correction for nitrification
      bname(32) = 'KHS_NITR_NH4N      ' 
        par(32) = 1.00   !32: Half saturation constant of nitrification for NH4N - mg/L N
      bname(33) = 'KHS_NITR_DOXY      ' 
        par(33) = 2.00   !33: Half saturation constant of nitrification for DOXY - mg/L O2
      bname(34) = 'K_DENITR           ' 
        par(34) = 0.08   !34: Denitrification rate constant at 20 C - 1/day
      bname(35) = 'THETA_DENITR       ' 
        par(35) = 1.04   !35: Temperature correction for denitrification
      bname(36) = 'KHS_DENITR_NO3N    ' 
        par(36) = 1.00   !36: Half saturation constant of denitrification for NO3N - mg/L N
      bname(37) = 'KHS_DENITR_DOC     ' 
        par(37) = 1.00   !37: Half saturation constant of denitrification for DOC - mg/L C
      bname(38) = 'KHS_DENITR_DOXY    ' 
        par(38) = 0.80   !38: Half saturation constant of denitrification for DOXY - mg/L O
      bname(39) = 'DENITR_YIELD       ' 
        par(39) = 10.0   !39: Denitrification yield
      bname(40) = 'DOXY_AT_ANOXIA     ' 
        par(40) = 0.50   !40: DOXY, under which anoxia begins - mg/L O2
      bname(41) = 'SOLID_PART_COEFF_NH4'
        par(41)= 3.50   !41: Solid part coeff for ammonium nitrogen (kg^-1)
      bname(42) = 'SOLID_PART_COEFF_PO4'
        par(42)= 3.50   !42: Solid part coeff for phosphate phosphorus (kg^-1) 

C Making all values zero to control reading
         do i=1, nc
          par(i) = 0.0          
        end do
        

c       --------------------------------------------------------
c       Get biocon file name
c       --------------------------------------------------------
c
        call getfnm('bioscon',file)
        write(*,*)
c
c ----------------- reading constants file--------------------
        if( file .ne. ' ' ) then
           iunit = ifileo(0,file,'form','old')
	   call trimline(file,ifile)
           write(*,*)'Constants initialized from file: ',file(1:ifile)
           if( iunit .le. 0 ) then
            write(6,'(a,a)') 'filename: ',file(1:ifile)
            stop 'PARAM_READ_SED: error stop:
     *        Cannot open parameters file'
           end if
c
c         --------------------------------------------------------
c         Read first line in file
c         --------------------------------------------------------

          read(iunit,'(a)') line
          
c          write(6,*) line
          
          ioff = 1
c
c         --------------------------------------------------------
c         Loop on lines
c         --------------------------------------------------------

    1     continue
            iw = nrdnls(name,value,text,line,ioff)

c   nrdnls     type of variable read :
c			-1 : error
c			 0 : end of line, nothing read
c			 1 : number variable with name
c			 2 : number variable without name
c			 3 : character variable with name
c			 4 : character variable without name
            
c            write(6,*)'iw= ', iw
            
            if( iw .le. 0 ) then
              read(iunit,'(a)',end=2) line
              
c              write(6,*) line
              
              ioff = 1
            else
c            print *,'name=',name,' value= ', value
       
              if( iw .eq. 1 ) text = ' '
              if( iw .eq. 2 ) value = 0.0D+0
              do j = 1,nc
              
              call triml(bname(j))
                
                if (name .eq. bname(j)) par(j) = value                
                 
!                 write(6,*) '%%%%%%%%%%%%%%%%%%%%%%%%%',
!     *                    'j=',j,name, bname(j),value, par(j)               
                
              end do

            end if
            goto 1
    2     continue
c-------end of loop on lines---------------------
 
        end if
c----------- end of reading constants file-------

c       --------------------------------------------------------
c       Writes constant value on the screen
c       --------------------------------------------------------

        write(6,*) 'Maximal number of parameters ',nc
        write(*,*)'Constants for bottom sediments module:'        

        do j = 1,nc
          if (bname(j). ne. ' ' )write(*,44)bname(j),par(j)
        end do
        write(*,*)

44      format(3x,a20,f12.7)
        end !param_sed_read
      

c********************************************************************
c********************************************************************





c********************************************************************
c********************************************************************

	subroutine aquabcini(par,par_sed,PHTIME,PHTAB,TEMPTIME,TEMPTAB,
     *                       NBIOTS, BIOTSNOD, BIOTSFUN,
     *                       NBIOTS_sed, BIOTSNOD_sed, BIOTSFUN_sed,	          
     *                       NDGTS, DGTSNOD, DGTSFUN,
     *                       NDGTS_sed,DGTSNOD_sed,DGTSFUN_sed)


c initializes eutro routines
      
	implicit none
	
	include 'aquabc.h'
	include 'aquabc_aout.h'
	
	integer ibsedim
	
      real par(800)
      real par_sed(100) 
      real  PHTAB(1000)
      INTEGER PHTIME(1000)

      real  TEMPTAB(1000)
      INTEGER TEMPTIME(1000)
      
      integer INFUN,  ifileo,nb, nb_sed 

      character*80 file
      character*80 file_sed

      
      
c    Reading of model parameters(constants) for WC

      call cur_param_read(par)
      ibsedim = par(790)
	
c     Reading of model parameters(constants) for BS	
	
      if (ibsedim.eq.1) then
       call cur_param_sed_read(par_sed)
      endif
    
C  Initialisation of state variables output
       write(6,*) '--------------------------------------------------'
     +          ,'---------------'
       write(6,*) ''

       write(6,*) 'INITIALIZING TIME SERIES AND DIAGNOSTIC ', 
     +          'ASCII OUTPUTS FOR'
     +          ,'AQUABC Water Column VARIABLES'
       write(6,*) '--------------------------------------------------'
     +          ,'---------------'
       write(6,*) ''

          
      call getfnm('bioaow',file)
      nb = ifileo(90,file,'f','old')

       print *,'Control file unit number for Water Column variables ',
     +       ' ASCII output',nb

       print *,'Control file name for Water Column variables ',
     +       ' ASCII output',file      
      
      
	call biotser_init(nb,nstate,
     *                        NBIOTSMX,NBIOTS,
     *                        BIOTSNOD,BIOTSFUN,
     *                        NDGTSMX,NDIAGVAR,
     *                        NDGTS,DGTSNOD,DGTSFUN)

       
     
     
        
       if (ibsedim.eq.1) then
       write(6,*) '--------------------------------------------------'
     +          ,'---------------'
       write(6,*) ''


       write(6,*) 'INITIALIZING TIME SERIES AND DIAGNOSTIC ', 
     +          'ASCII OUTPUTS FOR'
     +          ,'AQUABC Bottom Sediment VARIABLES'
       write(6,*) '--------------------------------------------------'
     +          ,'---------------'
       write(6,*) ''


         call getfnm('bioaos',file_sed)
         nb_sed = ifileo(90,file_sed,'f','old')
     
       print *,'Control file unit number for Botom Sediments',
     * '  variables  ASCII output',nb_sed
       print *,'Control file name for Botom Sediments variables ',
     * ' ASCII output',file_sed

 
	   call biotser_init(nb_sed,nsstate,
     *                        NBIOTSMX_sed,NBIOTS_sed,
     *                        BIOTSNOD_sed,BIOTSFUN_sed,
     *                        NDGTSMX_sed,NDIAGVAR_sed,
     *                        NDGTS_sed,DGTSNOD_sed,DGTSFUN_sed)
     
       endif 
     

	
C Initialisation of pH reading
		
      call getfnm('bioph',file)
      INFUN = ifileo(60,file,'f','u')
      
	Call READ_PH(INFUN, PHTIME, PHTAB)

	close(INFUN)
	
C Initialisation of temperature reading(temporary)
 
        call getfnm('biotemp',file)	
	INFUN = ifileo(60,file,'f','u')
	
	Call READ_TEMP(INFUN, TEMPTIME, TEMPTAB)	

	end !aquabcini

c********************************************************************
c********************************************************************
c********************************************************************
c********************************************************************     
          SUBROUTINE READ_PH(INFUN, PHTIME, PH)

          IMPLICIT NONE

          INTEGER PHSIZE
          PARAMETER(PHSIZE = 1000)

          INTEGER INFUN
          INTEGER PHTIME
          real PH


          DIMENSION PHTIME(PHSIZE)
          DIMENSION PH(PHSIZE)

          INTEGER I


          DO 1010, I = 1, PHSIZE
              PHTIME(I) = 0
              PH(I) = 0.0
 1010     CONTINUE


          I = 1


 1020     CONTINUE         

              READ(INFUN, *, ERR = 100, END = 200) 
     *             PHTIME(I), PH(I) 


              IF (I .GT. 1) THEN
                  IF (PHTIME(I) .LE. PHTIME(I - 1)) THEN
                      GOTO 101
                  END IF
              END IF


              IF ((PH(I) .LT. 0.0).OR.(PH(I) .GT. 14.0)) THEN
                  WRITE(6,*) "WRONG VALUE FOR pH"
                  GOTO 102
              END IF


              IF (I > 1000) THEN
                  GOTO 200
              END IF

              I = I + 1

          GOTO 1020
          
         

  100     CONTINUE

          WRITE(6,*) "READ_PH : File read error in pH time series file" 
          WRITE(6, FMT = "(A7, I5)") "LINE : ", I 
          STOP 


  101     CONTINUE

          WRITE(6,*) "READ_PH : Time must be greater than previous time" 
          WRITE(6, FMT = "(A7, I5)") "LINE : ", I 
          STOP


  102     CONTINUE

          WRITE(6,*) "READ_PH : pH must be between 0 and 14" 
          WRITE(6,*) "pH you entered is :", PH(I) 
          WRITE(6, FMT = "(A7, I5)") "LINE : ", I 
          STOP


  200     CONTINUE

          WRITE(6,*) "Finished reading the pH time series file"


          IF (I < PHSIZE) THEN

              PH(I) = -1

          END IF


      END

c********************************************************************
c********************************************************************
c********************************************************************
c********************************************************************      
 


          SUBROUTINE READ_TEMP(INFUN, TEMPTIME, TEMP)

          IMPLICIT NONE

          INTEGER TEMPSIZE
          PARAMETER(TEMPSIZE = 1000)

          INTEGER INFUN
          INTEGER TEMPTIME
          real TEMP


          DIMENSION TEMPTIME(TEMPSIZE)
          DIMENSION TEMP(TEMPSIZE)

          INTEGER I


          DO 1010, I = 1, TEMPSIZE
              TEMPTIME(I) = 0
              TEMP(I) = 0.0
 1010     CONTINUE


          I = 1


 1020     CONTINUE         

              READ(INFUN, *, ERR = 100, END = 200) 
     *             TEMPTIME(I), TEMP(I) 


              IF (I .GT. 1) THEN
                  IF (TEMPTIME(I) .LE. TEMPTIME(I - 1)) THEN
                      GOTO 101
                  END IF
              END IF


              IF ((TEMP(I) .LT. -2.0).OR.(TEMP(I) .GT. 65.0)) THEN
                  WRITE(6,*) "WRONG VALUE FOR TEMPERATURE"
                  GOTO 102
              END IF


              IF (I > 1000) THEN
                  GOTO 200
              END IF

              I = I + 1

          GOTO 1020
          
          

  100     CONTINUE

          WRITE(6,*) "READ_TEMP : File read error in temperature file" 
          WRITE(6, FMT = "(A7, I5)") "LINE : ", I 
          STOP 


  101     CONTINUE

          WRITE(6,*) "READ_TEMP : ",
     *    "Time must be greater than previous time" 
          
          WRITE(6, FMT = "(A7, I5)") "LINE : ", I 
          STOP


  102     CONTINUE

          WRITE(6,*) "READ_TEMP : ",
     *    "Temparature must be between -2.0 and 65.0" 
          
          WRITE(6,*) "Temperature you entered is :", TEMP(I) 
          WRITE(6, FMT = "(A7, I5)") "LINE : ", I 
          STOP


  200     CONTINUE

          WRITE(6,*) "Finished reading the temperature time file"

  
          IF (I < TEMPSIZE) THEN

              TEMP(I) = -1

          END IF


      END


c********************************************************************
c********************************************************************
c********************************************************************
c********************************************************************
     


      real FUNCTION GET_PH(TIME, PHTIME, PH)

          IMPLICIT NONE

          INTEGER PHSIZE
          PARAMETER(PHSIZE = 1000)

          INTEGER TIME
          INTEGER PHTIME
          real PH


          DIMENSION PHTIME(PHSIZE)
          DIMENSION PH(PHSIZE)

          INTEGER I
          INTEGER K
          INTEGER P_IND
          INTEGER N_IND
          INTEGER L_IND
          LOGICAL E_IND

          real PREV_H
          real NEXT_H
          real H_PLUS
          real H_INCR
          
          E_IND = .TRUE.


          DO 1010, I = 1, 1000
              IF (PH(I) .EQ. -1) THEN

                  L_IND = I - 1
                  E_IND = .FALSE.
                  EXIT
              END IF
 1010     CONTINUE

          IF (E_IND) THEN
              L_IND = 1000
          END IF


          IF (TIME.GE.PHTIME(L_IND)) THEN
              GET_PH = PH(L_IND)
          ELSE


              IF (TIME.LE.PHTIME(1)) THEN
                  GET_PH = PH(1)
              ELSE
                  DO 1020, I = 1, L_IND - 1

                      IF (PHTIME(I).GT.TIME) THEN
                          EXIT
                      END IF

 1020             CONTINUE


                  P_IND = I - 1
                  N_IND = I

                  PREV_H = 10.0**(-PH(P_IND))
                  NEXT_H = 10.0**(-PH(N_IND))

                  H_INCR = (NEXT_H - PREV_H) / 
     *                     (PHTIME(N_IND) - PHTIME(P_IND))

                  H_PLUS = PREV_H + 
     *                     (H_INCR * (TIME - PHTIME(P_IND)))

                  GET_PH = -LOG10(H_PLUS)

              END IF

          END IF 


      END
      
c********************************************************************
c********************************************************************
c********************************************************************
c********************************************************************
 
      real FUNCTION GET_TEMP(TIME, TEMPTIME, TEMP)

          IMPLICIT NONE

          INTEGER TEMPSIZE
          PARAMETER(TEMPSIZE = 1000)

          INTEGER TIME
          INTEGER TEMPTIME
          real TEMP


          DIMENSION TEMPTIME(TEMPSIZE)
          DIMENSION TEMP(TEMPSIZE)

          INTEGER I
          INTEGER K
          INTEGER P_IND
          INTEGER N_IND
          INTEGER L_IND
          LOGICAL E_IND

          real T_INCR

          E_IND = .TRUE.


          DO 1010, I = 1, 1000
              IF (TEMP(I) .EQ. -1) THEN

                  L_IND = I - 1
                  E_IND = .FALSE.
                  EXIT
              END IF
 1010     CONTINUE


          IF (E_IND) THEN
              L_IND = 1000
          END IF


          IF (TIME.GE.TEMPTIME(L_IND)) THEN
              GET_TEMP = TEMP(L_IND)
          ELSE


              IF (TIME.LE.TEMPTIME(1)) THEN
                  GET_TEMP = TEMP(1)
              ELSE
                  DO 1020, I = 1, L_IND - 1

                      IF (TEMPTIME(I).GT.TIME) THEN
                          EXIT
                      END IF

 1020             CONTINUE


                  P_IND = I - 1
                  N_IND = I


                  T_INCR = (TEMP(N_IND) - TEMP(P_IND)) / 
     *                     (TEMPTIME(N_IND) - TEMPTIME(P_IND))

                  GET_TEMP = TEMP(P_IND) + 
     *                       (T_INCR * (TIME - TEMPTIME(P_IND)))

              END IF

          END IF 


      END
c********************************************************************
c********************************************************************
c********************************************************************
c********************************************************************   

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
c*******************************************************************

        subroutine bio_av_shell(e)

c computes and writes average/min/max of bio variables
c
c id = 260
c
c e(1) average  == 261
c e(1) min      == 262
c e(1) max      == 263
c e(2) average  == 264
c ...

        implicit none

c parameter

        include 'param.h'
        include 'aquabc.h'

        !integer nstate
        !parameter( nstate = 9 )

        real e(nlvdim,nkndim,nstate)    !state vector

c local
        integer idtc,itmc,itsmed
        integer id,nvar
c function
        real getpar
c save
        double precision bioacu(nlvdim,nkndim,nstate)
        real biomin(nlvdim,nkndim,nstate)
        real biomax(nlvdim,nkndim,nstate)

        integer ivect(8)

        save bioacu,biomin,biomax
        save ivect

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

          idtc=nint(getpar('idtcon'))
          itmc=nint(getpar('itmcon'))

          nvar = nstate

          id = 260
          call cmed_init('bav',id,nvar,nlvdim,idtc,itmc
     +                          ,bioacu,biomin,biomax,ivect)

          icall = 1
        end if

        call cmed_accum(nlvdim,e,bioacu,biomin,biomax,ivect)

        end

c*************************************************************
c*************************************************************
c*************************************************************
c here we should fill in restart routines for aquabc state variables
c*************************************************************
c*************************************************************
c*************************************************************

        subroutine write_restart_eco(iunit)
        implicit none
        integer iunit
        integer nstate,nkn,i
        nstate = 0
        nkn = 0
        write(iunit) nstate,nkn
        end
        subroutine skip_restart_eco(iunit)
        implicit none
        integer iunit
        integer nstate,nkn,i
        read(iunit) nstate,nkn
        do i=1,nstate
          read(iunit)
        end do
        end
        subroutine read_restart_eco(iunit)
        implicit none
        integer iunit
        call skip_restart_eco(iunit)
        end

c*************************************************************

