c
c $Id: ht.f,v 1.76 2010-03-11 15:36:38 georg Exp $
c
c finite element model ht (version 3D)
c
c original version from march 1991
c
c revision log :
c
c revised 30.08.95      $$AUST - austausch coefficient introduced
c revised 11.10.95      $$BCLBND - boundary condition for baroclinic runs
c revised 04.08.97      $$ZEONV - new arrays for water level elementwise
c 19.03.1998	ggu	$$IPCCV close data items commented or deleted
c 03.04.1998	ggu	$$DESCRP BUG overwriting descrp with boundary vals.
c 30.04.1998    ggu     finally eliminated /semimp/, /trock/, /ffloat/
c 05.05.1998	ggu	izeit eliminated
c 28.05.1998	ggu	call to sp131g changed
c 14.08.1998	ggu	call to biocos (biological model) introduced
c 20.08.1998	ggu	iextpo finally eliminated, momentum arrays added
c 20.08.1998	ggu	spb11 -> sp111
c 21.08.1998    ggu     xv eliminated
c 03.09.1998    ggu     biological reactor integratated
c 06.11.1998    ggu     hv renamed into hkv
c 08.04.1999    ggu     equilibrium tide (tidal potential) introduced
c 19.04.1999    ggu     subroutine custom called
c 22.10.1999    ggu     igrv, ngrv eliminated (subst by ilinkv, lenkv)
c 20.01.2000    ggu     common block /dimdim/ eliminated
c 07.03.2000    ggu     arrays for vertical turbulent diffusion coefficient
c 20.06.2000    ggu     alv -> visv    slv -> difv
c 19.11.2001    ggu     h1hv eliminated, ulv, vlv eliminated
c 03.12.2001    ggu     new arrays [cst]difhv
c 21.08.2002    ggu     new arrays /kvolc/, /ivol/ copied from hp
c 31.07.2003    ggu     array winv eliminated
c 10.08.2003    ggu     big restructuring
c 13.08.2003    ggu     some more restructuring
c 14.08.2003    ggu     even more restructuring and cleaning up
c 04.12.2003    ggu     integration of wave and sediment module
c 06 03 2004    aac     lagrangian trajectories computation module
c 03.09.2004    ggu     restart now in ht
c 15.10.2004    ggu     gotm substituted by general turbulence closure
c 02.12.2004    ggu     variable time step implemented
c 17.01.2005    ggu     new routine diff_h_set, new difhv, ausv deleted
c 24.02.2005    ggu     new routine smagorinski and common smagv
c 03.03.2005    ggu     smagorinski uses difhv, inserted in subdif.f
c 03.03.2005    ggu     new 3d boundary arrays for C/T/S
c 12.08.2005    ggu     minimum level index, TVD arrays
c 07.11.2005    ggu     new array sed2dn introduced (file name for sediments)
c 16.02.2006    ggu     new routine atoxi3d (TOXI)
c 23.03.2006    ggu     changed time step to real
c 18.10.2006    ccf     radiation stress and waves included (with pipe)
c 10.11.2006    ggu     initialize depth values after restart
c 16.11.2006    ggu     turbulence values included
c 02.04.2007    ggu     changes in algorithm (ASYM)
c 31.05.2007    deb     new arrays bpresxv, bclevvar (debora)
c 26.09.2007    ggu     deleted arrays rcv,rtv,rsv
c 27.09.2007    ggu     deleted call to tstvol,tstvol1
c 20.03.2008    acc     new call for ERSEM ecological model (BFM MODULE)
c 07.04.2008    acc     new array bfm*bc introduced (file name for ersem)
c 10.04.2008    ggu&ccf	upro, waveov, stokes, z0bk
c 16.04.2008    ggu     evaporation mass flux (evapv)
c 22.04.2008    ggu     gradx/yv non global due to parallelization
c 23.04.2008    ggu     deleted r3{c|t|s}v
c 28.04.2008    ggu     new routine init_stability(), new conzm3sh()
c 29.04.2008    ggu&aac new BMF ERSEM model
c 12.11.2008    ggu	new sigma level initialization
c 10.12.2008    ggu	new array rfricv for bottom friction
c 18.12.2008    ggu	new routine debug_output(), might be removed later
c 04.03.2009    ggu	matrix amat into solver routines
c 11.03.2009    ggu	new arrays for meteo data
c 06.04.2009    ggu	renamed nlidim to nlkdim
c 30.11.2009    ggu	write output file for successfull completion
c 19.02.2010    ggu	init_stability() changed to reset_stability()
c 22.02.2010    ggu	new arrays wxv, wyv
c 26.02.2010    ggu	new arrays sauxe1/2
c 29.04.2010    ggu	write volumes (wrfvla)
c 26.01.2011    ggu	new arrays for observations and nudging
c 16.02.2011    ggu	new iarnv, call to aquabc
c 17.02.2011    ccf	new radiation stress in 3D
c 23.03.2011    ggu	new call to adjust_spherical()
c
c*****************************************************************

c----------------------------------------------------------------

	program ht

	include 'param.h'

c----------------------------------------------------------------

	parameter (ibndim=100)
	parameter (mardim=nlvdim*10)

c variables and coefficients

	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /mkonst/ eps1,eps2,pi,flag,high,hihi
	common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft
	common /femtim/ itanf,itend,idt,nits,niter,it

	common /level/ nlvdi,nlv

c run and basin description

	character*80 descrp,descrr
	common /descrp/ descrp
	common /descrr/ descrr

c boundary file names			!$$ST	!$$DESCRP

	character*80 boundn(nbcdim)
	character*80 conzn(nbcdim)
	character*80 saltn(nbcdim)
	character*80 tempn(nbcdim)
        character*80 bio2dn(nbcdim)
        character*80 sed2dn(nbcdim)
        character*80 tox3dn(nbcdim)
	character*80 bfm1bc(nbcdim)
        character*80 bfm2bc(nbcdim)
        character*80 bfm3bc(nbcdim)

	common /boundn/ boundn
        common /conzn/ conzn
        common /saltn/ saltn
        common /tempn/ tempn
        common /bio2dn/ bio2dn
        common /sed2dn/ sed2dn
        common /tox3dn/ tox3dn
        common /bfm1bc/bfm1bc
        common /bfm2bc/bfm2bc
        common /bfm3bc/bfm3bc
  
c various arrays

	common /knausc/ knausm,knaus(nexdim)

        common /kfluxc/ nsect,kfluxm,kflux(nfxdim)
        common /iflux/ iflux(3,nfxdim)

        integer nvols,kvold,kvolm,kvol(nfxdim)
        common /kvolc/ nvols,kvold,kvolm,kvol
        integer ivolm,ivol(nfxdim)
        common /ivol/ivolm,ivol

c basin arrays

	common /nen3v/nen3v(3,neldim)
	common /iarv/iarv(neldim)
	common /ipv/ipv(nkndim), /ipev/ipev(neldim)
	common /xgv/xgv(nkndim), /ygv/ygv(nkndim)
	common /hm3v/hm3v(3,neldim)

c static geometry information

	include 'evmain.h'

        common /ilinkv/ilinkv(nkndim+1)
        common /lenkv/lenkv(nlkdim)
        common /linkv/linkv(nlkdim)

	common /ieltv/ieltv(3,neldim)
	common /kantv/kantv(2,nkndim)
	common /dxv/dxv(nkndim), /dyv/dyv(nkndim)

	common /iarnv/iarnv(nkndim)	!area code information on nodes

c dynamic geometry information

	common /iwegv/iwegv(neldim)
        common /inodv/inodv(nkndim)

c boundary arrays

	common /bnd/bnd(ibndim,nbcdim)
	common /ierv/ierv(2,nrbdim)
	common /rhv/rhv(nrbdim), /rlv/rlv(nrbdim)
	common /rrv/rrv(nrbdim), /irv/irv(nrbdim)
	common /rzv/rzv(nkndim), /rqv/rqv(nkndim)
	common /iopbnd/iopbnd(nkndim)

        real rqpsv(nkndim), rqdsv(nkndim)
        common /rqpsv/rqpsv, /rqdsv/rqdsv
        real evapv(nkndim)
        common /evapv/evapv
        real mfluxv(nlvdim,nkndim)
        common /mfluxv/mfluxv

c depth structure of levels

	common /ilhv/ilhv(neldim)
	common /ilhkv/ilhkv(nkndim)
	common /hlv/hlv(nlvdim), /hldv/hldv(nlvdim)
	common /hlhv/hlhv(neldim)

        integer ilmv(neldim)
        common /ilmv/ilmv
        integer ilmkv(nkndim)
        common /ilmkv/ilmkv

	common /hkv/hkv(nkndim), /hev/hev(neldim)

c new depth and area arrays

	real hdknv(nlvdim,nkndim)
	common /hdknv/hdknv
	real hdkov(nlvdim,nkndim)
	common /hdkov/hdkov

	real hdenv(nlvdim,neldim)
	common /hdenv/hdenv
	real hdeov(nlvdim,neldim)
	common /hdeov/hdeov

        real areakv(nlvdim,nkndim)
        common /areakv/areakv

c water level and velocity arrays

	common /zov/zov(nkndim), /znv/znv(nkndim)
	common /zeov/zeov(3,neldim), /zenv/zenv(3,neldim)	!$$ZEONV

	common /uov/uov(neldim), /vov/vov(neldim)
	common /unv/unv(neldim), /vnv/vnv(neldim)

	common /ulov/ulov(nlvdim,neldim)
	common /ulnv/ulnv(nlvdim,neldim)
	common /vlov/vlov(nlvdim,neldim)
	common /vlnv/vlnv(nlvdim,neldim)
	common /wlov/wlov(0:nlvdim,nkndim)
	common /wlnv/wlnv(0:nlvdim,nkndim)

	common /utlov/utlov(nlvdim,neldim)
	common /utlnv/utlnv(nlvdim,neldim)
	common /vtlov/vtlov(nlvdim,neldim)
	common /vtlnv/vtlnv(nlvdim,neldim)

	common /uprv/uprv(nlvdim,nkndim)
	common /vprv/vprv(nlvdim,nkndim)
	common /upro/upro(nlvdim,nkndim)
	common /vpro/vpro(nlvdim,nkndim)
	common /wprv/wprv(0:nlvdim,nkndim)

	common /up0v/up0v(nkndim)
	common /vp0v/vp0v(nkndim)
        
        common /fxv/fxv(nlvdim,neldim)		!new HYDRO debora
        common /fyv/fyv(nlvdim,neldim)

	common /xv/xv(3,nkndim)

c concentration, salinity and temperature

	common /saltv/saltv(nlvdim,nkndim)
	common /tempv/tempv(nlvdim,nkndim)
	common /sobsv/sobsv(nlvdim,nkndim)
	common /tobsv/tobsv(nlvdim,nkndim)
	common /rtauv/rtauv(nlvdim,nkndim)	!relaxation time

	common /rhov/rhov(nlvdim,nkndim)
	common /bpresv/bpresv(nlvdim,nkndim)
        common /bpresxv/bpresxv(nlvdim,neldim)	!deb110407	debora
        common /bpresyv/bpresyv(nlvdim,neldim)	!deb110407

c coriolis parameter

	common /fcorv/fcorv(neldim)

c friction and diffusion

	common /rfricv/rfricv(neldim)
	common /czv/czv(neldim)
	common /austv/austv(neldim)			!$$AUST

	common /wdifhv/wdifhv(3,3,neldim)   	!weights for horizontal diff.

	common /difhv/difhv(nlvdim,neldim)   	!horizontal diffusion - 3D

        real rdistv(nkndim)
        common /rdistv/rdistv

	common /visv/visv(0:nlvdim,nkndim)	!viscosity (momentum)
	common /difv/difv(0:nlvdim,nkndim)	!diffusivity (scalars)

c special boundary arrays

	common /ruv/ruv(nkndim), /rvv/rvv(nkndim)	!momentum input (2D)
	common /crad/crad(neldim)			!$$GWI (radiation)

c meteo (wind and pressure)

	common /tauxnv/tauxnv(nkndim), /tauynv/tauynv(nkndim)
	common /wxv/wxv(nkndim), /wyv/wyv(nkndim)
	common /wxov/wxov(nkndim), /wyov/wyov(nkndim)
	common /wxnv/wxnv(nkndim), /wynv/wynv(nkndim)

	common /ppv/ppv(nkndim)
	common /pov/pov(nkndim), /pnv/pnv(nkndim)

        real metrad(nkndim),methum(nkndim)
        real mettair(nkndim),metcc(nkndim)
        real metwbt(nkndim),metws(nkndim)
        real metrain(nkndim)
        common /metrad/metrad, /methum/methum
        common /mettair/mettair, /metcc/metcc
        common /metwbt/metwbt, /metws/metws
        common /metrain/metrain

c tidal potential

        real xgeov(nkndim), ygeov(nkndim)
        common /xgeov/xgeov, /ygeov/ygeov
        real zeqv(nkndim)
        common /zeqv/zeqv

c wave sub-module

        real waveh(nkndim)      !wave height [m]
        real wavep(nkndim)      !wave period [s]
        real waved(nkndim)      !wave direction (same as wind direction)
        real waveov(nkndim)     !orbital velocity
        real stokesx(nkndim)    !stokes velocity x
        real stokesy(nkndim)    !stokes velocity y

        common /waveh/waveh, /wavep/wavep, /waved/waved, /waveov/waveov
        common /stokesx/stokesx, /stokesy/stokesy
	save /waveh/,/wavep/,/waved/,/waveov/,/stokesx/,/stokesy/

        real z0bk(nkndim)                   !bottom roughenss on nodes
        common /z0bk/z0bk
	save /z0bk/

c variables for pipe

	integer ipipe,idcoup
	save ipipe, idcoup

c radiation stress

        real radx(nlvdim,neldim),rady(nlvdim,neldim)
        common /radx/radx,/rady/rady
	save /radx/,/rady/

c global matrix

	common /rmat/rmat(mardim)
	common /rvec/rvec(2*nlvdim)

c auxiliary arrays

	common /v1v/v1v(nkndim), /v2v/v2v(nkndim)
	common /v3v/v3v(nkndim)
	common /ve1v/ve1v(neldim)
	common /saux1/saux1(nlvdim,nkndim)
	common /saux2/saux2(nlvdim,nkndim)
	common /saux3/saux3(nlvdim,nkndim)
	common /saux4/saux4(nlvdim,nkndim)
	common /sauxe1/sauxe1(nlvdim,neldim)
	common /sauxe2/sauxe2(nlvdim,neldim)

c deleted arrays

c	real*8 urv,vrv,zrv
c	common /urv/urv(neldim), /vrv/vrv(neldim), /zrv/zrv(neldim)

	common /tken/tken(0:nlvdim,nkndim)	!turbulent kinetic energy
	common /eps/eps(0:nlvdim,nkndim)	!dissipation rate
	common /rls/rls(0:nlvdim,nkndim)	!length scale

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%%%%%%% code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	call vers3d	!set version of model

c-----------------------------------------------------------
c dimensions
c-----------------------------------------------------------

	nlvdi=nlvdim

	call cstdim(nkndim,neldim,nrbdim,nbcdim
     +			,mbwdim,ngrdim,nardim,nexdim
     +			,nfxdim,nlkdim)

c-----------------------------------------------------------
c read STR file
c-----------------------------------------------------------

	call cstinit
	call cstfile(nkndim,neldim)

c-----------------------------------------------------------
c check dimensions
c-----------------------------------------------------------

	call sp131a(nkndim,neldim)
	call sp131b(nrbdim,nbcdim)
	call sp131m(mbwdim)
	call sp131g(mardim,2,2*nlvdim,2,0)

	call cstcheck	!this also sets time and coriolis

c	call pritst(1)

c-----------------------------------------------------------
c initialize triangles
c-----------------------------------------------------------

	call set_spherical
	call set_ev
	call adjust_spherical
	call print_spherical
	call set_geom

c-----------------------------------------------------------
c initialize barene data structures
c-----------------------------------------------------------

	call setweg(-1,n)
	call setnod
	call update_geom	!update ieltv - needs inodv

c-----------------------------------------------------------
c inititialize time independent vertical arrays
c-----------------------------------------------------------

	call init_vertical

c-----------------------------------------------------------
c initialize boundary conditions
c-----------------------------------------------------------

	call sp111(1)           !here zenv, utlnv, vtlnv are initialized

	call handle_sigma_init

c-----------------------------------------------------------
c initialize depth arrays and barene data structure
c-----------------------------------------------------------

	call setweg(0,n)
	call setznv		! -> change znv since zenv has changed

        call inirst             !restart

c----------------------------------------------- FIXME
	call makehev(hev)
	call makehkv(hkv,v1v)
	call set_last_layer
c----------------------------------------------- FIXME

	call setnod

	call setarea(nlvdim,areakv)

	call setdepth(nlvdim,hdknv,hdenv,zenv,areakv)	!FIXME (see above)
	call copydepth(nlvdim,hdknv,hdkov,hdenv,hdeov)
	call setdepth(nlvdim,hdknv,hdenv,zenv,areakv)

	!call check_depth

c-----------------------------------------------------------
c initialize open boundary routines
c-----------------------------------------------------------

	call bndo_init
	call bndo_info_file('bndo_info.dat')

c-----------------------------------------------------------
c initialize transports and velocities
c-----------------------------------------------------------

	call init_uv            !set vel, w, pr, ... from transports
	call barocl(0)

c-----------------------------------------------------------
c initialize modules
c-----------------------------------------------------------

	call init_others
	call init_chezy
	call init_nodal_area_code	!area codes on nodes
        call diffweight
        call diff_h_set
	call tideini
	call cstsetup
	call sp136(ic)
        call shdist(rdistv)

c-----------------------------------------------------------
c write input values to log file and perform check
c-----------------------------------------------------------

	call check_fem
	call check_values
c	call listdim
	call prilog

c        call bclevvar_ini       !chao debora
	call bclfix_ini

	call init_pipe(ipipe,idcoup)
	call read_pipe(ipipe,it,idcoup)
        call write_pipe(ipipe,it,idcoup)

	!call custom(it)		!call for initialization

	write(6,*) 'starting time loop'

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%%%%% time loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	do while( it .lt. itend )

	   !if( it .eq. 1429200 ) call test3d(66,100)

	   call reset_stability

           call set_timestep
           call get_timestep(dt)

	   call dobefor3d

	   call sp111(2)		!boundary conditions

           call read_pipe(ipipe,it,idcoup)
	   
	   call sp259f			!hydro

	   call wrfvla(it)		!write finite volume

	   call conz3sh			!concentration (for iconz == 1)
	   call conzm3sh		!multi concentration (for iconz > 1)
	   call barocl(1)		!baroclinic contribution

           call turb_closure

           call subwaves(it,dt)         !wave model
           call sedi(it,dt)             !sediment transport

	   call ecological_module(it,dt)	!ecological model
           call atoxi3d(it,dt)			!toxi

           call lagrange

	   call doafter3d

           if( niter .eq. 1 ) call useunit(200)

	   call pritime			!output to terminal

	   call total_energy
	   call check_all
	   !call check_special

           call write_pipe(ipipe,it,idcoup)

	   !call debug_output(it)
	   !if( it .gt. 50400 ) call test3d(66,100)
	   !if( it .gt. 3600 ) call test3d(66,100)
	end do

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%% end of time loop %%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	call endtime
        !call ht_finished

	!call pripar(15)
	!call prifnm(15)

        end

c*****************************************************************

        subroutine ht_finished

        implicit none

        open(1,file='SHYFEM_FINISHED',status='unknown',form='formatted')
        write(1,*) 'SHYFEM has finished successfully the simulation'
        close(1)

	end

c*****************************************************************

	subroutine debug_output(it)

	implicit none

	integer it

	include 'param.h'

        integer nen3v(3,neldim)
        common /nen3v/nen3v
        integer ilhv(neldim)
        common /ilhv/ilhv
        real zeov(3,neldim)
        common /zeov/zeov
        real hdeov(nlvdim,neldim)
        common /hdeov/hdeov
        real visv(0:nlvdim,nkndim)
        common /visv/visv
        real utlov(nlvdim,neldim),vtlov(nlvdim,neldim)
        common /utlov/utlov, /vtlov/vtlov
	real z0bk(nkndim)
	common /z0bk/z0bk

        real saltv(nlvdim,nkndim)
        real tempv(nlvdim,nkndim)
        common /saltv/saltv
        common /tempv/tempv

	write(66) it

	write(66) 3*neldim,3
	write(66) zeov

	write(66) nlvdim*neldim,nlvdim
	write(66) hdeov

	write(66) nlvdim*neldim,nlvdim
	write(66) utlov

	write(66) nlvdim*neldim,nlvdim
	write(66) vtlov

	write(66) (nlvdim+1)*nkndim,nlvdim+1
	write(66) visv

	write(66) nkndim,1
	write(66) z0bk

	write(66) 0,0

	end

c*****************************************************************

	subroutine check_depth

	implicit none

	include 'param.h'

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	integer ilhv(neldim)
	integer ilhkv(nkndim)
	common /ilhv/ilhv
	common /ilhkv/ilhkv

	real hev(neldim)
	common /hev/hev
	real hdknv(nlvdim,nkndim)
	common /hdknv/hdknv
	real hdkov(nlvdim,nkndim)
	common /hdkov/hdkov

	real hdenv(nlvdim,neldim)
	common /hdenv/hdenv
	real hdeov(nlvdim,neldim)
	common /hdeov/hdeov

	integer k,l,lmax,ie
	real hmax

	hmax = 0.
	do ie=1,nel
	  hmax = max(hmax,hev(ie))
	end do
	write(6,*) 'hmax (hev) = ',hmax

	hmax = 0.
	do k=1,nkn
	  lmax = ilhkv(k)
	  do l=1,lmax
	    hmax = max(hmax,hdknv(l,k))
	  end do
	end do
	write(6,*) 'hmax (hdknv) = ',hmax

	end

c*****************************************************************

	subroutine check_special

	include 'param.h'

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

	integer ilhkv(nkndim)
	common /ilhkv/ilhkv
        real saltv(nlvdim,nkndim)
        real tempv(nlvdim,nkndim)
        common /saltv/saltv
        common /tempv/tempv
	
	integer k,l,lmax

	do k=1,nkn
	  lmax = ilhkv(k)
	  do l=1,lmax
	    if( saltv(l,k) .gt. 40. ) then
		write(66,*) it,l,k,saltv(l,k)
	    end if
	  end do
	end do

	end

c*****************************************************************

