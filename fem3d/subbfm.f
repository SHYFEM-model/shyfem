c
c $Id: subbfm.f,v 1.2 2008-07-16 15:41:39 georg Exp $
c
c bfm module
c
c contents :
c
c subroutine bfm_module(it,dt)
c				administers bfm ecological model
c subroutine bfm_shell(it,dt)
c				computes ecological scalars with BFM  model
c subroutine bfm_init(nbfmv1,b1cn,nbfmv2,b2cn,nbfmv3,b3cn)
c				initializes bfm  arrays
c subroutine comp_out(ivs1,ivs2,ivs3,itmbfm,idtbfm)
c				outputs variables
c subroutine setlux
c				light routine ??
c
c revision log :
c
c 10.03.2008    aac     bfm ecological module from scratch
c 29.04.2008    ggu     bfm model integrated in main branch
c 30.04.2008    ggu     double to real (BUG)
c 31.05.2011    ggu     clean from useless common blocks
c
c**************************************************************
c
c DOCS  FILENAME        Boundary conditions
c
c Boundary conditions have to be given in a file in every
c section |bound|.
c
c |bfm2dn|      File name that contains boundary conditions for
c               concentration of the bfm ecological variables.
c               The format is the same as for the file |boundn|.
c               The unit of the values given in the second
c               and following columns (nbfmv data columns for ERSEM)
c               must the ones of the variable.
c
c DOCS  FILENAME        Initial conditions
c
c Initialization of variables are done by file. The files can be created
c by the progam |laplap|. They have to be given in section |name|.
c
c**************************************************************
c
c b1cn is real in main and double precision in init ???????????
c ligth -> light	dove viene utilizzato ???
c ddepth		a che cosa serve ????
c
c in bfm Makefile:
c
c	cleanall
c
c changed:
c
c	btanf -> ibtanf
c	btend -> itbend
c	lightflag -> blight
c
c cosa fanno tutte le variabili di gotm qui???
c 
c**************************************************************
c**************************************************************

        subroutine ecological_module(it,dt)

c general interface to ecological module

        implicit none

        integer it
        real dt

        call bfm_module(it,dt)

        end

c**************************************************************
	subroutine bfm_module(it,dt)

c administers bfm ecological model

	implicit none

	integer it
	real dt

	real getpar,areaele

	integer ibfm
	save ibfm
	data ibfm / 0 /
        
	integer ibtanf,ibtend
	save ibtanf,ibtend
     
	if( ibfm .lt. 0 ) return

c------------------------------------------------------------
c initialization
c------------------------------------------------------------

	if( ibfm .eq. 0 ) then
	  ibfm = nint(getpar('ibfm'))
	  if( ibfm .le. 0 ) ibfm = -1
	  if( ibfm .lt. 0 ) return

          ibtanf=getpar('ibtanf')
          ibtend=getpar('ibtend')

	  write(*,*) 'BFM ECOLOGICAL MODEL INCLUDED: ',ibfm
	end if

c------------------------------------------------------------
c is it time ?
c------------------------------------------------------------

        if( it .lt. ibtanf .or. it .gt. ibtend ) return

c------------------------------------------------------------
c call bfm module
c------------------------------------------------------------

	call bfm_shell(it,dt)

c------------------------------------------------------------
c end of routine
c------------------------------------------------------------

	end

c**************************************************************

	subroutine bfm_shell(it,dt)

c computes ecological scalars with BFM  model

	implicit none

	integer it
	real dt

	include 'param.h'
	include 'bfm_common.h'

	integer ndim
	parameter(ndim=nlvdim)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real grav,fcor,dcor,dirn,rowass,roluft
        common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft
        integer nlvdi,nlv
        common /level/ nlvdi,nlv

	integer nen3v(1)
	common /nen3v/nen3v
	real rhov(nlvdim,nkndim)
	common /rhov/rhov

	real hkv(nkndim)
	common /hkv/hkv
	
	real znv(nkndim)
        common /znv/znv

! OBC ARRAY AND VARIABLES

	character*80 bfm1bc(1)
        common /bfm1bc/bfm1bc
	character*80 bfm2bc(1)
        common /bfm2bc/bfm2bc
	character*80 bfm3bc(1)
        common /bfm3bc/bfm3bc

	character*10 what

	integer nbfmv1   	!number of solutes transported var 
        parameter (nbfmv1 = 7)
        integer nbfmv2   	!number of fitozoo transported var
        parameter (nbfmv2 = 9)
        integer nbfmv3   	!number of essudates transported var
        parameter (nbfmv3 = 4)
	

 	real b1cn(nlvdim,nkndim,nbfmv1)
        real b2cn(nlvdim,nkndim,nbfmv2)
	real b2cn_a(nlvdim,nkndim,nbfmv2)
	real b2cn_b(nlvdim,nkndim,nbfmv2)
	real b2cn_c(nlvdim,nkndim,nbfmv2)
	real b2cn_d(nlvdim,nkndim,nbfmv2)

        real b3cn(nlvdim,nkndim,nbfmv3)
	real b3cn_a(nlvdim,nkndim,nbfmv3)
	real b3cn_b(nlvdim,nkndim,nbfmv3)
	real b3cn_c(nlvdim,nkndim,nbfmv3)




 	real bfm1(nb3dim,nbcdim)   !boundary state for solutes
        real bfm2(nb3dim,nbcdim)   !boundary state for fitozoo
	real bfm3(nb3dim,nbcdim)   !boundary state for essudates

 	real b1bound(nbfmv1)       !boundary vector for solutes
 	real b2bound(nbfmv2)       !boundary vector for fitozoo
 	real b3bound(nbfmv3)       !boundary vector for essudates
	save b1bound,b2bound,b3bound

!O2o,N1p,N3n,N4n,O4n,N5s,N6r,B1c

        data b1bound /300.0 , 0.75 , 5.0 , 1.0 , 200.0 , 4.0 , 1.0 / 

!P1c,P2c,P3c,P4c,Z3c,Z4c,Z5c,Z6c

	data b2bound /15.7, 10.0, 10.0, 10.0, 10.0, 1.2, 1.2, 7.2, 7.2/ 

!R1c,R2c,R6c,R7c

	data b3bound /17.0 , 0.1 , 17.0 , 1.0/


	real fc2_a(nbfmv2),fc2_b(nbfmv2),fc2_c(nbfmv2)
     +  ,fc2_d(nbfmv2),fc3_a,fc3_b,fc3_c

	data fc2_a / 0.017,
     +             0.0126,
     +             0.0126,
     +             0.0126,
     +             0.0126,
     +             0.015,
     +             0.015,
     +             0.0167,
     + 	 	   0.0167 /
	data fc2_b / 0.0019,
     +             0.0007862,
     +             0.0007862,
     +             0.0007862,
     +             0.0007862,
     +             0.00167,
     +             0.00167,
     +             0.00185,
     +             0.00185 /
	data fc2_c / 0,
     +             0.05,
     +             0.03,
     +             0.07,
     +             0.02,
     +             0.,
     +             0.,
     +             0.,
     +             0. /
	data fc2_d / 0.,
     +             0.01,
     +             0.,
     +             0.,
     +             0.,
     +             0.,
     +             0.,
     +             0.,
     +             0. /

	

	parameter( fc3_a =0.0126,fc3_b =0.0007862,fc3_c =0.01 )

	real fct
!----------------------------------------------------------------------
!
!	Solutes & nutrients  ONP
!
! 41     O2o            Oxgen                                  mmol O2/m3
! 42     N1p            Phosphate    			       mmol P/m3        
! 43     N3n            Nitrate        			       mmol N/m3
! 44     N4n            Ammonium        		       mmol N/m3
! 45     O4n            NitrogenSink         		       mmol N/m3
! 46     N5s            Silicate         		       mmol Si/m3
! 47     N6r            Reduction Equivalents         	       mmol S--/m3
!
!	Bacteria Fito & Zoo plankton  BFZ
!
! 51     B1c            Pelagic Bacteria         	       mg C/m3
! 52     P1c            Diatoms         		       mg C/m3
! 53     P2c            Flagellates          		       mg C/m3
! 54     P3c            PicoPhytoPlankton      		       mg C/m3
! 55     P4c            Dinoflagellates     		       mg C/m3
! 56     Z3c            Carnivorous mesozooplankton            mg C/m3
! 57     Z4c            Omnivorous mesozooplankton             mg C/m3
! 58     Z5c            Microzooplankton         	       mg C/m3
! 59     Z6c            Heterotrophic nanoflagellates (HNAN)   mg C/m3
!
!       Bacteria Fito & Zoo plankton  NP -- nooutput ---
!         
!	
!	B1n	qnB1c*B1c		
!	P1n	fqnPc(1)*P1c	
!	P2n	fqnPc(2)*P2c	
!	P3n	fqnPc(3)*P3c	
!	P4n 	fqnPc(4)*P4c	
!	Z3n	fqnZc(1)*Z3c	
!	Z4n	fqnZc(2)*Z4c	
!	Z5n	fqnZc(3)*Z5c	
!	Z6n 	fqnZc(4)*Z6c	

!       B1p     qpB1c*B1c               
!       P1p     fqpPc(1)*P1c    
!       P2p     fqpPc(2)*P2c    
!       P3p     fqpPc(3)*P3c    
!       P4p     fqpPc(4)*P4c    
!       Z3p     fqpZc(1)*Z3c    
!       Z4p     fqpZc(2)*Z4c    
!       Z5p     fqpZc(3)*Z5c    
!       Z6p     fqpZc(4)*Z6c    

!       P1l     fqlPc(1)*P1c    
!       P2l     fqlPc(2)*P2c    
!       P3l     fqlPc(3)*P3c    
!       P4l     fqlPc(4)*P4c    

!       P1s     fqsPc(1)*P1c    

!	
!
!	Essudates and methabolic solutes POC
!
! 61     R1c            Labile Organic Carbon (LOC)            mg C/m3
! 62     R2c            CarboHydrates (sugars)                 mg C/m3
! 63     R6c            Particulate Organic Carbon (POC)       mg C/m3
! 64     R7c            Refractory DOC         		       mg C/m3
!
!
!	Essudates and methabolic solutes NP -- nooutput ---
!
!	R1n	?*R1c
!	R1p 	?*R1c
!	R6n	fqnR6c*R6c
!	R6p	fqpR6c*R6c
!
!	R6s	fqsR6c*R6c

!---------------------------------------------------------------------

!---------------OUTPUT--------------------------

	integer itmbfm,idtbfm,ivs1,ivs2,ivs3
	 
c---------------------------------------------------------------
c aux arrays superposed onto other aux arrays
c---------------------------------------------------------------

	real shearf2(nlvdim,nkndim)
	common /saux1/shearf2
	real buoyf2(nlvdim,nkndim)
	common /saux2/buoyf2
	real taub(1)
	common /v1v/taub
	real areaac(1)
	common /v2v/areaac

	integer ilhkv(1)
	common /ilhkv/ilhkv
	real uprv(nlvdim,nkndim)
	common /uprv/uprv
	real vprv(nlvdim,nkndim)
	common /vprv/vprv

	real visv(0:nlvdim,nkndim)
	common /visv/visv
	real difv(0:nlvdim,nkndim)
	common /difv/difv

        real difhv(nlvdim,1)
        common /difhv/difhv

	integer k,l
	integer laux
	integer nlev
	real czdef,taubot
	save czdef

	real ddepth(nkndim)
	common /ddepth/ddepth

	double precision z0s,z0b	!surface and bottom roughness length (m)

	real dtreal
	real areale

	real drr
        common /drr/drr


! function

	real getpar
	integer iround
 
	real difmol,rkpar
	real wsink
	real t

	integer namlst,ll,is1,is2,is3,i
	integer icall
	integer nintp

	save difmol,rkpar
	save icall
	save itmbfm,idtbfm
	save ivs1,ivs2,ivs3
	data icall / 0 /
	data wsink / 0. /

	real tempv(nlvdim,nkndim)
	common /tempv/tempv

	real saltv(nlvdim,nkndim)
        common /saltv/saltv

	! LIGTH DUMMY CALL TO BE CLEANED

	real r1,s1,ss1,tt1,t1

	real rtmp
        common /rtmp/rtmp
	
	real rsal
        common /rsal/rsal

        integer ilp,ivar


	integer ioi,iio
	save ioi,iio	

	real tto(1,nkndim)
	save tto

	real sto(1,nkndim)
        save sto


	real rtr
	save rtr
	real nkk,krf

	integer ok
c------------------------------------------------------
c documentation
c------------------------------------------------------

!------------------------------------------------------
! initialization
!------------------------------------------------------

	if( icall .lt. 0 ) return

	what = 'bfm'

	if( icall .eq. 0 ) then
	  ioi=1 ! DUMMY FOR LIGTH
	  iio=2
	  rkpar=getpar('chpar')
          difmol=getpar('difmol')
	  itmbfm = iround(getpar('itmbfm'))
          idtbfm = iround(getpar('idtbfm'))

!         --------------------------------------------------------
!         Initializes HYBRID HYDRO-BFM arrays
!         --------------------------------------------------------
       
	  call bfm_init(nbfmv1,b1cn,nbfmv2,b2cn
     +,b2cn_a,b2cn_b,b2cn_c,b2cn_d,nbfmv3,b3cn,b3cn_a,b3cn_b,b3cn_c)

!         --------------------------------------------------------
!         Initializes STANDALONE BFM arrays 
!         --------------------------------------------------------

          call feminit_bio(b1cn,nbfmv1,b2cn,nbfmv2,b2cn_a,
     +   b2cn_b,b2cn_c,b2cn_d,b3cn,nbfmv3,b3cn_a,b3cn_b,b3cn_c)
!	  --------------------------------------------------
!	  Initialize HYBRID HYDRO-BFM arrays from external file 
!	  ---------------------------------------------------

!	  call bfm_xy_init

!	  ----------------------------------------------------
!	  Initialize lux data DUMMY VERSION
!	  ----------------------------------------------------
!	  call setlux
!         ----------------------------------------------------
!         Initialize temp data DUMMY VERSION
!         ----------------------------------------------------

!         --------------------------------------------------
!         Initialize boundary conditions for all state variables
!         --------------------------------------------------

          nintp = 2

 	  call bnds_init(what,bfm1bc,nintp,nbfmv1,nb3dim,bfm1,b1bound)
   	  call bnds_init(what,bfm2bc,nintp,nbfmv2,nb3dim,bfm2,b2bound)
   	  call bnds_init(what,bfm3bc,nintp,nbfmv3,nb3dim,bfm3,b3bound)
!         ---------------------------------------------------------
!         INITIALIZES OUPUT
!         ---------------------------------------------------------

	  ivs1 = 0
          call confop(ivs1,itmbfm,idtbfm,1,7,'onp')

          ivs2 = 0
          call confop(ivs2,itmbfm,idtbfm,1,11,'bfz')

	  ivs3 = 0
	  call confop(ivs3,itmbfm,idtbfm,1,4,'poc')

	end if

!------------------------------------------------------
! normal call
!------------------------------------------------------

	icall = icall + 1
	t = it

!------------------------------------------------------
! light ?? (dove viene utilizzato?)
!------------------------------------------------------

!-------------------------------------------------------
! compute OBC for transported vars (HYBRID HYDRO-BFM arrays)
!-------------------------------------------------------
	
	call scal_bnd(what,t,bfm1)
	call scal_bnd(what,t,bfm2)
	call scal_bnd(what,t,bfm3)

!-----------------------------------------------------------
! advection and diffusion of hybrid hydro-bfm var
!-----------------------------------------------------------

!$OMP PARALLEL PRIVATE(is1,is2,is3)
!$OMP DO SCHEDULE(DYNAMIC)

 	do is1=1,nbfmv1
  	    call scal_adv(what,is1,b1cn(1,1,is1),bfm1
     +                          ,rkpar,wsink,difhv,difv,difmol)
 	end do

!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(DYNAMIC)

 	do is2=1,nbfmv2
  	    call scal_adv(what,is2,b2cn(1,1,is2),bfm2
     +                          ,rkpar,wsink,difhv,difv,difmol)
        end do

!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(DYNAMIC)

 	do is3=1,nbfmv3
  	    call scal_adv(what,is3,b3cn(1,1,is3),bfm3
     +                          ,rkpar,wsink,difhv,difv,difmol)
        end do

!$OMP END DO NOWAIT
!$OMP END PARALLEL

!------------------------------------------------------
! ASSIGN DEPTH TO NODE
!------------------------------------------------------

	do k=1,nkn
              ddepth(k) = hkv(k) + znv(k)
	end do


	

!------------------------------------------------------
! call BFM ecological routine
!------------------------------------------------------

	do k=1,nkn
	call do_BFM_ECO(k,b1cn,nbfmv1,b2cn,nbfmv2,b2cn_a,
     +   b2cn_b,b2cn_c,b2cn_d,b3cn,nbfmv3,b3cn_a,b3cn_b,b3cn_c)
	end do



!------------------------------------------------------
! compute output for ecological scalars
!------------------------------------------------------

	if( it .ge. itmbfm )then
      	  call comp_out(ivs1,ivs2,ivs3,itmbfm,idtbfm)
	end if

!-------------------------------------------------
! reset ecological scalars	  
!-------------------------------------------------

!	call reset_out

!------------------------------------------------------
! end of routine
!------------------------------------------------------

	end

c**************************************************************

	subroutine bfm_init(nbfmv1,b1cn,nbfmv2,b2cn
     +,b2cn_a,b2cn_b,b2cn_c,b2cn_d,nbfmv3,b3cn,b3cn_a,b3cn_b,b3cn_c)

c initializes bfm  arrays

	implicit none

	include 'param.h'
	include 'bfm_common.h'
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer nlvdi,nlv
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /level/ nlvdi,nlv
	
	integer nbfmv1
	real b1cn(nlvdim,nkndim,nbfmv1)
	integer nbfmv2
        real b2cn(nlvdim,nkndim,nbfmv2)
        real b2cn_a(nlvdim,nkndim,nbfmv2)
	real b2cn_b(nlvdim,nkndim,nbfmv2)
	real b2cn_c(nlvdim,nkndim,nbfmv2)
	real b2cn_d(nlvdim,nkndim,nbfmv2)
	integer nbfmv3
        real b3cn(nlvdim,nkndim,nbfmv3)
        real b3cn_a(nlvdim,nkndim,nbfmv3)
	real b3cn_b(nlvdim,nkndim,nbfmv3)
	real b3cn_c(nlvdim,nkndim,nbfmv3)



	save  /fO2o/,/fN1p/,/fN3n/,/fN4n/,/fO4n/,/fN5s/,/fN6r/,
     + /fB1c/,/fB1n/,/fB1p/,/fP1c/,/fP1n/,/fP1p/,/fP1l/,/fP1s/,
     + /fP2c/,/fP2n/,/fP2p/,/fP2l/,/fP3c/,/fP3n/,/fP3p/,/fP3l/,
     + /fP4c/,/fP4n/,/fP4p/,/fP4l/,/fZ3c/,/fZ3n/,/fZ3p/,/fZ4c/,
     + /fZ4n/,/fZ4p/,/fZ5c/,/fZ5n/,/fZ5p/,/fZ6c/,/fZ6n/,/fZ6p/,
     + /fR1c/,/fR1n/,/fR1p/,/fR2c/,/fR6c/,/fR6n/,/fR6p/,/fR6s/,
     + /fR7c/

	integer l,k,is

	 do k=1,nkn
          do l=0,nlv
	   do is=1,nbfmv1
	   b1cn(l,k,is)=0.
	   end do
	   do is=1,nbfmv2
           b2cn(l,k,is)=0.
           b2cn_a(l,k,is)=0.
	   b2cn_b(l,k,is)=0.
	   b2cn_c(l,k,is)=0.
	   b2cn_d(l,k,is)=0.
           end do
	   do is=1,nbfmv3
           b3cn(l,k,is)=0.
           b3cn_a(l,k,is)=0.
	   b3cn_b(l,k,is)=0.
	   b3cn_c(l,k,is)=0.
           end do
	  end do
	 end do
	
! init output var

        do k=1,nkn
	   fO2o(k)=0.
           fN1p(k)=0.
           fN3n(k)=0.
           fN4n(k)=0.
           fO4n(k)=0.
           fN5s(k)=0.
           fN6r(k)=0.
           fB1c(k)=0.
           fB1n(k)=0.
           fB1p(k)=0.
           fP1c(k)=0.
           fP1n(k)=0.
           fP1p(k)=0.
           fP1l(k)=0.
           fP1s(k)=0.
           fP2c(k)=0.
           fP2n(k)=0.
           fP2p(k)=0.
           fP2l(k)=0.
           fP3c(k)=0.
           fP3n(k)=0.
           fP3p(k)=0.
           fP3l(k)=0.
           fP4c(k)=0.
           fP4n(k)=0.
           fP4p(k)=0.
           fP4l(k)=0.
           fZ3c(k)=0.
           fZ3n(k)=0.
           fZ3p(k)=0.
           fZ4c(k)=0.
           fZ4n(k)=0.
           fZ4p(k)=0.
           fZ5c(k)=0.
           fZ5n(k)=0.
           fZ5p(k)=0.
           fZ6c(k)=0.
           fZ6n(k)=0.
           fZ6p(k)=0.
           fR1c(k)=0.
           fR1n(k)=0.
           fR1p(k)=0.
           fR2c(k)=0.
           fR6c(k)=0.
           fR6n(k)=0.
           fR6p(k)=0.
           fR6s(k)=0.
           fR7c(k)=0.
        end do
  
	end
  
c**************************************************************

	subroutine comp_out(ivs1,ivs2,ivs3,itmbfm,idtbfm)

c outputs variables

	implicit none

	integer ivs1,ivs2,ivs3
	integer itmbfm,idtbfm

	include 'param.h'
	include 'bfm_common.h'

	call confil(ivs1,itmbfm,idtbfm,41,1,fO2o)
!	call confil(ivs1,itmbfm,idtbfm,42,1,fN1p)
!        call confil(ivs1,itmbfm,idtbfm,43,1,fN3n)
!        call confil(ivs1,itmbfm,idtbfm,44,1,fN4n)
!        call confil(ivs1,itmbfm,idtbfm,45,1,fO4n)
!        call confil(ivs1,itmbfm,idtbfm,46,1,fN5s)
!        call confil(ivs1,itmbfm,idtbfm,47,1,fN6r)

!        call confil(ivs2,itmbfm,idtbfm,51,1,fB1c)
!        call confil(ivs2,itmbfm,idtbfm,52,1,fP1c)
!        call confil(ivs2,itmbfm,idtbfm,53,1,fP2c)
!        call confil(ivs2,itmbfm,idtbfm,54,1,fP3c)
!        call confil(ivs2,itmbfm,idtbfm,55,1,fP4c)
!        call confil(ivs2,itmbfm,idtbfm,56,1,fZ3c)
!        call confil(ivs2,itmbfm,idtbfm,57,1,fZ4c)
!        call confil(ivs2,itmbfm,idtbfm,58,1,fZ5c)
!        call confil(ivs2,itmbfm,idtbfm,59,1,fZ6c)
!	call confil(ivs2,itmbfm,idtbfm,72,1,fP1l)
!	call confil(ivs2,itmbfm,idtbfm,75,1,fP4l)
	


!        call confil(ivs3,itmbfm,idtbfm,61,1,fR1c)
!        call confil(ivs3,itmbfm,idtbfm,62,1,fR2c)
!        call confil(ivs3,itmbfm,idtbfm,63,1,fR6c)
!        call confil(ivs3,itmbfm,idtbfm,64,1,fR7c)
  
	end
  

