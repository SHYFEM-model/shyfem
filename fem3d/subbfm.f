
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
c 10.03.2008	aac	bfm ecological module from scratch
c 29.04.2008	ggu	bfm model integrated in main branch
c 30.04.2008	ggu	double to real (BUG)
c 23.03.2010	ggu	changed v6.1.1
c 17.02.2011	ggu	changed VERS_6_1_18
c 18.02.2011	ggu	changed VERS_6_1_19
c 31.05.2011	ggu	clean from useless common blocks
c 07.06.2011	ggu	changed VERS_6_1_25
c 15.01.2012	aac	advection for all BFM var introduced
c 17.02.2012	aac&ggu	restart for bfm
c 26.03.2012	ggu	bfm1-3 had wrong second dimension
c 01.06.2012	ggu	changed VERS_6_1_53
c 22.10.2012	ggu	saved some variables
c 17.06.2013	ggu	bug fix: wsinkv was not present in call
c 18.06.2014	ggu	changed VERS_6_1_77
c 18.07.2014	ggu	changed VERS_7_0_1
c 21.10.2014	ggu	converted to new boundary treatment
c 05.11.2014	ggu	changed VERS_7_0_5
c 23.12.2014	ggu	changed VERS_7_0_11
c 19.01.2015	ggu	changed VERS_7_1_3
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 30.07.2015	ggu	changed VERS_7_1_83
c 03.04.2018	ggu	changed VERS_7_5_43
c 16.02.2019	ggu	changed VERS_7_5_60
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
c ddepth		a che cosa serve ???? -> needed in standalone
c drr			-> needed in standalone
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

        subroutine ecological_module

c general interface to ecological module

        implicit none

        integer it
        real dt
	double precision dtime

	call get_timestep(dt)
	call get_act_dtime(dtime)
	it = dtime

        call bfm_module(it,dt)

        end

c**************************************************************

	subroutine bfm_module(it,dt)

c administers bfm ecological model

	implicit none

	integer it
	real dt

	real getpar,areaele

	integer, save :: ibfm = 0
        
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

	use mod_sinking
	use mod_depth
	use mod_ts
	use mod_diff_visc_fric
	use mod_hydro_print
	use mod_hydro
	use levels
	use basin

	implicit none

	integer it
	real dt

	include 'param.h'
	include 'bfm_common.h'

	integer ndim
	parameter(ndim=nlvdim)

	include 'pkonst.h'


	

! OBC ARRAY AND VARIABLES

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
	save b1cn,b2cn,b2cn_a,b2cn_b,b2cn_c,b2cn_d

        real b3cn(nlvdim,nkndim,nbfmv3)
	real b3cn_a(nlvdim,nkndim,nbfmv3)
	real b3cn_b(nlvdim,nkndim,nbfmv3)
	real b3cn_c(nlvdim,nkndim,nbfmv3)
	save b3cn,b3cn_a,b3cn_b,b3cn_c

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
	 



	real load(nlvdim,nkndim)
	common /load/load

	integer k,l
	integer laux
	integer nlev
	real czdef,taubot
	save czdef

	double precision dtime0
	integer idbfm1(nbcdim)
	integer idbfm2(nbcdim)
	integer idbfm3(nbcdim)
	save idbfm1,idbfm2,idbfm3

	real ddepth(nkndim)
	common /ddepth/ddepth		!-> needed in standalone

	double precision z0s,z0b	!surface and bottom roughness length (m)

	real dtreal
	real areale
	real rload

	real drr
        common /drr/drr			!-> needed in standalone

! function
	logical has_restart

	real getpar
	integer iround
 
	real difmol,rkpar
	real wsink
	real t

	integer namlst,ll,is1,is2,is3,is4,is5,is6,is7,is8,is9,i
	integer icall
	integer nintp

	save difmol,rkpar
	save icall
	save itmbfm,idtbfm
	save ivs1,ivs2,ivs3
	data icall / 0 /
	data wsink / 0. /

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

	  if( .not. has_restart(6) ) then       
	    call bfm_init(nbfmv1,b1cn,nbfmv2,b2cn
     +			,b2cn_a,b2cn_b,b2cn_c,b2cn_d,nbfmv3
     +			,b3cn,b3cn_a,b3cn_b,b3cn_c)
	  endif

!         --------------------------------------------------------
!         Initializes STANDALONE BFM arrays 
!         --------------------------------------------------------

          call feminit_bio(b1cn,nbfmv1,b2cn,nbfmv2,b2cn_a
     +			,b2cn_b,b2cn_c,b2cn_d,b3cn,nbfmv3
     +			,b3cn_a,b3cn_b,b3cn_c)

!         --------------------------------------------------
!         Initialize boundary conditions for all state variables
!         --------------------------------------------------

 	  !call bnds_init0('bfm1',bfm1bc,nintp,nbfmv1,nb3dim,bfm1,b1bound)
   	  !call bnds_init0('bfm2',bfm2bc,nintp,nbfmv2,nb3dim,bfm2,b2bound)
   	  !call bnds_init0('bfm3',bfm3bc,nintp,nbfmv3,nb3dim,bfm3,b3bound)

	  call get_first_dtime(dtime0)
          nintp = 2
          call bnds_init_new('bfm1',dtime0,nintp,nbfmv1,nkn,nlv
     +                          ,b1bound,idbfm1)
          call bnds_init_new('bfm2',dtime0,nintp,nbfmv2,nkn,nlv
     +                          ,b2bound,idbfm2)
          call bnds_init_new('bfm3',dtime0,nintp,nbfmv3,nkn,nlv
     +                          ,b3bound,idbfm3)

!         ---------------------------------------------------------
!         INITIALIZES load
!         ---------------------------------------------------------

	  load = 0.

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
	rload = 0.

!------------------------------------------------------
! light ?? (dove viene utilizzato?)
!------------------------------------------------------

!-------------------------------------------------------
! compute OBC for transported vars (HYBRID HYDRO-BFM arrays)
!-------------------------------------------------------
	
	!call scal_bnd('bfm1',t,bfm1)
	!call scal_bnd('bfm2',t,bfm2)
	!call scal_bnd('bfm3',t,bfm3)

!-----------------------------------------------------------
! advection and diffusion of hybrid hydro-bfm var
!-----------------------------------------------------------

!$OMP PARALLEL PRIVATE(is1,is2,is3,is4,is5,is6,is7,is8,is9,fct)
!$OMP DO SCHEDULE(DYNAMIC)

 	do is1=1,nbfmv1
  	    call scal_adv('bfm_1',is1,b1cn(1,1,is1),idbfm1
     +                          ,rkpar,wsink,difhv,difv,difmol)
 	end do

!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(DYNAMIC)

 	do is2=1,nbfmv2
  	    call scal_adv('bfm_2',is2,b2cn(1,1,is2),idbfm2
     +                          ,rkpar,wsink,difhv,difv,difmol)
        end do

!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(DYNAMIC)

 	do is3=1,nbfmv3
  	    call scal_adv('bfm_3',is3,b3cn(1,1,is3),idbfm3
     +                          ,rkpar,wsink,difhv,difv,difmol)
        end do

!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(DYNAMIC)

        do is4=1,nbfmv2
            fct=fc2_a(is4)
            call scal_adv_fact('bfm_4',is4,fct,b2cn_a(1,1,is4),idbfm2
     +                          ,rkpar,wsink,wsinkv,rload,load
     +				,difhv,difv,difmol)
        end do

!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(DYNAMIC)

        do is5=1,nbfmv2
            fct=fc2_b(is5)
            call scal_adv_fact('bfm_5',is5,fct,b2cn_b(1,1,is5),idbfm2
     +                          ,rkpar,wsink,wsinkv,rload,load
     +				,difhv,difv,difmol)
        end do

!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(DYNAMIC)

        do is6=2,5
            fct=fc2_c(is6)
            call scal_adv_fact('bfm_6',is6,fct,b2cn_a(1,1,is6),idbfm2
     +                          ,rkpar,wsink,wsinkv,rload,load
     +				,difhv,difv,difmol)
        end do

!$OMP END DO NOWAIT

!        do is=1,nbfmv2
            is7 = 2
            fct=fc2_d(is7)
            call scal_adv_fact('bfm_7',is7,fct,b2cn_d(1,1,is7),idbfm2
     +                          ,rkpar,wsink,wsinkv,rload,load
     +				,difhv,difv,difmol)
!        end do

!$OMP DO SCHEDULE(DYNAMIC)

         do is8=1,3,2
            fct=fc3_a
            call scal_adv_fact('bfm_8',is8,fct,b3cn_a(1,1,is8),idbfm3
     +                          ,rkpar,wsink,wsinkv,rload,load
     +				,difhv,difv,difmol)

            fct=fc3_b
            call scal_adv_fact('bfm_8',is8,fct,b3cn_b(1,1,is8),idbfm3
     +                          ,rkpar,wsink,wsinkv,rload,load
     +				,difhv,difv,difmol)
         end do

!$OMP END DO NOWAIT
!$OMP END PARALLEL

         is9 = 3
         fct=fc3_c
         call scal_adv_fact('bfm_9',is9,fct,b3cn_c(1,1,is9),idbfm3
     +                          ,rkpar,wsink,wsinkv,rload,load
     +				,difhv,difv,difmol)

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
 	  call do_BFM_ECO(k,b1cn,nbfmv1,b2cn,nbfmv2,b2cn_a
     +      ,b2cn_b,b2cn_c,b2cn_d,b3cn,nbfmv3,b3cn_a,b3cn_b,b3cn_c)
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
     +			,b2cn_a,b2cn_b,b2cn_c,b2cn_d,nbfmv3
     +			,b3cn,b3cn_a,b3cn_b,b3cn_c)

c initializes bfm  arrays

	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'
	include 'bfm_common.h'
	
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
 	call confil(ivs1,itmbfm,idtbfm,42,1,fN1p)
        call confil(ivs1,itmbfm,idtbfm,43,1,fN3n)
        call confil(ivs1,itmbfm,idtbfm,44,1,fN4n)
        call confil(ivs1,itmbfm,idtbfm,45,1,fO4n)
        call confil(ivs1,itmbfm,idtbfm,46,1,fN5s)
        call confil(ivs1,itmbfm,idtbfm,47,1,fN6r)

        call confil(ivs2,itmbfm,idtbfm,51,1,fB1c)
        call confil(ivs2,itmbfm,idtbfm,52,1,fP1c)
        call confil(ivs2,itmbfm,idtbfm,53,1,fP2c)
        call confil(ivs2,itmbfm,idtbfm,54,1,fP3c)
        call confil(ivs2,itmbfm,idtbfm,55,1,fP4c)
        call confil(ivs2,itmbfm,idtbfm,56,1,fZ3c)
        call confil(ivs2,itmbfm,idtbfm,57,1,fZ4c)
        call confil(ivs2,itmbfm,idtbfm,58,1,fZ5c)
        call confil(ivs2,itmbfm,idtbfm,59,1,fZ6c)
 	call confil(ivs2,itmbfm,idtbfm,72,1,fP1l)
 	call confil(ivs2,itmbfm,idtbfm,75,1,fP4l)
 
        call confil(ivs3,itmbfm,idtbfm,61,1,fR1c)
        call confil(ivs3,itmbfm,idtbfm,62,1,fR2c)
        call confil(ivs3,itmbfm,idtbfm,63,1,fR6c)
        call confil(ivs3,itmbfm,idtbfm,64,1,fR7c)
  
	end
  
c**************************************************************

	subroutine write_restart_eco(iunit)

	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'
        include 'bfm_common.h'

	integer iunit
	integer k,nstate

	nstate = 48

	write(6,*) 'write_restart_bfm: ',fO2o(1)

        write(iunit) nstate,nkn

        write(iunit)(fO2o(k),k=1,nkn)
        write(iunit)(fN1p(k),k=1,nkn)
        write(iunit)(fN3n(k),k=1,nkn)
        write(iunit)(fN4n(k),k=1,nkn)
        write(iunit)(fO4n(k),k=1,nkn)
        write(iunit)(fN5s(k),k=1,nkn)
        write(iunit)(fN6r(k),k=1,nkn)
        write(iunit)(fB1c(k),k=1,nkn)
        write(iunit)(fB1n(k),k=1,nkn)
        write(iunit)(fB1p(k),k=1,nkn)
        write(iunit)(fP1c(k),k=1,nkn)
        write(iunit)(fP1n(k),k=1,nkn)
        write(iunit)(fP1p(k),k=1,nkn)
        write(iunit)(fP1l(k),k=1,nkn)
        write(iunit)(fP1s(k),k=1,nkn)
        write(iunit)(fP2c(k),k=1,nkn)
        write(iunit)(fP2n(k),k=1,nkn)
        write(iunit)(fP2p(k),k=1,nkn)
        write(iunit)(fP2l(k),k=1,nkn)
        write(iunit)(fP3c(k),k=1,nkn)
        write(iunit)(fP3n(k),k=1,nkn)
        write(iunit)(fP3p(k),k=1,nkn)
        write(iunit)(fP3l(k),k=1,nkn)
        write(iunit)(fP4c(k),k=1,nkn)
        write(iunit)(fP4n(k),k=1,nkn)
        write(iunit)(fP4p(k),k=1,nkn)
        write(iunit)(fP4l(k),k=1,nkn)
        write(iunit)(fZ3c(k),k=1,nkn)
        write(iunit)(fZ3n(k),k=1,nkn)
        write(iunit)(fZ3p(k),k=1,nkn)
        write(iunit)(fZ4c(k),k=1,nkn)
        write(iunit)(fZ4n(k),k=1,nkn)
        write(iunit)(fZ4p(k),k=1,nkn)
        write(iunit)(fZ5c(k),k=1,nkn)
        write(iunit)(fZ5n(k),k=1,nkn)
        write(iunit)(fZ5p(k),k=1,nkn)
        write(iunit)(fZ6c(k),k=1,nkn)
        write(iunit)(fZ6n(k),k=1,nkn)
        write(iunit)(fZ6p(k),k=1,nkn)
        write(iunit)(fR1c(k),k=1,nkn)
        write(iunit)(fR1n(k),k=1,nkn)
        write(iunit)(fR1p(k),k=1,nkn)
        write(iunit)(fR2c(k),k=1,nkn)
        write(iunit)(fR6c(k),k=1,nkn)
        write(iunit)(fR6n(k),k=1,nkn)
        write(iunit)(fR6p(k),k=1,nkn)
        write(iunit)(fR6s(k),k=1,nkn)
        write(iunit)(fR7c(k),k=1,nkn)

	end

c**************************************************************

	subroutine skip_restart_eco(iunit)

	implicit none

	integer iunit

	integer nstate,i
	integer nstate_aux,nkn_aux

	nstate = 48

        read(iunit) nstate_aux,nkn_aux
	if( nstate_aux .ne. nstate ) then
	  write(6,*) 'nstate (file): ',nstate_aux
	  write(6,*) 'nstate (sim) : ',nstate
	  stop 'error stop read_restart_eco: incompatible parameters'
	end if

	do i=1,nstate
	  read(iunit)
	end do

	end

c**************************************************************

	subroutine read_restart_eco(iunit)

	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'
        include 'bfm_common.h'

	integer iunit

	integer k,nstate
	integer nstate_aux,nkn_aux

	nstate = 48

	write(6,*) 'read_restart_bfm: ',fO2o(1)

        read(iunit) nstate_aux,nkn_aux
	if( nstate_aux .ne. nstate .or. nkn_aux .ne. nkn ) then
	  write(6,*) 'nstate,nkn (file): ',nstate_aux,nkn_aux
	  write(6,*) 'nstate,nkn (sim) : ',nstate,nkn
	  stop 'error stop read_restart_eco: incompatible parameters'
	end if

        read(iunit)(fO2o(k),k=1,nkn)
        read(iunit)(fN1p(k),k=1,nkn)
        read(iunit)(fN3n(k),k=1,nkn)
        read(iunit)(fN4n(k),k=1,nkn)
        read(iunit)(fO4n(k),k=1,nkn)
        read(iunit)(fN5s(k),k=1,nkn)
        read(iunit)(fN6r(k),k=1,nkn)
        read(iunit)(fB1c(k),k=1,nkn)
        read(iunit)(fB1n(k),k=1,nkn)
        read(iunit)(fB1p(k),k=1,nkn)
        read(iunit)(fP1c(k),k=1,nkn)
        read(iunit)(fP1n(k),k=1,nkn)
        read(iunit)(fP1p(k),k=1,nkn)
        read(iunit)(fP1l(k),k=1,nkn)
        read(iunit)(fP1s(k),k=1,nkn)
        read(iunit)(fP2c(k),k=1,nkn)
        read(iunit)(fP2n(k),k=1,nkn)
        read(iunit)(fP2p(k),k=1,nkn)
        read(iunit)(fP2l(k),k=1,nkn)
        read(iunit)(fP3c(k),k=1,nkn)
        read(iunit)(fP3n(k),k=1,nkn)
        read(iunit)(fP3p(k),k=1,nkn)
        read(iunit)(fP3l(k),k=1,nkn)
        read(iunit)(fP4c(k),k=1,nkn)
        read(iunit)(fP4n(k),k=1,nkn)
        read(iunit)(fP4p(k),k=1,nkn)
        read(iunit)(fP4l(k),k=1,nkn)
        read(iunit)(fZ3c(k),k=1,nkn)
        read(iunit)(fZ3n(k),k=1,nkn)
        read(iunit)(fZ3p(k),k=1,nkn)
        read(iunit)(fZ4c(k),k=1,nkn)
        read(iunit)(fZ4n(k),k=1,nkn)
        read(iunit)(fZ4p(k),k=1,nkn)
        read(iunit)(fZ5c(k),k=1,nkn)
        read(iunit)(fZ5n(k),k=1,nkn)
        read(iunit)(fZ5p(k),k=1,nkn)
        read(iunit)(fZ6c(k),k=1,nkn)
        read(iunit)(fZ6n(k),k=1,nkn)
        read(iunit)(fZ6p(k),k=1,nkn)
        read(iunit)(fR1c(k),k=1,nkn)
        read(iunit)(fR1n(k),k=1,nkn)
        read(iunit)(fR1p(k),k=1,nkn)
        read(iunit)(fR2c(k),k=1,nkn)
        read(iunit)(fR6c(k),k=1,nkn)
        read(iunit)(fR6n(k),k=1,nkn)
        read(iunit)(fR6p(k),k=1,nkn)
        read(iunit)(fR6s(k),k=1,nkn)
        read(iunit)(fR7c(k),k=1,nkn)

	write(6,*) 'read_restart_bfm: ',fO2o(1)

	end

c**************************************************************

