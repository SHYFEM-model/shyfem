#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: Phyto
!
! DESCRIPTION
!   This process describes the dynamics of all phytoplankton
!    groups in the ERSEM model. The differences in behaviour
!    are expressed by differences in parameter-values only.
!    
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine PhytoDynamics(phyto, ppphytoc, ppphyton, ppphytop, ppphytos, &
    ppphytol)
!
! !USES:

  ! For the following Pelagic-states fluxes are defined: R1c, R6c, O2o, R2c, &
  ! N3n, N4n, N1p, R1n, R6n, R1p, R6p, N5s
  ! The following global scalar vars are used: SUNQ, ThereIsLight
  ! The following Pelagic 1-d global boxvars are modified : flP1R6s
  ! The following Pelagic 1-d global boxvars  are used: ETW, EIR, xEPS, Depth
  ! The following Pelagic 2-d global boxvars are modified : eiPI, sediPI
  ! The following Pelagic 2-d global boxvars got a value: sunPI
  ! The following Pelagic 2-d global boxvars  are used: qpPc, qnPc, qsPc, qlPc
  ! The following groupmember vars  are used: iiP1, iiP4
  ! The following 0-d global box parametes are used: p_small, &
  ! ChlLightFlag, LightForcingFlag
  ! The following 1-d global parameter vars are used: p_qchlc
  ! The following global constants are used: RLEN
  ! The following constants are used: SEC_PER_DAY, E2W, HOURS_PER_DAY

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
#ifdef NOPOINTERS
  use mem,  ONLY: D3STATE
#else
  use mem, ONLY: D3STATE, R1c, R6c, O2o, R2c, N3n, N4n, N1p, R1n, R6n, R1p, R6p, &
    N5s
#endif
  use mem, ONLY: ppR1c, ppR6c, ppO2o, ppR2c, ppN3n, ppN4n, ppN1p, ppR1n, &
    ppR6n, ppR1p, ppR6p, ppN5s, SUNQ, ThereIsLight, flP1R6s, ETW, EIR, xEPS, &
    Depth, eiPI, sediPI, sunPI, qpPc, qnPc, qsPc, qlPc, iiP1, iiP4, NO_BOXES, &
    iiBen, iiPel, flux_vector
  use constants,  ONLY: SEC_PER_DAY, E2W, HOURS_PER_DAY
  use mem_Param,  ONLY: p_small, ChlLightFlag, LightForcingFlag, p_qchlc
  use mem_Phyto


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following vector functions are used:eTq_vector, MM_vector, insw_vector
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,   ONLY: eTq_vector, MM_vector, insw_vector


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE

! !INPUT:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer,intent(IN)  :: phyto
  integer,intent(IN) :: ppphytoc
  integer,intent(IN) :: ppphyton
  integer,intent(IN) :: ppphytop
  integer,intent(IN) :: ppphytos
  integer,intent(IN) :: ppphytol

!  
!
! !AUTHORS
!   ERSEM group + J.G. Baretta-Bekker + W.Ebenhoeh
!     P. Ruardij (NIOZ)
!
!
!
! !REVISION_HISTORY
!   !
!
! COPYING
!   
!   Copyright (C) 2006 P. Ruardij, the mfstep group, the ERSEM team 
!   (rua@nioz.nl, vichi@bo.ingv.it)
!
!   This program is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation;
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTEABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!EOP
!-------------------------------------------------------------------------!
!BOC
!
!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Set up Local Variable for copy of state var. object
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN),dimension(NO_BOXES) :: phytoc
  real(RLEN),dimension(NO_BOXES) :: phyton
  real(RLEN),dimension(NO_BOXES) :: phytop
  real(RLEN),dimension(NO_BOXES) :: phytos
  real(RLEN),dimension(NO_BOXES) :: phytol
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer,dimension(NO_BOXES)  :: i
  real(RLEN),dimension(NO_BOXES)  :: r
  real(RLEN),dimension(NO_BOXES)  :: et
  real(RLEN),dimension(NO_BOXES)  :: sum
  real(RLEN),dimension(NO_BOXES)  :: sadap
  real(RLEN),dimension(NO_BOXES)  :: sea
  real(RLEN),dimension(NO_BOXES)  :: sdo
  real(RLEN),dimension(NO_BOXES)  :: rugc
  real(RLEN),dimension(NO_BOXES)  :: sra
  real(RLEN),dimension(NO_BOXES)  :: srs
  real(RLEN),dimension(NO_BOXES)  :: srt
  real(RLEN),dimension(NO_BOXES)  :: slc
  real(RLEN),dimension(NO_BOXES)  :: run
  real(RLEN),dimension(NO_BOXES)  :: pe_R6
  real(RLEN),dimension(NO_BOXES)  :: rupp
  real(RLEN),dimension(NO_BOXES)  :: rump
  real(RLEN),dimension(NO_BOXES)  :: misp
  real(RLEN),dimension(NO_BOXES)  :: rupn
  real(RLEN),dimension(NO_BOXES)  :: rumn3
  real(RLEN),dimension(NO_BOXES)  :: rumn4
  real(RLEN),dimension(NO_BOXES)  :: rumn
  real(RLEN),dimension(NO_BOXES)  :: netgrowth
  real(RLEN),dimension(NO_BOXES)  :: misn
  real(RLEN),dimension(NO_BOXES)  :: cqun3
  real(RLEN),dimension(NO_BOXES)  :: rums
  real(RLEN),dimension(NO_BOXES)  :: rups
  real(RLEN),dimension(NO_BOXES)  :: miss
  real(RLEN),dimension(NO_BOXES)  :: iN
  real(RLEN),dimension(NO_BOXES)  :: iN1p
  real(RLEN),dimension(NO_BOXES)  :: iNIn
  real(RLEN),dimension(NO_BOXES)  :: iN5s
  real(RLEN),dimension(NO_BOXES)  :: rrc
  real(RLEN),dimension(NO_BOXES)  :: rr1c
  real(RLEN),dimension(NO_BOXES)  :: rr1n
  real(RLEN),dimension(NO_BOXES)  :: rr1p
  real(RLEN),dimension(NO_BOXES)  :: rr6c
  real(RLEN),dimension(NO_BOXES)  :: rr6n
  real(RLEN),dimension(NO_BOXES)  :: rr6p
  real(RLEN),dimension(NO_BOXES)  :: rr6s
  real(RLEN),dimension(NO_BOXES)  :: runn
  real(RLEN),dimension(NO_BOXES)  :: runn3
  real(RLEN),dimension(NO_BOXES)  :: runn4
  real(RLEN),dimension(NO_BOXES)  :: runp
  real(RLEN),dimension(NO_BOXES)  :: runs
  real(RLEN),dimension(NO_BOXES)  :: Irr
  real(RLEN),dimension(NO_BOXES)  :: rho_Chl
  real(RLEN),dimension(NO_BOXES)  :: rate_Chl
  real(RLEN),dimension(NO_BOXES)  :: Photo_max
  real(RLEN),dimension(NO_BOXES)  :: flPIR2c
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Copy  state var. object in local var
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  phytoc = D3STATE(ppphytoc,:)
  phyton = D3STATE(ppphyton,:)
  phytop = D3STATE(ppphytop,:)
  phytol = D3STATE(ppphytol,:)
  if ( ppphytos > 0 )  phytos = D3STATE(ppphytos,:)



  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient limitation (intracellular) N, P
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  iN1p = min( 1.0D+00, max( p_small, ( qpPc(phyto, &
    :)- p_qplc(phyto))/( p_qpRc(phyto)- p_qplc(phyto))))
  iNIn = min( 1.0D+00, max( p_small, ( qnPc(phyto, &
    :)- p_qnlc(phyto))/( p_qnRc(phyto)- p_qnlc(phyto))))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient limitation due to intra- extracellular silicate.
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  select case ( phyto)
    case ( iiP1 )
      iN5s = min( 1.0D+00, max( p_small, ( qsPc(phyto, &
        :)- p_qslc(phyto))/(( p_qsRc(phyto)- p_qslc(phyto)))))



    case default
      iN5s  =   1.0D+00


  end select


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Phytoplankton growth is limited by nitrogen and phosphorus
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  select case ( p_limnut(phyto))
    case ( 0 )
      iN  =   (iN1p* iNIn)**(0.5D+00)  ! geometric mean


    case ( 1 )
      iN  =   min(  iN1p,  iNIn)  ! Liebig rule


    case ( 2 )
      iN  =   2.0D+00/( 1.0D+00/ iN1p+ 1.0D+00/ iNIn)  ! combined


  end select



  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Temperature response of Phytoplankton
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  et  =   eTq_vector(  ETW(:),  p_q10(phyto))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Photosynthesis (Irradiance EIR is in uE m-2 s-1, Irr is mid-layer EIR)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  if ( ChlLightFlag== 2) then
    !Irr = EIR*exp(-xEPS* 0.5 * Depth)*SEC_PER_DAY;
    Irr = max( p_small, EIR(:)/ xEPS(:)/ Depth(:)*( 1.0D+00- exp( &
      - xEPS(:)* Depth(:)))* SEC_PER_DAY)

    eiPI(phyto,:) = ( 1.0D+00- exp( - qlPc(phyto, :)* p_alpha_chl(phyto)/ &
      p_sum(phyto)* Irr))

    ! ( See CalcLightDistribution)
    ! if one use average light per day calculate sum according:
  end if

  select case ( LightForcingFlag)
    case ( 1 )
      sum  =   p_sum(phyto)* et* eiPI(phyto,:)* iN5s



    case ( 2 )
      sum  =   p_sum(phyto)* et* eiPI(phyto,:)* iN5s*( SUNQ/ HOURS_PER_DAY)



    case ( 3 )
      sum  =   p_sum(phyto)* et* eiPI(phyto,:)* iN5s* ThereIsLight


  end select

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Lysis and excretion
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  sdo  =  ( p_thdo(phyto)/( iN+ p_thdo(phyto)))* p_sdmo(phyto)  ! nutr. -stress lysis

  select case ( phyto)

    case ( iiP4 )  ! extra lysis for P4 only
      sdo  =   sdo+ p_seo(phyto)* MM_vector(  phytoc,  100.0D+00)



  end select


  sea  =   sum* p_pu_ea(phyto)  ! activity excretion

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Apportioning over R1 and R6:
  ! Cell lysis generates both DOM and POM.
  ! The nutr.-depleted cell has a nutrient-carbon ratio equal to p_q?lc.
  ! Assuming that this structural part is not easily degradable,
  ! at least a fraction equal to the minimum quota is released as POM.
  ! Therefore, nutrients (and C) in the structural part go to R6.
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  pe_R6 = min( p_qplc(phyto)/( qpPc(phyto, :)+ p_small), p_qnlc(phyto)/( &
    qnPc(phyto, :)+ p_small))
  pe_R6  =   min(  1.0D+00,  pe_R6)
  rr6c  =   pe_R6* sdo* phytoc
  rr1c  =  ( 1.0D+00- pe_R6)* sdo* phytoc

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Respiration rate
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  sra  =   p_pu_ra(phyto)*( sum- sea)  ! activity
  srs  =   et* p_srs(phyto)  ! rest
  srt  =   sra+ srs  ! total
  rrc  =   srt* phytoc  ! total actual respiration

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Production, productivity and C flows
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rugc  =   sum* phytoc  ! gross production
  slc  =   sea+ srt+ sdo  ! specific loss terms
  call flux_vector( iiPel, ppphytoc,ppphytoc,-(- rugc) )  ! source/sink.c
  call flux_vector( iiPel, ppphytoc,ppR1c, rr1c )  ! source/sink.c
  call flux_vector( iiPel, ppphytoc,ppR6c, rr6c )  ! source/sink.c

  !  These fluxes are first added before defining as flux!
  flPIR2c  =   sea* phytoc

  call flux_vector( iiPel, ppphytoc,ppphytoc,-( rrc) )  ! source/sink.c
  call flux_vector( iiPel, ppO2o,ppO2o,-( rrc/ 12.0D+00) )  ! source/sink.o
  call flux_vector( iiPel, ppO2o,ppO2o, rugc/ 12.0D+00 )  ! source/sink.o


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Potential-Net prim prod. (mgC /m3/d)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  sadap  =   max(  srs,  sum)
  run  =   max(  0.0D+00, ( sum- slc)* phytoc)  ! net production

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient Uptake: calculate maximal uptake of N, P
  ! Check if C-fixation is # larger to make of all C new biomass
  ! Assumed is that Si-depletion directly the growth rate in contradiction
  ! to N and P.
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  cqun3  =   p_lN4(phyto)/( p_lN4(phyto)+ N4n(:))
  rumn3  =   p_qun(phyto)* N3n(:)* phytoc* cqun3  ! max pot. uptake of N3
  rumn4  =   p_qun(phyto)* N4n(:)* phytoc  ! max pot. uptake of N4
  rumn  =   rumn3+ rumn4  ! max pot. uptake of NI

  rump  =   p_qup(phyto)* N1p(:)* phytoc  ! max pot. uptake

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Check if which fraction C-fixation can be used for new biomass
  ! by checking the potential nutrient avilability
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  netgrowth = min( run, ( rumn+ max( 0.0D+00, 0.05D+00* &
    rugc*( qnPc(phyto, :)- p_qnlc(phyto))))/ p_qnlc(phyto))
  netgrowth = min( netgrowth, ( rump+ max( 0.0D+00, &
    0.05D+00* rugc*( qpPc(phyto, :)- p_qplc(phyto))))/ p_qplc(phyto))


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Excrete C which can not be used for growth as carbo-hydrates:
  ! Correct net C-uptake
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  netgrowth  =   max(  netgrowth,  0.0D+00)
  flPIR2c  =   flPIR2c+ run- netgrowth
  call flux_vector( iiPel, ppphytoc,ppR2c, flPIR2c )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Apparent Net prim prod. (mgC /m3/d)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  run  =   netgrowth
  sunPI(phyto,:)  =   run/( p_small+ phytoc)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient dynamics: NITROGEN
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


  misn  =   sadap*( p_xqn(phyto)* p_qnRc(phyto)* phytoc- phyton)  ! Intracellular missing amount of N
  rupn  =   p_xqn(phyto)* p_qnRc(phyto)* run-( srs+ sdo)* phyton  ! N uptake based on net assimilat. C
  runn  =   min(  rumn,  rupn+ misn)  ! actual uptake of NI

  r  =   insw_vector(  runn)
  runn3  =   r* runn* rumn3/( p_small+ rumn)  ! actual uptake of Nn
  runn4  =   r* runn* rumn4/( p_small+ rumn)  ! actual uptake of Nn
  call flux_vector( iiPel, ppN3n,ppphyton, runn3 )  ! source/sink.n
  call flux_vector( iiPel, ppN4n,ppphyton, runn4 )  ! source/sink.n
  call flux_vector(iiPel, ppphyton,ppN4n,- runn*( 1.0D+00- r))  ! source/sink.n

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nuttrient dynamics: PHOSPHORUS
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  misp  =   sadap*( p_xqp(phyto)* p_qpRc(phyto)* phytoc- phytop)  ! intracellular missing amount of P
  rupp  =   p_xqp(phyto)* run* p_qpRc(phyto)-( sdo+ srs)* phytop  ! P uptake based on C uptake
  runp  =   min(  rump,  rupp+ misp)  ! actual uptake

  r  =   insw_vector(  runp)
  call flux_vector( iiPel, ppN1p,ppphytop, runp* r )  ! source/sink.p
  call flux_vector(iiPel, ppphytop,ppN1p,- runp*( 1.0D+00- r))  ! source/sink.p

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Excretion of N and P to PON and POP
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rr6n  =   pe_R6* sdo* phyton
  rr1n  =   sdo* phyton- rr6n

  rr6p  =   rr6c* p_qplc(phyto)
  rr1p  =   sdo* phytop- rr6p

  call flux_vector( iiPel, ppphyton,ppR1n, rr1n )  ! source/sink.n
  call flux_vector( iiPel, ppphyton,ppR6n, rr6n )  ! source/sink.n

  call flux_vector( iiPel, ppphytop,ppR1p, rr1p )  ! source/sink.p
  call flux_vector( iiPel, ppphytop,ppR6p, rr6p )  ! source/sink.p


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient dynamics: SILICATE
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  if ( phyto== iiP1) then

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    !  Nutrient uptake
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    rums  =   p_qus(phyto)* N5s(:)* phytoc  ! max pot uptake

    miss  =   sadap*( p_qsRc(phyto)* phytoc- phytos)  ! intracellular missing Si
    rups  =   run* p_qsRc(phyto)-( sdo+ srs)* phytos  ! Si uptake based on C uptake
    runs  =   min(  rums,  rups+ miss)  ! actual uptake
    r  =   insw_vector(  runs)
    call flux_vector( iiPel, ppN5s,ppphytos, runs* r )  ! source/sink.c
    call flux_vector(iiPel, ppphytos,ppN5s,- runs*( 1.0D+00- r))  ! source/sink.c

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Losses of Si
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    rr6s  =   sdo* phytos  ! Lysis, particulate

    ! Collect first all fluxes of P-->silica
    flP1R6s(:)  =   flP1R6s(:)+ rr6s

  end if



  if ( ChlLightFlag== 2) then
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Chl-a synthesis and photoacclimation
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


      WHERE ( ( phytoc> p_small) )
        rho_Chl = p_qchlc( phyto)* p_sum(phyto)* eiPI(phyto,:)* phytoc/( &
          p_alpha_chl(phyto)*( phytol+ p_small)* Irr)
        ! total synthesis, only when there is net production (run > 0)
!       rate_Chl = rho_Chl*( sum- sea- sra)* phytoc- sdo* phytol+ min( &
!         0.0D+00, sum- slc+ sdo)* max( 0.0D+00, phytol- p_qchlc( phyto)* phytoc)
         rate_Chl = rho_Chl*( max(srs,sum-slc) )* phytoc- sdo* phytol+ min( &
               -srs, sum- slc+ sdo)* max( 0.0D+00, phytol- p_qchlc( phyto)* phytoc)




      ELSEWHERE
        ! construction to overcome probelems with the GCC64bits compiler:
        rate_Chl = run* p_qchlc( phyto)+ phytoc*( p_qchlc( phyto)- &
          qlPc(phyto,:))


    END WHERE

    call flux_vector( iiPel, ppphytol,ppphytol, rate_Chl )
  end if



  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Sedimentation
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  if ( p_res(phyto)> 0.0D+00) then
    sediPI(phyto,:) = sediPI(phyto,:)+ p_res(phyto)* max( 0.0D+00, ( &
      p_esNI(phyto)- min( iN5s, iN)))
  end if


  ! End of computation section for process PhytoDynamics




  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
