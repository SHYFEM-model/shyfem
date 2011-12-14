#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: PhotoAvailableRadiation
!
! DESCRIPTION
!   This process computes the depth-integrated PAR and expresses
!	this in eiPI for each compartment (I). The change in EPLi
!	(optimal irradiance Iopt) due to daily variations is calculated in 
!	LightAdaptation.p.
!	The daily irradiance in W/m2 in each compartment is calculated
!	in the Light/Light.f and Light/VerticalDistribution.f and passed 
!       in EIR. The extinction-coefficient xEPS (/m) is calculated in 
!	in CalcVerticalExtinction.f
!
!   # Switch between diffferent light-production functions ("light type"):
!     iswLtyp = 0 : Steele (old ERSEM)  y*exp(1-y)
!     iswLtyp = 1 : Steele (Simpson)    y*exp(1-y)
!     iswLtyp = 2 : Ebenhoeh            2y/(1+y^2)
!     iswLtyp = 3 : ramp                min(1,y)
!     iswLtyp = 4 : step                1 if y>1 , 0 elsewhere
!	iswLtyp = 5 : Smith_average
!	iswLtyp = 6 : Smith II (actual_Irr)		
!
!     with y = irradiance/optimal
! 
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine PhotoAvailableRadiationDynamics(phyto, ppphytoc, ppphyton, &
    ppphytop, ppphytos, ppphytol)
!
! !USES:
  ! The following box states are used (NOT in fluxes): PhytoPlankton
  ! The following Pelagic 1-d global vars are used: D3STATETYPE
  ! The following Pelagic 1-d global boxvars  are used: EIR, xEPS, Depth
  ! The following Pelagic 2-d global boxvars are modified : EPLi
  ! The following Pelagic 2-d global boxvars got a value: eiPI
  ! The following constituent constants  are used: iiL
  ! The following 0-d global box parametes are used: LightForcingFlag
  ! The following global constants are used: RLEN,NOTRANSPORT

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,NOTRANSPORT
#ifdef NOPOINTERS
  use mem,  ONLY: D3STATE
#else
  use mem,  ONLY: D3STATE, PhytoPlankton
#endif
  use mem, ONLY: ppPhytoPlankton, D3STATETYPE, EIR, xEPS, Depth, EPLi, eiPI, &
    iiL, NO_BOXES, iiBen, iiPel, flux_vector
  use mem_Param,  ONLY: LightForcingFlag
  use mem_PhotoAvailableRadiation


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following vector functions are used:insw_vector
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,   ONLY: insw_vector


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
!   Original version by W. Ebenhoeh, Oldenburg University 
!                           Hanneke Baretta-Bekker, VKI
!       Translated to OpenSesame by Piet Ruardij, NIOZ
!	Dependency on phytoplankton species by M. Vichi, INGV
!
!
!
! !REVISION_HISTORY
!   File created on 8 feb. 1997
!	Modified by Daji Huang and JWB, 19/6/1998
!	Scaled the calculation of eiPI to a max. of 1 (at SUNQ of 16h) JWB040930
!	(Re)introduced the Smith formulations and removed the eiPI scaling JWB041006
!
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
  integer  :: i
  integer,dimension(NO_BOXES)  :: rampcontrol
  real(RLEN),dimension(NO_BOXES)  :: pIRRZ
  real(RLEN),dimension(NO_BOXES)  :: pIRR0
  real(RLEN),dimension(NO_BOXES)  :: xd
  real(RLEN),dimension(NO_BOXES)  :: exfac
  real(RLEN),dimension(NO_BOXES)  :: f_0_noon
  real(RLEN),dimension(NO_BOXES)  :: f_z_noon
  real(RLEN),dimension(NO_BOXES)  :: f_0_afternoon
  real(RLEN),dimension(NO_BOXES)  :: f_z_afternoon
  real(RLEN),dimension(NO_BOXES)  :: f_0_mean
  real(RLEN),dimension(NO_BOXES)  :: f_z_mean
  real(RLEN),dimension(NO_BOXES)  :: pirr0_noon
  real(RLEN),dimension(NO_BOXES)  :: pirrz_noon
  real(RLEN),dimension(NO_BOXES)  :: pirr0_afternoon
  real(RLEN),dimension(NO_BOXES)  :: pirrz_afternoon
  real(RLEN),dimension(NO_BOXES)  :: lx0
  real(RLEN),dimension(NO_BOXES)  :: lxz
  real(RLEN),dimension(NO_BOXES)  :: corr_mean
  real(RLEN),dimension(NO_BOXES)  :: corr_irra
  real(RLEN),dimension(NO_BOXES)  :: corr_irrb
  real(RLEN),dimension(NO_BOXES)  :: noon_light
  real(RLEN),dimension(NO_BOXES)  :: afternoon_light
  real(RLEN),dimension(NO_BOXES)  :: sum_state_phyto

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! user defined external functions
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer, external  :: rtoi
  integer, external  :: imin
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Copy  state var. object in local var
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  phytoc = D3STATE(ppphytoc,:)
  phyton = D3STATE(ppphyton,:)
  phytop = D3STATE(ppphytop,:)
  phytos = D3STATE(ppphytos,:)
  phytol = D3STATE(ppphytol,:)


  ! Recalculate Optimal light from the transported PLx.l

  i  =   ppPhytoPlankton(phyto,iiL)
  select case ( D3STATETYPE( i))

    case ( NOTRANSPORT )
      EPLi(phyto,:)  =   PhytoPlankton(phyto,iiL)



    case default
      EPLi(phyto,:)  =   PhytoPlankton(phyto,iiL)/ phytoc



  end select



  noon_light  =   EIR(:)* 1.7596D+00  ! magic number = 2 * 2PI/(4+PI)
  afternoon_light  =   EIR(:)* 1.0620D+00  ! magic number = (sqrt(2)+1)/2* 2PI/(4+PI)

  xd  =   xEPS(:)* Depth(:)
  exfac  =   exp( - xd)


  pIRR0  =   EIR(:)/ EPLi(phyto,:)
  pirr0_noon  =   noon_light/ EPLi(phyto,:)
  pirr0_afternoon  =   afternoon_light/ EPLi(phyto,:)

  pIRRZ  =   pIRR0* exfac
  pirrz_noon  =   pirr0_noon* exfac
  pirrz_afternoon  =   pirr0_afternoon* exfac

  select case ( p_iswLtyp(phyto))

    case ( 0 )
      ! Steele - Di Toro
      f_0_mean  =   exp(  1.0D+00- pIRR0)
      f_z_mean  =   exp(  1.0D+00- pIRRZ)
      corr_mean  =  ( f_z_mean- f_0_mean)/ xd
      corr_irra  =  ( f_z_mean- f_0_mean)/ xd* 6.0D+00
      corr_irrb  =   0.0D+00



    case ( 1 )
      ! Steele
      f_0_noon  =   exp(  1.0D+00- pirr0_noon)
      f_z_noon  =   exp(  1.0D+00- pirrz_noon)
      f_0_afternoon  =   exp(  1.0D+00- pirr0_afternoon)
      f_z_afternoon  =   exp(  1.0D+00- pirrz_afternoon)
      corr_irra  =  -( f_0_noon- f_z_noon)/ xd
      corr_irrb  =  -( f_0_afternoon- f_z_afternoon)/ xd



    case ( 2 )
      ! Ebenhoeh
      f_0_noon  =   atan(pirr0_noon)
      f_z_noon  =   atan(pirrz_noon)
      f_0_afternoon  =   atan(pirr0_afternoon)
      f_z_afternoon  =   atan(pirrz_afternoon)
      corr_irra  =   2.0D+00*( f_0_noon- f_z_noon)/ xd
      corr_irrb  =   2.0D+00*( f_0_afternoon- f_z_afternoon)/ xd



    case ( 3 )
      rampcontrol  =   2* int(insw_vector(  pirrz_noon- 1.0D+00))
      rampcontrol = min( 2, rampcontrol+ int(insw_vector( pirr0_noon- &
        1.0D+00)))

        WHERE (( rampcontrol)==2)
          corr_irra  =  ( log(  pirr0_noon)- log(  pirrz_noon))/ xd


        ELSEWHERE (( rampcontrol)==1)
          corr_irra  =  ( 1.0D+00+ log(  pirr0_noon)- pirrz_noon)/ xd


        ELSEWHERE (( rampcontrol)==0)
          corr_irra  =  ( pirr0_noon- pirrz_noon)/ xd


      END WHERE


      rampcontrol  =   2* int(insw_vector(  pirrz_afternoon- 1.0D+00))
      rampcontrol = min( 2, rampcontrol+ int(insw_vector( pirr0_afternoon- &
        1.0D+00)))


        WHERE (( rampcontrol)==2)
          corr_irrb  =   1.0D+00


        ELSEWHERE (( rampcontrol)==1)
          corr_irrb  =  ( 1.0D+00+ log(  pirr0_afternoon)- pirrz_afternoon)/ xd


        ELSEWHERE (( rampcontrol)==0)
          corr_irrb  =  ( pirr0_afternoon- pirrz_afternoon)/ xd


      END WHERE





    case ( 4 )
      !Step
      lx0  =   log(  pirr0_noon)
      lxz  =   log(  pirr0_afternoon)
      corr_irra  =  ( max(  0.0D+00,  lx0)- max(  0.0D+00,  lx0- xd))/ xd
      corr_irrb  =  ( max(  0.0D+00,  lxz)- max(  0.0D+00,  lxz- xd))/ xd



    case ( 5 )
      !Smith:
      f_0_mean = log( 1.0D+00+ sqrt( (1.0D+00/( 1.D-80+ pIRR0* &
        exp( 1.0D+00)))**(2.0D+00)+ 1.0D+00))
      f_z_mean = log( exfac+ sqrt( (1.0D+00/( 1.D-80+ pIRR0* &
        exp( 1.0D+00)))**(2.0D+00)+ exfac* exfac))
      corr_mean  =  ( f_0_mean- f_z_mean)/ xd
      corr_irra  =  ( f_0_mean- f_z_mean)/ xd* 6.0D+00
      corr_irrb  =   0.0D+00



    case ( 6 )
      ! Smith II
      f_0_noon = log( 1.0D+00+ sqrt( (pirr0_noon* exp( 1.0D+00))**(- &
        2.0D+00)+ 1.0D+00))
      f_z_noon = log( exfac+ sqrt( (pirr0_noon* exp( 1.0D+00))**(- &
        2.0D+00)+ exfac* exfac))
      f_0_afternoon = log( 1.0D+00+ sqrt( (pirr0_afternoon* exp( &
        1.0D+00))**(- 2.0D+00)+ 1.0D+00))
      f_z_afternoon = log( exfac+ sqrt( (pirr0_afternoon* exp( 1.0D+00))**(- &
        2.0D+00)+ exfac* exfac))
      corr_irra  =  ( f_0_noon- f_z_noon)/ xd
      corr_irrb  =  ( f_0_afternoon- f_z_afternoon)/ xd



  end select




  select case ( LightForcingFlag)

    case ( 1 )
      !rectangular integration:
      eiPI(phyto,:)  =   min(  1.0D+00,  corr_mean)



    case ( 3 )
      !   Simpson integration is used as default:
      eiPI(phyto,:)  =  ( corr_irra+ 4.0D+00* corr_irrb)/ 6.0D+00


  end select







  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
