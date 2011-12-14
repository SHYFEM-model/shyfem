#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Initialise BFM variables
!
! !INTERFACE:
   subroutine init_var_bfm(namlst,fname,unit,setup)
!
! !DESCRIPTION:
!  Allocate BFM variables and give initial values of
!  parameters and state variables
!
! !USES:
#ifndef NOT_STANDALONE
   use api_bfm
   use global_mem
#endif
#ifdef BFM_GOTM
   use bio_var
   use bio_bfm
#endif
#ifdef BFM_POM
   use api_pom
#endif
#ifdef BFM_OPA_OFFLINE
   use api_opa_offline
#endif
#ifdef BFM_OPA_PELAGOS
   use api_opa_pelagos
#endif
   use mem
   use mem_Phyto, ONLY: p_qnRc,p_qpRc,p_qsRc
   use mem_BenBac, ONLY: p_qnc,p_qpc
   use constants, ONLY: HOURS_PER_DAY
   use mem_Param, ONLY: CalcPelagicFlag,CalcBenthicFlag,p_small,p_qchlc, &
                        CalcPhytoPlankton,CalcMicroZooPlankton,          &
                        CalcMesoZooPlankton,CalcBenOrganisms,            &
                        CalcBenBacteria,CalcBacteria,CalcPelChemistry,   &
                        p_small
   use string_functions, ONLY: getseq_number,empty
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,          intent(in)        :: namlst
   character(len=*), intent(in)        :: fname
   integer,          intent(in)        :: unit
   integer,          intent(in)        :: setup
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
!
! !LOCAL VARIABLES:
   integer              :: icontrol,i,j,iiLastElement
   integer,parameter    :: NSAVE=100  ! Maximum no variables which can be saved
   character(len=64),dimension(NSAVE):: var_save
   character(len=64),dimension(NSAVE):: ave_save
   REALTYPE  :: n1p0,n3n0,n4n0,n5s0,n6r0,  &
                p1c0,p2c0,p3c0,p4c0,z3c0,  &
                z4c0,z5c0,z6c0,b1c0,r1c0,  &
                r2c0,r6c0,r7c0,o2o0,o4n0,  &
                p1l0,p2l0,p3l0,p4l0

   REALTYPE  :: y1c0, y2c0, y3c0, y4c0, y5c0, &
                q1c0, q11c0, q6c0, q1n0,      &
                q11n0, q6n0, q1p0, q11p0,     &
                q6p0, q6s0, k3n0, g4n0, &
                h1c0, h2c0, k1p0, k11p0, k21p0,     &
                k4n0, k14n0, k24n0, k6r0,k5s0, &
                d1m0, d2m0, d6m0, d7m0, d8m0, d9m0, g2o0
   namelist /bfm_init_nml/ surface_flux_method,       &
                           n_surface_fluxes,          &
                           n1p0,n3n0,n4n0,n5s0,n6r0,  &
                           p1c0,p2c0,p3c0,p4c0,z3c0,  &
                           z4c0,z5c0,z6c0,b1c0,r1c0,  &
                           r2c0,r6c0,r7c0,o2o0,o4n0,  &
                           p1l0,p2l0,p3l0,p4l0

   namelist /bfm_save_nml/ var_save, ave_save

   namelist /bfm_ben_init_nml/  &
                           y1c0, y2c0, y3c0, y4c0, y5c0,     &
                           q1c0, q11c0, q6c0, q1n0,          &
                           q11n0, q6n0, q1p0, q11p0,         &
                           q6p0,  q6s0, k3n0, g4n0,          &
                           h1c0, h2c0, k1p0, k11p0, k21p0,   &
                           k4n0, k14n0, k24n0, k6r0,k5s0,    &
                           d1m0, d2m0, d6m0, d7m0, d8m0, d9m0, g2o0
   interface
      subroutine init_cnps(c,n,p,s,nc,pc,sc)
         REALTYPE,dimension(:),intent(in)           :: c
         REALTYPE,intent(in),optional               :: nc,pc,sc
         REALTYPE,dimension(:),intent(out),optional :: n
         REALTYPE,dimension(:),intent(out),optional :: p
         REALTYPE,dimension(:),intent(out),optional :: s
      end subroutine init_cnps
   end interface
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
!EOP
!-----------------------------------------------------------------------
!BOC

   LEVEL2 'init_var_bfm'
   !---------------------------------------------
   ! Give reasonable initial values
   ! Overwritten by namelist parameters
   !---------------------------------------------
   surface_flux_method = -1
   n_surface_fluxes = 1

   !---------------------------------------------
   ! Pelagic variables
   !---------------------------------------------
   n1p0 = _ONE_
   n3n0 = _ONE_
   n4n0 = _ONE_
   n5s0 = _ONE_
   n6r0 = _ONE_
   o2o0 = 300.0
   o4n0 = _ONE_
   p1c0 = _ONE_
   p2c0 = _ONE_
   p3c0 = _ONE_
   p4c0 = _ONE_
   p1l0 = _ONE_
   p2l0 = _ONE_
   p3l0 = _ONE_
   p4l0 = _ONE_
   z3c0 = _ONE_
   z4c0 = _ONE_
   z5c0 = _ONE_
   z6c0 = _ONE_
   b1c0 = _ONE_
   r1c0 = _ONE_
   r2c0 = _ONE_
   r6c0 = _ONE_
   r7c0 = _ONE_

   !---------------------------------------------
   ! Benthic variables
   !---------------------------------------------
   y1c0  = _ONE_
   y2c0  = _ONE_
   y3c0  = _ONE_
   y4c0  = _ONE_
   y5c0  = _ONE_
   q1c0  = _ONE_
   q11c0 = _ONE_
   q6c0  = _ONE_
   q1n0  = _ONE_
   q11n0 = _ONE_
   q6n0  = _ONE_
   q1p0  = _ONE_
   q11p0 = _ONE_
   q6p0  = _ONE_
   q6s0  = _ONE_
   h1c0  = _ONE_
   h2c0  = _ONE_
   k1p0  = _ONE_
   k11p0 = _ONE_
   k21p0 = _ONE_
   k3n0  = _ONE_
   g4n0  = _ONE_
   k4n0  = _ONE_
   k14n0 = _ONE_
   k24n0 = _ONE_
   k6r0  = _ONE_
   k5s0  = _ONE_
   d1m0  = _ONE_
   d2m0  = _ONE_
   d6m0  = _ONE_
   d7m0  = _ONE_
   d8m0  = _ONE_
   d9m0  = _ONE_
   g2o0  = _ONE_

   !---------------------------------------------
   ! Open and read the namelist
   !---------------------------------------------
   icontrol=0
   open(namlst,file=fname,action='read',status='old',err=98)
   read(namlst,nml=bfm_init_nml,err=99)
   if (setup >=2 )  then
     read(namlst,nml=bfm_ben_init_nml,err=101)
   end if
   var_save=""
   ave_save=""
   var_ave=.false.
   read(namlst,nml=bfm_save_nml,err=100)
   close(namlst)
   icontrol=1
98 if ( icontrol == 0 ) then
     LEVEL3 'I could not open ',trim(fname)
     LEVEL3 'The initial values of the BFM variables are set to ONE'
     LEVEL3 'If thats not what you want you have to supply ',trim(fname)
   end if

   !---------------------------------------------
   ! Check variable to be saved and
   ! set the corresponding flag value in var_ids
   !---------------------------------------------
   do i=1,NSAVE
      if (.NOT.empty(var_save(i))) then
            j=getseq_number(var_save(i),var_names,stBenFluxE,.TRUE.)
            if ( j > 0 ) var_ids(j)=-1
      end if
      if ( .NOT.empty(var_save(i)) .AND. j==0 ) then
            STDERR 'Warning: variable ',trim(var_save(i)),' does not exist!'
      end if
   end do
   do i=1,NSAVE
      if (.NOT.empty(ave_save(i))) then
         j=getseq_number(ave_save(i),var_names,stBenFluxE,.TRUE.)
         if ( .NOT.empty(ave_save(i)) .AND. j==0 ) then
            STDERR 'Warning: variable ',trim(ave_save(i)),' does not exist!'
         else if ( var_ids(j) <0 ) then
            STDERR 'Warning: Variable ',trim(ave_save(i)), &
               ' is already selected for output in var_save'
         else if ( j > 0 ) then
            var_ids(j)=-1
            var_ave(j)=.true.
         end if
      end if
   end do

#ifdef BFM_GOTM
   !---------------------------------------------
   ! Create pointers
   !---------------------------------------------
    call pointers_gotm_bfm()
#endif

   !---------------------------------------------
   ! Initialize BFM parameters
   !---------------------------------------------
   call Initialize

   !---------------------------------------------
   ! Initially set the number of sun hours
   ! equal to the number of hours in a day.
   !---------------------------------------------
   SUNQ = HOURS_PER_DAY

   !---------------------------------------------
   ! Initialise pelagic state variables
   ! also if using a benthic-only setup
   ! (for boundary conditions)
   !---------------------------------------------
      N1p = n1p0
      N3n = n3n0
      N4n = n4n0
      N5s = n5s0
      N6r = n6r0
      O2o = o2o0
      O4n = o4n0
      P1c = p1c0
      if (p1l0 /= _ONE_) then
         P1l = p1l0
      else
         P1l = p1c0*p_qchlc(iiP1)
      end if
      P2c = p2c0
      if (p2l0 /= _ONE_) then
         P2l = p2l0
      else
         P2l = p2c0*p_qchlc(iiP2)
      end if
      P3c = p3c0
      if (p3l0 /= _ONE_) then
         P3l = p3l0
      else
         P3l = p3c*p_qchlc(iiP3)
      end if
      P4c = p4c0
      if (p4l0 /= _ONE_) then
         P4l = p4l0
      else
         P4l = p4c0*p_qchlc(iiP4)
      end if
      Z3c = z3c0
      Z4c = z4c0
      Z5c = z5c0
      Z6c = z6c0
      B1c = b1c0
      R1c = r1c0
      R2c = r2c0
      R6c = r6c0
      R7c = r7c0

      !---------------------------------------------
      ! Initialise other internal components
      ! with Redfield
      !---------------------------------------------
      call init_cnps(c=P1c,n=P1n,p=P1p,s=P1s,nc=p_qnRc(iiP1), &
           pc=p_qpRc(iiP1),sc=p_qsRc(iiP1))
      call init_cnps(c=P2c,n=P2n,p=P2p,nc=p_qnRc(iiP2), &
           pc=p_qpRc(iiP2))
      call init_cnps(c=P3c,n=P3n,p=P3p,nc=p_qnRc(iiP3), &
           pc=p_qpRc(iiP3))
      call init_cnps(c=P4c,n=P4n,p=P4p,nc=p_qnRc(iiP4), &
           pc=p_qpRc(iiP4))
      call init_cnps(c=Z3c,n=Z3n,p=Z3p)
      call init_cnps(c=Z4c,n=Z4n,p=Z4p)
      call init_cnps(c=Z5c,n=Z5n,p=Z5p)
      call init_cnps(c=Z6c,n=Z6n,p=Z6p)
      call init_cnps(c=B1c,n=B1n,p=B1p)
      call init_cnps(c=R1c,n=R1n,p=R1p)
      call init_cnps(c=R6c,n=R6n,p=R6p,s=R6s)

   !---------------------------------------------
   ! Initialise benthic state variables
   !---------------------------------------------
   !MAV: need to always give init non-zero values
   ! because there are still part of the
   ! benthic system which are computed when setup=1
!   if (setup >=2) then
      Y1c  = y1c0
      Y2c  = y2c0
      Y3c  = y3c0
      Y4c  = y4c0
      Y5c  = y5c0
      Q1c  = q1c0
      Q11c = q11c0
      Q6c  = q6c0
      Q1n  = q1n0
      Q11n = q11n0
      Q6n  = q6n0
      Q1p  = q1p0
      Q11p = q11p0
      Q6p  = q6p0
      Q6s  = q6s0
      H1c  = h1c0
      H2c  = h2c0
      K1p  = k1p0
      K11p = k11p0
      K21p = k21p0
      K3n  = k3n0
      G4n  = g4n0
      K4n  = k4n0
      K14n = k14n0
      K24n = k24n0
      K6r  = k6r0
      K5s  = k5s0
      D1m  = d1m0
      D2m  = d2m0
      D6m  = d6m0
      D7m  = d7m0
      D8m  = d8m0
      D9m  = d9m0
      G2o  = g2o0

      !---------------------------------------------
      ! Initialise organisms' internal components
      ! with Redfield
      !---------------------------------------------
      call init_cnps(c=Y1c,n=Y1n,p=Y1p)
      call init_cnps(c=Y2c,n=Y2n,p=Y2p)
      call init_cnps(c=Y3c,n=Y3n,p=Y3p)
      call init_cnps(c=Y4c,n=Y4n,p=Y4p)
      call init_cnps(c=Y5c,n=Y5n,p=Y5p)
      call init_cnps(c=H1c,n=H1n,p=H1p,nc=p_qnc(iiH1), &
           pc=p_qpc(iiH1))
      call init_cnps(c=H2c,n=H2n,p=H2p,nc=p_qnc(iiH2), &
           pc=p_qpc(iiH2))
!    end if

   !---------------------------------------------
   ! Check setup settings
   ! and finalize initialization
   !---------------------------------------------
   select case (setup)
      case (0)
      case (1) ! Pelagic only
         LEVEL2 "Pelagic-only setup (bio_setup=1), Switching off the benthic system"
         CalcBenthicFlag = 0
      case (2) ! Benthic only
         LEVEL2 "Benthic-only setup (bio_setup=2), Switching off the pelagic system"
         CalcPelagicFlag = .FALSE.
         CalcPhytoPlankton=.FALSE.
         CalcBacteria=.FALSE.
         CalcMesoZooPlankton=.FALSE.
         CalcMicroZooPlankton=.FALSE.
      case (3) ! Pelagic-Benthic coupling
         LEVEL2 "Pelagic-Benthic coupled setup (bio_setup=3)"
         if (CalcBenthicFlag == 0) &
            LEVEL3 'Warning, benthic system is switched off!'
         if (.NOT.CalcPelagicFlag) &
            LEVEL3 'Warning, pelagic system is switched off!'
   end select

   select case (CalcBenthicFlag)
     case (0)
        LEVEL3 "Benthic model is: not used"
     case (1)
        LEVEL3 "Benthic model is: simple nutrient return"
     case (2)
        LEVEL3 "Benthic model is: benthos + intermediate nutrient return"
     case (3)
        LEVEL3 "Benthic model is: benthos + Ruardij & Van Raaphorst"
   end select

   !---------------------------------------------
   ! Zeroing of the switched off state variables
   !---------------------------------------------
   do j = 1,iiPhytoPlankton
      iiLastElement = iiL
      if (.NOT.CalcPhytoPlankton(j)) then
         if (j==iiP1) iiLastElement=iiS
         do i = iiC,iiLastElement
            D3STATE(ppPhytoPlankton(j,i),:) = p_small
            D3STATETYPE(ppPhytoPlankton(j,i)) = NOTRANSPORT
         end do
      end if
   end do
   do j = 1,iiMesoZooPlankton
      iiLastElement = iiP
      if (.NOT.CalcMesoZooPlankton(j)) then
         do i = iiC,iiLastElement
            D3STATE(ppMesoZooPlankton(j,i),:) = p_small
            D3STATETYPE(ppMesoZooPlankton(j,i)) = NOTRANSPORT
         end do
      end if
   end do
   do j = 1,iiMicroZooPlankton
      iiLastElement = iiP
      if (.NOT.CalcMicroZooPlankton(j)) then
         do i = iiC,iiLastElement
            D3STATE(ppMicroZooPlankton(j,i),:) = p_small
            D3STATETYPE(ppMicroZooPlankton(j,i)) = NOTRANSPORT
         end do
      end if
   end do
   do j = 1,iiBenOrganisms
      iiLastElement = iiP
      if (.NOT.CalcBenOrganisms(j)) then
         do i = iiC,iiLastElement
            D2STATE(ppBenOrganisms(j,i),:) = p_small
         end do
      end if
   end do
   do j = 1,iiBenBacteria
      iiLastElement = iiP
      if (.NOT.CalcBenBacteria(j)) then
         do i = iiC,iiLastElement
            D2STATE(ppBenBacteria(j,i),:) = p_small
         end do
      end if
   end do
   if (.NOT.CalcBacteria) then
      B1c = p_small; B1n = p_small; B1p = p_small;
      D3STATETYPE(ppB1c) = NOTRANSPORT
      D3STATETYPE(ppB1n) = NOTRANSPORT
      D3STATETYPE(ppB1p) = NOTRANSPORT
   end if

   return

99  FATAL 'I could not read bfm_init_nml'
    stop 'init_var_bfm'
100 FATAL 'I could not read bfm_save_nml'
    stop 'init_var_bfm'
101 FATAL 'I could not read bfm_ben_init_nml'
    stop 'init_var_bfm'
102 FATAL 'I could not read bfm_ben_save_nml'
    stop 'init_var_bfm'

   end subroutine init_var_bfm
!EOC

