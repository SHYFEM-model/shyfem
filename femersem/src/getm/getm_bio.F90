!$Id$
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: getm_bio()
!
! !INTERFACE:
   module getm_bio
!
! !DESCRIPTION:
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,ioff,joff
   use domain, only: iimin,iimax,jjmin,jjmax,kmax
   use domain, only: az,au,av
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dxu,dxv,dyu,dyv,arcd1
#else
   use domain, only: dx,dy,ard1
#endif
   use variables_3d, only: uu,vv,ww,hun,hvn,ho,hn
   use variables_3d, only: nuh,T,S,light
   use variables_3d, only: cc3d,ccb3d, &
                           cc3d_out,ccb3d_out, counter_ave,flag_out, & 
                           p_poro_2d,p_p_ae_2d,p_poro_default,p_p_ae_default, &
                           ben_init_file 
   use advection_3d, only: do_advection_3d
   use meteo, only: swr
   use halo_zones, only: update_3d_halo,update_2d_halo, wait_halo,D_TAG
   use bio, only: init_bio, do_bio
   use bio, only: bio_calc
   use bio_var, only: numc,cc, &
                  numbc,ccb, &
                  numc_diag,diag, numbc_diag,diagb, bio_setup, &
                  var_ids,var_ave, &
                  stPelStateS,stPelDiagS,stPelFluxS,stBenStateS,stBenDiagS,stBenFluxS, & 
                  stPelStateE,stPelDiagE,stPelFluxE,stBenStateE,stBenDiagE,stBenFluxE      !BFM
   use gotm_error_msg, only: set_parallel_flag_for_gotm, output_gotm_error,get_warning_for_getm
   use exceptions, only:getm_error
   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
   public init_getm_bio, do_getm_bio, fill_diag_array
   logical, public           :: hotstart_bio=.true.
!
! !PRIVATE DATA MEMBERS:
   integer         :: bio_hor_adv=1
   integer         :: bio_ver_adv=1
   integer         :: bio_adv_split=1
   REALTYPE        :: bio_AH=10.
   REALTYPE        :: bio_missing=-9999.0
   logical         :: read_poro=.FALSE.
#ifdef STATIC
   REALTYPE        :: delxu(I2DFIELD),delxv(I2DFIELD)
   REALTYPE        :: delyu(I2DFIELD),delyv(I2DFIELD)
   REALTYPE        :: area_inv(I2DFIELD)
   REALTYPE        :: ffp(I3DFIELD)
   REALTYPE        :: ffb(I2DFIELD)
!  REALTYPE, dimension(:,:,:), pointer :: ff
#else
   REALTYPE, dimension(:,:), allocatable :: delxu,delxv
   REALTYPE, dimension(:,:), allocatable :: delyu,delyv
   REALTYPE, dimension(:,:), allocatable :: area_inv

#endif
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!  $Log$
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_getm_bio
!
! !INTERFACE:
   subroutine init_getm_bio(input_dir)

   IMPLICIT NONE

!
! !DESCRIPTION:! !INPUT PARAMETERS:
  character(len=*)                    :: input_dir

!  Reads the namelist and makes calls to the init functions of the
!  various model components.
!
! !REVISION HISTORY:
!  See the log for the module
!
!  !LOCAL VARIABLES
   integer, parameter                  :: namlst=10
   integer, parameter                  :: unit_bio=63
   REALTYPE                            :: h(0:kmax)
   character(len=80)                    :: brrbrr

   integer         :: rc
   integer         :: i,j,n

   namelist /getm_bio_nml/ hotstart_bio,bio_hor_adv,bio_ver_adv, &
                           bio_adv_split,bio_AH,read_poro,ben_init_file
!EOP
!-------------------------------------------------------------------------
!BOC
   LEVEL2 'init_getm_bio()'

   call init_bio(11,'bio.inp',unit_bio,kmax,h)

   if (bio_calc) then

      if ( bio_setup /=2 ) then
        if ( numc > 0 ) then 
          allocate(cc3d(numc,I3DFIELD),stat=rc)       ! pel. biological fields of states
          if (rc /= 0) stop 'init_getm_bio: Error allocating memory (cc3d)'
        endif
        n=count(var_ave(stPelStateS:stPelStateE) ==1) + &
                count( var_ids(stPelDiagS:stPelFluxE)/= 0 )
        if ( n > 0 )  then
          allocate(cc3d_out(n,I3DFIELD),stat=rc)   ! pel. biological fields of diagnos.
          if (rc /= 0) stop 'init_getm_bio: Error allocating memory (cc3d_out)'
          cc3d_out=bio_missing;
        endif
        allocate(counter_ave(I2DFIELD),stat=rc)   ! pel. biological fields of diagnos.
        if (rc /= 0) stop 'init_getm_bio: Error allocating memory (counter_ave)'
        counter_ave=0;
      endif

      if ( bio_setup >=2) then
        if ( numbc > 0 )  then
          allocate(ccb3d(numbc,I2DFIELD,0:1),stat=rc) ! bent, biological fields of states
          if (rc /= 0) stop 'init_getm_bio: Error allocating memory (ccb3d)'
        endif
        n=count( var_ave(stBenStateS:stBenStateE)==1) + &
                count( var_ids(stBenDiagS:stBenFluxE) /= 0)
        if ( n > 0 )  then
          allocate(ccb3d_out(n,I2DFIELD,0:1),stat=rc)   ! ben. biological fields of flux.
          if (rc /= 0) stop 'init_getm_bio: Error allocating memory (cc3d_out)'
          ccb3d_out=bio_missing;
        endif
      endif
      

      ben_init_file=''
      brrbrr="Reading from "//trim(input_dir)//"/getm_bio.inp"
      LEVEL2 brrbrr
      open(77,status='unknown',file=trim(input_dir) // "/getm_bio.inp")
!     open(77,status='unknown',file="/home/kbk/getm-cases/not_ready/ns_06nm/getm_bio.inp")
      read(77,NML=getm_bio_nml)

      LEVEL2 "Settings related to 3D biological calculations"
      LEVEL3 'bio_hor_adv=   ',bio_hor_adv
      LEVEL3 'bio_ver_adv=   ',bio_ver_adv
      LEVEL3 'bio_adv_split= ',bio_adv_split
      LEVEL3 'bio_AH=        ',bio_AH

      if (hotstart_bio) then
         LEVEL2 "Reading biological fields from hotstart file"
      else
         LEVEL2 "Initialise biological fields from namelist"
         do j=jjmin,jjmax
           do i=iimin,iimax
              if (az(i,j) .ne. 0 ) then
                if ( bio_setup /= 2 ) cc3d(1:numc,i,j,:)=cc(1:numc,:)
                if ( bio_setup >= 2 ) ccb3d(1:numbc,i,j,:)=ccb(1:numbc,:)
              endif
           enddo
         end do
         call init_2d_grid(ben_init_file)
      end if
      call set_2d_grid_prameters(0,0,0)

#ifndef STATIC
      allocate(delxu(I2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_getm_bio: Error allocating memory (delxu)'

      allocate(delxv(I2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_getm_bio: Error allocating memory (delxv)'

      allocate(delyu(I2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_getm_bio: Error allocating memory (delyu)'

      allocate(delyv(I2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_getm_bio: Error allocating memory (delyv)'

      allocate(area_inv(I2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_getm_bio: Error allocating memory (area_inv)'
#endif
#if defined(SPHERICAL) || defined(CURVILINEAR)
      delxu=dxu
      delxv=dxv
      delyu=dyu
      delyv=dyv
      area_inv=arcd1
#else
      delxu=dx
      delxv=dx
      delyu=dy
      delyv=dy
      area_inv=ard1
#endif

      do n=1,numc
         LEVEL3 'n=',n
         call update_3d_halo(cc3d(n,:,:,:),cc3d(n,:,:,:),az, &
                          iimin,jjmin,iimax,jjmax,kmax,D_TAG)
         call wait_halo(D_TAG)
      enddo

   end if
#ifdef PARALLEL
   call set_parallel_flag_for_gotm(.TRUE.)
#endif


   return
   end subroutine init_getm_bio
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  do_getm_bio()
!
! !INTERFACE:
   subroutine do_getm_bio(dt,write_3d)
!
! !DESCRIPTION:
!
! !USES:
   use meteo, only:swr,u10,v10                                     !BFM
!BFM  use meteo, only:swr
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
   logical, intent(in)                 :: write_3d
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  See the log for the module
!
! !LOCAL VARIABLES:
   integer         :: n
   integer         :: i,j,k
   REALTYPE        :: h1d(0:kmax),T1d(0:kmax),S1d(0:kmax)
   REALTYPE        :: nuh1d(0:kmax),light1d(0:kmax)
   REALTYPE        :: bioshade1d(0:kmax)
   REALTYPE        :: I_0,wind
   REALTYPE        :: r
   CHARACTER(LEN=80) :: msg,sub
   LOGICAL         :: error_flag,warning_flag
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL2 'do_getm_bio()'

   LEVEL3 'do advection processes'
!  then we do the advection of the biological variables
   if ( .not.bio_calc) return
   if ( bio_setup /= 2 )  then
     do n=1,numc
        ffp = cc3d(n,:,:,:)
        call do_advection_3d(dt,ffp,uu,vv,ww,hun,hvn,ho,hn, &
              delxu,delxv,delyu,delyv,area_inv,az,au,av, &
              bio_hor_adv,bio_ver_adv,bio_adv_split,bio_AH)
        call update_3d_halo(ffp,ffp,az, &
                         iimin,jjmin,iimax,jjmax,kmax,D_TAG)
        call wait_halo(D_TAG)
        cc3d(n,:,:,:) = ffp
     end do
   endif

!  First we do all the vertical processes
   LEVEL3 'do vertical processes'
   do j=jjmin,jjmax
      do i=iimin,iimax
         if (az(i,j) .ge. 1 ) then
          h1d=hn(i,j,:)
          if ( h1d(1) > 0.1D+00 ) then
            I_0=swr(i,j)
            T1d=T(i,j,:)
            S1d=S(i,j,:)
            nuh1d=nuh(i,j,:)
            light1d=light(i,j,:)
            call set_2d_grid_prameters(1,i,j)
            if ( bio_setup /= 2 ) cc(:,:)=cc3d(:,i,j,:)
            if ( bio_setup >= 2 ) ccb(:,:)=ccb3d(:,i,j,:)
            wind = sqrt(u10(i,j)*u10(i,j)+v10(i,j)*v10(i,j))                     !BFM
            call do_bio(kmax,I_0,wind,dt,h1d,T1d,S1d,nuh1d,light1d,bioshade1d)
            if ( bio_setup /= 2 ) cc3d(:,i,j,:)=cc(:,:)
            if ( bio_setup >= 2 ) ccb3d(:,i,j,:)=ccb(:,:)
            light(i,j,:)=bioshade1d
            call fill_diag( 1,write_3d, i,j, h1d,kmax)
            call output_gotm_error( error_flag, sub, msg)
            if ( error_flag ) then
               r=sum(h1d)
               STDERR 'i,j=',i+ioff,j+joff,'depth=',r
               call getm_error( sub, msg)
            endif
            call get_warning_for_getm(warning_flag )
            if ( warning_flag ) then
               r=sum(h1d)
               STDERR 'i,j=',i+ioff,j+joff,'depth=',r
            endif
          else
            STDERR 'shallow gridpoint i,j=',i+ioff,j+joff,'depth=',r
          endif
         end if
      end do
   end do

   LEVEL2 'end of vertical processes'


   return
   end subroutine do_getm_bio
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE:  fill_diag
! 
! !INTERFACE:
   subroutine fill_diag( mode,llwrite, ig,jg, h,nlev )
!
! !USES:
  use bio_var, only: stPelStateS,stPelDiagS,stPelFluxS,stBenStateS,stBenDiagS,stBenFluxS
  use bio_var, only: stPelStateE,stPelDiagE,stPelFluxE,stBenStateE,stBenDiagE,stBenFluxE

!
! !INPUT PARAMETERS: 
   implicit none 
   logical,intent(IN)                        :: llwrite
   integer,intent(IN)                        :: mode
   integer,intent(IN)                        :: ig
   integer,intent(IN)                        :: jg
   integer,intent(IN)                        :: nlev
   REALTYPE,intent(IN),dimension(0:nlev)     :: h
!
!
! !DESCRIPTION: 
!      With this routine all the output prepared for all non-state variables 
!      and for state variables of which the average value hato be collected
!      between 2 subsequent output steps.
!
!      All values are stored in the array cc3d_out :
!      The squenbces in the array is:
!       1, average values:
!                      a. state variables
!                      b. diagnostic variables which are calculated in the model.
!                      c. flux variables.
!      2. normal values:
!                      b. diagnostic variables
!                      c. normal variables. 
! !BUGS:  
!
! !SEE ALSO: 
!
! !SYSTEM ROUTINES: 
!
! !FILES USED:  
!
! !REVISION HISTORY: 
!
! 
!EOP
!-----------------------------------------------------------------------
!BOC
  integer                         ::k
  integer                         ::l

   select case (mode)
    case(0)
        counter_ave(iimin:iimax,jjmin:jjmax)=0
    case(1)
        k=0
        ! add value for average output:
        call fill_diag_array(1,1, stPelStateS,stPelStateE,ig,jg, h,nlev, k )
        call fill_diag_array(2,1, stPelDiagS,stPelDiagE,ig,jg, h,nlev, k )
        call fill_diag_array(3,1, stPelFluxS,stPelFluxE,ig,jg, h,nlev, k )
        if ( llwrite) then
          ! set value for normal output aftyer the average output:
          call fill_diag_array(2,0, stPelDiagS,stPelDiagE,ig,jg, h,nlev, k )
          call fill_diag_array(3,0, stPelFluxS,stPelFluxE,ig,jg, h,nlev, k )
        endif
        k=0
        ! add value for average output:
        call fill_diag_array(4,1, stBenStateS,stBenStateE,ig,jg, h,nlev, k )
        call fill_diag_array(5,1, stBenDiagS,stBenDiagE,ig,jg, h,nlev, k )
        call fill_diag_array(6,1, stBenFluxS,stBenFluxE,ig,jg, h,nlev, k )
        if ( llwrite) then
          ! set value for normal output after the average output:
          call fill_diag_array(5,0, stBenDiagS,stBenDiagE,ig,jg, h,nlev, k )
          call fill_diag_array(6,0, stBenFluxS,stBenFluxE,ig,jg, h,nlev, k )
        endif
        counter_ave(ig,jg)=counter_ave(ig,jg)+1
        flag_out=1
    case(2)
        ! calculate the acerage value: divide values for average output by countere_ave
        k=0
        call fill_diag_array(11,1, stPelStateS,stPelFluxE,0,0, h,nlev, k )
        call fill_diag_array(11,0, stPelDiagS,stPelFluxE,0,0, h,nlev, k )
        k=0
        call fill_diag_array(14,1, stBenStateS,stBenFluxE,0,0, h,nlev, k )
        call fill_diag_array(14,0, stBenDiagS,stBenFluxE,0,0, h,nlev, k )
    end select

   end subroutine fill_diag

!EOC
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: fill_diag_array 
! 
! !INTERFACE:
!
   subroutine fill_diag_array(mode,ave, from,to,ig,jg, h,nlev, k )
! !USES:
   use bio_var, only: diag, diagb,var_ids,var_ave

   use bio_var,only: c1dim
   use mem,only:make_flux_output
!
! !INPUT PARAMETERS: 
!
   implicit none 
   integer,intent(IN)                        :: mode
   integer,intent(IN)                        :: from
   integer                                   :: ave
   integer,intent(IN)                        :: to
   integer,intent(IN)                        :: ig
   integer,intent(IN)                        :: jg
   integer,intent(IN)                        :: nlev
   REALTYPE,intent(IN),dimension(0:nlev)     :: h
! !INPUT/OUTPUT PARAMETERS:
!
   integer,intent(INOUT)                     :: k

! !DESCRIPTION: 
!   fill 4d array with columns calculation in an efficient way
!
! !LOCAL VARIABLES:
!
   integer                        :: i
   integer                        :: j
   integer                        :: l
   integer                        :: m
   integer                        :: n
   integer                        :: k2

! !BUGS:  
!
! !SEE ALSO: 
!
! !SYSTEM ROUTINES: 
!
! !FILES USED:  
!
! !REVISION HISTORY: 
!
!  2006-05-22   Piet Ruardij  Initial code.
! 
!EOP
!-----------------------------------------------------------------------
!BOC
   select case (mode)
      case(10:)
        if ( flag_out ==0) return
        do i=from,to
          if ( ave==1 ) then
            if ( var_ids(i) /= 0 .and. var_ave(i) == 1) then 
              k=k+1
              select case(mode)
                case(11)
                  do n=jjmin,jjmax
                    do m=iimin,iimax
                      if (az(m,n) .ge. 1 ) then
                        if (counter_ave(m,n).gt.0 ) then
                           cc3d_out(k,m,n,0:nlev) = cc3d_out(k,m,n,0:nlev)/counter_ave(m,n)
                        else
                           cc3d_out(k,m,n,0:nlev)=bio_missing
                        endif
                      endif
                    enddo
                  enddo 
                case(14)
                  do n=jjmin,jjmax
                    do m=iimin,iimax
                      if (az(m,n) .ge. 1 ) then
                        if (counter_ave(m,n).gt.0 ) then
                           ccb3d_out(k,m,n,0:1) = ccb3d_out(k,m,n,0:1)/counter_ave(m,n)
                        else
                           ccb3d_out(k,m,n,0:1)=bio_missing
                        endif
                      endif
                    enddo
                  enddo 
              end select
            endif
          elseif ( var_ids(i) /= 0 .and. var_ave(i) == 0) then 
            k=k+1
            select case(mode)
              case(11)
                do n=jjmin,jjmax
                  do m=iimin,iimax
                    if (az(m,n) .ge. 1 ) then
                      if (counter_ave(m,n) ==0 ) cc3d_out(k,m,n,0:nlev)=bio_missing
                    endif
                  enddo
                enddo 
              case(14)
                do n=jjmin,jjmax
                  do m=iimin,iimax
                    if (az(m,n) .ge. 1 ) then
                      if (counter_ave(m,n) ==0 ) ccb3d_out(k,m,n,0:1)=bio_missing
                    endif
                  enddo
                enddo 
            end select
          endif
        enddo
      case default 
        j=0
        do i=from,to
          j=j+1
          if ( var_ids(i) /= 0.and. var_ave(i) == ave ) then
            k=k+1
            if ( counter_ave(ig,jg).gt.0  .and.  ave== 1) then
            ! add value for average output 
               select case(mode)
               case(1)
                    cc3d_out(k,ig,jg,0:nlev) = cc3d_out(k,ig,jg,0:nlev)+cc(j,0:nlev)
               case(2)
                    cc3d_out(k,ig,jg,0:nlev) = cc3d_out(k,ig,jg,:)+diag(j,0:nlev)
               case(3)
                   call make_flux_output(1,j,nlev, h, c1dim)
                   cc3d_out(k,ig,jg,0:nlev) = cc3d_out(k,ig,jg,0:nlev)+ c1dim(0:nlev)
               case(4)
                    ccb3d_out(k,ig,jg,0:1) = ccb3d_out(k,ig,jg,0:1)+ccb(j,0:1)
               case(5)
                    ccb3d_out(k,ig,jg,0:1) = ccb3d_out(k,ig,jg,0:1)+diagb(j,0:1)
               case(6)
                   call make_flux_output(2,j,nlev, h, c1dim)
                   ccb3d_out(k,ig,jg,0:1) = ccb3d_out(k,ig,jg,0:1)+ c1dim(0:1)
               end select
            else
            ! set value for normal output or reset value for average output on zero 
            ! when counter_ave ==0 . 
               select case(mode)
               case(1)
                    cc3d_out(k,ig,jg,0:nlev) = cc(j,0:nlev)
               case(2)
                    cc3d_out(k,ig,jg,0:nlev) = diag(j,0:nlev)
               case(3)
                   call make_flux_output(1,j,nlev, h, c1dim)
                   cc3d_out(k,ig,jg,0:nlev) = c1dim(0:nlev)
               case(4)
                    ccb3d_out(k,ig,jg,0:1) = ccb(j,0:1)
               case(5)
                    ccb3d_out(k,ig,jg,0:1) = diagb(j,0:1)
               case(6)
                   call make_flux_output(2,j,nlev, h, c1dim)
                   ccb3d_out(k,ig,jg,0:1) = c1dim(0:1)
               end select
            endif
          endif
        enddo
      end select
  end subroutine fill_diag_array
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: set_2d_grid_parameters 
! 
! !INTERFACE:
 subroutine set_2d_grid_prameters(mode,ig,jg)
!
! !IUSES:
 use mem_Param, ONLY:  p_p_ae,p_poro
!
! !INPUT PARAMETERS: 
      implicit none 
      integer,intent(IN)    ::mode
      integer,intent(IN)    ::ig
      integer,intent(IN)    ::jg
!
! !LOCAL VARIABLES 
      integer                ::rc
      integer                ::nc=0
      integer                ::i
      integer                ::j
      integer                ::status
      REALTYPE               :: kd
      REALTYPE               :: rho=2.65
 
!
! !DESCRIPTION: 
!         read porosities from netcdf-file
!         calculatest adsorption ceoof for Phosphate in aerobicx layer from
!         porosties
!
! !REVISION HISTORY: 
!
!  01072006 Created by Piet Ruardij
! 
!EOP
!-------------------------------------------------------------------------
!BOC
 
 if (bio_setup >=2) then
    if ( read_poro ) then
      select case (mode)
        case(0)
            p_poro_default=p_poro(1)
            p_p_ae_default=p_p_ae(1)
            allocate(p_poro_2d(E2DFIELD),stat=rc)   
            if (rc /= 0) stop 'read_2d_grid_parameters: Error allocating memory (p_poro_2d)'
            allocate(p_p_ae_2d(E2DFIELD),stat=rc)   
            p_p_ae_2d=bio_missing
            if (rc /= 0) stop 'read_2d_grid_parameters: Error allocating memory (p_ae_2d)'
 
            call init_2dbio_ncdf( 11,'Input/Ben_Sedprop.nc','Porosity',nc,status,p_poro_2d)
            if ( status.ne.0) call getm_error("set_2d_grid_prameters()",  &
                                  'Could not find name in Input/Ben_Sedprop.nc.' )
            do j=jjmin,jjmax
                do i=iimin,iimax
                    if (az(i,j) .ge. 1 ) then
                       if ( p_poro_2d(i,j) > 0 ) then
                         kd = 4.03087  * (p_poro_2d(i,j) - 0.38662 )/ 0.00415 
                         p_p_ae_2d(i,j) = kd*(1-p_poro_2d(i,j))/p_poro_2d(i,j)*rho
                       else
                         STDERR 'Warning No value for p_poro_2d at',i,j
                         p_p_ae_2d(i,j) = -9999;
                       endif
     
                    endif
                enddo
            enddo
        case(1)
           if ( p_poro_2d(ig,jg) > 0 ) then
             p_p_ae(1)=p_p_ae_2d(ig,jg)
             p_poro(1)=p_poro_2d(ig,jg)
           else
           endif
      end select
   endif
 endif
 
 end subroutine set_2d_grid_prameters 

!EOC
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: init_2d_states
! 
! !INTERFACE:
     subroutine init_2d_grid(file)
!
! !USES:
    use string_functions, only: empty
    use bio_var, only: var_names,stBenStateS,stBenStateE
!
! !INPUT PARAMETERS: 
     character(len=*)              ::file
!
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:  read 2d-field for initializarion of benthic states vars.
!
! !REVISION HISTORY: 
!
!  15072006 Piet Raurdij
! 
!EOP
!-------------------------------------------------------------------------
!BOC
     integer                    ::i
     integer                    ::l
     integer                    ::j
     integer                    ::nc
     integer                    ::status

     if ( empty(file) ) return
     l=1
     j=0
     do i=stBenStateS,stBenStateE
        j=j+1
        call init_2dbio_ncdf(l,file,var_names(i),nc,status,ccb3d(j,:,:,1))
        if ( status ==0 ) then
          if ( l== 1 ) LEVEL1 '---Benthic Initialization from: ',file,':'
           LEVEL2 var_names(i)
        endif
        l=0
     enddo
     call init_2dbio_ncdf(-1,file,'',nc,status,ccb3d(j,:,:,1))
     end subroutine init_2d_grid
!EOC
!-------------------------------------------------------------------------

   end module getm_bio

!-----------------------------------------------------------------------
! Copyright (C) 2004 - Karsten Bolding and Hans Burchard               !
!-----------------------------------------------------------------------
