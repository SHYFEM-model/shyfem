#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: standalone
!
! !INTERFACE:
   module standalone
!
! !DESCRIPTION:
! This module contains all the routines for the standard BFM
! simulation in standalone version, i.e. with a 0.5D setup
! (pelagic and benthic).
! It also contains all the ancillary functions for forcing functions.
!
! !USES:
!  default: all is private.
   use global_mem, only:RLEN
   use constants,  only:SEC_PER_DAY
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public envforcing_bfm,timestepping,init_standalone
   public temperature,salinity,light,lightAtTime,daylength,instLight
   public time_manager,bfm_to_hydro,hydro_to_bfm,rst_bfm_to_hydro
! !PUBLIC DATA MEMBERS:
   ! Note: all read from namelist
   !---------------------------------------------
   ! geographic and dimensional parameters
   !---------------------------------------------
   real(RLEN),public        :: latitude,longitude
   integer, public          :: nboxes
   real(RLEN),public        :: indepth
   !---------------------------------------------
   ! timestepping parameters
   ! real:
   ! mindelt: minimal time step allowed for computation
   ! endtim: endtime of integration
   ! time: actual time
   ! delt: actual time step of global integration
   ! maxdelt: maximal timestep
   !---------------------------------------------
   real(RLEN),public  :: maxdelt,mindelt,endtim, &
                         timesec,delt
   !---------------------------------------------
   ! integer:
   ! nmaxdelt: number of mindelts in maxdelts
   ! nendtim: number of total maxelts to endtim
   ! nmin: actual no. of mindelts in maxdelt intervall
   ! nstep: actual time step in mindelts
   ! ntime: actual time in maxdelts
   ! method: integration method
   !---------------------------------------------
   integer,public     :: nmaxdelt,nendtim,nmin,nstep,ntime, &
                         method
   !---------------------------------------------
   ! forcing function parameters
   !---------------------------------------------
   real(RLEN), public :: tw,ts,tde,sw,ss,lw,ls
   real(RLEN), public :: botdep_c,botdep_n,botdep_p,botdep_si,botox_o
   !---------------------------------------------
   ! arrays for integration routines
   !---------------------------------------------
   real(RLEN),public,dimension(:,:),allocatable :: bbccc3D,bccc3D,ccc_tmp3D
   real(RLEN),public,dimension(:,:),allocatable :: bbccc2D,bccc2D,ccc_tmp2D
   real(RLEN),public                            :: dtm1
   logical,public                               :: sspflag

!
! !PRIVATE DATA MEMBERS:
   real(RLEN),parameter :: PI=3.14159265,RFACTOR=PI/180.
   integer,parameter    :: namlst=10,unit=11

   REALTYPE,dimension(:),allocatable :: cdepth
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi (INGV), Momme Buthenschoen (UNIBO)
!
!  $Log: $
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the standalone BFM
!
! !INTERFACE:
   subroutine init_standalone()
!
! !DESCRIPTION:
!  Main communication of array dimensions between
!  BFM and the standalone setup
!  Read also integration limits
!
!
! !USES:
   use mem, only: NO_D3_BOX_STATES, NO_BOXES,          &
                  NO_BOXES_X, NO_BOXES_Y, NO_BOXES_Z,  &
                  NO_D2_BOX_STATES, NO_BOXES_XY,       &
                  NO_D2_BOX_DIAGNOSS, NO_D3_BOX_DIAGNOSS,&
                  NO_STATES,Depth,Depth_ben, D3STATE, D2STATE
   use global_mem, only:RLEN,LOGUNIT,NML_OPEN,NML_READ,error_msg_prn
   use api_bfm
   use netcdf_bfm, only: init_netcdf_bfm,init_save_bfm
   use time
   
   IMPLICIT NONE
   include '../../../fem3d/param.h'
!
! !INPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
!

!LOCAL VARIABLES:
   namelist /standalone_nml/ nboxes,indepth,maxdelt,    &
            mindelt,endtim,method,latitude,longitude
   namelist /anforcings_nml/ lw,ls,sw,ss,tw,ts,tde,     &
            botdep_c,botdep_n,botdep_p,botdep_si,botox_o
   namelist /time_nml/ timefmt,MaxN,start,stop,simdays
!
! !LOCAL VARIABLES:
!   namelist /standalone_nml/mindelt,method,latitude,longitude
!   namelist /anforcings_nml/tde,botdep_c,botdep_n,botdep_p,botdep_si,botox_o
!   namelist /time_nml/ timefmt,MaxN,start,stop,simdays
!
! !LOCAL VARIABLES:
   real(RLEN) :: tt
   integer    :: dtm1,i

! HYDROCODE DATA
  integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
  common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
  real drr
  common /drr/drr
  real hkv(nkndim)
  common /hkv/hkv

  integer itanf,itend,idt,nits,niter,it
  common /femtim/ itanf,itend,idt,nits,niter,it




!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL2 'init_standalone'
   !---------------------------------------------
   ! Give initial default values
   ! (overwritten with namelist)
   !---------------------------------------------
   nboxes      = 1 ! real data from model
   indepth     = 10.0  !fake data from model
   latitude    = 0.0 !fake from nml
   longitude   = 0.0 !fake from nml
   maxdelt     = 100.0 !fake 
   mindelt     = 1.0 !realdata from str
   endtim      = 360.0 !fake 
   method      = 1 !real from str
  !anforcings
   lw          = 300.0 !fake from model
   ls          = 600.0 !fake "
   sw          = 33.0 !fake "
   ss          = 37.0 !fake "
   tw          = 10.0 !fake "
   ts          = 25.0 !fake "
   tde         = 1.0 !fake "
   botdep_c    = 0.0 !real "
   botdep_n    = 0.0 !real "
   botdep_p    = 0.0 !real "
   botdep_si   = 0.0 !real "
   botox_o     = 0.0 !real  from model
	
   print*,' 2 step '
	
   open(namlst,file='standalone.nml',status='old',action='read',err=100)
   read(namlst,nml=standalone_nml,err=101)
   close(namlst)
   open(namlst,file='standalone.nml',status='old',action='read',err=100)
   read(namlst,nml=anforcings_nml,err=102)
   close(namlst)
   open(namlst,file='standalone.nml',status='old',action='read',err=100)
   read(namlst,nml=time_nml,err=103)
   close(namlst)

! initialize state variables from fem3d
 !  maxdelt =real(idt)
!   nboxes = 1 
!   mindelt = 1
!   print*,' drr',drr
   !---------------------------------------------
   ! set the dimensions
   !---------------------------------------------
   NO_BOXES_X  = nboxes
   NO_BOXES_Y  = 1
   NO_BOXES_Z  = 1
   NO_BOXES    = NO_BOXES_X * NO_BOXES_Y * NO_BOXES_Z
   NO_BOXES_XY = NO_BOXES_X * NO_BOXES_Y
   NO_STATES   = NO_D3_BOX_STATES * NO_BOXES +   &
                 NO_D2_BOX_STATES * NO_BOXES_XY

   LEVEL3 'Number of Boxes:',NO_BOXES
   !---------------------------------------------
   ! initialise the timestepping parameters
   ! Use the GOTM time-manager
   ! Time is given in Julian day and seconds of day
   ! (both integer values)
   !---------------------------------------------
   timestep = maxdelt
   call init_time(MinN,MaxN)
   if (HasRealTime .eqv. .true.) then
      timesec=julianday*SEC_PER_DAY+secondsofday
      simdays=nint(simtime/SEC_PER_DAY)
   else
      timesec=0.0
   end if
   nmaxdelt=1
   LEVEL3 'nmaxdelt: ',nmaxdelt
   tt=maxdelt/2.
   do while (tt.ge.mindelt)
      tt=tt/2.
      nmaxdelt=nmaxdelt*2
      LEVEL3 'nmaxdelt: ',nmaxdelt
   end do
   mindelt=maxdelt/nmaxdelt ! maxdelt = nmaxdelt*mindelt
   nendtim=MaxN
   nstep=nmaxdelt
   ntime=0
   nmin=0
   dtm1=maxdelt
   delt=maxdelt
   if (method.eq.3) delt=2*delt
   LEVEL3 'Integration method: ',method
   LEVEL3 'maxdelt (sec): ',maxdelt
   LEVEL3 'mindelt (sec): ',mindelt
   LEVEL3 'nmaxdelt: ',nmaxdelt
   LEVEL3 'Simulation time (days): ',simdays
   LEVEL3 'nendtim: ',nendtim
   LEVEL3 'Initial time (sec): ',timesec,tt

   !---------------------------------------------
   !---------------------------------------------
   ! initialise the timestepping parameters
   ! Use the GOTM time-manager
   ! Time is given in Julian day and seconds of day
   ! (both integer values)
   !---------------------------------------------
!   timestep = maxdelt
!   timesec=real(it)
!   nmaxdelt=1
!   LEVEL3 'nmaxdelt: ',nmaxdelt
!   tt=maxdelt/2.!

!   do while (tt.ge.mindelt)
!      tt=tt/2.
!      nmaxdelt=nmaxdelt*2
!      LEVEL3 'nmaxdelt: ',nmaxdelt
!   end do
!   mindelt=maxdelt/nmaxdelt ! maxdelt = nmaxdelt*mindelt
!   nstep=nmaxdelt
!   ntime=0
!   nmin=0
!   dtm1=maxdelt
!   delt=maxdelt
!   if (method.eq.3) delt=2*delt
!   LEVEL3 'Integration method: ',method
!   LEVEL3 'maxdelt (sec): ',maxdelt
!   LEVEL3 'mindelt (sec): ',mindelt
!   LEVEL3 'nmaxdelt: ',nmaxdelt
1   LEVEL3 'Initial time (sec): ',timesec

   !---------------------------------------------
   ! Initialise the BFM with standalone settings
   !---------------------------------------------
   call init_bfm(namlst)
   !---------------------------------------------
   ! Initialise state variable names and diagnostics
   !---------------------------------------------
   call set_var_info_bfm
   !---------------------------------------------
   ! Allocate memory and give initial values
   !---------------------------------------------
   ! the argument list is kept for compatibility with GOTM
   call init_var_bfm(namlst,'bfm.nml',unit,bio_setup)
   !---------------------------------------------
   ! Initialize depth 
   !---------------------------------------------
    Depth = indepth
    Depth_ben = Depth
!   LEVEL3 'Box Depth:',Depth,Depth_ben
   !---------------------------------------------
   ! initialise netcdf output
   !---------------------------------------------
   call calcmean_bfm(INIT)
   call init_netcdf_bfm(out_title,'01-01-0000',0,  &
             lat=latitude,lon=longitude,z=Depth,   &
             oceanpoint=(/(i,i=1,NO_BOXES)/),      &
             surfacepoint=(/(i,i=1,NO_BOXES_XY)/), &
             bottompoint=(/(i,i=1,NO_BOXES_XY)/))
   call init_save_bfm
   !---------------------------------------------
   ! allocate and initialise integration arrays
   !---------------------------------------------
   allocate(bbccc3D(NO_D3_BOX_STATES,NO_BOXES))
   allocate(bccc3D(NO_D3_BOX_STATES,NO_BOXES))
   allocate(ccc_tmp3D(NO_D3_BOX_STATES,NO_BOXES))
   allocate(bbccc2D(NO_D2_BOX_STATES,NO_BOXES))
   allocate(bccc2D(NO_D2_BOX_STATES,NO_BOXES))
   allocate(ccc_tmp2D(NO_D2_BOX_STATES,NO_BOXES))
   ! Initialize prior time step for leap-frog:
   if (method == 3) then
      bbccc3d = D3STATE
      bbccc2d = D2STATE
      ccc_tmp3D = D3STATE
      ccc_tmp2D = D3STATE
   end if
   return

100   call error_msg_prn(NML_OPEN,"standalone.f90","standalone.nml")
101   call error_msg_prn(NML_READ,"standalone.f90","standalone_nml")
102   call error_msg_prn(NML_READ,"standalone.f90","anforcings_nml")
103   call error_msg_prn(NML_READ,"standalone.f90","time_nml")

   end subroutine init_standalone
!EOC
!-----------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Light and other environmental forcings used in the BFM
!
! !INTERFACE
   subroutine envforcing_bfm()
!
! !DESCRIPTION
!
! !USES
   use api_bfm
   use global_mem, only: RLEN
   use mem,        only: ETW,ESW,EIR,ESS,SUNQ,ThereIsLight,Wind,&
                         rutQ6c,rutQ6n,rutQ6p,rutQ6s,R6c,R6n,R6p,R6s,O2o,&
			 Depth,Depth_ben
   use mem_Param,  only: LightForcingFlag,p_PAR
   IMPLICIT NONE
   include '../../../fem3d/param.h'
!
! !INPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Momme Butenschoen (UNIBO)
!
! !LOCAL VARIABLES:
   real(RLEN) :: dfrac,wlight,dtime
   integer    :: dyear
   real(RLEN),external :: GetDelta
   real(RLEN) :: biodelta
   real(RLEN) :: wx,wy
   integer    :: l
   real       :: tt,ss,ll

!  ENVFORCING FROM HYDRO

   integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
   common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
   integer itanf,itend,idt,nits,niter,it
   common /femtim/ itanf,itend,idt,nits,niter,it

   real ddepth(nkndim)
   common /ddepth/ddepth
	
   integer node
   common /node/node

!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   LEVEL1 'envforcing_bfm'
   LEVEL2 'time=',timesec
#endif
   timesec =it
   !---------------------------------------------
   ! Computes all the forcings
   !---------------------------------------------
   dtime = timesec/SEC_PER_DAY
   sunq=daylength(dtime,latitude)
   dfrac=(dtime-floor(dtime)) ! fraction of the day
   dyear=mod(dtime,360.) ! Day of the year
   wlight=light(dyear,dfrac)
   select case(LightForcingFlag)
    case (3) ! light on/off distribution for daylight average
      ThereIsLight=lightAtTime(dfrac,sunq)
      wlight=wlight*ThereIsLight
    case (1) ! instantaneous light distribution
      wlight=instLight(wlight,sunq,dfrac)
    case default ! light constant during the day
   end select
    ESS = 0.
!    ETW = temperature(dyear,dfrac)
!    ESW = salinity(dyear,dfrac)
!    EIR = wlight*p_PAR

#ifdef DEBUG
   LEVEL2 'ETW=',ETW
   LEVEL2 'ESW=',ESW
   LEVEL2 'EIR=',EIR
#endif
	
! feeding envforcing vector from HYDROcode

     call get_wind(node,wx,wy)
     Wind=sqrt(wx**2+wy**2)

     call get_light(node,ll)
     EIR = ll
	
     l = 1
     call getts(l,node,tt,ss)
     ETW = tt
     ESW = ss

!    print*, ETW,ESW,tempv(1,2),saltv(1,2)
!    stop
!    ETW=20
!    ESW=35
!    EIR =0
!    EIR=wlight*p_PAR   ! to be changed 
!   ESS=0.
	if(ddepth(node).le.20)then  !cucco mod
	 ddepth(node)=20
	endif
   Depth=ddepth(node)
   Depth_ben=Depth
!   print*,node,ETW,ESW,Wind,Depth
   call CalcVerticalExtinction

   if (bio_setup==2) then
      ! Bottom deposition and ventilation fluxes
      ! (mg C m^-2 d^-1 or mmol NUT m^-2 d^-1)
      ! currently constant deposition rates read from namelist
      ! (se to zero for no deposition)
      biodelta=GetDelta()
      R6c(:) = R6c(:)+botdep_c*biodelta
      R6n(:) = R6n(:)+botdep_n*biodelta
      R6p(:) = R6p(:)+botdep_p*biodelta
      R6s(:) = R6s(:)+botdep_si*biodelta
      O2o(:) = O2o(:)+botox_o*biodelta
   end if

   end subroutine envforcing_bfm
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:
!
! !INTERFACE:
   FUNCTION daylength(time,latitude)
!
! !DESCRIPTION:
! This function computes the length of the daylight period in hours
! as a function of time of the year (days) and latitude
!
! !USES:
   use global_mem, only:RLEN
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   real(RLEN),intent(in) :: time
   real(RLEN),intent(in) :: latitude
!
! !INPUT/OUTPUT PARAMETERS:

!
! !OUTPUT PARAMETERS:
   real(RLEN) :: daylength
!
! !REVISION HISTORY:
!  Author(s): Momme Butenschoen (UNIBO)
!
! !LOCAL VARIABLES:
   real(RLEN)           :: declination
   real(RLEN),parameter :: cycle=360.
!
!EOP
!-----------------------------------------------------------------------
!BOC
   declination = -0.406*cos(2.*PI*int(time)/cycle)
   daylength = acos(-tan(declination)*tan(latitude*RFACTOR))/PI*24.
   return

   END FUNCTION daylength
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:
!
! !INTERFACE:
   FUNCTION lightAtTime(df,dl)
!
! !DESCRIPTION:
!  This function determines whether there is light at a certain time
!  of the day. Returns an integer value 0 or 1
!
! !USES:
   use global_mem, only:RLEN
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
 real(RLEN),intent(in) :: df,dl

!
! !INPUT/OUTPUT PARAMETERS:

!
! !OUTPUT PARAMETERS:
 integer :: lightAtTime

!
! !REVISION HISTORY:
!  Author(s): Momme Butenschoen (UNIBO)
!
! !LOCAL VARIABLES:
 real(RLEN) :: daytime,daylength
!
!EOP
!-----------------------------------------------------------------------
!BOC
    daytime=df*24. ! time of the day = fraction of the day * 24
    daytime=abs(daytime-12.) ! distance from noon
    daylength=dl/2.
    if(daytime.lt.daylength) then
      lightAtTime=1
    else
      lightAtTime=0
    endif
    return
   END FUNCTION lightAtTime
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:
!
! !INTERFACE:
   FUNCTION instLight(l,dl,df)
!
! !DESCRIPTION:
!  This function computes the instantaneous light at a certain time of
!  the day
!
! !USES:
   use global_mem, only:RLEN
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   real(RLEN),intent(in) :: df,dl,l
!
! !INPUT/OUTPUT PARAMETERS:

!
! !OUTPUT PARAMETERS:
   real(RLEN) :: instLight
!
! !REVISION HISTORY:
!  Author(s): Momme Butenschoen (UNIBO)
!
! !LOCAL VARIABLES:
   real(RLEN) :: daylength,daytime
!
!EOP
!-----------------------------------------------------------------------
!BOC
     daytime=df*24. ! time of the day = fraction of the day * 24
     daytime=abs(daytime-12.) ! distance from noon
     daylength=dl/2.
     if(daytime.lt.daylength) then
       daytime=daytime/daylength*PI
       instLight=l*cos(daytime)+l
     else
       instLight=0.
     endif
     return
   END FUNCTION instLight
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:
!
! !INTERFACE:
   FUNCTION salinity(dy,df)
!
! !DESCRIPTION:
!  This function provides an articial salinity value given the
!  parameters in the standalone.nml namelist
!
! !USES:
   use global_mem, only:RLEN
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,intent(in)    :: dy
   real(RLEN),intent(in) :: df
!
!
! !INPUT/OUTPUT PARAMETERS:

!
! !OUTPUT PARAMETERS:
   real(RLEN) :: salinity

!
! !REVISION HISTORY:
!  Author(s): Momme Butenschoen (UNIBO)
!
! !LOCAL VARIABLES:

!
!EOP
!-----------------------------------------------------------------------
!BOC
     salinity=(ss+sw)/2.-(ss-sw)/2.*cos((dy+(df-.5))*RFACTOR)
   END FUNCTION salinity
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:
!
! !INTERFACE:
   FUNCTION temperature(dy,df)
!
! !DESCRIPTION:
!  This function provides an articial temperature value given the
!  parameters in the standalone.nml namelist
!
! !USES:
   use global_mem, only:RLEN
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,intent(in)    :: dy
   real(RLEN),intent(in) :: df
!
!
! !INPUT/OUTPUT PARAMETERS:

!
! !OUTPUT PARAMETERS:
   real(RLEN) :: temperature

!
! !REVISION HISTORY:
!  Author(s): Momme Butenschoen (UNIBO)
!
! !LOCAL VARIABLES:

!
!EOP
!-----------------------------------------------------------------------
!BOC
     temperature=(ts+tw)/2.-(ts-tw)/2.*cos((dy+(df-.5))*RFACTOR) &
                    -tde*.5*cos(2*Pi*df)
   END FUNCTION temperature
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:
!
! !INTERFACE:
   FUNCTION light(dy,df)
!
! !DESCRIPTION:
!  This function provides an articial light value given the
!  parameters in the standalone.nml namelist
!
! !USES:
   use global_mem, only:RLEN
   IMPLICIT NONE
! !INPUT PARAMETERS:
   integer,intent(in)    :: dy
   real(RLEN),intent(in) :: df
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   real(RLEN) :: light
! !REVISION HISTORY:
!  Author(s): Momme Butenschoen (UNIBO)
!
! !LOCAL VARIABLES:
!
!EOP
!-----------------------------------------------------------------------
!BOC
     light=(ls+lw)/2.-(ls-lw)/2.*cos(dy*RFACTOR)
   END FUNCTION light
!EOC


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:
!
! !INTERFACE:
   SUBROUTINE timestepping()
!
! !DESCRIPTION:
!
!
! !USES:
   use global_mem, only:RLEN
   use netcdf_bfm, only: save_bfm
   use mem
   use api_bfm, only: out_delta
   IMPLICIT NONE

   integer itanf,itend,idt,nits,niter,it
   common /femtim/ itanf,itend,idt,nits,niter,it

! !INPUT PARAMETERS:
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
! !REVISION HISTORY:
!  Author(s): Momme Butenschoen (UNIBO)
!
! !LOCAL VARIABLES:

!
!EOP
!-----------------------------------------------------------------------
!BOC

!   LEVEL1 'timesteppingntime='
      ntime=it/idt

!	   print*,ntime,it,idt
      call envforcing_bfm
      call EcologyDynamics
      select case (method)
         case (2)
            call integrationRK2
         case (3)
            call integrationLf
         case default
            call integrationEfw
      end select
      call calcmean_bfm(ACCUMULATE)
!      if (mod(ntime,out_delta).eq.0) then
!         LEVEL1 'OUTPUT' , timesec/SEC_PER_DAY
!         call calcmean_bfm(MEAN)
!         call save_bfm(timesec)
!      end if
      call ResetFluxes

   END SUBROUTINE timestepping
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Providing info on variables
!
! !INTERFACE:

   subroutine time_manager
  
   !---------------------------------------------
   ! manage the timestepping parameters
   !---------------------------------------------
  
 
   use mem
   use time
   use global_mem, only:RLEN
	
   IMPLICIT NONE
   
   integer itanf,itend,idt,nits,niter,it
   common /femtim/ itanf,itend,idt,nits,niter,it

   real drr
   common /drr/drr

   real(RLEN) :: tt
   integer    :: dtm1,i

   mindelt =1 
   maxdelt = real(idt)
   timestep = maxdelt
   timesec=real(it)-maxdelt
   nmaxdelt=1
!   LEVEL3 'nmaxdelt: ',nmaxdelt
   tt=maxdelt/2.

   do while (tt.ge.mindelt)
      tt=tt/2.
      nmaxdelt=nmaxdelt*2
!      LEVEL3 'nmaxdelt: ',nmaxdelt,tt
   end do
   mindelt=maxdelt/nmaxdelt ! maxdelt = nmaxdelt*mindelt
   nstep=nmaxdelt
   ntime=0
   nmin=0
   dtm1=maxdelt
   delt=maxdelt
!   LEVEL3 mindelt,nstep,dtm1,delt,timesec
   if (method.eq.3) delt=2*delt

!   maxdelt =drr
!   timestep = maxdelt
!   nmaxdelt=1
!   tt=maxdelt/2.

!    print*,tt ,mindelt
!   do while (tt.ge.mindelt)
!      tt=tt/2.
!      nmaxdelt=nmaxdelt*2
!       LEVEL3 'nmaxdelt: ',nmaxdelt,maxdelt
!	print*,'ok3'
!   end do
!   mindelt=maxdelt/nmaxdelt ! maxdelt = nmaxdelt*mindelt
!   nstep=nmaxdelt
!   nmin=0
!   dtm1=maxdelt
!   delt=maxdelt

!   timesec = it
   

   end subroutine time_manager

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Providing info on variables
!
! !INTERFACE:

   	subroutine bfm_to_hydro(k,b1cn,nbfmv1,b2cn,nbfmv2,b2cn_a,b2cn_b,b2cn_c,b2cn_d,b3cn,nbfmv3,b3cn_a,b3cn_b,b3cn_c)


   use mem
   implicit  none

   include '../../../fem3d/param.h'
   include '../../../fem3d/bfm_common.h'
   integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
   integer nlvdi,nlv
   common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
   common /level/ nlvdi,nlv


   integer l,k

	integer nbfmv1,nbfmv2,nbfmv3
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
   
                   fO2o(k)= 0.+ O2o(1)
	           fN1p(k)= 0.+ N1p(1)
                   fN3n(k)= 0.+ N3n(1)
                   fN4n(k)= 0.+ N4n(1)
                   fO4n(k)= 0.+ O4n(1)
                   fN5s(k)= 0.+ N5s(1)
                   fN6r(k)= 0.+ N6r(1)
                   fB1c(k)= 0.+ B1c(1)
                   fB1n(k)= 0.+ B1n(1)
                   fB1p(k)= 0.+ B1p(1)
                   fP1c(k)= 0.+ P1c(1)
                   fP1n(k)= 0.+ P1n(1)
                   fP1p(k)= 0.+ P1p(1)
                   fP1l(k)= 0.+ P1l(1)
                   fP1s(k)= 0.+ P1s(1)
                   fP2c(k)= 0.+ P2c(1)
                   fP2n(k)= 0.+ P2n(1)
                   fP2p(k)= 0.+ P2p(1)
                   fP2l(k)= 0.+ P2l(1)
                   fP3c(k)= 0.+ P3c(1)
                   fP3n(k)= 0.+ P3n(1)
                   fP3p(k)= 0.+ P3p(1)
                   fP3l(k)= 0.+ P3l(1)
                   fP4c(k)= 0.+ P4c(1)
                   fP4n(k)= 0.+ P4n(1)
                   fP4p(k)= 0.+ P4p(1)
                   fP4l(k)= 0.+ P4l(1)
                   fZ3c(k)= 0.+ Z3c(1)
                   fZ3n(k)= 0.+ Z3n(1)
                   fZ3p(k)= 0.+ Z3p(1)
                   fZ4c(k)= 0.+ Z4c(1)
                   fZ4n(k)= 0.+ Z4n(1)
                   fZ4p(k)= 0.+ Z4p(1)
                   fZ5c(k)= 0.+ Z5c(1)
                   fZ5n(k)= 0.+ Z5n(1)
                   fZ5p(k)= 0.+ Z5p(1)
                   fZ6c(k)= 0.+ Z6c(1)
                   fZ6n(k)= 0.+ Z6n(1)
                   fZ6p(k)= 0.+ Z6p(1)
                   fR1c(k)= 0.+ R1c(1)
                   fR1n(k)= 0.+ R1n(1)
                   fR1p(k)= 0.+ R1p(1)
                   fR2c(k)= 0.+ R2c(1)
                   fR6c(k)= 0.+ R6c(1)
                   fR6n(k)= 0.+ R6n(1)
                   fR6p(k)= 0.+ R6p(1)
                   fR6s(k)= 0.+ R6s(1)
                   fR7c(k)= 0.+ R7c(1)   
         do l=1,nlv
	          b1cn(l,k,1) = 0.+ O2o(1)
	          b1cn(l,k,2) = 0.+ N1p(1) 
                  b1cn(l,k,3) = 0.+ N3n(1) 
                  b1cn(l,k,4) = 0.+ N4n(1) 
                  b1cn(l,k,5) = 0.+ O4n(1) 
                  b1cn(l,k,6) = 0.+ N5s(1) 
                  b1cn(l,k,7) = 0.+ N6r(1) 
                  b2cn(l,k,1) = 0.+ B1c(1) 
                  b2cn(l,k,2) = 0.+ P1c(1) 
                  b2cn(l,k,3)= 0.+ P2c(1) 
                  b2cn(l,k,4)= 0.+ P3c(1) 
                  b2cn(l,k,5)= 0.+ P4c(1) 
                  b2cn(l,k,6)= 0.+ Z3c(1) 
                  b2cn(l,k,7)= 0.+ Z4c(1) 
                  b2cn(l,k,8)= 0.+ Z5c(1) 
                  b2cn(l,k,9)= 0.+ Z6c(1) 
                  b3cn(l,k,1)= 0.+ R1c(1) 
                  b3cn(l,k,2)= 0.+ R2c(1) 
                  b3cn(l,k,3)= 0.+ R6c(1) 
                  b3cn(l,k,4)= 0.+ R7c(1) 
 
        	b2cn_a(l,k,1) = 0. + B1n(1)
	        b2cn_b(l,k,1) = 0. + B1p(1)
                b2cn_a(l,k,2) = 0. + P1n(1)
                b2cn_b(l,k,2) = 0. + P1p(1)
                b2cn_c(l,k,2) = 0. + P1l(1)
                b2cn_d(l,k,2) = 0. + P1s(1)
                b2cn_a(l,k,3) = 0. + P2n(1)
                b2cn_b(l,k,3) = 0. + P2p(1)
                b2cn_c(l,k,3) = 0. + P2l(1)
                b2cn_a(l,k,4) = 0. + P3n(1)  
                b2cn_b(l,k,4) = 0. + P3p(1)
                b2cn_c(l,k,4) = 0. + P3l(1)
                b2cn_a(l,k,5) = 0. + P4n(1)
                b2cn_b(l,k,5) = 0. + P4p(1)
                b2cn_c(l,k,5) = 0. + P4l(1)
                b2cn_a(l,k,6) =  0. + Z3n(1)
                b2cn_b(l,k,6) = 0. + Z3p(1)
                b2cn_a(l,k,7) =  0. + Z4n(1)
                b2cn_b(l,k,7) = 0. + Z4p(1)
                b2cn_a(l,k,8) = 0. + Z5n(1)
                b2cn_b(l,k,8) = 0. + Z5p(1)
                b2cn_a(l,k,9) = 0. + Z6n(1) 
                b2cn_b(l,k,9) = 0. + Z6p(1)
 

                b3cn_a(l,k,1) = 0. + R1n(1)
                b3cn_b(l,k,1) = 0. + R1p(1)
                b3cn_a(l,k,3) = 0. + R6n(1)
                b3cn_b(l,k,3) = 0. + R6p(1)
                b3cn_c(l,k,3) = 0. + R6s(1) 
 
                                                     
       end do                        
 
                                      
                                      
                                      
   end subroutine bfm_to_hydro        
!EOC                                  
!-----------------------------------------------------------------------
!BOP                                  
!                                     
! !IROUTINE: Providing info on variables
!                                     
! !INTERFACE:                         

   	subroutine hydro_to_bfm(k,b1cn,nbfmv1,b2cn,nbfmv2,b2cn_a,b2cn_b,b2cn_c,b2cn_d,b3cn,nbfmv3,b3cn_a,b3cn_b,b3cn_c)

   use global_mem
   use mem

   implicit  none

   include '../../../fem3d/param.h'
   
   integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
   integer nlvdi,nlv
   common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
   common /level/ nlvdi,nlv

        integer nbfmv1,nbfmv2,nbfmv3
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


   integer l,k

   
     do l=1,nlv
   	  O2o(1) =b1cn(l,k,1)
          N1p(1) =b1cn(l,k,2)
          N3n(1) =b1cn(l,k,3)
          N4n(1) =b1cn(l,k,4)
          O4n(1) =b1cn(l,k,5)
          N5s(1) =b1cn(l,k,6)
          N6r(1) =b1cn(l,k,7)
          B1c(1) =b2cn(l,k,1)
          P1c(1) =b2cn(l,k,2)
          P2c(1) =b2cn(l,k,3)
          P3c(1) =b2cn(l,k,4)
          P4c(1) =b2cn(l,k,5)
          Z3c(1) =b2cn(l,k,6)
          Z4c(1) =b2cn(l,k,7)
          Z5c(1) =b2cn(l,k,8)
          Z6c(1) =b2cn(l,k,9)
          R1c(1) =b3cn(l,k,1)
          R2c(1) =b3cn(l,k,2)
          R6c(1) =b3cn(l,k,3)
          R7c(1) =b3cn(l,k,4)

       B1n(1)= b2cn_a(l,k,1)
       B1p(1)= b2cn_b(l,k,1)
       P1n(1)= b2cn_a(l,k,2)
       P1p(1)= b2cn_b(l,k,2)
       P1l(1)= b2cn_c(l,k,2)
       P1s(1)= b2cn_d(l,k,2)
       P2n(1)= b2cn_a(l,k,3)
       P2p(1)= b2cn_b(l,k,3)
       P2l(1)= b2cn_c(l,k,3)
       P3n(1)= b2cn_a(l,k,4)
       P3p(1)= b2cn_b(l,k,4) 
       P3l(1)= b2cn_c(l,k,4) 
        P4n(1)= b2cn_a(l,k,5) 
       P4p(1)= b2cn_b(l,k,5) 
       P4l(1)= b2cn_c(l,k,5) 
       Z3n(1)= b2cn_a(l,k,6) 
       Z3p(1)= b2cn_b(l,k,6) 
       Z4n(1)= b2cn_a(l,k,7) 
       Z4p(1)= b2cn_b(l,k,7) 
       Z5n(1)= b2cn_a(l,k,8) 
       Z5p(1)= b2cn_b(l,k,8) 
       Z6n(1)= b2cn_a(l,k,9) 
       Z6p(1)= b2cn_b(l,k,9) 

       R1n(1)= b3cn_a(l,k,1) 
       R1p(1)= b3cn_b(l,k,1) 
       R6n(1)= b3cn_a(l,k,3) 
       R6p(1)= b3cn_b(l,k,3) 
       R6s(1)= b3cn_c(l,k,3) 
 
 
        end do
 
 
   end subroutine hydro_to_bfm
!EOC                                  
!-----------------------------------------------------------------------
!BOP                                  
!                                     
! !IROUTINE: Providing info on variables
!                                     
! !INTERFACE:   


   	subroutine rst_bfm_to_hydro(k,b1cn,nbfmv1,b2cn,nbfmv2,b2cn_a,b2cn_b,b2cn_c,b2cn_d,b3cn,nbfmv3,b3cn_a,b3cn_b,b3cn_c)


   use mem
   implicit  none

   include '../../../fem3d/param.h'
   include '../../../fem3d/bfm_common.h'
   integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
   integer nlvdi,nlv
   common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
   common /level/ nlvdi,nlv


   integer l,k

	integer nbfmv1,nbfmv2,nbfmv3
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
   
         do l=1,nlv
	          b1cn(l,k,1) = 0.+ fO2o(k)
	          b1cn(l,k,2) = 0.+ fN1p(k) 
                  b1cn(l,k,3) = 0.+ fN3n(k) 
                  b1cn(l,k,4) = 0.+ fN4n(k) 
                  b1cn(l,k,5) = 0.+ fO4n(k) 
                  b1cn(l,k,6) = 0.+ fN5s(k) 
                  b1cn(l,k,7) = 0.+ fN6r(k) 
                  b2cn(l,k,1) = 0.+ fB1c(k) 
                  b2cn(l,k,2) = 0.+ fP1c(k) 
                  b2cn(l,k,3)= 0.+ fP2c(k) 
                  b2cn(l,k,4)= 0.+ fP3c(k) 
                  b2cn(l,k,5)= 0.+ fP4c(k) 
                  b2cn(l,k,6)= 0.+ fZ3c(k) 
                  b2cn(l,k,7)= 0.+ fZ4c(k) 
                  b2cn(l,k,8)= 0.+ fZ5c(k) 
                  b2cn(l,k,9)= 0.+ fZ6c(k) 
                  b3cn(l,k,1)= 0.+ fR1c(k) 
                  b3cn(l,k,2)= 0.+ fR2c(k) 
                  b3cn(l,k,3)= 0.+ fR6c(k) 
                  b3cn(l,k,4)= 0.+ fR7c(k) 
 
        	b2cn_a(l,k,1) = 0. + fB1n(k)
	        b2cn_b(l,k,1) = 0. + fB1p(k)
                b2cn_a(l,k,2) = 0. + fP1n(k)
                b2cn_b(l,k,2) = 0. + fP1p(k)
                b2cn_c(l,k,2) = 0. + fP1l(k)
                b2cn_d(l,k,2) = 0. + fP1s(k)
                b2cn_a(l,k,3) = 0. + fP2n(k)
                b2cn_b(l,k,3) = 0. + fP2p(k)
                b2cn_c(l,k,3) = 0. + fP2l(k)
                b2cn_a(l,k,4) = 0. + fP3n(k)  
                b2cn_b(l,k,4) = 0. + fP3p(k)
                b2cn_c(l,k,4) = 0. + fP3l(k)
                b2cn_a(l,k,5) = 0. + fP4n(k)
                b2cn_b(l,k,5) = 0. + fP4p(k)
                b2cn_c(l,k,5) = 0. + fP4l(k)
                b2cn_a(l,k,6) =  0. + fZ3n(k)
                b2cn_b(l,k,6) = 0. + fZ3p(k)
                b2cn_a(l,k,7) =  0. + fZ4n(k)
                b2cn_b(l,k,7) = 0. + fZ4p(k)
                b2cn_a(l,k,8) = 0. + fZ5n(k)
                b2cn_b(l,k,8) = 0. + fZ5p(k)
                b2cn_a(l,k,9) = 0. + fZ6n(k) 
                b2cn_b(l,k,9) = 0. + fZ6p(k)
 

                b3cn_a(l,k,1) = 0. + fR1n(k)
                b3cn_b(l,k,1) = 0. + fR1p(k)
                b3cn_a(l,k,3) = 0. + fR6n(k)
                b3cn_b(l,k,3) = 0. + fR6p(k)
                b3cn_c(l,k,3) = 0. + fR6s(k) 
 
                                                     
       end do                        
 
                                      
                                      
                                      
   end subroutine rst_bfm_to_hydro        
!EOC                                  
!-----------------------------------------------------------------------
!BOP                                  
!                                     
! !IROUTINE: Providing info on variables
!                                     
! !INTERFACE:                         
!
   END MODULE standalone
 
                                   
