!$Id: ode_solvers.F90,v 1.7 2005/11/18 10:59:35 kbk Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: General ODE solver
!
! !INTERFACE:
!BEGIN_BFM:MAV MODIFY
! addition of bio_setup and number of benthic variables
   subroutine ode_solver(solver,numc,nlev,dt,h,cc,t,bio_setup,numbc)
!END_BFM
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: solver,nlev,numc
   integer, intent(in)                 :: bio_setup,numbc
   REALTYPE, intent(in)                :: dt
   REALTYPE, intent(in)                :: t(0:nlev)
   REALTYPE, intent(in)                :: h(0:nlev)
!
! !INPUT/OUTPUT PARAMETER:
   REALTYPE, intent(inout)             :: cc(1:numc,0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC
   select case (solver)
      case (1)
         call euler_forward(dt,numc,nlev,h,t,bio_setup,numbc)
      case (2)
         call runge_kutta_2(dt,numc,nlev,h,t,bio_setup,numbc)
      case (3)
         call runge_kutta_4(dt,numc,nlev,h,t,bio_setup,numbc)
      case (4)
         call patankar(dt,numc,nlev,h,t,bio_setup,numbc)
      case (5)
         call patankar_runge_kutta_2(dt,numc,nlev,h,t,bio_setup,numbc)
      case (6)
         call patankar_runge_kutta_4(dt,numc,nlev,h,t,bio_setup,numbc)
      case (7)
         call modified_patankar(dt,numc,nlev,h,t,bio_setup,numbc)
      case (8)
         call modified_patankar_2(dt,numc,nlev,h,t,bio_setup,numbc)
      case (9)
         call modified_patankar_4(dt,numc,nlev,h,t,bio_setup,numbc)
      case (10)
         call emp_1(dt,numc,nlev,h,t,bio_setup,numbc)
      case (11)
         call emp_2(dt,numc,nlev,h,t,bio_setup,numbc)
      case default
         stop "bio: no valid solver method specified in bio.inp !"
   end select

   return
   end subroutine ode_solver
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Euler-forward scheme for geobiochemical models
!
! !INTERFACE:
   subroutine euler_forward(dt,numc,nlev,h,t,bio_setup,numbc)
!
! !DESCRIPTION:
!
! !USES
  use bio_var, only: cc,pp,dd,ccb,ppb,ddb
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
   integer, intent(in)                 :: numc,nlev
   integer, intent(in)                 :: numbc,bio_setup
   REALTYPE, intent(in)                :: h(0:nlev)
   REALTYPE, intent(in)                :: t(0:nlev)
!
! !INPUT/OUTPUT PARAMETER:

!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
  logical  :: first
  integer  :: i
!EOP
!-----------------------------------------------------------------------
!BOC
   first=.true.

   call process_model(first,numc,nlev,h,t)

   ! integrate pelagic variables 
   if (bio_setup/=2) then
     cc(:,:) = cc(:,:) + dt*sum(pp(:,:,:)-dd(:,:,:),2)
     do i=1,numc
       pp(i,i,1:nlev) = _ZERO_
       dd(i,i,1:nlev) = _ZERO_
     end do
   end if
   ! integrate benthic variables 
    if (bio_setup>1) then
      ccb(:,:) = ccb(:,:) + dt*sum(ppb(:,:,:)-ddb(:,:,:),2)
      ! reset diagonal terms only
      do i=1,numbc
        ppb(i,i,:) = _ZERO_
        ddb(i,i,:) = _ZERO_
      end do
    end if

   return
   end subroutine euler_forward
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Second-order Runge-Kutta scheme for geobiochemical models
!
! !INTERFACE:
   subroutine runge_kutta_2(dt,numc,nlev,h,t,bio_setup,numbc)
!
! !DESCRIPTION:
!
! !USES:
  use bio_var, only: cc,pp,dd,ccb,ppb,ddb
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
   integer, intent(in)                 :: numc,nlev
   integer, intent(in)                 :: numbc,bio_setup
   REALTYPE, intent(in)                :: h(0:nlev)
   REALTYPE, intent(in)                :: t(0:nlev)
!
! !INPUT/OUTPUT PARAMETER:

!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
  logical  :: first
  REALTYPE :: rhs(1:numc,0:nlev),rhs1(1:numc)
  REALTYPE :: cc0(1:numc,0:nlev)
  REALTYPE :: ccb0(1:numc,0:nlev)
  integer  :: i,j,ci
!EOP
!-----------------------------------------------------------------------
!BOC

   first=.true.

   if (bio_setup/=2) cc0 = cc
   if (bio_setup>1) ccb0 = ccb

   call process_model(first,numc,nlev,h,t)

   if (bio_setup/=2) then
     ! store the initial pelagic state
     cc0 = cc
     do ci=1,nlev
        do i=1,numc
           rhs(i,ci)=0.
           do j=1,numc
              rhs(i,ci)=rhs(i,ci)+pp(i,j,ci)-dd(i,j,ci)
           end do
           cc(i,ci)=cc0(i,ci)+dt*rhs(i,ci)
        end do
     end do
     do i=1,numc
       pp(i,i,1:nlev) = _ZERO_
       dd(i,i,1:nlev) = _ZERO_
     end do
   end if ! bio_setup/=2

   if (bio_setup>1) then
     ! store the initial benthic state
     ccb0 = ccb
     do ci=1,nlev
        do i=1,numc
           rhs(i,ci)=0.
           do j=1,numc
              rhs(i,ci)=rhs(i,ci)+pp(i,j,ci)-dd(i,j,ci)
           end do
           ccb(i,ci)=ccb0(i,ci)+dt*rhs(i,ci)
        end do
     end do
     do i=1,numc
       ppb(i,i,1:nlev) = _ZERO_
       ddb(i,i,1:nlev) = _ZERO_
     end do
   end if ! bio_setup>1
 
   
   call process_model(first,numc,nlev,h,t)

   if (bio_setup/=2) then
     do ci=1,nlev
        do i=1,numc
           rhs1(i)=0.
           do j=1,numc
              rhs1(i)=rhs1(i)+pp(i,j,ci)-dd(i,j,ci)
           end do
           cc(i,ci)=cc0(i,ci)+dt*0.5*(rhs(i,ci)+rhs1(i))
        end do
     end do
     do i=1,numc
       pp(i,i,1:nlev) = _ZERO_
       dd(i,i,1:nlev) = _ZERO_
     end do
   end if ! bio_setup/=2

   if (bio_setup>1) then
     do ci=1,nlev
        do i=1,numc
           rhs1(i)=0.
           do j=1,numc
              rhs1(i)=rhs1(i)+ppb(i,j,ci)-ddb(i,j,ci)
           end do
           ccb(i,ci)=ccb0(i,ci)+dt*0.5*(rhs(i,ci)+rhs1(i))
        end do
     end do
     do i=1,numc
       ppb(i,i,1:nlev) = _ZERO_
       ddb(i,i,1:nlev) = _ZERO_
     end do
   end if ! bio_setup>1


   return
   end subroutine runge_kutta_2
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Fourth-order Runge-Kutta scheme for geobiochemical models
!
! !INTERFACE:
   subroutine runge_kutta_4(dt,numc,nlev,h,t,bio_setup,numbc)
!
! !DESCRIPTION:
!
! !USES:
  use bio_var, only: cc,pp,dd,ccb,ppb,ddb
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
   integer, intent(in)                 :: numc,nlev
   integer, intent(in)                 :: numbc,bio_setup
   REALTYPE, intent(in)                :: h(0:nlev)
   REALTYPE, intent(in)                :: t(0:nlev)
!
! !INPUT/OUTPUT PARAMETER:
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
  logical  :: first
  REALTYPE :: rhs(1:numc,0:nlev),rhs1(1:numc,0:nlev)
  REALTYPE :: rhs2(1:numc,0:nlev),rhs3(1:numc,0:nlev)
  REALTYPE :: cc0(1:numc,0:nlev)
  integer  :: i,j,ci
!EOP
!-----------------------------------------------------------------------
!BOC
   ! store the initial state
   cc0 = cc

   first=.true.
   call process_model(first,numc,nlev,h,t)

   do ci=1,nlev
      do i=1,numc
         rhs(i,ci)=0.
         do j=1,numc
            rhs(i,ci)=rhs(i,ci)+pp(i,j,ci)-dd(i,j,ci)
         end do
         cc(i,ci)=cc0(i,ci)+dt*rhs(i,ci)
      end do
   end do

   do i=1,numc
     pp(i,i,1:nlev) = _ZERO_
     dd(i,i,1:nlev) = _ZERO_
   end do

   call process_model(first,numc,nlev,h,t)

   do ci=1,nlev
      do i=1,numc
         rhs1(i,ci)=0.
         do j=1,numc
            rhs1(i,ci)=rhs1(i,ci)+pp(i,j,ci)-dd(i,j,ci)
         end do
         cc(i,ci)=cc0(i,ci)+dt*rhs1(i,ci)
      end do
   end do

   do i=1,numc
     pp(i,i,1:nlev) = _ZERO_
     dd(i,i,1:nlev) = _ZERO_
   end do

   call process_model(first,numc,nlev,h,t)

   do ci=1,nlev
      do i=1,numc
         rhs2(i,ci)=0.
         do j=1,numc
            rhs2(i,ci)=rhs2(i,ci)+pp(i,j,ci)-dd(i,j,ci)
         end do
         cc(i,ci)=cc0(i,ci)+dt*rhs2(i,ci)
      end do
   end do

   do i=1,numc
     pp(i,i,1:nlev) = _ZERO_
     dd(i,i,1:nlev) = _ZERO_
   end do

   call process_model(first,numc,nlev,h,t)

   do ci=1,nlev
      do i=1,numc
         rhs3(i,ci)=0.
         do j=1,numc
            rhs3(i,ci)=rhs3(i,ci)+pp(i,j,ci)-dd(i,j,ci)
         end do
         cc(i,ci)=cc0(i,ci)+dt*1./3.*(0.5*rhs(i,ci)+rhs1(i,ci)+rhs2(i,ci)+0.5*rhs3(i,ci))
      end do
   end do

   do i=1,numc
     pp(i,i,1:nlev) = _ZERO_
     dd(i,i,1:nlev) = _ZERO_
   end do

   return
   end subroutine runge_kutta_4
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Patankar scheme for geobiochemical models
!
! !INTERFACE:
   subroutine patankar(dt,numc,nlev,h,t,bio_setup,numbc)
!
! !DESCRIPTION:
!
! !USES
  use bio_var, only: cc,pp,dd,ccb,ppb,ddb
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
   integer, intent(in)                 :: numc,nlev
   integer, intent(in)                 :: numbc,bio_setup
   REALTYPE, intent(in)                :: t(0:nlev)
   REALTYPE, intent(in)                :: h(0:nlev)
!
! !INPUT/OUTPUT PARAMETER:
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
  logical  :: first
  REALTYPE :: ppsum(numc,0:nlev),ddsum(numc,0:nlev)
  integer  :: i,j,ci
!EOP
!-----------------------------------------------------------------------
!BOC
   first=.true.
   call process_model(first,numc,nlev,h,t)

   if (bio_setup/=2) then
! scalar
!     do ci=1,nlev
!        do i=1,numc
!           ppsum=0.
!           ddsum=0.
!           do j=1,numc
!              ppsum=ppsum+pp(i,j,ci)
!              ddsum=ddsum+dd(i,j,ci)
!           end do
!           cc(i,ci)=(cc(i,ci)+dt*ppsum)/(1.+dt*ddsum/cc(i,ci))
!        end do
!     end do
     cc = (cc+dt*sum(pp(:,:,:),2))/(1.+dt*sum(dd(:,:,:),2)/cc)
     ! reset diagonal terms only
     do i=1,numc
       pp(i,i,1:nlev) = _ZERO_
       dd(i,i,1:nlev) = _ZERO_
     end do
   end if

   ! compute benthic variables 
    if (bio_setup>1) then
! scalar
!      do i=1,numbc
!         ppsum = _ZERO_
!         ddsum = _ZERO_
!         do j=1,numbc
!            ppsum=ppsum+ppb(i,j,1)
!            ddsum=ddsum+ddb(i,j,1)
!         end do
!         ccb(i,1)=(ccb(i,1)+dt*ppsum)/(1.+dt*ddsum/ccb(i,1))
!      end do
      ccb = (ccb+dt*sum(ppb(:,:,:),2))/(1.+dt*sum(ddb(:,:,:),2)/ccb)
      ! reset diagonal terms only
      do i=1,numbc
        ppb(i,i,:) = _ZERO_
        ddb(i,i,:) = _ZERO_
      end do
    end if
 
   return
   end subroutine patankar
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Patankar-Runge-Kutta (2nd-order) scheme for geobiochemical
!  models
!
! !INTERFACE:
   subroutine patankar_runge_kutta_2(dt,numc,nlev,h,t,bio_setup,numbc)
!
! !DESCRIPTION:
!
! !USES:
  use bio_var, only: cc,pp,dd,ccb,ppb,ddb
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
   integer, intent(in)                 :: numc,nlev
   integer, intent(in)                 :: numbc,bio_setup
   REALTYPE, intent(in)                :: h(0:nlev)
   REALTYPE, intent(in)                :: t(0:nlev)
!
! !INPUT/OUTPUT PARAMETER:
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
  logical  :: first
  REALTYPE :: ppsum(1:numc,0:nlev),ddsum(1:numc,0:nlev)
  REALTYPE :: cc0(1:numc,0:nlev)
  integer  :: i,j,ci
!EOP
!-----------------------------------------------------------------------
!BOC
   ! store the initial state
   cc0 = cc

   first=.true.
   call process_model(first,numc,nlev,h,t)

   do ci=1,nlev
      do i=1,numc
         ppsum(i,ci)=0.
         ddsum(i,ci)=0.
         do j=1,numc
            ppsum(i,ci)=ppsum(i,ci)+pp(i,j,ci)
            ddsum(i,ci)=ddsum(i,ci)+dd(i,j,ci)
         end do
         cc(i,ci)=(cc0(i,ci)+dt*ppsum(i,ci))/(1.+dt*ddsum(i,ci)/cc0(i,ci))
      end do
   end do

   do i=1,numc
     pp(i,i,1:nlev) = _ZERO_
     dd(i,i,1:nlev) = _ZERO_
   end do

   call process_model(first,numc,nlev,h,t)

   do ci=1,nlev
      do i=1,numc
         do j=1,numc
            ppsum(i,ci)=ppsum(i,ci)+pp(i,j,ci)
            ddsum(i,ci)=ddsum(i,ci)+dd(i,j,ci)
         end do
         cc(i,ci)=(cc0(i,ci)+0.5*dt*ppsum(i,ci))/(1.+0.5*dt*ddsum(i,ci)/cc(i,ci))
      end do
   end do

   do i=1,numc
     pp(i,i,1:nlev) = _ZERO_
     dd(i,i,1:nlev) = _ZERO_
   end do

   return
   end subroutine patankar_runge_kutta_2
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Patankar-Runge-Kutta (4th-order) scheme for geobiochemical
!  models
!
! !INTERFACE:
   subroutine patankar_runge_kutta_4(dt,numc,nlev,h,t,bio_setup,numbc)
!
! !DESCRIPTION:
!
! !USES:
  use bio_var, only: cc,pp,dd,ccb,ppb,ddb
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
   integer, intent(in)                 :: numc,nlev
   integer, intent(in)                 :: numbc,bio_setup
   REALTYPE, intent(in)                :: h(0:nlev)
   REALTYPE, intent(in)                :: t(0:nlev)
!
! !INPUT/OUTPUT PARAMETER:
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
  logical  :: first
  REALTYPE :: ppsum(1:numc,0:nlev),ddsum(1:numc,0:nlev)
  REALTYPE :: ppsum1(1:numc,0:nlev),ddsum1(1:numc,0:nlev)
  REALTYPE :: ppsum2(1:numc,0:nlev),ddsum2(1:numc,0:nlev)
  REALTYPE :: ppsum3(1:numc,0:nlev),ddsum3(1:numc,0:nlev)
  REALTYPE :: cc1(1:numc,0:nlev)
  integer  :: i,j,ci
!EOP
!-----------------------------------------------------------------------
!BOC
   first=.true.
   call process_model(first,numc,nlev,h,t)

   do ci=1,nlev
      do i=1,numc
         ppsum(i,ci)=0.
         ddsum(i,ci)=0.
         do j=1,numc
            ppsum(i,ci)=ppsum(i,ci)+pp(i,j,ci)
            ddsum(i,ci)=ddsum(i,ci)+dd(i,j,ci)
         end do
         cc1(i,ci)=(cc(i,ci)+dt*ppsum(i,ci))/(1.+dt*ddsum(i,ci)/cc(i,ci))
      end do
   end do

   call process_model(first,numc,nlev,cc1,pp,dd,h,t)

   do ci=1,nlev
      do i=1,numc
         ppsum1(i,ci)=0.
         ddsum1(i,ci)=0.
         do j=1,numc
            ppsum1(i,ci)=ppsum1(i,ci)+pp(i,j,ci)
            ddsum1(i,ci)=ddsum1(i,ci)+dd(i,j,ci)
         end do
         cc1(i,ci)=(cc(i,ci)+dt*ppsum1(i,ci))/(1.+dt*ddsum1(i,ci)/cc1(i,ci))
      end do
   end do

   call process_model(first,numc,nlev,cc1,pp,dd,h,t)

   do ci=1,nlev
      do i=1,numc
         ppsum2(i,ci)=0.
         ddsum2(i,ci)=0.
         do j=1,numc
            ppsum2(i,ci)=ppsum2(i,ci)+pp(i,j,ci)
            ddsum2(i,ci)=ddsum2(i,ci)+dd(i,j,ci)
         end do
         cc1(i,ci)=(cc(i,ci)+dt*ppsum2(i,ci))/(1.+dt*ddsum2(i,ci)/cc1(i,ci))
      end do
   end do

   call process_model(first,numc,nlev,cc1,pp,dd,h,t)

   do ci=1,nlev
      do i=1,numc
         ppsum3(i,ci)=0.
         ddsum3(i,ci)=0.
         do j=1,numc
            ppsum3(i,ci)=ppsum3(i,ci)+pp(i,j,ci)
            ddsum3(i,ci)=ddsum3(i,ci)+dd(i,j,ci)
         end do
         ppsum(i,ci)=1./3.*(0.5*ppsum(i,ci)+ppsum1(i,ci)+ppsum2(i,ci)+0.5*ppsum3(i,ci))
         ddsum(i,ci)=1./3.*(0.5*ddsum(i,ci)+ddsum1(i,ci)+ddsum2(i,ci)+0.5*ddsum3(i,ci))
         cc(i,ci)=(cc(i,ci)+dt*ppsum(i,ci))/(1.+dt*ddsum(i,ci)/cc1(i,ci))
      end do
   end do

   return
   end subroutine patankar_runge_kutta_4
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Modified Patankar scheme for geobiochemical models
!
! !INTERFACE:
   subroutine modified_patankar(dt,numc,nlev,h,t,bio_setup,numbc)
!
! !DESCRIPTION:
!
! !USES:
  use bio_var, only: cc,pp,dd,ccb,ppb,ddb
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
   integer, intent(in)                 :: numc,nlev
   integer, intent(in)                 :: numbc,bio_setup
   REALTYPE, intent(in)                :: h(0:nlev)
   REALTYPE, intent(in)                :: t(0:nlev)
!
! !INPUT/OUTPUT PARAMETER:
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
  logical  :: first
  REALTYPE :: a(1:numc,1:numc),r(1:numc)
  integer  :: i,j,ci
!EOP
!-----------------------------------------------------------------------
!BOC
   first=.true.
   call process_model(first,numc,nlev,h,t)

   do ci=1,nlev
      do i=1,numc
         a(i,i)=0.
         do j=1,numc
            a(i,i)=a(i,i)+dd(i,j,ci)
            if (i.ne.j) a(i,j)=-dt*pp(i,j,ci)/cc(j,ci)
         end do
         a(i,i)=dt*a(i,i)/cc(i,ci)
         a(i,i)=1.+a(i,i)
         r(i)=cc(i,ci)+dt*pp(i,i,ci)
      end do
      call matrix(numc,a,r,cc(:,ci))
   end do

   return
   end subroutine modified_patankar
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Modified Patankar-Runge-Kutta (2nd-order) scheme for
!  geobiochemical models
!
! !INTERFACE:
   subroutine modified_patankar_2(dt,numc,nlev,h,t,bio_setup,numbc)
!
! !DESCRIPTION:
!
! !USES:
  use bio_var, only: cc,pp,dd,ccb,ppb,ddb
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
   integer, intent(in)                 :: numc,nlev
   integer, intent(in)                 :: numbc,bio_setup
   REALTYPE, intent(in)                :: h(0:nlev)
   REALTYPE, intent(in)                :: t(0:nlev)
!
! !INPUT/OUTPUT PARAMETER:
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
  logical  :: first
  REALTYPE :: pp1(1:numc,1:numc,0:nlev),dd1(1:numc,1:numc,0:nlev)
  REALTYPE :: a(1:numc,1:numc),r(1:numc)
  REALTYPE :: cc1(1:numc,0:nlev)
  integer  :: i,j,ci
!EOP
!-----------------------------------------------------------------------
!BOC
   first=.true.
   call process_model(first,numc,nlev,h,t)

   do ci=1,nlev
      do i=1,numc
         a(i,i)=0.
         do j=1,numc
            a(i,i)=a(i,i)+dd(i,j,ci)
            if (i.ne.j) a(i,j)=-dt*pp(i,j,ci)/cc(j,ci)
         end do
         a(i,i)=dt*a(i,i)/cc(i,ci)
         a(i,i)=1.+a(i,i)
         r(i)=cc(i,ci)+dt*pp(i,i,ci)
      end do
      call matrix(numc,a,r,cc1(:,ci))
   end do

   call process_model(first,numc,nlev,h,t)

   pp=0.5*(pp+pp1)
   dd=0.5*(dd+dd1)

   do ci=1,nlev
      do i=1,numc
         a(i,i)=0.
         do j=1,numc
            a(i,i)=a(i,i)+dd(i,j,ci)
            if (i.ne.j) a(i,j)=-dt*pp(i,j,ci)/cc1(j,ci)
         end do
         a(i,i)=dt*a(i,i)/cc1(i,ci)
         a(i,i)=1.+a(i,i)
         r(i)=cc(i,ci)+dt*pp(i,i,ci)
      end do
      call matrix(numc,a,r,cc(:,ci))
   end do

   return
   end subroutine modified_patankar_2
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Modified Patankar-Runge-Kutta (4th-order) scheme for
!  geobiochemical models
!
! !INTERFACE:
   subroutine modified_patankar_4(dt,numc,nlev,h,t,bio_setup,numbc)
!
! !DESCRIPTION:
!
! !USES:
  use bio_var, only: cc,pp,dd,ccb,ppb,ddb
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
   integer, intent(in)                 :: numc,nlev
   integer, intent(in)                 :: numbc,bio_setup
   REALTYPE, intent(in)                :: h(0:nlev)
   REALTYPE, intent(in)                :: t(0:nlev)
!
! !INPUT/OUTPUT PARAMETER:
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
  logical  :: first
  REALTYPE :: pp1(1:numc,1:numc,0:nlev),dd1(1:numc,1:numc,0:nlev)
  REALTYPE :: pp2(1:numc,1:numc,0:nlev),dd2(1:numc,1:numc,0:nlev)
  REALTYPE :: pp3(1:numc,1:numc,0:nlev),dd3(1:numc,1:numc,0:nlev)
  REALTYPE :: a(1:numc,1:numc),r(1:numc)
  REALTYPE :: cc1(1:numc,0:nlev)
  integer  :: i,j,ci
!EOP
!-----------------------------------------------------------------------
!BOC
   first=.true.
   call process_model(first,numc,nlev,h,t)

   do ci=1,nlev
      do i=1,numc
         a(i,i)=0.
         do j=1,numc
            a(i,i)=a(i,i)+dd(i,j,ci)
            if (i.ne.j) a(i,j)=-dt*pp(i,j,ci)/cc(j,ci)
         end do
         a(i,i)=dt*a(i,i)/cc(i,ci)
         a(i,i)=1.+a(i,i)
         r(i)=cc(i,ci)+dt*pp(i,i,ci)
      end do
      call matrix(numc,a,r,cc1(:,ci))
   end do

   call process_model(first,numc,nlev,h,t)

   do ci=1,nlev
      do i=1,numc
         a(i,i)=0.
         do j=1,numc
            a(i,i)=a(i,i)+dd1(i,j,ci)
            if (i.ne.j) a(i,j)=-dt*pp1(i,j,ci)/cc1(j,ci)
         end do
         a(i,i)=dt*a(i,i)/cc1(i,ci)
         a(i,i)=1.+a(i,i)
         r(i)=cc(i,ci)+dt*pp1(i,i,ci)
      end do
      call matrix(numc,a,r,cc1(:,ci))
   end do

   call process_model(first,numc,nlev,h,t)

   do ci=1,nlev
      do i=1,numc
         a(i,i)=0.
         do j=1,numc
            a(i,i)=a(i,i)+dd2(i,j,ci)
            if (i.ne.j) a(i,j)=-dt*pp2(i,j,ci)/cc1(j,ci)
         end do
         a(i,i)=dt*a(i,i)/cc1(i,ci)
         a(i,i)=1.+a(i,i)
         r(i)=cc(i,ci)+dt*pp2(i,i,ci)
      end do
      call matrix(numc,a,r,cc1(:,ci))
   end do

   call process_model(first,numc,nlev,h,t)

   pp=1./3.*(0.5*pp+pp1+pp2+0.5*pp3)
   dd=1./3.*(0.5*dd+dd1+dd2+0.5*dd3)

   do ci=1,nlev
      do i=1,numc
         a(i,i)=0.
         do j=1,numc
            a(i,i)=a(i,i)+dd(i,j,ci)
            if (i.ne.j) a(i,j)=-dt*pp(i,j,ci)/cc1(j,ci)
         end do
         a(i,i)=dt*a(i,i)/cc1(i,ci)
         a(i,i)=1.+a(i,i)
         r(i)=cc(i,ci)+dt*pp(i,i,ci)
      end do
      call matrix(numc,a,r,cc(:,ci))
   end do

   return
   end subroutine modified_patankar_4
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: First-order extended modified Patankar scheme for
!  geobiochemical models. Submitted to Applied Numerical Mathematics
!  (2005), authors: Bruggeman, Burchard, Kooi, Sommeijer.
!
! !INTERFACE:
   subroutine emp_1(dt,numc,nlev,h,t,bio_setup,numbc)
!
! !DESCRIPTION:
!
! !USES:
  use bio_var, only: cc,pp,dd,ccb,ppb,ddb
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
   integer, intent(in)                 :: numc,nlev
   integer, intent(in)                 :: numbc,bio_setup
   REALTYPE, intent(in)                :: h(0:nlev)
   REALTYPE, intent(in)                :: t(0:nlev)
!
! !INPUT/OUTPUT PARAMETER:
!
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
  logical  :: first
  integer  :: i,ci
!BEGIN_BFM:MAV TMP solution for 1 benthic layer only
  integer, parameter  :: nlevb=1
!END_BFM
  REALTYPE :: pi, derivative(1:numc)
!EOP
!-----------------------------------------------------------------------
!BOC
   first=.true.
   call process_model(first,numc,nlev,h,t)

   ! integrate pelagic variables 
   if (bio_setup/=2) then
     do ci=1,nlev
       derivative(:) = sum(pp(:,:,ci),2)-sum(dd(:,:,ci),2)
       call findpi_bisection(numc, cc(:,ci), derivative(:), dt, 1.d-9, pi)
       cc(:,ci) = cc(:,ci) + dt*derivative(:)*pi
     end do
     ! reset diagonal terms only
     do i=1,numc
       pp(i,i,1:nlev) = _ZERO_
       dd(i,i,1:nlev) = _ZERO_
     end do
   end if
   ! integrate benthic variables 
   if (bio_setup>1) then
     do ci=1,nlevb
       derivative(:) = sum(ppb(:,:,ci),2)-sum(ddb(:,:,ci),2)
       call findpi_bisection(numbc, ccb(:,ci), derivative(:), dt, 1.d-9, pi)
       ccb(:,ci) = ccb(:,ci) + dt*derivative(:)*pi
     end do
     ! reset diagonal terms only
     do i=1,numbc
       ppb(i,i,:) = _ZERO_
       ddb(i,i,:) = _ZERO_
     end do
   end if

   return
   end subroutine emp_1
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Second-order extended modified Patankar scheme for geobiochemical
! models. Submitted to Applied Numerical Mathematics (2005),
! authors: Bruggeman, Burchard, Kooi, Sommeijer.
!
! !INTERFACE:
   subroutine emp_2(dt,numc,nlev,h,t,bio_setup,numbc)
!
! !DESCRIPTION:
!
! !USES:
  use bio_var, only: cc,pp,dd,ccb,ppb,ddb
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
   integer, intent(in)                 :: numc,nlev
   integer, intent(in)                 :: numbc,bio_setup
   REALTYPE, intent(in)                :: h(0:nlev)
   REALTYPE, intent(in)                :: t(0:nlev)
!
! !INPUT/OUTPUT PARAMETER:
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
!BEGIN_BFM:MAV TMP solution for 1 benthic layer only
  integer, parameter  :: nlevb=1
!END_BFM
  logical  :: first
  integer  :: i,ci
  REALTYPE :: pi, rhs(1:numc,0:nlev), cc_med(1:numc,0:nlev)
  REALTYPE :: rhsb(1:numbc,1:nlevb), ccb_med(1:numc,1:nlevb)
!EOP
!-----------------------------------------------------------------------
!BOC
   first=.true.

   call process_model(first,numc,nlev,h,t)

   ! integrate pelagic variables 
   if (bio_setup/=2) then
     do ci=1,nlev
       rhs(:,ci) = sum(pp(:,:,ci),2) - sum(dd(:,:,ci),2)
       call findpi_bisection(numc, cc(:,ci), rhs(:,ci), dt, 1.d-9, pi)
       cc_med(:,ci) = cc(:,ci) + dt*rhs(:,ci)*pi
     end do
     ! reset diagonal terms only
     do ci=1,numc
       pp(ci,ci,1:nlev) = _ZERO_
       dd(ci,ci,1:nlev) = _ZERO_
     end do
   end if
   ! integrate benthic variables 
   if (bio_setup>1) then
     do ci=1,nlevb
       rhsb(:,ci) = sum(ppb(:,:,ci),2) - sum(ddb(:,:,ci),2)
       call findpi_bisection(numbc, cc(:,ci), rhsb(:,ci), dt, 1.d-9, pi)
       ccb_med(:,ci) = ccb(:,ci) + dt*rhsb(:,ci)*pi
     end do
     ! reset diagonal terms only
     do ci=1,numbc
       ppb(ci,ci,1) = _ZERO_
       ddb(ci,ci,1) = _ZERO_
     end do
   end if

   call process_model(first,numc,nlev,h,t)

   ! integrate pelagic variables 
   if (bio_setup/=2) then
     do ci=1,nlev
       rhs(:,ci) = 0.5 * (rhs(:,ci) + sum(pp(:,:,ci),2) - sum(dd(:,:,ci),2))
       ! Correct for the state variables that will be included in 'pi'.
       do i=1,numc
         if (rhs(i,ci) .lt. 0.) rhs(:,ci) = rhs(:,ci) * cc(i,ci)/cc_med(i,ci)
       end do
       call findpi_bisection(numc, cc(:,ci), rhs(:,ci), dt, 1.d-9, pi)
       cc(:,ci) = cc(:,ci) + dt*rhs(:,ci)*pi
     end do ! ci (z-levels)
     ! reset diagonal terms only
     do ci=1,numc
       pp(ci,ci,1:nlev) = _ZERO_
       dd(ci,ci,1:nlev) = _ZERO_
     end do
   end if
   ! integrate benthic variables 
   if (bio_setup>1) then
     do ci=1,nlevb
       rhsb(:,ci) = 0.5 * (rhsb(:,ci) + sum(ppb(:,:,ci),2) - sum(ddb(:,:,ci),2))
       ! Correct for the state variables that will be included in 'pi'.
       do i=1,numbc
         if (rhsb(i,ci) .lt. 0.) rhsb(:,ci) = rhsb(:,ci) * ccb(i,ci)/ccb_med(i,ci)
       end do
       call findpi_bisection(numbc, ccb(:,ci), rhsb(:,ci), dt, 1.d-9, pi)
       ccb(:,ci) = ccb(:,ci) + dt*rhsb(:,ci)*pi
     end do ! ci (z-levels)
     ! reset diagonal terms only
     do ci=1,numbc
       ppb(ci,ci,1) = _ZERO_
       ddb(ci,ci,1) = _ZERO_
     end do
   end if

   return
   end subroutine emp_2
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Finds the Extended Modified Patankar product term 'pi'
!  with the bisection technique.
!
! !INTERFACE:
   subroutine findpi_bisection(numc, cc, derivative, dt, accuracy, pi)
!
! !DESCRIPTION:
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   integer, intent(in)   :: numc
   REALTYPE, intent(in)  :: cc(1:numc), derivative(1:numc)
   REALTYPE, intent(in)  :: dt, accuracy
!
! !OUTPUT PARAMETER:
   REALTYPE, intent(out) :: pi
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   REALTYPE :: pileft, piright, fnow
   REALTYPE :: relderivative(1:numc)
   integer  :: iter, i, potnegcount
!EOP
!-----------------------------------------------------------------------
!BOC
! Sort the supplied derivatives (find out which are negative).
   potnegcount = 0
   piright = 1.
   do i=1,numc

      if (derivative(i).lt.0.) then
!        State variable could become zero or less; include it in the
!        J set of the EMP scheme.
         if (cc(i).eq.0.) write (*,*) "Error: state variable ",i," is zero and has negative derivative!"
         potnegcount = potnegcount+1
         relderivative(potnegcount) = dt*derivative(i)/cc(i)

!        Derivative is negative, and therefore places an upper bound on pi.
         if (-1./relderivative(potnegcount).lt.piright) piright = -1./relderivative(potnegcount)
     end if

   end do

   if (potnegcount.eq.0) then
!     All derivatives are positive, just do Euler.
      pi = 1.0
      return
   end if

   pileft = 0.      ! polynomial(0) = 1

!  Determine maximum number of bisection iterations from
!  requested accuracy.
!  maxiter = -int(ceiling(dlog10(accuracy)/dlog10(2.D0)))

   do iter=1,20
!     New pi to test is middle of current pi-domain.
      pi = 0.5*(piright+pileft)

!     Calculate polynomial value.
      fnow = 1.
      do i=1,potnegcount
         fnow = fnow*(1.+relderivative(i)*pi)
      end do

      if (fnow>pi) then
!        Polynomial(pi)>0; we have a new left bound for pi.
         pileft = pi
      elseif (fnow<pi) then
!       Polynomial(pi)<0; we have a new right bound for pi.
        piright = pi
      else
!       Freak occurrence: polynomial(pi)=0, we happened to pinpoint
!       the exact pi.
        exit
      end if
!     Check if we now pi accurately enough (accuracy refers to the
!     number of decimals we know).
      if ((piright-pileft)/pi<accuracy) exit
   end do

!  Low pi values imply very large negative relative derivative. This happens
!  for stiff systems (or very high delta_t), and for non-positive systems.
!  Then EMP is not suitable (it will stall state variable values), so warn user.
   if (pi.lt.1.d-4) then
     write (*,*) "Warning: small pi=",pi," in Extended Modified Patankar slows down system!"
!    write (*,*) "relative derivatives: ",derivative(:)*dt/cc(:)
!    write (*,*) "You system may be stiff or non-positive, or you time step is too large."
!    stop "ode_solvers::findpi_bisection"
   end if

   return

   end subroutine findpi_bisection
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Matrix solver
!
! !INTERFACE:
   subroutine matrix(n,a,r,c)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: n
!
! INPUT/OUTPUT PARAMETERS:
  REALTYPE                             :: a(1:n,1:n),r(1:n)
!
! OUTPUT PARAMETERS:
  REALTYPE, intent(out)                :: c(1:n)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
  integer  :: i,j,k
!EOP
!-----------------------------------------------------------------------
!BOC
   do i=1,n
      r(i)=r(i)/a(i,i)
      do j=n,i,-1
         a(i,j)=a(i,j)/a(i,i)
      end do
      do k=i+1,n
         r(k)=r(k)-a(k,i)*r(i)
         do j=i+1,n
            a(k,j)=a(k,j)-a(k,i)*a(i,j)
         end do
      end do
   end do

   do i=n,1,-1
      c(i)=r(i)
      do j=i+1,n
         c(i)=c(i)-a(i,j)*c(j)
      end do
   end do

   return
   end subroutine matrix
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
