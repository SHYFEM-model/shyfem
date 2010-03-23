#include"cppdefs.h"
!-----------------------------------------------------------------------
! !ROUTINE: to call do_turbulence from SHYFEM

   subroutine do_gotm_turb(nlev,dt,depth,u_taus,u_taub,z0s,z0b,h,       &
			   NN,SS,num1d,nuh1d,tke1d,eps1d,L1d)

! !USES:
   use turbulence,  only: do_turbulence
   use turbulence,  only: num,nuh,tke,eps,L

   IMPLICIT NONE

! !INPUT PARAMETERS:

!  number of vertical layers
   integer,  intent(in)                :: nlev

!  time step (s)
   REALTYPE, intent(in)                :: dt

!  distance between surface
!  and bottom(m)
   REALTYPE, intent(in)      	       :: depth

!  surface and bottom
!  friction velocity (m/s)
   REALTYPE, intent(in)                :: u_taus,u_taub

!  surface and bottom
!  roughness length (m)
   REALTYPE, intent(in)                :: z0s,z0b

!  layer thickness (m)
   REALTYPE, intent(in)                :: h(0:nlev)

!  boyancy frequency squared (1/s^2)
   REALTYPE, intent(in)                :: NN(0:nlev)

!  shear-frequency squared (1/s^2)
   REALTYPE, intent(in)                :: SS(0:nlev)

! !INPUT,OUTPUT PARAMETERS:

!  turbulent diffusivities
!  of momentum
   REALTYPE			:: num1d(0:nlev)	

!  turbulent diffusivities
!  of temperature
   REALTYPE			:: nuh1d(0:nlev)

!  TKE
   REALTYPE			:: tke1d(0:nlev)

!  rate of dissipation
   REALTYPE			:: eps1d(0:nlev)

!  turbulent length-scale
   REALTYPE			:: L1d(0:nlev)

! !Update GOTM variables

   num = num1d
   nuh = nuh1d
   tke = tke1d
   eps = eps1d
   L = L1d

! !Call GOTM turbulence routine

   call do_turbulence(nlev,dt,depth,u_taus,u_taub,z0s,z0b,h,      &
                            NN,SS)

! !Update variables to pass back to SHYFEM

   num1d = num
   nuh1d = nuh
   tke1d = tke
   eps1d = eps
   L1d = L

   end

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! !ROUTINE: to call init_turbulence from SHYFEM

   subroutine init_gotm_turb(namlst,fn,nlev)

   use turbulence, only: init_turbulence
   use mtridiagonal,only: init_tridiagonal

! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,          intent(in)        :: namlst
   character(len=*), intent(in)        :: fn
   integer,          intent(in)        :: nlev

   call init_tridiagonal(nlev)

   call init_turbulence(namlst,fn,nlev)

   end

!-----------------------------------------------------------------------

   subroutine report_gotm(iu)

   !use turbulence, only: report_model
   use turbulence

   implicit none

   integer iu

   write(iu,*) ' '
   write(iu,*) '--------------------------------------------------------'
   write(iu,*) 'You are using the k-epsilon model'
   write(iu,*) 'with the following properties:'
   write(iu,*) ' '
   write(iu,*) 'ce1                                  =', ce1
   write(iu,*) 'ce2                                  =', ce2
   write(iu,*) 'ce3minus                             =', ce3minus
   write(iu,*) 'ce3plus                              =', ce3plus
   write(iu,*) 'sig_k                                =', sig_k
   write(iu,*) 'sig_e                                =', sig_e
   write(iu,*) ' '
   write(iu,*) 'Value of the stability function'
   write(iu,*) 'in the log-law,                   cm0 =', cm0
   write(iu,*) 'in shear-free turbulence,        cmsf =', cmsf
   write(iu,*) ' '
   write(iu,*) 'von Karman constant,           kappa =', kappa
!   write(iu,*) 'homogeneous decay rate,            d =', gen_d
!   write(iu,*) 'spatial decay rate (no shear), alpha =', gen_alpha
!   write(iu,*) 'length-scale slope (no shear),     L =', gen_l
!   write(iu,*) 'steady-state Richardson-number, Ri_st=', ri_st
   write(iu,*) '--------------------------------------------------------'
   write(iu,*) ' '

   end

! changed to public:

!   REALTYPE, public                              :: ri_st=0.25

!   REALTYPE, public                              :: gen_d=-1.2
!   REALTYPE, public                              :: gen_alpha=-2.0
!   REALTYPE, public                              :: gen_l=0.2

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! next routines are used for debug
!
! at start call save_gotm_init (time step, do loop, whatever)
! data is saved to unit 44
! if wanted, after gotm call, with save_gotm_write the data is written to 45
! save data (in gotm) with save_gotm_array etc..
!
!-----------------------------------------------------------------------

	subroutine save_gotm_init()

	implicit none

	rewind(44)

	end

!-----------------------------------------------------------------------

	subroutine save_gotm_write()

	implicit none

	character*80 line

	close(44)
	open(2,file='fort.44',status='old',form='formatted',err=3)

    1	continue
	  read(2,'(a)',end=2) line
	  write(45,'(a)') line
	  goto 1
    2	continue

	close(2)

    3	continue
	end

!-----------------------------------------------------------------------

	subroutine save_gotm_array(text,nlev,val)

	implicit none

	character*(*) text
	integer nlev
	double precision val(0:nlev)

	integer i

	write(44,*) text,(val(i),i=0,nlev)

	end

!-----------------------------------------------------------------------

	subroutine save_gotm_val(text,val)

	implicit none

	character*(*) text
	double precision val

	write(44,*) text,val

	end

!-----------------------------------------------------------------------

	subroutine save_gotm_int(text,ival)

	implicit none

	character*(*) text
	integer ival

	write(44,*) text,ival

	end

!-----------------------------------------------------------------------












