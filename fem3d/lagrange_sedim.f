c
c $Id: lagrange_larve.f,v 1.2 2008-11-03 10:42:26 georg Exp $
c
c simulates lagrangian sediment transport
c
c revision log :
c
c 06.05.2015    ccf	written from scratch
c
!*******************************************************************
!================================================================
        module lgr_sedim_module
!================================================================

        implicit none
        save

        integer                       :: nls       !number of sediment classes
        real, allocatable             :: gs(:)     !particle grain size [mm]
        real, allocatable             :: fr(:)     !fraction of sediment type [0-1]
        double precision, allocatable :: ws(:)     !particle settling velocity [m/s]

!================================================================
        contains
!================================================================
!*******************************************************************
! initialize sediment classes, settling velocities and grainsizes

        subroutine lgr_sedim_init

        implicit none

        nls = 3

        allocate (gs(nls))
        allocate (ws(nls))
        allocate (fr(nls))

               ! clay  silt  sand
        gs = (/0.004, 0.015, 0.1/)
        ws = (/0.0001, 0.0010, 0.01/)
        fr = (/0.2, 0.3, 0.5/)

        end subroutine lgr_sedim_init

!*******************************************************************

	subroutine lgr_sediment(it)

	use mod_lagrange
	use levels

	implicit none

	include 'param.h'

	integer, intent(in)   :: it
	integer               :: i,ie,lb,lmax
	real                  :: z
	integer, save         :: icall =  0 

!-----------------------------------------------------
! initialization
!-----------------------------------------------------

	if( icall .eq. 0 ) then
	  write(6,*)'Initialization on lagrangian sediment model'
	  icall = 1
	end if

!-----------------------------------------------------
! normal call
!-----------------------------------------------------

!-----------------------------------------------------
! do not advect particles on bottom of last layer
!-----------------------------------------------------
	
	do i = 1,nbdy
	  ie   = lgr_ar(i)%ie
          lmax = ilhv(ie)
	  lb   = lgr_ar(i)%l
          z    = lgr_ar(i)%z		!rel. depth   0=surface  1=bottom

	  if ( lb.eq.lmax .and. z.eq.1. ) lgr_ar(i)%l = -1	!particle on bottom
	end do

!----------------------------------------------------------------
! end of routine
!----------------------------------------------------------------

	end subroutine lgr_sediment

!*******************************************************************
! assign sediment class, settlign velocity and grainsize to particle
! called by subroutine lgr_set_properties

        subroutine lgr_set_sedim(sc,sv,sg)

        implicit none

        integer, intent(out)          :: sc	!sediment class
        double precision, intent(out) :: sv	!settling velocity [m/s]
	real, intent(out)             :: sg	!sediment grainsize

	real                          :: r
	integer                       :: nr,na
	integer, parameter            :: nperc=100
	integer, dimension(nperc)     :: array
	integer                       :: n,n1,n2
	
!----------------------------------------------------------------
! Create array with index values according to fractions
!----------------------------------------------------------------

	n1 = 0
	n2 = 0
	do n = 1,nls
	  if (n .eq. 1) then
	    n1 = 1
	  else
	    n1 = n1 + fr(n-1) * nperc
	  end if
	  n2 = n2 + fr(n) * nperc
	  array(n1:n2) = n
	end do

!----------------------------------------------------------------
! Get random weighted index from array
!----------------------------------------------------------------

        call random_number(r)
	na = r*nperc + 1
	na = min(nperc,na)
	nr = array(na)

!----------------------------------------------------------------
! Set properties
!----------------------------------------------------------------

	sv = ws(nr)
	sg = gs(nr)
	sc = nr

	end subroutine lgr_set_sedim

!*******************************************************************
!================================================================
        end module lgr_sedim_module
!================================================================

