c
c $Id: lagrange_larve.f,v 1.2 2008-11-03 10:42:26 georg Exp $
c
c simulates lagrangian sediment transport
c
c revision log :
c
c 06.05.2015    ccf	written from scratch
c 19.05.2017    ccf	include deposition and resuspension dynamics
c
!*******************************************************************
!================================================================
        module lgr_sedim_module
!================================================================

        implicit none
        save

        integer                       :: nls     !number of sediment classes
        real, allocatable             :: gs(:)   !particle grain size [mm]
        real, allocatable             :: fr(:)   !fraction of sediment type [0-1]
        double precision, allocatable :: ws(:)   !particle settling velocity [mm/s]
        real, allocatable             :: tcd(:)  !critical stress deposition [N/m**2]
        real, parameter            :: tce = 0.1  !critical threshold for erosion [N/m**2]

!================================================================
        contains
!================================================================
!*******************************************************************
! initialize sediment classes, settling velocities and grainsizes

        subroutine lgr_sedim_init

        implicit none

        nls = 3
        !nls = 1

        allocate (gs(nls))
        allocate (ws(nls))
        allocate (fr(nls))
        allocate (tcd(nls))

               ! clay  silt  sand
        !gs = (/0.004, 0.015, 0.1/)
        !ws = (/0.0001, 0.0010, 0.01/)
        !fr = (/0.2, 0.3, 0.5/)

        !gs = (/0.000, 0.005, 0.015/)
        !ws = (/0.0, 0.5, 5.0/)
        !fr = (/0.34,0.33,0.33/)
        !tcd = (/0.00,0.001,0.05/)

        gs = (/0.001, 0.005, 0.010/)
        ws = (/0.1, 0.5, 2.0/)
        fr = (/0.34,0.33,0.33/)
        tcd = (/0.001,0.005,0.05/)

        end subroutine lgr_sedim_init

!*******************************************************************
! Sediment dynamics (must be that tcd < tce)
! if tau < tcd       particle sink and could deposit  
! if tcd < tau < tce particle stay in suspension 
! if tau > tce       particle on bottom could be resuspended 
!
! Probablility of deposition/erosion is function of tau

	subroutine lgr_sediment

        use basin
        use mod_diff_visc_fric
	use mod_lagrange
	use levels

	implicit none

	include 'param.h'

	integer                 :: i,ie,ii,k,l
	integer			:: lb,lmax,nlev,lbe
	real                    :: z
	integer, save           :: icall =  0 
        integer                 :: nr
        integer                 :: mode
	double precision	:: sv
	real, dimension(nel)    :: tauev     !bottom current stress for elements
	real, dimension(nkn)    :: tauw      !wave bottom stress for node
	integer, dimension(nlv) :: sdyn	     !flag for erosion/deposition/suspension
        real                    :: taue(nel) !bottom shear stress in element
	real		        :: tc,tw,tau
	real                    :: r,perc,maxtcd,tcdc
	integer			:: id

!-----------------------------------------------------
! initialization
!-----------------------------------------------------

	if( icall .eq. 0 ) then
	  write(6,*)'Initialization on lagrangian sediment model'
	  bcompress = .true.
	  icall = 1
	end if

!-----------------------------------------------------
! normal call
!-----------------------------------------------------

!-----------------------------------------------------
! compute stresses and sediment dynamics at elements
!-----------------------------------------------------

	!call current_bottom_stress_el(tauev)
        !call wave_bottom_stress(tauw)

	maxtcd = maxval(tcd)
	do ie = 1,nel

!         -----------------------------------------------------
!         compute combined current wave stress at elements
!         -----------------------------------------------------
	   tc = tauev(ie)
	   tw = 0.
           do ii=1,3
             k = nen3v(ii,ie)
	     tw = tw + tauw(k)
	   end do
           tw = tw / 3.
           if( tc+tw == 0. ) then
             tau = 0.
           else
             tau = tc * ( 1. + 1.2 * ( tw/(tc+tw) )**3.2 )
           end if

	   taue(ie) = tau

!         -----------------------------------------------------
!         case sensitive for sediment dynamics at elements
!         -----------------------------------------------------

   	  if ( tau < maxtcd ) then
!           -----------------------------------------------------
!           Deposition allowed
!           -----------------------------------------------------
	    sdyn(ie) = 2
          else if( tau > tce ) then
!           -----------------------------------------------------
!           Erosion allowed
!           -----------------------------------------------------
	    sdyn(ie) = 1
	  else
!           -----------------------------------------------------
!           Particle in suspension
!           -----------------------------------------------------
	    sdyn(ie) = 0
	  end if 
	end do

!-----------------------------------------------------
! Loop over particles
!-----------------------------------------------------

	do i = 1,nbdy
          id = lgr_ar(i)%id
	  ie   = lgr_ar(i)%actual%ie
          if ( ie < 1 ) cycle
          lmax = ilhv(ie)
	  lb   = lgr_ar(i)%actual%l		!layer were particle is located
	  if ( lb > 0 .and. lb /= lmax ) cycle	!skip if particle is not on bottom
	  sv   = lgr_ar(i)%sinking
	  if ( sv == 0.d0 ) cycle	!skip if it has settling velocity = 0
          z    = lgr_ar(i)%actual%z	!rel. depth   0=surface  1=bottom
	  nr   = lgr_ar(i)%type		!sediment class id
	  tcdc = tcd(nr)
          mode = sdyn(ie) 
          tau  = taue(ie)
          call random_number(r)

!         -----------------------------------------------------
!         Set vertical velocity 
!         -----------------------------------------------------
	  !lgr_ar(i)%sv = 0. 

!         -----------------------------------------------------
!         Case sensitive for sediment dynamics
!         -----------------------------------------------------
          select case ( mode )
            case ( 0 )
!             -----------------------------------------------------
!             Do nothing, keep particle in suspension
!             -----------------------------------------------------
            case ( 1 )
!             -----------------------------------------------------
!             Erosion if particle on bottom, otherwise keep in suspension
!             -----------------------------------------------------
              perc = r**(tce/tau)
	      if ( lb == -1 .and. perc > 0.5 ) then
		lgr_ar(i)%actual%l = lmax
		lgr_ar(i)%actual%z = 0.5**((perc-0.5)/0.5) !between 1 and 0.5
	      end if
            case ( 2 )
!             -----------------------------------------------------
!             Deposition of particle on bottom
!             -----------------------------------------------------
	      if ( z == 1. .and. tau < tcdc ) then 	!deposition for that specific class
	        lgr_ar(i)%actual%l = -1    		!particle reached bottom
	      end if
            case default
	      write(*,*) 'mode: ', mode
              stop 'error stop lgr_sediment: internal error'
          end select
	end do

!----------------------------------------------------------------
! end of routine
!----------------------------------------------------------------

	end subroutine lgr_sediment

!*******************************************************************
! assign sediment class, settlign velocity and grainsize to particle
! called by subroutine lgr_set_properties

        subroutine lgr_set_sedim(sc,sv,sg,nc)

        implicit none

        integer, intent(out)          :: sc	!sediment class
        double precision, intent(out) :: sv	!settling velocity [m/s]
	real, intent(out)             :: sg	!sediment grainsize
	integer, intent(out)          :: nc	!number of custom properties

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

	sv = ws(nr)*0.001           !from mm/s to m/s
	sg = ws(nr)
	!sg = gs(nr)
	sc = nr

	nc = 1

	end subroutine lgr_set_sedim

!*******************************************************************
!================================================================
        end module lgr_sedim_module
!================================================================

