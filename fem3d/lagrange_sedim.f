
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

c simulates lagrangian sediment transport
c
c revision log :
c
c 06.05.2015	ccf	written from scratch
c 21.05.2015	ggu	changed VERS_7_1_11
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 16.11.2015	ggu	changed VERS_7_3_14
c 19.05.2017	ccf	include deposition and resuspension dynamics
c 25.10.2018	ggu	changed VERS_7_5_51
c 14.02.2019	ggu	changed VERS_7_5_56
c 16.02.2019	ggu	changed VERS_7_5_60
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
        real                          :: tce     !critical threshold for erosion [N/m**2]

!================================================================
        contains
!================================================================
!*******************************************************************
! initialize sediment classes, settling velocities and grainsizes
! still set by hand

        subroutine lgr_sedim_init

        implicit none

	real                    :: maxtcd
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

        gs = (/0.001, 0.005, 0.010/)
        fr = (/0.34,0.33,0.33/)
        ws = (/0.003, 0.01, 0.03/)	!mm/s 
        ws = ws*0.001           	!from mm/s to m/s

        tce = 0.2
	tcd = 0.5*2800.*ws**1.03	!from sedtrans05

        maxtcd = maxval(tcd)
        if ( maxtcd > tce ) then
          write(6,*) 'maxtcd,tce: ',maxtcd,tce
          write(6,*) 'maxtcd must be < tce'
          stop 'error stop lgr_sedim_init'
	end if

        end subroutine lgr_sedim_init

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
	integer, allocatable          :: array(:)
	integer                       :: n,n1,n2
	
!----------------------------------------------------------------
! Create array with index values according to fractions
!----------------------------------------------------------------

        allocate(array(nperc))

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
	sg = ws(nr)
	sg = gs(nr)
	sc = nr

	nc = 1

	end subroutine lgr_set_sedim

!================================================================
        end module lgr_sedim_module
!================================================================

!*******************************************************************
! Sediment dynamics (must be that tcd < tce)
! if tau < tcd       particle sinks and could deposit  
! if tcd < tau < tce particle stays in suspension 
! if tau > tce       particle on bottom could be resuspended 
!
! Probablility of deposition/erosion is function of tau

	subroutine lgr_sediment

        use basin
        use mod_diff_visc_fric
	use mod_lagrange
	use lgr_sedim_module
	use levels
	use mod_bstress, only : taubot
	use mod_sediment, only : isedi

	implicit none

	integer                 :: i,ie,l
	integer			:: lb,lmax
	real                    :: z
	integer, save           :: icall = 0 
        integer                 :: nr
	double precision	:: sv
	real		        :: tau,tcdp
	real                    :: r,perc
	integer			:: id
	real			:: rfa = 0.5   !probabilistic factor [0.3-0.5]
        double precision 	:: xi(3)
        integer          	:: ii(3)

!-----------------------------------------------------
! initialization
!-----------------------------------------------------

	if( icall .eq. 0 ) then
	  write(6,*)'Initialization on lagrangian sediment model'
	  bcompress = .false.
	  icall = 1
	end if

!-----------------------------------------------------
! compute bottom shear stresses 
!-----------------------------------------------------

	if ( isedi == 0 ) call bottom_stress(taubot)

!-----------------------------------------------------
! Loop over particles
!-----------------------------------------------------

	do i = 1,nbdy
          id = lgr_ar(i)%id
	  ie   = lgr_ar(i)%actual%ie
          if ( ie < 1 ) cycle
          lmax = ilhv(ie)
	  lb   = lgr_ar(i)%actual%l		!layer were particle is located
	  if ( lb > 0 .and. lb /= lmax ) cycle	!skip if particle is not on bottom layer 
	  sv   = lgr_ar(i)%sinking
	  if ( sv == 0.d0 ) cycle	!skip if it has settling velocity = 0
          z    = lgr_ar(i)%actual%z	!rel. depth   0=surface  1=bottom
	  nr   = lgr_ar(i)%type		!sediment class id
	  tcdp = tcd(nr)		!deposition threshold
          xi   = lgr_ar(i)%actual%xi(:)
          ii   = nen3v(:,ie)
          tau  = sum(taubot(ii)*xi)	!bottom stress at particle position
          call random_number(r)		!random number

!         -----------------------------------------------------
!         Deposition of particle 
!         -----------------------------------------------------
          if ( tau < tcdp ) then
            perc = r**(tau/tcdp)
            if ( lb /= -1 .and. z == 1. .and. perc >= rfa ) then 
               lgr_ar(i)%actual%l = -1    		!particle reached bottom
            end if

!         -----------------------------------------------------
!         Erosion of particle, otherwise keep in suspension
!         -----------------------------------------------------
          else if ( tau > tce ) then
            perc = r**(tce/tau)
            if ( lb == -1 .and. perc >= rfa ) then
	      lgr_ar(i)%actual%l = lmax
	      lgr_ar(i)%actual%z = 0.7**((perc-rfa)/(1.-rfa)) !between 0.7 and 1
            end if

          end if

        end do

!----------------------------------------------------------------
! end of routine
!----------------------------------------------------------------

	end subroutine lgr_sediment

!*******************************************************************
