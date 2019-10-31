
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

! bottom shear stress subroutines
!
! revision log :
!
! 30.11.2018	ggu&ccf	moved from subwaves.f, subn35.f, subssed.f
! 30.11.2018	ccf	create module mod_bstress
! 14.02.2019	ggu	changed VERS_7_5_56
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	changed VERS_7_5_61
! 31.10.2019	ggu	introduced bdry, bugfix for wave bottom stress (depth<0)
!
!*****************************************************************

!==================================================================
        module mod_bstress
!==================================================================

        implicit none

        integer, private, save  :: nkn_stress = 0
        real, allocatable, save :: taubot(:)     !bottom shear stress [N/m2]

        logical, save  		:: bsdry = .true.  !compute stress in dry areas

        integer, save  		:: ibstress  = 0 !parameter for stress module
        double precision, save  :: da_str(4) = 0 !for output

!==================================================================
        contains
!==================================================================

        subroutine mod_bstress_init(nkn)

        integer  :: nkn

        if( nkn == nkn_stress ) return

        if( nkn > 0 ) then
          if( nkn == 0 ) then
            write(6,*) 'nkn: ',nkn
            stop 'error stop mod_bstress_init: incompatible parameters'
          end if
        end if

        if( nkn_stress > 0 ) then
          deallocate(taubot)
        end if

        nkn_stress = nkn

        if( nkn == 0 ) return

        allocate(taubot(nkn))

        taubot = 0.

        end subroutine mod_bstress_init

!==================================================================
        end module mod_bstress
!==================================================================

        subroutine init_bstress

! initialize output file 

        use mod_bstress

        implicit none

        real 		:: getpar             !get parameter function
        logical 	:: has_output_d
        integer 	:: id,nvar
	integer		:: isedi

        ibstress = nint(getpar('ibstrs'))
        isedi    = nint(getpar('isedi'))

        if ( isedi > 0 ) ibstress = 0	!stress already written in .sed file

        if( ibstress <= 0 ) return

        nvar = 1
        call init_output_d('itmout','idtout',da_str)
        if( has_output_d(da_str) ) then
          call shyfem_init_scalar_file('bstress',nvar,.true.,id)
          da_str(4) = id
        end if

        end subroutine init_bstress

!**************************************************************

        subroutine bstress

! compute bottom stress and write to file 

        use mod_bstress

        implicit none

        logical 	  :: next_output_d
        integer 	  :: id,nvar
        double precision  :: dtime

        if( ibstress <= 0 ) return

        call bottom_stress(taubot)

        if( next_output_d(da_str) ) then
          id = nint(da_str(4))
          call get_act_dtime(dtime)
          call shy_write_scalar_record2d(id,dtime,893,taubot)
        end if

        end subroutine bstress

!**************************************************************

        subroutine bottom_stress(taubot)

! computes bottom stress at nodes (current + waves)

        use basin

        implicit none

        real, intent(out) 		:: taubot(nkn)

        real, allocatable		:: taucur(:)
        real, allocatable		:: tauwav(:)
        integer 			:: k
        real 				:: tc,tw,tm

        allocate(taucur(nkn))
        allocate(tauwav(nkn))

        call current_bottom_stress(taucur)
        call wave_bottom_stress(tauwav)

        do k=1,nkn
          tc = taucur(k)
          tw = tauwav(k)
          if( tc+tw == 0. ) then
            tm = 0.
          else
            tm = tc * ( 1. + 1.2 * ( tw/(tc+tw) )**3.2 )
            tm = max(tm,tw)
          end if
          taubot(k) = tm
        end do

        end

!***********************************************************

        subroutine current_bottom_stress(taucur)

! computes bottom stress at nodes due to currents

        use basin
        use levels
        use evgeom
        use mod_diff_visc_fric
	use mod_geom_dynamic
        use mod_ts
        use mod_bstress

        implicit none

        include 'pkonst.h'

        real, intent(out)	:: taucur(nkn)

        integer 		:: ie,k,ii,lmax
        real 			:: area,rho,bnstress
        real, allocatable	:: aux(:)

	allocate(aux(nkn))

        call bottom_friction    !bottom stress on elements (bnstressv)

        taucur = 0.
        aux = 0.

        do ie=1,nel
          lmax = ilhv(ie)
          area = ev(10,ie)
          bnstress = bnstressv(ie)
          do ii=1,3
            k = nen3v(ii,ie)
            rho = rowass + rhov(lmax,k)
	    if( bsdry .or. iwegv(ie) == 0 ) then
              taucur(k) = taucur(k) + rho * bnstress * area
              aux(k) = aux(k) + area
	    end if
          end do
        end do

        where( aux > 0. ) taucur = taucur / aux

        end

!*********************************************************************

	subroutine wave_bottom_stress(tauwav)

! computes bottom stress from waves (on nodes)

	use basin
	use mod_parwaves
	use mod_waves
        use mod_depth
        use mod_hydro, only : znv

	implicit none

	real, intent(out)	:: tauwav(nkn)

	integer			:: k
	real 			:: h,p,depth

        tauwav = 0.
	if( iwave == 0 ) return

	do k = 1,nkn
	  h = waveh(k)
	  p = wavep(k)
	  depth = hkv(k) + znv(k)
	  call compute_wave_bottom_stress(h,p,depth,z0,tauwav(k))
	end do

	end subroutine wave_bottom_stress

!*****************************************************************

	subroutine compute_wave_bottom_stress(h,p,depth,z0,tauw)

	implicit none

	real h,p	!wave height and period
	real depth	!depth of water column
	real z0		!bottom roughness
	real tauw	!stress at bottom (return)

	include 'pkonst.h'

	real omega,zeta,a,eta,fw,uw
	real, parameter :: pi = 3.14159

	tauw = 0.
	if( p == 0. ) return
	if( depth <= 0. ) return		!bug fix

        omega = 2.*pi/p
        zeta = omega * omega * depth / grav
        if( zeta .lt. 1. ) then
          eta = sqrt(zeta) * ( 1. + 0.2 * zeta )
        else
          eta = zeta * ( 1. + 0.2 * exp(2.-2.*zeta) )
        end if

	if( eta > 80. ) eta = 80.		!GGUZ0

        uw = pi * h / ( p * sinh(eta) )
        a = uw * p / (2.*pi)
        if( a .gt. 1.e-5 ) then			!GGUZ0
          fw = 1.39 * (z0/a)**0.52
        else
          fw = 0.
        end if

        tauw = 0.5 * rowass * fw * uw * uw

	end subroutine compute_wave_bottom_stress

c******************************************************************

	subroutine compute_bstress_dry(bset)

        use mod_bstress

	implicit none

	logical bset

	bsdry = bset

	end

c******************************************************************

