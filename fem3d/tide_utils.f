
!--------------------------------------------------------------------------
!
!    Copyright (C) 2018-2019  Christian Ferrarin
!    Copyright (C) 2019  Georg Umgiesser
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

! subroutines for computing tidal potential and tidal analysis/prediction
! Considered constituents:
!    - semi-diurnal: M2,S2,N2,K2,NU2,MU2,L2,T2
!    - diurnal:      K1,O1,P1,Q1,J1,OO1,S1
!    - long-period:  MF,MM,MSM,MSF,SSA,SA
!
!********************************************************************
!
! revision log :
!
! 29.03.2018	ccf	converted previous routines to tide_utils
! 13.02.2019	ccf	introduced more constituents
! 13.03.2019	ggu	changed VERS_7_5_61
! 15.07.2020	ccf&ggu	new default value for itidana=-1 (no analysis)
! 22.09.2020	ggu	correct warnings for PGI compiler
!
!********************************************************************

!==================================================================
        module tide
!==================================================================

        implicit none

        integer, private, save  :: nkn_tide = 0
        integer, save           :: rtide    !parameter for the tidal potential
        real, save           	:: ltidec   !calibration coefficient for load tide
        real, allocatable, save :: zeqv(:)  !tidal equilibrium
        integer, save           :: itidana = -1 !parameter for calling tidal analysis

        double precision, allocatable :: vvk(:,:) !astronomical argument adjustment for phase
        double precision, allocatable :: uvk(:,:) !nodal modulation factor for phase
        double precision, allocatable :: fvk(:,:) !nodal modulation factor for amplitude

	integer, parameter	:: ntsemi=8 !number of semi-diurnal tidal const
	integer, parameter	:: ntdiur=7 !number of diurnal tidal const
	integer, parameter	:: ntlong=6 !number of long-term tidal const
        integer, parameter      :: ntd=ntsemi+ntdiur+ntlong !num of total tidal const

        real, parameter		:: rayl=0.95!Rayleigh criterion 

        type  :: const_entry
          character(5)          :: name     !constituent name
	  integer               :: nj       !number of satellites for the const.
          integer               :: dood(6)  !Doodson numbers
          integer	        :: mask	    !Mask: =1, tide included in tidef
          real		        :: amp	    !Tidal constituent amplitude [m]
          real		        :: loven    !Elasticy factor
          double precision      :: semi     !phase correction
          double precision      :: freq     !frequency (cy/day)
        end type const_entry
        type(const_entry), save, allocatable :: const_ar(:)

	double precision, private, allocatable :: astc(:,:) !Coefficients of the ephermides

	integer, private, save 	:: nsat	    !total number of satellites
        type  :: sat_entry
          sequence
	  integer          	:: ddel(3)  !changes in the last three Doodson numb
	  integer          	:: ir       !flag for amplitude corr
	  double precision 	:: ph       !phase correction
	  double precision 	:: ee       !amp ration of sat tidal pot to main const
        end type sat_entry
        type(sat_entry), private, save, allocatable :: sat_ar(:)

	integer, private, save 	:: nshal    !total number of shallow consts
        type  :: shall_entry
          character(5)     	:: name     !constituent name
	  double precision 	:: coef	    !phase correction
        end type shall_entry
        type(shall_entry), private, save, allocatable :: shall_ar(:)

	double precision        :: astr(6) !five ephermides + tau
	double precision        :: ader(6) !time derivative of astro + dtau

!==================================================================
        contains
!==================================================================
! Allocate tidal array

        subroutine tide_init(nkn)

        integer	:: nkn

        if( nkn == nkn_tide ) return

        if( nkn_tide > 0 ) then
            deallocate(zeqv)
            deallocate(vvk)
            deallocate(uvk)
            deallocate(fvk)
        end if

        nkn_tide = nkn

        if( nkn == 0 ) return

        allocate(zeqv(nkn))
        allocate(vvk(ntd,nkn))
        allocate(uvk(ntd,nkn))
        allocate(fvk(ntd,nkn))

        end subroutine tide_init

!********************************************************************
! Set tidal parameters

	subroutine tidepar_init

	implicit none

! intitalize coefficients of the ephermides
        double precision, parameter, dimension(4) :: hcp =
     +      (/279.696678, 0.9856473354, 0.00002267, 0.000000000/)
        double precision, parameter, dimension(4) :: ppc =
     +      (/281.220844, 0.0000470684, 0.0000339, 0.000000070/)
        double precision, parameter, dimension(4) :: scp =
     +      (/270.434164, 13.1763965268, -0.0000850, 0.000000039/)
        double precision, parameter, dimension(4) :: pcp =
     +      (/334.329556, 0.1114040803, -0.0007739, -0.00000026/)
        double precision, parameter, dimension(4) :: npc =
     +      (/-259.183275, 0.0529539222, -0.0001557, -0.000000050/)

	allocate(astc(4,5))

        astc(:,1) = hcp
        astc(:,2) = ppc
        astc(:,3) = scp
        astc(:,4) = pcp
        astc(:,5) = npc

! intitalize tidal constituent parameters
	allocate(const_ar(ntd))

        const_ar%name = (/
     +		'M2 ','S2 ','N2 ','K2 ','NU2','MU2','L2 ','T2 ',
     +    	'K1 ','O1 ','P1 ','Q1 ','J1 ','OO1','S1 ',
     +  	'MF ','MM ','MSM','MSF','SSA','SA '/)

        const_ar%freq = (/
     +		1.93227362d0, 2.00000000d0, 1.89598197d0, 2.00547582d0,
     +    	1.90083888d0, 1.86454723d0, 1.96856526d0, 1.99726222d0,
     +  	1.00273791d0, 0.92953571d0, 0.99726209d0, 0.89324406d0,
     +          1.03902956d0, 1.07594011d0, 1.00000013d0,
     +          0.07320220d0, 0.03629165d0, 0.06772638d0, 0.03143470d0,
     +          0.00547582d0, 0.00273778d0/)

        const_ar%nj = (/
     +		 9, 3, 4, 5, 4, 3, 5, 0,
     +          10, 8, 6,10,10, 8, 2,
     +		 0, 0, 0, 0, 0, 0/)

	!semi-diurnal (8)
        const_ar(1)%dood  = (/2, 0, 0, 0, 0, 0/)
        const_ar(2)%dood  = (/2, 2,-2, 0, 0, 0/)
        const_ar(3)%dood  = (/2,-1, 0, 1, 0, 0/)
        const_ar(4)%dood  = (/2, 2, 0, 0, 0, 0/)
        const_ar(5)%dood  = (/2,-1, 2,-1, 0, 0/)
        const_ar(6)%dood  = (/2,-2, 2, 0, 0, 0/)
        const_ar(7)%dood  = (/2, 1, 0,-1, 0, 0/)
        const_ar(8)%dood  = (/2, 2,-3, 0, 0, 1/)
	!diurnal (7)
        const_ar(9)%dood  = (/1, 1, 0, 0, 0, 0/)
        const_ar(10)%dood = (/1,-1, 0, 0, 0, 0/)
        const_ar(11)%dood = (/1, 1,-2, 0, 0, 0/)
        const_ar(12)%dood = (/1,-2, 0, 1, 0, 0/)
        const_ar(13)%dood = (/1, 2, 0,-1, 0, 0/)
        const_ar(14)%dood = (/1, 3, 0, 0, 0, 0/)
        const_ar(15)%dood = (/1, 1,-1, 0, 0, 1/)
	!long-period (6)
        const_ar(16)%dood = (/0, 2, 0, 0, 0, 0/)
        const_ar(17)%dood = (/0, 1, 0,-1, 0, 0/)
        const_ar(18)%dood = (/0, 1,-2, 1, 0, 0/)
        const_ar(19)%dood = (/0, 2,-2, 0, 0, 0/)
        const_ar(20)%dood = (/0, 0, 2, 0, 0, 0/)
        const_ar(21)%dood = (/0, 0, 1, 0, 0,-1/)

        const_ar%mask = (/1, 1, 1, 1, 1, 1, 1, 1,
     +                    1, 1, 1, 1, 1, 1, 1,
     +                    1, 1, 1, 1, 1, 1/)

!        const_ar%amp = (/
!     +         0.242334,0.112841,0.046398,0.030704,0.008877,0.007463,
!     +         0.006899,0.007563,
!     +         0.141465,0.100514,0.046843,0.019256,0.007965,0.004361,
!     +         0.001116,
!     +         0.042017,0.022191,0.004239,0.003678,0.019542,0.003104/)

!Kanta 1995, Table 6.1.1
        const_ar%amp = (/
     +	       0.244102,0.113572,0.046735,0.030875,0.008877,0.007463,
     +	       0.006899,0.006636,
     +         0.142408,0.101266,0.047129,0.019387,0.007965,0.004361,
     +	       0.001116,
     +         0.042017,0.022191,0.004239,0.003678,0.019542,0.003104/)

        const_ar%loven = (/
     +	       0.6930,0.6930,0.6930,0.6930,0.6930,0.6930,0.6930,0.6930,
     +         0.7364,0.6950,0.7059,0.6946,0.6911,0.6925,0.7126,
     +         0.6930,0.6930,0.6930,0.6930,0.6930,0.6930/)

        const_ar%semi = (/
     +         0.d0,0.d0,0d0,0.d0,0.d0,0.d0,-0.50d0,0.d0,
     +        -0.75d0,-0.25d0,-0.25d0,-0.25d0,-0.75d0,-0.75d0,-0.75d0,
     +         0.d0,0.d0,0d0,0.d0,0.d0,0.d0/)

! intitalize satellite parameters
	nsat = sum(const_ar%nj)
	allocate(sat_ar(nsat))

        sat_ar(:)%ddel(1) = (/
     +		-1,-1, 0, 0, 1, 1, 1, 2, 2,
     +		 0, 1, 2,
     +		-2,-1, 0, 0,
     +		-1,-1, 0, 0, 0,
     +		 0, 1, 2, 2,
     +          -1,-1, 0,
     +           0, 2, 2, 2, 2,
     +		-2,-1,-1,-1, 0, 0, 0, 0, 1, 1,
     +		-1, 0, 0, 1, 1, 1, 2, 2,
     +		 0, 0, 0, 1, 2, 2,
     +		-2,-2,-1,-1,-1, 0,-1, 0, 1, 2,
     +		 0, 0, 0, 1, 1, 1, 1, 2, 2, 2,
     +		-2,-2,-2,-1,-1, 0, 0, 0,
     +		 0, 0/)

        sat_ar(:)%ddel(2) = (/
     +          -1, 0,-2,-1,-1, 0, 1, 0, 1,
     +		-1, 0, 0,
     +          -2, 0,-2,-1,
     +	  	 0, 1,-1, 1, 2,
     +          -1, 0, 0, 1,
     +          -1, 0,-1,
     +          -1,-1, 0, 1, 2,
     +		-1,-1, 0, 1,-2,-1, 1, 2, 0, 1,
     + 		 0,-2,-1,-1, 0, 1, 0, 1,
     +		-2,-1, 0, 0, 0, 1,
     +		-3,-2,-2,-1, 0,-2, 0,-1, 0, 0,
     +		-1, 1, 2,-1, 0, 1, 2, 0, 1, 2,
     +          -1, 0, 1, 0, 1, 1, 2, 3,
     +		 0, 1/)

        sat_ar(:)%ddel(3) = (/
     +           0, 0, 0, 0, 0, 0, 0, 0, 0, 
     +		 0, 0, 0, 
     +           0, 1, 0, 0,
     +		 0, 0, 0, 0, 0,
     +           0, 0, 0, 0,
     +		 0, 0, 0, 
     +		 0, 0, 0, 0, 0,
     +		 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     +		 0, 0, 0, 0, 0, 0, 0, 0,
     +		 0, 0, 2, 0, 0, 0, 
     +		 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
     +		 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     +		 0, 0, 0, 0, 0, 0, 0, 0,
     +          -2, 0/)

        sat_ar%ir = (/
     +           2, 2, 0, 0, 2, 2, 2, 0, 0, 
     +		 0, 2, 0,
     +           0, 0, 0, 0, 
     +		 2, 2, 0, 0, 0,
     +		 0, 2, 0, 0,
     +		 2, 2, 0,
     +		 0, 0, 0, 0, 0,
     +		 0, 1, 1, 1, 0, 0, 0, 0, 1, 1,
     +		 1, 0, 0, 1, 1, 1, 0, 0, 
     +		 0, 0, 0, 1, 0, 0,
     +		 0, 0, 1, 1, 1, 0, 0, 0, 1, 0,
     +		 0, 0, 0, 1, 1, 1, 1, 0, 0, 0,
     +		 0, 0, 0, 1, 1, 0, 0, 0,
     +		 0, 0/)

        sat_ar%ph = (/
     +      0.75d0,0.75d0,0.00d0,0.50d0,0.25d0,0.75d0,0.75d0,0.00d0,
     +      0.00d0,
     +      0.00d0,0.75d0,0.00d0,
     +      0.50d0,0.00d0,0.00d0,0.50d0,
     +      0.75d0,0.75d0,0.50d0,0.00d0,0.00d0,
     +	    0.50d0,0.75d0,0.00d0,0.50d0,
     +	    0.25d0,0.25d0,0.50d0,
     +	    0.50d0,0.00d0,0.50d0,0.50d0,0.50d0,
     +      0.00d0,0.75d0,0.25d0,0.75d0,0.00d0,0.50d0,0.00d0,0.50d0,
     +	    0.25d0,0.25d0,
     +      0.25d0,0.50d0,0.00d0,0.25d0,0.75d0,0.25d0,0.50d0,0.50d0,
     +      0.00d0,0.50d0,0.50d0,0.75d0,0.50d0,0.50d0,
     +      0.50d0,0.50d0,0.75d0,0.75d0,0.75d0,0.50d0,0.00d0,0.00d0,
     +      0.75d0,0.50d0,
     +	    0.50d0,0.00d0,0.50d0,0.75d0,0.25d0,0.25d0,0.25d0,0.50d0,
     +	    0.50d0,0.50d0,
     +      0.50d0,0.00d0,0.00d0,0.25d0,0.25d0,0.00d0,0.00d0,0.00d0,
     +	    0.00d0,0.50d0/)

        sat_ar%ee = (/
     +      1.d-4,4.d-4,5.d-4,3.73d-2,1.0d-4,9.d-4,2.d-4,6.d-4,2.d-4,
     +	    2.2d-3,1.d-4,1.d-4,
     +      3.9d-3,8.d-4,5.d-4,3.73d-2,
     +      2.4d-3,4.d-4,1.28d-2,0.298d0,3.24d-2,
     +	    3.73d-2,4.2d-3,4.2d-3,3.6d-3,
     +	    1.8d-3,1.04d-2,3.75d-2,
     +	    3.66d-2,4.7d-3,2.505d-1,1.102d-1,1.56d-2,
     +      2.d-4,1.d-4,7.d-4,1.d-4,1.d-4,1.98d-2,0.1356d0,2.9d-3,
     +      2.d-4,1.d-4,
     +      3.d-4,5.8d-3,0.1885d0,4.d-4,2.9d-3,4.d-4,6.4d-3,1.d-3,
     +      8.d-4,1.12d-2,4.d-4,4.d-4,1.5d-3,3.d-4,
     +      7.d-4,3.9d-3,1.d-3,1.15d-2,2.92d-2,5.7d-3,8.d-4,0.1884d0,
     +      1.8d-3,2.8d-3,
     +	    2.94d-2,1.98d-1,4.7d-3,2.7d-3,8.16d-2,3.31d-2,2.7d-3,
     +      1.52d-2,9.8d-3,5.7d-3,
     +	    3.7d-3,0.1496d0,2.96d-2,2.4d-2,9.9d-3,0.6398d0,0.1342d0,
     +      8.6d-3,
     +	    3.534d-1,2.64d-2/)

! shallow constituents - NOT YET IMPLEMENTED
        nshal = 0
	allocate(shall_ar(nshal))

        end subroutine tidepar_init

!**********************************************************************
! This subroutine calculates the five ephermides of the sun and moon
! and their derivative. Units for the ephermides are cycles and for 
! their derivatives are cycles/365 days.
!   astr(1) = h	  !mean longitude of the sum
!   astr(2) = pp  !mean longitude of the solar perigee
!   astr(3) = s	  !mean longitude of the moon
!   astr(4) = p	  !mean longitude of the lunar perigee
!   astr(5) = np  !neg. of the long. of the mean ascending node
!   astr(6) = tau !lunar time tau, based on fractional part of solar day
!   ader(:) = ... !time derivative of astr

        subroutine astro(d)

        implicit none

        double precision, intent(in)  :: d       !Universal time (days since 1899-12-31::12)

        double precision, parameter    :: f=360.d0
	double precision	       :: d2,d22,d23,f2
	integer			       :: i

        d2  = d*1.d-4
        d22 = d2*d2
        d23 = d2**3
        f2  = f/365.d0

	do i = 1,5
          astr(i) = astc(1,i)+astc(2,i)*d+astc(3,i)*d22+astc(4,i)*d23
          astr(i) = dmod(astr(i)/f, 1.d0)
          ader(i) = astc(2,i)+2.d-4*astc(3,i)*d2+3.d-4*astc(4,i)*d22
          ader(i) = dmod(ader(i)/f2, 1.d0)
        end do

        astr(6) = dmod((d+0.5d0),1.d0) + astr(1) - astr(3)
        ader(6) = 1.d0 + ader(1) - ader(3)

        return

        end subroutine astro

!**********************************************************************
! setvuf calculates the V,U,F values at time d for all constituents
!  ntd is the number of main constituents
!  ntotal is the number of constituents (main + shallow water)
!  F is the nodal modulation adjustment factor for amplitude
!  U is the nodal modulation adjustment factor for phase
!  V is the astronomical argument adjustment for phase.

        subroutine setvuf(lat,v,u,f)

        implicit none

        double precision, intent(in)  :: lat    !latitude
        double precision, intent(out) :: v(ntd) !ASTRONOMICAL ARGUMENT ADJUSTMENT FOR PHASE
        double precision, intent(out) :: u(ntd) !NODAL MODULATION FACTOR FOR PHASE
        double precision, intent(out) :: f(ntd) !NODAL MODULATION FACTOR FOR AMPLITUDE
  
        integer 		      :: ntotal
        integer 		      :: jbase,k,k1,km1,l,lk,j,j1,jl
        double precision 	      :: vdbl
        double precision 	      :: slat,slat2,sumc,sums,rr,uudbl
        double precision, parameter   :: pi=4.d0*atan(1.d0)
        double precision, parameter   :: twopi=2.d0*pi
        double precision, parameter   :: eps = 1.d-20

        slat = dsin(pi*lat/180.)
        slat2 = slat*slat

! ONLY THE FRACTIONAL PART OF A SOLAR DAY NEED BE RETAINED FOR COMPUTING
! THE LUNAR TIME tau.

        jbase = 0
        do k = 1,ntd
          vdbl=const_ar(k)%dood(1)*astr(6)+const_ar(k)%dood(2)*astr(3)+
     +         const_ar(k)%dood(3)*astr(1)+const_ar(k)%dood(4)*astr(4)+ 
     +         const_ar(k)%dood(5)*astr(5)+const_ar(k)%dood(6)*astr(2)+ 
     +    const_ar(k)%semi
          vdbl = mod(vdbl,1.D+0)
          j1 = jbase+1
          jl = jbase + const_ar(k)%nj
  
! SATELLITE AMPLITUDE RATIO ADJUSTMENT FOR LATITUDE
          sumc = 1.
          sums = 0.
          do j = j1,jl
            rr = sat_ar(j)%ee
            l  = sat_ar(j)%ir
  	    select case (l)
  	      case (1)	!latitude correction for diurnal constituents
                rr = sat_ar(j)%ee*0.36309*(1.-5.*slat2)/(slat+eps)
  	      case (2)	!latitude correction for semi-diurnal constituents
                rr = sat_ar(j)%ee*2.59808*slat
  	    end select
            uudbl=sat_ar(j)%ddel(1)*astr(4)+sat_ar(j)%ddel(2)*astr(5)+ 
     +            sat_ar(j)%ddel(3)*astr(2)+sat_ar(j)%ph 
  	    uudbl = dmod(uudbl,1.d0) * twopi 
            sumc = sumc + rr*cos(uudbl)
            sums = sums + rr*sin(uudbl)
      	  end do
          f(k) = sqrt(sumc*sumc + sums*sums)
          v(k) = vdbl
          u(k) = atan2(sums,sumc)/twopi
          jbase = jl
        end do

!  F AND V+U OF THE SHALLOW WATER CONSTITUENTS ARE COMPUTED FROM
!  THE VALUES OF THE MAIN CONSTITUENT FROM WHICH THEY ARE DERIVED.
!!!!!  SHALLOW WATER CONSTITUENTS NOT YET IMPLEMENTED !!!!!!

!        jbase = 0
!        ntotal = ntd + nshal
!        k1 = ntd+1
!        if(k1 > ntotal) return
!
!        do k = k1,ntotal
!          f(k) = 1.0
!          v(k) = 0.0
!          u(k) = 0.
!          do lk = 1,mf
!            if(const_ar(k)%name == name(lk)) go to 268
!          end do
!268       j1=jbase+1
!          jl=jbase+const_ar(k)%nj
!          do j=j1,jl
!            km1 = k-1
!            do l = 1,km1
!              if(const_ar(l)%name /= shall_ar(j)%name) then
!                write(6,*)shall_ar(j)%name
!                stop 'error stop SETVUF shallow'
!              end if
!            end do
!            f(k) = f(k)*f(l)**abs(shall_ar(j)%coef)
!            v(k) = v(k) + shall_ar(j)%coef*v(l)
!            u(k) = u(k) + shall_ar(j)%coef*u(l)
!          end do
!          jbase = jl
!        end do
  
        return

        end subroutine setvuf

!**********************************************************************
! Check Rayleigh criterion. 
! In case of violation set mask = 0 overwritting previous setting

        subroutine rayleigh_crit(dtday)

        implicit none

        real, intent(in)                :: dtday  !time series lenght [day]
        real                            :: minres !min freq [ch/day]

        minres = rayl / dtday

        where( const_ar%freq < minres ) const_ar%mask=const_ar%mask*0

        end subroutine rayleigh_crit

!==================================================================
        end module tide
!==================================================================#

!**********************************************************************
! Computete Universal time in days from absolute time

        subroutine tide_time(atime,d)

        implicit none

        double precision, intent(in)  :: atime !absolute time [s]

        double precision, intent(out) :: d     !Universal time [days since 1899-12-31::12]

        double precision              :: secs  !seconds in the day [s]
        integer                       :: jd    !julian day [d]
        integer                       :: year,month,day
        integer                       :: days
        double precision, parameter   :: secs_in_day = 86400.d0

        days = atime/secs_in_day
        secs = atime - secs_in_day*days
        call days2date(days,year,month,day)
        call date2j(year,month,day,jd)

        d = jd + 365.d0*(year - 1900.d0) + int((year-1901.d0)/4.d0)
        d = d - 0.5d0 + secs/secs_in_day

        end subroutine tide_time

!********************************************************************
! Performs tidal analysis and compute amplitude and phase of selected
! tidal constants

        subroutine tda_solve(nequ,ndat,var,acov,bcov,h,g)

        use tide, only : ntd

        implicit none

!  input variables
        integer, intent(in)           :: nequ
	integer, intent(in)           :: ndat
	double precision, intent(in)  :: var
        double precision, intent(in)  :: acov(nequ,nequ)
        double precision, intent(in)  :: bcov(nequ)

!  output variables
        real, intent(out)	      :: h(ntd)	!amplitude [m]
        real, intent(out)	      :: g(ntd)	!phase [deg]

!  local variables
	integer			      :: nequm1
	integer 		      :: icon,k,iequ1,iequ2,ic,is
	integer			      :: job,info
	double precision	      :: arg
	double precision 	      :: dndat,lvar
	double precision, allocatable :: a(:,:),acovl(:,:)
	double precision, allocatable :: usvd(:,:)
	double precision, allocatable :: vsvd(:,:)
	double precision, allocatable :: wsvd(:)
	double precision, allocatable :: e(:)
	double precision, allocatable :: tmpsvd(:)
	double precision, allocatable :: b(:)
	double precision, allocatable :: b1(:)
	double precision, allocatable :: c1(:)
	double precision, allocatable :: c(:)
        double precision, parameter   :: d2r = 4.0*atan(1.0)/180.0  !Degrees to radians

	allocate(a(nequ,nequ))
	allocate(acovl(nequ,nequ))
	allocate(usvd(nequ,nequ))
	allocate(vsvd(nequ,nequ))
	allocate(wsvd(nequ))
	allocate(e(nequ))
	allocate(tmpsvd(nequ))
	allocate(b(nequ))
	allocate(b1(nequ))
	allocate(c1(nequ))
	allocate(c(nequ))

	dndat  = dble(ndat)
        nequm1 = nequ-1

	acovl = acov

!  Divide VAR, ACOV and BCOV by number of data points and
!  generate A and B:
        lvar = var/dndat
        a = acovl/dndat
        b = bcov/dndat

!  Solve A*C=B, using SVD:
        job=11
        call dsvdc(a,nequ,nequ,nequ,wsvd,e,usvd,nequ,vsvd,nequ,
     +           tmpsvd,job,info)

!  Back substitute:
        call dsvbksb2(usvd,wsvd,vsvd,nequ,nequ,nequ,
     +                nequ,nequ,b,c,tmpsvd)

!  Compute h and g
        do icon=1,ntd
          ic = 2*icon
          is = ic+1
          h(icon) = sqrt(c(ic)**2+c(is)**2)
          if(c(is).ne.0.0.or.c(ic).ne.0.0) then
            arg = datan2(c(is),c(ic))/d2r
            if(arg.lt.0.0) arg = arg + 360.0
          else
            arg = 0.0
          endif
          g(icon) = arg
        end do

	end subroutine tda_solve

! ******************************************************************************
! Singular value decomposition solution

      subroutine dsvbksb2(u,w,v,n,p,minnp,nmax,pmax,b,x,tmp)

      implicit none

      integer		:: n,p,minnp,nmax,pmax
      double precision	:: u(nmax,minnp),w(minnp),v(pmax,minnp)
      double precision	:: b(nmax),x(minnp),tmp(minnp)

      double precision  :: s
      integer           :: j

      do j=1,minnp
        s=0.d0
        if(w(j).ne.0.)then
          s = sum(u(:,j)*b(:))/w(j)
        endif
        tmp(j)=s
      end do

      do j=1,p
        x(j) = sum(v(j,:)*tmp(:))
      end do

      return

      end subroutine dsvbksb2

!*********************************************************************
! Computes overdetermined array

        subroutine get_overarray(nequ,v,u,f,x)

        use tide

        implicit none

!  input variables
        integer, intent(in)           :: nequ
        double precision, intent(in)  :: v(ntd)  !astronomical argument adjustment for phase
        double precision, intent(in)  :: u(ntd)  !nodal modulation factor for phase
        double precision, intent(in)  :: f(ntd)  !nodal modulation factor for amplitude

!  output variables
        double precision, intent(out) :: x(nequ) !overdetermined array

!  local variables
        double precision              :: arg
        double precision              :: ux,vx,fx
        integer                       :: icon,ic,is
        double precision, parameter   :: twopi=8.d0*atan(1.d0)

        x(1) = 1.d0
        do icon = 1,ntd
          ic = 2*icon
          is = ic+1
          ux = u(icon)
          vx = v(icon)
          fx = f(icon)
          arg = (vx + ux)*twopi
          x(ic) = fx*dcos(arg)
          x(is) = fx*dsin(arg)
        end do

        end subroutine get_overarray

!**********************************************************************
