
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

! revision log :
!
! 10.07.2015	ggu	changed VERS_7_1_50
! 28.04.2016	ggu	changed VERS_7_5_9
! 09.09.2016	ggu	changed VERS_7_5_17
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62

!**************************************************************************

        module mod_meteo

        implicit none

c metrain and evapv are in [m/s]
c metrain is read from file in [mm/day] and converted to [m/s]


        integer, private, save  :: nkn_meteo = 0

        real, allocatable, save :: wxv(:)	! wind velocity in x [m/s]
        real, allocatable, save :: wyv(:)	! wind velocity in y [m/s]
        real, allocatable, save :: ppv(:)	! pressure (atmos) [Pa,mbar]

        real, allocatable, save :: tauxnv(:)	! wind stress in x [N/m**2]
        real, allocatable, save :: tauynv(:)	! wind stress in y [N/m**2]

        real, allocatable, save :: metrad(:)	! downward sw solar rad [W/m**2]
        real, allocatable, save :: methum(:)	! humidity [%]
        real, allocatable, save :: metdew(:)    ! dew point temperature [C]  
        real, allocatable, save :: mettair(:)	! 10 m air temperature [C]
        real, allocatable, save :: metcc(:)	! cloud cover [0-1]
        real, allocatable, save :: metrain(:)	! precipitation [m/s]
        real, allocatable, save :: metwbt(:)	! wet bulb temperature [C] 

        real, allocatable, save :: metice(:)	! ice cover [0-1]
        real, allocatable, save :: metws(:)	! wind speed [m/s]

        real, allocatable, save :: windcd(:)	! wind drag coefficient
        real, allocatable, save :: evapv(:)	! evaporation [m/s]

        contains

!************************************************************

        subroutine mod_meteo_init(nkn)

        integer  :: nkn

        if( nkn == nkn_meteo ) return

        if( nkn_meteo > 0 ) then
          deallocate(wxv)
          deallocate(wyv)
          deallocate(ppv)
          deallocate(tauxnv)
          deallocate(tauynv)
          deallocate(metrad)
          deallocate(methum)
          deallocate(metdew)  
          deallocate(mettair)
          deallocate(metcc)
          deallocate(metrain)
          deallocate(metwbt)
          deallocate(metice)
          deallocate(metws)
          deallocate(windcd)
          deallocate(evapv)
        end if

        nkn_meteo = nkn

        if( nkn == 0 ) return

        allocate(wxv(nkn))
        allocate(wyv(nkn))
        allocate(ppv(nkn))
        allocate(tauxnv(nkn))
        allocate(tauynv(nkn))
        allocate(metrad(nkn))
        allocate(methum(nkn))
        allocate(metdew(nkn))  
        allocate(mettair(nkn))
        allocate(metcc(nkn))
        allocate(metrain(nkn))
        allocate(metwbt(nkn))
        allocate(metice(nkn))
        allocate(metws(nkn))
        allocate(windcd(nkn))
        allocate(evapv(nkn))

        end subroutine mod_meteo_init

!************************************************************

        end module mod_meteo


