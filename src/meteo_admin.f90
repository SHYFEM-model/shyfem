!
! $Id$
!
! meteo files
!
! revision log :
!
! 16.02.2011    ggu     created by copying mainly from subn11.f
! 30.05.2014    ggu     new imreg == 3
! 10.07.2014    ggu     only new file format allowed
! 30.04.2015    ggu     ice cover implemented
!
!*********************************************************************
!---------------------------------------------------------------------
        module meteo_admin
!---------------------------------------------------------------------
        contains
!---------------------------------------------------------------------

	subroutine convert_distributed

! converts distributed source from [m/s] to [m**3/s]

        use bnd_dynamic
        use evgeom
        use basin
        use shympi

        implicit none


        include 'param.h'

        integer k,ie,ii
        double precision area3
#ifdef DEBUGON
        integer e
        double precision v1v(nkn_local)
#else
        double precision v1v(nkn)
#endif

        do k=1,nkn
          v1v(k) = 0.d0
        end do

#ifdef DEBUGON
	do e=1,nel_local
          if(bmpi) then
            ie=domain%elems%mapID(e)
          else
            ie=e
          endif
#else
	do ie=1,nel
#endif
          area3 = 4. * ev(10,ie)
          do ii=1,3
            k = nen3v(ii,ie)
            v1v(k) = v1v(k) + area3
          end do
        end do

#ifndef DEBUGON
        call shympi_exchange_and_sum_2d_nodes(v1v)      !gmica
#endif

        do k=1,nkn
          rqdsv(k) = rqdsv(k) * v1v(k)
        end do

        end

!*******************************************************************
!*******************************************************************
!*******************************************************************

        subroutine evap_init

! initializes evaporation mass flux

        use meteo
        use basin, only : nkn,nel,ngr,mbw

        implicit none

        include 'param.h'


        integer k

        do k=1,nkn
          evapv(k) = 0.d0
        end do

        end

!*******************************************************************

        subroutine rain_evap_set

! adds evaporation mass flux to distributed source

        use para
        use meteo
        use bnd_dynamic
        use basin, only : nkn,nel,ngr,mbw

        implicit none

        include 'param.h'
        include 'femtime.h'

        integer k,ievap

        ievap = nint(getpar('ievap'))

        do k=1,nkn
          rqdsv(k) = rqdsv(k) + metrain(k) - ievap*evapv(k)
        end do

        end

!*******************************************************************
!*******************************************************************
!*******************************************************************

        subroutine windcd_init

! initializes evaporation mass flux

        use para
        use meteo
        use basin, only : nkn,nel,ngr,mbw

        implicit none

        include 'param.h'


        integer k
        double precision dragco

        dragco = getpar('dragco')

        do k=1,nkn
          windcd(k) = dragco
        end do

        end

!*******************************************************************
!*******************************************************************
!*******************************************************************

        subroutine meteo_init

! initializes meteo variables

        use meteo_forcing

        implicit none

        write(6,*) 'initializing meteo'

        call evap_init
        call windcd_init        !windcd(:) = dragco

        call meteo_forcing_fem  !new file format

        end

!*******************************************************************

        subroutine meteo_force

! update meteo variables and admin rain/evaporation

        use meteo_forcing

        implicit none

        call meteo_forcing_fem  !new file format

	!call compute_heat_flux

        call rain_evap_set              !add evaporation
        call convert_distributed        !convert from [m/s] to [m**3/s]

        end

!*******************************************************************
!*******************************************************************
!*******************************************************************

	subroutine get_meteo_forcing(k,wx,wy,tauxn,tauyn,p)

! returns wind (wx/y), normalized stress (taux/yn) and pressure (p)

        use meteo

        implicit none

        integer k		!node number
        double precision wx,wy		!wind velocity [m/s]
        double precision tauxn,tauyn	!normalized stress [m**2/s**2]
        double precision p			!pressure [Pa]

        include 'param.h'

        wx = wxv(k)
        wy = wyv(k)
        tauxn = tauxnv(k)
        tauyn = tauynv(k)
        p = ppv(k)

        end

!*******************************************************************

	subroutine get_light(k,rad_light)

! returns light intensity [W/m**2]

        use meteo_forcing

	implicit none

        integer k               !node number
        double precision rad_light          !watt/m**2

	call meteo_get_solar_radiation(k,rad_light)

        end

!*******************************************************************

	subroutine get_ice(k,ice_cover)

! returns ice cover [fraction 0-1]

	use meteo

	implicit none

	include 'param.h'

        integer k               !node number
        double precision ice_cover          ![0-1]

	ice_cover = metice(k)

        end

!*******************************************************************

	subroutine get_ice_all(ice_cover)

! returns ice cover array [fraction 0-1]

	use meteo
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'

        double precision ice_cover(nkn)          ![0-1]

	integer k

	do k=1,nkn
	  ice_cover(k) = metice(k)
	end do

        end

!*******************************************************************

!---------------------------------------------------------------------
        end module meteo_admin
!---------------------------------------------------------------------
