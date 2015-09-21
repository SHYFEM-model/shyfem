c
c $Id$
c
c meteo files
c
c revision log :
c
c 16.02.2011    ggu     created by copying mainly from subn11.f
c 30.05.2014    ggu     new imreg == 3
c 10.07.2014    ggu     only new file format allowed
c 30.04.2015    ggu     ice cover implemented
c
c*********************************************************************

	subroutine convert_distributed

c converts distributed source from [m/s] to [m**3/s]

	use mod_bound_dynamic
	use evgeom
	use basin

	implicit none


	include 'param.h'

	integer k,ie,ii
	real area3
	real v1v(nkn)

	do k=1,nkn
	  v1v(k) = 0.
	end do

	do ie=1,nel
	  area3 = 4. * ev(10,ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    v1v(k) = v1v(k) + area3
	  end do
	end do

	do k=1,nkn
	  rqdsv(k) = rqdsv(k) * v1v(k)
	end do

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine evap_init

c initializes evaporation mass flux

	use mod_meteo
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'


	integer k

	do k=1,nkn
	  evapv(k) = 0.	
	end do

	end

c*******************************************************************

	subroutine rain_evap_set

c adds evaporation mass flux to distributed source

	use mod_meteo
	use mod_bound_dynamic
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'


	integer k,ievap
	real getpar

	ievap = nint(getpar('ievap'))

	do k=1,nkn
	  rqdsv(k) = rqdsv(k) + metrain(k) - ievap*evapv(k)
	end do

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine windcd_init

c initializes evaporation mass flux

	use mod_meteo
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'


	integer k
	real dragco

	real getpar

	dragco = getpar('dragco')

	do k=1,nkn
	  windcd(k) = dragco
	end do

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine meteo_init

c initializes meteo variables

	use meteo_forcing_module

	implicit none

	write(6,*) 'initializing meteo'

	call evap_init
	call windcd_init

	call meteo_forcing_fem	!new file format

	end

c*******************************************************************

	subroutine meteo_force

c update meteo variables and admin rain/evaporation

	use meteo_forcing_module

	implicit none

	call meteo_forcing_fem	!new file format

	!call compute_heat_flux

	call rain_evap_set		!add evaporation
	call convert_distributed	!convert from [m/s] to [m**3/s]

	end

c*******************************************************************

	subroutine compute_heat_flux

c computes heat flux through bulk formulas

	use mod_ts
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'

	include 'femtime.h'

	double precision dq
	real dt

        call get_timestep(dt)
        call qflux3d(it,dt,nkn,nlvdi,tempv,dq)	!compute heat flux

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine get_meteo_forcing(k,wx,wy,tauxn,tauyn,p)

c returns wind (wx/y), normalized stress (taux/yn) and pressure (p)

	use mod_meteo

	implicit none

	integer k		!node number
	real wx,wy		!wind velocity [m/s]
	real tauxn,tauyn	!normalized stress [m**2/s**2]
	real p			!pressure [Pa]

	include 'param.h'

	wx = wxv(k)
	wy = wyv(k)
	tauxn = tauxnv(k)
	tauyn = tauynv(k)
	p = ppv(k)

	end

c*******************************************************************

	subroutine get_light(k,rad_light)

c returns light intensity [W/m**2]

	implicit none

        integer k               !node number
        real rad_light          !watt/m**2

	call meteo_get_solar_radiation(k,rad_light)

        end

c*******************************************************************

	subroutine get_ice(k,ice_cover)

c returns ice cover [fraction 0-1]

	use mod_meteo

	implicit none

	include 'param.h'

        integer k               !node number
        real ice_cover          ![0-1]

	ice_cover = metice(k)

        end

c*******************************************************************

	subroutine get_ice_all(ice_cover)

c returns ice cover array [fraction 0-1]

	use mod_meteo
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'

        real ice_cover(nkn)          ![0-1]

	integer k

	do k=1,nkn
	  ice_cover(k) = metice(k)
	end do

        end

c*******************************************************************

