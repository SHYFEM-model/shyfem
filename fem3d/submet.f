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
c
c*********************************************************************

	subroutine convert_distributed

c converts distributed source from [m/s] to [m**3/s]

	implicit none

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	real rqpsv(1), rqdsv(1)
	common /rqpsv/rqpsv, /rqdsv/rqdsv
	integer nen3v(3,1)
	common /nen3v/nen3v
	include 'ev.h'
	real v1v(1)
	common /v1v/v1v

	integer k,ie,ii
	real area3

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

	implicit none

	include 'meteo.h'

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	integer k

	do k=1,nkn
	  evapv(k) = 0.	
	end do

	end

c*******************************************************************

	subroutine rain_evap_set

c adds evaporation mass flux to distributed source

	implicit none

	include 'meteo.h'

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real rqdsv(1)
	common /rqdsv/rqdsv

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

	implicit none

	include 'meteo.h'

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

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

	implicit none

	include 'param.h'

	integer itanf,itend,idt,nits,niter,it
	common /femtim/ itanf,itend,idt,nits,niter,it
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real tempv(nlvdim,nkndim)
	common /tempv/tempv
	double precision dq
	real dt

        call get_timestep(dt)
        call qflux3d(it,dt,nkn,nlvdim,tempv,dq)	!compute heat flux

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine get_meteo_forcing(k,wx,wy,tauxn,tauyn,p)

c returns wind (wx/y), normalized stress (taux/yn) and pressure (p)

	implicit none

	integer k		!node number
	real wx,wy		!wind velocity [m/s]
	real tauxn,tauyn	!normalized stress [m**2/s**2]
	real p			!pressure [Pa]

	include 'meteo.h'

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

