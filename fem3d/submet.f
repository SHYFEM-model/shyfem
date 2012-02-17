c
c $Id$
c
c meteo files
c
c revision log :
c
c 16.02.2011    ggu     created by copying mainly from subn11.f
c
c*********************************************************************

	subroutine rain0d( it )

c simulates rain (0D) increasing the water level

	implicit none

	integer it

	integer ndim
	parameter (ndim=100)
	real barray(ndim)
	save barray

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real metrain(1)
	common /metrain/metrain

	character*(80)	file
	integer k
	real rw
	real t,zdist

	real getpar

	integer icall
	save icall
	data icall / 0 /

	if( icall .eq. -1 ) return

c---------------------------------------------------------------
c initialization
c---------------------------------------------------------------

	if( icall .eq. 0 ) then		!initialize
	   
	   zdist = getpar('zdist')
	   call getfnm('rain',file)

	   if( zdist .eq. 0. .and. file .eq. ' ' ) icall = -1
	   if( icall .eq. -1 ) return

	   call exffild(file,ndim,barray,zdist)

	   icall = 1

	end if

c---------------------------------------------------------------
c normal call
c---------------------------------------------------------------

	t = it
	call exfintp(barray,t,rw)

	do k=1,nkn
	  metrain(k) = rw
	end do

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

	end

c*******************************************************************

	subroutine rain2distributed

c adjustes distributed source with rain

	implicit none

	integer it

	real zconv
	parameter( zconv = 1.e-3 / 86400. )

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real metrain(1)
	common /metrain/metrain
	real rqdsv(1)
	common /rqdsv/rqdsv

	integer k

c---------------------------------------------------------------
c rain is in mm/day -> convert to m/s (water level change)
c---------------------------------------------------------------

	do k=1,nkn
	  rqdsv(k) = rqdsv(k) + metrain(k) * zconv
	end do

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

	end

c*******************************************************************

	subroutine convert_distributed

c converts distributed source from [m/s] to [m**3/s]

	implicit none

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	real rqpsv(1), rqdsv(1)
	common /rqpsv/rqpsv, /rqdsv/rqdsv
	real evapv(1)
	common /evapv/evapv
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
	  evapv(k) = evapv(k) * v1v(k)		!only for output purpose
	end do

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine evap_init

c initializes evaporation mass flux

	implicit none

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real evapv(1)
	common /evapv/evapv

	integer k

	do k=1,nkn
	  evapv(k) = 0.	
	end do

	end

c*******************************************************************

	subroutine evap_set

c adds evaporation mass flux to distributed source

	implicit none

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real evapv(1)
	common /evapv/evapv
	real rqdsv(1)
	common /rqdsv/rqdsv

	integer k,ievap
	real getpar

	ievap = nint(getpar('ievap'))
	if( ievap .le. 0 ) return

	do k=1,nkn
	  rqdsv(k) = rqdsv(k) - evapv(k)
	end do

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine meteo_init

c initializes meteo variables

	integer imreg
	real getpar

	imreg = nint(getpar('imreg'))

	if( imreg .eq. 1 ) then
	  call meteo_regular
	else
	  call windad
	  call qflux_init
	end if

	call evap_init

	end

c*******************************************************************

	subroutine meteo_force

c update meteo variables and admin rain/evaporation

	integer itanf,itend,idt,nits,niter,it
	common /femtim/ itanf,itend,idt,nits,niter,it

	integer imreg
	real getpar

	imreg = nint(getpar('imreg'))

	if( imreg .eq. 1 ) then
	  call meteo_regular
	else
	  call windad			!wind
	  call qflux_read(it)		!heat flux
	  call rain0d(it)			!rain
	end if

	call rain2distributed		!copy rain to distributed source
	call evap_set			!add evaporation
	call convert_distributed	!convert from [m/s] to [m**3/s]

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine get_meteo_forcing(k,wx,wy,tauxn,tauyn,p)

c returns wind (wx/y), normalized stress (taux/yn) and pressure (p)

	integer k		!node number
	real wx,wy		!wind velocity [m/s]
	real tauxn,tauyn	!normalized stress [m**2/s**2]
	real p			!pressure [Pa]

        real tauxnv(1),tauynv(1)
        common /tauxnv/tauxnv,/tauynv/tauynv
        real wxv(1),wyv(1)
        common /wxv/wxv,/wyv/wyv
        real ppv(1)
        common /ppv/ppv

	wx = wxv(k)
	wy = wyv(k)
	tauxn = tauxnv(k)
	tauyn = tauynv(k)
	p = ppv(k)

	end

c*******************************************************************

	subroutine get_light(k,rad_light)

c returns light intensity [W/m**2]

        integer k               !node number
        real rad_light          !watt/m**2

	call meteo_get_solar_radiation(k,rad_light)

        end

c*******************************************************************

