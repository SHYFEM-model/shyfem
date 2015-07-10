
        module mod_meteo

        implicit none

        !real wxv(nkndim),wyv(nkndim)
        !common /wxv/wxv,/wyv/wyv
        !real ppv(nkndim)
        !common /ppv/ppv

        !real metrad(nkndim),methum(nkndim)
        !real mettair(nkndim),metcc(nkndim)
        !common /metrad/metrad, /methum/methum
        !common /mettair/mettair, /metcc/metcc

        !real metrain(nkndim)
        !common /metrain/metrain

        !real tauxnv(nkndim),tauynv(nkndim)
        !common /tauxnv/tauxnv,/tauynv/tauynv

        !real metwbt(nkndim),metws(nkndim)
        !common /metwbt/metwbt, /metws/metws

        !real windcd(nkndim)          !wave drag coefficient
        !common /windcd/windcd

	!real evapv(nkndim)		!evaporation
	!common /evapv/evapv

c metrain and evapv are in [m/s]
c metrain is read from file in [mm/day] and converted to [m/s]


        integer, private, save  :: nkn_meteo = 0

        real, allocatable, save :: wxv(:)	! wind velocity in x [m/s]
        real, allocatable, save :: wyv(:)	! wind velocity in y [m/s]
        real, allocatable, save :: ppv(:)	! pressure (atmospheric) [Pa,mbar]

        real, allocatable, save :: tauxnv(:)	! wind stress in x [N/m**2]
        real, allocatable, save :: tauynv(:)	! wind stress in y [N/m**2]

        real, allocatable, save :: metrad(:)	! downward shortwave solar radiation [W/m**2]
        real, allocatable, save :: methum(:)	! humidity [%]
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


