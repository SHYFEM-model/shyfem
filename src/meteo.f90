
        module meteo

        implicit none

! metrain and evapv are in [m/s]
! metrain is read from file in [mm/day] and converted to [m/s]


        integer, private, save  :: nkn_meteo = 0

        double precision, allocatable, save :: wxv(:)       ! wind velocity in x [m/s]
        double precision, allocatable, save :: wyv(:)       ! wind velocity in y [m/s]
        double precision, allocatable, save :: ppv(:)       ! pressure (atmos) [Pa,mbar]

        double precision, allocatable, save :: tauxnv(:)    ! wind stress in x [N/m**2]
        double precision, allocatable, save :: tauynv(:)    ! wind stress in y [N/m**2]

        double precision, allocatable, save :: metrad(:)    ! downward sw solar rad [W/m**2]
        double precision, allocatable, save :: methum(:)    ! humidity [%]
        double precision, allocatable, save :: metdew(:)    ! dew point temperature [C]  
        double precision, allocatable, save :: mettair(:)   ! 10 m air temperature [C]
        double precision, allocatable, save :: metcc(:)     ! cloud cover [0-1]
        double precision, allocatable, save :: metrain(:)   ! precipitation [m/s]
        double precision, allocatable, save :: metwbt(:)    ! wet bulb temperature [C] 

        double precision, allocatable, save :: metice(:)    ! ice cover [0-1]
        double precision, allocatable, save :: metws(:)     ! wind speed [m/s]

        double precision, allocatable, save :: windcd(:)    ! wind drag coefficient
        double precision, allocatable, save :: evapv(:)     ! evaporation [m/s]

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

        end module meteo


