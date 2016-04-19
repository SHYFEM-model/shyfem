
	module mod_hydro_baro

	implicit none

	integer, private, save :: nel_hydro_baro = 0

	real, allocatable, save :: uov(:), vov(:)
	real, allocatable, save :: unv(:), vnv(:)

	contains

!************************************************************

        subroutine mod_hydro_baro_init(nel)

        integer nel

        if( nel == nel_hydro_baro ) return

        if( nel_hydro_baro > 0 ) then
          deallocate(uov)
          deallocate(vov)

          deallocate(unv)
          deallocate(vnv)
        end if

        nel_hydro_baro = nel

        if( nel == 0 ) return

        allocate(uov(nel))
        allocate(vov(nel))

        allocate(unv(nel))
        allocate(vnv(nel))

        end subroutine mod_hydro_baro_init

!************************************************************

        end module mod_hydro_baro
