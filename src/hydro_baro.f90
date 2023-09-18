
	module hydro_baro

	implicit none

        !double precision uov(neldim), vov(neldim)
        !common /uov/uov, /vov/vov
        !double precision unv(neldim), vnv(neldim)
        !common /unv/unv, /vnv/vnv
        !save /uov/,/vov/,/unv/,/vnv/

	integer, private, save :: nel_hydro_baro = 0

	double precision, allocatable, save :: uov(:), vov(:)
	double precision, allocatable, save :: unv(:), vnv(:)

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

        end module hydro_baro
