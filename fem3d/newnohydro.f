c
c routines for non hydrostatic terms
c
c revision log :
c
c 10.05.2013    dbf     written from scratch
c 31.05.2013    dbf     written from scratch
c 17.06.2016    ggu&wjm adapted to new version
c
c********************************************************************

	subroutine nonhydro_init

	implicit none

	integer inohyd
	real getpar

	inohyd = nint(getpar('inohyd'))
	if( inohyd /= 0 ) then
	  write(6,*) 'inohyd = ',inohyd
	  stop 'error stop nonhydro_init: cannot run non-hydrostatic'
	end if

	end

c********************************************************************

	subroutine nonhydro_get_flag(bnohyd)

	implicit none

	logical bnohyd

	bnohyd = .false.

	end

c********************************************************************

	subroutine nonhydro_adjust

	end

c********************************************************************

	subroutine nonhydro_copy

        end

c********************************************************************

	subroutine sp256wnh 

	end

c********************************************************************

	subroutine nonhydro_set_explicit 

	end 

c********************************************************************

	subroutine nonhydro_prepare_matrix

	end

c********************************************************************

        subroutine nonhydro_adjust_value

	end

c********************************************************************

	subroutine nonhydro_correct_uveta

	end

c********************************************************************

	subroutine nh_handle_output(dtime)

	implicit none

	double precision dtime

	end

c********************************************************************

