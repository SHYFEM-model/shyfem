
c*****************************************************************

        subroutine mod_bnd_init(nb)
	implicit none
	integer nb
	end

        subroutine mod_bnd_adjust(nb)
	implicit none
	integer nb
	end

        subroutine mod_bnd_reinit(nb)
	implicit none
	integer nb
	end

        subroutine mod_irv_init(nrb)
        implicit none
        integer nrb
        end

c*****************************************************************

        subroutine ev_init(nel)
	implicit none
	integer nel
	end

        subroutine mod_geom_init(nkn,nel,ngr)
	implicit none
	integer nkn,nel,ngr
	end

        subroutine mod_depth_init(nkn,nel,nlv)
	implicit none
	integer nkn,nel,nlv
	end

        subroutine mod_plot2d_init(nkn,nel)
	implicit none
	integer nkn,nel
	end

        subroutine mod_hydro_plot_init(nkn,nel)
	implicit none
	integer nkn,nel
	end

c*****************************************************************

        subroutine levels_init(nkn,nel,nlv)
	implicit none
	integer nkn,nel,nlv
	end

        subroutine mod_hydro_init(nkn,nel,nlv)
	implicit none
	integer nkn,nel,nlv
	end

        subroutine mod_hydro_vel_init(nkn,nel,nlv)
	implicit none
	integer nkn,nel,nlv
	end

        subroutine mod_hydro_print_init(nkn,nlv)
	implicit none
	integer nkn,nlv
	end

        subroutine mod_plot3d_init(nkn,nel,nlv)
	implicit none
	integer nkn,nel,nlv
	end

c*****************************************************************

c*****************************************************************

        subroutine levels_get_dimension(nl)

! returns vertical dimension (static)

        integer nl

	include 'param.h'

        nl = nlvdim

        end subroutine levels_get_dimension

c*****************************************************************

