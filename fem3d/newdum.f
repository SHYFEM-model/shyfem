c
c $Id: newdum.f,v 1.5 2010-02-22 15:38:36 georg Exp $
c
c dummy routines for compatibility with 2D version
c
c contents :
c
c subroutine inclos             initializes closing sections
c subroutine rdclos(isc)        reads closing sections
c subroutine ckclos             post-processes closing sections
c subroutine prclos             prints info on closing sections
c subroutine tsclos             tests closing sections
c
c subroutine inoxy              initializes oxygen module
c subroutine rdoxy              reads oxygen module
c subroutine ckoxy              post-processes oxygen module
c subroutine proxy              prints info on oxygen module
c subroutine tsoxy              tests oxygen module
c
c revision log :
c
c 22.05.1998	ggu	created for closing sections
c 22.01.1999	ggu	oxygen module added
c 19.02.2010	ggu	massconc eliminated
c
c*****************************************************************

c dummies for closing section

c        subroutine inclos
c	end
c
c        subroutine rdclos(isc)
c	integer isc
c	write(6,*) 'Closing sections not yet supported in 3d version'
c	write(6,*) 'closing section : ',isc
c	stop 'error stop rdclos'
c	end
c
c        subroutine ckclos
c	end
c
c        subroutine prclos
c	end
c
c        subroutine tsclos
c	end

c*****************************************************************

c dummies for oxygen module

        subroutine inoxy
	implicit none
	end

        subroutine rdoxy
	implicit none
	write(6,*) 'Oxygen module not yet supported in 3d version'
	stop 'error stop rdoxy'
	end

        subroutine ckoxy
	implicit none
	end

        subroutine proxy
	implicit none
	end

        subroutine tsoxy
	implicit none
	end

c*****************************************************************

c*****************************************************************

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

	subroutine mod_hydro_baro_init(nel)
	implicit none
	integer nel
	end

	subroutine mod_diff_visc_fric_init(nkn,nel,nlv)
	implicit none
	integer nkn,nel,nlv
	end

	subroutine mod_roughness_init(nkn)
	implicit none
	integer nkn
	end

	subroutine mod_ts_init(nkn,nlv)
	implicit none
	integer nkn,nlv
	end

	subroutine mod_area_init(nkn,nlv)
	implicit none
	integer nkn,nlv
	end

	subroutine mod_aux_array_init(nkn,nel,nlv)
	implicit none
	integer nkn,nel,nlv
	end

	subroutine mod_bound_dynamic_init(nkn,nlv)
	implicit none
	integer nkn,nlv
	end

	subroutine mod_gotm_aux_init(nkn,nlv)
	implicit none
	integer nkn,nlv
	end

	subroutine mod_diff_aux_init(nel)
	implicit none
	integer nel
	end

        subroutine mod_bnd_aux_init(nkn,nel)
	implicit none
	integer nkn,nel
	end

        subroutine mod_geom_dynamic_init(nkn,nel)
	implicit none
	integer nkn,nel
	end

        subroutine mod_nudging_init(nkn)
	implicit none
	integer nkn
	end

        subroutine mod_depth_init(nkn,nel)
	implicit none
	integer nkn,nel
	end

        subroutine mod_layer_thickness_init(nkn,nel,nlv)
	implicit none
	integer nkn,nel,nlv
	end

        subroutine mod_internal_init(nkn,nel,nlv)
	implicit none
	integer nkn,nel,nlv
	end

        subroutine mod_nohyd_init(nkn,nlv)
	implicit none
	integer nkn,nlv
	end

        subroutine mod_bclfix_init(nkn,nel,nlv)
	implicit none
	integer nkn,nel,nlv
	end

        subroutine mod_fluidmud_init(nkn,nlv)
	implicit none
	integer nkn,nlv
	end

        subroutine mod_sinking_init(nkn,nlv)
	implicit none
	integer nkn,nlv
	end

        subroutine mod_turbulence_init(nkn,nlv)
	implicit none
	integer nkn,nlv
	end

        subroutine mod_waves_init(nkn,nel,nlv)
	implicit none
	integer nkn,nel,nlv
	end

        subroutine mod_meteo_init(nkn)
	implicit none
	integer nkn
	end

	subroutine mod_geom_init(nkn,nel,ngr)
	implicit none
	integer nkn,nel,ngr
	end

	subroutine ev_init(nel)
	implicit none
	integer nel
	end

	subroutine mod_conz_init(ncs,nkn,nlv)
	implicit none
	integer ncs,nkn,nlv
	end

	subroutine mod_bound_geom_init(nkn,nrb)
	implicit none
	integer nkn,nrb
	end

	subroutine mod_bound_geom_reinit(nkn,nrb)
	implicit none
	integer nkn,nrb
	end

	subroutine mod_irv_init(nrb)
	implicit none
	integer nrb
	end

c*****************************************************************

        subroutine mod_bound_geom_info

	use mod_bnd
	use mod_bound_geom
	use basin, only : nkn,nel,ngr,mbw

	include 'param.h'

        integer iu,i

        iu = 88

        write(iu,*) 'mod_bound_geom_info: ',nkn,nrb
        write(iu,*) 'irv: ',nrb,(irv(i),i=1,nrb)
        write(iu,*) 'ierv: ',(ierv(1,i),i=1,nrb)
        write(iu,*) 'ierv: ',(ierv(2,i),i=1,nrb)
        write(iu,*) 'rhv: ',(rhv(i),i=1,nrb)
        write(iu,*) 'rlv: ',(rlv(i),i=1,nrb)
        write(iu,*) 'rrv: ',(rrv(i),i=1,nrb)
        write(iu,*) 'iopbnd: '
        do i=1,nkn
          if( iopbnd(i) .ne. 0 ) write(iu,*) i,iopbnd(i)
        end do
        write(iu,*) 'mod_bound_geom_info end'

        end subroutine mod_bound_geom_info

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

        subroutine mod_tvd_init(nel)
	implicit none
	integer nel
	end

        subroutine mod_tides_init(nkn)
	implicit none
	integer nkn
	end

        subroutine mod_bndo_init(ngr,nrb)
	implicit none
	integer ngr,nrb
	end

        subroutine mod_nudge_init(nkn)
	implicit none
	integer nkn
	end

c*****************************************************************

	subroutine mod_system_init(nkn,nel,mbw)
	implicit none
	integer nkn,nel,mbw
	end

	subroutine mod_system_amat_init(nkn,mbw)
	implicit none
	integer nkn,mbw
	end

c*****************************************************************


