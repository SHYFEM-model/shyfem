	module mod_diff_visc_fric

	implicit none

        !real rfricv(neldim)
        !common /rfricv/rfricv
        !real czv(neldim)
        !common /czv/czv
        !save /rfricv/,/czv/
        !real austv(neldim)                     !$$AUST
        !common /austv/austv
        !real difhv(nlvdim,neldim)      !horizontal diffusion - 3D
        !common /difhv/difhv
        !save /austv/,/difhv/
	!--- Commented in the old code
        !!!real wdifhv(3,3,neldim)       !weights for horizontal diff.
        !!!common /wdifhv/wdifhv
        !!!save /wdifhv/
	!---
        !real visv(0:nlvdim,nkndim)      !viscosity (momentum)
        !common /visv/visv
        !real difv(0:nlvdim,nkndim)      !diffusivity (scalars)
        !common /difv/difv
        !save /visv/,/difv/

	integer, private, save :: nkn_diff_visc_fric = 0
	integer, private, save :: nel_diff_visc_fric = 0
	integer, private, save :: nlv_diff_visc_fric = 0

	real, allocatable, save :: rfricv(:)
	real, allocatable, save :: czv(:)
	real, allocatable, save :: austv(:)
	real, allocatable, save :: difhv(:,:)
	!!!real, allocatable, save :: wdifhv(3,3,:)
	real, allocatable, save :: visv(:,:)
	real, allocatable, save :: difv(:,:)

	contains

!************************************************************

        subroutine mod_diff_visc_fric_init(nkn,nel,nlv)

        integer nkn, nel, nlv

        if( nkn == nkn_diff_visc_fric .and. nel == nel_diff_visc_fric
     +      .and. nlv == nlv_diff_visc_fric ) return

        if( nkn_diff_visc_fric > 0 .or. nel_diff_visc_fric > 0 .or.
     +      nlv_diff_visc_fric > 0 ) then
          deallocate(rfricv)
          deallocate(czv)
          deallocate(austv)
          deallocate(difhv)
          !!!deallocate(wdifhv)
          deallocate(visv)
          deallocate(difv)
        end if

        if( nkn == 0 ) return
        if( nel == 0 ) return
        if( nlv == 0 ) return

        nkn_diff_visc_fric = nkn
        nel_diff_visc_fric = nel
        nlv_diff_visc_fric = nlv

        allocate(rfricv(nel))
        allocate(czv(nel))
        allocate(austv(nel))
        allocate(difhv(nlv,nel))
        !!!allocate(wdifhv(3,3,nel))
        allocate(visv(0:nlv,nkn))
        allocate(difv(0:nlv,nkn))

        end subroutine mod_diff_visc_fric_init

!************************************************************

        end module mod_diff_visc_fric
