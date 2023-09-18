
        module fluidmud

        implicit none

        !double precision rhosed		!Mud primary particle density (kg/m3)
        !common /rhosed/ rhosed
	!save /rhosed/

 	!double precision dm0,nf
 	!common /dm0/ dm0
	!common /nf/ nf
	!save /dm0/,/nf/

        !double precision z0bkmud(nkndim)       !bottom roughenss on nodes for mud
        !common /z0bkmud/z0bkmud
        !save /z0bkmud/

        !double precision mudc(nlvdim,nkndim)        !Fluid mud concentrationarray (kg/m3)
        !common /mudc/mudc
        !double precision rhomud(nlvdim,nkndim) !Mud floc part. density (kg/m3)
        !common /rhomud/rhomud
        !save /mudc/,/rhomud/

        !double precision visv_yield(0:nlvdim,nkndim) !viscosity (mud)
        !common /visv_yield/visv_yield
        !double precision difv_yield(0:nlvdim,nkndim) !diffusivity (mud)
        !common /difv_yield/difv_yield
	!save /visv_yield/,/difv_yield/

        !double precision lambda(nlvdim,nkndim)    ! Structural parameter
        !common /lambda/lambda
        !double precision vts(0:nlvdim,nkndim)        ! Rheological Viscosity [m2/s]
        !common /vts/vts
        !double precision dmf_mud(nlvdim,nkndim)  ! Floc size array.
        !common /dmf_mud/dmf_mud
	!save /lambda/,/vts/,/dmf_mud/

        !double precision wprvs(0:nlvdim,nkndim)    !water density
        !common /wprvs/wprvs
	!save /wprvs/


        integer, private, save  :: nkn_fmud = 0
        integer, private, save  :: nlv_fmud = 0

        double precision, private, save     :: rhosed  ! Mud primary particle density (kg/m3)
        double precision, private, save     :: dm0     ! ccf viene definita anche in submud.f ???
        double precision, private, save     :: nf      ! ccf scalare o array???

        double precision, allocatable, save :: z0bkmud(:)          ! bottom roughenss on nodes for mud
        double precision, allocatable, save :: mudc(:,:)           ! Fluid mud concentrationarray (kg/m3)
        double precision, allocatable, save :: rhomud(:,:) ! Mud floc part. density (kg/m3)
        double precision, allocatable, save :: visv_yield(:,:)     ! viscosity (mud)                     
        double precision, allocatable, save :: diff_yield(:,:)     ! diffusivity (mud)                     
        double precision, allocatable, save :: lambda(:,:)         ! Structural parameter                  
        double precision, allocatable, save :: vts(:,:)            ! Rheological Viscosity [m2/s]          
        double precision, allocatable, save :: dmf_mud(:,:)        ! Floc size array.                     
        double precision, allocatable, save :: wprvs(:,:)          ! Water density (kg/m3)              

        contains

!************************************************************

        subroutine mod_fluidmud_init(nkn,nlv)

        integer  :: nkn
        integer  :: nlv

        if( nkn == nkn_fmud .and. nlv == nlv_fmud ) return

        if( nkn > 0 .or. nlv > 0 ) then
          if( nkn == 0 .or. nlv == 0 ) then
            write(6,*) 'nkn,nlv: ',nkn,nlv
            stop 'error stop mod_fluidmud_init: incompatible parms'
          end if
        end if

        if( nkn_fmud > 0 ) then
          deallocate(z0bkmud)
          deallocate(mudc)
          deallocate(rhomud)
          deallocate(visv_yield)
          deallocate(diff_yield)
          deallocate(lambda)
          deallocate(vts)
          deallocate(dmf_mud)
          deallocate(wprvs)
        end if

        nkn_fmud = nkn
        nlv_fmud = nlv

        if( nkn == 0 ) return

        allocate(z0bkmud(nkn))
        allocate(mudc(nlv,nkn))
        allocate(rhomud(nlv,nkn))
        allocate(visv_yield(0:nlv,nkn))
        allocate(diff_yield(0:nlv,nkn))
        allocate(lambda(nlv,nkn))
        allocate(vts(0:nlv,nkn))
        allocate(dmf_mud(nlv,nkn))
        allocate(wprvs(0:nlv,nkn))

        end subroutine mod_fluidmud_init

!************************************************************

        subroutine mod_fluidmud_dummy_init(nkn,nlv)

	integer nkn,nlv

        allocate(vts(0:nlv,nkn))
	vts = 0.

        end subroutine mod_fluidmud_dummy_init

!************************************************************

        end module fluidmud

