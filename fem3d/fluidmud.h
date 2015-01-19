
        double precision rhosed		!Mud primary particle density (kg/m3)
        common /rhosed/ rhosed
	save /rhosed/

 	double precision dm0,nf
 	common /dm0/ dm0
	common /nf/ nf
	save /dm0/,/nf/

        real z0bkmud(nkndim)       !bottom roughenss on nodes for mud
        common /z0bkmud/z0bkmud
        save /z0bkmud/

        real mudc(nlvdim,nkndim)        !Fluid mud concentrationarray (kg/m3)
        common /mudc/mudc
        double precision rhomud(nlvdim,nkndim) !Mud floc part. density (kg/m3)
        common /rhomud/rhomud
        save /mudc/,/rhomud/

        real visv_yield(0:nlvdim,nkndim) !viscosity (mud)
        common /visv_yield/visv_yield
        real difv_yield(0:nlvdim,nkndim) !diffusivity (mud)
        common /difv_yield/difv_yield
	save /visv_yield/,/difv_yield/

        real lambda(nlvdim,nkndim)    ! Structural parameter
        common /lambda/lambda
        real vts(0:nlvdim,nkndim)        ! Rheological Viscosity [m2/s]
        common /vts/vts
        real dmf_mud(nlvdim,nkndim)  ! Floc size array.
        common /dmf_mud/dmf_mud
	save /lambda/,/vts/,/dmf_mud/

        real wprvs(0:nlvdim,nkndim)    !water density
        common /wprvs/wprvs
	save /wprvs/

