
        real rfricv(neldim)
        common /rfricv/rfricv
        real czv(neldim)
        common /czv/czv

	save /rfricv/,/czv/

        real austv(neldim)                     !$$AUST
        common /austv/austv
        real difhv(nlvdim,neldim)      !horizontal diffusion - 3D
        common /difhv/difhv

	save /austv/,/difhv/

        !real wdifhv(3,3,neldim)       !weights for horizontal diff.
        !common /wdifhv/wdifhv
	!save /wdifhv/

        real visv(0:nlvdim,nkndim)      !viscosity (momentum)
        common /visv/visv
        real difv(0:nlvdim,nkndim)      !diffusivity (scalars)
        common /difv/difv

	save /visv/,/difv/

