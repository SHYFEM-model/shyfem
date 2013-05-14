
	integer nfxboxdim
	parameter(nfxboxdim=2000)	!total number of nodes in sections

        integer nscboxdim
        parameter(nscboxdim=140)        !maximum number of sections

        integer nbxdim
        parameter(nbxdim=50)            !maximum number of boxes

        integer nbox,nsect,nbc_ob,kfluxm,kflux(nfxboxdim)
        common /kfluxboxc/ nbox,nsect,nbc_ob,kfluxm,kflux
	save /kfluxboxc/

        integer iflux(3,nfxboxdim)
        common /ifluxbox/iflux
	save /ifluxbox/

	integer iboxes(neldim)
	common /iboxes/iboxes
	save /iboxes/

	integer ikboxes(nkndim)
	common /ikboxes/ikboxes
	save /ikboxes/

	integer isects(4,nscboxdim)	!internal section description
	common /isects/isects
	save /isects/

	integer iscbnd(4,nbcdim)	!boundary condition descriptions
	common /iscbnd/iscbnd
	save /iscbnd/

