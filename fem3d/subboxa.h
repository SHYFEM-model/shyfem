
	integer nfxboxdim
	parameter(nfxboxdim=2000)	!total number of nodes in sections

        integer nscboxdim
        parameter(nscboxdim=140)        !maximum number of sections

        integer nbxdim
        parameter(nbxdim=50)            !maximum number of boxes

        integer nbox,nsect,kfluxm,kflux(nfxboxdim)
        common /kfluxboxc/ nbox,nsect,kfluxm,kflux

        integer iflux(3,nfxboxdim)
        common /ifluxbox/iflux

	integer iboxes(neldim)
	common /iboxes/iboxes

	integer isects(4,nscboxdim)
	common /isects/isects

        save /kfluxboxc/,/ifluxbox/,/iboxes/,/isects/

