
	integer nfxflxdim
	parameter(nfxflxdim=100)	!total number of nodes in sections

        integer nscflxdim
        parameter(nscflxdim=20)         !maximum number of sections

        integer nsect,kfluxm,kflux(nfxflxdim)
        common /kfluxflxc/ nsect,kfluxm,kflux

        integer iflux(3,nfxflxdim)
        common /ifluxflx/iflux

        save /kfluxflxc/,/ifluxflx/

