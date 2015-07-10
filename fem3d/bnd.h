
	integer nbvdim
	!parameter(nbvdim=100)
	parameter(nbvdim=25)

        real bnd(nbvdim,nbcdim)
        common /bnd/bnd
	save /bnd/

