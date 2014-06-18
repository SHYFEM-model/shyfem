
	integer nbndim
	parameter(nbndim=100)

        real bnd(nbndim,nbcdim)
        common /bnd/bnd
	save /bnd/

