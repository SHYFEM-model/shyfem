
        real fvlv(nlvdim,nkndim)
        common /fvlv/ fvlv

        real wauxv(0:nlvdim,nkndim)
        common /wauxv/ wauxv

        real het3v(nlvdim,neldim)
        common /het3v/ het3v

        real p3(nlvdim,2*neldim)        !is good for nodes, elements & arrays
        common /p3/ p3

	save /fvlv/,/wauxv/
	save /het3v/,/p3/

