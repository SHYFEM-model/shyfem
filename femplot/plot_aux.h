
        real fvlv(nlvdim,nkndim)
        common /fvlv/ fvlv
        real arfvlv(nkndim)
        common /arfvlv/ arfvlv

        real wauxv(0:nlvdim,nkndim)
        common /wauxv/ wauxv

        logical bwater(neldim)
        common /bwater/ bwater
        logical bkwater(nkndim)
        common /bkwater/ bkwater

        real hetv(neldim)
        common /hetv/ hetv
        real het3v(nlvdim,neldim)
        common /het3v/ het3v
        !real hl(nlvdim)
        !common /hl/ hl

        real parray(neldim)             !is good for nodes and elements
        common /parray/ parray
        real p3(nlvdim,2*neldim)        !is good for nodes, elements & arrays
        common /p3/ p3

	save /fvlv/,/arfvlv/,/wauxv/,/bwater/,/bkwater/
	save /hetv/,/het3v/,/parray/,/p3/
	!save /hetv/,/het3v/,/hl/,/parray/,/p3/

