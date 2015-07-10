
        real arfvlv(nkndim)
        common /arfvlv/ arfvlv

        logical bwater(neldim)
        common /bwater/ bwater
        logical bkwater(nkndim)
        common /bkwater/ bkwater

        real hetv(neldim)
        common /hetv/ hetv

        real parray(neldim)             !is good for nodes and elements
        common /parray/ parray

	save /arfvlv/,/bwater/,/bkwater/
	save /hetv/,/parray/

