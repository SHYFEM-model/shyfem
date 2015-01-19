
        real uprv(nlvdim,nkndim)
        common /uprv/uprv
        real vprv(nlvdim,nkndim)
        common /vprv/vprv
        real upro(nlvdim,nkndim)
        common /upro/upro
        real vpro(nlvdim,nkndim)
        common /vpro/vpro
        real wprv(0:nlvdim,nkndim)
        common /wprv/wprv

	save /uprv/,/vprv/,/upro/,/vpro/,/wprv/

        real up0v(nkndim)
        common /up0v/up0v
        real vp0v(nkndim)
        common /vp0v/vp0v

	save /up0v/,/vp0v/

        real xv(3,nkndim)
        common /xv/xv

	save /xv/

