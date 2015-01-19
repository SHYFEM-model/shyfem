
        real rzv(nkndim), rqv(nkndim)
        common /rzv/rzv, /rqv/rqv

        real rqpsv(nkndim), rqdsv(nkndim)
        common /rqpsv/rqpsv, /rqdsv/rqdsv
        !real evapv(nkndim)
        !common /evapv/evapv
        real mfluxv(nlvdim,nkndim)
        common /mfluxv/mfluxv

	save /rzv/,/rqv/
	save /rqpsv/,/rqdsv/,/mfluxv/

