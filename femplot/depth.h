
        real hkv(nkndim)
        common /hkv/hkv
        real hev(neldim)
        common /hev/hev

        real hdknv(nlvdim,nkndim)
        common /hdknv/hdknv
        real hdkov(nlvdim,nkndim)
        common /hdkov/hdkov

        real hdenv(nlvdim,neldim)
        common /hdenv/hdenv
        real hdeov(nlvdim,neldim)
        common /hdeov/hdeov

	save /hkv/,/hev/,/hdknv/,/hdkov/,/hdenv/,/hdeov/

	real hkv_min(nkndim), hkv_max(nkndim)
	common /hkv_min/hkv_min, /hkv_max/hkv_max
	save /hkv_min/,/hkv_max/

