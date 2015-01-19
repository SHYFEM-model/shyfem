
        real rdistv(nkndim)
        common /rdistv/rdistv

        real fcorv(neldim)
        common /fcorv/fcorv
        real fxv(nlvdim,neldim)          !new HYDRO debora
        common /fxv/fxv
        real fyv(nlvdim,neldim)
        common /fyv/fyv

	save /rdistv/,/fcorv/,/fxv/,/fyv/

        integer iuvfix(neldim)
        common /iuvfix/iuvfix
	save /iuvfix/

        double precision ddxv(2*nlvdim,neldim)  !ASYM
        double precision ddyv(2*nlvdim,neldim)  !ASYM
        common /ddxv/ddxv, /ddyv/ddyv
        save /ddxv/, /ddyv/

