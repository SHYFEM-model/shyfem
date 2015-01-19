
        real ulov(nlvdim,neldim)
        common /ulov/ulov
        real ulnv(nlvdim,neldim)
        common /ulnv/ulnv
        real vlov(nlvdim,neldim)
        common /vlov/vlov
        real vlnv(nlvdim,neldim)
        common /vlnv/vlnv
        real wlov(0:nlvdim,nkndim)
        common /wlov/wlov
        real wlnv(0:nlvdim,nkndim)
        common /wlnv/wlnv

	save /ulov/,/ulnv/,/vlov/,/vlnv/,/wlov/,/wlnv/

