
        real v1v(nkndim)
        common /v1v/v1v
        real v2v(nkndim)
        common /v2v/v2v
        real v3v(nkndim)
        common /v3v/v3v
        real ve1v(neldim)
        common /ve1v/ve1v
        real saux1(nlvdim,nkndim)
        common /saux1/saux1
        real saux2(nlvdim,nkndim)
        common /saux2/saux2
        real saux3(nlvdim,nkndim)
        common /saux3/saux3
        real saux4(nlvdim,nkndim)
        common /saux4/saux4
        real sauxe1(nlvdim,neldim)
        common /sauxe1/sauxe1
        real sauxe2(nlvdim,neldim)
        common /sauxe2/sauxe2

	save /v1v/,/v2v/,/v3v/,/ve1v/
	save /saux1/,/saux2/,/saux3/,/saux4/
	save /sauxe1/,/sauxe2/

