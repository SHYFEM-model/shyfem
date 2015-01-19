
	double precision numv_gotm(0:nlvdim,nkndim)
	double precision nuhv_gotm(0:nlvdim,nkndim)
	double precision tken_gotm(0:nlvdim,nkndim)
	double precision eps_gotm(0:nlvdim,nkndim)
	double precision rls_gotm(0:nlvdim,nkndim)

        common /numv_gotm/ numv_gotm
        common /nuhv_gotm/ nuhv_gotm
        common /tken_gotm/ tken_gotm
        common /eps_gotm/ eps_gotm
        common /rls_gotm/ rls_gotm

	save /numv_gotm/,/nuhv_gotm/
	save /tken_gotm/,/eps_gotm/,/rls_gotm/

        real shearf2(nlvdim,nkndim)
        common /shearf2/shearf2
        real buoyf2(nlvdim,nkndim)
        common /buoyf2/buoyf2
	save /shearf2/,/buoyf2/

