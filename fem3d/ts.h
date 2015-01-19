
	real rhov(nlvdim,nkndim)
	common /rhov/rhov

        real saltv(nlvdim,nkndim)
        common /saltv/saltv
        real tempv(nlvdim,nkndim)
        common /tempv/tempv

        real sobsv(nlvdim,nkndim)
        common /sobsv/sobsv
        real tobsv(nlvdim,nkndim)
        common /tobsv/tobsv
        real rtauv(nlvdim,nkndim)      !relaxation time
        common /rtauv/rtauv

        real bpresv(nlvdim,nkndim)
        common /bpresv/bpresv
        real bpresxv(nlvdim,neldim)
        common /bpresxv/bpresxv
        real bpresyv(nlvdim,neldim)
        common /bpresyv/bpresyv

	save /rhov/,/saltv/,/tempv/
	save /sobsv/,/tobsv/,/rtauv/
	save /bpresv/,/bpresxv/,/bpresyv/

