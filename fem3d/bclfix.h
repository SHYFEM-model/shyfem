
	integer nbadim			!dimension of boundary array
	parameter (nbadim = 20000)

	integer nlxdim			!this must be bigger than bound. data
	parameter (nlxdim = nlvdim)

	character*80 fixfile
	logical bfix,bosigma
	integer nbfix,lbmax
	real tramp,tnudge,vfact

	common /fixfile/fixfile
	common /fixlog/bfix,bosigma
	common /fixint/nbfix,lbmax
	common /fixreal/tramp,tnudge,vfact

	save /fixfile/, /fixlog/, /fixint/, /fixreal/

        integer ielfix(0:3,neldim)
        common /ielfix/ielfix
        save /ielfix/

