
	integer ip_iunit		!unit number, -1 for not open
	integer ip_nintp		!interpolation, 2: linear, 4: cubic
	integer ip_nvar			!number of variables stored
	integer ip_nsize		!total size for one variable
	integer ip_ndata		!total size for all variables
	integer ip_ndim			!dimension of array
	integer ip_nextra		!extra header information
	integer ip_ires			!pointer to where results are stored
	integer ip_nspace		!space needed for all info in array
	integer ip_np			!number of horizontal points
	integer ip_lmax			!max number of levels
	integer ip_iformat		!is formatted?

c np = 0, nsize = 0	-> time series
c np > 0, nsize > 0	-> fem file format

	parameter( ip_iunit   =  1 )
	parameter( ip_nintp   =  2 )
	parameter( ip_nvar    =  3 )
	parameter( ip_nsize   =  4 )
	parameter( ip_ndata   =  5 )
	parameter( ip_ndim    =  6 )
	parameter( ip_nextra  =  7 )
	parameter( ip_ires    =  8 )
	parameter( ip_nspace  =  9 )
	parameter( ip_np      = 10 )
	parameter( ip_lmax    = 11 )
	parameter( ip_iformat = 12 )

        integer nextra
        parameter( nextra = 13 )

        real rguard
        parameter( rguard = 1.234543e+20 )

