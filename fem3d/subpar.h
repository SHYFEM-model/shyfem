
c-------------------------------------------

	integer nnamdi		!total number of parameters
	integer nichdi		!total number of strings
	integer niardi		!total number of arrays
	integer narrdi		!total number of values in arrays

	parameter (nnamdi=300)
	parameter (niardi=20)
	parameter (nichdi=100)
	parameter (narrdi=1000)

c-------------------------------------------

	character*6 nampar(nnamdi)      	!names of parameters
	character*6 secpar(nnamdi)      	!section names
	integer itypar(nnamdi)			!type of parameters

c itypar:  1: numeric  2: numeric array  3: string

	double precision valpar(nnamdi) 	!parameter values

	character*80 chapar(nichdi)		!strings

	integer ip_arrpar(2,niardi)		!filling of arrays
	double precision arrpar(narrdi)		!values in arrays

	common /d_par/valpar,arrpar
	common /i_par/itypar,ip_arrpar
	common /c6_par/nampar,secpar
	common /c80_par/chapar

	save /d_par/, /i_par/, /c6_par/, /c80_par/

c-------------------------------------------

	character*6 actsec,defsec,auxsec	!default section
	integer nentry				!total number of values
	integer incha				!total number of strings
	integer inarr				!total number of arrays
	integer infarr				!filling of arrays

	common /parsco/actsec,defsec,auxsec
	common /iac_par/nentry,incha,inarr,infarr

	save /parsco/, /iac_par/

c-------------------------------------------


