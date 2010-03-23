
c-------------------------------------------

	integer nnamdi		!total number of parameters
	integer niardi		!total number of arrays
	integer nichdi		!total number of text strings
	integer narrdi		!total number of values in arrays
	integer nchadi		!total number of characters in text strings

	parameter (nnamdi=200)
	parameter (niardi=20)
	parameter (nichdi=10)
	parameter (narrdi=1000)
	parameter (nchadi=1000)

c-------------------------------------------

	character*6 nampar(nnamdi)      	!names of parameters
	character*6 secpar(nnamdi)      	!section names
	integer itypar(nnamdi)			!type of parameters

	double precision valpar(nnamdi) 	!parameter values

	integer ip_arrpar(2,niardi)		!filling of arrays
	double precision arrpar(narrdi)		!values in arrays

	integer ip_chapar(2,nichdi)		!filling of characters
	character*1 chapar(nchadi)		!characters in arrays

	common /d_par/valpar,arrpar
	common /i_par/itypar,ip_arrpar,ip_chapar
	common /c6_par/nampar,secpar
	common /c_par/chapar

	save /d_par/, /i_par/, /c6_par/, /c_par/

c-------------------------------------------

	character*6 actsec,defsec,auxsec	!default section
	integer nentry				!total number of values
	integer inarr				!total number of arrays
	integer incha				!total number of char strings
	integer infarr				!filling of arrays
	integer infcha				!filling of chars

	common /parsco/actsec,defsec,auxsec
	common /iac_par/nentry,inarr,incha,infarr,infcha

	save /parsco/, /iac_par/

c-------------------------------------------


