
	program nc_include

	implicit none

	include 'netcdf.inc'

	integer ncid,retval
	character*80 file

	file='test.nc'
	retval = nf_open(file, nf_nowrite, ncid)

	end
