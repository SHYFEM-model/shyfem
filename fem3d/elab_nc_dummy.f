
!********************************************************************

	subroutine nc_output_init(ncid,title,nvar,ivars)

	implicit none

	integer ncid			!id of file (return)
	character*(*) title		!name of simulation
	integer nvar			!total number of variables to be written
	integer ivars(nvar)		!variable id of SHYFEM

	write(6,*) 'Cannot initialize NETCDF module'
	write(6,*) 'NETCDF support has not been enabled'
	write(6,*) 'Please enable NETCDF support in Rules.make'
	write(6,*) '(set "NETCDF=true")'

	stop 'error stop nc_output_init: no netcdf support'

	end

!********************************************************************

	subroutine nc_output_record(ncid,var_id,cv3)

	use basin
	use levels

	implicit none

	integer ncid
	integer var_id
	real cv3(nlvdi,nkn)

	end

!********************************************************************

	subroutine nc_output_hydro(ncid,znv,uprv,vprv)

	use basin
	use levels

	implicit none

	integer ncid
	real znv(nkn)
	real uprv(nlvdi,nkn)
	real vprv(nlvdi,nkn)

	end

!********************************************************************

	subroutine nc_output_time(ncid,dtime)

	use shyelab_out

	implicit none

	integer ncid
	double precision dtime

	end

!********************************************************************

	subroutine nc_output_final(ncid)

	implicit none

	integer ncid

	end

!********************************************************************

