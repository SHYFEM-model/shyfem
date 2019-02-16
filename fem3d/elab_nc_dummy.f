
!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2018  Georg Umgiesser
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

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

        subroutine nc_output_record_reg(ncid,var_id,nlvd,np,cv3)


        implicit none

        integer ncid
        integer var_id
        integer nlvd
        integer np
        real cv3(nlvd,np)

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

