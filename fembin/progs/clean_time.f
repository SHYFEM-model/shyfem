
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

	program clean_time

c cleans time series from syncronization time records

	implicit none

	integer ndim	!number of values to read
	parameter(ndim=6)

	integer idtsyn,time
	integer ival(2)
	real rval(3)

	read(5,*) idtsyn

    1	continue
	  read(5,*,end=2) time,rval,ival
	  if( mod(time,idtsyn) .eq. 0 ) goto 1
	  write(6,*) time,rval,ival
	  goto 1
    2	continue

	end

