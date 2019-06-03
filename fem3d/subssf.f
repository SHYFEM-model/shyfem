
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

! revision log :
!
! 23.03.2010	ggu	changed v6.1.1
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62

!***********************************************************************

	subroutine read_ts(file,ndim,nvar,n,itime,array)

! reads time series file

	implicit none

	character*(*) file		!file name
	integer ndim			!dimension of arrays
	integer nvar			!number of expected variables in file
	integer n			!number of records read (return)
	integer itime(ndim)		!time column
	real array(ndim,nvar)		!value column(s)

	integer naux
	parameter(naux=100)

	integer iunit,i,it
	real aux(naux)
	integer ifileo

	if( nvar .gt. naux ) stop 'error stop read_ts: naux'

	iunit = 55
        iunit = ifileo(iunit,file,'form','old')
        if( iunit .le. 0 ) stop 'error stop read_ts: no such file'

	n = 0
    1	continue
	  read(iunit,*,end=2) it,(aux(i),i=1,nvar)
	  n = n + 1
	  if( n .gt. ndim ) stop 'error stop read_ts: ndim'
	  itime(n) = it
	  do i=1,nvar
	    array(n,i) = aux(i)
	  end do
	goto 1
    2	continue

	close(iunit)

	end

!***********************************************************************

