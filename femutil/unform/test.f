
!--------------------------------------------------------------------------
!
!    Copyright (C) 2019  Georg Umgiesser
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
! 13.03.2019	ggu	changed VERS_7_5_61

	program test_unformatted

	implicit none

	integer iunit,nl
	integer fields(3)

	open(1,file='unform.dat',form='unformatted',status='unknown')
	write(1) 1
	close(1)

	open(1,file='stream.dat',form='unformatted',status='unknown'
     +			,access='stream')
	write(1) 1
	close(1)

	iunit = 1
	open(iunit,file='unform.dat',form='unformatted',status='old'
     +			,access='stream')
	call skip_record(iunit,nl)
	close(iunit)

	write(6,*) nl,' bytes contained in record'

	end

!**************************************************************

	subroutine skip_record(iunit,nl)

	implicit none

	integer iunit
	integer nl

	integer nle

	call read_length(iunit,nl)
	call skip_bytes(iunit,nl)
	call read_length(iunit,nle)

	if( nl /= nle ) then
	  write(6,*) 'cannot read record: ',nl,nle
	  stop 'error stop read_record: nl/=nle'
	end if

	end

!**************************************************************

	subroutine read_length(iunit,nl)

	implicit none

	integer iunit
	integer nl

	read(iunit) nl

	end

!**************************************************************

	subroutine skip_bytes(iunit,nb)

	implicit none

	integer iunit
	integer nb

	integer (kind=1), allocatable :: aux(:)

	allocate(aux(nb))

	read(iunit) aux(1:nb)

	end

!**************************************************************

	subroutine read_words(iunit,nb,words)

	implicit none

	integer iunit
	integer nb
	integer words(*)

	integer nw

	nw = nb / 4
	if( 4*nw /= nb ) then
	  write(6,*) 'nb = ',nb,'  nw = ',nw
	  stop 'error stop read_words: number of bytes no multiple of 4'
	end if

	read(iunit) words(1:nw)

	end

!**************************************************************

