
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
! 31.03.2017	ggu	changed VERS_7_5_24
! 06.07.2018	ggu	changed VERS_7_5_48
! 14.02.2019	ggu	changed VERS_7_5_56
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62

!***************************************************************
!
! these routines read line(s) defining area where averaging will be done
!
! old file format: xy
! new file format: grd, bnd, xy
!
! ieflag and ikflag are set and used in elab routines
!
!***************************************************************

        subroutine handle_area

! this works for more than one line

        use basin
        use elabutil

        implicit none

	integer is,ie,nn,il
        integer		  :: np
	real, allocatable :: xx(:),yy(:)
	integer, allocatable :: ifl(:)
	integer, allocatable :: ikf(:),ief(:)

        ieflag = 1
        ikflag = 1

	if( .not. barea ) return

	np = 0
        call read_all_lines(areafile,np,xx,yy,ifl)
        if( np <= 0 ) goto 99
        allocate(xx(np),yy(np),ifl(np))
        call read_all_lines(areafile,np,xx,yy,ifl)

        ieflag = -1
        ikflag = 0
        allocate(ikf(nkn),ief(nel))

	ie = 0
	il = 0
	do
	  is = ie + 1
	  call get_next_line(np,ifl,is,ie,nn)
	  if( nn == 0 ) exit
	  il = il + 1
	  write(6,*) 'reading line: ',il,nn,is,ie
          call check_elements(nn,xx(is:ie),yy(is:ie),ief,ikf)
	  where( ikf == 1 ) ikflag = 1
	  where( ief == 1 ) ieflag = 1
	  where( ief == 0 .and. ieflag < 0 ) ieflag = 0
	end do

	call write_grd_with_flags('flags.grd',ikflag,ieflag)

	baverbas = .true.

	return
   99	continue
	write(6,*) 'error reading area file: ',trim(areafile)
	stop 'error stop handle_area: area file'
	end subroutine handle_area

!***************************************************************
!
! Read coordinates of points in line
! for n == 0 only checks how many nodes to read
! for n > 0 reads nodes into nodes() (error if n is too small)
!
!***************************************************************

	subroutine get_next_line(n,ifl,is,ie,nn)

	implicit none

	integer istart,n,is,ie,nn
	integer ifl(n)

	integer i

	nn = 0
	if( is > n ) return
	if( ifl(is) /= 1 ) goto 99
	ie = is

	do i=is+1,n
	  if( ifl(i) == 1 ) exit
	  ie = i
	end do

	nn = ie - is + 1

	return
   99	continue
	write(6,*) 'error in line file: ',is,ifl(is)
	stop 'error stop get_next_line: file structure'
	end

!***************************************************************

	subroutine write_grd_with_flags(file,ikf,ief)

	use basin

	implicit none

	character*(*) file
	integer ikf(nkn)
	integer ief(nel)

	integer k,ie,ii
	real x,y

	open(1,file=file,status='unknown',form='formatted')

	do k=1,nkn
	  x = xgv(k)
	  y = ygv(k)
	  write(1,1000) 1,k,ikf(k),x,y
	end do

	do ie=1,nel
	  x = xgv(k)
	  y = ygv(k)
	  write(1,2000) 2,ie,ief(ie),3,(nen3v(ii,ie),ii=1,3)
	end do

	close(1)

	return
 1000	format(i1,2i10,2f14.6)
 2000	format(i1,6i10)
	end

!***************************************************************

