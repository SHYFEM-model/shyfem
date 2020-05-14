
!--------------------------------------------------------------------------
!
!    Copyright (C) 2018-2019  Georg Umgiesser
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

! utility routines to convert between 3D arrays and linear array
!
! contents :
!
! revision log :
!
! 09.11.2018	ggu	written from scartch
! 11.11.2018	ggu	some bug fixes, extension to double
! 12.11.2018	ggu	utility routines for easy reading/writing
! 18.12.2018	ggu	changed VERS_7_5_52
! 14.02.2019	ggu	changed VERS_7_5_56
! 16.02.2019	ggu	changed VERS_7_5_60
!
! usage :
!
!	call count_linear
!	allocate(rlin)
!	read rlin
!	call linear2vals
!
!	call count_linear
!	allocate(rlin)
!	call vals2linear
!	write rlin
!
!************************************************************************

        subroutine count_linear(nlvddi,n,m,il,nlin)

        implicit none

        integer nlvddi,n,m
        integer il(n)
	integer nlin

	integer j,lmax

        nlin = 0

        do j=1,n
          lmax = min(nlvddi,il(j))
	  nlin = nlin + lmax
	end do

	nlin = nlin * m

	end

!************************************************************************

        subroutine vals2linear(nlvddi,n,m,il,vals,rlin,nlin)

        implicit none

        integer nlvddi,n,m
        integer il(n)
        real vals(nlvddi,n,m)
        real rlin(nlin)
        integer nlin

        integer i,j,lmax,nl,ne

        nl = 0

        do i=1,m
          do j=1,n
            lmax = min(nlvddi,il(j))
	    ne = nl + lmax
	    if( ne > nlin ) goto 99
            rlin(nl+1:ne) = vals(1:lmax,j,i)
	    nl = ne
          end do
        end do

	nlin = nl

	return
   99	continue
	write(6,*) nl,ne,nlin
	stop 'error stop vals2linear: nl>nlin'
        end

!************************************************************************

        subroutine linear2vals(nlvddi,n,m,il,vals,rlin,nlin)

        implicit none

        integer nlvddi,n,m
        integer il(n)
        real vals(nlvddi,n,m)
        real rlin(nlvddi*n*m)
        integer nlin

        integer i,j,lmax,nl,ne

	ne = 0
        nl = 0

        do i=1,m
          do j=1,n
            lmax = min(nlvddi,il(j))
	    ne = nl + lmax
	    if( ne > nlin ) goto 99
            vals(1:lmax,j,i) = rlin(nl+1:ne)
	    nl = ne
          end do
        end do

	if( nl /= nlin ) stop 'error stop linear2vals: nl/=nlin'

	return
   99	continue
	write(6,*) nl,ne,nlin
	stop 'error stop linear2vals: nl>nlin'
        end

!************************************************************************

        subroutine dvals2linear(nlvddi,n,m,il,vals,rlin,nlin)

        implicit none

        integer nlvddi,n,m
        integer il(n)
        double precision vals(nlvddi,n,m)
        double precision rlin(nlin)
        integer nlin

        integer i,j,lmax,nl,ne

        nl = 0

        do i=1,m
          do j=1,n
            lmax = min(nlvddi,il(j))
	    ne = nl + lmax
	    if( ne > nlin ) goto 99
            rlin(nl+1:ne) = vals(1:lmax,j,i)
	    nl = ne
          end do
        end do

	nlin = nl

	return
   99	continue
	write(6,*) nl,ne,nlin
	stop 'error stop dvals2linear: nl>nlin'
        end

!************************************************************************

        subroutine dlinear2vals(nlvddi,n,m,il,vals,rlin,nlin)

        implicit none

        integer nlvddi,n,m
        integer il(n)
        double precision vals(nlvddi,n,m)
        double precision rlin(nlvddi*n*m)
        integer nlin

        integer i,j,lmax,nl,ne

	ne = 0
        nl = 0

        do i=1,m
          do j=1,n
            lmax = min(nlvddi,il(j))
	    ne = nl + lmax
	    if( ne > nlin ) goto 99
            vals(1:lmax,j,i) = rlin(nl+1:ne)
	    nl = ne
          end do
        end do

	if( nl /= nlin ) stop 'error stop dlinear2vals: nl/=nlin'

	return
   99	continue
	write(6,*) nl,ne,nlin
	stop 'error stop dlinear2vals: nl>nlin'
        end

!************************************************************************

        subroutine linear2read(iunit,nlvddi,n,il,vals,ierr)

	implicit none

	integer iunit
	integer nlvddi,n
	integer il(n)
	real vals(nlvddi,n)
	integer ierr

	integer m,nlin,i
	real, allocatable :: rlin(:)

	m = 1

        call count_linear(nlvddi,n,m,il,nlin)
        allocate(rlin(nlin))
        read(iunit,iostat=ierr) (rlin(i),i=1,nlin)
        call linear2vals(nlvddi,n,m,il,vals,rlin,nlin)

	end

!************************************************************************

        subroutine linear2write(iunit,nlvddi,n,il,vals,ierr)

	implicit none

	integer iunit
	integer nlvddi,n
	integer il(n)
	real vals(nlvddi,n)
	integer ierr

	integer m,nlin,i
	real, allocatable :: rlin(:)

	m = 1

        call count_linear(nlvddi,n,m,il,nlin)
        allocate(rlin(nlin))
        call vals2linear(nlvddi,n,m,il,vals,rlin,nlin)
        write(iunit,iostat=ierr) (rlin(i),i=1,nlin)

	end

!************************************************************************

