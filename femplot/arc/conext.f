
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

	programm conext

c extracts special point from con file

	implicit none

c parameters
	integer nkndim, neldim, ngrdim
	parameter (nkndim=11000,neldim=22000)
	parameter (ngrdim=11)
	integer matdim
	parameter (matdim=nkndim*ngrdim)
c description
	character*80 descrr,descrp
	common /descrr/ descrr
	common /descrp/ descrp
c FEM parameters
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real eps,eps2,pi,flag,high
	real grav,fcor,dcor,dirn,rowass,roluft
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /mkonst/ eps,eps2,pi,flag,high
	common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft
c basin
	real xgv(nkndim), ygv(nkndim)
	real hm3v(3,neldim)
	integer nen3v(3,neldim)
	integer ipev(neldim), ipv(nkndim)
	integer iarv(neldim)
	common /xgv/xgv, /ygv/ygv
	common /hm3v/hm3v
	common /nen3v/nen3v
	common /ipev/ipev, /ipv/ipv
	common /iarv/iarv
c boundary etc.
	integer kantv(2,nkndim)
	real hv(nkndim), hev(neldim)
	real hetv(neldim)
	real parray(nkndim)
	common /kantv/kantv
	common /hv/hv, /hev/hev
	common /hetv/hetv
	common /parray/parray
c auxiliary
	real v1v(nkndim)
	common /v1v/v1v
	real v2v(nkndim)
	common /v2v/v2v
	real vev(neldim)
	common /vev/vev
	real amat(matdim)
	common /amat/amat
	real uvnv(neldim), vvnv(neldim)
	common /uvnv/uvnv, /vvnv/vvnv
	real uv(nkndim), vv(nkndim)
	common /uv/uv, /vv/vv
c local
	integer mode
	integer iapini
	real getpar

	eps=1.e-5
	eps2=1.e-6
	pi=3.141592653
	flag=-9988765.0
	high=1.e+35
	grav=9.81

c======================================================================

c read basin

	if(iapini(7,nkndim,neldim,matdim).eq.0) then
		stop 'error stop : iapini'
	end if

c make depth on nodes and elements

	call mkhv(hv,v1v,nkn,nel)
	call mkhev(hev,nel)

	call ichoice(mode)

	if( mode .eq. 1 ) call plobas
	if( mode .eq. 2 ) call plosim(.true.)
	if( mode .eq. 3 ) call plosim(.false.)
	if( mode .eq. 4 ) call plozet
	if( mode .eq. 5 ) call plocon('.con')
	if( mode .eq. 6 ) call plocon('.tem')
	if( mode .eq. 7 ) call plocon('.sal')
	if( mode .eq. 8 ) call plocon('.rms')

	write(6,*) 'data written to file 99'

	end

c**********************************************************

	subroutine plocon(type)

	implicit none

	character*(*) type

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real parray(1)
	common /parray/parray

	integer it
	logical connext,oktime

	call conopen(type)

	do while( connext(it,nkn,1,parray) )
	  write(99,*) it,parray(1594)
	end do

	call conclose

	end

c**********************************************************

	subroutine plobas
	end

	subroutine plozet
	end

	subroutine plosim(bdum)
	logical bdum
	end

c**********************************************************

	subroutine ichoice(mode)

	implicit none

	integer mode

	write(6,*)
	write(6,*) ' basin .............. 1'
	write(6,*) ' velocity ........... 2'
	write(6,*) ' transport .......... 3'
	write(6,*) ' water level ........ 4'
	write(6,*) ' concentration ...... 5'
	write(6,*) ' temperature ........ 6'
	write(6,*) ' salinity ........... 7'
	write(6,*) ' rms ................ 8'
	write(6,*)
	write(6,*) ' Enter choice : '
	write(6,*)

	read(5,'(i10)') mode

	end

c**********************************************************

