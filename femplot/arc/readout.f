
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

	program readout

c reads file OUT

        implicit none

	integer nkndim,neldim
	parameter(nkndim=10000,neldim=2*nkndim)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

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

c OUT
        real xv(3,nkndim)
        real zenv(3,neldim)
        real usnv(neldim), vsnv(neldim)
        common /xv/xv
        common /zenv/zenv
        common /usnv/usnv, /vsnv/vsnv

c aux
        real hv(nkndim), hev(neldim)
        common /hv/hv, /hev/hev
        real v1v(nkndim)
        common /v1v/v1v

        integer nrec,it,ie
        logical outnext
	integer iapini
	real href
	real getpar

        nrec = 0

        if( iapini(3,nkndim,neldim,0) .eq. 0 ) then
                stop 'error stop : iapini'
        end if

        call mkhv(hv,v1v,nkn,nel)
        call mkhev(hev,nel)

        call outopen

	href = getpar('href')
	write(6,*) 'adjusting depth with href = ',href
	do ie=1,nel
	  hev(ie) = hev(ie) - href
	end do

        do while( outnext(it) )
          nrec = nrec + 1
          write(6,*) nrec,it
	  call check
        end do

        call outclose

        end

c**********************************************************

	subroutine check

	implicit none

	real drittl
	parameter(drittl=1./3.)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        real xv(3,1)
        real zenv(3,1)
        real usnv(1), vsnv(1)
        common /xv/xv
        common /zenv/zenv
        common /usnv/usnv, /vsnv/vsnv

	real hev(1)
	common /hev/hev

	integer ie,iemax,iemin,ii
	real u,v,uv,uvmax
	real htmax,htmin,htot,rh,zm

	uvmax = 0.
	iemax = 0
	iemin = 0
	htmax = 0.
	htmin = 999.

	do ie=1,nel
	  u = usnv(ie)
	  v = vsnv(ie)

	  zm = 0.
	  do ii=1,3
	    zm = zm + zenv(ii,ie)
	  end do
	  zm = zm * drittl
	  htot = hev(ie) + zm
	  rh = 1. / htot

	  u = u * rh
	  v = v * rh

	  uv = u*u + v*v
	  if( uv .gt. uvmax ) then
	    uvmax = uv
	    htmax = htot
	    iemax = ie
	  end if

	  if( htot .lt. htmin ) then
	    htmin = htot
	    iemin = ie
	  end if
	end do

	uvmax = sqrt(uvmax)
	write(6,*) iemax,uvmax,htmax,htmin,iemin

	end
