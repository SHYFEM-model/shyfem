
!--------------------------------------------------------------------------
!
!    Copyright (C) 2009  Georg Umgiesser
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

c utilitiy routines for section plot (velocities)
c
c revision log :
c
c 13.10.2009    ggu     routines written from scratch
c
c*******************************************************************

	subroutine prepare_vel(pp3)

	use mod_hydro_print
	use mod_hydro_plot
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	real pp3(nlvdi,nkn)

	integer k,l
	real u,v,w
	real href

	call make_vertical_velocity	!compute wlnv from utlnv,vtlnv

	href = 0.
	call mkht3(nlvdi,het3v,href)

	call make_vel_from_tra(het3v)
	call vel_to_node
	
	do k=1,nkn
	  do l=1,nlv
	    u = uprv(l,k)
	    v = vprv(l,k)
	    w = wprv(l,k)
	    pp3(l,k) = sqrt( u*u + v*v + w*w )
	  end do
	end do

	end

c*******************************************************************

	subroutine make_vel_from_tra(het3v)

c from transports to velocities (on elements)

	use mod_hydro_vel
	use mod_hydro
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

        real het3v(nlvdi,nel)		!layer depth at elements

	integer ie,l,lmax
	real h,rh

	do ie=1,nel
	  lmax = ilhv(ie)
	  do l=1,lmax
	    h = het3v(l,ie)
	    if( h .le. 0. ) goto 99
	    rh = 1. / h
	    ulnv(l,ie) = utlnv(l,ie) * rh
	    vlnv(l,ie) = vtlnv(l,ie) * rh
	  end do
	end do

	return
   99	continue
	write(6,*) 'ie,l,h: ',ie,l,h
	stop 'error stop make_vel_from_tra: zero depth'
	end

c*******************************************************************

	subroutine vel_to_node

c transfers velocities at elements to nodes 
c and vertical velocities to center of layer

	use mod_hydro_print
	use mod_hydro_vel
	use evgeom
	use levels
	use basin

	implicit none

	integer ie,ii,k,l,lmax
	real aj

        do k=1,nkn
          do l=1,nlv
            uprv(l,k)=0.
            vprv(l,k)=0.
            wprv(l,k)=0.
          end do
        end do

        do ie=1,nel
          aj=ev(10,ie)
	  lmax = ilhv(ie)
          do ii=1,3
            k=nen3v(ii,ie)
	    do l=1,lmax
              wprv(l,k)=wprv(l,k)+aj
              uprv(l,k)=uprv(l,k)+aj*ulnv(l,ie)
              vprv(l,k)=vprv(l,k)+aj*vlnv(l,ie)
            end do
	  end do
	end do

        do k=1,nkn
	  do l=1,nlv
            if(wprv(l,k).gt.0.) then
              uprv(l,k)=uprv(l,k)/wprv(l,k)
              vprv(l,k)=vprv(l,k)/wprv(l,k)
            end if
          end do
        end do

        do k=1,nkn
          do l=1,nlv
            wprv(l,k)=0.5*(wlnv(l,k)+wlnv(l-1,k))
          end do
        end do

	end

c*******************************************************************

