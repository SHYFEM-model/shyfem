
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

c laplacian interpolation routines
c
c contents :
c
c revision log :
c
c 20.08.2003    ggu     new routines lapint and lapl_assemble
c 20.08.2003    ggu     new routine prepare_bc
c 30.10.2003    ggu     routine prepare_bc moved to laplap.f file
c
c******************************************************************

	subroutine lapint(rmat,zv,rzv,flag)

c solves linear system matrix for laplacian
c
c uses the following variables and arrays:
c
c	nkn,nel,mbw
c	nen3v
c	ev

	use basin, only : nkn,nel,ngr,mbw

	implicit none

c arguments
	real rmat(1)		!band matrix to be assembled
	real zv(1)		!solution to system on return
	real rzv(1)		!boundary conditions
	real flag		!flag to distinguish boundary condition

c common
c local
	integer i,k
	integer matdim
	integer ier
	real epseps

	epseps = 1.e-6
	matdim = nkn*(mbw+1)	!dimension of rmat - only for symmetric matrix

        do i=1,matdim
            rmat(i)=0.
        end do

        do k=1,nkn
            zv(k)=0.
        end do

	call lapl_assemble(rmat,zv,rzv,flag)

        call mchb(zv,rmat,nkn,1,mbw,1,epseps,ier)
        if(ier.ne.0) goto 99

	return
   99	continue
	write(6,*) 'Error in inverting matrix solving laplacian'
	write(6,*) 'ier : ',ier
	stop 'error stop : lapint'
	end

c******************************************************************

	subroutine lapl_assemble(rmat,zv,rzv,flag)

c assembles linear system matrix for laplacian

	use evgeom
	use basin

	implicit none

c arguments
	real rmat(1)		!band matrix to be assembled
	real zv(1)		!constant vector for z system
	real rzv(1)		!boundary conditions
	real flag		!flag to distinguish boundary condition
c common

	include 'param.h'
c local
	integer kn(3)
	integer ie,i,j,j1,j2,n,m,kk
	integer ngl
	real area,rw
	real hia(3,3),hik(3)
	real b(3),c(3)
c fucntion
	integer locsps

	ngl=nkn

	do ie=1,nel

	  do i=1,3
		kn(i)=nen3v(i,ie)
		b(i)=ev(i+3,ie)
		c(i)=ev(i+6,ie)
	  end do

	  area = 12. * ev(10,ie)

c set up local matrix

	  do n=1,3
	    do m=1,3
	      hia(n,m) = area * ( b(n)*b(m)+c(n)*c(m) )
	    end do
	    hik(n) = 0.
	  end do

c implement boundary conditions

	  do i=1,3
	    if(rzv(kn(i)).ne.flag) then
		rw=rzv(kn(i))
		j1=mod(i,3)+1
		j2=mod(i+1,3)+1
		hik(j1)=hik(j1)-rw*hia(j1,i)
		hik(j2)=hik(j2)-rw*hia(j2,i)
		hia(i,j1)=0.
		hia(i,j2)=0.
		hia(j1,i)=0.
		hia(j2,i)=0.
		hik(i)=rw*hia(i,i)
	    end if
	  end do

c in hia(i,j),hik(i),i,j=1,3 is system

          do i=1,3
            do j=1,3
              kk=locsps(kn(i),kn(j),ngl,mbw)
              if(kk.gt.0) then
                rmat(kk)=rmat(kk)+hia(i,j)
              end if
            end do
	    zv(kn(i))=zv(kn(i))+hik(i)
          end do

	end do

	end

c******************************************************************

