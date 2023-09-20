
!--------------------------------------------------------------------------
!
!    Copyright (C) 2017,2019-2020  Georg Umgiesser
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
! 04.11.2017	ggu	changed VERS_7_5_34
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62
! 28.01.2020	ggu	new code to transfrom trans2vel on element, vorticity
! 02.02.2020	ggu	update for vorticity
! 09.05.2023    lrp     introduce top layer index variable
! 05.06.2023    lrp     introduce z-star

!***************************************************************

        subroutine prepare_hydro(bvel,nndim,cv3all,znv,uprv,vprv)

        use basin
        use levels
        use mod_depth

        implicit none

        logical bvel
        integer nndim
        real cv3all(nlvdi,nndim,0:4)
        real znv(nkn)
        real uprv(nlvdi,nkn)
        real vprv(nlvdi,nkn)

        real, allocatable :: zenv(:)
        real, allocatable :: uv(:,:)
        real, allocatable :: vv(:,:)

        allocate(zenv(3*nel))
        allocate(uv(nlvdi,nel))
        allocate(vv(nlvdi,nel))

        znv(1:nkn)     = cv3all(1,1:nkn,1)
        zenv(1:3*nel)  = cv3all(1,1:3*nel,2)
        uv(:,1:nel)    = cv3all(:,1:nel,3)
        vv(:,1:nel)    = cv3all(:,1:nel,4)

        call shy_transp2vel(bvel,nel,nkn,nlv,nlvdi,hev,zenv,nen3v
     +                          ,ilhv,hlv,uv,vv
     +                          ,uprv,vprv)

        deallocate(zenv,uv,vv)

        end

!***************************************************************

        subroutine prepare_hydro_e(bvel,nndim,cv3all,zev,ue3v,ve3v)

        use basin
        use levels
        use mod_depth

        implicit none

        logical bvel
        integer nndim
        real cv3all(nlvdi,nndim,0:4)
        real zev(nel)
        real ue3v(nlvdi,nel)
        real ve3v(nlvdi,nel)

        real, allocatable :: znv(:)
        real, allocatable :: zenv(:)
        real, allocatable :: uv(:,:)
        real, allocatable :: vv(:,:)

        allocate(znv(nkn))
        allocate(zenv(3*nel))
        allocate(uv(nlvdi,nel))
        allocate(vv(nlvdi,nel))

        znv(1:nkn)     = cv3all(1,1:nkn,1)
        zenv(1:3*nel)  = cv3all(1,1:3*nel,2)
        uv(:,1:nel)    = cv3all(:,1:nel,3)
        vv(:,1:nel)    = cv3all(:,1:nel,4)

        call shy_transp2vel_e(bvel,nel,nkn,nlv,nlvdi,hev,zenv,nen3v
     +                          ,ilhv,hlv,uv,vv
     +                          ,zev,ue3v,ve3v)

        deallocate(znv,zenv,uv,vv)

        end

!***************************************************************

        subroutine convert_to_speed(uprv,vprv,sv,dv)

        use basin
        use levels

        implicit none

        real uprv(nlvdi,nkn)
        real vprv(nlvdi,nkn)
        real sv(nlvdi,nkn)
        real dv(nlvdi,nkn)

        integer k,lmax,l
        real u,v,s,d

        do k=1,nkn
          lmax = ilhkv(k)
          do l=1,lmax
            u = uprv(l,k)
            v = vprv(l,k)
            call c2p_ocean(u,v,s,d)   !d is ocean convention
            sv(l,k) = s
            dv(l,k) = d
          end do
        end do

        end

!***************************************************************

        subroutine shy_transp2vel(bvel,nel,nkn,nlv,nlvddi
     +				,hev,zenv,nen3v
     +                          ,ilhv,hlv,utlnv,vtlnv
     +                          ,uprv,vprv)

! transforms transports at elements to velocities at nodes

        implicit none

	logical bvel			!if true compute velocities
        integer nel
        integer nkn
        integer nlv
        integer nlvddi
        real hev(nel)
        real zenv(3,nel)
        integer nen3v(3,nel)
        integer ilhv(nel)
        real hlv(nlvddi)
        real utlnv(nlvddi,nel)
        real vtlnv(nlvddi,nel)
        real uprv(nlvddi,nkn)
        real vprv(nlvddi,nkn)

        real weight(nlvddi,nkn)         !aux variable for weights
        real hl(nlvddi)                 !aux variable for real level thickness

        logical bsigma
        integer ie,ii,k,l,lmax,lmin,nsigma,nadapt,nlvaux
        real hmed,u,v,area,zeta,zmin
        real hsigma,hadapt

        real area_elem

        call get_sigma_info(nlvaux,nsigma,hsigma)
        if( nlvaux .gt. nlvddi ) stop 'error stop transp2vel: nlvddi'
        bsigma = nsigma .gt. 0

	weight = 0.
	uprv = 0.
	vprv = 0.
	hl = 1.		!in case of transports

        do ie=1,nel

          area = area_elem(ie)
          lmax = ilhv(ie)
	  lmin = 1
	  if( bvel ) then
	    zeta = sum(zenv(:,ie)) / 3.	!average of zeta on element
	    zmin = minval(zenv(:,ie))   !min: z-adapt coords works with zmin
	    call compute_zadapt_info(zmin,hlv,nsigma,lmax,
     +			             lmin,nadapt,hadapt)
	    call get_layer_thickness(lmax,lmin,nsigma,nadapt,
     +				     hsigma,hadapt,zeta,hev(ie),hlv,hl)
	  end if

          do l=1,lmax
            hmed = hl(l)
	    if (hmed .gt. 0.) then
              u = utlnv(l,ie) / hmed
              v = vtlnv(l,ie) / hmed
	    else
	      u = 0.
	      v = 0.
      	    end if
            do ii=1,3
              k = nen3v(ii,ie)
              uprv(l,k) = uprv(l,k) + area * u
              vprv(l,k) = vprv(l,k) + area * v
              weight(l,k) = weight(l,k) + area
            end do
          end do
        end do

	where( weight > 0. )
	  uprv = uprv / weight
	  vprv = vprv / weight
	end where

        end

!***************************************************************

        subroutine shy_transp2vel_e(bvel,nel,nkn,nlv,nlvddi
     +				,hev,zenv,nen3v
     +                          ,ilhv,hlv,utlnv,vtlnv
     +                          ,zev,ue3v,ve3v)

! transforms transports at elements to velocities at elements

        implicit none

	logical bvel			!if true compute velocities
        integer nel
        integer nkn
        integer nlv
        integer nlvddi
        real hev(nel)
        real zenv(3,nel)
        integer nen3v(3,nel)
        integer ilhv(nel)
        real hlv(nlvddi)
        real utlnv(nlvddi,nel)
        real vtlnv(nlvddi,nel)
        real zev(nel)
        real ue3v(nlvddi,nel)
        real ve3v(nlvddi,nel)

        real hl(nlvddi)                 !aux variable for real level thickness

        logical bsigma
        integer ie,ii,k,l,lmin,lmax,nsigma,nadapt,nlvaux
        real hmed,zeta,zmin
        real hsigma,hadapt

        call get_sigma_info(nlvaux,nsigma,hsigma)
        if( nlvaux .gt. nlvddi ) stop 'error stop transp2vel: nlvddi'
        bsigma = nsigma .gt. 0

	ue3v = 0.
	ve3v = 0.
	hl = 1.		!in case of transports

        do ie=1,nel

          lmax = ilhv(ie)
	  lmin = 1
	  zeta = sum(zenv(:,ie)) / 3.	!average of zeta on element
          zmin = minval(zenv(:,ie))     !min: z-adapt coords works with zmin	  
	  if( bvel ) then
	    call compute_zadapt_info(zmin,hlv,nsigma,lmax,
     +			             lmin,nadapt,hadapt)
	    call get_layer_thickness(lmax,lmin,nsigma,nadapt,
     +				     hsigma,hadapt,zeta,hev(ie),hlv,hl)
	  end if

	  zev(ie) = zeta
          do l=1,lmax
            hmed = hl(l)
            if (hmed .gt. 0.) then
              ue3v(l,ie) = utlnv(l,ie) / hmed
              ve3v(l,ie) = vtlnv(l,ie) / hmed
            else
              ue3v(l,ie) = 0.
              ve3v(l,ie) = 0.		    
            end if
          end do
        end do

        end

!******************************************************************

	subroutine compute_vorticity(nndim,cv3all,cv3)

        use basin
        use levels
        use mod_depth
        use evgeom

	implicit none

	integer, parameter :: nvar = 4
	logical, parameter :: bnode = .true.	!use nodal values to compute

	integer nndim
	real cv3all(nlvdi,nndim,0:nvar)
	real cv3(nlvdi,nkn)

	integer ie,ii,lmax,l,k
	real b(3),c(3),u,v,aj,cm

        real, allocatable :: aux(:,:)
        real, allocatable :: zv(:)
        real, allocatable :: u3v(:,:)
        real, allocatable :: v3v(:,:)

        allocate(aux(nlvdi,nkn))

	if( bnode ) then
          allocate(zv(nkn))
          allocate(u3v(nlvdi,nkn))
          allocate(v3v(nlvdi,nkn))
	  call prepare_hydro(.true.,nndim,cv3all,zv,u3v,v3v)
	else
          allocate(zv(nel))
          allocate(u3v(nlvdi,nel))
          allocate(v3v(nlvdi,nel))
          call prepare_hydro_e(.true.,nndim,cv3all,zv,u3v,v3v)
	end if

	cv3 = 0.
	aux = 0.

        do ie=1,nel
          aj = ev(10,ie)
	  b(:) = ev(4:6,ie)
	  c(:) = ev(7:9,ie)
          lmax = ilhv(ie)
	  if( bnode ) then
           do l=1,lmax
	    cm = 0.
	    do ii=1,3
              k = nen3v(ii,ie)
	      cm = cm + v3v(l,k)*b(ii) - u3v(l,k)*c(ii)
	    end do
	    do ii=1,3
              k = nen3v(ii,ie)
              cv3(l,k) = cv3(l,k) + aj*cm
              aux(l,k) = aux(l,k) + aj
	    end do
	   end do
	  else
           do l=1,lmax
            u = u3v(l,ie)
            v = v3v(l,ie)
            do ii=1,3
              k = nen3v(ii,ie)
              cv3(l,k) = cv3(l,k) + 3.*aj*(u*c(ii)-v*b(ii))
              aux(l,k) = aux(l,k) + aj
            end do
           end do
	  end if
        end do

        where ( aux > 0. ) cv3 = cv3 / aux

	end

!******************************************************************

