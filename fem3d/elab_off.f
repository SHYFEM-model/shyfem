
!--------------------------------------------------------------------------
!
!    Copyright (C) 1991,1998,2018-2019  Georg Umgiesser
!    Copyright (C) 2018  Christian Ferrarin
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

! utility routines for shyelab: elab_off
!
! revision log :
!
! 27.08.1991	ggu	(from scratch) (make_vertical_velocity_off)
! 14.08.1998	ggu	w=0 at open boundary nodes (make_vertical_velocity_off)
! 20.08.1998	ggu	some documentation (make_vertical_velocity_off)
! 24.05.2018	ccf	written from scratch
! 06.07.2018	ggu	changed VERS_7_5_48
! 12.11.2018	ggu	linear arrays introduced
! 18.12.2018	ggu	changed VERS_7_5_52
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	changed VERS_7_5_61
! 21.05.2019	ggu	changed VERS_7_5_62
!
!***************************************************************
!***************************************************************
!***************************************************************

        subroutine off_output_hydro(idout,dtime,nndim,cv3all)

! write off records from shy (hydro only)

        use basin
        use levels

        implicit none

        integer, intent(in)		:: idout
        double precision, intent(in)	:: dtime
        integer, intent(in)		:: nndim
        real, intent(in) 		:: cv3all(nlvdi,nndim,0:4)

        integer it
	integer ii,ie,l,k
	integer i,nlin,nlink,nline

        double precision, allocatable :: ut(:,:)
        double precision, allocatable :: vt(:,:)
        double precision, allocatable :: ze(:)
        double precision, allocatable :: wn(:,:)
        double precision, allocatable :: wnaux(:,:)
        double precision, allocatable :: zn(:)
        double precision, allocatable :: sn(:,:)
        double precision, allocatable :: tn(:,:)
        double precision, allocatable :: rlin(:)

! allocate arrays
        allocate(ut(nlvdi,nel))
        allocate(vt(nlvdi,nel))
        allocate(ze(3*nel))
        allocate(wn(0:nlvdi,nkn))
        allocate(wnaux(nlvdi,nkn))
        allocate(zn(nkn))
        allocate(sn(nlvdi,nkn))
        allocate(tn(nlvdi,nkn))

! assign values
 	it = int(dtime)
        zn(1:nkn)    = cv3all(1,1:nkn,1)
        ze(1:3*nel)  = cv3all(1,1:3*nel,2)
        ut(:,1:nel)  = cv3all(:,1:nel,3)
        vt(:,1:nel)  = cv3all(:,1:nel,4)
	call make_vertical_velocity_off(ut,vt,wn)
	wnaux(1:nlvdi,:) = wn(1:nlvdi,:)
	sn = 0.
	tn = 0.

! set up linear arrays
	call count_linear(nlvdi,nkn,1,ilhkv,nlink)
	nlink = nlink + nkn				!account for wn(0:...)
	call count_linear(nlvdi,nel,1,ilhv,nline)
	allocate(rlin(max(nlink,nline)))

! write to output file
        write(idout) it,nkn,nel,3
        write(idout) (ilhv(ie),ie=1,nel)
        write(idout) (ilhkv(k),k=1,nkn)

	nlin = nline
	call dvals2linear(nlvdi,nel,1,ilhv,ut,rlin,nlin)
        write(idout) (rlin(i),i=1,nlin)
        !write(idout) ((ut(l,ie),l=1,ilhv(ie)),ie=1,nel)
	nlin = nline
	call dvals2linear(nlvdi,nel,1,ilhv,vt,rlin,nlin)
        write(idout) (rlin(i),i=1,nlin)
        !write(idout) ((vt(l,ie),l=1,ilhv(ie)),ie=1,nel)

        write(idout) (ze(ii),ii=1,3*nel)
	nlin = nlink
	call dvals2linear(nlvdi,nkn,1,ilhkv,wnaux,rlin,nlin)
        write(idout) (rlin(i),i=1,nlin)
        !write(idout) ((wn(l,k),l=1,ilhkv(k)),k=1,nkn)
        write(idout) (zn(k),k=1,nkn)

	nlin = nlink
	call dvals2linear(nlvdi,nkn,1,ilhkv,sn,rlin,nlin)
        write(idout) (rlin(i),i=1,nlin)
        !write(idout) ((sn(l,k),l=1,ilhkv(k)),k=1,nkn)
	nlin = nlink
	call dvals2linear(nlvdi,nkn,1,ilhkv,tn,rlin,nlin)
        write(idout) (rlin(i),i=1,nlin)
        !write(idout) ((tn(l,k),l=1,ilhkv(k)),k=1,nkn)

	end subroutine off_output_hydro

c******************************************************************

        subroutine make_vertical_velocity_off(utlnv,vtlnv,wlnv)

c computes vertical velocities
c
c from sp256w in new3di.F
c
c velocities are computed on S/T points (top and bottom of layer)
c bottom velocity of the whole column is assumed to be 0
c -> maybe change this
c
c computes volume differences and from these computes vertical
c velocities at every time step so that the computed velocities
c satisfy the continuity equation for every single time step
c
c wlnv is computed horizontally at a node and vertically
c it is at the center of the layer -> there are nlv velocities
c computed
c
c b,c are 1/m, (phi is dimensionless)
c aj is m**2
c utlnv... is m**2/s
c dvol is in m**3/s
c vv is m**2 (area)
c
c wlnv is first used to accumulate volume difference -> dvol
c at the end it receives the vertical velocity
c
c wlnv (dvol)   aux array for volume difference
c vv            aux array for area

        use evgeom
        use levels
        use basin

        implicit none

        double precision utlnv(nlvdi,nel)
        double precision vtlnv(nlvdi,nel)
        double precision wlnv(0:nlvdi,nkn)

        integer k,ie,ii,kk,l
        integer ilevel
        integer inwater
        double precision aj,wbot,wtop,ff

	double precision, allocatable :: wauxv(:,:)

	allocate(wauxv(0:nlvdi,nkn))

        wauxv = 0.
        wlnv = 0.

c compute difference of velocities for each layer
c
c f(ii) > 0 ==> flux into node ii

        inwater = 0

        do ie=1,nel
         !if( bwater(ie) ) then           !FIXME        !not working
          inwater = inwater + 1
          aj=4.*ev(10,ie)               !area of triangle / 3
          ilevel = ilhv(ie)
          do l=1,ilevel
            do ii=1,3
                kk=nen3v(ii,ie)
                ff = utlnv(l,ie)*ev(ii+3,ie) + vtlnv(l,ie)*ev(ii+6,ie)
                wlnv(l,kk) = wlnv(l,kk) + 3. * aj * ff
                wauxv(l,kk)=wauxv(l,kk)+aj
            end do
          end do
         !end if
        end do

c from vel difference get absolute velocity (w_bottom = 0)
c       -> wlnv(nlv,k) is already in place !
c       -> wlnv(nlv,k) = 0 + wlnv(nlv,k)
c w of bottom of last layer must be 0 ! -> shift everything up
c wlnv(nlv,k) is always 0
c
c dividing dvol(m**3/s) by area (wauxv) gives vertical velocity

        do k=1,nkn
          wbot = 0.
          do l=nlv,1,-1
            wtop = wlnv(l,k)
            wlnv(l,k) = wbot
            wbot = wbot + wtop
          end do
          wlnv(0,k) = wbot
        end do

        do k=1,nkn
          do l=1,nlv
            if( wauxv(l,k) .gt. 0. ) then
              wlnv(l-1,k) = wlnv(l-1,k) / wauxv(l,k)
            end if
          end do
        end do

c set w to zero at open boundary nodes (new 14.08.1998)
c
c FIXME -> only for ibtyp = 1,2 !!!!

c        do k=1,nkn
c          if( inodv(k) .gt. 0 ) then    !open boundary node
c            do l=0,nlv
c               wlnv(l,k) = 0.
c            end do
c          end if
c        end do

        return
        end

c******************************************************************
