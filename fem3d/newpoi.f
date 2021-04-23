
!--------------------------------------------------------------------------
!
!    Copyright (C) 1991-1992,2001,2015,2018-2020  Georg Umgiesser
!    Copyright (C) 2015  William McKiver
!    Copyright (C) 2015  Debora Bellafiore
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

! module for solving poisson equation in 2d and 3d
!
! revision log :
!
! 18.02.1991	ggu	(from scratch)
! 04.06.1991	ggu	(c=(1) : friction term has been corrected)
! 01.10.1992	ggu	(staggered FE - completely restructured)
! 12.01.2001	ggu	solve for znv and not level difference (ZNEW)
! 15.12.2015	ggu&dbf&wmk	new version written from scratch
! 18.12.2015	ggu	changed VERS_7_3_17
! 26.04.2018	ggu	changed VERS_7_5_46
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	changed VERS_7_5_61
! 21.05.2019	ggu	changed VERS_7_5_62
! 16.02.2020    ggu     femtime eliminated
! 13.03.2021    clr&ggu adapted for petsc solver
c 23.04.2021    clr     alternative implementation to replace pragma directives
!
! notes :
!
! depth values are not used
! the equation is scaled vertically to be compatible with horizontal scale
!
!========================================================
	module poisson
!========================================================

	implicit none

	integer, save :: ipoiss = 0
	logical, save :: bpoi3d = .false.

!========================================================
	end module poisson
!========================================================

!******************************************************************

	subroutine poisson_init

	use mod_system
	use poisson

	implicit none

	real getpar

	ipoiss = nint(getpar('ipoiss'))

	bsys3d = .false.
	if( ipoiss /= 0 ) bsys3d = .true.

	end

!******************************************************************

	subroutine poisson_compute

	use basin
	use poisson

	implicit none

	if( ipoiss == 0 ) return

	call poisson_2d
	call poisson_3d

	stop 'end of poisson case'
	end

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine poisson_2d

	use basin
	use poisson

	implicit none

	real pvar(nkn)

	bpoi3d = .false.

	call poisson_set_obc(1,pvar)
	call poisson_2d_solve(pvar)
	call poisson_write(1,pvar)

	end

!******************************************************************

	subroutine poisson_2d_solve(pvar)

	use basin

	implicit none

	real pvar(nkn)
	real p0(nkn)

	p0 = 0.

	call system_init

	call poisson_2d_assemble(pvar)

	call system_solve(nkn,p0) 	!solves system matrix for pvar
	call system_get(nkn,pvar)	!copies solution to new pvar

	end

!******************************************************************

	subroutine poisson_2d_assemble(pvar)

! assembles linear system matrix
!
! vqv		flux boundary condition vector
!
! semi-implicit scheme for 3d model

	use mod_internal
	use mod_depth
	use evgeom
	use levels
	use basin
        use mod_zeta_system, only : kn,hia,hik

	implicit none

	real pvar(nkn)

	real drittl
	parameter (drittl=1./3.)

	include 'mkonst.h'
	include 'pkonst.h'
 
	integer ie,i,j,j1,j2,n,m,kk,l,k
	integer ngl
	integer ilevel
	integer lmax
	real aj,aj4,aj12
	real ht
	real h11,hh999
	real delta
	real amatr(3,3)
	real b(3),c(3),z(3)

	integer locsps,loclp,iround
	real getpar

!-------------------------------------------------------------
! loop over elements
!-------------------------------------------------------------

	do ie=1,nel

!	  ------------------------------------------------------
!	  initialize element values
!	  ------------------------------------------------------

	  aj=ev(10,ie)
          aj4=4.*aj
          aj12=12.*aj
	  do i=1,3
	    kk=nen3v(i,ie)
	    kn(i)=kk
	    b(i)=ev(i+3,ie)
	    c(i)=ev(i+6,ie)
	  end do

!	  ------------------------------------------------------
!	  set element matrix and RHS
!	  ------------------------------------------------------

	  do n=1,3
	    do m=1,3
	      hia(n,m) = -aj12 * ( b(n) * b(m) + c(n) * c(m) )
	    end do
	    hik(n) = 0.
	  end do

!	  ------------------------------------------------------
!	  boundary conditions
!	  ------------------------------------------------------

	  do i=1,3
	    if( pvar(kn(i)) .ne. flag ) then
	      j1=mod(i,3)+1
	      j2=mod(i+1,3)+1
              hia(i,i)=1.
              hia(i,j1)=0.
              hia(i,j2)=0.
              hik(i)=pvar(kn(i))
	    end if
	  end do

!	  ------------------------------------------------------
!	  in hia(i,j),hik(i),i,j=1,3 is system
!	  ------------------------------------------------------

	  !call system_assemble(ie,nkn,mbw,kn,hia,hik)
	  call system_assemble(ie)

	end do

!-------------------------------------------------------------
! end of loop over elements
!-------------------------------------------------------------

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------

	end

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine poisson_3d

	use basin
	use levels
	use poisson

	implicit none

	real pvar(nlvdi,nkn)

	bpoi3d = .true.

	call poisson_set_obc(nlvdi,pvar)
	call poisson_3d_solve(nlvdi,pvar)
	call poisson_write(nlvdi,pvar)

	end

!******************************************************************

	subroutine poisson_3d_solve(nlvdi,pvar)

	use basin
	use levels, only : nlv

	implicit none

	integer nlvdi
	real pvar(nlvdi,nkn)
	real p0(nlvdi,nkn)

	p0 = 0.

	call system_init

	call poisson_3d_assemble(pvar)

	call system_solve_3d(nkn,nlvdi,nlv,p0) !solves system matrix for pvar
	call system_get_3d(nkn,nlvdi,nlv,pvar) !copies solution to new pvar

	end

!******************************************************************

	subroutine poisson_3d_assemble(pvar)

! assembles linear system matrix
!
! vqv		flux boundary condition vector
!
! semi-implicit scheme for 3d model

	use mod_internal
	use mod_depth
	use mod_layer_thickness
	use evgeom
	use levels
	use basin

	implicit none

	real pvar(nlvdi,nkn)

	real drittl
	parameter (drittl=1./3.)

	include 'mkonst.h'
	include 'pkonst.h'
 
	integer kn(3)
	integer ie,i,j,j1,j2,n,m,kk,l,k,iii
	integer ngl
	integer ilevel
	integer lmax
	real aj,aj4,aj12
	real ht
	real h11,hh999
	real delta,r,hh
	real hia(3,3),hik(3),amatr(3,3)
	real hia3d(-1:+1,3,3)
	real b(3),c(3),z(3)
	double precision darea
	real dist

	real hd,hdm,hdp,rhm,rhp,rhc
	real hldaux(0:nlvdi+1)

	integer locsps,loclp,iround

	hldaux = 0.

!-------------------------------------------------------------
! scaling of vertical distances
!-------------------------------------------------------------

	darea = 0.
	do ie=1,nel
	  darea = darea + 12 * ev(10,ie)
	end do
	dist = sqrt(darea/nel)
	
!-------------------------------------------------------------
! loop over elements
!-------------------------------------------------------------

	do ie=1,nel

!	  ------------------------------------------------------
!	  initialize element values
!	  ------------------------------------------------------

	  aj=ev(10,ie)
          aj4=4.*aj
          aj12=12.*aj
	  do i=1,3
	    kk=nen3v(i,ie)
	    kn(i)=kk
	    b(i)=ev(i+3,ie)
	    c(i)=ev(i+6,ie)
	  end do

!	  ------------------------------------------------------
!	  set element matrix and RHS
!	  ------------------------------------------------------

	  lmax = ilhv(ie)
	  hldaux(1:lmax) = hdenv(1:lmax,ie)
	  hldaux(lmax+1:nlv) = 0.

	  do l=1,lmax

	  hia3d = 0.

	  hh = dist
	  hd = hh
	  hdm = hh
	  hdp = hh
	  !hd = hldaux(l)
	  !hdm = hldaux(l-1)
	  !hdp = hldaux(l+1)

	  rhm = 2. / ( hd + hdm )
	  if( l == 1 ) rhm = 0.
	  rhp = 2. / ( hd + hdp )
	  if( l == lmax ) rhp = 0.
	  rhc = rhm + rhp

	  do n=1,3
	    do m=1,3
	      hia3d(0,n,m) = -aj12 * hd * ( b(n) * b(m) + c(n) * c(m) )
	      if ( n .eq. m ) then
	        hia3d(0,n,m) = hia3d(0,n,m) - aj4 * rhc
	        hia3d(-1,n,m) = aj4 * rhm
	        hia3d(+1,n,m) = aj4 * rhp
              end if
	    end do
	    hik(n) = 0.
	  end do

!	  ------------------------------------------------------
!	  boundary conditions
!	  ------------------------------------------------------

	  do i=1,3
	    if( pvar(l,kn(i)) .ne. flag ) then
	      r = 1.
	      r = hia3d(0,i,i)
	      j1=mod(i,3)+1
	      j2=mod(i+1,3)+1
              hia3d(-1,i,i)=0.
              hia3d(+1,i,i)=0.
              hia3d(0,i,i)=r
              hia3d(0,i,j1)=0.
              hia3d(0,i,j2)=0.
              hik(i)=pvar(l,kn(i)) * r
	    end if
	  end do

!	  ------------------------------------------------------
!	  in hia(i,j),hik(i),i,j=1,3 is system
!	  ------------------------------------------------------

	  call system_assemble_3d(ie,l,nlv,kn,hia3d,hik)

	  end do

	end do

!-------------------------------------------------------------
! end of loop over elements
!-------------------------------------------------------------

	call system_adjust_matrix_3d

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------

	end

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine poisson_set_obc(nlvdi,pvar)

	use basin
	use poisson

	implicit none

	integer nlvdi
	real pvar(nlvdi,nkn)

	include 'mkonst.h'

	integer k,l
	real bnd

	if( ipoiss < 1 .or. ipoiss > 6 ) then
	  write(6,*) 'ipoiss = ',ipoiss
	  stop 'error stop poisson_set_obc'
	end if

	pvar = flag

	do k=1,nkn
	  bnd = 0.
	  if( ipoiss == 1 ) then
            if( ipv(k) .ge. 1 .and. ipv(k) .le. 9 ) then
		bnd = 10.
            else if( ipv(k) .ge. 3601 .and. ipv(k) .le. 3609 ) then
		bnd = 20.
	    end if
	  else if( ipoiss == 2 ) then
            if( ipv(k) .ge. 1 .and. ipv(k) .le. 9 ) then
		bnd = 10.
            else if( ipv(k) .ge. 217 .and. ipv(k) .le. 225 ) then
		bnd = 20.
	    end if
	  else if( ipoiss == 3 ) then
            if( ipv(k) .ge. 1 .and. ipv(k) .le. 2 ) then
		bnd = 10.
            else if( ipv(k) .ge. 5 .and. ipv(k) .le. 6 ) then
		bnd = 20.
	    end if
	  else if( ipoiss == 4 ) then
            if( ipv(k) .ge. 1 .and. ipv(k) .le. 2 ) then
		bnd = 10.
            else if( ipv(k) .ge. 7 .and. ipv(k) .le. 8 ) then
		bnd = 20.
	    end if
	  end if
	  if( bnd /= 0. ) pvar(:,k) = bnd
	end do

	do l=1,nlvdi
	  bnd = 0.
	  if( ipoiss == 5 ) then
	    if( l == 1 ) then
	      bnd = 10.
	      pvar(l,:) = bnd
	    else if( l == nlvdi ) then
	      bnd = 20.
	      pvar(l,:) = bnd
	    end if
	  else if( ipoiss == 6 ) then
	    if( l == 1 ) then
	      bnd = 10.
	      pvar(l,:) = bnd
	    else if( l == nlvdi ) then
	      bnd = 20.
	      pvar(l,:) = bnd
	    else
	      do k=1,nkn
	        bnd = 0.
                if( ipv(k) .ge. 1 .and. ipv(k) .le. 9 ) then
	  	  bnd = 10.
                else if( ipv(k) .ge. 217 .and. ipv(k) .le. 225 ) then
		  bnd = 20.
	        end if
	        if( bnd /= 0. ) pvar(l,k) = bnd
	      end do
	    end if
	  end if
	end do

	!write(6,*) 'poisson boundary: ',pvar
	do k=1,nkn
	  write(6,*) k
	  write(6,'(5f15.4)') pvar(:,k)
	end do

	end

!******************************************************************

	subroutine poisson_write(nlvdi,pvar)

	use basin
	use poisson

	implicit none

	integer nlvdi
	real pvar(nlvdi,nkn)

	integer k,l,iu
	integer, save :: iu2d = 0
	integer, save :: iu3d = 0

	l = (nlvdi+1)/2
	iu = 635
	if( bpoi3d ) iu = 666

	write(6,*) 'poisson writing: ',bpoi3d,nlvdi,l,iu,nkn

        do k = 1,nkn
          if ( xgv(k) .eq. 0.0 ) then
            write(iu,*) ygv(k),pvar(l,k),xgv(k) 
          end if
        end do
	flush(iu)

	do l=1,nlvdi
	  write(iu+1,'(i6,5f12.4)') l,(pvar(l,k),k=1,nkn,nkn/5)
	end do

	if( nlvdi == 0 ) then
	  call conwrite(iu2d,'.p2d',1,10,nlvdi,pvar)
	  flush(iu2d)
	else
	  call conwrite(iu3d,'.p3d',1,10,nlvdi,pvar)
	  flush(iu3d)
	end if

	end

!******************************************************************
!******************************************************************
!******************************************************************

