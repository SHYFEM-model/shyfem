
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

c utility routines for flux computations
c
c contents :
c
c subroutine flxscs(n,kflux,iflux,az,fluxes)	flux through sections
c subroutine flxsec(n,kflux,iflux,az,fluxes)	flux through section
c
c subroutine flxini				initializes flux routines
c subroutine flx_init(kfluxm,kflux,nsect,iflux)	sets up array iflux
c subroutine flxinf(m,kflux,iflux)		sets up one info structure
c function igtnsc(k1,k2)			gets number of internal section
c
c revision log :
c
c 09.05.2013	ggu	separated from subflxa.f
c 14.05.2013	ggu	deleted error check between 2d and 3d computation
c 13.06.2013	ggu	changed VERS_6_1_65
c 19.01.2015	ggu	changed VERS_7_1_2
c 19.01.2015	ggu	changed VERS_7_1_3
c 10.07.2015	ggu	changed VERS_7_1_50
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 16.12.2015	ggu	changed VERS_7_3_16
c 18.12.2015	ggu	changed VERS_7_3_17
c 15.04.2016	ggu	changed VERS_7_5_8
c 26.10.2016	ccf	bug fix in flxsec
c 12.01.2017	ggu	changed VERS_7_5_21
c 30.03.2017	ggu	changed accumulator to time step dt, not number of calls
c 04.02.2018	ggu	new routines with accumulator in double
c 22.02.2018	ggu	changed VERS_7_5_42
c 03.04.2018	ggu	changed VERS_7_5_43
c 16.02.2019	ggu	changed VERS_7_5_60
c
c notes :
c
c These routines can also be used internally to compute the flux
c over various sections. The following calling sequence must be respected:
c
c call flx_init(kfluxm,kflux,nsect,iflux)		initializes iflux
c
c call flxscs(kfluxm,kflux,iflux,az,fluxes) computes fluxes 
c
c Initialization can be done anytime.
c
c******************************************************************
c******************************************************************
c******************************************************************
c******************************************************************
c******************************************************************


c******************************************************************
c******************************************************************
c******************************************************************
c******************************************************************
c******************************************************************


c******************************************************************
c******************************************************************
c******************************************************************
c******************************************************************
c******************************************************************

	subroutine fluxes_init(nlvddi,nsect,nlayers,tr,masst)

c initializes nr and masst

	implicit none

	integer nlvddi,nsect
	integer nlayers(nsect)
	real tr
	real masst(0:nlvddi,3,nsect)

	integer i,l,lmax

        tr = 0.
	masst = 0.

	end

c******************************************************************

	subroutine fluxes_accum(nlvddi,nsect,nlayers,dt,tr,masst,fluxes)

c accumulates fluxes into masst

	implicit none

	integer nlvddi,nsect
	integer nlayers(nsect)
	real dt
	real tr
	real masst(0:nlvddi,3,nsect)
	real fluxes(0:nlvddi,3,nsect)

        tr = tr + dt
	masst = masst + fluxes * dt

	end

c******************************************************************

	subroutine fluxes_aver(nlvddi,nsect,nlayers,tr,masst,fluxes)

c averages masst and puts result into fluxes

	implicit none

	integer nlvddi,nsect
	integer nlayers(nsect)
	real tr
	real masst(0:nlvddi,3,nsect)
	real fluxes(0:nlvddi,3,nsect)

	if( tr == 0. ) return
        fluxes = masst / tr

	end

c******************************************************************
c******************************************************************
c******************************************************************

	subroutine fluxes_init_d(nlvddi,nsect,nlayers,tr,masst)

c initializes nr and masst

	implicit none

	integer nlvddi,nsect
	integer nlayers(nsect)
	double precision tr
	double precision masst(0:nlvddi,3,nsect)

	integer i,l,lmax

        tr = 0.
	masst = 0.

	end

c******************************************************************

	subroutine fluxes_accum_d(nlvddi,nsect,nlayers,dt,tr,masst,fluxes)

c accumulates fluxes into masst

	implicit none

	integer nlvddi,nsect
	integer nlayers(nsect)
	real dt
	double precision tr
	double precision masst(0:nlvddi,3,nsect)
	real fluxes(0:nlvddi,3,nsect)

        tr = tr + dt
	masst = masst + fluxes * dt

	end

c******************************************************************

	subroutine fluxes_aver_d(nlvddi,nsect,nlayers,tr,masst,fluxes)

c averages masst and puts result into fluxes

	implicit none

	integer nlvddi,nsect
	integer nlayers(nsect)
	double precision tr
	double precision masst(0:nlvddi,3,nsect)
	real fluxes(0:nlvddi,3,nsect)

	if( tr == 0. ) return
        fluxes = masst / tr

	end

c******************************************************************
c******************************************************************
c******************************************************************
c******************************************************************
c******************************************************************

	subroutine flxscs(kfluxm,kflux,iflux,az,fluxes,is,scalar)

c computes flux through all sections and returns them in fluxes
c
c flux are divided into total, positive and negative

	use levels, only : nlvdi,nlv

	implicit none

	integer kfluxm
	integer kflux(kfluxm)
	integer iflux(3,kfluxm)
	real az
	real fluxes(0:nlvdi,3,*)	!computed fluxes (return)
	integer is			!type of scalar (0=mass)
	real scalar(nlvdi,*)

	integer nnode,ifirst,ilast,ntotal
	integer ns
	logical nextline

	nnode = 0
	ns = 0

	do while( nextline(kflux,kfluxm,nnode,ifirst,ilast) )
	  ns = ns + 1
	  ntotal = ilast - ifirst + 1
	  !write(66,*) 'section ',ns,ntotal,is
	  call flxsec(ntotal,kflux(ifirst),iflux(1,ifirst),az
     +				,fluxes(0,1,ns),is,scalar)
	end do

	end

c******************************************************************

	subroutine flxsec(n,kflux,iflux,az,fluxes,is,scalar)

c computes flux through one section and returns it in fluxes

	use levels, only : nlvdi,nlv

	implicit none

	integer n
	integer kflux(n)
	integer iflux(3,n)
	real az
	real fluxes(0:nlvdi,3)		!computed fluxes (return)
	integer is			!type of scalar (0=mass)
	real scalar(nlvdi,*)

	integer i,k,l,lkmax
	integer istype,iafter,ibefor
	real port,port2d,sport
	real flux(nlvdi)

	fluxes = 0.

	do i=1,n
		k = kflux(i)
		istype = iflux(1,i)
		ibefor = iflux(2,i)
		iafter = iflux(3,i)

		call flx2d(k,ibefor,iafter,istype,az,port2d)

		call flx3d(k,ibefor,iafter,istype,az,lkmax,flux)

		do l=1,lkmax
		  port = flux(l)
		  sport = port
		  if( is .gt. 0 ) sport = port * scalar(l,k)	!not mass

		  fluxes(l,1) = fluxes(l,1) + sport
		  if( port .gt. 0. ) then
		    fluxes(l,2) = fluxes(l,2) + sport
		  else
		    fluxes(l,3) = fluxes(l,3) - sport
		  end if
		end do
	end do

! compute vertical integrated fluxes

        fluxes(0,1) = sum(fluxes(1:nlvdi,1))
        fluxes(0,2) = sum(fluxes(1:nlvdi,2))
        fluxes(0,3) = sum(fluxes(1:nlvdi,3))

	end
	  
c******************************************************************
c******************************************************************
c******************************************************************

	subroutine flx_init(kfluxm,kflux,nsect,iflux)

c sets up array iflux

	implicit none

        integer kfluxm		!total number of nodes in kflux
        integer kflux(kfluxm)	!nodes in sections
	integer nsect		!number of section (return)
	integer iflux(3,*)	!internal array for flux computation (return)

	integer ifirst,ilast,nnode,ntotal

	integer klineck
	logical nextline

c----------------------------------------------------------
c check nodes for compatibility
c----------------------------------------------------------

	nsect = klineck(kfluxm,kflux)

	if( nsect .lt. 0 ) then
	  write(6,*) 'errors setting up fluxes ',nsect
	  stop 'error stop : flx_init'
	end if

c----------------------------------------------------------
c now set info structure for sections
c----------------------------------------------------------

	nnode = 0

	do while( nextline(kflux,kfluxm,nnode,ifirst,ilast) )
	  ntotal = ilast - ifirst + 1
c	  write(6,*) kfluxm,nnode,ifirst,ilast,ntotal
	  call flxinf(ntotal,kflux(ifirst),iflux(1,ifirst))
	end do

c----------------------------------------------------------
c end of routine
c----------------------------------------------------------

	end

c******************************************************************

	subroutine flxinf(n,kflux,iflux)

c sets up info structure iflux(3,1) for one section

	implicit none

	integer n
	integer kflux(n)
	integer iflux(3,n)

	integer i,k
	integer ktype
	integer kafter,kbefor

	integer igtnsc,flxtype

	do i=1,n
	  k = kflux(i)
	  ktype = flxtype(k)

	  iflux(1,i) = ktype

	  kbefor = 0
	  if( i .ne. 1 ) kbefor = kflux(i-1)
	  kafter = 0
	  if( i .ne. n ) kafter = kflux(i+1)

	  iflux(2,i) = igtnsc(k,kbefor)
	  iflux(3,i) = igtnsc(k,kafter)
	end do

	end

c******************************************************************

	function igtnsc(k1,k2)

c gets number of internal section in link index of k1

	use mod_geom

	implicit none

	integer igtnsc
	integer k1,k2
	integer elems(maxlnk)

	integer k,i,n,ie

	integer knext,kbhnd

	call get_elems_around(k1,maxlnk,n,elems)

c	deal with no section given

	igtnsc = 0
	if( k2 .le. 0 ) return

c	section in the middle

	do i=1,n
	  igtnsc = i
	  ie = elems(i)
	  k = knext(k1,ie)
	  if( k .eq. k2 ) return
	end do

c	in case we are on the boundary

	i = n
	igtnsc = igtnsc + 1
	ie = elems(i)
	k = kbhnd(k1,ie)
	if( k .eq. k2 ) return

c	no node found

	write(6,*) k1,k2
	write(6,*) k1,n
	write(6,*) (elems(i),i=1,n)
	call get_elems_around(k2,maxlnk,n,elems)
	write(6,*) k2,n
	write(6,*) (elems(i),i=1,n)
	stop 'error stop igtnsc: internal error (2)'
	end
	      
c**********************************************************************
c**********************************************************************
c**********************************************************************
c**********************************************************************
c**********************************************************************

	subroutine get_nlayers(kfluxm,kflux,nlayers,nlmax)

c computes maximum numer of layers for sections

	use levels

	implicit none

	integer kfluxm
	integer kflux(*)
	integer nlayers(*)	!total number of layers for sections (return)
	integer nlmax		!maximum layers for all sections (return)

	integer ns
	integer nnode,ifirst,ilast
	integer i,k,l,lmax

	logical nextline

	ns = 0
	nnode = 0
	nlmax = 0

	do while( nextline(kflux,kfluxm,nnode,ifirst,ilast) )
	  ns = ns + 1
	  lmax = 0
	  do i=ifirst,ilast
	    k = kflux(i)
	    l = ilhkv(k)
	    lmax = max(lmax,l)
	  end do
	  nlayers(ns) = lmax
	  nlmax = max(nlmax,lmax)
	end do

	end

c**********************************************************************
c**********************************************************************
c**********************************************************************
c**********************************************************************
c**********************************************************************

