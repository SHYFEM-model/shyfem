
!--------------------------------------------------------------------------
!
!    Copyright (C) 2010,2019  Georg Umgiesser
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

c header for closing sections
c
c---------------------------------------------------------------------

c	dimension of vectors

! revision log :
!
! 23.03.2010	ggu	changed v6.1.1
! 16.02.2019	ggu	changed VERS_7_5_60

c---------------------------------------------------------------------


	integer ipcdim,ivcdim
	parameter (ipcdim=500,ivcdim=500)

c---------------------------------------------------------------------

c	declaration of vectors

	integer ipccv(ipcdim)
	integer ivccv(ivcdim)
	common /ipccv/ipccv
	common /ivccv/ivccv
	real rpccv(ipcdim)
	real rvccv(ivcdim)
	equivalence (ipccv(1),rpccv(1))
	equivalence (ivccv(1),rvccv(1))

	save /ipccv/ , /ivccv/

c---------------------------------------------------------------------

c	size of header and section

	integer lhead,lsect

	parameter (lhead=7)		!size of header
	parameter (lsect=25)		!size of one section

c	pointer into header (must be called with isect = 0)

	integer lipdim,livdim,lipful,livful,lnsect

	parameter (lipdim=lsect-lhead+1)	!dimension of ipccv
	parameter (livdim=lsect-lhead+2)	!dimension ov ivccv
	parameter (lipful=lsect-lhead+3)	!filling of ipccv
	parameter (livful=lsect-lhead+4)	!filling of ivccv
	parameter (lnsect=lsect-lhead+7)	!number of sections

c	pointer into sections (must be called with isect > 0)

	integer lnkboc,lkboc,lniboc,liboc
	integer lhboc,lnitb,litb,lkout,lkin
	integer lkref,lkdir,lisoft,lzdate,lvdate
	integer lmnstp,liclos,listp
	integer lscal,lhref,liact,limode
	integer lzdiff,lflux,libnd,libndz

	parameter (lnkboc=1)
	parameter (lkboc=2)
	parameter (lniboc=3)
	parameter (liboc=4)
	parameter (lhboc=5)
	parameter (lnitb=6)
	parameter (litb=7)
	parameter (lkout=8)
	parameter (lkin=9)
	parameter (lkref=10)
	parameter (lkdir=11)
	parameter (lisoft=12)
	parameter (lzdate=13)
	parameter (lvdate=14)
	parameter (lmnstp=15)
	parameter (liclos=16)
	parameter (listp=17)
	parameter (lscal=18)
	parameter (lhref=19)
	parameter (liact=20)
	parameter (limode=21)
	parameter (lzdiff=22)
	parameter (lflux=23)
	parameter (libnd=24)
	parameter (libndz=25)

c---------------------------------------------------------------------
c
c	ipnt(id,isect) = lhead + lsect*(isect-1) + id
c
c---------------------------------------------------------------------

