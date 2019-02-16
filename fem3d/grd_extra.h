
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

c extra data structure for grd files

	integer nlidim
	integer nlndim

	parameter ( nlidim = 100 )	!maximum number of lines
	parameter ( nlndim = nkndim )	!maximum number of nodes in lines

	integer nli			!number of lines read
	integer iplv(nlidim)		!external line numbers

	!integer iarnv(nkndim)		!types for nodes
	integer iarlv(nkndim)		!types for lines

	real hllv(nlidim)		!depth of lines

	integer ipntlv(0:nlidim)	!pointer into line structure
	integer inodlv(nlndim)		!nodes of lines

	common /nli/nli
	common /iplv/iplv
	!common /iarnv/iarnv
	common /iarlv/iarlv
	common /hllv/hllv
	common /ipntlv/ipntlv
	common /inodlv/inodlv

	!save /nli/,/iplv/,/iarnv/,/iarlv/,/hllv/,/ipntlv/,/inodlv/
	save /nli/,/iplv/,/iarlv/,/hllv/,/ipntlv/,/inodlv/

