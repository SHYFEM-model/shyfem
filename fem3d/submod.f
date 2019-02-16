
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

c***************************************************************
c
c look in subcst for initialization
c
c cstinit	M_INIT		start of hp/ht
c cstcheck	M_CHECK		after read of STR file
c cstsetup	M_SETUP		after setup of basic arrays (ev,...)
c
c look in subnsh for read, write and others
c
c prilog	M_PRINT		after everything has been setup (before loop)
c pritst	M_TEST		only for debug
c
c dobefor	M_BEFOR		beginning of each time loop
c doafter	M_AFTER		end of each time loop
c
c nlsh2d	M_READ		read in STR file
c
c***************************************************************
c
c modules still to transform:
c
c wrouta
c resid
c rmsvel
c adjust_chezy (bottom friction)
c prwnds
c prarea (checy)
c prclos
c proxy, prlgr
c prbnds
c pripar
c biocos
c locous, locspc
c
c***************************************************************

	subroutine modules(mode)

c handles module framework

	implicit none

	integer mode		!mode of call

c------------------------------------------------------------
c modules
c------------------------------------------------------------

	call mod_ext(mode)		!extra nodes
c	call mod_ets(mode)		!extra time series nodes
c	call mod_box(mode)		!boxes
c	call mod_flx(mode)		!fluxes through sections
c	call mod_vol(mode)		!volumes in areas

c------------------------------------------------------------
c end of routine
c------------------------------------------------------------

	end

c***************************************************************

