
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

	integer nbc_dim
	parameter ( nbc_dim = 100 )

        character*80 boundn(nbc_dim)
        character*80 conzn(nbc_dim)
        character*80 saltn(nbc_dim)
        character*80 tempn(nbc_dim)
        character*80 bio2dn(nbc_dim)
        character*80 sed2dn(nbc_dim)
        character*80 mud2dn(nbc_dim)
        character*80 lam2dn(nbc_dim)
        character*80 dmf2dn(nbc_dim)
        character*80 tox3dn(nbc_dim)
        character*80 bfm1bc(nbc_dim)
        character*80 bfm2bc(nbc_dim)
        character*80 bfm3bc(nbc_dim)
        character*80 vel3dn(nbc_dim)
        character*80 bfmbcn(nbc_dim)

        common /boundn/ boundn
        common /conzn/ conzn
        common /saltn/ saltn
        common /tempn/ tempn
        common /bio2dn/ bio2dn
        common /sed2dn/ sed2dn
        common /mud2dn/ mud2dn
        common /lam2dn/ lam2dn  !!!!!!!!!!!!!!!!! BUG
        common /dmf2dn/ dmf2dn
        common /tox3dn/ tox3dn
        common /bfm1bc/ bfm1bc
        common /bfm2bc/ bfm2bc
        common /bfm3bc/ bfm3bc
        common /vel3dn/ vel3dn
        common /bfmbcn/ bfmbcn

	save /boundn/,/conzn/,/saltn/,/tempn/
	save /bio2dn/,/sed2dn/,/mud2dn/,/lam2dn/,/dmf2dn/,/tox3dn/
	save /bfm1bc/,/bfm2bc/,/bfm3bc/
	save /vel3dn/,/bfmbcn/

