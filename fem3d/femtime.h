
!--------------------------------------------------------------------------
!
!    Copyright (C) 2014-2015,2017-2019  Georg Umgiesser
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
! 26.11.2014	ggu	changed VERS_7_0_7
! 23.12.2014	ggu	changed VERS_7_0_11
! 23.09.2015	ggu	changed VERS_7_2_4
! 12.10.2015	ggu	changed VERS_7_3_3
! 04.11.2017	ggu	changed VERS_7_5_34
! 03.04.2018	ggu	changed VERS_7_5_43
! 16.02.2019	ggu	changed VERS_7_5_60
! 18.05.2022	ggu	new arrays for accumulation of cpu timing info

!--------------------------------------------------------------------------

	integer itanf,itend,idt,nits,niter,it
	common /femtim/ itanf,itend,idt,nits,niter,it

        integer itunit,idtorig
        common /femtimu/ itunit,idtorig

	double precision t_act,dt_act,dt_orig,atime0,dtanf,dtend
        common /femtimd/ t_act,dt_act,dt_orig,atime0,dtanf,dtend

	logical bsync
        common /femtiml/ bsync

	character*20 aline_act
        common /femtimc/ aline_act

	integer, parameter :: ncpu = 10
	double precision cputime(ncpu)
	double precision acutime(ncpu)
        common /femtimp/ cputime,acutime

        save /femtim/,/femtimu/,/femtimd/,/femtiml/,/femtimc/,/femtimp/

