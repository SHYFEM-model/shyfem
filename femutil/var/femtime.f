
!--------------------------------------------------------------------------
!
!    Copyright (C) 2020  Georg Umgiesser
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
! 06.03.2020	ggu	changed VERS_7_5_69

	program test

	call set_include
	call set_module
	call write_include
	call write_module

	end

	module femtime
	include 'femtime.h'
	end module

	module femtime1
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
        save /femtim/,/femtimu/,/femtimd/,/femtiml/,/femtimc/
	end module

	subroutine set_include
	include 'femtime.h'
	it = 99
	end

	subroutine set_module
	use femtime
	dt = 88
	end

	subroutine write_include
	include 'femtime.h'
	write(6,*) 'include: ',it,dt
	end

	subroutine write_module
	use femtime1
	write(6,*) 'module: ',it,dt
	end

