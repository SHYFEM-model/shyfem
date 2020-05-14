
!--------------------------------------------------------------------------
!
!    Copyright (C) 2010,2018  Georg Umgiesser
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
! 23.03.2010	ggu	changed v6.1.1
! 22.02.2018	ggu	changed VERS_7_5_42
! 25.10.2018	ggu	changed VERS_7_5_51
! 18.12.2018	ggu	changed VERS_7_5_52

        integer legdim			!maximum number of legend entries
        parameter(legdim=400)

        integer nleg,nlegdi,iplotleg
        common /nleg/ nleg,nlegdi,iplotleg

        real xleg(2,legdim), yleg(2,legdim)
        common /xyleg/ xleg,yleg

        character*80 legleg(legdim)
        common /legleg/ legleg

        character*4 whatleg(legdim)
        common /whatleg/ whatleg

        integer legsiz(legdim)
        common /legsiz/ legsiz

        real aleg(legdim), cleg(legdim)
        common /aleg/ aleg
        common /cleg/ cleg

	save /nleg/,/xyleg/,/legsiz/
	save /legleg/
	save /whatleg/
	save /aleg/,/cleg/

