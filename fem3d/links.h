
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

	integer maxlnk
	parameter( maxlnk = 100 )

        integer ilinkv(nkndim+1)
        common /ilinkv/ ilinkv
        integer lenkv(nlkdim)
        common /lenkv/ lenkv
        integer lenkiiv(nlkdim)
        common /lenkiiv/ lenkiiv
        integer linkv(nlkdim)
        common /linkv/ linkv

        !integer ilinkv(1)
        !common /ilinkv/ ilinkv
        !integer lenkv(1)
        !common /lenkv/ lenkv
        !integer linkv(1)
        !common /linkv/ linkv

	integer lnk_nodes(maxlnk)
	common /lnk_nodes/ lnk_nodes

	integer lnk_elems(maxlnk)
	common /lnk_elems/ lnk_elems

	save /ilinkv/,/lenkv/,/linkv/,/lenkiiv/
	save /lnk_nodes/,/lnk_elems/

