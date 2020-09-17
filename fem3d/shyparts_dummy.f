
!--------------------------------------------------------------------------
!
!    Copyright (C) 2017-2019  Georg Umgiesser
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

c revision log :
c
c 19.05.2020	ccf	started from scratch
c
c****************************************************************

        program shyparts

c shyparts dummy routine

	implicit none

        write(6,*)' For automatic partitioning of a grid install' 
        write(6,*)' one of the following libraries and set the'
        write(6,*)' parameters PARTS and PARTSDIR in the'
        write(6,*)' Rules.make configuration file' 
        write(6,*)'   - METIS'
        write(6,*)' Then recompile: "make fem"'

	end

c*******************************************************************

        subroutine node_test
        end

c*******************************************************************


