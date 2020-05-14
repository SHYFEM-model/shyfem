
!--------------------------------------------------------------------------
!
!    Copyright (C) 2018-2019  Georg Umgiesser
!    Copyright (C) 2018  Petras Zemlys
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

! custtomization for various aspects of the plot
!
! revision log :
!
! 07.06.2018	ggu	new file for customization from Petras
! 07.06.2018	pzy	font size definitions
! 06.07.2018	ggu	changed VERS_7_5_48
! 18.12.2018	ggu	changed VERS_7_5_52
! 21.05.2019	ggu	changed VERS_7_5_62
!
!**************************************************************

!==================================================================
	module plot_fonts
!==================================================================
!
! Module for font size definition of different map components
!
! Written by Petras May-June 2018
!
!------------------------------------------------------------------

	implicit none
           
	! font size for colorbar
	integer, parameter :: fs_cbar_text = 24 ! for legend text
	integer, parameter :: fs_cbar_num  = 15 ! for numbers
		
	! font size for black-white frame		
	integer, parameter :: fs_bw_frame = 15 ! for coordinate values
		
!==================================================================
	end module plot_fonts
!==================================================================

