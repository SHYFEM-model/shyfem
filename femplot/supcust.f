!
! custtomization for various aspects of the plot
!
! revision log :
!
! 07.06.2018    ggu     new file for customization from Petras
! 07.06.2018    pzy     font size definitions
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

