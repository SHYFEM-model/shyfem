!
! shyelab_average.f: utility for averaging
!
! revision log :
!
! 07.10.2017    ggu     started
!
!***********************************************************

	subroutine average_vertical_node(lmax,hlv,z,htot,values,aver)

! averages vertically a profile of values

	implicit none

	integer lmax
	real hlv(lmax)
	real z,htot
	real values(lmax)
	real aver		!return

	integer l
	integer nlvaux,nsigma
	real hsigma
	real h
	real hd(lmax)
	double precision vaccum,haccum

	aver = 0.
	if( lmax == 1 ) aver = values(1)
	if( lmax <= 1 ) return

        call get_sigma_info(nlvaux,nsigma,hsigma)
        call get_layer_thickness(lmax,nsigma,hsigma,z,htot,hlv,hd)

	vaccum = 0.
	haccum = 0.
	do l=1,lmax
	  h = hd(l)
	  vaccum = vaccum + values(l) * h
	  haccum = haccum + h
	end do

	aver = vaccum / haccum

	end

!***********************************************************

