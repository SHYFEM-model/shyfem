!
! $Id: bastreat.f,v 1.14 2010-02-26 17:35:06 georg Exp $
!
! revision log :
!
! 06.04.1999    ggu     completely restructured
! 04.06.1999    ggu     new statistics are computed
! 08.09.2003    ggu     mode 5 -> write depth values from elements
! 23.09.2004    ggu     interpolq() changed for bathy interpolation
! 02.10.2004    ggu     interpole() for exponential interpolation
! 01.11.2004    ggu     whole program simplyfied
! 06.12.2008    ggu     smoothing introduced
! 06.04.2009    ggu     read param.h
! 29.05.2009    ggu     does only depth limiting and smoothing
! 20.11.2009    ggu     possibility to smooth only on specific areas
! 30.03.2011    ggu     new routines to delete elements
! 13.06.2013    ggu     copy_depth() renamed to transfer_depth()
!
!****************************************************************

        subroutine bas_partition

! performs modifications on basin
!
! takes care of lat/lon coordinates

	use mod_geom
	use mod_depth
	use evgeom
	use basin
	use grd
	use basutil

	implicit none

	integer k,i,nl,il,n,ib,in,node,ilext,np,ic
	integer, allocatable :: nc(:)
	real x,y,perc
	real xx(nkn)
	real yy(nkn)

	logical inconvex,inpoly,filex

!-----------------------------------------------------------------
! open and read file containing lines
!-----------------------------------------------------------------

	if( .not. filex(lfile) ) then
	  write(6,*) 'Cannot open file ',trim(lfile)
	  stop 'error stop bas_partition: no such file'
	end if

	call grd_read(lfile)

	iarnv = 0
	nl = nl_grd

!-----------------------------------------------------------------
! loop over lines
!-----------------------------------------------------------------

        do i=1,nl
          il = i
          ilext = ipplv(il)
          n = ipntlv(il) - ipntlv(il-1)
          ib = ipntlv(il-1)
	  !write(6,*) 'line: ',i,ilext,n
	  do in=1,n
	    node = inodlv(ib+in)
	    x = xv(node)
	    y = yv(node)
	    xx(in) = x
	    yy(in) = y
	    !write(6,*) in,node,x,y
	  end do
	  if( xx(1) == xx(n) .and. yy(1) == yy(n) ) n = n - 1
	  np = 0
	  do k=1,nkn
	    x = xgv(k)
	    y = ygv(k)
	    !write(6,*) 'looking for... ',k,x,y
	    if( inpoly(n,xx,yy,x,y) ) then
	      iarnv(k) = il
	      np = np + 1
	      !write(6,*) 'inside... ',il,n,k,x,y
	    end if
	  end do
	  perc = (100.*np)/nkn
	  !write(6,*) 'inside points found... ',i,il,n,np,perc
        end do

	nl = nl + 1
	where( iarnv == 0 ) iarnv = nl

!-----------------------------------------------------------------
! write information to terminal
!-----------------------------------------------------------------

	allocate(nc(0:nl))
	nc = 0
	do k=1,nkn
	  ic = iarnv(k)
	  if( ic < 1 .or. ic > nl ) then
	    write(6,*) 'ic,nl: ',ic,nl
	    stop 'error stop bas_partition: internal error (1)'
	  end if
	  nc(ic) = nc(ic) + 1
	end do
	write(6,*) 'Information on domains: '
	write(6,*) '   domain     nodes   percent'
	do ic=1,nl
	  write(6,'(2i10,f10.2)') ic,nc(ic),(100.*nc(ic))/nkn
	end do

!-----------------------------------------------------------------
! write file
!-----------------------------------------------------------------

        call basin_to_grd
        call grd_write('bas_partition.grd')
        write(6,*) 'The partition has been written to bas_partition.grd'

!-----------------------------------------------------------------
! end of routine
!-----------------------------------------------------------------

	return
	end

!*******************************************************************

