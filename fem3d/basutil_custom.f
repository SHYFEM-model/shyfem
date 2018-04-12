c
c $Id: bastreat.f,v 1.14 2010-02-26 17:35:06 georg Exp $
c
c revision log :
c
c 06.04.1999    ggu     completely restructured
c 04.06.1999    ggu     new statistics are computed
c 08.09.2003    ggu     mode 5 -> write depth values from elements
c 23.09.2004    ggu     interpolq() changed for bathy interpolation
c 02.10.2004    ggu     interpole() for exponential interpolation
c 01.11.2004    ggu     whole program simplyfied
c 06.12.2008    ggu     smoothing introduced
c 06.04.2009    ggu     read param.h
c 29.05.2009    ggu     does only depth limiting and smoothing
c 20.11.2009    ggu     possibility to smooth only on specific areas
c 30.03.2011    ggu     new routines to delete elements
c 13.06.2013    ggu     copy_depth() renamed to transfer_depth()
c
c****************************************************************

        subroutine bas_custom

c performs modifications on basin
c
c takes care of lat/lon coordinates

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
	character*80 file

	logical inconvex,inpoly,filex

c-----------------------------------------------------------------
c what to do
c-----------------------------------------------------------------

	file = 'lines.grd'
	if( .not. filex(file) ) return

	call grd_read(file)

	iarnv = 0

	nl = nl_grd

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

	allocate(nc(0:nl))
	nc = 0
	do k=1,nkn
	  ic = iarnv(k)
	  if( ic > nl ) stop 'error...'
	  nc(ic) = nc(ic) + 1
	end do
	do ic=0,nl
	  write(6,*) ic,nc(ic),(100.*nc(ic))/nkn
	end do

c-----------------------------------------------------------------
c write
c-----------------------------------------------------------------

        call basin_to_grd
        call grd_write('basnodes.grd')
        write(6,*) 'The basin has been written to basnodes.grd'

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	return
	end

c*******************************************************************

	subroutine color_areai_0(iestart,icol,icon)

	use basin
	use mod_geom

	implicit none

	integer iestart,icol
	integer icon(nel)

	integer ip,ien,ii,ie
	integer list(nel)

	ip = 1
	list(ip) = iestart

	do while( ip .gt. 0 ) 
	  ie = list(ip)
	  icon(ie) = icol
	  !write(6,*) icol,ip,ie
	  ip = ip -1
	  do ii=1,3
	    ien = ieltv(ii,ie)
	    if( ien .gt. 0 ) then
	      if( icon(ien) .eq. 0 ) then
	        ip = ip + 1
	        list(ip) = ien
	      else if( icon(ien) .ne. icol ) then
	        goto 99
	      end if
	    end if
	  end do
	end do

	return
   99	continue
	write(6,*) ip,ie,ien,icol,icon(ien)
	stop 'error stop color_area: internal error (1)'
	end

c*******************************************************************

