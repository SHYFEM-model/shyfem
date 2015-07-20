c
c $Id: nos_elab.f,v 1.2 2009-04-07 10:43:57 georg Exp $
c
c revision log :
c
c 24.01.2011    ggu     written from scratch
c
c****************************************************************

	program basres

c computes resolution in basin

	use basin

	implicit none

	include 'param.h'

c--------------------------------------------------
c--------------------------------------------------

	character*50 file
	character*80 title
	real cvres(nkndim)
	real cvaux(nkndim)

	integer ilhkv(nkndim)
	real hlv(nlvdim)
	real hev(neldim)

	integer nb,ierr
	real rmin,rmax

	integer iapini,ifileo

c---------------------------------------------------------------
c open basin
c---------------------------------------------------------------

	if(iapini(1,nkndim,neldim,0).eq.0) then
		stop 'error stop : iapini'
	end if

c---------------------------------------------------------------
c initializing units and nodes to be extracted
c---------------------------------------------------------------

	call makehev(hev)

        call mkname(' ','basres','.nos',file)
        write(6,*) 'writing file ',file(1:50)
        nb = ifileo(0,file,'unform','new')
        if( nb .le. 0 ) goto 98

	title = 'basin resolution'
        call wfnos(nb,3,nkn,nel,1,1,title,ierr)
        if( ierr .ne. 0 ) goto 97
        call wsnos(nb,ilhkv,hlv,hev,ierr)
        if( ierr .ne. 0 ) goto 97


c---------------------------------------------------------------
c compute and write resolution
c---------------------------------------------------------------

	call compute_resolution(cvres,cvaux)

	call mima(cvres,nkn,rmin,rmax)
	write(6,*) 'min/max resolution: ',rmin,rmax

        call wrnos(nb,0,334,1,ilhkv,cvres,ierr)    !aver
        if( ierr .ne. 0 ) goto 99

	write(6,*)
	write(6,*) 'data written to file basres.nos'
	write(6,*)

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

	stop
   97   continue
        write(6,*) 'error writing header'
        stop 'error stop nosaver'
   98   continue
        write(6,*) 'error opening file'
        stop 'error stop nosaver'
   99   continue
        write(6,*) 'error writing data'
        stop 'error stop nosaver'
	end

c***************************************************************

	subroutine compute_resolution(cvres,cvaux)

c computes horizontal resolution

	use basin

	implicit none

	real cvres(1)
	real cvaux(1)

	include 'param.h'

	integer k,ie,ii,iii,k1,k2
	real dd

	real dist

	do k=1,nkn
	  cvres(k) = 0.
	  cvaux(k) = 0.
	end do
	
	do ie=1,nel
	  do ii=1,3
	    iii = 1 + mod(ii,3)
	    k1 = nen3v(ii,ie)
	    k2 = nen3v(iii,ie)
	    dd = dist(k1,k2)
	    cvres(k1) = cvres(k1) + dd
	    cvaux(k1) = cvaux(k1) + 1.
	    cvres(k2) = cvres(k2) + dd
	    cvaux(k2) = cvaux(k2) + 1.
	  end do
	end do

	do k=1,nkn
	  if( cvaux(k) .gt. 0. ) cvres(k) = cvres(k) / cvaux(k)
	end do

	end

c***************************************************************


	function dist(k1,k2)

c computes distance between nodes

	use basin

	implicit none

	real dist
	integer k1,k2

	include 'param.h'

	real x1,x2,y1,y2,dx,dy

	x1 = xgv(k1)
	y1 = ygv(k1)
	x2 = xgv(k2)
	y2 = ygv(k2)

	dx = x1 - x2
	dy = y1 - y2

	dist = sqrt( dx*dx + dy*dy )

	end

c***************************************************************
