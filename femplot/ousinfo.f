c
c $Id: ousinfo.f,v 1.2 2009-04-07 09:33:35 georg Exp $
c
	programm ousinfo

c reads ous file
c
c this program demonstrates how to use the routines ousopen etc..

        implicit none

	include 'param.h'

	character*80 descrp
	common /descrp/ descrp
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real  grav,fcor,dcor,dirn,rowass,roluft
        common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft
        integer nlvdi,nlv
        common /level/ nlvdi,nlv

        real hm3v(neldim,1)
        common /hm3v/hm3v
        real xgv(nkndim), ygv(nkndim)
        common /xgv/xgv, /ygv/ygv
        integer nen3v(neldim,1), iarv(neldim)
        common /nen3v/nen3v, /iarv/iarv
        integer ipv(nkndim), ipev(neldim)
        common /ipv/ipv, /ipev/ipev

        real hlv(nlvdim), hev(neldim)
        common /hlv/hlv, /hev/hev
        integer ilhv(neldim)
        common /ilhv/ilhv
        real zenv(3,neldim)
        common /zenv/zenv
        real znv(nkndim)
        common /znv/znv
        real utlnv(nlvdim,neldim)
        common /utlnv/utlnv
        real vtlnv(nlvdim,neldim)
        common /vtlnv/vtlnv

	logical ousnext
        integer iapini
	integer it,nvers,l
        real rmin,rmax

c----------------------------------------------------

c	initialize dummy params

	nkn = nkndim
	nel = neldim
	nlv = nlvdim
	nlvdi = nlvdim

c	open simulation with basin

        if( iapini(3,nkndim,neldim,0) .eq. 0 ) then
                stop 'error stop : iapini'
        end if

c	initialize simulation

	call ousopen

        call ousinfo(nvers,nkn,nel,nlv)
        write(6,*) nvers,nkn,nel,nlv

c	loop on output

	do while( ousnext(it) )

	  write(6,*) it
          do l=1,nlv
            call minmax3d(nlvdim,nel,l,ilhv,utlnv,rmin,rmax)
            write(6,*) 'utmin/max: ',rmin,rmax
            call minmax3d(nlvdim,nel,l,ilhv,vtlnv,rmin,rmax)
            write(6,*) 'vtmin/max: ',rmin,rmax
          end do

	end do

c----------------------------------------------------

	end

c********************************************************************

        subroutine minmax3d(nlvdim,n,l,ilhv,array,rmin,rmax)

        implicit none

        integer nlvdim,n,l
        integer ilhv(1)
        real array(nlvdim,1)
        real rmin,rmax

        integer i

        rmin = 1.e+30
        rmax = -1.e+30

        do i=1,n
          if( l .le. ilhv(i) ) then
            rmin=min(rmin,array(l,i))
            rmax=max(rmax,array(l,i))
          end if
        end do

        end

c********************************************************************

