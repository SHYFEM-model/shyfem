c
c $Id: ousinfo.f,v 1.2 2009-04-07 09:33:35 georg Exp $
c
	program ousinfo

c reads ous file
c
c this program demonstrates how to use the routines ousopen etc..

        implicit none

	include 'param.h'
	include 'basin.h'

	include 'nlevel.h'

	include 'levels.h'
	include 'hydro.h'

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

