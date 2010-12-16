c
c $Id: ousutil.f,v 1.3 2009-09-14 08:20:58 georg Exp $
c
c utilities for OUS files
c
c revision log :
c
c 16.12.2010	ggu	copied from ousextr_gis.f
c
c******************************************************************

        subroutine get_records_from_stdin(ndim,irec)

c gets records to extract from stdin

        implicit none

        integer ndim
        integer irec(ndim)

        integer i,ir

        do i=1,ndim
          irec(i) = 0
        end do

        write(6,*) 'Please enter the record numbers to be extracted.'
        write(6,*) 'Enter every record on a single line.'
        write(6,*) 'Finish with 0 on the last line.'
        write(6,*) 'example:'
        write(6,*) '   5'
        write(6,*) '  10'
        write(6,*) '  15'
        write(6,*) '  0'
        write(6,*) ' '

        do while(.true.)
          write(6,*) 'Enter record to extract (0 to end): '
          ir = 0
          read(5,'(i10)') ir

          if( ir .le. 0 ) then
            return
          else if( ir .gt. ndim ) then
            write(6,*) 'Cannot extract records higher than ',ndim
            write(6,*) 'Please change ndim and recompile.'
          else
            irec(ir) = 1
          end if
        end do

        end

c******************************************************************

        subroutine transp2vel(nel,nkn,nlv,nlvdim,hev,zenv,nen3v
     +				,ilhv,hlv,utlnv,vtlnv
     +                          ,uprv,vprv,weight)

c transforms transports at elements to velocities at nodes

        implicit none

        integer nel
        integer nkn
	integer nlv
        integer nlvdim
        real hev(1)
        real zenv(3,1)
	integer nen3v(3,1)
	integer ilhv(1)
	real hlv(1)
        real utlnv(nlvdim,1)
        real vtlnv(nlvdim,1)
        real uprv(nlvdim,1)
        real vprv(nlvdim,1)
        real weight(nlvdim,1)

	logical bsigma
        integer ie,ii,k,l,lmax
        real zmed,hmed,u,v,w
	real hbot,htop,htot,hfact,hadd

	bsigma = hlv(nlv) .eq. -1.

	do k=1,nkn
	  do l=1,nlv
	    weight(l,k) = 0.
	    uprv(l,k) = 0.
	    vprv(l,k) = 0.
	  end do
	end do
	      
        do ie=1,nel
          zmed = 0.
          do ii=1,3
            zmed = zmed + zenv(ii,ie)
          end do
          zmed = zmed / 3.

	  htot = hev(ie)
	  htop = -zmed
          hfact = 1.
	  hadd = 0.
          if( bsigma ) htop = 0.
          if( bsigma ) then
	    hfact = -(htot+zmed)
	    hadd = -zmed
	  end if
	  lmax = ilhv(ie)
	  do l=1,lmax
	    if( l .eq. lmax ) hbot = hev(ie)
	    hbot = hadd + hfact * hlv(l)
	    hmed = hbot - htop
	    u = utlnv(l,ie) / hmed
	    v = vtlnv(l,ie) / hmed
	    do ii=1,3
	      k = nen3v(ii,ie)
	      uprv(l,k) = uprv(l,k) + u
	      vprv(l,k) = vprv(l,k) + v
	      weight(l,k) = weight(l,k) + 1.
	    end do
	    htop = hbot
	  end do
	end do

	do k=1,nkn
	  do l=1,nlv
	    w = weight(l,k)
	    if( w .gt. 0. ) then
	      uprv(l,k) = uprv(l,k) / w
	      vprv(l,k) = vprv(l,k) / w
	    end if
	  end do
	end do
	      
	end

c******************************************************************

	subroutine level_e2k(nkn,nel,nen3v,ilhv,ilhkv)

c computes ilhkv from ilhv

	implicit none

	integer nkn,nel
	integer nen3v(3,1)
	integer ilhv(1)
	integer ilhkv(1)

	integer k,ie,ii,lmax

	do k=1,nkn
	  ilhkv(k) = 1
	end do

	do ie=1,nel
	  lmax = ilhv(ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    if( ilhkv(k) .lt. lmax ) ilhkv(k) = lmax
	  end do
	end do

	end

c******************************************************************

