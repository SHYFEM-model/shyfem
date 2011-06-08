c
c $Id: ousutil.f,v 1.3 2009-09-14 08:20:58 georg Exp $
c
c utilities for OUS files
c
c revision log :
c
c 16.12.2010	ggu	copied from ousextr_gis.f
c 03.06.2011	ggu	some routines transfered to genutil.f
c 08.06.2011	ggu	new routine transp2nodes()
c
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

        subroutine transp2nodes(nel,nkn,nlv,nlvdim,hev,zenv,nen3v
     +				,ilhv,hlv,utlnv,vtlnv
     +                          ,utprv,vtprv,weight)

c transforms transports at elements to transports at nodes

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
        real utprv(nlvdim,1)
        real vtprv(nlvdim,1)
        real weight(nlvdim,1)

	logical bsigma
        integer ie,ii,k,l,lmax
        real u,v,w

	bsigma = hlv(nlv) .eq. -1.

	do k=1,nkn
	  do l=1,nlv
	    weight(l,k) = 0.
	    utprv(l,k) = 0.
	    vtprv(l,k) = 0.
	  end do
	end do
	      
        do ie=1,nel
	  lmax = ilhv(ie)
	  do l=1,lmax
	    u = utlnv(l,ie)
	    v = vtlnv(l,ie)
	    do ii=1,3
	      k = nen3v(ii,ie)
	      utprv(l,k) = utprv(l,k) + u
	      vtprv(l,k) = vtprv(l,k) + v
	      weight(l,k) = weight(l,k) + 1.
	    end do
	  end do
	end do

	do k=1,nkn
	  do l=1,nlv
	    w = weight(l,k)
	    if( w .gt. 0. ) then
	      utprv(l,k) = utprv(l,k) / w
	      vtprv(l,k) = vtprv(l,k) / w
	    end if
	  end do
	end do
	      
	end

c***************************************************************
c***************************************************************
c***************************************************************

        subroutine open_ous_file(name,status,nunit)

        implicit none

        character*(*) name,status
        integer nunit

        integer nb
        character*80 file

        integer ifileo

        call mkname(' ',name,'.ous',file)
        nb = ifileo(0,file,'unform',status)
        if( nb .le. 0 ) stop 'error stop open_ous_file: opening file'

        nunit = nb

        end

c***************************************************************

        subroutine qopen_ous_file(text,status,nunit)

c asks for name and opens ous file

        implicit none

        character*(*) text,status
        integer nunit

        character*80 name

        write(6,*) text
        read(5,'(a)') name
        write(6,*) name

        call open_ous_file(name,status,nunit)

        end

c***************************************************************
c***************************************************************
c***************************************************************

