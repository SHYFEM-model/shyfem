c
c $Id: newbsig.f,v 1.3 2009-03-24 17:49:19 georg Exp $
c
c routines for handling initialization from sigma levels
c
c contents :
c
c revision log :
c
c 07.11.2008    ggu     program written from scratch
c 06.12.2008    ggu     bbsig set from STR (nbsig)
c 24.03.2009    ggu     use nbsig to check if initialized
c 16.12.2010    ggu     adjusted for sigma levels, renamed sigma.h to bsig.h
c 25.10.2011    ggu     some routines renamed
c 03.11.2011    ggu     use routine set_hybrid_depth() for hybrid levels
c 02.12.2011    ggu     check in set_bsig_depths()
c 02.12.2011    ggu     use of intp_aver() for lmax = 1
c 27.01.2012    deb&ggu adapted to hybrid levels
c 05.10.2012    ggu     intp_vert() and intp_aver() transfered to subvintp.f
c
c*****************************************************************

	subroutine handle_bsig_init

c initializes u/v/z from sigma level data

	implicit none

	include 'param.h'
	include 'bsig.h'

	character*(60) sigma_l,sigma_d,sigma_u,sigma_v,sigma_z,sigma_b
	logical bbsig
        real getpar

	bbsig = nint(getpar('nbsig')) .gt. 0	!if true read in sigma data

	nbsig = 0			!number of sigma levels read (init)

	if( .not. bbsig ) return	!no initialization with sigma data

	sigma_l = 'RosSigmas.txt'
	sigma_d = 'PtsBathyAll.txt'
	sigma_b = 'PtsBathyBoundary.txt'

	sigma_u = 'Uvelocities_init.txt'
	sigma_v = 'Vvelocities_init.txt'
	sigma_z = 'elevation_init.txt'

	call set_bsig_levels(sigma_l)
	call set_bsig_depths(sigma_d)
	call init_vel_from_sigma(sigma_u,sigma_v,sigma_z)

	end

c*****************************************************************

	subroutine set_bsig_levels(file)

c read in sigma levels

	implicit none

	character*(*) file

	include 'param.h'
	include 'bsig.h'

	integer i

	nbsig = 0

	open(1,file=file)
	read(1,*) nbsig
        nbsig = nbsig - 1	!number of layers, not interfaces

	if( nbsig .gt. nsidim ) goto 98

	do i=0,nbsig
	  read(1,*) siglev(i)
	end do

	close(1)

	if( siglev(0) .ne. 0. ) goto 99
	if( siglev(nbsig) .ne. 1. ) goto 99

	write(6,*) nbsig,' sigma layers read from file ',file
	write(6,*) 'sigma levels read: ',(siglev(i),i=0,nbsig)

	return
   98	continue
	write(6,*) 'too many sigma levels'
	write(6,*) 'nbsig = ',nbsig,'    nsidim = ',nsidim
	stop 'error stop set_bsig_levels: nsidim'
   99	continue
	write(6,*) 'problems in sigma levels...'
	write(6,*) nbsig
	write(6,*) (siglev(i),i=0,nbsig)
	stop 'error stop set_bsig_levels: sigma'
	end

c*****************************************************************

	subroutine set_bsig_depths(file)

c read in depth file for sigma layers

	implicit none

	include 'param.h'
	include 'bsig.h'

	character*(*) file

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        logical berror
	integer k,l
	integer kint,kext
	integer nknaux
	real x,y
	real hfem,hfd

	do k=1,nkn
	  sigdep(k) = -1.
	end do

	open(1,file=file)
	read(1,*) nknaux

	if( nknaux .ne. nkn ) goto 98

	do k=1,nkn
	  !read(1,*) kint,kext,l,x,y,hfem,hfd
	  read(1,*) kext,x,y,hfem,hfd
	  sigdep(k) = hfd
	end do

	close(1)

	write(6,*) 'sigma layer depth read from file ',file

        berror = .false.
	do k=1,nkn
          if( sigdep(k) .eq. -1. ) then
            berror = .true.
            write(6,*) k
          end if
	end do
        if( berror ) then
          write(6,*) 'errors in reading sigma depth'
	  stop 'error stop set_bsig_depths: depth'
        end if

	return
   98	continue
	write(6,*) 'not compatible number of total nodes'
	write(6,*) 'nknaux = ',nknaux,'    nkn = ',nkn
	stop 'error stop set_bsig_depths: nkn'
	end

c*****************************************************************

	subroutine init_vel_from_sigma(fileu,filev,filez)

c initializes u/v/z from sigma level data

	implicit none

	include 'param.h'
	include 'bsig.h'

	character*(*) fileu,filev,filez

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        integer nen3v(3,1)
        common /nen3v/nen3v
        integer ilhv(1)
        common /ilhv/ilhv
        real znv(1)
        common /znv/znv
        real zenv(3,1)
        common /zenv/zenv
        real ulnv(nlvdim,1),vlnv(nlvdim,1)
        common /ulnv/ulnv, /vlnv/vlnv

	real sigu(nsidim,nkndim)
	real sigv(nsidim,nkndim)
	real sigz(nkndim)

	real femu(nlvdim)
	real femv(nlvdim)

	integer k,l,lmax,ie,ii,ifact
	integer kint,kext
	integer nknaux
	integer itime
	real x,y,h

	integer ipext

	if( nbsig .gt. nsidim ) goto 99

c----------------------------------------------------------
c read in water levels
c----------------------------------------------------------

	write(6,*) 'start reading water levels from file ',filez
	open(1,file=filez)
	read(1,*) itime,nknaux
	if( nknaux .ne. nkn ) goto 98
	do k=1,nkn
	  !read(1,*) kint,kext,lmax,x,y,sigz(k)
	  read(1,*) kext,x,y,h,lmax,sigz(k)
	  if( ipext(k) .ne. kext ) goto 97
	end do
	close(1)

c----------------------------------------------------------
c read in velocities on nodes and sigma levels
c----------------------------------------------------------

	write(6,*) 'start reading u velocities from file ',fileu
	open(1,file=fileu)
	read(1,*) itime,nknaux
	if( nknaux .ne. nkn ) goto 98
	do k=1,nkn
	  !read(1,*) kint,kext,lmax,x,y,(sigu(l,k),l=1,lmax)
	  read(1,*) kext,x,y,h,lmax,(sigu(l,k),l=1,lmax)
	  if( ipext(k) .ne. kext ) goto 97
	end do
	close(1)

	write(6,*) 'start reading v velocities from file ',filev
	open(1,file=filev)
	read(1,*) itime,nknaux
	if( nknaux .ne. nkn ) goto 98
	do k=1,nkn
	  !read(1,*) kint,kext,lmax,x,y,(sigv(l,k),l=1,lmax)
	  read(1,*) kext,x,y,h,lmax,(sigv(l,k),l=1,lmax)
	  if( ipext(k) .ne. kext ) goto 97
	end do
	close(1)

	write(6,*) 'finished reading initial fields'

c----------------------------------------------------------
c initialize water levels and adjust depth
c----------------------------------------------------------

	do k=1,nkn
	  znv(k) = sigz(k)
	end do

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    zenv(ii,ie) = znv(k)
	  end do
	end do

c----------------------------------------------------------
c loop over elements and do interpolation
c----------------------------------------------------------

	ifact = 0		!do not compute fact in interpolation routine

	do ie=1,nel

	  lmax = ilhv(ie)

	  do l=1,lmax
	    femu(l) = 0.
	    femv(l) = 0.
	  end do

	  do ii=1,3
	    k = nen3v(ii,ie)
	    call add_sigma_to_element(ie,k,ifact,nsidim
     +                          ,sigu,sigv,femu,femv)
	  end do

	  do l=1,lmax
	    ulnv(l,ie) = femu(l) / 3.
	    vlnv(l,ie) = femv(l) / 3.
	  end do

	end do

c----------------------------------------------------------
c set all other variables
c----------------------------------------------------------

	call make_new_depth	!makes layer thickness
	call vtot		!makes utlnv, vtlnv
	call init_uv		!here we need zenv, utlnv, vtlnv, hdenv

c----------------------------------------------------------
c end of routine
c----------------------------------------------------------

	return
   97	continue
	write(6,*) 'incompatible nodes'
	write(6,*) k,ipext(k),kext
	stop 'error stop init_vel_from_sigma: kext'
   98	continue
	write(6,*) 'not compatible number of total nodes'
	write(6,*) 'nknaux = ',nknaux,'    nkn = ',nkn
	stop 'error stop init_vel_from_sigma: nkn'
   99	continue
	write(6,*) 'too many sigma levels'
	write(6,*) 'nbsig = ',nbsig,'    nsidim = ',nsidim
	stop 'error stop init_vel_from_sigma: nsidim'
	end

c*****************************************************************

	subroutine add_sigma_to_element(ie,k,ifact,ndim,uval,vval,u,v)

	implicit none

	include 'param.h'
	include 'bsig.h'

	integer ie,k,ifact,ndim
	real uval(ndim,nkndim)
	real vval(ndim,nkndim)
	real u(1),v(1)

        integer ilhv(1)
        common /ilhv/ilhv
        real znv(1),hev(1)
        common /znv/znv, /hev/hev
        real hlv(1)
	common /hlv/hlv

	real hlfem(0:nlvdim+1)
	real femuval(nlvdim)
	real femvval(nlvdim)

	real hsig(0:nsidim+1)
	real siguval(nsidim+1)
	real sigvval(nsidim+1)

	logical bsigma,bcons
	integer nsigma
	integer l,lmax
	real hfem,hfd,zfem
	real fact
	real hsigma

	if( nbsig .le. 0 ) goto 98

	bcons = .false.

	call get_sigma(nsigma,hsigma)
	bsigma = nsigma .gt. 0

	hfd = sigdep(k)
	zfem = znv(k)		!this might not be set at the boundary
	hfem = hev(ie)
	lmax = ilhv(ie)

	hlfem(0) = -zfem
	call set_hybrid_depth(lmax,zfem,hfem,hlv,nsigma,hsigma,hlfem(1))

	fact = 1.
	if( ifact .ne. 0 ) fact = (hfd+zfem) / (hfem+zfem)

	call make_sigma_layer_bottom(nbsig,siglev,hfem,zfem,hsig)

	do l=1,nbsig
	  siguval(l) = uval(l,k)
	  sigvval(l) = vval(l,k)
	end do

        if( lmax .eq. 1 ) then
          call intp_aver(nbsig,hsig,siguval,femuval(1))
          call intp_aver(nbsig,hsig,sigvval,femvval(1))
        else
	  call intp_vert(bcons,nbsig,hsig,siguval,lmax,hlfem,femuval)
	  call intp_vert(bcons,nbsig,hsig,sigvval,lmax,hlfem,femvval)
        end if

	do l=1,lmax
	  u(l) = u(l) + fact * femuval(l)
	  v(l) = v(l) + fact * femvval(l)
	end do

	return
   98	continue
	stop 'error stop add_sigma_to_element: sigma layers not init'
	end

c*****************************************************************

	subroutine make_sigma_layer_thickness(nbsig,siglev,h,z,hlayer)

c compute layer thickness for depth h and water level z

	implicit none

	integer nbsig
	real siglev(0:nbsig)
	real h,z
	real hlayer(nbsig)

	integer l

	do l=1,nbsig
	  hlayer(l) = (h+z) * (siglev(l) - siglev(l-1))
	end do

	end

c*****************************************************************

	subroutine make_sigma_layer_bottom(nbsig,siglev,h,z,hbottom)

c compute bottom of layers for depth h and water level z

	implicit none

	integer nbsig
	real siglev(0:nbsig)
	real h,z
	real hbottom(0:nbsig)

	integer l

	do l=0,nbsig
	  hbottom(l) = - z - (h+z) * siglev(l)
	end do

	end

c*****************************************************************

