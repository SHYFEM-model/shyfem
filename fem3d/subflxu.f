c
c $Id: subflxa.f,v 1.25 2009-05-21 09:24:00 georg Exp $
c
c utility routines for flux computations
c
c contents :
c
c subroutine flxscs(n,kflux,iflux,az,fluxes)	flux through sections
c subroutine flxsec(n,kflux,iflux,az,fluxes)	flux through section
c
c subroutine flxini				initializes flux routines
c subroutine flx_init(kfluxm,kflux,nsect,iflux)	sets up array iflux
c subroutine flxinf(m,kflux,iflux)		sets up one info structure
c function igtnsc(k1,k2)			gets number of internal section
c
c revision log :
c
c 09.05.2013    ggu     separated from subflxa.f
c 14.05.2013    ggu     deleted error check between 2d and 3d computation
c
c notes :
c
c These routines can also be used internally to compute the flux
c over various sections. The following calling sequence must be respected:
c
c call flx_init(kfluxm,kflux,nsect,iflux)		initializes iflux
c
c call flxscs(kfluxm,kflux,iflux,az,fluxes) computes fluxes 
c
c Initialization can be done anytime.
c
c******************************************************************
c******************************************************************
c******************************************************************
c******************************************************************
c******************************************************************


c******************************************************************
c******************************************************************
c******************************************************************
c******************************************************************
c******************************************************************


c******************************************************************
c******************************************************************
c******************************************************************
c******************************************************************
c******************************************************************

	subroutine fluxes_init(nlvddi,nsect,nlayers,nr,masst)

c initializes nr and masst

	implicit none

	integer nlvddi,nsect
	integer nlayers(nsect)
	integer nr
	real masst(0:nlvddi,3,nsect)

	integer i,l,lmax

        nr = 0
	masst = 0.

	end

c******************************************************************

	subroutine fluxes_accum(nlvddi,nsect,nlayers,nr,masst,fluxes)

c accumulates fluxes into masst

	implicit none

	integer nlvddi,nsect
	integer nlayers(nsect)
	integer nr
	real masst(0:nlvddi,3,nsect)
	real fluxes(0:nlvddi,3,nsect)

	integer i,l,lmax

        nr = nr + 1
	do i=1,nsect
	  lmax = nlayers(i)
	  do l=0,lmax
	    masst(l,1,i) = masst(l,1,i) + fluxes(l,1,i)
	    masst(l,2,i) = masst(l,2,i) + fluxes(l,2,i)
	    masst(l,3,i) = masst(l,3,i) + fluxes(l,3,i)
	  end do
	end do

	end

c******************************************************************

	subroutine fluxes_aver(nlvddi,nsect,nlayers,nr,masst,fluxes)

c averages masst and puts result into fluxes

	implicit none

	integer nlvddi,nsect
	integer nlayers(nsect)
	integer nr
	real masst(0:nlvddi,3,nsect)
	real fluxes(0:nlvddi,3,nsect)

	integer i,l,lmax
	real rr

        rr=1./nr
        do i=1,nsect
	  lmax = nlayers(i)
	  do l=0,lmax
            fluxes(l,1,i) = masst(l,1,i) * rr
            fluxes(l,2,i) = masst(l,2,i) * rr
            fluxes(l,3,i) = masst(l,3,i) * rr
	  end do
	end do

	end

c******************************************************************
c******************************************************************
c******************************************************************
c******************************************************************
c******************************************************************

	subroutine flxscs(kfluxm,kflux,iflux,az,fluxes,is,scalar)

c computes flux through all sections and returns them in fluxes
c
c flux are divided into total, positive and negative

	use levels, only : nlvdi,nlv

	implicit none

	integer kfluxm
	integer kflux(kfluxm)
	integer iflux(3,kfluxm)
	real az
	real fluxes(0:nlvdi,3,*)	!computed fluxes (return)
	integer is			!type of scalar (0=mass)
	real scalar(nlvdi,*)

	integer nnode,ifirst,ilast,ntotal
	integer ns
	logical nextline

	nnode = 0
	ns = 0

	do while( nextline(kflux,kfluxm,nnode,ifirst,ilast) )
	  ns = ns + 1
	  ntotal = ilast - ifirst + 1
	  !write(66,*) 'section ',ns,ntotal,is
	  call flxsec(ntotal,kflux(ifirst),iflux(1,ifirst),az
     +				,fluxes(0,1,ns),is,scalar)
	end do

	end

c******************************************************************

	subroutine flxsec(n,kflux,iflux,az,fluxes,is,scalar)

c computes flux through one section and returns it in fluxes

	use levels, only : nlvdi,nlv

	implicit none

	integer n
	integer kflux(n)
	integer iflux(3,n)
	real az
	real fluxes(0:nlvdi,3)		!computed fluxes (return)
	integer is			!type of scalar (0=mass)
	real scalar(nlvdi,*)

	integer i,k,l,lkmax
	integer istype,iafter,ibefor
	real port,ptot,port2d,sport,sptot
	real flux(nlvdi)

	fluxes = 0.

	do i=1,n
		k = kflux(i)
		istype = iflux(1,i)
		ibefor = iflux(2,i)
		iafter = iflux(3,i)

		call flx2d(k,ibefor,iafter,istype,az,port2d)

		call flx3d(k,ibefor,iafter,istype,az,lkmax,flux)

		ptot = 0.
		sptot = 0.
		do l=1,lkmax
		  port = flux(l)
		  sport = port
		  if( is .gt. 0 ) sport = port * scalar(l,k)	!not mass

		  fluxes(l,1) = fluxes(l,1) + sport
		  if( port .gt. 0. ) then
		    fluxes(l,2) = fluxes(l,2) + sport
		  else
		    fluxes(l,3) = fluxes(l,3) - sport
		  end if
		  ptot = ptot + port
		  sptot = sptot + sport
		end do
		!write(66,*) i,k,istype,ibefor,iafter,flux

		!the next will give an error on partially dry nodes
		!if( abs(port2d-ptot) .gt. 1. ) then
		!  write(6,*) '***** integrated fluxes: ',k
		!  write(6,*) '   ',port2d,ptot,abs(port2d-ptot)
		!end if

		fluxes(0,1) = fluxes(0,1) + sptot
		if( ptot .gt. 0. ) then
		  fluxes(0,2) = fluxes(0,2) + sptot
		else
		  fluxes(0,3) = fluxes(0,3) - sptot
		end if
	end do

	end
	  
c******************************************************************
c******************************************************************
c******************************************************************

	subroutine flx_init(kfluxm,kflux,nsect,iflux)

c sets up array iflux

	implicit none

        integer kfluxm		!total number of nodes in kflux
        integer kflux(kfluxm)	!nodes in sections
	integer nsect		!number of section (return)
	integer iflux(3,*)	!internal array for flux computation (return)

	integer ifirst,ilast,nnode,ntotal

	integer klineck
	logical nextline

c----------------------------------------------------------
c check nodes for compatibility
c----------------------------------------------------------

	nsect = klineck(kfluxm,kflux)

	if( nsect .lt. 0 ) then
	  write(6,*) 'errors in section $FLUX'
	  stop 'error stop : flxini'
	end if

c----------------------------------------------------------
c now set info structure for sections
c----------------------------------------------------------

	nnode = 0

	do while( nextline(kflux,kfluxm,nnode,ifirst,ilast) )
	  ntotal = ilast - ifirst + 1
c	  write(6,*) kfluxm,nnode,ifirst,ilast,ntotal
	  call flxinf(ntotal,kflux(ifirst),iflux(1,ifirst))
	end do

c----------------------------------------------------------
c end of routine
c----------------------------------------------------------

	end

c******************************************************************

	subroutine flxinf(n,kflux,iflux)

c sets up info structure iflux(3,1) for one section

	implicit none

	integer n
	integer kflux(n)
	integer iflux(3,n)

	integer i,k
	integer ktype
	integer kafter,kbefor

	integer igtnsc,flxtype

	do i=1,n
	  k = kflux(i)
	  ktype = flxtype(k)

	  iflux(1,i) = ktype

	  kbefor = 0
	  if( i .ne. 1 ) kbefor = kflux(i-1)
	  kafter = 0
	  if( i .ne. n ) kafter = kflux(i+1)

	  iflux(2,i) = igtnsc(k,kbefor)
	  iflux(3,i) = igtnsc(k,kafter)
	end do

	end

c******************************************************************

	function igtnsc(k1,k2)

c gets number of internal section in link index of k1

	use mod_geom

	implicit none

	integer igtnsc
	integer k1,k2
	integer elems(maxlnk)

	integer k,i,n,ie

	integer knext,kbhnd

	call get_elems_around(k1,maxlnk,n,elems)

c	deal with no section given

	igtnsc = 0
	if( k2 .le. 0 ) return

c	section in the middle

	do i=1,n
	  igtnsc = i
	  ie = elems(i)
	  k = knext(k1,ie)
	  if( k .eq. k2 ) return
	end do

c	in case we are on the boundary

	i = n
	igtnsc = igtnsc + 1
	ie = elems(i)
	k = kbhnd(k1,ie)
	if( k .eq. k2 ) return

c	no node found

	write(6,*) k1,k2
	write(6,*) k1,n
	write(6,*) (elems(i),i=1,n)
	call get_elems_around(k2,maxlnk,n,elems)
	write(6,*) k2,n
	write(6,*) (elems(i),i=1,n)
	stop 'error stop igtnsc: internal error (2)'
	end
	      
c**********************************************************************
c**********************************************************************
c**********************************************************************
c**********************************************************************
c**********************************************************************

	subroutine get_nlayers(kfluxm,kflux,nlayers,nlmax)

c computes maximum numer of layers for sections

	use levels

	implicit none

	integer kfluxm
	integer kflux(*)
	integer nlayers(*)	!total number of layers for sections (return)
	integer nlmax		!maximum layers for all sections (return)

	integer ns
	integer nnode,ifirst,ilast
	integer i,k,l,lmax

	logical nextline

	ns = 0
	nnode = 0
	nlmax = 0

	do while( nextline(kflux,kfluxm,nnode,ifirst,ilast) )
	  ns = ns + 1
	  lmax = 0
	  do i=ifirst,ilast
	    k = kflux(i)
	    l = ilhkv(k)
	    lmax = max(lmax,l)
	  end do
	  nlayers(ns) = lmax
	  nlmax = max(nlmax,lmax)
	end do

	end

c**********************************************************************
c**********************************************************************
c**********************************************************************
c**********************************************************************
c**********************************************************************

