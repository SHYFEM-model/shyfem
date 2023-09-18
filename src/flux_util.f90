!
! $Id: subflxa.f,v 1.25 2009-05-21 09:24:00 georg Exp $
!
! utility routines for flux computations
!
! contents :
!
! subroutine flxscs(n,kflux,iflux,az,fluxes)	flux through sections
! subroutine flxsec(n,kflux,iflux,az,fluxes)	flux through section
!
! subroutine flxini				initializes flux routines
! subroutine flx_init(kfluxm,kflux,nsect,iflux)	sets up array iflux
! subroutine flxinf(m,kflux,iflux)		sets up one info structure
! function igtnsc(k1,k2)			gets number of internal section
!
! revision log :
!
! 09.05.2013    ggu     separated from subflxa.f
! 14.05.2013    ggu     deleted error check between 2d and 3d computation
!
! notes :
!
! These routines can also be used internally to compute the flux
! over various sections. The following calling sequence must be respected:
!
! call flx_init(kfluxm,kflux,nsect,iflux)		initializes iflux
!
! call flxscs(kfluxm,kflux,iflux,az,fluxes) computes fluxes 
!
! Initialization can be done anytime.
!
!******************************************************************
!******************************************************************
!******************************************************************
!------------------------------------------------------------------
        module flux_util
!------------------------------------------------------------------
        contains
!------------------------------------------------------------------

	subroutine fluxes_init(nlvddi,nsect,nlayers,nr,masst)

	implicit none

	integer nlvddi,nsect
	integer nlayers(nsect)
	integer nr
	double precision masst(0:nlvddi,3,nsect)

	integer i,l,lmax

        nr = 0
        do i=1,nsect
	  lmax = nlayers(i)
	  do l=0,lmax
            masst(l,1,i) = 0.
            masst(l,2,i) = 0.
            masst(l,3,i) = 0.
          end do
        end do

	end

!******************************************************************

	subroutine fluxes_accum(nlvddi,nsect,nlayers,nr,masst,fluxes)

	implicit none

	integer nlvddi,nsect
	integer nlayers(nsect)
	integer nr
	double precision masst(0:nlvddi,3,nsect)
	double precision fluxes(0:nlvddi,3,nsect)

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

!******************************************************************

	subroutine fluxes_aver(nlvddi,nsect,nlayers,nr,masst,fluxes)

	implicit none

	integer nlvddi,nsect
	integer nlayers(nsect)
	integer nr
	double precision masst(0:nlvddi,3,nsect)
	double precision fluxes(0:nlvddi,3,nsect)

	integer i,l,lmax
	double precision rr

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

!******************************************************************
!******************************************************************
!******************************************************************
!******************************************************************
!******************************************************************

	subroutine flxscs(kfluxm,kflux,iflux,az,fluxes,is,scalar)

! computes flux through all sections and returns them in flux
!
! flux are divided into total, positive and negative

	use levels, only : nlvdi,nlv
        use line_admin

	implicit none

	integer kfluxm
	integer kflux(kfluxm)
	integer iflux(3,kfluxm)
	double precision az
	double precision fluxes(0:nlvdi,3,*)	!computed fluxes (return)
	integer is			!type of scalar (0=mass)
	double precision scalar(nlvdi,*)

	integer nnode,ifirst,ilast,ntotal
	integer ns

	nnode = 0
	ns = 0

	do while( nextline(kflux,kfluxm,nnode,ifirst,ilast) )
	  ns = ns + 1
	  ntotal = ilast - ifirst + 1
	  !write(66,*) 'section ',ns,ntotal,is
	  call flxsec(ntotal,kflux(ifirst),iflux(1,ifirst),az,fluxes(0,1,ns),is,scalar)
	end do

	end

!******************************************************************

	subroutine flxsec(n,kflux,iflux,az,fluxes,is,scalar)

! computes flux through one section

	use levels, only : nlvdi,nlv
        use discharge_flux

	implicit none

	integer n
	integer kflux(n)
	integer iflux(3,n)
	double precision az
	double precision fluxes(0:nlvdi,3)		!computed fluxes (return)
	integer is			!type of scalar (0=mass)
	double precision scalar(nlvdi,*)

	integer i,k,l,lkmax
	integer istype,iafter,ibefor
	double precision port,ptot,port2d,sport,sptot
	double precision flux(nlvdi)

	do l=0,nlvdi
	  fluxes(l,1) = 0.
	  fluxes(l,2) = 0.
	  fluxes(l,3) = 0.
	end do

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
	  
!******************************************************************
!******************************************************************
!******************************************************************

	subroutine flx_init(kfluxm,kflux,nsect,iflux)

! sets up array iflux

        use line_admin

	implicit none

        integer kfluxm		!total number of nodes in kflux
        integer kflux(kfluxm)	!nodes in sections
	integer nsect		!number of section (return)
	integer iflux(3,*)	!internal array for flux computation (return)

	integer ifirst,ilast,nnode,ntotal

!----------------------------------------------------------
! check nodes for compatibility
!----------------------------------------------------------

	nsect = klineck(kfluxm,kflux)

	if( nsect .lt. 0 ) then
	  write(6,*) 'errors in section $FLUX'
	  stop 'error stop : flxini'
	end if

!----------------------------------------------------------
! now set info structure for sections
!----------------------------------------------------------

	nnode = 0

	do while( nextline(kflux,kfluxm,nnode,ifirst,ilast) )
	  ntotal = ilast - ifirst + 1
!	  write(6,*) kfluxm,nnode,ifirst,ilast,ntotal
	  call flxinf(ntotal,kflux(ifirst),iflux(1,ifirst))
	end do

!----------------------------------------------------------
! end of routine
!----------------------------------------------------------

	end

!******************************************************************

	subroutine flxinf(n,kflux,iflux)

! sets up info structure iflux(3,1) for one section

        use discharge_flux

	implicit none

	integer n
	integer kflux(n)
	integer iflux(3,n)

	integer i,k
	integer ktype
	integer kafter,kbefor

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

!******************************************************************

	function igtnsc(k1,k2)

! gets number of internal section in link index of k1

	use geom
        use lnku

	implicit none

	integer igtnsc
	integer k1,k2
	integer elems(maxlnk)

	integer k,i,n,ie

	call get_elems_around(k1,maxlnk,n,elems)

!	deal with no section given

	igtnsc = 0
	if( k2 .le. 0 ) return

!	section in the middle

	do i=1,n
	  igtnsc = i
	  ie = elems(i)
	  k = knext(k1,ie)
	  if( k .eq. k2 ) return
	end do

!	in case we are on the boundary

	i = n
	igtnsc = igtnsc + 1
	ie = elems(i)
	k = kbhnd(k1,ie)
	if( k .eq. k2 ) return

!	no node found

	write(6,*) k1,k2
	write(6,*) k1,n
	write(6,*) (elems(i),i=1,n)
	call get_elems_around(k2,maxlnk,n,elems)
	write(6,*) k2,n
	write(6,*) (elems(i),i=1,n)
	stop 'error stop igtnsc: internal error (2)'
	end
	      
!**********************************************************************
!**********************************************************************
!**********************************************************************
!**********************************************************************
!**********************************************************************

	subroutine get_nlayers(kfluxm,kflux,nlayers,nlmax)

	use levels
        use line_admin

	implicit none

	integer kfluxm
	integer kflux(*)
	integer nlayers(*)
	integer nlmax

	integer ns
	integer nnode,ifirst,ilast
	integer i,k,l,lmax

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

!**********************************************************************
!**********************************************************************
!**********************************************************************

!------------------------------------------------------------------
        end module flux_util
!------------------------------------------------------------------
