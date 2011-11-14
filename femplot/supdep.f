c
c $Id: supdep.f,v 1.5 2009-11-18 17:16:00 georg Exp $
c
c routines for averaging depth
c
c contents :
c
c subroutine mkhkv(hkv,auxv,nkn,nel)	makes hkv (nodewise depth)
c subroutine mkhev(hev,nel)		makes hev (elementwise depth)
c subroutine mkht(hetv,href)		makes hetv (elem depth of actual layer)
c subroutine mkht3(nlvdim,het3v,href)	makes het3v (3D depth structure)
c function hlthick(l,lmax,hl)		layer thickness
c
c revision log :
c
c 26.05.2000    ggu     routines written from scratch
c 17.09.2008    ggu     routine mkht changed for layer = -1
c 13.10.2009    ggu     new routine mkht3
c 13.10.2009    ggu     new routine mkht3
c 17.12.2010    ggu     substituted hv with hkv, new routine hlthick()
c 30.03.2011    ggu     new routine mkareafvl()
c 14.11.2011    ggu     use get_layer_thickness() for layer structure
c
c******************************************************************

	subroutine mkhkv(hkv,auxv,nkn,nel)

c makes hkv (nodewise depth)

	implicit none

c arguments
	real hkv(1)
	real auxv(1)
	integer nkn,nel
c common
	integer nen3v(3,1)
	real hm3v(3,1)
	common /nen3v/nen3v, /hm3v/hm3v
c local
	integer ie,ii,k,kn

        do k=1,nkn
          hkv(k) = 0.
        end do

        do ie=1,nel
          do ii=1,3
            kn=nen3v(ii,ie)
            hkv(kn)=hkv(kn)+hm3v(ii,ie)
	    auxv(kn)=auxv(kn)+1.
          end do
        end do

        do k=1,nkn
          hkv(k) = hkv(k) / auxv(k)
        end do

	return
	end

c******************************************************************

	subroutine mkhev(hev,nel)

c makes hev (elementwise depth)

	implicit none

c arguments
	real hev(1)
	integer nel
c common
	real hm3v(3,1)
	common /hm3v/hm3v
c local
	integer ie,ii
	real h

        do ie=1,nel
	  h=0.
          do ii=1,3
            h=h+hm3v(ii,ie)
          end do
	  hev(ie) = h / 3.
        end do

	return
	end

c******************************************************************

	subroutine mkht(hetv,href)

c makes hetv (elementwise depth of actual layer)
c
c uses level to decide what to do

	implicit none

c arguments
	real hetv(1)
	real href
c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer nlvdi,nlv
        common /level/ nlvdi,nlv
        real hlv(1), hev(1)
        common /hlv/hlv, /hev/hev
        integer ilhv(1)
        common /ilhv/ilhv
	real zenv(3,1)
	common /zenv/zenv
        real vev(1)
        common /vev/vev
        real hl(1)
        common /hl/hl
c local
	logical bdebug,bzeta
	integer ie,ii
	integer level,lmax
	real z,h
c functions
	integer getlev
	real hlthick

c-------------------------------------------------------------------
c initialization
c-------------------------------------------------------------------

        bdebug = .true.
        bdebug = .false.
	bzeta = .true.

	level = getlev()

c-------------------------------------------------------------------
c handle different kind of levels
c-------------------------------------------------------------------

	do ie=1,nel
	  lmax = ilhv(ie)
	  call get_layer_thickness(ie,lmax,bzeta,hl)
	  hetv(ie) = hlthick(level,lmax,hl)
	end do

c-------------------------------------------------------------------
c debug output
c-------------------------------------------------------------------

	if( bdebug ) then
	  ie = 1644
	  write(6,*) 'debugging mkht...'
	  write(6,*) ilhv(ie),hetv(ie),hev(ie)
	  do ii=1,ilhv(ie)
	    write(6,*) hlv(ii)
	  end do
	end if

c-------------------------------------------------------------------
c end of routine
c-------------------------------------------------------------------

	end

c******************************************************************

	subroutine mkht3(nlvdim,het3v,href)

c makes het3v (3D depth structure)

	implicit none

c arguments
	integer nlvdim
	real het3v(nlvdim,1)
	real href
c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer nlvdi,nlv
        common /level/ nlvdi,nlv
        integer ilhv(1)
        common /ilhv/ilhv
c local
	logical bdebug,bzeta
	integer ie,ii
	integer l,lmax
	real z,h

c-------------------------------------------------------------------
c initialization
c-------------------------------------------------------------------

        bdebug = .true.
        bdebug = .false.
	bzeta = .true.

c-------------------------------------------------------------------
c compute layer thickness
c-------------------------------------------------------------------

	do ie=1,nel
	  lmax = ilhv(ie)
	  call get_layer_thickness(ie,lmax,bzeta,het3v(1,ie))
	end do

c-------------------------------------------------------------------
c end of routine
c-------------------------------------------------------------------

	end

c******************************************************************

	function hlthick(l,lmax,hl)

c computes thickness of layer l
c
c works also for sigma layers

	implicit none

	real hlthick		! layer thickness (return)
	integer l		! layer to compute thickness
	integer lmax		! maximum layers
	real hl(1)		! level thickness

	integer ll
	real hm

	if( l .eq. 0 ) then			! total layer
	  hm = 0.
	  do ll=1,lmax
	    hm = hm + hl(ll)
	  end do
	  hlthick = hm
	else if( l .ge. 1 .and. l .le. lmax ) then
	  hlthick = hl(l)
	else if( l .eq. -1 ) then		! bottom layer
	  hlthick = hl(lmax)
	else
	  hlthick = 0.
	end if

	end

c******************************************************************

	subroutine mkareafvl

c makes area of finite volume

	implicit none

c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nen3v(3,1)
	common /nen3v/nen3v
	real arfvlv(1)
	common /arfvlv/arfvlv
	include 'ev.h'
c local
	integer ie,ii,k
	real afvl

	do k=1,nkn
	  arfvlv(k) = 0.
	end do

	do ie=1,nel
	  afvl = 4. * ev(10,ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    arfvlv(k) = arfvlv(k) + afvl
	  end do
	end do

	end

c******************************************************************

