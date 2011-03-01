c
c $Id: newini.f,v 1.45 2009-11-18 16:50:37 georg Exp $
c
c routines for initialization
c
c contents :
c
c subroutine inicfil(name,var,nvar)
c				initializes nodal value variable from file
c subroutine inic2fil(name,var,nvar)
c				initializes nodal value variable from file (2D)
c
c revision log :
c
c 20.08.1998	ggu	initialization routines copied from subn11/new11
c 20.08.1998	ggu	new routine inicfil to init nodal values from file
c 06.11.1998	ggu	computation of hkv/hev taken out of sp211 -> cstset
c 22.10.1999	ggu	igrv, ngrv eliminated (subst by ilinkv, lenkv)
c 07.03.2000    ggu     arrays for vertical turbulent diffusion coefficient
c 07.03.2000    ggu     useless parts commented
c 20.06.2000    ggu     useless parts deleted
c 21.06.2000    ggu     dzreg introduced
c 08.08.2000    ggu     hlvmin is now percentage of last layer thickness
c 25.03.2003    ggu     voldist adjusted for 3D, new routine surinl
c 07.08.2003    ggu     deleted sp212
c 09.08.2003    ggu     cleaned up (sp211 is completely different)
c 14.08.2003    ggu     error in check_ilevels (integer instead of real)
c 14.08.2003    ggu     error in check_ilevels (must test .lt.)
c 14.08.2003    ggu     new routine set_depth -> transferred from cstcheck
c 14.08.2003    ggu     renamed sp211 into init_vertical and init_others
c 14.08.2003    ggu     new routines init_z, init_const_z, init_file_z
c 14.08.2003    ggu     new routines init_uvt
c 02.09.2003    ggu     routine inicfil restructured for lmax=0
c 03.09.2003    ggu     some more checks for check_ilevels
c 04.03.2004    ggu     change in inicfil() for more variables
c 05.03.2004    ggu     LMAX - changed bound from 0 to 1
c 02.09.2004    ggu     rdrst transfered to subrst.f
c 15.03.2005    ggu     set_austausch() and check_austausch() eliminated
c 05.04.2005    ggu     minimum levels (for baroclinic terms) -> set_min_levels
c 28.11.2005    ggu     in set_depth use makehkv to compute hkv
c 16.02.2006    ggu     describe file format in inicfil()
c 15.11.2006    ggu     new parameters to construct streched vert. coordinates
c 27.08.2007    ccf     variable coriolis from spherical coordinates (isphe)
c 17.03.2008    ggu     better error message if missing levels
c 07.04.2008    ggu     deleted surinl, volinl, voldist
c 12.11.2008    ggu     set_last_layer() rewritten, new adjust_k_depth()
c 06.12.2008    ggu&ccf small bug fix in set_last_layer()
c 21.01.2009    ggu	cleaned up, new read_scalar()
c 27.01.2009    ggu	mzreg and hlvmax deleted (not used)
c 24.03.2009    ggu	bug fix: in set_last_layer() do not adjust for 2D
c 21.04.2009    ggu	new routine inic2fil()
c 28.09.2010    ggu	bug fix in init_coriolis() for isphe=1
c 08.10.2010    ggu	bug fix in init_coriolis() -> ym not set for isphe=1
c 16.12.2010    ggu	big restructering for sigma levels (look for bsigma)
c 21.02.2011    ggu	error check for dzreg, nsigma and levels
c
c**************************************************************

	subroutine init_vertical

c set up time independent vertical vectors

	implicit none

	write(6,*) 'setting up vertical structure'

c------------------------------------------------------------------
c sanity check
c------------------------------------------------------------------

	call check_nlv

c------------------------------------------------------------------
c levels read in from $levels section
c------------------------------------------------------------------

	call adjust_levels	!sets hlv, hldv

c------------------------------------------------------------------
c adjust depth values
c------------------------------------------------------------------

	call set_depth	!computes hev, hkv and reassigns to hm3v -> href

c------------------------------------------------------------------
c set up depth vectors
c------------------------------------------------------------------

	call set_ilhv		!sets ilhv (elemental)
	call set_last_layer	!adjusts hev, hm3v, ilhv - sets up hlhv
	call set_ilhkv		!sets ilhkv (nodal)
	call set_min_levels	!sets ilmkv and ilmv

	call adjust_k_depth	!adjusts nodal depth values (hkv)

c------------------------------------------------------------------
c check data structure
c------------------------------------------------------------------

	call check_vertical

c------------------------------------------------------------------
c end of routine
c------------------------------------------------------------------

	end

c**************************************************************

	subroutine init_others

c set up various arrays (coriolis, eddy)

	implicit none

c------------------------------------------------------------------
c set others
c------------------------------------------------------------------

	call init_coriolis
	call set_eddy

c------------------------------------------------------------------
c end of routine
c------------------------------------------------------------------

	end

c*****************************************************************

	subroutine set_eddy

c sets vertical eddy coefficient

	implicit none

	include 'param.h'

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nlvdi,nlv
	common /level/ nlvdi,nlv

	real visv(0:nlvdim,1)
	real difv(0:nlvdim,1)
	common /visv/visv, /difv/difv

	integer k,l
	real vistur,diftur

	real getpar

c------------------------------------------------------------------
c get parameters
c------------------------------------------------------------------

	vistur=getpar('vistur')
	diftur=getpar('diftur')

c------------------------------------------------------------------
c set eddy coefficient
c------------------------------------------------------------------

	do k=1,nkn
	  do l=0,nlv
	    visv(l,k)=vistur
	    difv(l,k)=diftur
	  end do
	end do

c------------------------------------------------------------------
c end of routine
c------------------------------------------------------------------

	end

c*****************************************************************

	subroutine check_eddy

c checks vertical eddy coefficient

	implicit none

	include 'param.h'

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nlvdi,nlv
	common /level/ nlvdi,nlv

	real visv(0:nlvdim,1)
	real difv(0:nlvdim,1)
	common /visv/visv, /difv/difv

	integer k,l
	real v,d

c------------------------------------------------------------------
c check eddy coefficient
c------------------------------------------------------------------

	do k=1,nkn
	  do l=0,nlv
	    v = visv(l,k)
	    d = difv(l,k)
	    if( v .lt. 0. .or. v .gt. 1.e+5 ) goto 99
	    if( d .lt. 0. .or. d .gt. 1.e+5 ) goto 99
	  end do
	end do

c------------------------------------------------------------------
c end of routine
c------------------------------------------------------------------

	return
   99	continue
	write(6,*) k,l,v,d
	stop 'error stop check_eddy: error in values'
	end

c*****************************************************************

	subroutine set_depth

c sets up depth arrays and adjusts depth values

	implicit none

        real hkv(1), hev(1)
        common /hkv/hkv, /hev/hev
	real v1v(1)
	common /v1v/v1v
	real hldv(1)
	common /hldv/hldv

	logical bsigma
	real hmin,hmax,href
	real getpar

c       call bocche     !FIXME

	bsigma = hldv(1) .lt. 0.

        hmin=getpar('hmin')
        hmax=getpar('hmax')
        href=getpar('href')

        call depadj(hmin,hmax,href)	!adjusts h=h-href and hmax<h<hmin

	if( bsigma ) then
	  call set_sigma_hkhe
	else
          call huniqu(hev,hkv)		!computes hev and reassignes to hm3v
          call makehkv(hkv,v1v)		!computes hkv as average
	end if

	end

c*****************************************************************

	subroutine adjust_k_depth

c adjusts nodal depth values

	implicit none

        real hkv(1)
        common /hkv/hkv
	real v1v(1)
	common /v1v/v1v

        call makehkv(hkv,v1v)		!computes hkv as average

	end

c*****************************************************************

	subroutine check_vertical

c checks arrays containing vertical structure

	implicit none

	call check_nlv
	call check_levels
	call check_ilevels
	call check_hlhv

	end

c*****************************************************************

	subroutine check_nlv

c checks nlv and associated parameters

	implicit none

	include 'param.h'

	integer nlvdi,nlv
	common /level/ nlvdi,nlv

	write(6,*) 'check_nlv : ',nlvdim,nlvdi,nlv

	if(nlvdim.ne.nlvdi) stop 'error stop check_nlv: level dimension'
	if(nlv.gt.nlvdim) stop 'error stop check_nlv: level dimension'

	end

c*****************************************************************

	subroutine adjust_levels

c adjusts levels read in from $levels section
c
c creates hlv if not set
c from hlv creates hldv

	implicit none

c common
	integer nlvdi,nlv
	common /level/ nlvdi,nlv
	real hldv(1),hlv(1)
	common /hldv/hldv, /hlv/hlv

c local
	logical bstop,bsigma
	integer l,nsigma
	real dzreg,hbot,hl,fact
	real getpar

c--------------------------------------------------------------
c create hlv values
c--------------------------------------------------------------

	dzreg = getpar('dzreg')
	nsigma = nint(getpar('nsigma'))

	if( nlv .le. 0 ) then		!no level section given
	  if( nsigma .gt. 0 ) then	!sigma layers
	    nlv = nsigma
	    hl = -1. / nlv
	    do l=1,nlv
	      hbot = hbot + dzreg
	      hlv(l) = l * hl
	    end do
	  else if( dzreg .gt. 0. ) then	!regular grid
	    nlv = nlvdi
	    hbot = 0.
	    do l=1,nlv
	      hbot = hbot + dzreg
	      hlv(l) = hbot
	    end do
	  else				!just one layer
	    nlv = 1
	    hlv(1) = 10000.
	  end if
	else
	  if( dzreg .gt. 0. ) then
	    write(6,*) 'nlv,dzreg: ',nlv,dzreg
	    write(6,*) 'You cannot give both levels and dzreg.'
	    stop 'error stop adjust_levels: dzreg'
	  end if
	  if( nsigma .gt. 0 ) then
	    write(6,*) 'nlv,nsigma: ',nlv,nsigma
	    write(6,*) 'You cannot give both levels and nsigma.'
	    stop 'error stop adjust_levels: nsigma'
	  end if
	end if

	write(6,*) 'adjust_levels: ',nlv,(hlv(l),l=1,nlv)

c--------------------------------------------------------------
c check hlv values
c--------------------------------------------------------------

	bstop = .false.
	bsigma = hlv(nlv) .eq. -1.	!sigma levels

	fact = 1.
	if( bsigma ) fact = -1.

	do l=2,nlv
	  if( fact*hlv(l) .le. fact*hlv(l-1) ) then
	    write(6,*) 'Error in level values for level : ',l
	    write(6,*) '   hlv(l-1),hlv(l) :',hlv(l-1),hlv(l)
	    bstop = .true.
	  end if
	end do
	if( bstop ) stop 'error stop adjust_levels: level error'

c--------------------------------------------------------------
c create hldv values
c--------------------------------------------------------------

	hldv(1)=hlv(1)
	do l=2,nlv
	  hldv(l)=hlv(l)-hlv(l-1)
	end do

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	end

c*****************************************************************

	subroutine check_levels

c checks arrays hlv and hldv

	implicit none

c common
	integer nlvdi,nlv
	common /level/ nlvdi,nlv
	real hldv(1),hlv(1)
	common /hldv/hldv, /hlv/hlv

c local
	logical bstop,bsigma
	integer l
	real h,hd,fact

	bsigma = hldv(1) .lt. 0.
	fact = 1.
	if( bsigma ) fact = -1.		!sigma levels

c--------------------------------------------------------------
c check hlv values
c--------------------------------------------------------------

	bstop = .false.
	do l=2,nlv
	  if( fact*hlv(l) .le. fact*hlv(l-1) ) then
	    write(6,*) 'Error in level values for level : ',l
	    write(6,*) '   hlv(l-1),hlv(l) :',hlv(l-1),hlv(l)
	    bstop = .true.
	  end if
	end do
	if( bstop ) stop 'error stop check_levels: level error'

c--------------------------------------------------------------
c check hldv values
c--------------------------------------------------------------

	h = 0.
	do l=1,nlv
	  hd=hlv(l)-h
	  if( hd .ne. hldv(l) ) then
	    write(6,*) 'Error in dlevel values for level : ',l
	    write(6,*) '   hd,hldv(l) :',hd,hldv(l)
	    bstop = .true.
	  end if
	  h = hlv(l)
	end do
	if( bstop ) stop 'error stop check_levels: dlevel error'

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	end

c*****************************************************************

	subroutine set_ilhv

c sets array ilhv

	implicit none

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nlvdi,nlv
	common /level/ nlvdi,nlv

	integer ilhv(1)
	common /ilhv/ilhv
	real hlv(1)
	common /hlv/hlv
	real hev(1)
	common /hev/hev

c local
	logical bsigma
	integer ie,l,lmax
	real h,hmax

	lmax=0
	hmax = 0.

	bsigma = hlv(nlv) .eq. -1.

	do ie=1,nel
	  hmax = max(hmax,hev(ie))
	end do

	do ie=1,nel

	  h=hev(ie)

	  if( bsigma ) then
	    l = nlv
	  else
	    do l=1,nlv
	      if(hlv(l).ge.h) goto 1
	    end do
    1	    continue
	    if( l .gt. nlv ) goto 99
	  end if

	  ilhv(ie)=l
	  lmax = max(lmax,l)

	end do

	nlv = lmax

	return
   99	continue
	write(6,*) ie,l,nlv,h,hlv(nlv)
	write(6,*) 'maximum basin depth: ',hmax
	write(6,*) 'maximum layer depth: ',hlv(nlv)
	stop 'error stop set_ilhv: not enough layers'
	end

c*****************************************************************

	subroutine set_ilhkv

c set ilhkv array

	implicit none

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	integer ilhv(1)
	common /ilhv/ilhv
	integer ilhkv(1)
	common /ilhkv/ilhkv
	integer nen3v(3,1)
	common /nen3v/nen3v

	integer ie,ii,k,l

	do k=1,nkn
	  ilhkv(k)=0
	end do

	do ie=1,nel
	  l=ilhv(ie)
	  do ii=1,3
	    k=nen3v(ii,ie)
	    if(l.gt.ilhkv(k)) ilhkv(k)=l
	  end do
	end do

	end

c*****************************************************************

	subroutine set_min_levels

c set minimum number of levels for node and element

	implicit none

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nlvdi,nlv
	common /level/ nlvdi,nlv

	integer ilhv(1)
	common /ilhv/ilhv
	integer ilhkv(1)
	common /ilhkv/ilhkv

	integer ilmv(1)
	common /ilmv/ilmv
	integer ilmkv(1)
	common /ilmkv/ilmkv

	integer nen3v(3,1)
	common /nen3v/nen3v

	integer ie,ii,k,l
	integer lmin,lmax

	do k=1,nkn
	  ilmkv(k) = 99999
	end do

	do ie=1,nel
	  l=ilhv(ie)
	  do ii=1,3
	    k=nen3v(ii,ie)
	    if(l.lt.ilmkv(k)) ilmkv(k)=l
	  end do
	end do

	do k=1,nkn
	  lmin = ilmkv(k)
	  lmax = ilhkv(k)
	  if( lmin .gt. lmax .or. lmin .lt. 0 ) then
	    stop 'error stop set_min_levels: lmin'
	  end if
	  !write(6,*) 'set_min_levels (nodes): ',k,lmin,lmax
	end do

	do ie=1,nel
	  lmin = 99999
	  lmax = ilhv(ie)
	  do ii=1,3
	    k=nen3v(ii,ie)
	    if(lmin.gt.ilmkv(k)) lmin = ilmkv(k)
	  end do
	  ilmv(ie) = lmin
	  !write(6,*) 'set_min_levels (elems): ',ie,lmin,lmax
	end do

	end

c*****************************************************************

	subroutine check_ilevels

c checks arrays ilhv and ilhkv

	implicit none

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	integer nen3v(3,1)
	common /nen3v/nen3v
	integer ilhv(1)
	common /ilhv/ilhv
	integer ilhkv(1)
	common /ilhkv/ilhkv
	real hev(1)
	common /hev/hev
	real hlv(1)
	common /hlv/hlv
	real hldv(1)
	common /hldv/hldv
	real hlhv(1)
	common /hlhv/hlhv

	logical bsigma
	integer nsigma
	integer ie,ii,k,lmax,lk
	real h,hlev

	bsigma = hldv(1) .lt. 0.
	nsigma = ilhv(1)

	do ie=1,nel
	  lmax = ilhv(ie)
	  hlev = hlv(lmax)
	  h = hev(ie)
	  if( lmax .le. 0 ) goto 99
	  !if( h .gt. hlev ) goto 99		!might have changed
	  if( .not. bsigma .and. lmax .gt. 1 ) then
	    if( hlv(lmax-1) + hlhv(ie) - h .ne. 0. ) goto 99
	  end if
	end do

	do ie=1,nel
	  lmax=ilhv(ie)
	  do ii=1,3
	    k=nen3v(ii,ie)
	    lk = ilhkv(k)
	    if( lk .le. 0 ) goto 98
	    if( lk .lt. lmax ) goto 98
	  end do
	  if( bsigma .and. lmax .ne. nsigma ) goto 96
	end do

	do k=1,nkn
	  lmax=ilhkv(k)
	  if( lmax .le. 0 ) goto 97
	  if( bsigma .and. lmax .ne. nsigma ) goto 96
	end do

	return
   96	continue
	write(6,*) ie,k,lmax,nsigma
	stop 'error stop check_ilevels: error in vertical structure (4)'
   97	continue
	write(6,*) k,lmax
	stop 'error stop check_ilevels: error in vertical structure (3)'
   98	continue
	write(6,*) ie,lmax,k,lk
	stop 'error stop check_ilevels: error in vertical structure (2)'
   99	continue
	write(6,*) ie,lmax,hlev,h,hlhv(ie)
	stop 'error stop check_ilevels: error in vertical structure (1)'
	end

c*****************************************************************

	subroutine set_last_layer

c sets last layer thickness
c
c adjusts hev, hm3v, ilhv - sets up hlhv

	implicit none

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nlvdi,nlv
	common /level/ nlvdi,nlv

	integer ilhv(1)
	common /ilhv/ilhv
	real hlv(1)
	common /hlv/hlv
	real hldv(1)
	common /hldv/hldv
	real hlhv(1)
	common /hlhv/hlhv
	real hev(1)
	common /hev/hev
	real hm3v(3,1)
	common /hm3v/hm3v

c local
	logical bwrite
	logical badjust,b2d,bsigma
	integer ie,l,ii
	integer ihtot,lmax
	integer ilytyp
        integer ic1,ic2,ic3
	real h,hold,hnew,hlast
	real hlvmin
	real hmin,hmax

	integer ieext
	real getpar

	bwrite = .false.

c------------------------------------------------------------
c see if 2D application
c------------------------------------------------------------

	b2d = nlv .eq. 1			!2D application

c------------------------------------------------------------
c check if sigma levels
c------------------------------------------------------------

	bsigma = hldv(1) .lt. 0.
	if( bsigma ) then		!sigma layers
	  do ie=1,nel
	    hlhv(ie) = 0.
	  end do
	  write(6,*) 'set_last_layer (nlv used) : ',nlv
	  write(6,*) 'sigma layers detected'
	  return
	end if

c------------------------------------------------------------
c set hlvmin
c------------------------------------------------------------

c	ilytyp: 
c	  0=no adjustment  
c	  1=adjust to full layers (change depth)
c         2=adjust only if h<hlvmin (change depth)
c         3=add to last layer (keep depth but change layer)

	hlvmin = getpar('hlvmin')		!min percentage of last layer
	ilytyp = nint(getpar('ilytyp'))

	if( hlvmin .gt. 1. .or. hlvmin .lt. 0. ) goto 98

c------------------------------------------------------------
c adjust last layer thickness and store it in hlhv
c------------------------------------------------------------

	ic1 = 0
	ic2 = 0
	ic3 = 0

	ihtot = 0
	lmax = 0

	do ie=1,nel

	  l = ilhv(ie)
	  hold = hev(ie)			!original depth of element
	  h = hold
	  if( l .gt. 1 ) h = h - hlv(l-1)	!actual layer thickness
	  hlast = hldv(l)			!regular layer thickness
	  hmin = hlvmin * hlast

	  if( l .gt. 1 .and. h .le. 0. ) goto 99

	  badjust = .false.

	  if( l .gt. 1 ) then			!only for more than 1 layer
	    if( ilytyp .eq. 0 ) then
c		no adjustment
	    else if( h .le. hmin ) then		!take away layer
		badjust = .true.
		l = l-1
	    else if( ilytyp .eq. 1 ) then
		badjust = .true.
	    end if
	  end if

	  if( ilytyp .eq. 1 ) badjust = .true.	!adjust also first layer
	  if( b2d ) badjust = .false.		!but never for 2D application

	  if( badjust ) then
	    if( ilytyp .le. 2 ) then		!adjust depth
		h = hldv(l)
		hnew = hlv(l)
		hev(ie) = hnew
		do ii=1,3
		  hm3v(ii,ie) = hnew
		end do
	    else if( ilytyp .eq. 3 ) then
		h = hldv(l) + h
	    else
		stop 'error stop set_last_layer: internal error (7)'
	    end if

	    ilhv(ie) = l

	    if( bwrite ) then
	      if( ihtot .eq. 0 ) then
		  write(6,*) 'set_last_layer: Adjustment of depth values'
	      end if
	      write(6,'(2i8,3f12.2)') ie,ieext(ie),l,hev(ie),hold,h
	    end if
	    ihtot = ihtot + 1
	  end if

	  hlhv(ie) = h			!finally set last layer thickness
	  lmax = max(lmax,l)

	  if( hlv(l) .ne. hev(ie) ) ic1 = ic1 + 1
	  if( hldv(l) .ne. hlhv(ie) ) ic2 = ic2 + 1
	  if( l .eq. 1 ) ic3 = ic3 + 1
	end do

c------------------------------------------------------------
c adjust nlv and write final statistics
c------------------------------------------------------------

	nlv=lmax

	write(6,*) 'set_last_layer (nlv used) : ',nlv
	write(6,*) 'Total number of elements with adjusted depth: '
     +			,ihtot
	write(6,*) 'Incomplete depth:     ',ic1
	write(6,*) 'Differing last layer: ',ic2
	write(6,*) 'One layer elements:   ',ic3

c------------------------------------------------------------
c end of routine
c------------------------------------------------------------

	return
   98	continue
	write(6,*) 'hlvmin is given in fraction of last layer'
	write(6,*) 'therefore:  0 <= hlvmin <= 1'
	write(6,*) 'hlvmin = ',hlvmin
	stop 'error stop set_last_layer: hlvmin'
   99	continue
	stop 'error stop set_last_layer: layer depth 0 (internal error)'
	end

c*****************************************************************

	subroutine check_hlhv

c checks depth structure 

	implicit none

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	integer ilhv(1)
	common /ilhv/ilhv
	real hldv(1)
	common /hldv/hldv
	real hev(1)
	common /hev/hev
	real hlhv(1)
	common /hlhv/hlhv

c local
	logical bstop
	logical bsigma
	integer ie,l
	integer lmax
	real h

	bstop = .false.
	bsigma = hldv(1) .lt. 0.

	do ie=1,nel

	  lmax = ilhv(ie)
	  h=hlhv(ie)

	  do l=1,lmax-1
	    h=h+hldv(l)
	  end do

	  if( bsigma ) then
	    h = h + hldv(lmax)
	    h = - h * hev(ie)
	  end if

	  if(abs(h-hev(ie)).gt.1.e-4) then
c	  if(abs(h-hev(ie)).gt.0.) then
		write(6,*) 'ie,htot,hsum,lmax : ',ie,hev(ie),h,lmax
		bstop = .true.
	  end if

	end do

	if( bstop ) then
	  stop 'error stop check_hlhv: error in depth structure'
	end if

	end

c*****************************************************************

	subroutine init_coriolis

c sets coriolis parameter

	implicit none

	real omega2	!double frequency of earth rotation
	parameter ( omega2 = 2.0 * 0.729E-4 )
	real rearth	!radius of earth
	parameter ( rearth = 6371000. )

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real eps1,eps2,pi,flag,high
	real grav,fcor,dcor,dirn,rowass,roluft
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /mkonst/ eps1,eps2,pi,flag,high
	common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft

	integer nen3v(3,1)
	real ygv(1)
	real fcorv(1)
	common /nen3v/nen3v
	common /ygv/ygv
	common /fcorv/fcorv 

	integer k,ie,ii
	integer icor
	integer isphe
	real yc,ym,y,ymin,ymax
	real aux1,aux2,rad
	real getpar

c if coordinates are cartesian (isphe=0) then icor determines 
c how Coriolis is used:
c
c 0:	no Coriolis (default)
c 1:	beta-plane (linear varying Coriolis parameter)
c 2:	f-plane (constant Coriolis parameter)
c
c please note that in case of icor=1 or 2 also the parameter dlat (the average 
c latitude of the basin) has to be set
c
c with spherical coordinates (isphe=1) the default is to always use
c Coriolis. If you really do not want to use Coriolis, then please set
c icor = -1. The parameter dlat is not needed.

	icor=nint(getpar('icor'))	!flag how to use Coriolis
	isphe=nint(getpar('isphe'))	!flag if coordinates are spherical

	rad = pi / 180.

	yc=0.
	ymin=ygv(1)
	ymax=ygv(1)
	do k=1,nkn
	  y = ygv(k)
	  yc = yc + y
	  ymin = min(y,ymin)
	  ymax = max(y,ymax)
	end do
	yc=yc/nkn

	aux1 = 0.
	aux2 = 0.

	if(icor.eq.1) then	!beta plane
		aux1 = omega2 * sin(dcor*rad)
		aux2 = omega2 * cos(dcor*rad) / rearth
	else if(icor.eq.2) then	!f plane
		aux1 = omega2 * sin(dcor*rad)
		aux2 = 0.
	end if

	write(6,*) 'pi , dcor    : ',pi,dcor
	write(6,*) 'f_0, beta    : ',aux1,aux2
	write(6,*) 'yc,ymin,ymax : ',yc,ymin,ymax

	do ie=1,nel
	  ym=0.
	  do ii=1,3
	    ym=ym+ygv(nen3v(ii,ie))
	  end do
	  ym=ym/3.
	  if( isphe .eq. 0 ) then	!cartesian
	    fcorv(ie)=aux1+aux2*(ym-yc)
	  else				!spherical
	    if( icor .lt. 0 ) then	! -> do not use
	      fcorv(ie) = 0.
	    else
	      !fcorv(ie) = omega2*cos(ym*rad)	!BUG
	      fcorv(ie) = omega2*sin(ym*rad)
	    end if
	  end if
	end do

	end

c*****************************************************************

	subroutine check_coriolis

c checks coriolis parameter

	implicit none

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	real fcorv(1)
	common /fcorv/fcorv 

	integer ie
	real f,fmax

	fmax = 2.0 * 0.729E-4

	do ie=1,nel
	  f = fcorv(ie)
	  if( fmax - f .lt. 0. ) then
	    write(6,*) ie,f,fmax
	    stop 'error stop check_coriolis: error in array'
	  end if
	end do

	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine inicfil(name,var,nvar)

c initializes nodal value variable from file

	implicit none

	include 'param.h'

	character*(*) name		!name of variable
	real var(nlvdim,nkndim,1)	!variable to set
        integer nvar

	integer it
	integer nlvdi
	character*80 file

	call getfnm(name,file)
	if( file .eq. ' ' ) return	!nothing to initialize

	it = -1
	nlvdi = nlvdim
	call read_scalar(file,it,nlvdi,var,nvar)

	end

c*****************************************************************

	subroutine inic2fil(name,var,nvar)

c initializes nodal value variable from file (2D version)

	implicit none

	include 'param.h'

	character*(*) name		!name of variable
	real var(nkndim,1)		!variable to set
        integer nvar

	integer it
	integer nlvdi
	character*80 file

	call getfnm(name,file)
	if( file .eq. ' ' ) return	!nothing to initialize

	it = -1
	nlvdi = 1
	call read_scalar(file,it,nlvdi,var,nvar)

	end

c*****************************************************************

	subroutine read_scalar(file,it,nlvdi,var,nvar)

c reads scalar value from file

c format of file:
c
c read(nb) it,nkn,lmax,ivars
c read(nb) (hlv(l),l=1,lmax)	!only if lmax > 1
c do ivar=1,ivars
c    read(nb) ((var(l,k,ivar),l=1,lmax),k=1,nkn)
c end do
c
c it is the time of the data record - should be equal to the parameter
c	passed into the routine. If there are more records, it is used to
c	choose the right one. If you just want to take the first record
c	then pass -1 to the subroutine
c nkn is number of nodal values - must be the same as in FEM simulation
c ivars is number of variables  - must be the same as nvar (parameter)
c lmax is number of layers given
c	lmax = 0 or 1	=>	only one layer given, second record not read
c	lmax > 1	=>	1 or more layers are given
c				in this case the hlv(l) read must be identical
c				to the one in the FEM simulation
c if lmax is less than the number of available layers 
c in the simulation (nlvdi), the value of var(lmax,k,ivar) is 
c used to initialize the underlying layers

	implicit none

	include 'param.h'

	character*(*) file		!file to read
	integer it			!time to look for, time found (return)
	integer nlvdi			!vertical dimension of var
	real var(nlvdi,nkndim,1)	!variable to contain scalar
        integer nvar			!how many variables

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real hlv(1)
	common /hlv/hlv

	integer nb,irec,itsearch
	integer nkk,lmax
	integer l,k
        integer ivars,ivar
	real val
	real rlaux(nlvdim)

	integer ifileo

c-------------------------------------------------------
c open file
c-------------------------------------------------------

	nb = ifileo(0,file,'unform','old')
	if( nb .le. 0 ) goto 97

c-------------------------------------------------------
c loop on file
c-------------------------------------------------------

        ivar = 0
	itsearch = it
	it = it + 1			!make different

	do while( itsearch .ne. it )

c	  -------------------------------------------------------
c	  read first record
c	  -------------------------------------------------------

	  irec = 1
	  read(nb,err=90) it,nkk,lmax,ivars

	  if( nkk .ne. nkn .or. lmax .gt. nlvdi ) goto 99
	  if( ivars .ne. nvar ) goto 96

c	  -------------------------------------------------------
c	  read second record (only if lmax > 0)
c	  -------------------------------------------------------

	  irec = 2
	  if( lmax .gt. 1 ) then          !changed from 0 to 1 (5.3.2004) LMAX
	    read(nb,err=90) (rlaux(l),l=1,lmax)
	  else
	    lmax = 1
	  end if

c	  -------------------------------------------------------
c	  read data records
c	  -------------------------------------------------------

	  irec = 3
          do ivar=1,nvar
	    read(nb,err=90) ((var(l,k,ivar),l=1,lmax),k=1,nkn)
          end do

	  if( itsearch .eq. -1 ) itsearch = it		!take first record
	end do

c-------------------------------------------------------
c done -> close file
c-------------------------------------------------------

	close(nb)

c-------------------------------------------------------
c initialize the other levels if only surface is given
c-------------------------------------------------------

	if( lmax .gt. 1 ) then
	  do l=1,lmax
	    if( hlv(l) .ne. rlaux(l) ) goto 98
	  end do
	end if

	if( lmax .lt. nlvdi ) then
          do ivar=1,nvar
	    do k=1,nkn
	      val = var(lmax,k,ivar)
	      do l=lmax+1,nlvdi
	        var(l,k,ivar) = val
	      end do
	    end do
          end do
	end if

c-------------------------------------------------------
c end of routine
c-------------------------------------------------------

	write(6,*) 'Variable succesfully read from file '
	write(6,*) file
	write(6,*) it,nkk,lmax,ivars

	return
   90	continue
	write(6,*) 'read error in record = ',irec,' ivar = ',ivar
	write(6,*) '... reading file',file
	stop 'error stop inicfil'
   96	continue
	write(6,*) 'ivars not compatible with nvar: ',ivars,nvar
	stop 'error stop inicfil'
   97	continue
	write(6,*) 'Cannot open file ',file
	stop 'error stop inicfil'
   98	continue
	write(6,*) 'levels are not the same from init file ',file
	write(6,*) (hlv(l),l=1,lmax)
	write(6,*) (rlaux(l),l=1,lmax)
	stop 'error stop inicfil'
   99	continue
	write(6,*) 'parameters are not the same from init file ',file
	write(6,*) 'nkn, lmax from file  : ',nkk,lmax
	write(6,*) 'nkn, lmax from model : ',nkn,nlvdim
	stop 'error stop inicfil'
	end

c*******************************************************************

	subroutine init_z(const)

c initializes variables

	implicit none

	real const		!constant z value to impose

	include 'param.h'

	character*(80) name

c--------------------------------------------------------
c get name of file
c--------------------------------------------------------

        call getfnm('bound',name)

c--------------------------------------------------------
c initialize from file or with constant
c--------------------------------------------------------

	if(name.ne.' ') then
	  call init_file_z(name)
	else
	  call init_const_z(const)
	end if

c--------------------------------------------------------
c transfer nodal values to element values
c--------------------------------------------------------

	call setzev

c--------------------------------------------------------
c end of routine
c--------------------------------------------------------

	end

c*****************************************************************

	subroutine init_const_z(const)

c initializes water level with constant

	implicit none

	real const		!constant z value to impose

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	real znv(1),zov(1)
	common /znv/znv, /zov/zov
        real zenv(3,1), zeov(3,1)
        common /zenv/zenv, /zeov/zeov

	integer k,ie,ii

	do k=1,nkn
	  znv(k) = const
	  !zov(k) = const                !FIXME -> not needed
	end do

	do ie=1,nel
	  do ii=1,3
	    zenv(ii,ie) = const
	    !zeov(ii,ie) = const         !FIXME -> not needed
	  end do
	end do

        write(6,*) 'Water levels initialized with constant z = ',const

	end

c*******************************************************************

	subroutine init_file_z(name)

c initializes water level from file

	implicit none

	character*(*) name	!file name

	include 'param.h'

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	real znv(1)
	common /znv/znv

	integer nb,k
	integer ifileo

	nb=ifileo(10,name,'unform','old')
	if(nb.le.0) goto 81
	read(nb,err=82,end=82) (znv(k),k=1,nkn)
	close(nb)

	write(6,*) 'Initial water levels read from file : '
	write(6,*) name

	return
   81	continue
	write(6,*) 'Error opening initial water level file :'
	write(6,*) name
	write(6,*) 'on unit ',nb
	stop 'error stop init_file_z'
   82	continue
	write(6,*) 'Error reading from initial water level file'
	write(6,*) name
	stop 'error stop init_file_z'
	end

c*******************************************************************

	subroutine init_uvt

c initializes transport in levels

	implicit none

	include 'param.h'

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	real utlnv(nlvdim,1), vtlnv(nlvdim,1)
	common /utlnv/utlnv, /vtlnv/vtlnv

	integer ie,l

	do ie=1,nel
	  do l=1,nlvdim
	    utlnv(l,ie) = 0.
	    vtlnv(l,ie) = 0.
	  end do
	end do

	end

c*******************************************************************

	subroutine init_z0b

c initializes bottom roughness z0bn(k)

	implicit none

	include 'param.h'

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        real z0bk(nkndim)                   !bottom roughenss on nodes
        common /z0bk/z0bk

	integer k

	do k = 1,nkn
	  z0bk(k) = 0.03 * 0.03       !ggu
	end do

	end

c*******************************************************************

	subroutine set_spherical

	implicit none

	integer isphe
	real getpar

	isphe = nint(getpar('isphe'))
	call set_coords_ev(isphe)

	end

c*******************************************************************

