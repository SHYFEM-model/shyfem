!
! $Id: newini.f,v 1.45 2009-11-18 16:50:37 georg Exp $
!
! routines for initialization
!
! contents :
!
! subroutine inicfil(name,var,nvar)
!				initializes nodal value variable from file
! subroutine inic2fil(name,var,nvar)
!				initializes nodal value variable from file (2D)
!
! revision log :
!
! 20.08.1998	ggu	initialization routines copied from subn11/new11
! 20.08.1998	ggu	new routine inicfil to init nodal values from file
! 06.11.1998	ggu	computation of hkv/hev taken out of sp211 -> cstset
! 22.10.1999	ggu	igrv, ngrv eliminated (subst by ilinkv, lenkv)
! 07.03.2000    ggu     arrays for vertical turbulent diffusion coefficient
! 07.03.2000    ggu     useless parts commented
! 20.06.2000    ggu     useless parts deleted
! 21.06.2000    ggu     dzreg introduced
! 08.08.2000    ggu     hlvmin is now percentage of last layer thickness
! 25.03.2003    ggu     voldist adjusted for 3D, new routine surinl
! 07.08.2003    ggu     deleted sp212
! 09.08.2003    ggu     cleaned up (sp211 is completely different)
! 14.08.2003    ggu     error in check_ilevels (integer instead of real)
! 14.08.2003    ggu     error in check_ilevels (must test .lt.)
! 14.08.2003    ggu     new routine set_depth -> transferred from cstcheck
! 14.08.2003    ggu     renamed sp211 into init_vertical and init_others
! 14.08.2003    ggu     new routines init_z, init_const_z, init_file_z
! 14.08.2003    ggu     new routines init_uvt
! 02.09.2003    ggu     routine inicfil restructured for lmax=0
! 03.09.2003    ggu     some more checks for check_ilevels
! 04.03.2004    ggu     change in inicfil() for more variables
! 05.03.2004    ggu     LMAX - changed bound from 0 to 1
! 02.09.2004    ggu     rdrst transfered to subrst.f
! 15.03.2005    ggu     set_austausch() and check_austausch() eliminated
! 05.04.2005    ggu     minimum levels (for baroclinic terms) -> set_min_levels
! 28.11.2005    ggu     in set_depth use makehkv to compute hkv
! 16.02.2006    ggu     describe file format in inicfil()
! 15.11.2006    ggu     new parameters to construct streched vert. coordinates
! 27.08.2007    ccf     variable coriolis from spherical coordinates (isphe)
! 17.03.2008    ggu     better error message if missing levels
! 07.04.2008    ggu     deleted surinl, volinl, voldist
! 12.11.2008    ggu     set_last_layer() rewritten, new adjust_k_depth()
! 06.12.2008    ggu&ccf small bug fix in set_last_layer()
! 21.01.2009    ggu	cleaned up, new read_scalar()
! 27.01.2009    ggu	mzreg and hlvmax deleted (not used)
! 24.03.2009    ggu	bug fix: in set_last_layer() do not adjust for 2D
! 21.04.2009    ggu	new routine inic2fil()
! 28.09.2010    ggu	bug fix in init_coriolis() for isphe=1
! 08.10.2010    ggu	bug fix in init_coriolis() -> ym not set for isphe=1
! 16.12.2010    ggu	big restructering for sigma levels (look for bsigma)
! 21.02.2011    ggu	error check for dzreg, nsigma and levels
! 23.03.2011    ggu	new routine adjust_spherical(), error check isphe
! 24.08.2011    ggu	eliminated hbot for sigma due to run time error
! 25.10.2011    ggu	hlhv eliminated
! 04.11.2011    ggu	adapted for hybrid coordinates
! 11.11.2011    ggu	bug fix in adjust_levels: for zeta levels set nlv
! 23.08.2013    ggu	renamed adjust_k_depth() to make_hkv() -> to subdep
! 23.08.2013    ggu	depth routines to subdep.f
! 05.09.2013    ggu	in adjust_levels() allow for nlv==1
! 31.10.2014    ccf	initi_z0 for zos and zob
! 25.05.2015    ggu	file cleaned and prepared for module
! 05.11.2015    ggu	can now initialize z,u,v from file
! 28.06.2016    ggu	coriolis computation changed -> yc=(ymax-ymin)/2.
!
! notes :
!
! for information on hybrid levels see adjust_levels()
!
!**************************************************************
!--------------------------------------------------------------
        module initialize_par
!--------------------------------------------------------------
        contains
!--------------------------------------------------------------

	subroutine check_nlv

! checks nlv and associated parameters

	use levels, only : nlvdi,nlv

	implicit none

	write(6,*) 'check_nlv : ',nlvdi,nlv

	if(nlv.gt.nlvdi) stop 'error stop check_nlv: level dimension'

	end

!*****************************************************************

	subroutine estimate_nlv(nlv_est,hmax)

! estimates maximum value for nlv

	use basin
        use para
        use sigma_admin

	implicit none

	integer nlv_est		!nlv_read on entry, estimate on return
	real hmax		!maximum depth

	include 'param.h'

	integer ie,ii
	integer nsigma,nreg
        double precision hsigma,dzreg

	dzreg = getpar('dzreg')
	call get_sigma(nsigma,hsigma)

	nreg = 0
	if( dzreg > 0 ) nreg = hmax/dzreg

	!nlv_est = nlv_est + nsigma + nreg + 1
	nlv_est = nlv_est + nsigma + nreg
	nlv_est = max(nlv_est,1)

	end

!*****************************************************************

	subroutine adjust_levels(hmax)

! adjusts levels read in from $levels section
!
! creates hlv if not set
! from hlv creates hldv
! needs hm3v to compute hmax
!
! strategy:
!
! set hsigma < hmax to ask for hybrid levels
! set nsigma > 1 to ask for sigma coordinates
! set dzreg > 0. to ask for regular zeta levels
!
! only nsigma				only sigma levels
! only dzreg				only regular zeta levels
! nsigma, dzreg and hsigma		hybrid levels
!
! only sigma levels			only sigma levels
! only zeta levels			only zeta levels
! sigma levels, dzreg and hsigma	hybrid levels
! zeta levels, nsigma and hsigma	hybrid levels
! sigma and zeta levels and hsigma	hybrid levels

	use levels
	use basin
        use mpi_common_struct
        use para
        use sigma
        use sigma_admin

	implicit none

	real hmax

	logical bsigma,bhybrid,bzeta
	integer l,ie,ii,nsigma,second
        double precision dzreg,h,hl,fact,hsigma
        double precision hbot,htop

        if(my_id.eq.0) write(6,*) 'adjust layer structure'

!--------------------------------------------------------------
! create hlv values
!--------------------------------------------------------------

	second = 2
	dzreg = getpar('dzreg')
	call get_sigma(nsigma,hsigma)
        if(my_id.eq.0) write(6,*)'nlv,nsigma,hsigma: ',nlv,nsigma,hsigma

	bsigma = nsigma .gt. 0
	bhybrid = hsigma .lt. hmax
	bzeta = dzreg .gt. 0.

	if( nsigma .eq. 1 ) goto 92

	if( nlv .le. 0 ) then		! hlv not set -> must set

	  if( bsigma ) then		!sigma layers
	    if( nsigma > nlvdi ) goto 86
	    call make_sigma_levels(nsigma,hlv)
	    nlv = nsigma
	  end if

	  if( bhybrid .or. bzeta ) then	!zeta levels
	    if( bhybrid .and. .not. bzeta ) goto 99
	    if( bhybrid .and. .not. bsigma ) goto 90
	    nlv = nlvdi
	    if( bhybrid ) then
	      call make_zeta_levels(nsigma,hsigma,dzreg,nlv,hlv)
	    else if( .not. bsigma ) then	!no sigma layers given
	      call make_zeta_levels(0,0.d0,dzreg,nlv,hlv)
	    else
	      goto 89
	    end if
	  else if( .not. bsigma ) then		!just one layer
	    nlv = 1
	    hlv(1) = hsigma
	  end if

	else				!level section present

	  if( nlv .eq. 1 ) then		!just one level (barotropic)
	    if( bsigma .or. bhybrid ) goto 98
	  else if( hlv(second) .gt. hlv(1) ) then	!zeta layers given
	    if( bzeta ) goto 97
	    if( bsigma ) then		!put sigma layers on top
	      if( .not. bhybrid ) goto 88
	      if( hlv(1) .ne. hsigma ) goto 96
	      if( nsigma + nlv - 1 .gt. nlvdi ) goto 95
	      do l=1,nlv
	        hlv(nsigma+l-1) = hlv(l)
	      end do
	      call make_sigma_levels(nsigma,hlv)
	      hlv(nsigma) = hsigma
	      nlv = nlv + nsigma - 1
	    end if
	  else if( hlv(nlv) .eq. -1. ) then	!only sigma layers given
	    if( bsigma ) goto 91
	    nsigma = nlv
	    if( bhybrid ) then		!put zeta levels below
	      if( .not. bzeta ) goto 94
	      nlv = nlvdi
	      call make_zeta_levels(nsigma,hsigma,dzreg,nlv,hlv)
	    end if
	  else				!both sigma and zeta levels given
	    if( bsigma ) goto 91
	    if( bzeta ) goto 97
	    if( .not. bhybrid ) goto 88
	    do l=2,nlv
	      if( hlv(l) .gt. hlv(l-1) ) exit
	    end do
	    if( l > nlv ) goto 93	!not found
	    if( hlv(l-1) .eq. -1. ) then
	      nsigma = l-1
	      hlv(l-1) = hsigma
	    else if( hlv(l) .eq. hsigma ) then
	      nsigma = l
	    else
	      goto 93
	    end if
	    if( hlv(nsigma) .eq. hlv(nsigma+1) ) then	!double entry
	      do l=nsigma+1,nlv
		hlv(l-1) = hlv(l)
	      end do
	      nlv = nlv - 1
	    else if( hlv(nsigma) .ge. hlv(nsigma+1) ) then
	      goto 93
	    end if
	  end if

	end if

	call set_sigma(nsigma,hsigma)		!uses putpar
	call set_sigma_info(nlv,nsigma,hsigma)	!sets internal data structure

!--------------------------------------------------------------
! create hldv values
!--------------------------------------------------------------

	hldv(1)=hlv(1)
	do l=2,nlv
	  htop = hlv(l-1)
	  hbot = hlv(l)
	  if( l .eq. nsigma ) hbot = -1.
	  hldv(l) = hbot - htop
	end do

        if(my_id.eq.0) then
	  write(6,*) 'adjust_levels: '
	  write(6,*) 'nlv,nsigma,hsigma: ',nlv,nsigma,hsigma
	  write(6,'(5g14.6)') 'hlv:  ',(hlv(l),l=1,nlv)
	  write(6,'(5g14.6)') 'hldv: ',(hldv(l),l=1,nlv)
        end if

!--------------------------------------------------------------
! check hlv and hldv values
!--------------------------------------------------------------

	hbot = hlv(nlv)
	if( hbot /= -1. .and. hmax > hbot ) goto 87

	call check_levels

	write(6,*) 'finished adjusting layer structure ',nlv

!--------------------------------------------------------------
! end of routine
!--------------------------------------------------------------

	return
   86	continue
	write(6,*) 'nsigma,nlvdi: ',nsigma,nlvdi
	stop 'error stop adjust_levels: dimension too small'
   87	continue
	write(6,*) 'nlv,hlv(nlv),hmax: ',nlv,hlv(nlv),hmax
	stop 'error stop adjust_levels: hlv too low'
   88	continue
	write(6,*) 'hsigma: ',hsigma
	write(6,*) 'for hybrid levels hsigma must be set'
	stop 'error stop adjust_levels: hsigma'
   89	continue
	write(6,*) 'nsigma,dzreg: ',nsigma,dzreg
	write(6,*) 'with pure sigma layers cannot give dzreg'
	stop 'error stop adjust_levels: bsigma and dzreg'
   90	continue
	write(6,*) 'nsigma,hsigma: ',nsigma,hsigma
	write(6,*) 'with hsigma set nsigma must be > 0'
	stop 'error stop adjust_levels: nsigma and hsigma'
   91	continue
	write(6,*) 'nlv,nsigma: ',nlv,nsigma
	write(6,*) 'You cannot give both sigma levels and nsigma'
	stop 'error stop adjust_levels: nsigma'
   92	continue
	write(6,*) 'nsigma must be > 1'
	stop 'error stop adjust_levels: nsigma = 1'
   93	continue
	write(6,*) 'error in level structure'
	write(6,*) (hlv(l),l=1,nlv)
	stop 'error stop adjust_levels: hlv'
   94	continue
	write(6,*) 'for hybrid levels dzreg must be > 0'
	stop 'error stop adjust_levels: dzreg'
   95	continue
	write(6,*) 'nlv,nsigma,nlvdi: ',nlv,nsigma,nlvdi
	write(6,*) 'not enough space to add zeta levels'
	stop 'error stop adjust_levels: nlvdi'
   96	continue
	write(6,*) 'hlv(1),hsigma: ',hlv(1),hsigma
	write(6,*) 'for hybrid levels first zeta level must be hsigma'
	stop 'error stop adjust_levels: hsigma'
   97	continue
	write(6,*) 'nlv,dzreg: ',nlv,dzreg
	write(6,*) 'You cannot give both zeta levels and dzreg'
	stop 'error stop adjust_levels: dzreg'
   98	continue
	write(6,*) 'only one level given in level section'
	write(6,*) hlv(1)
	write(6,*) 'cannot have sigma levels...'
	stop 'error stop adjust_levels: nlv = 1'
   99	continue
	write(6,*) 'for hybrid levels dzreg must be > 0'
	stop 'error stop adjust_levels: dzreg = 0'
	end

!*****************************************************************

	subroutine check_levels

! checks arrays hlv and hldv

	use levels
        use sigma_admin

	implicit none

! common
	include 'param.h'

! local
	logical bstop,bsigma
	integer l,nsigma,levmin
        double precision h,hd,fact,hsigma
        double precision hbot,htop

	call get_sigma(nsigma,hsigma)
	bsigma = nsigma .gt. 0

	bstop = .false.

!--------------------------------------------------------------
! check hlv values
!--------------------------------------------------------------

	if( nsigma .gt. 0 ) then
	  h = hlv(nsigma)
	  if( h .ne. -1. .and. h .ne. hsigma ) then
	    write(6,*) h,hsigma
	    stop 'error stop check_levels: hsigma'
	  end if
	end if

	hbot = -hlv(1)
	do l=2,nsigma
	  htop = hbot
	  hbot = -hlv(l)
	  if( l .eq. nsigma ) hbot = 1.
	  if( hbot .le. htop ) then
	    write(6,*) 'Error in level values for level : ',l
	    write(6,*) '   hlv(l-1),hlv(l) :',hlv(l-1),hlv(l)
	    bstop = .true.
	  end if
	end do

	levmin = nsigma + 1
	if( levmin .eq. 1 ) levmin = 2
	do l=levmin,nlv
	  if( hlv(l) .le. hlv(l-1) ) then
	    write(6,*) 'Error in level values for level : ',l
	    write(6,*) '   hlv(l-1),hlv(l) :',hlv(l-1),hlv(l)
	    bstop = .true.
	  end if
	end do

	if( bstop ) stop 'error stop check_levels: level error'

!--------------------------------------------------------------
! check hldv values
!--------------------------------------------------------------

	hbot = hlv(1)
	do l=2,nlv
	  htop = hbot
	  hbot = hlv(l)
	  if( l .eq. nsigma ) hbot = -1.
	  hd=hbot-htop
	  if( hd .ne. hldv(l) ) then
	    write(6,*) 'Error in dlevel values for level : ',l
	    write(6,*) '   hd,hldv(l) :',hd,hldv(l)
	    bstop = .true.
	  end if
	  hbot = hlv(l)
	end do
	if( bstop ) stop 'error stop check_levels: dlevel error'

!--------------------------------------------------------------
! end of routine
!--------------------------------------------------------------

	end

!*****************************************************************

	subroutine set_ilhv

! sets nlv and ilhv - only needs hm3v and hlv, hev is temporary array

	use levels
	use basin
        use sigma_admin

	implicit none

	logical bsigma
	integer ie,ii,l,lmax,nsigma
        double precision hsigma

        double precision h,hmax,hm
        double precision hev(nel)		!local

	lmax=0
	hmax = 0.

	call get_sigma(nsigma,hsigma)
	bsigma = nsigma .gt. 0

	do ie=1,nel
	  hm = 0.
	  do ii=1,3
	    hm = hm + hm3v(ii,ie)
	  end do
	  hm = hm / 3.
	  hev(ie) = hm
	  hmax = max(hmax,hm)
	end do

	do ie=1,nel

	  h=hev(ie)

	  if( bsigma .and. h .le. hsigma ) then		!only sigma levels
	    l = nsigma
	  else
	    do l=nsigma+1,nlv
	      if(hlv(l).ge.h) exit
	    end do
	    if( l .gt. nlv ) goto 99
	  end if

	  ilhv(ie)=l
	  lmax = max(lmax,l)

	end do

	nlv = lmax

	write(6,*) 'finished setting ilhv and nlv'
	write(6,*) 'nsigma,hsigma: ',nsigma,hsigma
	write(6,*) 'nlv,lmax,hmax: ',nlv,lmax,hmax
	write(6,'(5g14.6)') (hlv(l),l=1,nlv)

	return
   99	continue
	write(6,*) ie,l,nlv,h,hlv(nlv)
	write(6,*) 'maximum basin depth: ',hmax
	write(6,*) 'maximum layer depth: ',hlv(nlv)
	write(6,'(5g14.6)') (hlv(l),l=1,nlv)
	stop 'error stop set_ilhv: not enough layers'
	end

!*****************************************************************

	subroutine set_ilhkv

! set ilhkv array - only needs ilhv

	use levels
	use basin

	implicit none

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

!*****************************************************************

	subroutine set_last_layer

! sets last layer thickness
!
! adjusts nlv, hm3v, ilhv

	use levels
	use basin
        use para
        use sigma_admin
        use fem_util

	implicit none

	include 'param.h'

	logical bwrite
	logical badjust,b2d,bsigma,binsigma
	integer ie,l,ii
	integer ihtot,lmax
	integer ilytyp
        integer ic1,ic2,ic3
	integer nsigma
        double precision h,hold,hnew,hlast
        double precision hlvmin
        double precision hmin,hmax,hsigma
        double precision hm
	!double precision hm

	bwrite = .false.

!------------------------------------------------------------
! see if 2D application
!------------------------------------------------------------

	b2d = nlv .eq. 1			!2D application

!------------------------------------------------------------
! check if sigma levels
!------------------------------------------------------------

	call get_sigma(nsigma,hsigma)
	bsigma = nsigma .gt. 0

!	if( bsigma ) then		!sigma layers
!	  write(6,*) 'set_last_layer (nlv used) : ',nlv
!	  write(6,*) 'sigma layers detected'
!	  return
!	end if

!------------------------------------------------------------
! set hlvmin
!------------------------------------------------------------

!	ilytyp: 
!	  0=no adjustment  
!	  1=adjust to full layers (change depth)
!         2=adjust to full layers only if h<hlvmin (change depth)
!         3=add to last layer if h<hlvmin (keep depth but change layer)

	hlvmin = getpar('hlvmin')		!min percentage of last layer
	ilytyp = nint(getpar('ilytyp'))

	if( hlvmin .gt. 1. .or. hlvmin .lt. 0. ) goto 98

!------------------------------------------------------------
! adjust last layer thickness
!------------------------------------------------------------

	ic1 = 0
	ic2 = 0
	ic3 = 0

	ihtot = 0
	lmax = 0

	do ie=1,nel

	  hm = 0.
	  do ii=1,3
	    hm = hm + hm3v(ii,ie)
	  end do
	  hm = hm / 3.

	  l = ilhv(ie)
	  binsigma = l .le. nsigma
	  hold = hm				!original depth of element
	  hnew = hold
	  h = hold
	  if( l .gt. 1 ) h = h - hlv(l-1)	!actual layer thickness
	  if( binsigma ) h = -hold*hldv(l)	!sigma layer
	  hlast = hldv(l)			!regular layer thickness
	  hmin = hlvmin * hlast

	  if( l .gt. 1 .and. h .le. 0. ) goto 99

	  badjust = .false.

	  if( l .gt. 1 ) then			!only for more than 1 layer
	    if( ilytyp .eq. 0 ) then
!		no adjustment
	    else if( h .le. hmin ) then		!take away layer
		badjust = .true.
		l = l-1
	    else if( ilytyp .eq. 1 ) then
		badjust = .true.
	    end if
	  end if

	  if( ilytyp .eq. 1 ) badjust = .true.	!adjust in any case
	  if( b2d ) badjust = .false.		!but never for 2D application
	  if( l .le. nsigma ) badjust = .false.	!we are in sigma layer

	  if( badjust ) then
	    if( ilytyp .le. 2 ) then		!adjust depth
		h = hldv(l)
		hnew = hlv(l)
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
	      write(6,'(2i8,3f12.2)') ie,ieext(ie),l,hnew,hold,h
	    end if
	    ihtot = ihtot + 1
	  end if

	  hlast = h			!finally set last layer thickness
	  lmax = max(lmax,l)

	  if( hlv(l) .ne. hnew ) ic1 = ic1 + 1
	  if( hldv(l) .ne. hlast ) ic2 = ic2 + 1
	  if( l .eq. 1 ) ic3 = ic3 + 1
	end do

!------------------------------------------------------------
! adjust nlv and write final statistics
!------------------------------------------------------------

	nlv=lmax

	write(6,*) 'set_last_layer (nlv used) : ',nlv
	write(6,*) 'Total number of elements with adjusted depth: ',ihtot
	write(6,*) 'Incomplete depth:     ',ic1
	write(6,*) 'Differing last layer: ',ic2
	write(6,*) 'One layer elements:   ',ic3

!------------------------------------------------------------
! end of routine
!------------------------------------------------------------

	return
   98	continue
	write(6,*) 'hlvmin is given in fraction of last layer'
	write(6,*) 'therefore:  0 <= hlvmin <= 1'
	write(6,*) 'hlvmin = ',hlvmin
	stop 'error stop set_last_layer: hlvmin'
   99	continue
	write(6,*) 'ie,l,hold,h: ',ie,l,hold,h
	write(6,*) 'nsigma,hsigma: ',nsigma,hsigma
	write(6,*) 'hlv: ',nlv,(hlv(l),l=1,nlv)
	if( nsigma > 0 .and. hold < 0. ) then
	  write(6,*) 'cannot yet handle salt marshes with sigma layers'
	end if
	stop 'error stop set_last_layer: layer depth 0 (internal error)'
	end

!*******************************************************************

	subroutine get_hmax_global(hmax)

	use basin
	use mpi_common_struct

	implicit none
        integer ierr
	real hmax,hout

        hmax = maxval(hm3v)
        if(bmpi) then
          call MPI_ALLREDUCE(hmax, hout, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
          hmax=hout
        end if

	end

!*******************************************************************
!--------------------------------------------------------------
        end module initialize_par
!--------------------------------------------------------------
