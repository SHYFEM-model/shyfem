
!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2018  Georg Umgiesser
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

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
c 07.03.2000	ggu	arrays for vertical turbulent diffusion coefficient
c 07.03.2000	ggu	useless parts commented
c 20.06.2000	ggu	useless parts deleted
c 21.06.2000	ggu	dzreg introduced
c 08.08.2000	ggu	hlvmin is now percentage of last layer thickness
c 25.03.2003	ggu	voldist adjusted for 3D, new routine surinl
c 07.08.2003	ggu	deleted sp212
c 09.08.2003	ggu	cleaned up (sp211 is completely different)
c 14.08.2003	ggu	error in check_ilevels (integer instead of real)
c 14.08.2003	ggu	error in check_ilevels (must test .lt.)
c 14.08.2003	ggu	new routine set_depth -> transferred from cstcheck
c 14.08.2003	ggu	renamed sp211 into init_vertical and init_others
c 14.08.2003	ggu	new routines init_z, init_const_z, init_file_z
c 14.08.2003	ggu	new routines init_uvt
c 02.09.2003	ggu	routine inicfil restructured for lmax=0
c 03.09.2003	ggu	some more checks for check_ilevels
c 04.03.2004	ggu	change in inicfil() for more variables
c 05.03.2004	ggu	LMAX - changed bound from 0 to 1
c 02.09.2004	ggu	rdrst transfered to subrst.f
c 15.03.2005	ggu	set_austausch() and check_austausch() eliminated
c 05.04.2005	ggu	minimum levels (for baroclinic terms) -> set_min_levels
c 28.11.2005	ggu	in set_depth use makehkv to compute hkv
c 16.02.2006	ggu	describe file format in inicfil()
c 15.11.2006	ggu	new parameters to construct streched vert. coordinates
c 27.08.2007	ccf	variable coriolis from spherical coordinates (isphe)
c 17.03.2008	ggu	better error message if missing levels
c 07.04.2008	ggu	deleted surinl, volinl, voldist
c 12.11.2008	ggu	set_last_layer() rewritten, new adjust_k_depth()
c 06.12.2008	ggu&ccf	small bug fix in set_last_layer()
c 21.01.2009	ggu	cleaned up, new read_scalar()
c 27.01.2009	ggu	mzreg and hlvmax deleted (not used)
c 24.03.2009	ggu	bug fix: in set_last_layer() do not adjust for 2D
c 21.04.2009	ggu	new routine inic2fil()
c 23.03.2010	ggu	changed v6.1.1
c 09.04.2010	ggu	changed v6.1.3
c 28.09.2010	ggu	bug fix in init_coriolis() for isphe=1
c 08.10.2010	ggu	bug fix in init_coriolis() -> ym not set for isphe=1
c 16.12.2010	ggu	big restructering for sigma levels (look for bsigma)
c 27.01.2011	ggu	changed VERS_6_1_17
c 21.02.2011	ggu	error check for dzreg, nsigma and levels
c 01.03.2011	ggu	changed VERS_6_1_20
c 23.03.2011	ggu	new routine adjust_spherical(), error check isphe
c 14.04.2011	ggu	changed VERS_6_1_22
c 24.08.2011	ggu	eliminated hbot for sigma due to run time error
c 25.10.2011	ggu	hlhv eliminated
c 04.11.2011	ggu	adapted for hybrid coordinates
c 10.11.2011	ggu	changed VERS_6_1_36
c 11.11.2011	ggu	bug fix in adjust_levels: for zeta levels set nlv
c 22.11.2011	ggu	changed VERS_6_1_37
c 12.12.2011	ggu	changed VERS_6_1_39
c 29.08.2012	ggu	changed VERS_6_1_56
c 10.05.2013	ggu	changed VERS_6_1_64
c 13.06.2013	ggu	changed VERS_6_1_65
c 23.08.2013	ggu	renamed adjust_k_depth() to make_hkv() -> to subdep
c 23.08.2013	ggu	depth routines to subdep.f
c 05.09.2013	ggu	in adjust_levels() allow for nlv==1
c 12.09.2013	ggu	changed VERS_6_1_67
c 25.10.2013	ggu	changed VERS_6_1_68
c 31.10.2014	ccf	initi_z0 for zos and zob
c 05.11.2014	ggu	changed VERS_7_0_5
c 05.12.2014	ggu	changed VERS_7_0_8
c 12.12.2014	ggu	changed VERS_7_0_9
c 23.12.2014	ggu	changed VERS_7_0_11
c 19.01.2015	ggu	changed VERS_7_1_2
c 19.01.2015	ggu	changed VERS_7_1_3
c 30.04.2015	ggu	changed VERS_7_1_9
c 21.05.2015	ggu	changed VERS_7_1_11
c 25.05.2015	ggu	file cleaned and prepared for module
c 05.06.2015	ggu	changed VERS_7_1_12
c 10.07.2015	ggu	changed VERS_7_1_50
c 13.07.2015	ggu	changed VERS_7_1_51
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 10.10.2015	ggu	changed VERS_7_3_2
c 05.11.2015	ggu	can now initialize z,u,v from file
c 28.04.2016	ggu	changed VERS_7_5_9
c 30.05.2016	ggu	changes in set_last_layer(), possible bug fix ilytyp==1
c 28.06.2016	ggu	coriolis computation changed -> yc=(ymax-ymin)/2.
c 05.10.2016	ggu	changed VERS_7_5_19
c 05.12.2017	ggu	changed VERS_7_5_39
c 24.01.2018	ggu	changed VERS_7_5_41
c 03.04.2018	ggu	changed VERS_7_5_43
c 19.04.2018	ggu	changed VERS_7_5_45
c 11.05.2018	ggu	new routine exchange_vertical() to compute global nlv,hlv
c 06.07.2018	ggu	changed VERS_7_5_48
c 14.02.2019	ggu	changed VERS_7_5_56
c 16.02.2019	ggu	changed VERS_7_5_60
c 13.03.2019	ggu	changed VERS_7_5_61
c
c notes :
c
c for information on hybrid levels see adjust_levels()
c
c**************************************************************

	subroutine init_vertical

c set up time independent vertical vectors

	use levels
	use basin, only : nkn,nel,ngr,mbw
	use shympi

	implicit none

	integer nlv_est,nlv_read,nlv_final
	integer nlv_e,nlv_k
	real hmax
	real, allocatable :: hlv_aux(:)

	write(6,*) 'setting up vertical structure'

c------------------------------------------------------------------
c sanity check
c------------------------------------------------------------------

	call get_hmax_global(hmax)
	write(6,*) 'maximum depth: ',hmax

	nlv_est = nlv
	call estimate_nlv(nlv_est,hmax)
	write(6,*) 'nlv,nlv_est,nlvdi: ',nlv,nlv_est,nlvdi

	call check_nlv

	if( nlv > 0 ) then
	  allocate(hlv_aux(nlv))
	  hlv_aux(1:nlv) = hlv(1:nlv)
	  call levels_hlv_init(0)
	end if
	call levels_init(nkn,nel,nlv_est)
	if( nlv > 0 ) then
	  hlv(1:nlv) = hlv_aux(1:nlv)
	  deallocate(hlv_aux)
	end if

c------------------------------------------------------------------
c levels read in from $levels section
c------------------------------------------------------------------

	call adjust_levels(hmax)	!sets hlv, hldv, nlv, sigma_info, etc.

c------------------------------------------------------------------
c set up layer vectors
c------------------------------------------------------------------

	call set_ilhv		!sets nlv, ilhv (elemental)
	call set_last_layer	!adjusts nlv, ilhv, hm3v
	call set_ilhkv		!sets ilhkv (nodal)
	call set_ilmkv		!sets ilmkv (nodal)
	call exchange_levels	!copies from other domains and sets nlv
	call set_ilmv		!sets ilmv (elemental)

c------------------------------------------------------------------
c compute final nlv
c------------------------------------------------------------------

	nlv_e = maxval(ilhv)
	nlv_k = maxval(ilhkv)
	nlv = max(nlv_e,nlv_k)

c------------------------------------------------------------------
c check data structure
c------------------------------------------------------------------

	nlv_final = nlv
	call levels_reinit(nlv_final)

	call check_vertical
	call shympi_set_hlv(nlv,hlv)

	write(6,*) 'init_vertical: nlvdi = ',nlvdi,'  nlv = ',nlv

c------------------------------------------------------------------
c end of routine
c------------------------------------------------------------------

	end

c**************************************************************

	subroutine exchange_vertical(nlv,hlv)

! exchanges vertical structure between domains
!
! to be deleted...

	use shympi

	implicit none

	integer nlv
	real hlv(nlv)

	integer ia
	real, parameter :: flag = -1.35472E+10
	real, allocatable :: hlvs(:,:)

        nlv_global = shympi_max(nlv)

	allocate(hlv_global(nlv_global))
	allocate(hlvs(nlv_global,n_threads))

	hlv_global = flag
	hlv_global(1:nlv) = hlv
	call shympi_gather(hlv_global,hlvs)

	do ia=1,n_threads
	  if( hlvs(nlv_global,ia) /= flag ) exit
	end do
	if( ia > n_threads ) then
	  write(6,*) 'error setting global nlv'
	  write(6,*) nlv_global,nlv
	  write(6,*) hlvs
	  stop 'error stop exchange_vertical: global nlv'
	end if

	hlv_global = hlvs(:,ia)

	!write(6,*) 'exchange_vertical: global nlv set: ',nlv,nlv_global
	!write(6,*) hlv_global

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

	use mod_diff_visc_fric
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

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

	visv=vistur
	difv=diftur

c------------------------------------------------------------------
c end of routine
c------------------------------------------------------------------

	end

c*****************************************************************

	subroutine check_eddy

c checks vertical eddy coefficient

	use mod_diff_visc_fric
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

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

	subroutine check_vertical

c checks arrays containing vertical structure

	use levels
	use shympi

	implicit none

	call check_nlv
	call check_hlv
	call check_levels
	call check_ilevels

        if(bmpi_debug) then
	  call shympi_check_2d_node(ilhkv,'ilhkv')
	  call shympi_check_2d_node(ilmkv,'ilmkv')
	  call shympi_check_2d_elem(ilhv,'ilhv')
	  call shympi_check_2d_elem(ilmv,'ilmv')
        end if

	end

c*****************************************************************

	subroutine check_hlv

	use levels

	implicit none

	integer l

	write(6,*) 'check_hlv: ',nlv,nlvdi
	write(6,'(5g14.6)') (hlv(l),l=1,nlv)

	end

c*****************************************************************

	subroutine check_nlv

c checks nlv and associated parameters

	use levels, only : nlvdi,nlv

	implicit none

	write(6,*) 'check_nlv : ',nlvdi,nlv

	if(nlv.gt.nlvdi) stop 'error stop check_nlv: level dimension'

	end

c*****************************************************************

	subroutine estimate_nlv(nlv_est,hmax)

c estimates maximum value for nlv

	use basin

	implicit none

	integer nlv_est		!nlv_read on entry, estimate on return
	real hmax

	integer ie,ii
	integer nsigma,nreg
	real hsigma,dzreg

	real getpar

	dzreg = getpar('dzreg')
	call get_sigma(nsigma,hsigma)

	nreg = 0
	if( dzreg > 0 ) nreg = hmax/dzreg

	!nlv_est = nlv_est + nsigma + nreg + 1
	nlv_est = nlv_est + nsigma + nreg
	nlv_est = max(nlv_est,1)

	end

c*****************************************************************

	subroutine adjust_levels(hmax)

c adjusts levels read in from $levels section
c
c creates hlv if not set
c from hlv creates hldv
c needs hm3v to compute hmax
c
c strategy:
c
c set hsigma < hmax to ask for hybrid levels
c set nsigma > 1 to ask for sigma coordinates
c set dzreg > 0. to ask for regular zeta levels
c
c only nsigma				only sigma levels
c only dzreg				only regular zeta levels
c nsigma, dzreg and hsigma		hybrid levels
c
c only sigma levels			only sigma levels
c only zeta levels			only zeta levels
c sigma levels, dzreg and hsigma	hybrid levels
c zeta levels, nsigma and hsigma	hybrid levels
c sigma and zeta levels and hsigma	hybrid levels

	use levels
	use basin

	implicit none

	real hmax

	logical bsigma,bhybrid,bzeta
	integer l,ie,ii,nsigma,second
	real dzreg,hl,fact,hsigma
	real hbot,htop

	real getpar

	write(6,*) 'adjust layer structure'

c--------------------------------------------------------------
c create hlv values
c--------------------------------------------------------------

	second = 2
	dzreg = getpar('dzreg')
	call get_sigma(nsigma,hsigma)
	write(6,*) 'nlv,nsigma,hsigma: ',nlv,nsigma,hsigma

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
	      call make_zeta_levels(0,0.,dzreg,nlv,hlv)
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

c--------------------------------------------------------------
c create hldv values
c--------------------------------------------------------------

	hldv(1)=hlv(1)
	do l=2,nlv
	  htop = hlv(l-1)
	  hbot = hlv(l)
	  if( l .eq. nsigma ) hbot = -1.
	  hldv(l) = hbot - htop
	end do

	write(6,*) 'adjust_levels: '
	write(6,*) 'nlv,nsigma,hsigma: ',nlv,nsigma,hsigma
	write(6,'(5g14.6)') 'hlv:  ',(hlv(l),l=1,nlv)
	write(6,'(5g14.6)') 'hldv: ',(hldv(l),l=1,nlv)

c--------------------------------------------------------------
c check hlv and hldv values
c--------------------------------------------------------------

	hbot = hlv(nlv)
	if( hbot /= -1. .and. hmax > hbot ) goto 87

	call check_levels

	write(6,*) 'finished adjusting layer structure ',nlv

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

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

c*****************************************************************

	subroutine check_levels

c checks arrays hlv and hldv

	use levels

	implicit none

c local
	logical bstop,bsigma
	integer l,nsigma,levmin
	real h,hd,fact,hsigma
	real hbot,htop

	call get_sigma(nsigma,hsigma)
	bsigma = nsigma .gt. 0

	bstop = .false.

c--------------------------------------------------------------
c check hlv values
c--------------------------------------------------------------

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

c--------------------------------------------------------------
c check hldv values
c--------------------------------------------------------------

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

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	end

c*****************************************************************

	subroutine exchange_levels

c exchanges level info with other domains - sets nlv

	use levels
	use shympi

	implicit none

	!call shympi_comment('exchanging ilhkv, ilmkv')
	call shympi_exchange_2d_node(ilhkv)
	call shympi_exchange_2d_node(ilmkv)
	!call shympi_barrier

	end

c*****************************************************************

	subroutine set_ilhv

c sets nlv and ilhv - only needs hm3v and hlv, hev is temporary array

	use levels
	use basin

	implicit none

c local
	logical bsigma
	integer ie,ii,l,lmax,nsigma
	real hsigma

	real h,hmax,hm
	real hev(nel)		!local

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

	  if( bsigma .and. h .le. hsigma ) then	!only sigma levels
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

c*****************************************************************

	subroutine set_ilhkv

c set ilhkv array - only needs ilhv

	use levels
	use basin
	use shympi

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

	if( shympi_partition_on_elements() ) then
          !call shympi_comment('shympi_elem: exchange ilhkv - max')
          call shympi_exchange_2d_nodes_max(ilhkv)
	else
          call shympi_exchange_2d_node(ilhkv)
	end if

	end

c*****************************************************************

	subroutine set_ilmkv

c set minimum number of levels for node

	use levels
	use basin
	use shympi

	implicit none

	integer ie,ii,k,l
	integer lmin,lmax

	ilmkv = huge(1)

	do ie=1,nel
	  l=ilhv(ie)
	  do ii=1,3
	    k=nen3v(ii,ie)
	    if(l.lt.ilmkv(k)) ilmkv(k)=l
	  end do
	end do

	if( shympi_partition_on_elements() ) then
          !call shympi_comment('shympi_elem: exchange ilmkv - min')
          call shympi_exchange_2d_nodes_min(ilmkv)
	else
          call shympi_exchange_2d_node(ilmkv)
	end if

	do k=1,nkn
	  lmin = ilmkv(k)
	  lmax = ilhkv(k)
	  if( lmin .gt. lmax .or. lmin .le. 0 ) then
	    stop 'error stop set_min_levels: lmin'
	  end if
	end do

	end

c*****************************************************************

	subroutine set_ilmv

c set minimum number of levels for elems

	use levels
	use basin

	implicit none

	integer ie,ii,k,lmin

	do ie=1,nel
	  lmin = huge(1)
	  do ii=1,3
	    k=nen3v(ii,ie)
	    if(lmin.gt.ilmkv(k)) lmin = ilmkv(k)
	  end do
	  ilmv(ie) = lmin
	end do

	end

c*****************************************************************

	subroutine check_ilevels

c checks arrays ilhv and ilhkv

	use levels
	use basin

	implicit none

	logical bsigma,bspure
	integer nsigma
	integer ie,ii,k,lmax,lk
	real hmax,hsigma

	call get_sigma(nsigma,hsigma)
	bsigma = nsigma .gt. 0

	hmax = 0.
	do ie=1,nel
	  lmax = ilhv(ie)
	  if( lmax .le. 0 ) goto 99
	  do ii=1,3
	    hmax = max(hmax,hm3v(ii,ie))
	  end do
	end do

	bspure = bsigma .and. hmax .le. hsigma	!pure sigma coordinates

	do ie=1,nel
	  lmax=ilhv(ie)
	  do ii=1,3
	    k=nen3v(ii,ie)
	    lk = ilhkv(k)
	    if( lk .le. 0 ) goto 98
	    if( lk .lt. lmax ) goto 98
	  end do
	  if( bspure .and. lmax .ne. nsigma ) goto 96
	  if( bsigma .and. lmax .lt. nsigma ) goto 96
	end do

	do k=1,nkn
	  lmax=ilhkv(k)
	  if( lmax .le. 0 ) goto 97
	  if( bspure .and. lmax .ne. nsigma ) goto 96
	  if( bsigma .and. lmax .lt. nsigma ) goto 96
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
	write(6,*) ie,lmax
	stop 'error stop check_ilevels: error in vertical structure (1)'
	end

c*****************************************************************

	subroutine set_last_layer

c sets last layer thickness
c
c adjusts nlv, hm3v, ilhv

	use levels
	use basin

	implicit none

	logical bwrite
	logical b2d,bsigma,binsigma
	logical bdepth,blayer
	integer ie,l,ii
	integer ihtot,lmax
	integer ilytyp
        integer ic1,ic2,ic3
	integer nsigma
	real h,hold,hnew,hlast
	real hlvmin
	real hmin,hmax,hsigma
	real hm
	!double precision hm

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

	call get_sigma(nsigma,hsigma)
	bsigma = nsigma .gt. 0

c	if( bsigma ) then		!sigma layers
c	  write(6,*) 'set_last_layer (nlv used) : ',nlv
c	  write(6,*) 'sigma layers detected'
c	  return
c	end if

c------------------------------------------------------------
c set hlvmin
c------------------------------------------------------------

c	ilytyp: 
c	  0=no adjustment  
c	  1=adjust to full layers (change depth)
c         2=adjust to full layers only if h<hlvmin (change depth)
c         3=add to last layer if h<hlvmin (keep depth but change layer)

	hlvmin = getpar('hlvmin')		!min percentage of last layer
	ilytyp = nint(getpar('ilytyp'))

	if( hlvmin .gt. 1. .or. hlvmin .lt. 0. ) goto 98

c------------------------------------------------------------
c adjust last layer thickness
c------------------------------------------------------------

	ic1 = 0
	ic2 = 0
	ic3 = 0

	ihtot = 0
	lmax = 0

	do ie=1,nel

	  !hm = 0.
	  !do ii=1,3
	  !  hm = hm + hm3v(ii,ie)
	  !end do
	  !hm = hm / 3.
	  hm = sum(hm3v(:,ie))/3.

	  l = ilhv(ie)
	  binsigma = l .le. nsigma
	  hold = hm				!original depth of element
	  hnew = hold
	  h = hold
	  if( l .gt. 1 ) h = h - hlv(l-1)	!actual last layer thickness
	  if( binsigma ) h = -hold*hldv(l)	!sigma layer
	  hlast = hldv(l)			!regular last layer thickness
	  hmin = hlvmin * hlast

	  if( l .gt. 1 .and. h .le. 0. ) goto 99

	  bdepth = .false.			!adjust depth? (implies blayer)
	  blayer = .false.			!adjust layer?

	  if( l .gt. 1 ) then			!only for more than 1 layer
	    if( ilytyp .eq. 0 ) then
c		no adjustment
	    else if( ilytyp .eq. 1 ) then
	      if( h < hlast ) then		!last layer not full layer
		bdepth = .true.
	      end if
	    else if( ilytyp .eq. 2 ) then
	      if( h < hmin ) then		!last layer too small
		bdepth = .true.
	      end if
	    else if( ilytyp .eq. 3 ) then
	      if( h < hmin ) then		!add to layer above
		blayer = .true.
	      end if
	    else
	      write(6,*) 'ilytyp = ',ilytyp
	      stop 'error stop set_last_layer: internal error (9)'
	    end if
	  end if

	  if( b2d .or. binsigma ) then		!dont for 2D or in sigma layer
	    bdepth = .false.
	    blayer = .false.
	  end if

	  if( blayer .or. bdepth ) l = l -1

	  if( bdepth ) then		!depth and layer have changed
	    h = hldv(l)
	    hnew = hlv(l)
	    hm3v(:,ie) = hnew
	  else if( blayer ) then	!only layer has changed
	    h = hldv(l) + h
	  end if

	  ilhv(ie) = l

	  if( blayer .or. bdepth ) then		!write to terminal
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
	write(6,*) 'ie,l,hold,h: ',ie,l,hold,h
	write(6,*) 'nsigma,hsigma: ',nsigma,hsigma
	write(6,*) 'hlv: ',nlv,(hlv(l),l=1,nlv)
	if( nsigma > 0 .and. hold < 0. ) then
	  write(6,*) 'cannot yet handle salt marshes with sigma layers'
	end if
	stop 'error stop set_last_layer: layer depth 0 (internal error)'
	end

c*****************************************************************

	subroutine init_coriolis

c sets coriolis parameter

	use mod_internal
	use basin
        use coordinates
	use shympi

	implicit none

	real omega2	!double frequency of earth rotation
	parameter ( omega2 = 2.0 * 0.729E-4 )
	real rearth	!radius of earth
	parameter ( rearth = 6371000. )

	include 'mkonst.h'
	include 'pkonst.h'

	logical bgeo
	integer k,ie,ii
	integer icor
	integer isphe
	real yc,ym,y,ymin,ymax,dlat
	real aux1,aux2,rad
	real getpar
	real, dimension(nkn)	:: yaux

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
c Coriolis with latitude read from basin. If you really do not want to 
c use Coriolis, then please set icor = -1. The parameter dlat is not needed.

c with cartesian coordinates (isphe=0) the default is to use a constant 
c latitude (dlat) for Coriolis (icor > 0). If you want a spatially varying 
c Coriolis parameter you have to convert the cartesian coordinates to 
c spherical setting the basin projection (iproj > 0)

	icor=nint(getpar('icor'))	!flag how to use Coriolis
	call get_coords_ev(isphe)

	rad = pi / 180.
	dlat = dcor			! average latitude
	fcorv = 0.

	if( icor < 0 ) return		! no coriolis

	bgeo = ( isphe .eq. 1 .or. iproj .ne. 0 )  !use geographical coords

	if( bgeo ) then
	  yaux = ygeov
	else
	  yaux = ygv
	end if
       
! next is handled differently - shympi FIXME - might break compatibility

	ymin = shympi_min(yaux)
	ymax = shympi_max(yaux)
	yc = (ymax+ymin)/2.
	!yc   = sum(yaux)/nkn

	if( bgeo ) dlat = yc		! get directly from basin

	aux1 = 0.
	aux2 = 0.

	if(icor.eq.1) then	!beta plane
		aux1 = omega2 * sin(dlat*rad)
		aux2 = omega2 * cos(dlat*rad) / rearth
	else if(icor.eq.2) then	!f plane
		aux1 = omega2 * sin(dlat*rad)
		aux2 = 0.
	end if

	fcor = aux1		! coriolis value of average latitude

	write(6,*) 'pi, dlat     : ',pi,dlat
	write(6,*) 'icor, fcor   : ',icor,fcor
	write(6,*) 'f_0, beta    : ',aux1,aux2
	write(6,*) 'yc,ymin,ymax : ',yc,ymin,ymax

	do ie=1,nel
	  ym=0.
	  do ii=1,3
	    ym=ym+yaux(nen3v(ii,ie))
	  end do
	  ym=ym/3.
	  if( bgeo ) then			!spherical
	    !fcorv(ie) = omega2*cos(ym*rad)	!BUG
	    fcorv(ie) = omega2*sin(ym*rad)
	  else if( isphe .eq. 0 ) then		!cartesian
  	    fcorv(ie)=aux1+aux2*(ym-yc)
	  else
	    write(6,*) 'isphe = ',isphe
	    if( isphe .eq. -1 ) write(6,*) '...not initialized'
	    stop 'error stop init_coriolis: value for isphe not allowed'
	  end if
	end do

	end

c*****************************************************************

	subroutine check_coriolis

c checks coriolis parameter

	use mod_internal
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer ie
	real f,fmax

	fmax = 2.0 * 0.729E-4

	do ie=1,nel
	  f = fcorv(ie)
	  if( fmax - abs(f) .lt. 0. ) then
	    write(6,*) ie,f,fmax
	    stop 'error stop check_coriolis: f too big'
	  end if
	end do

	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine inicfil(name,var,nvar)

c initializes nodal value variable from file

	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	character*(*) name		!name of variable
	real var(nlvdi,nkn,1)		!variable to set
        integer nvar

	integer it
	character*80 file

	call getfnm(name,file)
	if( file .eq. ' ' ) return	!nothing to initialize

	write(6,*) 'this call is not supported anymore...'
	write(6,*) name
	write(6,*) file
	write(6,*) 'please use tracer_file_init()'
	stop 'error stop inicfil: internal error - unsupported call'

	end

c*****************************************************************

	subroutine inic2fil(name,var,nvar)

c initializes nodal value variable from file (2D version)

	use basin, only : nkn,nel,ngr,mbw

	implicit none

	character*(*) name		!name of variable
	real var(nkn,1)			!variable to set
        integer nvar

	integer it
	integer nlvdi
	character*80 file

	call getfnm(name,file)
	if( file .eq. ' ' ) return	!nothing to initialize

	write(6,*) 'this call is not supported anymore...'
	write(6,*) name
	write(6,*) file
	write(6,*) 'please use tracer_file_init()'
	stop 'error stop inic2fil: internal error - unsupported call'

	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine init_z(zconst)

c initializes water levels (znv and zenv)

	implicit none

	real zconst		!constant z value to impose

	character*80 name
	logical rst_use_restart

c--------------------------------------------------------
c see if we already have hydro restart data
c--------------------------------------------------------

	if( rst_use_restart(1) ) return

c--------------------------------------------------------
c get name of file
c--------------------------------------------------------

        call getfnm('zinit',name)

c--------------------------------------------------------
c initialize from file or with constant
c--------------------------------------------------------

	if(name.ne.' ') then
	  call init_file_z(name)
	else
	  call init_const_z(zconst)
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

	subroutine init_const_z(zconst)

c initializes water level with constant

	use mod_hydro
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	real zconst		!constant z value to impose

	integer k,ie,ii

	znv = zconst
	zenv = zconst

        write(6,*) 'Water levels initialized with constant z = ',zconst

	end

c*******************************************************************

	subroutine init_file_z(name)

c initializes water level from file

	use mod_hydro
	use basin, only : nkn,nel,ngr,mbw
	use intp_fem_file

	implicit none

	character*(*) name	!file name

	integer nb,k
	integer nvar,nintp,np,lmax,ibc
	integer idzeta
	double precision dtime
	integer nodes(1)
	real zconst(1)
	character*10 what

	call get_first_dtime(dtime)
        nodes = 0
        nvar = 1
        nintp = 2
        np = nkn
        lmax = 0
        ibc = 0                         !no lateral boundary
        what = 'zeta init'
        zconst = 0.

	write(6,*) 'Initializing water levels...'
        call iff_init(dtime,name,nvar,np,lmax,nintp
     +                          ,nodes,zconst,idzeta)
        call iff_set_description(idzeta,ibc,what)

        lmax = 1
        call iff_read_and_interpolate(idzeta,dtime)
        call iff_time_interpolate(idzeta,dtime,1,np,lmax,znv)

	call iff_forget_file(idzeta)

	write(6,*) 'Initial water levels read from file : '
	write(6,*) trim(name)

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine init_uvt

c initializes transport in levels

	use mod_hydro

	implicit none

	character*80 name
	logical rst_use_restart

c--------------------------------------------------------
c see if we already have hydro restart data
c--------------------------------------------------------

	if( rst_use_restart(1) ) return

c--------------------------------------------------------
c get name of file
c--------------------------------------------------------

        call getfnm('uvinit',name)

c--------------------------------------------------------
c initialize from file or with constant
c--------------------------------------------------------

	if(name.ne.' ') then
	  call init_file_uv(name)
	  call vtot
	else
	  utlnv = 0.
	  vtlnv = 0.
	  call ttov
	end if

c--------------------------------------------------------
c end of routine
c--------------------------------------------------------

	end

c*******************************************************************

	subroutine init_file_uv(name)

c initializes water level from file

	use mod_hydro
	use mod_hydro_vel
	use mod_hydro_print
	use basin, only : nkn,nel,ngr,mbw
	use intp_fem_file
	use levels

	implicit none

	character*(*) name	!file name

	integer nb,k
	integer nvar,nintp,np,lmax,ibc
	integer ntype,iformat
	integer idvel
	double precision dtime
	integer nodes(1)
	real uvconst(2)
	character*10 what

	integer fem_file_regular

	call get_first_dtime(dtime)
        nodes = 0
        nvar = 2
        nintp = 2
        np = nkn			!velocities must be on node - relax later
        lmax = nlvdi
        ibc = 0                         !no lateral boundary
        what = 'uv init'
        uvconst = 0.

	write(6,*) 'Initializing velocities...'

	call fem_file_test_formatted(name,np,nvar,ntype,iformat)
	if( nvar /= 2 ) then
	  write(6,*) 'expecting 2 variables, found ',nvar
	  write(6,*) 'read error from file ',trim(name)
	  stop 'error stop init_file_uv: nvar'
	end if
	if( fem_file_regular(ntype) > 0 ) then
	  np = nel
	else if( np /= nkn .and. np /= nel ) then
	  write(6,*) 'unexpected value for np found ',np
	  write(6,*) 'possible values: ',nkn,nel
	  write(6,*) 'read error from file ',trim(name)
	  stop 'error stop init_file_uv: np'
	end if

        call iff_init(dtime,name,nvar,np,lmax,nintp
     +                          ,nodes,uvconst,idvel)
        call iff_set_description(idvel,ibc,what)

        call iff_read_and_interpolate(idvel,dtime)
	if( np == nkn ) then
          call iff_time_interpolate(idvel,dtime,1,np,lmax,uprv)
          call iff_time_interpolate(idvel,dtime,2,np,lmax,vprv)
	  call prtouv
	else
          call iff_time_interpolate(idvel,dtime,1,np,lmax,ulnv)
          call iff_time_interpolate(idvel,dtime,2,np,lmax,vlnv)
	end if

	call make_prvel

	call iff_forget_file(idvel)

	write(6,*) 'Initial velocities read from file : '
	write(6,*) trim(name)

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine init_z0

c initializes surface z0sk(k) and bottom z0bn(k) roughness 

	use mod_roughness

	implicit none

	z0sk = z0sini
	z0bk = z0bini

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine set_spherical

	implicit none

	integer isphe
	real getpar

	isphe = nint(getpar('isphe'))
	call set_coords_ev(isphe)

	end

c*******************************************************************

	subroutine adjust_spherical

	implicit none

	integer isphe
	real rsphe

	call get_coords_ev(isphe)
	rsphe = isphe
	call putpar('isphe',rsphe)

	end

c*******************************************************************

        subroutine print_spherical

	use shympi

        implicit none

        integer isphe

	if( .not. shympi_output() ) return

	call get_coords_ev(isphe)

        write(6,*) 'setting for coordinates: isphe = ',isphe
        if( isphe .eq. 0 ) then
          write(6,*) 'using cartesian coordinates'
        else
          write(6,*) 'using lat/lon coordinates'
        end if

        end

c*******************************************************************

