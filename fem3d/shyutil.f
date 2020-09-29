
!--------------------------------------------------------------------------
!
!    Copyright (C) 2010-2019  Georg Umgiesser
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

! utilities for SHY files
!
! revision log :
!
! 29.04.2010	ggu	new file from nosaver
! 07.05.2010	ggu	new routines qopen_nos_file(), checks ev initialization
! 15.12.2010	ggu	volume computation also for sigma layers
! 10.11.2011	ggu	new routines for hybrid levels, init_volume() changed
! 02.12.2011	ggu	bug fix for call to get_sigma_info() (missing argument)
! 10.02.2012	ggu	new routines to get initial/final time of records
! 25.01.2013	ggu	new routines nos_get_vars()
! 05.09.2013	ggu	new call to get_layer_thickness()
! 20.01.2014	ggu	new helper routines
! 23.09.2015	ggu	close files in nos_get_it_start() nos_get_it_end()
! 19.02.2016	ggu	changed VERS_7_5_2
! 28.04.2016	ggu	changed VERS_7_5_9
! 25.05.2016	ggu	changed VERS_7_5_10
! 30.05.2016	ggu	changed VERS_7_5_11
! 08.09.2016	ggu	new flag bforce to force output
! 05.10.2016	ggu	use zeps for computing volumes
! 21.03.2017	ggu	use new values for initialization of accum
! 21.03.2017	ggu	computation of std revised (no NaNs)
! 31.03.2017	ggu	changed VERS_7_5_24
! 11.07.2017	ggu	changed VERS_7_5_30
! 04.11.2017	ggu	changed VERS_7_5_34
! 14.11.2017	ggu	changed VERS_7_5_36
! 12.06.2018	ggu	bug fix for area in make_aver_3d()
! 06.07.2018	ggu	changed VERS_7_5_48
! 16.02.2019	ggu	changed VERS_7_5_60
! 22.09.2020    ggu     correct warnings for PGI compiler
!
!***************************************************************

!==================================================================
	module shyutil
!==================================================================

	implicit none

        real, save, allocatable :: vol3e(:,:)
        real, save, allocatable :: vol3k(:,:)
        real, save, allocatable :: vol2e(:)
        real, save, allocatable :: vol2k(:)
        real, save, allocatable :: areae(:)
        real, save, allocatable :: areak(:)

        real, save, allocatable :: shy_znv(:)
        real, save, allocatable :: shy_zenv(:,:)
        real, save, allocatable :: shy_zeta(:)

        integer, save, allocatable :: naccu(:,:)
        double precision, save, allocatable :: accum(:,:,:,:)
        double precision, save, allocatable :: std(:,:,:,:)
        double precision, save, allocatable :: dir(:,:,:,:,:)

        integer, parameter :: idir = 16
        real, parameter :: adir = 360. / idir

	double precision, parameter :: accum_high = 1.e+30
	double precision, save :: accum_init = 0.

	logical, save :: bzeta = .false.

!==================================================================
	contains
!==================================================================

	subroutine shyutil_init(nkn,nel,nlv)

	integer nkn,nel,nlv

	allocate(vol3e(nlv,nel))
	allocate(vol3k(nlv,nkn))
	allocate(vol2e(nel))
	allocate(vol2k(nkn))
	allocate(areae(nel))
	allocate(areak(nkn))

	allocate(shy_znv(nkn))
	allocate(shy_zenv(3,nel))
	allocate(shy_zeta(nel))

	vol3e = 1.
	vol3k = 1.
	vol2e = 1.
	vol2k = 1.
	areae = 1.
	areak = 1.

	shy_znv = 0.
	shy_zenv = 0.
	shy_zeta = 0.

	end subroutine shyutil_init

!***************************************************************

	subroutine shyutil_init_accum(avermode,nlv,nn,nvar,istep)

	integer avermode
	integer nlv,nn,nvar,istep

        allocate(naccu(0:nvar,istep))
        allocate(accum(nlv,nn,0:nvar,istep))
        allocate(std(nlv,nn,0:nvar,istep))

	if( avermode == 3 ) then	!min
	  accum_init = accum_high
	else if( avermode == 4 ) then	!max
	  accum_init = -accum_high
	else				!aver, etc.
	  accum_init = 0.
	end if

	naccu = 0.
	std = 0.
	accum = accum_init

	end subroutine shyutil_init_accum

!***************************************************************

	subroutine shyutil_init_dir(nlv,nn,nvar,istep)

	integer nlv,nn,nvar,istep

        allocate(dir(nlv,nn,idir,0:nvar,istep))

	dir = 0.

	end subroutine shyutil_init_dir

!==================================================================
	end module shyutil
!==================================================================

	subroutine shy_make_vert_aver(idims,nndim,cv3,cv2)

	use basin
	use levels
	use shyutil

	implicit none

	integer idims(4)
	integer nndim
	real cv3(nlvdi,nndim)
	real cv2(nndim)

	integer ivar,lmax,nn,ie,i
	real cmin,cmax,cmed,cstd,atot,vtot
	integer iflag(nndim)

        ivar = idims(4)
        lmax = idims(3)
        nn = idims(1) * idims(2)

	if( lmax == 1 ) then
	  cv2(:) = cv3(1,:)
	else if( ivar == 3 ) then	!transports - must sum
	  do ie=1,nel
	    cv2(ie) = sum( cv3(1:ilhv(ie),ie) )
	  end do
	else
	  iflag = 1
	  if( nn == nkn ) then
	    call make_aver_3d(nlvdi,nn,cv3,areak,vol3k,ilhkv,iflag
     +				,cmin,cmax,cmed,cstd,atot,vtot,cv2)
	  else if( nn == nel ) then
	    call make_aver_3d(nlvdi,nn,cv3,areae,vol3e,ilhv,iflag
     +				,cmin,cmax,cmed,cstd,atot,vtot,cv2)
	  else
	    write(6,*) ivar,nn,nkn,nel
	    stop 'error stop shy_make_vert_aver: not possible'
	  end if
	end if

	end

!***************************************************************

	subroutine shy_make_basin_aver(idims,nlvddi,nndim,cv3,iflag
     +				,cmin,cmax,cmed,cstd,atot,vtot)

	use basin
	use levels
	use shyutil

	implicit none

	integer idims(4)
	integer nlvddi
	integer nndim
	real cv3(nlvddi,nndim)
	integer iflag(nndim)
	real cmin,cmax,cmed,cstd,atot,vtot

	integer ivar,lmax,nn
	real ze(3*nel)
	real zaux(nel)
	real cv2(nndim)

        ivar = idims(4)
        lmax = idims(3)
        nn = idims(1) * idims(2)

	if( abs(ivar) == 1 ) then		! water level - 2D
	  if( lmax /= 1 .or. nlvddi /= 1 ) then
	    write(6,*) ivar,nn,lmax,nkn,nel,nlvddi
	    stop 'error stop shy_make_basin_aver: level must have lmax=1'
	  end if
	  if( nn == nkn ) then
	    call make_aver_2d(nlvdi,nn,cv3,areak,iflag
     +				,cmin,cmax,cmed,cstd,atot,vtot)
	    vtot = sum(vol2k,mask=iflag>0)
	  else if( nn == nel ) then
	    call make_aver_2d(nlvdi,nn,cv3,areae,iflag
     +				,cmin,cmax,cmed,cstd,atot,vtot)
	    vtot = sum(vol2e,mask=iflag>0)
	  else if( nn == 3*nel ) then			!zenv
	    ze(:) = cv3(1,1)
	    call shy_make_zeta_from_elem(nel,ze,zaux)
	    call make_aver_2d(nlvdi,nel,zaux,areae,iflag
     +				,cmin,cmax,cmed,cstd,atot,vtot)
	    vtot = sum(vol2e,mask=iflag>0)
	  else
	    write(6,*) ivar,nn,nkn,nel
	    stop 'error stop shy_make_basin_aver: not possible'
	  end if
	else
	  if( nlvddi /= nlvdi ) then
	    write(6,*) ivar,nn,lmax,nkn,nel,nlvddi,nlvdi
	    stop 'error stop shy_make_basin_aver: nlvdi/=nlvddi'
	  end if
	  if( nn == nkn ) then
	    call make_aver_3d(nlvddi,nn,cv3,areak,vol3k,ilhkv,iflag
     +				,cmin,cmax,cmed,cstd,atot,vtot,cv2)
	  else if( nn == nel ) then
	    call make_aver_3d(nlvddi,nn,cv3,areae,vol3e,ilhv,iflag
     +				,cmin,cmax,cmed,cstd,atot,vtot,cv2)
	  else
	    write(6,*) ivar,nn,nkn,nel
	    stop 'error stop shy_make_basin_aver: not possible'
	  end if
	end if

	end

!***************************************************************

	subroutine make_aver_2d(nlvddi,nn,cv3,area,iflag
     +				,cmin,cmax,cmed,cstd,atot,vtot)

	implicit none

	integer nlvddi
	integer nn
	real cv3(nlvddi,nn)
	real area(nn)
	integer iflag(nn)
	real cmin,cmax,cmed,cstd,atot,vtot

	integer i
	double precision c,v
	double precision cctot,vvtot,c2tot,ccmin,ccmax
	real, parameter :: high = 1.e+30

	ccmin = +high
	ccmax = -high
	cctot = 0.
	vvtot = 0.
	c2tot = 0.

	do i=1,nn
	    if( iflag(i) <= 0 ) cycle
	    c = cv3(1,i)
	    v = area(i)
	    ccmin = min(ccmin,c)
	    ccmax = max(ccmax,c)
	    cctot = cctot + v*c
	    c2tot = c2tot + v*c*c
	    vvtot = vvtot + v
	end do

	cmin = ccmin
	cmax = ccmax
	vtot = vvtot
	atot = vvtot
	if( vtot == 0 ) vvtot = 1.	!avoid Nan
	cmed = cctot / vvtot
	cstd = c2tot / vvtot - cmed**2
	if( cstd < 0 ) cstd = 0.
	cstd = sqrt( cstd )

	end

!***************************************************************

	subroutine make_aver_3d(nlvddi,nn,cv3,area,vol,il,iflag
     +				,cmin,cmax,cmed,cstd,atot,vtot,cv2)

	implicit none

	integer nlvddi
	integer nn
	real cv3(nlvddi,nn)
	real area(nn)
	real vol(nlvddi,nn)
	integer il(nn)
	integer iflag(nn)
	real cmin,cmax,cmed,cstd,atot,vtot
	real cv2(nn)

	integer i,l,lmax
	double precision c,v,a
	double precision cctot,vvtot,c2tot,aatot,ccmin,ccmax
	double precision c2,v2
	real, parameter :: high = 1.e+30
        integer :: ks = 0
        logical bdebug
        logical bcompute

	ccmin = +high
	ccmax = -high
	cctot = 0.
	c2tot = 0.
	vvtot = 0.
	aatot = 0.

	do i=1,nn
	  bcompute = ( iflag(i) > 0 )
	  lmax = il(i)
	  bdebug = i == ks
	  c2 = 0.
	  v2 = 0.
	  a = area(i)
	  do l=1,lmax
	    c = cv3(l,i)
	    v = vol(l,i)
	    c2 = c2 + v*c
	    v2 = v2 + v
	    if( bcompute ) then
	      ccmin = min(ccmin,c)
	      ccmax = max(ccmax,c)
	      cctot = cctot + v*c
	      c2tot = c2tot + v*c*c
	      vvtot = vvtot + v
	    end if
	    if( bdebug ) write(41,*) l,v,c
	  end do
	  cv2(i) = c2 / v2
	  if( bcompute ) aatot = aatot + a
	end do

	cmin = ccmin
	cmax = ccmax
	vtot = vvtot
	atot = aatot
	if( vtot == 0 ) vvtot = 1.	!avoid Nan
	cmed = cctot / vvtot
	cstd = c2tot / vvtot - cmed**2
	if( cstd < 0 ) cstd = 0.
	cstd = sqrt( cstd )

	end

!***************************************************************
!***************************************************************
!***************************************************************

	subroutine shy_make_zeta(ftype)

	use basin
	use shyutil

	implicit none

	integer ftype

	if( ftype == 1 ) then
	  call shy_make_zeta_from_elem(nel,shy_zenv,shy_zeta)
	else if( ftype == 2 ) then
	  call shy_make_zeta_from_node(nkn,nel,nen3v,shy_znv,shy_zeta)
	else
	  shy_zeta = 0.
	end if

	end

!***************************************************************

	subroutine shy_make_zeta_from_node(nkn,nel,nen3v,znv,zeta)

	implicit none

	integer nkn,nel
	integer nen3v(3,nel)
	real znv(nkn)
	real zeta(nel)

	integer ie,k,ii
	real z

	do ie=1,nel
	  z = 0.
	  do ii=1,3
	    k = nen3v(ii,ie)
	    z = z + znv(k)
	  end do
	  zeta(ie) = z / 3.
	end do

	end

!***************************************************************

	subroutine shy_make_zeta_from_elem(nel,zenv,zeta)

	implicit none

	integer nel
	real zenv(3,nel)
	real zeta(nel)

	integer ie,ii
	real z

	do ie=1,nel
	  z = 0.
	  do ii=1,3
	    z = z + zenv(ii,ie)
	  end do
	  zeta(ie) = z / 3.
	end do

	end

!***************************************************************
!***************************************************************
!***************************************************************

	subroutine shy_make_area

! computes area of elements and nodes

	use basin
	use evgeom
	use shyutil

	implicit none

	logical binit
	integer ie,ii,k
	real area

	call is_init_ev(binit)
	if( .not. binit ) then
	  stop 'error stop shy_make_area: ev not initialized'
	end if

	areak = 0.

	do ie=1,nel
	  area = 4. * ev(10,ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    areak(k) = areak(k) + area
	  end do
	  areae(ie) = 3. * area
	end do

	end

!***************************************************************

	subroutine shy_make_volume

! initializes volumes

	use basin
	use shyutil
	use levels
	use mod_depth

	implicit none

	logical bvolwrite,bdebug
	integer ie,ii,k,l,lmax,nsigma,nlvaux,ks
	real z,h,hsigma,zeps
	double precision ak,vk,ve
	double precision, allocatable :: vole(:,:),volk(:,:)
	real hl(nlv)		!aux vector for layer thickness

	zeps = 0.01

	bvolwrite = .true.
	bvolwrite = .false.
	ks = 2985
	ks = 3096
	ks = 0

        call get_sigma_info(nlvaux,nsigma,hsigma)

	vol3e = 0.
	vol3k = 0.
	allocate(vole(nlvdi,nel))
	allocate(volk(nlvdi,nkn))
	vole = 0.
	volk = 0.

	do ie=1,nel
	  ak = areae(ie) / 3.	!area of vertex
	  h = hev(ie)
	  z = shy_zeta(ie)
	  if( h+z < zeps ) z = zeps - h	!make volume positive
	  lmax = ilhv(ie)
	  !write(6,*) ie,lmax,nlv,nlvdi
	  call get_layer_thickness(lmax,nsigma,hsigma,z,h,hlv,hl)
	  do l=1,lmax
	    vk = ak * hl(l)
	    ve = 3. * vk
	    vole(l,ie) = vole(l,ie) + ve
	    do ii=1,3
	      k = nen3v(ii,ie)
	      bdebug = k == ks
	      volk(l,k) = volk(l,k) + vk
	      if( bdebug ) write(81,*) ie,ii,lmax,l,vk,ak,hl(l),h
	    end do
	  end do
	  if( bvolwrite .and. mod(ie,nel/nel) == -1 ) then
	    write(61,*) ie,h
	    write(61,*) hl
	  end if
	end do

	vol3e = vole
	vol3k = volk
	vol2e = sum(vol3e,dim=1)
	vol2k = sum(vol3k,dim=1)

	deallocate(vole,volk)

	if( nlv <= 1 ) bvolwrite = .false.
	if( bvolwrite ) then
	write(71,*) 'volume from shy_make_volume'
	do k=1,nkn,nkn/10
	  write(71,*) k
	  write(71,*) vol3k(:,k)
	end do
	end if

	end

!***************************************************************
!***************************************************************
!***************************************************************

	subroutine shy_time_aver(bforce,avermode,iv,nread,ifreq,istep
     +				,nndim,idims,threshold,cv3,bout,bverb)

! avermode:  1:aver  2:sum  3:min  4:max  5:std  6:rms  7:thres  8:averdir
!
! avermode negative: only transform, do not accumulate

	use basin
	use shyutil
	use levels

	implicit none

	logical bforce
	integer avermode
	integer iv
	integer nread,ifreq,istep
	integer nndim
	integer idims(4)
	double precision threshold
	real cv3(nlvdi,nndim)	!input for accumulation, output if transformed
	logical bout		!.true. if in cv3 are transformed results
	logical bverb

	integer ip,naccum,mmode
	integer k,l,id,idmax
	integer ivar,lmax,nn
	double precision dmax

!---------------------------------------------------------------
! see if we have to do something
!---------------------------------------------------------------

	if( avermode .eq. 0 ) return

!---------------------------------------------------------------
! prepare for accumulation
!---------------------------------------------------------------

	! either istep == 1 or ifreq == 0

	if( istep < 1 .or. istep > 1 .and. ifreq > 0 ) then
	  write(6,*) 'istep,ifreq: ',istep,ifreq
	  stop 'error stop shy_time_aver: error in parameters'
	end if

	ip = mod(nread,istep)
	if( ip .eq. 0 ) ip = istep

	!write(6,*) 'ip: ',ip,istep,nread,avermode
	!write(6,*) '... ',ifreq,avermode,ip,istep,naccu(iv,ip)

        ivar = idims(4)
        lmax = idims(3)
        nn = idims(1) * idims(2)

!---------------------------------------------------------------
! accumulate
!---------------------------------------------------------------

	if( avermode > 0 ) naccu(iv,ip) = naccu(iv,ip) + 1

	if( avermode == 1 .or. avermode == 2 ) then
	  accum(:,:,iv,ip) = accum(:,:,iv,ip) + cv3(:,:)
	else if( avermode == 3 ) then
	  !accum(:,:,iv,ip) = min(accum(:,:,iv,ip),cv3(:,:))
	  where ( cv3(:,:) < accum(:,:,iv,ip) ) 
	    accum(:,:,iv,ip) = cv3(:,:)
	  end where
	else if( avermode == 4 ) then
	  !accum(:,:,iv,ip) = max(accum(:,:,iv,ip),cv3(:,:))
	  where ( cv3(:,:) > accum(:,:,iv,ip) )
	    accum(:,:,iv,ip) = cv3(:,:)
	  end where
	else if( avermode == 5 ) then
	  accum(:,:,iv,ip) = accum(:,:,iv,ip) + cv3(:,:)
	  std(:,:,iv,ip) = std(:,:,iv,ip) + cv3(:,:)**2
	else if( avermode == 6 ) then
	  accum(:,:,iv,ip) = accum(:,:,iv,ip) + cv3(:,:)**2
	else if( avermode == 7 ) then
	  where( cv3(:,:) >= threshold )
	    accum(:,:,iv,ip) = accum(:,:,iv,ip) + 1.
	  end where
	else if( avermode == 8 ) then
	  call shy_accum_dir(nlvdi,lmax,nn,iv,ip,cv3)
	end if

!---------------------------------------------------------------
! transform (average)
!---------------------------------------------------------------

	bout = ( naccu(iv,ip) == ifreq .or. avermode < 0 .or. bforce )

	if( .not. bout ) return

	naccum = max(1,naccu(iv,ip))
	mmode = abs(avermode)
	if( mmode == 2 ) naccum = 1			!sum
	if( mmode == 3 ) naccum = 1			!min
	if( mmode == 4 ) naccum = 1			!max
	if( mmode == 7 ) naccum = 1			!threshold
	if( naccum > 0 ) cv3(:,:) = accum(:,:,iv,ip)/naccum
	if( mmode == 5 ) then
	  !cv3(:,:) = sqrt( std(:,:,iv,ip)/naccum - cv3(:,:)**2 )
	  cv3(:,:) = std(:,:,iv,ip)/naccum - cv3(:,:)**2
	  where( cv3 > 0 )
	    cv3 = sqrt( cv3 )
	  else where
	    cv3 = 0.
	  end where
	else if( mmode == 6 ) then
	  cv3(:,:) = sqrt( cv3(:,:) )
	else if( mmode == 8 ) then
	  call shy_elab_dir(nlvdi,lmax,nn,iv,ip,cv3)
	end if

	if( bverb ) then
	  write(6,*) 'averaging: ',ip,naccum,naccu(iv,ip)
	end if

	naccu(iv,ip) = 0
	std(:,:,iv,ip) = 0.
	accum(:,:,iv,ip) = accum_init

!---------------------------------------------------------------
! end of routine
!---------------------------------------------------------------

	end

!***************************************************************

	subroutine shy_accum_dir(nlvddi,lmax,nn,iv,ip,cv3)

	use shyutil

	implicit none

	integer nlvddi
	integer nn,lmax
	integer iv,ip
	real cv3(nlvddi,nn)

	integer k,l,id

	do k=1,nn
	  do l=1,lmax
	    id = nint( cv3(l,k)/adir )
	    if( id == 0 ) id = idir
	    if( id < 0 .or. id > idir ) stop 'error stop: direction'
	    dir(l,k,id,iv,ip) = dir(l,k,id,iv,ip) + 1.
	  end do
	end do

	end

!***************************************************************

	subroutine shy_elab_dir(nlvddi,lmax,nn,iv,ip,cv3)

	use shyutil

	implicit none

	integer nlvddi
	integer nn,lmax
	integer iv,ip
	real cv3(nlvddi,nn)

	integer k,l,id,idmax
	double precision dmax

	do k=1,nn
	  do l=1,lmax
	    dmax = 0.
	    idmax = 0
	    do id=1,idir
	      if( dir(l,k,id,iv,ip) > dmax ) then
	        idmax = id
	        dmax = dir(l,k,id,iv,ip)
	      end if
	    end do
		if( idmax == idir ) idmax = 0
		cv3(l,k) = idmax * adir
	  end do
	end do

	dir(:,:,:,iv,ip) = 0.

	end

!***************************************************************

