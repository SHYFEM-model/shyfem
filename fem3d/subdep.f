
!--------------------------------------------------------------------------
!
!    Copyright (C) 1994,1997-1999,2002,2005-2006,2010-2011  Georg Umgiesser
!    Copyright (C) 2013-2019  Georg Umgiesser
!    Copyright (C) 2006  Christian Ferrarin
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

c depth utility routines
c
c contents :
c
c function igtdep(k,f)                  get depth for node
c function igtdpa(mode,h)               gets unique depth for all nodes
c subroutine huniqu(hev,hkv)            make depth unique for every node
c subroutine makehev(hev)		makes hev (elementwise depth)
c subroutine makehkv(hkv)		makes hkv (nodewise depth)
c subroutine depadj(hmin,hmax,href)	adjusts depth to ref/min/max values
c
c revision log :
c
c 02.02.1994	ggu	$$nmax - check error condtion nmax
c 29.06.1997	ggu	$$ndim - dimension of f is passed
c 29.06.1997	ggu	depth routines in one file
c 06.11.1998	ggu	new huniqu to compute hev and hkv
c 19.10.1999	ggu	new routine makehv from subutl
c 25.03.2002	ggu	new routines makehkv (before makehv) and makehev
c 28.11.2005	ggu	makehkv changed (uses real aux value, area weight)
c 24.02.2006	ggu	bug in makehkv -> haux was integer
c 18.10.2006	ccf	bug in makehkv -> no area multiplication
c 23.03.2010	ggu	changed v6.1.1
c 16.12.2010	ggu	in depadj() do not set hm3v to constant
c 17.05.2011	ggu	new routines to adjourn depth
c 31.05.2011	ggu	changed VERS_6_1_23
c 18.10.2011	ggu	changed VERS_6_1_33
c 18.11.2011	ggu	new routine makehkv_minmax()
c 05.09.2013	ggu	new routine set_sigma_hkv_and_hev() from newsig.f
c 12.09.2013	ggu	changed VERS_6_1_67
c 25.06.2014	ggu	computa also hkv_min and hkv_max
c 12.12.2014	ggu	changed VERS_7_0_9
c 23.12.2014	ggu	changed VERS_7_0_11
c 19.01.2015	ggu	changed VERS_7_1_3
c 05.05.2015	ggu	changed VERS_7_1_10
c 21.05.2015	ggu	changed VERS_7_1_11
c 25.05.2015	ggu	some changes in depth computation
c 05.06.2015	ggu	changed VERS_7_1_12
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 18.12.2015	ggu	changed VERS_7_3_17
c 25.05.2016	ggu	changed VERS_7_5_10
c 05.12.2017	ggu	changed VERS_7_5_39
c 19.04.2018	ggu	changed VERS_7_5_45
c 16.02.2019	ggu	changed VERS_7_5_60
c 13.03.2019	ggu	changed VERS_7_5_61
c 21.05.2019	ggu	changed VERS_7_5_62
c 05.12.2022	ggu	use ie_mpi for accumulation and computation of hkv
c
c********************************************************************

	function igtdep(k,f,ndim)

c gets depth given a node number
c if different values of depth are associated
c with a node, it returns all these values
c
c k             node number (internal)
c f             vector in which the depth values are
c               ...stored at return
c ndim		dimension of f
c igtdep        number of different values found

	use basin
	use shympi

	implicit none

	integer igtdep
	integer k,ndim
	real f(ndim)

	integer iact,ie,i,ii,ie_mpi

	iact=0
        do ie_mpi=1,nel
	  ie = ip_sort_elem(ie_mpi)
	  do ii=1,3
	    if(nen3v(ii,ie).eq.k) then
		do i=1,iact
		    if(f(i).eq.hm3v(ii,ie)) goto 1
		end do

		iact=iact+1     		!new depth
		if(iact.gt.ndim) goto 99	!$$nmax !$$ndim
		f(iact)=hm3v(ii,ie)

    1           continue       			!old depth
	    end if
	  end do
	end do

	igtdep=iact

	return
   99	continue
	stop 'error stop igtdep : nmax'		!$$nmax
	end

c********************************************************************

	function igtdpa(mode,h)

c gets unique depth for all nodes
c
c mode          switch
c               1       deepest value is returned
c               -1      most shallow value
c h             vector in which the depth values are
c               ...stored (return value)
c igtdpa        return status
c               1       unique depth
c               0       no unique depth
c               -1      error

	use basin

	implicit none

	integer igtdpa
	integer mode
	real h(nkn)

	real high
	parameter(high=1.e+30)

	logical buniq
	integer ie,ii,i,k
	real hh,hhh,hflag

	buniq=.true.

	if(mode.eq.1) then
		hflag=-high
	else if(mode.eq.-1) then
		hflag=high
	else
		write(6,*) 'Value for mode not allowed :',mode
		igtdpa=-1
		return
	end if

	do i=1,nkn
	   h(i)=hflag
	end do

	do ie=1,nel
	 do ii=1,3
	   k=nen3v(ii,ie)
	   hh=hm3v(ii,ie)
	   hhh=h(k)
	   if(mode.eq.1) then
		if(hh.gt.h(k)) h(k)=hh
	   else
		if(hh.lt.h(k)) h(k)=hh
	   end if
	   if(hhh.ne.hflag.and.hhh.ne.hh) buniq=.false.
	 end do
	end do

	do i=1,nkn
	   if(h(i).eq.hflag) then
		write(6,*) 'igtdpa : Nodes without depth'
		igtdpa=-1
		return
	   end if
	end do

	if(buniq) then
		igtdpa=1
	else
		igtdpa=0
	end if

	end

c********************************************************************

	subroutine huniqu(hev,hkv)

c make depth unique for every node (changes hm3v)
c nodal values are the highest (deepest) value
c
c hev		element averaged depth values
c hkv            array with unique depth values

	use basin

	implicit none

	real hev(nel)
	real hkv(nkn)

	integer ie,ii,k
	logical bstop
	real h,flag,hm

	integer ipext

c flag nodal values

	flag = -999.

	do k=1,nkn
	  hkv(k) = flag
	end do

c create element averaged depth values and assign to nodal values

	do ie=1,nel
	  hm = 0.
	  do ii=1,3
	    hm = hm + hm3v(ii,ie)
	  end do

	  hev(ie) = hm / 3.

	  do ii=1,3
	    hm3v(ii,ie) = hev(ie)
	  end do

	  h = hev(ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    if( k .le. 0 ) write(6,*) 'huniqu: ',ie,ii,k,hev(ie)
	    if( h .gt. hkv(k) ) hkv(k) = h
	  end do
	end do

c check if all depth values are available

	bstop = .false.

	do k=1,nkn
	  if( hkv(k) .eq. flag ) then
		write(6,*) 'No depth for node ',ipext(k)
		bstop = .true.
	  end if
	end do

	if( bstop ) stop 'error stop huniqu'

	end

c********************************************************************
c********************************************************************
c********************************************************************

        subroutine makehev(hev)

c makes hev (elementwise depth)

	use basin

        implicit none

c arguments
        real hev(nel)
c local
        integer ie,ii
	double precision hm

        do ie=1,nel
	  hm = 0.
          do ii=1,3
	    hm = hm + hm3v(ii,ie)
          end do
	  hev(ie) = hm / 3.
        end do

        end

c********************************************************************

        subroutine makehkv(hkv)

c makes hkv (nodewise depth)

	use basin
	use shympi

        implicit none

c arguments
        real hkv(nkn)
c local
        integer ie,ii,k,kn,ie_mpi
	real weight
        real haux(nkn)   !aux array -> bug - was integer

	real weight_elem

        do k=1,nkn
          hkv(k) = 0.
          haux(k) = 0.
        end do

        do ie_mpi=1,nel
	  ie = ip_sort_elem(ie_mpi)
	  weight = weight_elem(ie)
          do ii=1,3
            kn=nen3v(ii,ie)
            hkv(kn)=hkv(kn)+hm3v(ii,ie)*weight	!ccf
            haux(kn)=haux(kn)+weight
          end do
        end do

	if( shympi_partition_on_elements() ) then
          !call shympi_comment('shympi_elem: exchange hkv, haux')
          call shympi_exchange_and_sum_2d_nodes(hkv)
          call shympi_exchange_and_sum_2d_nodes(haux)
	end if

        hkv = hkv / haux

	if( shympi_partition_on_nodes() ) then
          call shympi_exchange_2d_node(hkv)
	end if

        end

c********************************************************************

        subroutine makehkv_minmax(hkv,itype)

c makes hkv (nodewise depth)
c
c itype:  -1: min  0: aver  +1: max

	use basin
	use shympi

        implicit none

        real hkv(nkn)
        integer itype

        integer k,ie,ii
        real hinit,h

c-------------------------------------------------------
c initialize
c-------------------------------------------------------

	if( itype .eq. 0 ) then
          call makehkv(hkv)
	  return
	end if

        if( itype .lt. 0. ) then
          hinit = +1.e+30
        else
          hinit = -1.e+30
        end if

        do k=1,nkn
          hkv(k) = hinit
        end do

c-------------------------------------------------------
c set hkv
c-------------------------------------------------------

        do ie=1,nel
          do ii=1,3
            h = hm3v(ii,ie)
            k = nen3v(ii,ie)
            if( itype .lt. 0 ) then
              hkv(k) = min(hkv(k),h)
            else
              hkv(k) = max(hkv(k),h)
            end if
          end do
        end do

	if( shympi_partition_on_elements() ) then
          if( itype .lt. 0 ) then
            call shympi_exchange_2d_nodes_min(hkv)
            !call shympi_comment('shympi_elem: exchange hkv_min')
          else
            call shympi_exchange_2d_nodes_max(hkv)
            !call shympi_comment('shympi_elem: exchange hkv_max')
          end if
	else
          call shympi_exchange_2d_node(hkv)
	end if

c-------------------------------------------------------
c end of routine
c-------------------------------------------------------

        end

c********************************************************************
c********************************************************************
c********************************************************************

	subroutine depadj(hmin,hmax,href)

c adjusts depth to reference and min/max values - only hm3v is changed

	use basin

	implicit none

	real hmin,hmax,href

	integer iaux,ie,ii
	real hmed

c adjust depth to constant in element %%%%%%%%%%%%%%%%%%%%%%


c adjust depth to minimum depth %%%%%%%%%%%%%%%%%%%%%%%%%%%%

	iaux=0
	do ie=1,nel
	 do ii=1,3
	  if(hm3v(ii,ie).lt.hmin) then
	    hm3v(ii,ie)=hmin
	    iaux=iaux+1
	  end if
	 end do
	end do

	if(iaux.gt.0) then
	  write(6,*) '***********************************************'
	  write(6,*) '***********************************************'
	  write(6,*) '***********************************************'
	  write(6,*) '*********** hmin = ',     hmin  ,' ************'
	  write(6,*) '*********** changes = ',  iaux  ,' ************'
	  write(6,*) '***********************************************'
	  write(6,*) '***********************************************'
	  write(6,*) '***********************************************'
	end if

c adjust depth to maximum depth %%%%%%%%%%%%%%%%%%%%%%%%%%%%

	iaux=0
	do ie=1,nel
	 do ii=1,3
	  if(hm3v(ii,ie).gt.hmax) then
	    hm3v(ii,ie)=hmax
	    iaux=iaux+1
	  end if
	 end do
	end do

	if(iaux.gt.0) then
	  write(6,*) '***********************************************'
	  write(6,*) '***********************************************'
	  write(6,*) '***********************************************'
	  write(6,*) '*********** hmax = ',     hmax  ,' ************'
	  write(6,*) '*********** changes = ',  iaux  ,' ************'
	  write(6,*) '***********************************************'
	  write(6,*) '***********************************************'
	  write(6,*) '***********************************************'
	end if

c adjust depth to reference level %%%%%%%%%%%%%%%%%%%%%%%%%%%%

	do ie=1,nel
	 do ii=1,3
	  hm3v(ii,ie)=hm3v(ii,ie)-href
	 end do
	end do

	end

c********************************************************************

	subroutine adjust_depth

c adjusts depth values - only hm3v is changed

	implicit none

	real hmin,hmax,href
	real getpar

c       call bocche     !FIXME

        hmin=getpar('hmin')
        hmax=getpar('hmax')
        href=getpar('href')

        call depadj(hmin,hmax,href)	!adjusts h=h-href and hmax<h<hmin

	end

c********************************************************************
c********************************************************************
c********************************************************************

	subroutine get_hmax_global(hmax)

	use basin
	use shympi

	implicit none

	real hmax
	real h

	h = maxval(hm3v)
	hmax = shympi_max(h)

	write(6,*) 'hmax: ',my_id,h,hmax

	end

c********************************************************************
c********************************************************************
c********************************************************************

	subroutine set_depth

c sets up depth arrays

	implicit none

	logical bsigma
	integer nlv,nsigma
	real hsigma

	call get_sigma_info(nlv,nsigma,hsigma)
	bsigma = nsigma .gt. 0

	if( bsigma ) then	!sigma or hybrid layers
	  call check_sigma_hsigma
	  call flatten_hm3v(hsigma)
	else
	  call flatten_hm3v(-999.)
	end if

	call make_hev
	call make_hkv

	call exchange_depth

	end

c********************************************************************

	subroutine exchange_depth

c exchanges nodal depth values between domains

	use mod_depth
	use shympi
	!use basin
	!use mod_random

	implicit none

	!call shympi_comment('exchanging hkv, hdk_min, hdk_max')
	call shympi_exchange_2d_node(hkv)
	call shympi_exchange_2d_node(hkv_min)
	call shympi_exchange_2d_node(hkv_max)
	!call shympi_barrier

	call shympi_check_2d_elem(hev,'hev')
	call shympi_check_2d_node(hkv,'hkv')
	call shympi_check_2d_node(hkv_min,'hkv_min')
	call shympi_check_2d_node(hkv_max,'hkv_max')

	end


c********************************************************************

	subroutine flatten_hm3v(hsigma)

c flattens bottom below hsigma - only hsigma is used, hm3v is changed

	use basin

	implicit none

	real hsigma

	integer ie,ii
	real hm

	do ie=1,nel
	  hm = 0.
	  do ii=1,3
	    hm = hm + hm3v(ii,ie)
	  end do
	  hm = hm / 3.
	  if( hm .gt. hsigma ) then
	    do ii=1,3
	      hm3v(ii,ie) = hm
	    end do
	  end if
	end do

	end

c********************************************************************

	subroutine make_hkv

c adjusts nodal depth values

	use mod_depth
	use basin, only : nkn,nel,ngr,mbw

	implicit none

        call makehkv(hkv)		!computes hkv as average
        call makehkv_minmax(hkv_min,-1)
        call makehkv_minmax(hkv_max,+1)

	end

c********************************************************************

	subroutine make_hev

c adjusts elemental depth values

	use mod_depth

	implicit none

        call makehev(hev)

	end

c********************************************************************
c********************************************************************
c********************************************************************

	subroutine adjourne_depth_from_hm3v

c adjourns hev and hkv from hm3v (if it has been changed)

	use mod_depth
	use shympi

	implicit none

        call make_hev
        call make_hkv
        !call set_last_layer		!FIXME

	call shympi_exchange_2d_node(hkv)

	end

c********************************************************************

	subroutine adjourn_depth_from_hev

c adjourns hev and hkv from hm3v (if it has been changed)

	use mod_depth
	use basin

	implicit none

	integer ie,ii

	do ie=1,nel
	  do ii=1,3
	    hm3v(ii,ie) = hev(ie)
	  end do
	end do

        call make_hkv
        !call set_last_layer		!FIXME

	end

c********************************************************************
c********************************************************************
c********************************************************************

	subroutine check_sigma_hsigma

c checks hkv and hsigma
c uses information about sigma layers and hsigma (hybrid)

	use basin

	implicit none

	logical berror
	integer k,ie,ii
	integer inc,ihmin,ihmax
	integer nlv,nsigma
	real flag
	real hm,h
	real hsigma
	real hkv(nkn)		!local

c-------------------------------------------------------
c initialize
c-------------------------------------------------------

	call get_sigma_info(nlv,nsigma,hsigma)

	if( nsigma == 0 ) return

	flag = -999.

	do k=1,nkn
	  hkv(k) = flag
	end do

	berror = .false.
	inc = 0

c-------------------------------------------------------
c set if hkv is continuous in sigma layers
c-------------------------------------------------------

	do ie=1,nel

	  hm = 0.
	  do ii=1,3
	    h = hm3v(ii,ie)
	    hm = hm + h
	  end do
	  hm = hm / 3.

	  if( hm > hsigma ) cycle

	  do ii=1,3
	    h = hm3v(ii,ie)
	    k = nen3v(ii,ie)
	    if( hkv(k) .eq. flag ) then
	      hkv(k) = h
	    else
	      if( h .ne. hkv(k) ) then
		write(6,*) 'depth of node not unique: ',ie,k,h,hkv(k)
	        inc = inc + 1
	      end if
	    end if
	  end do

	end do

	if( inc .gt. 0 ) then
	  write(6,*) 'number of occurences found: ',inc
	  stop 'error stop set_hkv_and_hev: depth not unique'
	end if

c-------------------------------------------------------
c check hsigma crossing
c-------------------------------------------------------

        do ie=1,nel
          ihmin = 0
          ihmax = 0
          do ii=1,3
            h = hm3v(ii,ie)
            if( h .lt. hsigma ) then
              ihmin = ihmin + 1
            else if( h .gt. hsigma ) then
              ihmax = ihmax + 1
            end if
          end do
          if( ihmin .gt. 0 .and. ihmax .gt. 0 ) then
	    write(6,*) 'hsigma crossing: ',ie,(hm3v(ii,ie),ii=1,3)
	    berror = .true.
	  end if
	end do

	if( berror ) then
	  write(6,*) 'elements with hsigma crossing depths'
	  stop 'error stop set_hkv_and_hev: hsigma crossing'
	end if

c-------------------------------------------------------
c end of routine
c-------------------------------------------------------

	end

c********************************************************************
c********************************************************************
c********************************************************************

	subroutine read_in_hev(file)

	use mod_depth
	use basin, only : nkn,nel,ngr,mbw

	character*(*) file

	integer ie,ios
	logical, save :: berror = .true.	!throw error if not found

	open(1,file=file,status='old',form='formatted',iostat=ios)

	if( ios /= 0 ) then
	  write(6,*) '*** cannot open file: ',trim(file)
	  write(6,*) '...not initializing hev'
	  if( berror ) stop 'error stop read_in_hev'
	  return
	end if

	read(1,*) nelaux
	if( nel .ne. nelaux ) stop 'error stop read_in_hev: nel'
	read(1,*) (hev(ie),ie=1,nel)
	close(1)

	write(6,*) '======================================'
	write(6,*) 'hev data read from file: ',file
	write(6,*) '======================================'

	end

c********************************************************************

	subroutine write_out_hev(file)

	use mod_depth
	use basin, only : nkn,nel,ngr,mbw

	character*(*) file

	integer ie

	open(1,file=file,status='unknown',form='formatted')
	write(1,*) nel
	write(1,*) (hev(ie),ie=1,nel)
	close(1)

	end

c********************************************************************

