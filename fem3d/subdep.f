c
c $Id: subdep.f,v 1.10 2008-07-16 15:41:39 georg Exp $
c
c depth utility routines
c
c contents :
c
c function igtdep(k,f)                  get depth for node
c function igtdpa(mode,h)               gets unique depth for all nodes
c subroutine huniqu(hev,hkv)            make depth unique for every node
c subroutine makehev(hev)		makes hev (elementwise depth)
c subroutine makehkv(hkv,haux)		makes hkv (nodewise depth)
c subroutine depadj(hmin,hmax,href)	adjusts depth to ref/min/max values
c
c revision log :
c
c 29.06.1997	ggu	depth routines in one file
c 06.11.1998	ggu	new huniqu to compute hev and hkv
c 19.10.1999	ggu	new routine makehv from subutl
c 25.03.2002	ggu	new routines makehkv (before makehv) and makehev
c 28.11.2005	ggu	makehkv changed (uses real aux value, area weight)
c 24.02.2006	ggu	bug in makehkv -> haux was integer
c 18.10.2006	ccf	bug in makehkv -> no area multiplication
c 16.12.2010	ggu	in depadj() do not set hm3v to constant
c 17.05.2011	ggu	new routines to adjourn depth
c 18.11.2011	ggu	new routine makehkv_minmax()
c 05.09.2013	ggu	new routine set_sigma_hkv_and_hev() from newsig.f
c 25.06.2014	ggu	computa also hkv_min and hkv_max
c 25.05.2015	ggu	some changes in depth computation
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
c
c revised 02.02.94 by ggu	$$nmax - check error condtion nmax
c revised 29.06.97 by ggu	$$ndim - dimension of f is passed

	use basin

	implicit none

	integer igtdep
	integer k,ndim
	real f(1)

	include 'param.h'

	integer iact,ie,i,ii

	iact=0
	do ie=1,nel
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
	real h(1)

	real high
	parameter(high=1.e+30)

	include 'param.h'

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

	real hev(1)
	real hkv(1)


	include 'param.h'

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

        subroutine makehev(hev)

c makes hev (elementwise depth)

	use basin

        implicit none

c arguments
        real hev(1)
c common
	include 'param.h'
c local
        integer ie,ii
	real hm

        do ie=1,nel
	  hm = 0.
          do ii=1,3
	    hm = hm + hm3v(ii,ie)
          end do
	  hev(ie) = hm / 3.
        end do

        end

c********************************************************************

        subroutine makehkv(hkv,haux)

c makes hkv (nodewise depth)

	use basin

        implicit none

c arguments
        real hkv(1)
        real haux(1)   !aux array -> bug - was integer
c common
	include 'param.h'
c local
        integer ie,ii,k,kn
	real weight

	real weight_elem

        do k=1,nkn
          hkv(k) = 0.
          haux(k) = 0.
        end do

        do ie=1,nel
	  weight = weight_elem(ie)
          do ii=1,3
            kn=nen3v(ii,ie)
            hkv(kn)=hkv(kn)+hm3v(ii,ie)*weight	!ccf
            haux(kn)=haux(kn)+weight
          end do
        end do

        do k=1,nkn
          hkv(k) = hkv(k) / haux(k)
        end do

        end

c********************************************************************

        subroutine makehkv_minmax(hkv,haux,itype)

c makes hkv (nodewise depth)
c
c itype:  -1: min  0: aver  +1: max

	use basin

        implicit none

        real hkv(1)
        real haux(1)
        integer itype


	include 'param.h'

        integer k,ie,ii
        real hinit,h

c-------------------------------------------------------
c initialize
c-------------------------------------------------------

	if( itype .eq. 0 ) then
          call makehkv(hkv,haux)
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

c-------------------------------------------------------
c end of routine
c-------------------------------------------------------

        end

c********************************************************************

	subroutine depadj(hmin,hmax,href)

c adjusts depth to reference and min/max values - only hm3v is changed

	use basin

	implicit none

	real hmin,hmax,href

	include 'param.h'

	integer iaux,ie,ii
	real hmed

c adjust depth to constant in element %%%%%%%%%%%%%%%%%%%%%%

c        do ie=1,nel
c          hmed=0.
c          do ii=1,3
c            hmed=hmed+hm3v(ii,ie)
c          end do
c          hmed=hmed/3.
c          do ii=1,3
c            hm3v(ii,ie)=hmed
c          end do
c        end do

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

	end

c********************************************************************

	subroutine flatten_hm3v(hsigma)

	use basin

	implicit none

	real hsigma

	include 'param.h'

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

	include 'param.h'

	real v1v(nkn)

        call makehkv(hkv,v1v)		!computes hkv as average
        call makehkv_minmax(hkv_min,v1v,-1)
        call makehkv_minmax(hkv_max,v1v,+1)

	end

c********************************************************************

	subroutine make_hev

c adjusts elemental depth values

	use mod_depth

	implicit none

	include 'param.h'

        call makehev(hev)

	end

c********************************************************************

	subroutine adjourne_depth_from_hm3v

c adjourns hev and hkv from hm3v (if it has been changed)

	implicit none

        call make_hev
        call make_hkv
        !call set_last_layer		!FIXME

	end

c********************************************************************

	subroutine adjourn_depth_from_hev

c adjourns hev and hkv from hm3v (if it has been changed)

	use mod_depth
	use basin

	implicit none


	include 'param.h'

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

	include 'param.h'

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

	include 'param.h'

	integer ie

	open(1,file=file,status='old',form='formatted')
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

	include 'param.h'

	integer ie

	open(1,file=file,status='unknown',form='formatted')
	write(1,*) nel
	write(1,*) (hev(ie),ie=1,nel)
	close(1)

	end

c********************************************************************

