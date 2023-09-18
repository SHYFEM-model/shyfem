!
! $Id: subdep.f,v 1.10 2008-07-16 15:41:39 georg Exp $
!
! depth utility routines
!
! contents :
!
! function igtdep(k,f)                  get depth for node
! function igtdpa(mode,h)               gets unique depth for all nodes
! subroutine huniqu(hev,hkv)            make depth unique for every node
! subroutine makehev(hev)		makes hev (elementwise depth)
! subroutine makehkv(hkv)		makes hkv (nodewise depth)
! subroutine depadj(hmin,hmax,href)	adjusts depth to ref/min/max values
!
! revision log :
!
! 29.06.1997	ggu	depth routines in one file
! 06.11.1998	ggu	new huniqu to compute hev and hkv
! 19.10.1999	ggu	new routine makehv from subutl
! 25.03.2002	ggu	new routines makehkv (before makehv) and makehev
! 28.11.2005	ggu	makehkv changed (uses double precision aux value, area weight)
! 24.02.2006	ggu	bug in makehkv -> haux was integer
! 18.10.2006	ccf	bug in makehkv -> no area multiplication
! 16.12.2010	ggu	in depadj() do not set hm3v to constant
! 17.05.2011	ggu	new routines to adjourn depth
! 18.11.2011	ggu	new routine makehkv_minmax()
! 05.09.2013	ggu	new routine set_sigma_hkv_and_hev() from newsig.f
! 25.06.2014	ggu	computa also hkv_min and hkv_max
! 25.05.2015	ggu	some changes in depth computation
!
!********************************************************************
!--------------------------------------------------------------------
        module depth_util_2nc
!--------------------------------------------------------------------

        use basin
        use para
        use depth
        use fem_util
        use sigma

!--------------------------------------------------------------------
        contains
!--------------------------------------------------------------------

	function igtdep(k,f,ndim)

! gets depth given a node number
! if different values of depth are associated
! with a node, it returns all these values
!
! k             node number (internal)
! f             vector in which the depth values are
!               ...stored at return
! ndim		dimension of f
! igtdep        number of different values found
!
! revised 02.02.94 by ggu	$$nmax - check error condtion nmax
! revised 29.06.97 by ggu	$$ndim - dimension of f is passed

	implicit none

	integer igtdep
	integer k,ndim
	double precision f(ndim)

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

!********************************************************************

	function igtdpa(mode,h)

! gets unique depth for all nodes
!
! mode          switch
!               1       deepest value is returned
!               -1      most shallow value
! h             vector in which the depth values are
!               ...stored (return value)
! igtdpa        return status
!               1       unique depth
!               0       no unique depth
!               -1      error

	implicit none

	integer igtdpa
	integer mode
	double precision h(nkn)

	double precision high
	parameter(high=1.e+30)

	logical buniq
	integer ie,ii,i,k
	double precision hh,hhh,hflag

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

!********************************************************************

	subroutine huniqu(hev,hkv)

! make depth unique for every node (changes hm3v)
! nodal values are the highest (deepest) value
!
! hev		element averaged depth values
! hkv            array with unique depth values

	implicit none

	double precision hev(nel)
	double precision hkv(nkn)

	integer ie,ii,k
	logical bstop
	double precision h,flag,hm

! flag nodal values

	flag = -999.

	do k=1,nkn
	  hkv(k) = flag
	end do

! create element averaged depth values and assign to nodal values

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

! check if all depth values are available

	bstop = .false.

	do k=1,nkn
	  if( hkv(k) .eq. flag ) then
		write(6,*) 'No depth for node ',ipext(k)
		bstop = .true.
	  end if
	end do

	if( bstop ) stop 'error stop huniqu'

	end

!********************************************************************
!********************************************************************
!********************************************************************

        subroutine makehev(hev)

! makes hev (elementwise depth)

        implicit none

! arguments
        double precision hev(nel)
! local
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

!********************************************************************

        subroutine makehkv(hkv)

! makes hkv (nodewise depth)

        use evgeom_2nc

        implicit none

! arguments
        double precision hkv(nkn)
! local
        integer ie,ii,k,kn
	double precision weight
        double precision haux(nkn)   !aux array -> bug - was integer

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

!********************************************************************

        subroutine makehkv_minmax(hkv,itype)

! makes hkv (nodewise depth)
!
! itype:  -1: min  0: aver  +1: max

        implicit none

        double precision hkv(nkn)
        integer itype

        integer k,ie,ii
        double precision hinit,h

!-------------------------------------------------------
! initialize
!-------------------------------------------------------

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

!-------------------------------------------------------
! set hkv
!-------------------------------------------------------

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

!-------------------------------------------------------
! end of routine
!-------------------------------------------------------

        end

!********************************************************************

	subroutine depadj(hmin,hmax,href)

! adjusts depth to reference and min/max values - only hm3v is changed

	implicit none

	double precision hmin,hmax,href

	integer iaux,ie,ii
	double precision hmed

! adjust depth to constant in element %%%%%%%%%%%%%%%%%%%%%%

!        do ie=1,nel
!          hmed=0.
!          do ii=1,3
!            hmed=hmed+hm3v(ii,ie)
!          end do
!          hmed=hmed/3.
!          do ii=1,3
!            hm3v(ii,ie)=hmed
!          end do
!        end do

! adjust depth to minimum depth %%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

! adjust depth to maximum depth %%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

! adjust depth to reference level %%%%%%%%%%%%%%%%%%%%%%%%%%%%

	do ie=1,nel
	 do ii=1,3
	  hm3v(ii,ie)=hm3v(ii,ie)-href
	 end do
	end do

	end

!********************************************************************

	subroutine adjust_depth

! adjusts depth values - only hm3v is changed

	implicit none

	double precision hmin,hmax,href

!       call bocche     !FIXME

        hmin=getpar('hmin')
        hmax=getpar('hmax')
        href=getpar('href')

        call depadj(hmin,hmax,href)	!adjusts h=h-href and hmax<h<hmin

	end

!********************************************************************

	subroutine set_depth

! sets up depth arrays

	implicit none

	logical bsigma
	integer nlv,nsigma
	double precision hsigma

	call get_sigma_info(nlv,nsigma,hsigma)
	bsigma = nsigma .gt. 0

	if( bsigma ) then	!sigma or hybrid layers
	  call check_sigma_hsigma
	  call flatten_hm3v(hsigma)
	else
	  call flatten_hm3v(-999.d0)
	end if

	call make_hev
	call make_hkv

	end

!********************************************************************

	subroutine flatten_hm3v(hsigma)

	implicit none

	double precision hsigma

	integer ie,ii
	double precision hm

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

!********************************************************************

	subroutine make_hkv

! adjusts nodal depth values

	use basin, only : nkn,nel,ngr,mbw

	implicit none

        call makehkv(hkv)		!computes hkv as average
        call makehkv_minmax(hkv_min,-1)
        call makehkv_minmax(hkv_max,+1)

	end

!********************************************************************

	subroutine make_hev

! adjusts elemental depth values

	implicit none

        call makehev(hev)

	end

!********************************************************************

	subroutine adjourne_depth_from_hm3v

! adjourns hev and hkv from hm3v (if it has been changed)

	implicit none

        call make_hev
        call make_hkv
        !call set_last_layer		!FIXME

	end

!********************************************************************

	subroutine adjourn_depth_from_hev

! adjourns hev and hkv from hm3v (if it has been changed)

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

!********************************************************************
!********************************************************************
!********************************************************************

	subroutine check_sigma_hsigma

! checks hkv and hsigma
! uses information about sigma layers and hsigma (hybrid)

	implicit none

	logical berror
	integer k,ie,ii
	integer inc,ihmin,ihmax
	integer nlv,nsigma
	double precision flag
	double precision hm,h
	double precision hsigma
	double precision hkv(nkn)		!local

!-------------------------------------------------------
! initialize
!-------------------------------------------------------

	call get_sigma_info(nlv,nsigma,hsigma)

	if( nsigma == 0 ) return

	flag = -999.

	do k=1,nkn
	  hkv(k) = flag
	end do

	berror = .false.
	inc = 0

!-------------------------------------------------------
! set if hkv is continuous in sigma layers
!-------------------------------------------------------

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

!-------------------------------------------------------
! check hsigma crossing
!-------------------------------------------------------

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

!-------------------------------------------------------
! end of routine
!-------------------------------------------------------

	end

!********************************************************************
!********************************************************************
!********************************************************************

	subroutine read_in_hev(file)

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

!********************************************************************

	subroutine write_out_hev(file)

	use basin, only : nkn,nel,ngr,mbw

	character*(*) file

	integer ie

	open(1,file=file,status='unknown',form='formatted')
	write(1,*) nel
	write(1,*) (hev(ie),ie=1,nel)
	close(1)

	end

!********************************************************************

!--------------------------------------------------------------------
        end module depth_util_2nc
!--------------------------------------------------------------------
