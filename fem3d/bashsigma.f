c
c $Id: bashsigma.f,v 1.20 2010-03-22 15:29:31 georg Exp $
c
c revision log :
c
c 16.11.2011    ggu     copied from basinf.f and newsig.f
c
c****************************************************************

        program bashsigma

c changes basin to conform with hsigma value
c
c needs continuous depth and produces continuous depth

	use mod_depth
	use evgeom
	use basin

	implicit none

	include 'param.h'


	!include 'aux_array.h'
	real haux(nkndim)

	logical bnode,belem
	integer idm,idp,iadjust
	real hflag,hsigma

	integer iapini

c-----------------------------------------------------------------
c read in basin
c-----------------------------------------------------------------

	hflag = -999.

        if( iapini(1,nkndim,neldim,0) .le. 0 ) stop

	call bas_info

	write(6,*) 'Enter value for hsigma: '
	read(5,*) hsigma
	write(6,*) 'Vaue used for hsigma: ',hsigma

	call set_ev

c-----------------------------------------------------------------
c check if depth values are continuous
c-----------------------------------------------------------------

	call check_depth(hflag)

        call makehev(hev)
        call makehkv(hkv,haux)

	call check_continuous(hsigma,idm,idp)

        if( idm .gt. 0 .or. idp .gt. 0 ) then	!this might be relaxed later
          write(6,*) 'basin has non continuous depth values: ',idm,idp
          stop 'error stop hsigma: depth values not on nodes'
        end if

c-------------------------------------------------------
c eliminate hsigma crossing
c-------------------------------------------------------

        iadjust = 0
        call check_hsigma_crossing(hsigma,iadjust)
        write(6,*) 'must adjust elements for hybrid levels: ',iadjust
        call check_hsigma_crossing(hsigma,iadjust)

	if( iadjust .gt. 0 ) then
          write(6,*) 'still hsigma crossing: ',iadjust
          stop 'error stop hsigma: still hsigma crossing'
	end if

c-------------------------------------------------------
c finalize depth values
c-------------------------------------------------------

        call makehev(hev)
        call makehkv_minmax(hkv,haux,-1)	!use minimum depth
	call check_depth(hflag)

c-----------------------------------------------------------------
c write to file
c-----------------------------------------------------------------

        open(1,file='hsigma.grd',status='unknown',form='formatted')
        call wrgrd_hsigma(1,hkv,hev,10000.)
        close(1)
	write(6,*) 'Basin written to file hsigma.grd'

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	end

c*******************************************************************

	subroutine check_continuous(hsigma,idm,idp)

c checks for continuous depth

	use mod_depth
	use basin

	implicit none

	real hsigma
	integer idm
	integer idp

	include 'param.h'


	integer ie,ii,k
	real h,hm

        idm = 0
        idp = 0

        do ie=1,nel
          hm = hev(ie)
          do ii=1,3
            h = hm3v(ii,ie)
            k = nen3v(ii,ie)
            if( abs(hkv(k)-h) .gt. 1.e-3 ) then
              if( hm .lt. hsigma .and. h .ne. hsigma ) then
                idm = idm + 1
              else
                idp = idp + 1
              end if
            end if
          end do
        end do

	end

c*******************************************************************

	subroutine check_depth(hflag)

c checks for flag values

	use basin

	implicit none

	real hflag

	include 'param.h'

	integer ie,ii
	integer iflag
	real h,hmin,hmax

	iflag = 0
	hmin = hm3v(1,1)
	hmax = hm3v(1,1)

        do ie=1,nel
          do ii=1,3
            h = hm3v(ii,ie)
	    if( h .eq. hflag ) iflag = iflag + 1
	    hmin = min(hmin,h)
	    hmax = max(hmax,h)
          end do
        end do

	write(6,*) 'hm3v iflag/hmin,hmax: ',iflag,hmin,hmax

	end

c*******************************************************************

        subroutine check_hsigma_crossing(hsigma,iadjust)

c checks and adjusts hsigma crossing

	use basin

        implicit none

        real hsigma
        integer iadjust

	include 'param.h'

	!include 'aux_array.h'

        logical berror,bdebug
        integer k,ie,ii
        integer ihmin,ihmax
        real h,hm
        real f(3)
        real v1v(nkn)

        bdebug = .true.
        bdebug = .false.
        iadjust = 0

        do k=1,nkn
          v1v(k) = 0.
        end do

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
            iadjust = iadjust + 1
            if(bdebug) write(6,*) 'before ',ie,(hm3v(ii,ie),ii=1,3)
            call adjust_for_hybrid(hsigma,hm3v(1,ie),f)
            if(bdebug) write(6,*) 'after  ',ie,(hm3v(ii,ie),ii=1,3)
            do ii=1,3
              k = nen3v(ii,ie)
              v1v(k) = v1v(k) + f(ii)
            end do
          end if
        end do

        do ie=1,nel
          do ii=1,3
            h = hm3v(ii,ie)
            k = nen3v(ii,ie)
	    if( v1v(k) .ne. 0. ) hm3v(ii,ie) = hsigma
          end do
        end do

        end

c*******************************************************************

        subroutine adjust_for_hybrid(hsigma,h,f)

c adjusts depth values for hybrid levels

        implicit none

        real hsigma
        real h(3),f(3)

        logical baver
        integer ii,iii,ip
        real hm
        real dh1,dh2

        baver = .false.

        hm = 0.
        do ii=1,3
          hm = hm + h(ii)
          f(ii) = 0.
        end do
        hm = hm / 3.

        if( baver ) then                !decide based on average

        if( hm .gt. hsigma ) then
          do ii=1,3
            if( h(ii) .lt. hsigma ) then
              h(ii) = hsigma
              f(ii) = 1.
            end if
          end do
        else
          do ii=1,3
            if( h(ii) .gt. hsigma ) then
              h(ii) = hsigma
              f(ii) = 1.
            end if
          end do
        end if

        else                            !decide based on min/max

        do ii=1,3
          iii = mod(ii,3) + 1
          dh1 = h(ii)-hsigma
          dh2 = h(iii)-hsigma
          if( dh1*dh2 .lt. 0. ) then
            if( abs(dh1) .lt. abs(dh2) ) then
              ip = ii
            else
              ip = iii
            end if
            h(ip) = hsigma
            f(ip) = 1.
          end if
        end do

        end if

        end

c*******************************************************************
c*******************************************************************
c*******************************************************************

        function areatr(ie)

c determination of area of element
c
c ie            number of element (internal)
c areatr        element area (return value)

	use basin

	real areatr
	integer ie

	include 'param.h'

	integer ii,i1,i2,k1,k2
	double precision f,x(3),y(3)

        do ii=1,3
          k=nen3v(ii,ie)
	  x(ii) = xgv(k)
	  y(ii) = ygv(k)
        end do

	f = (x(2)-x(1))*(y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1))

        areatr = f / 2.D0

        end

c*******************************************************************

	subroutine write_grd_from_bas

c writes grd file extracting info from bas file

	use basin

	implicit none

	include 'param.h'

	integer k,ie,ii,ia
	real x,y,h

	open(8,file='bas.grd',status='unknown',form='formatted')

	do k=1,nkn
	  x = xgv(k)
	  y = ygv(k)
	  !write(8,2000) 1,k,0,x,y
	  write(8,2000) 1,ipv(k),0,x,y
	end do

	do ie=1,nel
	  h = hm3v(1,ie)
	  ia = iarv(ie)
	  !write(8,2100) 2,ie,0,3,(nen3v(ii,ie),ii=1,3),h
	  write(8,2100) 2,ipev(ie),ia,3,(ipv(nen3v(ii,ie)),ii=1,3),h
	end do
	  
	close(8)

	return
 2000	format(i1,i10,i5,2e14.6)
 2100	format(i1,i10,2i5,3i10,e14.6)
	end

c*******************************************************************

        subroutine wrgrd_hsigma(iunit,hkv,hev,hsigma)

c writes grd file from bas
c
c hev and hkv must be set

	use basin

        implicit none

        integer iunit
        real hsigma
        real hkv(1)
        real hev(1)

        include 'param.h'

        integer k,ie,ii
        real h

        do k=1,nkn
          h = hkv(k)
          if( h .gt. hsigma ) then
            write(iunit,1000) 1,ipv(k),0,xgv(k),ygv(k)
          else
            write(iunit,1000) 1,ipv(k),0,xgv(k),ygv(k),hkv(k)
          end if
        end do

        write(iunit,*)

        do ie=1,nel
          h = hev(ie)
          if( h .gt. hsigma ) then
            write(iunit,1100) 2,ipev(ie),iarv(ie)
     +          ,3,(ipv(nen3v(ii,ie)),ii=1,3),hev(ie)
          else
            write(iunit,1100) 2,ipev(ie),iarv(ie)
     +          ,3,(ipv(nen3v(ii,ie)),ii=1,3)
          end if
        end do

        return
 1000   format(i1,2i10,3e16.8)
 1100   format(i1,2i10,i4,3i10,e16.8)
        end

c********************************************************************

