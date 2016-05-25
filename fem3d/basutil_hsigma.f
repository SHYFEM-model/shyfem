
c*******************************************************************

	subroutine bashsigma(hsigma)

	use basin
	use mod_depth

	implicit none

	real hsigma

	integer idm,idp,iadjust
	real :: hflag = -999.

c-----------------------------------------------------------------
c check if depth values are continuous
c-----------------------------------------------------------------

	call check_depth(hflag)

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
        call makehkv_minmax(hkv,-1)	!use minimum depth
	call check_depth(hflag)

c-----------------------------------------------------------------
c write to file
c-----------------------------------------------------------------

        call basin_to_grd
        call grd_write('bashsigma.grd')
	write(6,*) 'Basin written to file bashsigma.grd'

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

