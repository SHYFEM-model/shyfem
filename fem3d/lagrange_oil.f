
c*******************************************************************

	subroutine init_diff_oil

	use mod_lagrange
	use basin, only : nkn,nel,ngr,mbw

	implicit none

        include 'param.h'


	integer ie

	do ie=1,nel
	  rwhvar(ie) = rwhpar
	end do

	end

c*******************************************************************

	subroutine set_diff_oil

	use mod_lagrange
	use evgeom
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'femtime.h'

	real oil(nel)
	real zfilm(nel)
	real v1v(nkn)
	real v2v(nkn)

	integer ie,i
	integer ioil,ifreq
	real qoil,rhooil,ztresh,fact
	real area,roil,zoil,dz,rfact
	real rwhnew

	integer icall
	save icall
	data icall /0/

	ioil = 1
	ifreq = 6	!frequency of output

	qoil = 10.	!quantity of oil for one particle
	rhooil = 1000.
	ztresh = 0.001	!in meters
	fact = 1.

	if( ioil .le. 0 ) return

	icall = icall + 1

	do ie=1,nel
	  oil(ie) = 0.
	end do

        do i=1,nbdy
          ie=ie_body(i)
	  oil(ie) = oil(ie) + 1.
        end do

	do ie=1,nel
	  area = 12. * ev(10,ie)
	  roil = qoil * oil(ie)
	  zoil = roil / ( rhooil * area )

	  zfilm(ie) = zoil

	  rwhnew = rwhpar
	  if( zoil .gt. ztresh ) then
	    dz = zoil - ztresh
	    rfact = dz / ztresh
	    rwhnew = ( 1. + rfact*fact ) * rwhpar
	  end if
	  rwhvar(ie) = rwhnew

	  if( zoil .gt. 0. ) then
	    write(6,*) ie,zoil,ztresh,rwhnew
	  end if

	end do

        if( ifreq .gt. 0 .and. mod(icall,ifreq) .eq. 0 ) then
	  call e2n2d(zfilm,v1v,v2v)
          call wrnos2d_index(it,icall,'film','film thickness',v1v)
        end if

	end

c*******************************************************************

