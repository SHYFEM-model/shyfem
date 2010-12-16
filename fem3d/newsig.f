c
c routines for sigma levels
c
c revision log :
c
c 16.12.2010    ggu     program partially finished
c
c********************************************************************

	subroutine init_sigma_vertical

	implicit none

	call check_nlv

	call set_sigma_hkhe

	end

c********************************************************************

	subroutine set_sigma_hkhe

	implicit none

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        integer nen3v(3,1)
        common /nen3v/nen3v
        real hm3v(3,1)
        common /hm3v/hm3v
        real hkv(1)
        common /hkv/hkv
        real hev(1)
        common /hev/hev

	logical berror
	integer k,ie,ii
	real flag,h,hm

c-------------------------------------------------------
c initialize with flag
c-------------------------------------------------------

	berror = .false.
	flag = -999.

	do k=1,nkn
	  hkv(k) = flag
	end do

c-------------------------------------------------------
c set hkv and hev
c-------------------------------------------------------

	do ie=1,nel
	  hm = 0.
	  do ii=1,3
	    h = hm3v(ii,ie)
	    k = nen3v(ii,ie)
	    if( hkv(k) .eq. flag ) then
	      hkv(k) = h
	    else if( hkv(k) .ne. h ) then
	      write(6,*) 'different depth: ',k,hkv(k),h
	      berror = .true.
	    end if
	    hm = hm + h
	  end do
	  hev(ie) = hm / 3.
	end do

	if( berror ) stop 'error stop set_sigma_hkhe: depth values'

c-------------------------------------------------------
c check hkv
c-------------------------------------------------------

	do k=1,nkn
	  if( hkv(k) .eq. flag ) then
	    write(6,*) 'flag in depth: ',k
	    berror = .true.
	  end if
	end do

	if( berror ) stop 'error stop set_sigma_hkhe: flag in depth'

c-------------------------------------------------------
c end of routine
c-------------------------------------------------------

	end

c********************************************************************

	subroutine set_sigma_levels

	implicit none

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer nlvdi,nlv
        common /level/ nlvdi,nlv

        integer ilhv(1)
        common /ilhv/ilhv
        integer ilhkv(1)
        common /ilhkv/ilhkv

        real hlv(1)
        common /hlv/hlv


	integer i,k,ie
	integer nsiglv
	real ds

	real getpar

	nsiglv = nint(getpar('nsiglv'))

	if( nsiglv .gt. 0 .and. nlv .gt. 0 ) then
	  write(6,*) 'nsiglv,nlv: ',nsiglv,nlv
	  write(6,*) 'Either nsiglv or single levels can be given'
	  stop 'error stop set_sigma_levels: nsiglv'
	else if( nsiglv .gt. nlvdi ) then
	  write(6,*) 'nsiglv,nlvdi: ',nsiglv,nlvdi
	  write(6,*) 'Number of requested sigma levels too high'
	  write(6,*) 'Please change dimension of layers: nlvdim'
	  stop 'error stop set_sigma_levels: nlvdi'
	end if

	if( nsiglv .gt. 0 ) then
	  nlv = nsiglv
	  ds = 1. / nlv
	  do i=1,nlv
	    hlv(i) = -i*ds
	  end do
	  hlv(nlv) = -1.		!just to be sure
	else
	  
	end if

	do k=1,nkn
	  ilhkv(k) = nlv
	end do

	do ie=1,nel
	  ilhv(ie) = nlv
	end do

	end

c********************************************************************

