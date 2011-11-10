
c utilities for NOS files
c
c revision log :
c
c 29.04.2010    ggu     new file from nosaver
c 07.05.2010    ggu     new routines qopen_nos_file(), checks ev initialization
c 15.12.2010    ggu     volume computation also for sigma layers
c 10.11.2011    ggu     new routines for hybrid levels, init_volume() changed
c
c***************************************************************

	subroutine make_vert_aver(nlvdim,nkn,ilhkv,cv3,vol3,cv2)

	implicit none

	integer nlvdim
	integer nkn
	integer ilhkv(1)
	real cv3(nlvdim,1)
	real vol3(nlvdim,1)
	real cv2(1)

	integer k,l,lmax
	double precision c,v
	double precision cctot,vvtot

	do k=1,nkn
	  cctot = 0.
	  vvtot = 0.
	  lmax = ilhkv(k)
	  do l=1,lmax
	    c = cv3(l,k)
	    v = vol3(l,k)
	    cctot = cctot + c*v
	    vvtot = vvtot + v
	  end do
	  cctot = cctot / vvtot
	  cv2(k) = cctot
	end do

	end

c***************************************************************

	subroutine make_aver(nlvdim,nkn,ilhkv,cv3,vol3
     +				,cmin,cmax,cmed,vtot)

	implicit none

	integer nlvdim
	integer nkn
	integer ilhkv(1)
	real cv3(nlvdim,1)
	real vol3(nlvdim,1)
	real cmin,cmax,cmed,vtot

	integer k,l,lmax
	real c,v
	double precision cctot,vvtot

	cmin = cv3(1,1)
	cmax = cv3(1,1)
	cctot = 0.
	vvtot = 0.

	do k=1,nkn
	  lmax = ilhkv(k)
	  do l=1,lmax
	    c = cv3(l,k)
	    v = vol3(l,k)
	    cmin = min(cmin,c)
	    cmax = max(cmax,c)
	    cctot = cctot + c*v
	    vvtot = vvtot + v
	  end do
	end do

	cmed = cctot / vvtot
	vtot = vvtot

	end

c***************************************************************

	subroutine make_acumulate(nlvdim,nkn,ilhkv,cv3,cvacu)

	implicit none

	integer nlvdim
	integer nkn
	integer ilhkv(1)
	real cv3(nlvdim,1)
	real cvacu(nlvdim,1)

	integer k,l,lmax

        do k=1,nkn
          lmax = ilhkv(k)
          do l=1,lmax
            cvacu(l,k) = cvacu(l,k) + cv3(l,k)
          end do
        end do

	end

c***************************************************************
c***************************************************************
c***************************************************************

	subroutine open_nos_file(name,status,nunit)

	implicit none

	character*(*) name,status
	integer nunit

	integer nb
	character*80 file

        integer ifileo

	call mkname(' ',name,'.nos',file)
	nb = ifileo(0,file,'unform',status)
	if( nb .le. 0 ) stop 'error stop open_nos_file: opening file'

	nunit = nb

	end

c***************************************************************

        subroutine qopen_nos_file(text,status,nunit)

c asks for name and opens nos file

        implicit none

        character*(*) text,status
        integer nunit

        character*80 name

        write(6,*) text
        read(5,'(a)') name
        write(6,*) name

        call open_nos_file(name,status,nunit)

        end

c***************************************************************
c***************************************************************
c***************************************************************

	subroutine init_volume(nlvdim,nkn,nel,nlv,nen3v,ilhkv
     +				,hlv,hev,hl,vol3)

c initializes volumes just in case no volume file is found
c
c we just set everything to 1.
c we could do better using information on node area and depth structure

	implicit none

	integer nlvdim
	integer nkn,nel,nlv
	integer nen3v(3,1)
	integer ilhkv(1)
	real hlv(1)
	real hev(1)
	real hl(1)		!aux vector for layer thickness
	real vol3(nlvdim,1)

        include 'ev.h'

	logical belem,bsigma,bzeta
	integer ie,ii,k,l,lmax,nsigma
	real area,hsigma

	call is_init_ev(belem)		!use ev if initialized
        call get_sigma_info(nsigma,hsigma)
        bsigma = nsigma .gt. 0
	bzeta = .false.

	do k=1,nkn
	  lmax = ilhkv(k)
	  do l=1,nlv
	    vol3(l,k) = 0.
	  end do
	end do

	do ie=1,nel
	  area = 1.
	  if( belem ) area = 4. * ev(10,ie)
	  call get_layer_thickness(ie,nlv,bzeta,nsigma,hsigma,hl)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    lmax = ilhkv(k)
	    do l=1,lmax
	      vol3(l,k) = vol3(l,k) + area * hl(l)
	    end do
	  end do
	end do

	end

c***************************************************************

	subroutine get_volume(nvol,it,nlvdim,ilhkv,vol3)

c reads volumes

	implicit none

	integer nvol
	integer it
	integer nlvdim
	integer ilhkv(1)
	real vol3(nlvdim,1)

	integer ivar,ierr

	integer icall,itold
	save icall,itold
	data icall,itold /0,0/

	if( icall .gt. 0 .and. it .eq. itold ) return	!already read
	if( icall .le. 0 ) itold = it - 1		!force read
	if( it .lt. itold ) goto 95			!error - it too small

	do while( itold .lt. it )
	  call rdnos(nvol,itold,ivar,nlvdim,ilhkv,vol3,ierr)
          if(ierr.gt.0) goto 94			!read error
          if(ierr.ne.0) goto 93			!EOF
          if(ivar.ne.66) goto 92		!ivar should be 66
	end do

	icall = 1
	if( it .lt. itold ) goto 95		!HACK - just temporary
	if( itold .ne. it ) goto 96

	return
   92	continue
	write(6,*) ivar,66
	stop 'error stop get_volume: wrong variable'
   93	continue
	stop 'error stop get_volume: EOF found reading nos file'
   94	continue
	stop 'error stop get_volume: error reading nos file'
   95	continue
	write(6,*) 'it in vol file is higher than requested: ',itold,it
	return		!FIXME -> should signal error
	stop 'error stop get_volume: it too small'
   96	continue
	write(6,*) it,itold
	stop 'error stop get_volume: no vol record for it'
	end

c***************************************************************
c***************************************************************
c***************************************************************

	subroutine check_equal_i(text,n,a1,a2)

c tests array to be equal

	implicit none

	character*(*) text
	integer n
	integer a1(1)
	integer a2(1)

	integer i,imin,imax

	do i=1,n
	  if( a1(i) .ne. a2(i) ) goto 99
	end do

	return
   99	continue
	imin = max(1,i-5)
	imax = min(n,i-5)
	write(6,*) 'first array: ',imin,i,imax
	write(6,*) (a1(i),i=imin,imax)
	write(6,*) 'second array: ',imin,i,imax
	write(6,*) (a2(i),i=imin,imax)
	write(6,*) 'arrays are not equal: ',text
	stop 'error stop check_iqual_i: arrays differ'
	end

c***************************************************************

	subroutine check_equal_r(text,n,a1,a2)

c tests array to be equal

	implicit none

	character*(*) text
	integer n
	real a1(1)
	real a2(1)

	integer i,imin,imax

	do i=1,n
	  if( a1(i) .ne. a2(i) ) goto 99
	end do

	return
   99	continue
	imin = max(1,i-5)
	imax = min(n,i-5)
	write(6,*) 'first array: ',imin,i,imax
	write(6,*) (a1(i),i=imin,imax)
	write(6,*) 'second array: ',imin,i,imax
	write(6,*) (a2(i),i=imin,imax)
	write(6,*) 'arrays are not equal: ',text
	stop 'error stop check_iqual_r: arrays differ'
	end

c***************************************************************

