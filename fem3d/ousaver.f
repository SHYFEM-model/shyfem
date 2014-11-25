c
c $Id: ousaver.f,v 1.4 2009-04-07 10:43:57 georg Exp $
c
c interpolation of velocities onto nodes
c
c revision log :
c
c 02.09.2003	ggu	adapted to new OUS format
c 24.01.2005	ggu	computes maximum velocities for 3D (only first level)
c 04.03.2005	ggu	computes 3D velocities
c 16.10.2007	ggu	problems in z averaging shows blank areas -> zaver
c 14.04.2008	ggu	commented, blank area problems resolved
c 12.07.2011	ggu	some changes on how to treat partially dry areas
c
c***************************************************************

	program ousaver

c reads ous file and averages with frequency nfreq

	implicit none

	include 'param.h'

	character*80 descrr,descrp
	common /descrr/ descrr
	common /descrp/ descrp
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	real xgv(nkndim), ygv(nkndim)
	real hm3v(3,neldim)
	integer nen3v(3,neldim)
	integer ipev(neldim), ipv(nkndim)
	integer iarv(neldim)
	common /xgv/xgv, /ygv/ygv
	common /hm3v/hm3v
	common /nen3v/nen3v
	common /ipev/ipev, /ipv/ipv
	common /iarv/iarv

	integer ilhv(neldim)
	real hlv(nlvdim)
        real utlnv(nlvdim,neldim)
        real vtlnv(nlvdim,neldim)
	common /ilhv/ilhv
	common /hlv/hlv
        common /utlnv/utlnv
        common /vtlnv/vtlnv

        real uprv(nlvdim,nkndim)
        real vprv(nlvdim,nkndim)
        common /uprv/uprv
        common /vprv/vprv

	logical bwater(neldim)
	integer iwet(neldim)
	integer iwetv(neldim)
	real hev(neldim)
	real weight(nlvdim,nkndim)

	real znv(nkndim)
	real zenv(3,neldim)
	common /zenv/zenv

	real aux(nkndim)

	integer ndim
	parameter(ndim=100)
	real xpn(ndim), ypn(ndim)
	integer ielv(ndim)

	integer n,nx,ny
	integer nfreq,nrec
	integer ii,l,lmax,nbout
	integer icwet
	integer ifileo

        integer nvers,nin,nlv
        integer itanf,itend,idt,idtous
	integer it,ie,i
        integer ierr,nread,ndry
        integer nknous,nelous,nlvous
	integer ks
        real href,hzoff,hlvmin
	real zmin,zmax
	real umin,umax
	real vmin,vmax
	real cmin,cmax

c	integer rdous,rfous
	integer iapini,ideffi

	double precision zacum(3,neldim)
	double precision uacum(nlvdim,neldim)
	double precision vacum(nlvdim,neldim)

        real flag
        parameter ( flag = -999.0 )

c-------------------------------------------------------------------
c initialize parameters
c-------------------------------------------------------------------

	nfreq = 5

	nread = 0
	nrec  = 0

	ks = 6068	!special node for debug
	ks = 0

c-------------------------------------------------------------------
c read in basin and header of simulation
c-------------------------------------------------------------------

	if(iapini(3,nkndim,neldim,0).eq.0) then
		stop 'error stop : iapini'
	end if

	nin=ideffi('datdir','runnam','.ous','unform','old')
	if(nin.le.0) goto 100

	nvers=1
        call rfous(nin
     +			,nvers
     +			,nknous,nelous,nlvous
     +			,href,hzoff
     +			,descrp
     +			,ierr)

	nlv=nlvous
	call dimous(nin,nkndim,neldim,nlvdim)

        write(6,*)
        write(6,*)   descrp
        write(6,*)
        write(6,*) ' nvers        : ',nvers
        write(6,*) ' href,hzoff   : ',href,hzoff
        write(6,*) ' nkn,nel      : ',nknous,nelous
        write(6,*) ' nlv          : ',nlvous
        write(6,*)

	call rsous(nin,ilhv,hlv,hev,ierr)

c-------------------------------------------------------------------
c prepare file for output
c-------------------------------------------------------------------

	nbout = 55
        nbout=ifileo(nbout,'aver.ous','unform','new')
        if(nbout.le.0) then
	  stop 'error stop: Cannot open OUS file for writing'
	end if
	call wfous(nbout,nvers,nkn,nel,nlv,href,hzoff,descrp,ierr)
	if( ierr .ne. 0 ) goto 95
	call wsous(nbout,ilhv,hlv,hev,ierr)
	if( ierr .ne. 0 ) goto 95

c-------------------------------------------------------------------
c initialize arrays
c-------------------------------------------------------------------

	call aver(nread,nlvdim,neldim,nel
     +			,zacum,uacum,vacum,zenv,utlnv,vtlnv,iwet,iwetv)

c-------------------------------------------------------------------
c loop over records
c-------------------------------------------------------------------

  300   continue

c	---------------------------------
c	read record
c	---------------------------------

        call rdous(nin,it,nlvdim,ilhv,znv,zenv,utlnv,vtlnv,ierr)

        if( ierr .gt. 0 ) goto 96
	if( ierr .ne. 0 ) goto 100

	nread=nread+1

	write(6,*) 
	write(6,*) 'time = ',it,'    nread = ',nread
	write(6,*) 

c	---------------------------------
c	accumulate
c	---------------------------------

	call init_dry_mask(bwater)
	call set_dry_mask(bwater,znv,href,hzoff)

	icwet = 0
	do ie=1,nel
	  if( bwater(ie) ) iwet(ie) = iwet(ie) + 1
	  if( bwater(ie) ) icwet = icwet + 1
	  do ii=1,3
	    zacum(ii,ie) = zacum(ii,ie) + zenv(ii,ie)
	  end do
	  lmax = ilhv(ie)
	  do l=1,lmax
	    uacum(l,ie) = uacum(l,ie) + utlnv(l,ie)
	    vacum(l,ie) = vacum(l,ie) + vtlnv(l,ie)
	  end do
	end do

	!write(6,*) 'wet are ... ',icwet,nel

c	---------------------------------
c	elaborate if time is right
c	---------------------------------

	if( nfreq .gt. 0 .and. nread .eq. nfreq ) then
	  nrec = nrec + 1
	  call aver(nread,nlvdim,neldim,nel
     +			,zacum,uacum,vacum,zenv,utlnv,vtlnv,iwet,iwetv)
	  call zaver(nkndim,neldim,nkn,nel,nen3v,zenv,znv,iwetv,aux)
	  call debug_write_node(ks,it,nrec
     +		,nkndim,neldim,nlvdim,nkn,nel,nlv
     +		,nen3v,zenv,znv,utlnv,vtlnv)
          call wrous(nbout,it,nlvdim,ilhv,znv,zenv,utlnv,vtlnv,ierr)
	  write(6,*) 'averaged record written ',nread,nrec,it
          if(ierr.ne.0.) goto 95
	  nread = 0
	end if

	goto 300
  100	continue

c-------------------------------------------------------------------
c end of loop -> finish
c-------------------------------------------------------------------

	if( nread .ne. 0 ) then
	  nrec = nrec + 1
	  call aver(nread,nlvdim,neldim,nel
     +			,zacum,uacum,vacum,zenv,utlnv,vtlnv,iwet,iwetv)
	  call zaver(nkndim,neldim,nkn,nel,nen3v,zenv,znv,iwetv,aux)
	  call debug_write_node(ks,it,nrec
     +		,nkndim,neldim,nlvdim,nkn,nel,nlv
     +		,nen3v,zenv,znv,utlnv,vtlnv)
          call wrous(nbout,it,nlvdim,ilhv,znv,zenv,utlnv,vtlnv,ierr)
	  write(6,*) 'averaged record written ',nread,nrec,it
          if(ierr.ne.0.) goto 95
	end if

c-------------------------------------------------------------------
c end of routine
c-------------------------------------------------------------------

	stop
   96	continue
	write(6,*) 'ierr = ',ierr
	stop 'error stop: error in reading file'
   95	continue
	stop 'error stop: Cannot write output file OUS'
	end

c******************************************************************

	subroutine zaver(nkndim,neldim,nkn,nel,nen3v,zenv,znv,iwetv,aux)

c averages nodal z values -> must be done better, otherwise dry areas FIXME

	implicit none

	integer nkndim,neldim,nkn,nel
	integer nen3v(3,neldim)
	real zenv(3,neldim)
	real znv(nkndim)
	integer iwetv(neldim)
	real aux(nkndim)

	integer ie,ii,k
	integer icwet

	icwet = 0

	do k=1,nkn
	  znv(k) = 0.
	  aux(k) = 0.
	end do

	do ie=1,nel
	 if( iwetv(ie) .ne. 0 ) then
	  icwet = icwet + 1
	  do ii=1,3
	    k = nen3v(ii,ie)
	    znv(k) = znv(k) + zenv(ii,ie)
	    aux(k) = aux(k) + 1.
	  end do
	 end if
	end do

	do k=1,nkn
	  if( aux(k) .gt. 0. ) then
	    znv(k) = znv(k) / aux(k)
	  end if
	end do

	do ie=1,nel
	 if( iwetv(ie) .ne. 0 ) then
	  do ii=1,3
	    k = nen3v(ii,ie)
	    zenv(ii,ie) = znv(k)
	  end do
	 end if
	end do

	!write(6,*) 'icwet: ',icwet

	end

c******************************************************************

	subroutine aver(nread,nlvdim,neldim,nel
     +			,zacum,uacum,vacum,zenv,utlnv,vtlnv,iwet,iwetv)

c averages variables using accumulated values

	implicit none

	integer nread,nlvdim,neldim,nel
	double precision zacum(3,neldim)
	double precision uacum(nlvdim,neldim)
	double precision vacum(nlvdim,neldim)
	real zenv(3,neldim)
	real utlnv(nlvdim,neldim)
	real vtlnv(nlvdim,neldim)
	integer iwet(neldim)
	integer iwetv(neldim)

	logical bplotdry
	integer ie,ii,l,nwet,nfact
	integer icwet
	double precision rr

c	-------------------------------
c	initialize parameter
c	-------------------------------

c	set to true if you want averaging be done in partially dry areas
c	if set to true only in always wet areas the averaging is performed

	bplotdry = .false.	!average in partially dry areas

c	-------------------------------
c	initialize factor
c	-------------------------------

	if( nread .le. 0 ) then
	  rr = 0.
	  nwet = 1
	else
	  rr = 1./nread
	  nwet = nread
	end if

	nfact = nread		  !this makes iwetv=1 if always wet, 0 else
	if( bplotdry ) nfact = 1  !this makes iwetv=0 only if always dry
	if( nfact .le. 0 ) nfact = 1

c	-------------------------------
c	compute average
c	-------------------------------

	icwet = 0
	do ie=1,nel
	  iwetv(ie) = iwet(ie) / nfact		!is 1 if always wet and 0 else
	  if( iwet(ie) .eq. nwet ) icwet = icwet + 1
	  do ii=1,3
	    zenv(ii,ie) = zacum(ii,ie) * rr
	  end do
	  do l=1,nlvdim
	    utlnv(l,ie) = uacum(l,ie) * rr
	    vtlnv(l,ie) = vacum(l,ie) * rr
	  end do
	end do

	!write(6,*) 'icwet = ',icwet

c	-------------------------------
c	initialize arrays
c	-------------------------------

	do ie=1,nel
	  iwet(ie) = 0
	  do ii=1,3
	    zacum(ii,ie) = 0.
	  end do
	  do l=1,nlvdim
	    uacum(l,ie) = 0.
	    vacum(l,ie) = 0.
	  end do
	end do

c	-------------------------------
c	end of routine
c	-------------------------------

	end

c******************************************************************

	subroutine debug_write_node_to_be_deleted(it,nrec
     +		,nkndim,neldim,nlvdim,nkn,nel,nlv
     +		,nen3v,zenv,znv,utlnv,vtlnv)

c debug write ...

	implicit none

	integer it,nrec
	integer nkndim,neldim,nlvdim,nkn,nel,nlv
	integer nen3v(3,neldim)
	real znv(nkndim)
	real zenv(3,neldim)
	real utlnv(nlvdim,neldim)
	real vtlnv(nlvdim,neldim)

	integer ie,ii,k,l,ks
	logical bk

	ks = 6068

	write(66,*) 'time: ',it,nrec
	write(66,*) 'kkk: ',znv(ks)

	do ie=1,nel
	  bk = .false.
	  do ii=1,3
	    k = nen3v(ii,ie)
	    if( k .eq. ks ) then
	      write(66,*) 'ii: ',ii,ie,zenv(ii,ie)
	      bk = .true.
	    end if
	  end do
	  if( bk ) then
	  do l=1,nlv
	    write(66,*) 'ie: ',ie,l,utlnv(l,ie),vtlnv(l,ie)
	  end do
	  end if
	end do

	end

c******************************************************************

