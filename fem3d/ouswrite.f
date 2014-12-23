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

	program ouswrite

c writes ous file

	implicit none

	include 'param.h'

	character*80 descrr,descrp
	common /descrr/ descrr
	common /descrp/ descrp
	include 'nbasin.h'

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
	integer it,ie,i,k,iu
        integer ierr,nread,ndry
        integer nknous,nelous,nlvous
	integer nknaux,nelaux
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

c-------------------------------------------------------------------
c read in basin and header of simulation
c-------------------------------------------------------------------

	if(iapini(1,nkndim,neldim,0).eq.0) then
		stop 'error stop : iapini'
	end if

	call makehev(hev)

	iu = 1
	open(iu,file='zuv.dat',status='old',form='formatted')

c-------------------------------------------------------------------
c prepare file for output
c-------------------------------------------------------------------

	nbout = 55
	nvers = 1
	nlv = 1
	href = 0.
	hzoff = 5.
	descrp = 'created from zuv file'

	hlv(1) = 10000.
	do ie=1,nel
	  ilhv(ie) = 1
	end do

        nbout=ifileo(nbout,'zuv.ous','unform','new')
        if(nbout.le.0) then
	  stop 'error stop: Cannot open OUS file for writing'
	end if
	call wfous(nbout,nvers,nkn,nel,nlv,href,hzoff,descrp,ierr)
	if( ierr .ne. 0 ) goto 95
	call wsous(nbout,ilhv,hlv,hev,ierr)
	if( ierr .ne. 0 ) goto 95

c-------------------------------------------------------------------
c loop over records
c-------------------------------------------------------------------

  300   continue

c	---------------------------------
c	read record
c	---------------------------------

	call read_data(iu,it,nlvdim,nknaux,nelaux,znv,utlnv,vtlnv,ierr)

	if( nkn .ne. nknaux ) goto 99
	if( nel .ne. nelaux ) goto 99
	if( nlv .ne. 1 ) goto 98

        if( ierr .gt. 0 ) goto 96
	if( ierr .ne. 0 ) goto 100

	nread=nread+1

	write(6,*) 
	write(6,*) 'time = ',it,'    nread = ',nread
	write(6,*) 

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    zenv(ii,ie) = znv(k)
	  end do
	end do

        call wrous(nbout,it,nlvdim,ilhv,znv,zenv,utlnv,vtlnv,ierr)

	goto 300
  100	continue

c-------------------------------------------------------------------
c end of routine
c-------------------------------------------------------------------

	stop
   99	continue
	write(6,*) 'nkn,nel: ',nkn,nel
	write(6,*) 'nknaux,nelaux: ',nknaux,nelaux
	stop 'error stop: parameters'
   98	continue
	write(6,*) 'nlv = ',nlv
	stop 'error stop: only nlv = 1 allowed'
   96	continue
	write(6,*) 'ierr = ',ierr
	stop 'error stop: error in reading file'
   95	continue
	stop 'error stop: Cannot write output file OUS'
	end

c******************************************************************

	subroutine read_data(iu,it,nlvdim,nknaux,nelaux
     +			,znv,utlnv,vtlnv,ierr)

c averages nodal z values -> must be done better, otherwise dry areas FIXME

	implicit none

	integer iu
	integer it
	integer nlvdim
	integer nknaux,nelaux
	real znv(1)
	real utlnv(nlvdim,1)
	real vtlnv(nlvdim,1)
	integer ierr

	integer k,ie

	read(iu,*,end=99) it
	read(iu,*) nknaux,nelaux
	read(iu,*) (znv(k),k=1,nknaux)
	read(iu,*) (utlnv(1,ie),ie=1,nelaux)
	read(iu,*) (vtlnv(1,ie),ie=1,nelaux)

	ierr = 0
	return
   99	continue
	ierr = -1
	return
	end

c******************************************************************

