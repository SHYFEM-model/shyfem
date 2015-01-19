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
c 05.08.2012	ggu	new name noswrite
c
c***************************************************************

	program noswrite

c writes nos file from data file

	implicit none

	include 'param.h'

	character title
	include 'basin.h'
	include 'simul.h'


	include 'levels.h'

	include 'hydro_print.h'

	logical bwater(neldim)
	integer iwet(neldim)
	integer iwetv(neldim)
	real hev(neldim)
	real weight(nlvdim,nkndim)

	real znv(nkndim)
	include 'hydro.h'

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
	integer it,ie,i,k,iu,ivar
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

	real cv3(nlvdim,nkndim)

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
	open(iu,file='node.dat',status='old',form='formatted')

c-------------------------------------------------------------------
c prepare file for output
c-------------------------------------------------------------------

	nbout = 55
	nvers = 3
	nlv = 1
	title = 'created from node file'
	ivar = 222

	hlv(1) = 10000.
	do k=1,nkn
	  ilhkv(k) = 1
	end do

        nbout=ifileo(nbout,'zuv.nos','unform','new')
        if(nbout.le.0) then
	  stop 'error stop: Cannot open NOS file for writing'
	end if

	call whnos(nbout,nvers,nkn,nel,nlv,1,ilhkv,hlv,hev,title)

c-------------------------------------------------------------------
c loop over records
c-------------------------------------------------------------------

  300   continue

c	---------------------------------
c	read record
c	---------------------------------

	call read_data(iu,it,nlvdim,nknaux,cv3,ierr)

	if( nkn .ne. nknaux ) goto 99

        if( ierr .gt. 0 ) goto 96
	if( ierr .ne. 0 ) goto 100

	nread=nread+1

	write(6,*) 
	write(6,*) 'time = ',it,'    nread = ',nread
	write(6,*) 

        call wrnos(nbout,it,ivar,nlvdim,ilhkv,cv3,ierr)

	goto 300
  100	continue

c-------------------------------------------------------------------
c end of routine
c-------------------------------------------------------------------

	stop
   99	continue
	write(6,*) 'nkn,nel: ',nkn,nel
	write(6,*) 'nknaux: ',nknaux
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

	subroutine read_data(iu,it,nlvdim,nknaux
     +			,cv3,ierr)

c reads file with values

	implicit none

	integer iu
	integer it
	integer nlvdim
	integer nknaux
	real cv3(nlvdim,1)
	integer ierr

	integer k,ie

	read(iu,*,end=99) it
	read(iu,*) nknaux
	read(iu,*) (cv3(1,k),k=1,nknaux)

	ierr = 0
	return
   99	continue
	ierr = -1
	return
	end

c******************************************************************

