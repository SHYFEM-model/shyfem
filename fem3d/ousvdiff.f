c
c $Id: ousvdiff.f,v 1.2 2009-04-07 10:43:57 georg Exp $
c
c interpolation of velocities onto nodes
c
c revision log :
c
c 02.09.2003	ggu	adapted to new OUS format
c 24.01.2005	ggu	computes maximum velocities for 3D (only first level)
c 04.03.2005	ggu	computes 3D velocities
c 03.11.2010	ggu	prepared to use for water level differences
c
c***************************************************************

	program ousdiff

c reads ous files and computes difference

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
	common /ilhv/ilhv
	real hlv(nlvdim)
	common /hlv/hlv

        real utln1v(nlvdim,neldim)
        real vtln1v(nlvdim,neldim)
        common /utln1v/utln1v
        common /vtln1v/vtln1v
        real utln2v(nlvdim,neldim)
        real vtln2v(nlvdim,neldim)
        common /utln2v/utln2v
        common /vtln2v/vtln2v

        real uprv(nlvdim,nkndim)
        real vprv(nlvdim,nkndim)
        common /uprv/uprv
        common /vprv/vprv

	real hev(neldim)
        common /hev/hev

	real weight(nlvdim,nkndim)

	real zn1v(nkndim)
	real zen1v(3,neldim)
	real zn2v(nkndim)
	real zen2v(3,neldim)

	integer ndim
	parameter(ndim=100)
	real xpn(ndim), ypn(ndim)
	integer ielv(ndim)

	character*80 file1,file2
	integer n,nx,ny
	integer nfreq
	integer ii,l,lmax,nbout,nb1,nb2
	integer ifileo

	logical bnew
        integer nvers,nin,nlv
        integer itanf,itend,idt,idtous
	integer it,ie,i,k
	integer it1,it2
        integer ierr,nread,ndry
        integer nknous,nelous,nlvous
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

c-------------------------------------

	write(6,*) 'This program computes the difference of two ous files.'
	write(6,*) 'It needs the name of the basin and the two ous files.'
	write(6,*) 'The first is subtracted from the second.'

c-------------------------------------

	nread=0

	if(iapini(1,nkndim,neldim,0).eq.0) then
		stop 'error stop : iapini'
	end if

c-------------------------------------

	write(6,*) 'Enter name of first OUS file: '
	read (5,'(a)') file1

	call init_ousfile(nb1,file1
     +				,nkndim,neldim,nlvdim
     +				,nknous,nelous,nlvous
     +				,ilhv,hlv,hev
     +				,href,hzoff)

	if( nkn .ne. nknous ) goto 88
	if( nel .ne. nelous ) goto 88

c-------------------------------------

	write(6,*) 'Enter name of second OUS file: '
	read (5,'(a)') file2

	call init_ousfile(nb2,file2
     +				,nkndim,neldim,nlvdim
     +				,nknous,nelous,nlvous
     +				,ilhv,hlv,hev
     +				,href,hzoff)

	if( nkn .ne. nknous ) goto 88
	if( nel .ne. nelous ) goto 88

c-------------------------------------

	nbout = 55
	nvers = 1
        nbout=ifileo(nbout,'diff.ous','unform','new')
        if(nbout.le.0) then
	  stop 'error stop: Cannot open OUS file for writing'
	end if
	call wfous(nbout,nvers,nkn,nel,nlv,href,hzoff,descrp,ierr)
	if( ierr .ne. 0 ) goto 95
	call wsous(nbout,ilhv,hlv,hev,ierr)
	if( ierr .ne. 0 ) goto 95

c-------------------------------------

	it1 = -1024*1024*1024
	it2 = it1

  300   continue

	bnew = .true.
	do while( bnew .or. it1 .ne. it2 )
	  do while( bnew .or. it1 .lt. it2 )
            call rdous(nb1,it1,nlvdim,ilhv
     +			,zn1v,zen1v,utln1v,vtln1v,ierr)
	    if( ierr .gt. 0 ) goto 96	!error
	    if( ierr .lt. 0 ) goto 100	!EOF
	    bnew = .false.
	  end do
	  do while( bnew .or. it2 .lt. it1 )
            call rdous(nb2,it2,nlvdim,ilhv
     +			,zn2v,zen2v,utln2v,vtln2v,ierr)
	    if( ierr .gt. 0 ) goto 96	!error
	    if( ierr .lt. 0 ) goto 100	!EOF
	    bnew = .false.
	  end do
	end do

	nread=nread+1
	write(6,*) nread,it1,it2

	do k=1,nkn
	  !zn1v(k) = 0.5 * ( zn2v(k) + zn1v(k) )
	  zn1v(k) = zn2v(k) - zn1v(k)
	end do

	do ie=1,nel
	  do ii=1,3
	    !zen1v(ii,ie) = 0.5 * ( zen2v(ii,ie) + zen1v(ii,ie) )
	    zen1v(ii,ie) = zen2v(ii,ie) - zen1v(ii,ie)
	  end do
	  lmax = ilhv(ie)
	  do l=1,lmax
	    utln1v(l,ie) = utln2v(l,ie) - utln1v(l,ie)
	    vtln1v(l,ie) = vtln2v(l,ie) - vtln1v(l,ie)
	  end do
	end do

        call wrous(nbout,it1,nlvdim,ilhv,zn1v,zen1v,utln1v,vtln1v,ierr)
        if(ierr.ne.0.) goto 95

	write(6,*) 
	write(6,*) 'time : ',it
	write(6,*) 

	goto 300

  100	continue

	write(6,*)
	write(6,*) nread,' records elaborated'
	write(6,*)

	stop
   88	continue
	write(6,*) 'nkn, nel from basin: ',nkn,nel
	write(6,*) 'nkn, nel from simul: ',nknous,nelous
	stop 'error stop: simul parameters'
   96	continue
	write(6,*) 'error in reading file : ',ierr
	stop 'error stop: Cannot read file OUS'
   95	continue
	stop 'error stop: Cannot write output file OUS'
	end

c******************************************************************

	subroutine init_ousfile(nb,file
     +				,nkndim,neldim,nlvdim,nkn,nel,nlv
     +				,ilhv,hlv,hev
     +				,href,hzoff)

c initializes OUS file

	integer nb
	character*(*) file
	integer nkndim,neldim,nlvdim
	integer nkn,nel,nlv
	real href,hzoff

	character*80 descrp
	common /descrp/ descrp
	integer ilhv(neldim)
	real hlv(nlvdim)
	real hev(neldim)

	integer ierr,nvers

        nb=ifileo(55,file,'unform','old')
	if(nb.le.0) goto 100

	nvers=1
        call rfous(nb
     +			,nvers
     +			,nkn,nel,nlv
     +			,href,hzoff
     +			,descrp
     +			,ierr)
	if( ierr .ne. 0 ) goto 200

	call dimous(nb,nkndim,neldim,nlvdim)

	call rsous(nb,ilhv,hlv,hev,ierr)
	if( ierr .ne. 0 ) goto 200

        write(6,*)
        write(6,*)   descrp
        write(6,*)
        write(6,*) ' nvers        : ',nvers
        write(6,*) ' href,hzoff   : ',href,hzoff
        write(6,*) ' nkn,nel      : ',nkn,nel
        write(6,*) ' nlv          : ',nlv
        write(6,*)

	return
  100	continue
	write(6,*) file
	stop 'error stop init_ousfile: opening file'
  200	continue
	write(6,*) 'ierr = ',ierr
	stop 'error stop init_ousfile: reading header'
	end

c******************************************************************

