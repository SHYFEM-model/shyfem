c
c $Id: readext.f,v 1.9 2009-04-07 10:43:57 georg Exp $
c
c reads EXT file and writes data to several files
c
c notes :
c
c The header information is given as
c
c n    href  hzmin
c 1    k1    d1
c 2    k2    d2
c ...
c n    kn    dn
c
c where n is total number of nodes written, href is the reference
c value of the water levels and hzmin is the minimum depth for
c a node not to be dry. The following n lines give a consecutive
c number, the node number (k1,k2,...) as given in the STR file
c and the depth (d1,d2,...) at each node.
c
c The data records are written as
c
c it
c 1 u v z
c 2 u v z
c ...
c n u v z
c
c where it is the time stamp of the record and the
c next n lines give (u,v,z) for every node written.
c
c revision log :
c
c 16.09.1999    ggu	writes u/v/z/m.dat files + fort.76
c
c*****************************************************************

	program readext

	use basin

	implicit none

	integer noddim
	parameter (noddim=1000)

	include 'param.h'

	integer nvers,knausm
	integer nrec,it,i,nin,nout,nout2,kn
	integer nu,nv,nm,nz
	integer maxout
	real href,hzmin
	real err
	character*80 descrp

	integer knaus(noddim)
	real xv(3*noddim)
	real h(noddim)
	real u(noddim),v(noddim),z(noddim)
	real uv(noddim)
	integer iz(noddim)

c global data ------------------------------------
c ------------------------------------------------

	integer iapini, ideffi, ifileo
	real read7,rdrc7

	nrec = 0
	nout = 76
	nout2 = 77

c open files ---------------------------------------------------

        if(iapini(2,0,0,0).eq.0) then
                stop 'error stop : iapini'
        end if

        nin=ideffi('datdir','runnam','.ext','unform','old')
        if(nin.le.0) goto 97

	open(nout,file='fort.76',status='unknown',form='formatted')
	nu = ifileo(81,'u.dat','formatted','unknown')
	nv = ifileo(82,'v.dat','formatted','unknown')
	nm = ifileo(83,'m.dat','formatted','unknown')
	nz = ifileo(84,'z.dat','formatted','unknown')

c read and write header information ----------------------------

        err = read7(nin,noddim,nvers,knausm,knaus,h
     +                          ,href,hzmin,descrp)

	if( err .ne. 0. ) goto 99

	write(6,*) 'Title of simulation     : '
	write(6,*)
	write(6,*) descrp
	write(6,*)
	write(6,*) 'Version of file         : ',nvers
	write(6,*) 'Total number of nodes   : ',knausm
	write(6,*) 'Reference depth         : ',href
	write(6,*) 'Minimum depth           : ',hzmin
	write(6,*)
	write(6,*) '...reading data'

	write(nout,'(i10,2f20.8)') knausm,href,hzmin
	write(6,'(i10,2f20.8)') knausm,href,hzmin
	do i=1,knausm
	  kn = knaus(i)			!internal node number
	  write(nout,'(2i10,3f19.7)') i,ipv(kn),h(i),xgv(kn),ygv(kn)
	  write(6,'(2i10,3f19.7)') i,ipv(kn),h(i),xgv(kn),ygv(kn)
	end do

c loop over data -----------------------------------------------

   10	continue

        err = rdrc7(nin,nvers,it,knausm,xv)
	if( err .eq. -1. ) goto 100
	if( err .ne. 0. ) goto 98
	nrec = nrec + 1

	do i=1,knausm
	  u(i) = xv(i)
	  v(i) = xv(i+knausm)
	  z(i) = xv(i+2*knausm)
	  uv(i) = sqrt( u(i)*u(i) + v(i)*v(i) )
	end do

	if(mod(nrec,100).eq.0) write(6,*) nrec,' data records read'

	write(nout,*) it
	do i=1,knausm
	  write(nout,'(i10,3f20.8)') i,u(i),v(i),z(i)
	end do

	maxout = 20
	maxout=min(maxout,knausm)
	write(nu,1000) it,(u(i),i=1,maxout)
	write(nv,1000) it,(v(i),i=1,maxout)
	write(nm,1000) it,(uv(i),i=1,maxout)
	write(nz,1000) it,(z(i),i=1,maxout)

	goto 10
  100	continue

c end of data ------------------------------------------------

c	write(6,'(i10,5f10.4)') it,(v(i),i=1,knausm)

	write(6,*)
	write(6,*) 'Total number of data records read : ',nrec
	write(6,*)
	write(6,*) 'Data written to unit number 76'
	write(6,*)

	stop
 1000	format(i10,20f7.3)
   99	continue
	stop 'Error reading header of EXT file'
   98	continue
	stop 'Error reading data record of EXT file'
   97	continue
	stop 'Cannot open EXT file'
	end

c*****************************************************************

