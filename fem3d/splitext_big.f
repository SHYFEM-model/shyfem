c
c $Id: splitext.f,v 1.11 2009-02-13 17:22:44 georg Exp $
c
c revision log :
c
c 09.04.1999	ggu	restructured from readext
c 28.09.1999	ggu	reads now all data and then writes it
c 13.02.2009	ggu	big arrays in common for parallel environment
c
c***************************************************************

	program splitext

c This routine reads an EXT file and writes the data to unit 76
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

	implicit none

	integer noddim			!total number of nodes
	!parameter (noddim=27)
	!parameter (noddim=35)
	parameter (noddim=60)
	!parameter (noddim=20)

	character*80 line,file
	integer nvers,knausm
	integer nrec,it,i,nin,kn,in,nout,itold
	real href,hzmin
	real err,rm
	character*80 descrp

	integer knaus(noddim)
	real xv(3*noddim)
	real h(noddim)

	integer iapini, ifemop, ifileo
	real read7,rdrc7

c---------------------------------------------------------------
c open files
c---------------------------------------------------------------

        if(iapini(2,0,0,0).eq.0) then
                stop 'error stop : iapini'
        end if

	nin = ifemop('.ext','unform','old')
        if(nin.le.0) goto 97

c---------------------------------------------------------------
c read and write header information
c---------------------------------------------------------------

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

	if( knausm .gt. noddim ) goto 88

c---------------------------------------------------------------
c loop over data
c---------------------------------------------------------------

	nrec = 0

   10	continue

c	write(6,*) nin,nvers,it,knausm
        err = rdrc7(nin,nvers,it,knausm,xv)
	if( err .eq. -1. ) goto 100
	if( err .ne. 0. ) goto 98
	nrec = nrec + 1

	if(mod(nrec,100).eq.0) write(6,*) nrec,' data records read'

	itold = it

	do i=1,knausm
	  write(100+i,*) it,xv(i)
	  write(200+i,*) it,xv(i+knausm)
	  write(300+i,*) it,xv(i+2*knausm)
	  rm = sqrt( xv(i+knausm)**2 + xv(i+2*knausm)**2 )
	  write(400+i,*) it,rm
	end do

	goto 10
  100	continue

c---------------------------------------------------------------
c end of data
c---------------------------------------------------------------

	write(6,*)
	write(6,*) 'Total number of data records read : ',nrec
	write(6,*) 'Last time value read: ',itold
	write(6,*)

c---------------------------------------------------------------
c writing files
c---------------------------------------------------------------

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

	stop
   88	continue
	stop 'error stop: noddim'
   99	continue
	stop 'Error reading header of EXT file'
   98	continue
	write(6,*) err
	stop 'Error reading data record of EXT file'
   97	continue
	stop 'Cannot open EXT file'
	end

c*******************************************************************

	subroutine wrts(n,it,data,name,number)

c writes data to file name.number

	implicit none

	integer n
	integer it(1)
	real data(1)
	character*(*) name
	integer number

	integer i
	integer in,nout
	character*80 numlin,file
	integer ialfa

	nout = 1
	in = ialfa(float(number),numlin,-1,-1)
	file = name // '.' // numlin(1:in)

	open(nout,file=file,status='unknown',form='formatted')

	do i=1,n
	  write(nout,*) it(i),data(i)
	end do

	close(nout)

	end

c*******************************************************************

