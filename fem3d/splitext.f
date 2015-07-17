c
c $Id: splitext.f,v 1.11 2009-02-13 17:22:44 georg Exp $
c
c revision log :
c
c 09.04.1999	ggu	restructured from readext
c 28.09.1999	ggu	reads now all data and then writes it
c 13.02.2009	ggu	big arrays in common for parallel environment
c 16.11.2012	ggu	allow for direct writing values to files (bdirect)
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
	!parameter (noddim=60)
	parameter (noddim=90)
	!parameter (noddim=20)
	integer datdim			!total number of data records
	!parameter (datdim=120000)
	!parameter (datdim=80000)
	parameter (datdim=100000)

	character*80 line,file
	integer nvers,knausm,nb
	integer nrec,it,i,nin,kn,in,nout
	logical bdirect,bspecial
	real href,hzmin
	real err
	character*80 descrp

	integer knaus(noddim)
	real xv(3*noddim)
	real h(noddim)

	integer itime(datdim)
	real udata(datdim,noddim)
	real vdata(datdim,noddim)
	real zdata(datdim,noddim)
	real mdata(datdim,noddim)

	integer iapini, ifemop
	real read7,rdrc7

c---------------------------------------------------------------
c set parameter
c---------------------------------------------------------------

	bdirect = .false.
	bdirect = .true.	!write directly after read (for big files)
	bspecial = .true.	!special output
	bspecial = .false.

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
c open output files
c---------------------------------------------------------------

	if( bdirect ) then
	  do i=1,knausm
	    call opents(100+i,'u',i)
	    call opents(200+i,'v',i)
	    call opents(300+i,'z',i)
	    call opents(400+i,'m',i)
	  end do
	end if
c---------------------------------------------------------------
c loop over data
c---------------------------------------------------------------

	nb = knausm
	nrec = 0

   10	continue

c	write(6,*) nin,nvers,it,knausm
        err = rdrc7(nin,nvers,it,knausm,xv)
	if( err .eq. -1. ) goto 100
	if( err .ne. 0. ) goto 98
	nrec = nrec + 1

	if( .not. bdirect .and. nrec .gt. datdim ) then
	  write(6,*) 'Cannot read more than ',datdim,' data records'
	  nrec = nrec - 1
	  goto 100
	end if
	if(mod(nrec,100).eq.0) write(6,*) nrec,' data records read'

	if( bdirect ) then
	  do i=1,knausm
	    write(100+i,*) it,xv(i)		!u
	    write(200+i,*) it,xv(i+nb)		!v
	    write(300+i,*) it,xv(i+2*nb)	!z
	    write(400+i,*) it,sqrt( xv(i)**2 + xv(i+nb)**2 )
	  end do
	  if( bspecial .and. mod(it,3600) .eq. 0 ) then	!special output
	    write(10,*) it,xv(4+2*nb)
	  end if
	else
	  itime(nrec) = it
	  do i=1,knausm
	    udata(nrec,i) = xv(i)
	    vdata(nrec,i) = xv(i+knausm)
	    zdata(nrec,i) = xv(i+2*knausm)
	    mdata(nrec,i) = sqrt( udata(nrec,i)**2 + vdata(nrec,i)**2 )
	  end do
	end if

	goto 10
  100	continue

c---------------------------------------------------------------
c end of data
c---------------------------------------------------------------

	write(6,*)
	write(6,*) 'Total number of data records read : ',nrec
	if( .not. bdirect ) then
	  write(6,*) 'Last time value read: ',itime(nrec)
	end if
	write(6,*)

c---------------------------------------------------------------
c writing files
c---------------------------------------------------------------

	if( .not. bdirect ) then
	  do i=1,knausm
	    call wrts(nrec,itime,udata(1,i),'u',i)
	    call wrts(nrec,itime,vdata(1,i),'v',i)
	    call wrts(nrec,itime,zdata(1,i),'z',i)
	    call wrts(nrec,itime,mdata(1,i),'m',i)
	  end do
	end if

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

	stop
   88	continue
	stop 'error stop: noddim'
   89	continue
	stop 'error stop: datdim'
   99	continue
	stop 'Error reading header of EXT file'
   98	continue
	write(6,*) err
	stop 'Error reading data record of EXT file'
   97	continue
	stop 'Cannot open EXT file'
	end

c*******************************************************************

	subroutine opents(iunit,name,number)

c opens file name.number on unit iunit

	implicit none

	integer iunit
	character*(*) name
	integer number

	integer in,nout
	character*80 numlin,file
	integer ialfa

	nout = iunit
	in = ialfa(float(number),numlin,-1,-1)
	file = name // '.' // numlin(1:in)

	open(nout,file=file,status='unknown',form='formatted')

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
	call opents(nout,name,number)

	do i=1,n
	  write(nout,*) it(i),data(i)
	end do

	close(nout)

	end

c*******************************************************************

