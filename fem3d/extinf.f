c
c $Id: extinf.f,v 1.2 2003/03/25 14:08:54 georg Exp $
c
c revision log :
c
c 09.04.1999	ggu	restructured from readext
c 28.09.1999	ggu	reads now all data and then writes it
c
c***************************************************************

	program extinf

c This routine reads an EXT file and writes information

	implicit none

	integer noddim			!total number of nodes
	parameter (noddim=200)

	character*80 line,file
	integer nvers,knausm
	integer nrec,it,i,nin,kn,in,nout
	real href,hzmin
	real err
	character*80 descrp

	integer itanf,itend
	integer knaus(noddim)
	real xv(3*noddim)
	real h(noddim)
	integer iznv(noddim)

	integer iapini,ifemop
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

	if( knausm .gt. noddim ) goto 89

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

	itend = it
	if( nrec .eq. 1 ) itanf = it

	!if(mod(nrec,100).eq.0) write(6,*) nrec,' data records read'
	do i=1,knausm
	  iznv(i) = nint(100.*xv(i+2*knausm))
	end do
	write(6,2000) it,(iznv(i),i=1,knausm)

	goto 10
  100	continue

c---------------------------------------------------------------
c end of data
c---------------------------------------------------------------

	write(6,*)
	write(6,*) 'Total number of data records read : ',nrec
	write(6,*)
	write(6,*) 'Start/End time : ',itanf,itend
	write(6,*)

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

	stop
 2000	format(i10,13i5)
   89	continue
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

