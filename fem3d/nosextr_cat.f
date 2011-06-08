c
c $Id: nosextr_records.f,v 1.1 2008-07-16 15:41:39 georg Exp $
c
c concatenates NOS files
c
c revision log :
c
c 18.11.1998    ggu     check dimensions with dimnos
c 03.06.2011    ggu     routine adjourned
c 08.06.2011    ggu     routine rewritten
c
c**********************************************************

	program nosextr_records

c extracts whole records from nos file
c
c records have to be specified on stdin

        implicit none

	include 'param.h'

        integer ndim
	parameter ( ndim = 10000 )

	integer iu(ndim)
	character*80 name,file

        integer nrdim
	parameter ( nrdim = 2000 )

	integer irec(nrdim)

	character*80 title
	real cv(nkndim)
	real cv3(nlvdim,nkndim)

        double precision accum(nlvdim,nkndim)
	real amin(nlvdim,nkndim)
	real amax(nlvdim,nkndim)

	integer ilhkv(nkndim)
	real hlv(nlvdim)
	real hev(neldim)

	logical ball,bwrite
        integer nread,ivarold,nextr
	integer nread1,nread2
        integer l,k,nin,nb
        integer nkn,nel,nlv,nvar
        integer it,ivar
	integer itfirst,itsecond
        integer ierr
        integer nvers
        real r,rnull
	real conz,high

        integer iapini,ideffi,ifileo

c-------------------------------------------------------------------
c initialize params
c-------------------------------------------------------------------

	itfirst = -1
	itsecond = -1
	nread=0
	nread1=0
	nread2=0
	nextr=0
	rnull=0.
        ivarold = 0
	high = 1.e+30

c-------------------------------------------------------------------
c open NOS file and read header
c-------------------------------------------------------------------

	call qopen_nos_file('Enter first file: ','old',nin)

        nvers=3
	call rfnos(nin,nvers,nkn,nel,nlv,nvar,title,ierr)
        if(ierr.ne.0) goto 100

        write(6,*) 'nvers    : ',nvers
        write(6,*) 'nkn,nel  : ',nkn,nel
        write(6,*) 'nlv,nvar : ',nlv,nvar
        write(6,*) 'title    : ',title

        call dimnos(nin,nkndim,neldim,nlvdim)

	call rsnos(nin,ilhkv,hlv,hev,ierr)
        if(ierr.ne.0) goto 100

	write(6,*) 'Available levels: ',nlv
	write(6,*) (hlv(l),l=1,nlv)

c-------------------------------------------------------------------
c open NOS output file
c-------------------------------------------------------------------

        call mkname(' ','nos_cat','.nos',file)
        write(6,*) 'writing file ',file(1:50)
        nb = ifileo(55,file,'unform','new')
        if( nb .le. 0 ) goto 98
        call wfnos(nb,3,nkn,nel,nlv,1,title,ierr)
        if( ierr .ne. 0 ) goto 99
        call wsnos(nb,ilhkv,hlv,hev,ierr)
        if( ierr .ne. 0 ) goto 99

c-------------------------------------------------------------------
c loop on input records
c-------------------------------------------------------------------

  300   continue

	call rdnos(nin,it,ivar,nlvdim,ilhkv,cv3,ierr)

        if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
        if(ierr.ne.0) goto 100

	itfirst = it
	nread=nread+1
	nread1=nread1+1
	write(6,*) 'time : ',nread,it,ivar

	call wrnos(nb,it,ivar,nlvdim,ilhkv,cv3,ierr)
	if( ierr .ne. 0 ) goto 99

	goto 300
  100	continue

	close(nin)
	call delnos(nin)

c-------------------------------------------------------------------
c open NOS file and read header
c-------------------------------------------------------------------

	call qopen_nos_file('Enter second file: ','old',nin)

        nvers=3
	call rfnos(nin,nvers,nkn,nel,nlv,nvar,title,ierr)
        if(ierr.ne.0) goto 100

        write(6,*) 'nvers    : ',nvers
        write(6,*) 'nkn,nel  : ',nkn,nel
        write(6,*) 'nlv,nvar : ',nlv,nvar
        write(6,*) 'title    : ',title

        call dimnos(nin,nkndim,neldim,nlvdim)

	call rsnos(nin,ilhkv,hlv,hev,ierr)
        if(ierr.ne.0) goto 100

	write(6,*) 'Available levels: ',nlv
	write(6,*) (hlv(l),l=1,nlv)

c-------------------------------------------------------------------
c loop on input records
c-------------------------------------------------------------------

  301   continue

	call rdnos(nin,it,ivar,nlvdim,ilhkv,cv3,ierr)

        if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
        if(ierr.ne.0) goto 101
	if( it .le. itfirst ) goto 301	!to make records unique

	if( itsecond .eq. -1 ) itsecond = it
	nread=nread+1
	nread2=nread2+1
	write(6,*) 'time : ',nread,it,ivar

	call wrnos(nb,it,ivar,nlvdim,ilhkv,cv3,ierr)
	if( ierr .ne. 0 ) goto 99

	goto 301
  101	continue

	close(nin)

c-------------------------------------------------------------------
c end of loop
c-------------------------------------------------------------------

	write(6,*)
	write(6,*) nread,' records read'
	write(6,*)
	write(6,*) 'nread1/2: ',nread1,nread2
	write(6,*) 'it: ',itfirst,itsecond
	write(6,*)

c-------------------------------------------------------------------
c end of routine
c-------------------------------------------------------------------

	stop
   91	continue
	write(6,*) 'file may have only one type of variable'
	write(6,*) 'error ivar : ',ivar,ivarold
	write(6,*) 'You should use nossplit to extract scalars first'
	stop 'error stop nosextr_records: ivar'
   98	continue
	write(6,*) 'error opening outout file'
	stop 'error stop nosextr_records'
   99	continue
	write(6,*) 'error writing file'
	stop 'error stop nosextr_records'
	end

c***************************************************************

