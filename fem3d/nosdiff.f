c
c $Id: nosdiff.f,v 1.3 2008-07-16 15:41:39 georg Exp $
c
c revision log :
c
c 18.11.1998    ggu     check dimensions with dimnos
c 06.04.1999    ggu     some cosmetic changes
c 03.12.2001    ggu     some extra output -> place of min/max
c 09.12.2003    ggu     check for NaN introduced
c 10.04.2008    ggu     finished program
c 07.05.2010    ggu     some utils transfered to nosutil.f
c
c**************************************************************

	program nosdiff

c computes difference between two NOS files - same grid

	implicit none

	include 'param.h'

	real cv3(nlvdim,nkndim)
	real cv3_new(nlvdim,nkndim)

	integer ilhkv(nkndim)
	real hlv(nlvdim)
	real hev(neldim)

	integer nread,nin1,nin2,nout
	integer nvers
	integer nkn,nel,nlv,nvar
	integer nkn2,nel2,nlv2,nvar2
	integer ierr
	integer it,ivar
	integer it1,ivar1
	integer it2,ivar2
	integer l,k
	character*80 title
	real cmin,cmax

	integer iapini
	integer ifem_open_file

c--------------------------------------------------------------
c initialize params
c--------------------------------------------------------------

	nread=0

c--------------------------------------------------------------
c open basin and simulations
c--------------------------------------------------------------

	!if(iapini(2,nkndim,neldim,0).eq.0) then
	!	stop 'error stop : iapini'
	!end if

	nread=0
	write(6,*) 'before...'

        call qopen_nos_file('Enter name of simulation 1: ','old',nin1)
        nvers=3
	call rhnos(nin1,nvers,nkndim,neldim,nlvdim,nkn,nel,nlv,nvar
     +				,ilhkv,hlv,hev,title)

        call qopen_nos_file('Enter name of simulation 2: ','old',nin2)
        nvers=3
	call rhnos(nin2,nvers,nkndim,neldim,nlvdim,nkn2,nel2,nlv2,nvar2
     +				,ilhkv,hlv,hev,title)

	if( nkn .ne. nkn2 ) goto 96
	if( nel .ne. nel2 ) goto 96
	if( nlv .ne. nlv2 ) goto 96
	if( nvar .ne. nvar2 ) goto 96

c-----------------------------------------------------------------
c open output files
c-----------------------------------------------------------------

        call open_nos_file('diff','new',nout)
        call whnos(nout,3,nkn,nel,nlv,1,ilhkv,hlv,hev,title)

c--------------------------------------------------------------
c loop on data
c--------------------------------------------------------------

	do while(.true.)

c	   ----------------------------------------------------
c	   read next record
c	   ----------------------------------------------------

	   call rdnos(nin1,it1,ivar1,nlvdim,ilhkv,cv3,ierr)
           if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
           if(ierr.ne.0) goto 100

	   call rdnos(nin2,it2,ivar2,nlvdim,ilhkv,cv3_new,ierr)
           if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
           if(ierr.ne.0) goto 100

	   if( it1 .ne. it2 ) goto 95
	   if( ivar1 .ne. ivar2 ) goto 95
	   it = it1
	   ivar = ivar1

	   nread=nread+1
	   write(6,*) 'time : ',it,'   ivar : ',ivar

c	   ----------------------------------------------------
c	   make difference
c	   ----------------------------------------------------

	   do l=1,nlv
	     do k=1,nkn
	       cv3(l,k) = cv3_new(l,k) - cv3(l,k)
	     end do
	   end do

c	   ----------------------------------------------------
c	   write to file
c	   ----------------------------------------------------

	  call wrnos(nout,it,ivar,nlvdim,ilhkv,cv3,ierr)

	end do	!do while

c--------------------------------------------------------------
c end of loop on data
c--------------------------------------------------------------

  100	continue

	write(6,*)
	write(6,*) 'output written to file diff.nos'
	write(6,*)
	write(6,*) nread,' records read'
	write(6,*)

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	stop
   95	continue
	write(6,*) 'it1,it2: ',it1,it2
	write(6,*) 'ivar1,ivar2: ',ivar1,ivar2
	stop 'error stop nosdiff: it or ivar differing'
   96	continue
	write(6,*) 'nkn,nkn2: ',nkn,nkn2
	write(6,*) 'nel,nel2: ',nel,nel2
	write(6,*) 'nlv,nlv2: ',nlv,nlv2
	write(6,*) 'nvar,nvar2: ',nvar,nvar2
	stop 'error stop nosdiff: simulation are different'
	end

c***************************************************************

