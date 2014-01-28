c
c $Id: nos_elab.f,v 1.2 2009-04-07 10:43:57 georg Exp $
c
c revision log :
c
c 18.11.1998    ggu     check dimensions with dimnos
c 24.02.1999    ggu     use n2int for node number translation
c 03.12.2001    ggu     cleaned up, hakata bay
c
c****************************************************************

	program nos_elab

c extracts single nodes from nos file -> creates time series
c
c interactive version

	implicit none

	include 'param.h'

	integer ntrdim
	parameter (ntrdim=3)

c--------------------------------------------------
        character*80 descrr
        common /descrr/descrr
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        real xgv(nkndim), ygv(nkndim)
        real hm3v(3,neldim)
        integer nen3v(3,neldim)
        integer ipv(nkndim), ipev(neldim)
        integer iarv(neldim)

        common /xgv/xgv, /ygv/ygv
        common /hm3v/hm3v
        common /nen3v/nen3v
        common /ipv/ipv, /ipev/ipev
        common /iarv/iarv
c--------------------------------------------------

	character*80 title
	real cv(nkndim)
	real cv3(nlvdim,nkndim)
	real count(nlvdim,nkndim,ntrdim)
	real threshold(ntrdim)

	integer ilhkv(nkndim)
	real hlv(nlvdim)
	real hev(neldim)

	logical berror
	integer i,n,k,ke,l
	integer nread,nunit
	integer nvers
	integer nlv,nvar,ivar,ierr
	integer nin,it,it1,it2,nin1,nin2,nb
	integer ivar1,ivar2
	integer idt,itold
	integer ntres,nt
	real thres,tunit
	real c1,c2,c3
	real soglia

	integer iapini,ideffi,ialfa
	integer open_file,ifileo

c--------------------------------------------------
	character*50 file

c---------------------------------------------------------------
c nodes for extraction
c---------------------------------------------------------------

	integer ndim
	integer nnodes
	parameter( ndim = nkndim )
	integer nodes(ndim)	!node numbers
	integer nodese(ndim)	!external node numbers

	integer ifem_open_file

c---------------------------------------------------------------
c open simulation and basin
c---------------------------------------------------------------

	if(iapini(2,nkndim,neldim,0).eq.0) then
		stop 'error stop : iapini'
	end if

c---------------------------------------------------------------
c open NOS file and read header
c---------------------------------------------------------------

	nin = ifem_open_file('.nos','old')

        nvers=3
	call rfnos(nin,nvers,nkn,nel,nlv,nvar,title,ierr)
        if(ierr.ne.0) goto 100
        call dimnos(nin,nkndim,neldim,nlvdim)

        write(6,*) 'nvers    : ',nvers
        write(6,*) 'nkn,nel  : ',nkn,nel
        write(6,*) 'nlv,nvar : ',nlv,nvar
        write(6,*) 'title    : ',title

	call rsnos(nin,ilhkv,hlv,hev,ierr)
        if(ierr.ne.0) goto 100

	write(6,*) 'Available levels: ',nlv
	write(6,*) (hlv(l),l=1,nlv)

c---------------------------------------------------------------
c initializing units and nodes to be extracted
c---------------------------------------------------------------

        call mkname(' ','elab','.nos',file)
        write(6,*) 'writing file ',file(1:50)
        nb = ifileo(0,file,'unform','new')
        if( nb .le. 0 ) goto 98

        call wfnos(nb,3,nkn,nel,nlv,1,title,ierr)
        if( ierr .ne. 0 ) goto 97
        call wsnos(nb,ilhkv,hlv,hev,ierr)
        if( ierr .ne. 0 ) goto 97

c---------------------------------------------------------------
c initializing arrays and parameters
c---------------------------------------------------------------

	tunit = 3600.		!time unit (3600=hours)

	nread=0
	idt = 0
	itold = 0
	ntres = 3
	if( ntres .gt. ntrdim ) stop 'error stop: ntrdim'
	threshold(1) = 100
	threshold(2) = 200
	threshold(3) = 500

	do nt=1,ntres
          do k=1,nkn
            do l=1,nlv
	      count(l,k,nt) = 0.
	    end do
	  end do
	end do
	
c---------------------------------------------------------------
c loop on input records
c---------------------------------------------------------------

  300   continue

	call rdnos(nin,it,ivar1,nlvdim,ilhkv,cv3,ierr)

        if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
        if(ierr.ne.0) goto 100

	nread=nread+1
	write(6,*) 'time : ',it,'   record = ',nread
	if( nread .eq. 2 ) idt = it - itold
	if( idt .gt. 0 .and. idt .ne. it-itold ) then
	  write(6,*) it,itold,idt
	  stop 'error stop: idt not uniform'
	end if

	do nt=1,ntres
	  thres = threshold(nt)
          do k=1,nkn
            do l=1,nlv
	      if( cv3(l,k) .ge. thres ) then
	        count(l,k,nt) = count(l,k,nt) + 1. 
	      end if
	    end do
	  end do
	end do

	itold = it

	goto 300

  100	continue

c---------------------------------------------------------------
c end of loop
c---------------------------------------------------------------

	do nt=1,ntres
          do k=1,nkn
            do l=1,nlv
	      count(l,k,nt) = count(l,k,nt) * idt / tunit
	    end do
	  end do
	end do

	do nt=1,ntres
	  it = nint(threshold(nt))
          call wrnos(nb,it,335,nlvdim,ilhkv,count(1,1,nt),ierr)
          if( ierr .ne. 0 ) goto 99
	end do

	write(6,*)
	write(6,*) nread,' records read'
	write(6,*) ' time step = ',idt
	write(6,*)
	write(6,*) 'data written to file elab.nos'
	write(6,*)

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

	stop
   97   continue
        write(6,*) 'error writing header'
        stop 'error stop nosaver'
   98   continue
        write(6,*) 'error opening file'
        stop 'error stop nosaver'
   99   continue
        write(6,*) 'error writing data'
        stop 'error stop nosaver'
	end

c***************************************************************

	function open_file(text,ext)

	implicit none

	integer open_file
	character*(*) text
	character*(*) ext

	integer nin
	character*80 name,file
	integer ifileo

	write(6,*) text
	read(5,'(a)',end=100) name
	if( name .eq. ' ' ) goto 100

	call mkname(' ',name,ext,file)

	nin = ifileo(0,file,'unform','old')
	if( nin .le. 0 ) then
	  write(6,*) 'Cannot open file: ',file
	  stop 'error stop open_file: opening file'
	end if

	open_file = nin

	return
  100	continue
	open_file = 0
	return
	end

c***************************************************************

