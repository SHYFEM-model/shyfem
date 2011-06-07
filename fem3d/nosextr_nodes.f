c
c $Id: nosextr_nodes.f,v 1.3 2009-04-07 10:43:57 georg Exp $
c
c extract nodes from NOS file
c
c revision log :
c
c 18.11.1998    ggu     check dimensions with dimnos
c 24.02.1999    ggu     use n2int for node number translation
c 03.12.2001    ggu     cleaned up, hakata bay
c 03.06.2011    ggu     routine adjourned
c
c****************************************************************

	program nosextr_nodes

c extracts single nodes from nos file -> creates time series

	implicit none

	include 'param.h'

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

	integer ilhkv(nkndim)
	real hlv(nlvdim)
	real hev(neldim)

	logical berror
	integer i,n,k,ke,l,lmax
	integer nread,nunit
	integer nvers
	integer nlv,nvar,ivar,ierr
	integer nin,it

	integer iapini,ideffi,ialfa

c--------------------------------------------------
	integer nudim
	parameter(nudim=300)
	integer iunit(nudim)
	character*50 file

c---------------------------------------------------------------
c nodes for extraction
c---------------------------------------------------------------

	integer ndim
	integer nnodes
	parameter( ndim = nkndim )
	integer nodes(ndim)	!node numbers
	integer nodese(ndim)	!external node numbers

c---------------------------------------------------------------
c open simulation and basin
c---------------------------------------------------------------

	if(iapini(3,nkndim,neldim,0).eq.0) then
		stop 'error stop : iapini'
	end if

c---------------------------------------------------------------
c open NOS file and read header
c---------------------------------------------------------------

	nin=ideffi('datdir','runnam','.nos','unform','old')
	if(nin.le.0) goto 100

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

c---------------------------------------------------------------
c initializing units and nodes to be extracted
c---------------------------------------------------------------

	nread=0

	nunit = 60
	do i=1,nudim
	  iunit(i) = 0
	end do

        call get_nodes_from_stdin(ndim,nnodes,nodes,nodese)

	if( nnodes .le. 0 ) goto 100

c---------------------------------------------------------------
c loop on input records
c---------------------------------------------------------------

  300   continue

	call rdnos(nin,it,ivar,nlvdim,ilhkv,cv3,ierr)

        if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
        if(ierr.ne.0) goto 100

	nread=nread+1
	write(6,*) 'time : ',it,ivar

	if( ivar .gt. nudim ) then              !ivar too high
          write(6,*) 'nudim is too low ',ivar,nudim
          stop 'error stop: nudim'
        else if( iunit(ivar) .eq. 0 ) then      !not yet initialized
	  nunit = nunit + 1
	  file = ' '
	  n = ialfa(float(ivar),file,-1,-1)
	  file(n+1:) = '.dat'
	  write(6,*) 'opening file ',file
	  open(nunit,file=file,status='unknown',form='formatted')
	  iunit(ivar) = nunit
	else                                    !already initialized
	  nunit = iunit(ivar)
	end if

c	---------------------------------------------------------
c	write to file
c	---------------------------------------------------------

        write(nunit,'(i10,30e12.4)') it,(cv3(1,nodes(i)),i=1,nnodes)

	do i=1,nnodes
	  k = nodes(i)
	  ke = nodese(i)
	  lmax = ilhkv(k)
	  write(4,*) it,i,ke,k,lmax,ivar
	  write(4,*) (cv3(l,k),l=1,lmax)
	  write(3,*) it,i,ke,k,lmax,ivar
	  write(3,'((6f10.2))') (cv3(l,k),l=1,lmax)
	end do

	goto 300

  100	continue

c---------------------------------------------------------------
c end of loop
c---------------------------------------------------------------

	write(6,*)
	write(6,*) nread,' records read'
	write(6,*)
	write(6,*) 'data written to file 3 and 4'
	write(6,*)

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

	end

c***************************************************************

