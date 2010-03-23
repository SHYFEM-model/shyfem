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

	program nosextr_nodes

c extracts single nodes from nos file -> creates time series
c
c interactive version

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
	real cv1(nlvdim,nkndim)
	real cv2(nlvdim,nkndim)
	real cv3(nlvdim,nkndim)

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

c---------------------------------------------------------------
c open simulation and basin
c---------------------------------------------------------------

	!if(iapini(3,nkndim,neldim,0).eq.0) then
	!	stop 'error stop : iapini'
	!end if

c---------------------------------------------------------------
c open NOS file and read header
c---------------------------------------------------------------

	nin1 = open_file('Enter first file: ','.nos')
	if(nin1.le.0) goto 100
	nin2 = open_file('Enter second file: ','.nos')
	if(nin2.le.0) goto 100

        nvers=3
	call rfnos(nin1,nvers,nkn,nel,nlv,nvar,title,ierr)
        if(ierr.ne.0) goto 100
        call dimnos(nin1,nkndim,neldim,nlvdim)

	call rfnos(nin2,nvers,nkn,nel,nlv,nvar,title,ierr)
        if(ierr.ne.0) goto 100
        call dimnos(nin1,nkndim,neldim,nlvdim)

        write(6,*) 'nvers    : ',nvers
        write(6,*) 'nkn,nel  : ',nkn,nel
        write(6,*) 'nlv,nvar : ',nlv,nvar
        write(6,*) 'title    : ',title

	call rsnos(nin1,ilhkv,hlv,hev,ierr)
        if(ierr.ne.0) goto 100
	call rsnos(nin2,ilhkv,hlv,hev,ierr)
        if(ierr.ne.0) goto 100

	write(6,*) 'Available levels: ',nlv
	write(6,*) (hlv(l),l=1,nlv)

c---------------------------------------------------------------
c initializing units and nodes to be extracted
c---------------------------------------------------------------

	nread=0

        call mkname(' ','elab','.nos',file)
        write(6,*) 'writing file ',file(1:50)
        nb = ifileo(0,file,'unform','new')
        if( nb .le. 0 ) goto 98

        call wfnos(nb,3,nkn,nel,nlv,1,title,ierr)
        if( ierr .ne. 0 ) goto 97
        call wsnos(nb,ilhkv,hlv,hev,ierr)
        if( ierr .ne. 0 ) goto 97

c---------------------------------------------------------------
c loop on input records
c---------------------------------------------------------------

  300   continue

	call rdnos(nin1,it1,ivar1,nlvdim,ilhkv,cv1,ierr)

        if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
        if(ierr.ne.0) goto 100

	call rdnos(nin2,it2,ivar2,nlvdim,ilhkv,cv2,ierr)

        if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
        if(ierr.ne.0) goto 100

	if( it1 .ne. it2 ) then
	  write(6,*) 'time is different: ',it1,it2
	  stop 'error stop: time'
	end if

	it = it1

	nread=nread+1
	write(6,*) 'time : ',it,ivar1,ivar2

c	---------------------------------------------------------
c	write to file
c	---------------------------------------------------------

	soglia = 0.1

        do l=1,nlv
          do k=1,nkn
            c1 = cv1(l,k)
            c2 = cv2(l,k)

	    c3 = 0.
	    if( c1 .gt. soglia .and. c2 .gt. soglia ) then
	      c3 = c2 / c1
	    end if
	    !if( c1 .ne. 0. ) c3 = c2 / c1

            cv3(l,k) = c3
          end do
        end do

        call wrnos(nb,it,333,nlvdim,ilhkv,cv3,ierr)    !aver
        if( ierr .ne. 0 ) goto 99

	goto 300

  100	continue

c---------------------------------------------------------------
c end of loop
c---------------------------------------------------------------

	write(6,*)
	write(6,*) nread,' records read'
	write(6,*)
	write(6,*) 'data written to file elab.nos'
	write(6,*)

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

	return
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

        subroutine get_nodes_from_stdin(ndim,nnodes,nodes,nodese)

c gets records to extract from stdin

        implicit none

        integer ndim		!dimension of nodes
        integer nnodes		!total number of nodes read
        integer nodes(ndim)	!array with node numbers (nnodes in total)
        integer nodese(ndim)	!array with external node numbers

        integer ir
	integer ipint

	nnodes = 0

        write(6,*) 'Please enter the node numbers to be extracted.'
        write(6,*) 'Enter every node on a single line.'
        write(6,*) 'Finish with 0 on the last line.'
        write(6,*) 'example:'
        write(6,*) '  5'
        write(6,*) '  100'
        write(6,*) '  1505'
        write(6,*) '  0'
        write(6,*) ' '

        do while(.true.)
          write(6,*) 'Enter node to extract (0 to end): '
          ir = 0
          read(5,'(i10)') ir

          if( ir .le. 0 ) return

	  nnodes = nnodes + 1

          if( nnodes .gt. ndim ) then
            write(6,*) 'Cannot extract more than ',ndim,' nodes'
	    stop 'error stop get_nodes_from_stdin: ndim'
          else
            nodese(nnodes) = ir
            nodes(nnodes) = ipint(ir)
	    if( nodes(nnodes) .le. 0 ) then
	      write(6,*) 'No such node ',ir,' ... ignoring'
	      nnodes = nnodes - 1
	    end if
          end if
        end do

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

