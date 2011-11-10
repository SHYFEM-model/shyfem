c
c $Id: ousutil.f,v 1.3 2009-09-14 08:20:58 georg Exp $
c
c general utilities for output files
c
c revision log :
c
c 03.06.2011	ggu	copied from ousutil.f
c 07.11.2011	ggu	layer thickness for hybrid coordinates
c
c******************************************************************

        subroutine get_records_from_stdin(ndim,irec,ball)

c gets records to extract from stdin

        implicit none

        integer ndim
        integer irec(ndim)
	logical ball

        integer i,ir

	ball = .false.
        do i=1,ndim
          irec(i) = 0
        end do

        write(6,*) 'Please enter the record numbers to be extracted.'
        write(6,*) 'Enter every record number on a single line.'
        write(6,*) 'Enter -1 if you want all records.'
        write(6,*) 'Finish with 0 on the last line.'
        write(6,*) 'example:'
        write(6,*) '   5'
        write(6,*) '  10'
        write(6,*) '  15'
        write(6,*) '  0'
        write(6,*) ' '

        do while(.true.)
          write(6,*) 'Enter record to extract (0 to end): '
          ir = 0
          read(5,'(i10)') ir

          if( ir .le. 0 ) then
	    if( ir .eq. -1 ) ball = .true.
            return
          else if( ir .gt. ndim ) then
            write(6,*) 'Cannot extract records higher than ',ndim
            write(6,*) 'Please change ndim and recompile.'
            stop 'error stop get_records_from_stdin: ndim'
          else
            irec(ir) = 1
          end if
        end do

        end

c******************************************************************

        subroutine get_nodes_from_stdin(ndim,nnodes,nodes,nodese)

c gets records to extract from stdin

        implicit none

        integer ndim            !dimension of nodes
        integer nnodes          !total number of nodes read
        integer nodes(ndim)     !array with node numbers (nnodes in total)
        integer nodese(ndim)    !array with external node numbers

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

c******************************************************************

	subroutine level_e2k(nkn,nel,nen3v,ilhv,ilhkv)

c computes ilhkv from ilhv

	implicit none

	integer nkn,nel
	integer nen3v(3,1)
	integer ilhv(1)
	integer ilhkv(1)

	integer k,ie,ii,lmax

	do k=1,nkn
	  ilhkv(k) = 1
	end do

	do ie=1,nel
	  lmax = ilhv(ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    if( ilhkv(k) .lt. lmax ) ilhkv(k) = lmax
	  end do
	end do

	end

c******************************************************************

        subroutine write_connections(nkn,nel,nen3v,xgv,ygv)

c writes FEM connectivity
c
c there are nkn nodes and nel elements
c elements are always made out of 3 nodes
c
c legend:
c
c nkn                   total number of nodes
c nel                   total number of elements
c k                     node number [1...nkn]
c ie                    element number [1...nel]
c xgv(k),ygv(k)         coordinates of node k
c (nen3v(ii,ie),ii=1,3) element connectivity - nodes of element ie

        implicit none

        integer nkn,nel
        integer nen3v(3,1)
        real xgv(1), ygv(1)

        integer k,ie,ii

        open(1,file='connections.gis')

        write(1,*) nkn                          !total number of nodes
        do k=1,nkn
          write(1,*) k,xgv(k),ygv(k)            !node number, x, y
        end do

        write(1,*) nel                          !total number of elements
        do ie=1,nel
          write(1,*) ie,(nen3v(ii,ie),ii=1,3)   !element number, n1, n2, n3
        end do

        close(1)

        end

c******************************************************************

        subroutine make_name_with_time(pre,post,it,name)

        implicit none

        character*(*) pre,post,name
        integer it

        integer i
        character*10 naux

        write(naux,'(i10)') it

        do i=1,10
          if( naux(i:i) .eq. ' ' ) naux(i:i) = '_'
        end do

        name = pre // naux // post

        end

c******************************************************************
c******************************************************************
c******************************************************************

	subroutine get_sigma_info(nsigma,hsigma)

	implicit none

	integer nsigma
	real hsigma

        integer nlvdi,nlv
        common /level/ nlvdi,nlv

	real hlv(1)
	common /hlv/hlv

	integer l

c---------------------------------------------------------
c scan depth structure
c---------------------------------------------------------

	hsigma = 10000.

	do l=2,nlv
	  if( hlv(l) .gt. hlv(l-1) ) goto 1
	end do

c---------------------------------------------------------
c only sigma layers found
c---------------------------------------------------------

	if( hlv(nlv) .ne. -1 ) then
	  stop 'error stop get_sigma_info: internal error (1)'
	end if
	nsigma = nlv
	return

c---------------------------------------------------------
c zeta or hybrid levels found
c
c this algorithm cannot handle hybrid levels with only 2 sigma layers
c---------------------------------------------------------

    1	continue
	if( l .eq. 2 ) then	!only zeta levels
	  nsigma = 0
	else			!hybrid levels
	  nsigma = l
	  hsigma = hlv(l)
	end if

c---------------------------------------------------------
c end of routine
c---------------------------------------------------------

	end

c******************************************************************

	subroutine get_layer_thickness(ie,lmax,bzeta,nsigma,hsigma,hl)

c returns layer thickness - works also for lmax higher than actual layers
c
c works also for lmax higher than actual layers
c in this case the last values for hl are 0

	implicit none

	integer ie
	integer lmax
	logical bzeta
	integer nsigma
	real hsigma
	real hl(1)

	real hlv(1)
	common /hlv/hlv
	real hev(1)
	common /hev/hev
	real zenv(3,1)
	common /zenv/zenv

	integer ii,l
	real zmed
	real htot,hsig,htop,hbot

c---------------------------------------------------------
c compute contribution of zeta
c---------------------------------------------------------

	zmed = 0.
	if( bzeta ) then
          zmed = 0.
          do ii=1,3
            zmed = zmed + zenv(ii,ie)
          end do
          zmed = zmed / 3.
	end if

c---------------------------------------------------------
c compute level structure of sigma levels
c---------------------------------------------------------

	htot = hev(ie)
	hsig = min(htot,hsigma) + zmed

	hbot = 0.
	do l=1,nsigma
	  htop = hbot
	  hbot = hlv(l)
	  if( l .eq. nsigma ) hbot = -1.
	  hl(l) = -hsig * (hbot-htop)
	end do

c---------------------------------------------------------
c compute level structure of zeta and/or hybrid levels
c---------------------------------------------------------

	if( lmax .gt. nsigma ) then		!also zeta coordinates
	  if( lmax .eq. 1 ) then		!just one layer
	    hl(1) = htot + zmed
	  else
	    hbot = hsigma
	    if( nsigma .eq. 0 ) hbot = -zmed
	    do l=nsigma+1,lmax
	      htop = hbot
	      hbot = hlv(l)
	      if( hbot .gt. htot ) hbot = htot	!last layer
	      hl(l) = hbot - htop
	    end do
	  end if
	end if

c---------------------------------------------------------
c end of routine
c---------------------------------------------------------

	end

c******************************************************************

