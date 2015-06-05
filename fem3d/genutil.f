c
c $Id: ousutil.f,v 1.3 2009-09-14 08:20:58 georg Exp $
c
c general utilities for output files
c
c revision log :
c
c 03.06.2011	ggu	copied from ousutil.f
c 07.11.2011	ggu	layer thickness for hybrid coordinates
c 14.11.2011	ggu	sigma layer routines copied to sigmautil.f
c 25.01.2013	ggu	new routine level_k2e()
c 05.09.2013	ggu	new routine compute_levels_on_element()
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

          write(6,'(i10)') ir

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

	subroutine level_k2e(nkn,nel,nen3v,ilhkv,ilhv)

c computes ilhv from ilhkv - is not exact (need hev to do better)

	implicit none

	integer nkn,nel
	integer nen3v(3,1)
	integer ilhkv(1)
	integer ilhv(1)

	integer k,ie,ii,lmax,lmin

	lmax = 0
	do k=1,nkn
	  lmax = max(lmax,ilhkv(k))
	end do

	do ie=1,nel
	  lmin = lmax
	  do ii=1,3
	    k = nen3v(ii,ie)
	    lmin = min(lmin,ilhkv(k))
	  end do
	  ilhv(ie) = lmin
	end do

	end

c******************************************************************

	subroutine compute_levels_on_element(ie,zenv,zeta)

	implicit none

	integer ie
	real zenv(3,1)
	real zeta

	integer ii
	real z

	z = 0.
	do ii=1,3
	  z = z + zenv(ii,ie)
	end do

	zeta = z / 3.

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
        character*20 line

	logical dts_has_date

	if( dts_has_date() ) then
	  call dtsgf(it,line)
          name = pre // line // post
	else
          write(line,'(i10)') it
          do i=1,10
            if( line(i:i) .eq. ' ' ) line(i:i) = '_'
          end do
          name = pre // line(1:10) // post
	end if

        end

c******************************************************************
c******************************************************************
c******************************************************************

