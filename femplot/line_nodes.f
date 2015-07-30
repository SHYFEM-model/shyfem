c
c $Id: line_nodes.f,v 1.2 2010-03-11 15:32:16 georg Exp $
c
c routines to find nodes close to line
c
c revision log :
c
c 14.09.2009    ggu     routines written from scratch
c
c************************************************************************
c
c given a line defined by 2 or more points (coordinates x,y)
c the program returns the nodes that closest aproximate the line
c
c main routine to call: find_line_nodes
c
c needs basin information and ilinkv(1),lenkv(1),linkv(1) data structure
c to be link with links.f sublnku.f
c
c****************************************************************

	subroutine find_line_nodes(nl,x,y,ndim,n,nodes)

c find nodes close to line
c
c given a line defined by 2 or more points (coordinates x,y)
c the program returns the nodes that closest aproximate the line

	implicit none

	integer nl		!total number of coordinates defining line
	real x(nl), y(nl)	!coordinates of line
	integer ndim		!dimension of nodes()
	integer n		!total number of nodes found (return)
	integer nodes(ndim)	!nodes found (return)

	integer nt,ns,i

	nt = 1

	do i=2,nl
	  nt = nt - 1	!get rid of last node - otherwise we have it twice
	  call find_segment_nodes(x(i-1),y(i-1),x(i),y(i)
     +					,ndim,ns,nodes(nt+1))
	  nt = nt + ns	!add new nodes found
	end do

	n = nt

	end

c****************************************************************

	subroutine find_segment_nodes(x1,y1,x2,y2,ndim,n,nodes)

c find nodes close to one segment of line

	implicit none

	real x1,y1		!coordinates of initial point
	real x2,y2		!coordinates of final point
	integer ndim		!dimension of nodes()
	integer n		!total number of nodes found (return)
	integer nodes(ndim)	!nodes found (return)

	integer k1,k2,ka,kn

	call get_closest_node(x1,y1,k1)
	call get_closest_node(x2,y2,k2)

	call write_node_info(k1)
	call write_node_info(k2)

	ka = k1
	n = 1
	nodes(n) = ka

	do while( ka .ne. k2 )

	  call get_next_node(k1,k2,ka,kn)

	  n = n + 1
	  if( n .gt. ndim ) goto 99
	  nodes(n) = kn
	  ka = kn

	end do

	write(6,*) 'find_segment_nodes finish: ',n

	return
   99	continue
	stop 'error stop find_segment_nodes: ndim'
	end

c****************************************************************

	subroutine get_next_node(k1,k2,ka,kn)

c gets next node close to line segment starting from ka

	implicit none

	integer k1,k2		!start/end node of line segment
	integer ka		!last node found close to line segment
	integer kn		!next node found close to line segment (return)

	integer ndim
	parameter(ndim=100)	!must be at least maximum grade
	integer nodes(ndim)
	integer n,i,k
	real dist,distc,t,t0

	call get_nodes_around(ka,ndim,n,nodes)
	
	call get_distance_from_line(k1,k2,ka,dist,t0)
	kn = 0
	distc = 0.

	do i=1,n
	  k = nodes(i)
	  call get_distance_from_line(k1,k2,k,dist,t)
	  if( t .gt. t0 ) then				!only nodes ahead
	    if( kn .eq. 0 .or. dist .lt. distc ) then   !choose closest node
		kn = k
		distc = dist
	    end if
	  end if
	end do

	if( kn .eq. 0 ) then
	  stop 'error stop get_next_node: no node found'
	end if

	end

c****************************************************************

	subroutine get_distance_from_line(k1,k2,kp,dist,t)

c computes distance of kp from line given by k1,k2

	use basin

	implicit none

	include 'param.h'

	integer k1,k2		!start/end node of line segment
	integer kp		!node to compute distance to line segment
	real dist		!distance of kp from line segment (return)
	real t			!position of intersection on segment (return)

c (xs,ys) is intersection point on segment
c (xa,ya) is vector of segment : (xa,ya) = (x2,y2) - (x1,y1)
c (xs,ys) = (x1,y1) + t * (xa,ya)

	real x1,y1,x2,y2,xp,yp
	real xa,ya,xs,ys


	x1 = xgv(k1)
	y1 = ygv(k1)
	x2 = xgv(k2)
	y2 = ygv(k2)
	xp = xgv(kp)
	yp = ygv(kp)

	xa = x2 - x1
	ya = y2 - y1

	t = ( xa*(xp-x1) + ya*(yp-y1) ) / ( xa*xa + ya*ya )

	xs = x1 + t * xa
	ys = y1 + t * ya

	dist = sqrt( (xp-xs)**2 + (yp-ys)**2 )

	end

c****************************************************************
c****************************************************************
c****************************************************************

	subroutine get_closest_node(x0,y0,kc)

c finds closest node to coordinate (x0,y0)

	use basin

	implicit none

	include 'param.h'

	real x0,y0		!coordinates of point
	integer kc		!closest node to point (return)


	integer k
	real dist,distc

	kc = 1
	distc = (xgv(1)-x0)**2 + (ygv(1)-y0)**2

	do k=2,nkn
	  dist = (xgv(k)-x0)**2 + (ygv(k)-y0)**2
	  if( dist .lt. distc ) then
	    distc = dist
	    kc = k
	  end if
	end do

	end

c****************************************************************

	subroutine write_node_info(k)

c writes node info for node k

	use basin

	implicit none

	include 'param.h'

	integer k


	write(6,*) 'node = ',k,xgv(k),ygv(k)

	end

c****************************************************************
c****************************************************************
c****************************************************************
c****************************************************************
c****************************************************************

	subroutine line_points_test

	use mod_geom
	use basin

	implicit none

	include 'param.h'

	integer ndim
	parameter (ndim=500)

	character*80 basfil
	integer nodes(ndim)
	integer n,nb,nl
	integer nlkdi
	real x(ndim), y(ndim)

	integer iapini

	call shyfem_copyright('line_nodes - find node along line')

	call read_line(ndim,nl,x,y)

	if( iapini(1,0,0,0) .le. 0 ) stop

	call mod_geom_init(nkn,nel,ngr)

	nlkdi = 3*nel + 2*nkn
        call mklenk(nlkdi,nkn,nel,nen3v,ilinkv,lenkv)
        call mklink(nkn,ilinkv,lenkv,linkv)

	call find_line_nodes(nl,x,y,ndim,n,nodes)

	write(6,*) n,' nodes found in line'
	call write_line(n,nodes)
	write(6,*) 'output written to file 66 and 77'

	end

c****************************************************************

	subroutine read_memory_basin(basin)

	implicit none

	character*(*) basin

	integer i

	basin = ' '

	open(1,file='.memory',err=99)
	read(1,*) 
	read(1,'(a)') basin 
	close(1)

	do i=len(basin),1,-1
	  if( basin(i:i) .ne. ' ' ) goto 1
	end do
    1	continue
	basin(i+1:) = '.bas'

	write(6,*) 'using basin file: ',basin(1:i+4)
	return
   99	continue
	write(6,*) 'Cannot read memory file .memory'
	stop 'error stop read_memory_basin: memory'
	end

c****************************************************************

	subroutine write_line(n,nodes)

	use basin

	implicit none

	include 'param.h'

	integer n
	integer nodes(n)

	integer i,k,kext
	real x,y

	integer ipext

c--------------------------------------------------------------
c write grd file
c--------------------------------------------------------------

        open(66,file='line_nodes.grd',status='unknown',form='formatted')

	write(66,*)
	do i=1,n
	  k = nodes(i)
	  kext = ipext(k)
	  x = xgv(k)
	  y = ygv(k)
	  write(66,1000) 1,kext,0,x,y
	  nodes(i) = kext
	end do

	write(66,*)
	write(66,2000) 3,k,0,n
	write(66,3000) (nodes(i),i=1,n)
	write(66,*)

        close(66)

        write(6,*) 'grd file written to line_nodes.grd'

c--------------------------------------------------------------
c write txt file
c--------------------------------------------------------------

        open(77,file='line_nodes.txt',status='unknown',form='formatted')

	do i=1,n
	  kext = nodes(i)
	  write(77,*) kext
	end do
	write(77,*) 0

        close(77)

        write(6,*) 'txt file written to line_nodes.txt'

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	return
 1000	format(i1,2i10,2f14.4)
 2000	format(i1,3i10)
 3000	format((5i10))
	end

c****************************************************************

	subroutine read_line(ndim,nl,x,y)

	implicit none

	integer ndim,nl
	real x(ndim), y(ndim)

	integer i,k,ityp
	real xp,yp

	write(6,*) 'reading line from STDIN'

	nl = 0
    1	continue
	  !read(5,*,end=2) k,ityp,xp,yp
	  read(5,*,end=2) xp,yp,ityp
	  nl = nl + 1
	  if( nl .gt. ndim ) stop 'error stop read_line: ndim'
	  x(nl) = xp
	  y(nl) = yp
	  goto 1
    2	continue

	write(6,*) 'points read: ',nl
	do i=1,nl
	  write(6,*) x(i),y(i)
	end do

	end

c****************************************************************

	program line_points_main
	implicit none
	call line_points_test
	end

c****************************************************************

