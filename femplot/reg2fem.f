
	program reg2fem

	implicit none

	integer ndim
	parameter (ndim=2000)

	logical bdata
	integer iunit,it,nvar,nx,ny,n,i
	integer nrec,ifreq,ivar,ibas,nmax
	integer ndata
	real flag
	real x0,y0,dx,dy
	real data(ndim)
	real newdata(2*ndim)
	real hev(4*ndim)
	real xgv(2*ndim), ygv(2*ndim)

	ifreq = 12
	nmax = 300
	ivar = 2
	ibas = 0
	x0 = 11.
	y0 = 38.
	dx = 0.5
	dy = 0.5
	flag = -999.

	do i=1,4*ndim
	  hev(i) = 0.
	end do

	nvar = 0
	nx = 0
	ny = 0
	nrec = 0
	iunit = 5
	write(6,*) 'reading from STDIN'

    1	continue
	  call rgf_read(iunit,ndim,it,nvar,nx,ny,data,bdata)
	  if( .not. bdata ) goto 2
	  !nx = 21
	  !ny = 19
	  !nx = 5
	  !ny = 7

	  nrec = nrec + 1
	  write(6,*) nrec,it,nvar,nx,ny
	  if( nrec .gt. nmax ) goto 2

	  if( ifreq .eq. 0 .or. mod(nrec,ifreq) .eq. 0 ) then
	    if( ibas .eq. 0 ) then
	      write(6,*) 'creating basin'
	      call create_basin(nx,ny,x0,y0,dx,dy,xgv,ygv)
	      ibas = 1
	    end if
	    write(6,*) 'creating data for ',it
	    n = nx*ny
	    call create_data(nx,ny,flag,data(1+(ivar-1)*n),newdata,ndata)
	    write(6,*) ndata,' have been created'
	    call write_data(it,nx,ny,newdata,hev)

	    write(67,*) it,ndata
	    do i=1,ndata
	      write(67,*) i,xgv(i),ygv(i),newdata(i)
	    end do
	  end if

	  goto 1
    2	continue

	write(6,*) nrec,' records read'

	end

c*****************************************************************

	subroutine create_basin(nx,ny,x0,y0,dx,dy,xgv,ygv)

	implicit none

	integer nx,ny
	real x0,y0,dx,dy
	real xgv(1),ygv(1)

	integer n,i,j
	integer ie,nl,nc,nu
	real x,y

	open(1,file='regular.grd',form='formatted',status='unknown')

	write(1,*) 
	write(1,'(a,2i8)') '0 regular grid ',nx,ny
	write(1,*) 

	n = 0
	do j=0,ny-1

	  y = y0 + j*dy
	  do i=0,nx-1
	    x = x0 + i*dx
	    n = n + 1
	    write(1,1000) 1,n,0,x,y
	    xgv(n) = x
	    ygv(n) = y
	  end do

	  if( j .lt. ny-1 ) then
	    y = y0 + j*dy + dy/2.
	    do i=1,nx-1
	      x = x0 + i*dx - dx/2.
	      n = n + 1
	      write(1,1000) 1,n,0,x,y
	      xgv(n) = x
	      ygv(n) = y
	    end do
	  end if

	end do

	write(1,*) 

	ie = 0
	do j=1,ny-1
	  do i=1,nx-1
	    nl = i + (j-1)*(2*nx-1)
	    nc = nl + nx
	    nu = nc + nx
	    write(1,2000) 2,ie+1,0,3,nl,nl+1,nc
	    write(1,2000) 2,ie+2,0,3,nl+1,nu,nc
	    write(1,2000) 2,ie+3,0,3,nu,nu-1,nc
	    write(1,2000) 2,ie+4,0,3,nu-1,nl,nc
	    ie = ie + 4
	  end do
	end do

	write(1,*) 

	close(1)

	return
 1000	format(i1,i10,i5,2f12.4)
 2000	format(i1,i10,i5,i5,3i10)
	end

c*****************************************************************

	subroutine create_data(nx,ny,flag,data,newdata,ndata)

	implicit none

	integer nx,ny
	real flag
	real data(nx,ny)
	real newdata(1)
	integer ndata

	integer n,i,j,k
	integer ie,nl,nc,nu
	real val,valtot
	integer ip(4)

	ip(1) = 0
	ip(2) = 1
	ip(3) = 2*nx-1
	ip(4) = 2*nx

	n = 0
	do j=1,ny
	  do i=1,nx
	    n = n + 1
	    newdata(n) = data(i,j)
	    !newdata(n) = n
	  end do
	  n = n + nx - 1
	end do

	do j=1,ny-1
	  do i=1,nx-1
	    nl = i + (j-1)*(2*nx-1)
	    nc = nl + nx

	    n = 0
	    valtot = 0.
	    do k=1,4
	      val = newdata(nl+ip(k))
	      if( val .ne. flag ) then
		valtot = valtot + val
		!valtot = valtot + nl+ip(k)
		n = n + 1
	      end if
	    end do

	    if( n .gt. 0 ) then
	      newdata(nc) = valtot / n
	    else
	      newdata(nc) = flag
	    end if

	  end do
	end do

	ndata = nx*ny + (nx-1)*(ny-1)

	end

c*****************************************************************

	subroutine write_data(it,nx,ny,newdata,hev)

	implicit none

	integer it,nx,ny
	real newdata(1)
	real hev(1)

	integer nvers,i
	integer nkn,nel,nlv,nvar
	integer ivar,nlvdim,ierr
	integer ilhkv(1)
	real hlv(1)
	character*80 title

	integer iunit
	save iunit
	data iunit / 0 /

	nvers = 3
	nkn = nx*ny + (nx-1)*(ny-1)
	nel = 4 * (nx-1)*(ny-1)
	nlv = 1
	nvar = 1
	title = 'data from regular grid'

	write(6,*) 'regular data: ',it,nx,ny,nkn,nel
	write(66,*) 'regular data: ',it,nx,ny,nkn,nel
	do i=1,nkn+5
	  write(66,*) i,newdata(i)
	end do
	
	if( iunit .eq. 0 ) then
	  iunit = 1
	  open(1,file='regdata.nos',form='unformatted',status='unknown')

          call whnos        (iunit,nvers
     +                          ,nkn,nel,nlv,nvar
     +                          ,ilhkv,hlv,hev
     +                          ,title
     +                          )
	end if

	ivar = 999
	nlvdim = 1
	call wrnos(iunit,it,ivar,nlvdim,ilhkv,newdata,ierr)

	end

c*****************************************************************

