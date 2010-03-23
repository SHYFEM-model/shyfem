c
c smoothing program
c
c contents :
c
c revision log :
c 
c 01.02.2002    ggu     handle periodic line
c 06.08.2002    ggu     read also depth value
c 23.09.2004    ggu     adapted for malta, bug fix
c 19.10.2005    ggu     documentation and description
c
c description :
c
c This program smothes and reduces the number of points on a line
c
c The program ask for the GRD file that contains the lines
c It also asks for the smoothing parameter (sigma) and reduction
c	parameter (reduct)
c The units of these values must be in the units of the line coordinates.
c 
c Smoothing:
c
c The program uses a gaussian curve as the smoothing operator, with
c its standard deviation given by sigma. Therefore, if (x,y) coordinates
c of the lines are in meters and sigma=100, the lengthscale for the 
c smoothing is about 100 meters, which means that points closer than 100
c meters to the point will have much more effect than points further away
c than 100 m.
c
c Reduction:
c
c The program uses the parameter reduct to decide if a point should
c be kept or not. If reduct=100, than all points closer than 100 m
c to the reference point are discarded. At the end a line is left with
c a distance between the points of minimum 100 m.
c
c Adaption of resolution:
c
c If the points of the line have no depth values, then the 
c smoothing and reduction is working as described above. If, however,
c depth values do exists, then they are used to change the resolution
c of the line points. For a value of 2, twice the points in the reduction
c algorithm are retained. In other words, the value for reduct is changed
c locally to a value of 100/2=50, and only points closer as 50 are
c eliminated. Line points without depth receive an implicit value of 1.
c The value -1 indicates that this point should never be eliminated 
c nor should it be smoothed.
c
c Please note that the change in resolution is applicable only for the
c reduction algorithm, and that the smoothing algorithm actually does
c only use -1 for not smoothing.
c
c***********************************************************************

	program grids

c smoothing program
c
c same as gridf but smooth on length of line, not points

	implicit none

	integer nkndim,neldim,nlidim,nlndim,ndim
	parameter( nkndim = 50000 )
	parameter( neldim = 10000 )
	parameter( nlidim = 1000 )
	parameter( nlndim = 50000 )
	parameter( ndim = 10000 )

	integer ner
	integer nco,nkn,nknh,nel,nelh,nli
	logical bstop
	character*80 file

	integer ipv(nkndim)
	integer ipev(neldim)
	integer iarv(neldim)
	integer nen3v(3,neldim)

	real xgv(nkndim)
	real ygv(nkndim)
	real hkv(nkndim)
	real hev(neldim)

	integer iplv(nlidim)
	integer ialrv(nlidim)
	integer ipntlv(0:nlidim)
	integer inodlv(nlndim)

	logical bperiod
	integer nl,nt
	real xt(ndim)
	real yt(ndim)
	real ht(ndim)
	real aux(ndim)

	integer l
	integer nline
	integer nnode
	real sigma
	real reduct

c-----------------------------------
c sigma		smoothing parameter (size of gaussian kernel)
c reduct	reduction of points in line (in meters)
c-----------------------------------
	sigma = 2.0
	sigma = 1.5
	sigma = 500.
	sigma = 200.
	reduct = 4.0
	reduct = 200.
c----------------------------------- hakata
	sigma = 200.
	reduct = 200.
c----------------------------------- circle
	sigma = 0.0
	reduct = 0.1
c----------------------------------- lido
	sigma = 0.
	reduct = 100.
c----------------------------------- curonian lagoon
	sigma = 200.
	reduct = 200.
c----------------------------------- brasile
	sigma = 0.005
	reduct = 0.005
c----------------------------------- malta
	sigma = 0.1
	reduct = 0.2
	sigma = 0.1 * 245.
	reduct = 0.2 * 245.
c----------------------------------- hue
	sigma = 30.
	reduct = 300.
c-----------------------------------

c------------------------------------------------------
c read in parameters and file name
c------------------------------------------------------

	file = ' '
	write(6,*) 'Enter file name with line: '
	read(5,'(a)') file
	if( file .eq. ' ' ) stop 'no file given'

	write(6,*) 'Enter sigma and reduct: '
	read(5,*) sigma,reduct

	write(6,*) 'file   :',file
	write(6,*) 'sigma  :',sigma
	write(6,*) 'reduct :',reduct

c------------------------------------------------------

	nline = 0
	nnode = 0

	ner = 6
	bstop = .false.

        call rdgrd(file,ner,bstop,nco,nkn,nknh,nel,nelh,nli
     +                  ,nkndim,neldim,nlidim
     +                  ,ipv,ipev,iarv,nen3v,xgv,ygv,hev,hkv
     +                  ,nlndim,iplv,ialrv,ipntlv,inodlv)

	write(6,*) 'comments : ',nco
	write(6,*) 'nodes    : ',nkn,nknh
	write(6,*) 'elements : ',nel,nelh
	write(6,*) 'lines    : ',nli

	call insnod(nkn,ipv)
c	call prline(nli,iplv,ialrv,ipntlv,inodlv)
c	call prlixy(nli,iplv,ialrv,ipntlv,inodlv,xgv,ygv)

	open(99,file='smooth.grd',status='unknown',form='formatted')
	open(98,file='reduce.grd',status='unknown',form='formatted')

	do l=1,nli
	  nl = ndim
	  call extrli(l,nli,iplv,ialrv,ipntlv,inodlv,xgv,ygv,hkv
     +				,xt,yt,ht,nl,nt)
	  call prlisxy(l,nli,iplv,ialrv,ipntlv,inodlv,xgv,ygv)
	  call rotate_line(nl,xt,yt,ht,aux)
	  call prlisxy(l,nli,iplv,ialrv,ipntlv,inodlv,xgv,ygv)
	  nline = iplv(l)

	  call mkperiod(xt,yt,nl,bperiod)
	  call intpdep(ht,nl,bperiod)

	  call smooth(sigma,xt,yt,ht,nl,bperiod)
	  call wrline(99,nline,nnode,nl,xt,yt,ht,nt,bperiod)

	  if( bperiod ) then
	    nl = nl + 1
	    ht(nl) = ht(1)
	  end if

	  !call reduce(reduct,xt,yt,ht,nl)
	  call reduce_line(reduct,xt,yt,ht,nl)
	  call wrline(98,nline,nnode,nl,xt,yt,ht,nt,bperiod)
	end do

	write(6,*) 'routine finished...'

	close(99)
	close(98)

	write(6,*) 'files smooth.grd and reduce.grd written'

	end

c**********************************************************

	subroutine distxy(ndim,nl,xt,yt,ht,dxy)

c sets up distance between nodes

c dxy(i) -> distance to next node (i,i+1)
c dxy(i-1) -> distance to previous node (i,i-1)

	implicit none

	integer ndim
	integer nl
	real xt(1)
	real yt(1)
	real ht(1)
	real dxy(-ndim:2*ndim)

	integer i
	real xo,yo,xn,yn
	real dx,dy,dist
	real ho,hn,h

	xn = xt(nl)
	yn = yt(nl)
	hn = ht(nl)

	do i=1,nl
	  xo = xn
	  yo = yn
	  ho = hn
	  xn = xt(i)
	  yn = yt(i)
	  hn = ht(i)

	  dx = xn - xo
	  dy = yn - yo
	  dist = sqrt( dx*dx + dy*dy )

	  dxy(i-1) = dist
	  dxy(i-1+nl) = dist
	  dxy(i-1-nl) = dist
	end do

	dxy(2*nl) = dxy(0)

	end
	  
c********************************************************

	subroutine smooth(sigma,xt,yt,ht,nl,bperiod)

c smoothing

	real sigma
	real xt(1)
	real yt(1)
	real ht(1)
	integer nl
	logical bperiod

	integer ndim
	parameter(ndim=10000)

	integer i
	integer ngk
	real dist,dx,dy,distot
	real gk(-ndim:ndim)
	real raux(-ndim:2*ndim)
	real dxy(-ndim:2*ndim)

	if( nl .gt. ndim ) stop 'error stop smooth: ndim'

        if( sigma .le. 0. ) return

c set up dxy

	call distxy(ndim,nl,xt,yt,ht,dxy)

	distot = 0.
	do i=1,nl-1
	  distot = distot + dxy(i)
	end do

	dx = xt(1) - xt(nl)
	dy = yt(1) - yt(nl)
	dist = sqrt ( dx**2 + dy**2 )
	write(6,*) 'distance : ',nl,dist,distot,dist/distot,bperiod

c create smoothing kernel

c	ngk = ndim
c	call gkernel(sigma,ndim,ngk,gk)
c	write(6,*) 'Gaussian kernel : ',sigma,ngk

c smooth x and y

	call grsmooth(ndim,nl,sigma,xt,raux,dxy,ht,bperiod)
	call grsmooth(ndim,nl,sigma,yt,raux,dxy,ht,bperiod)

	end

c********************************************************

	subroutine grsmooth(ndim,nl,sigma,rt,raux,dxy,ht,bperiod)

c gaussian smoothing

	integer ndim
	integer nl
	real sigma
	real rt(1)
	real raux(-ndim:2*ndim)
	real dxy(-ndim:2*ndim)
	real ht(1)
	logical bperiod

	integer i,j
	real val
	real rkern,tkern
	real rintv

c set up circular array

	do i=1,nl
	  val = rt(i)
	  raux(i) = val
	  raux(i+nl) = val
	  raux(i-nl) = val
	end do

c set up gaussian kernel parameters

        pi = 4.*atan(1.)
        eps = 1.e-7

        a = 1. / ( sqrt(2.*pi) * sigma )
        b = - 1. / ( 2 * sigma * sigma )

c convolution

	do i=1,nl

	  rintv = 0.5 * ( dxy(i-1) + dxy(i) )
	  rkern = a
	  weight = rkern * rintv
	  tkern = weight
	  val = raux(i) * weight
	  jmin = i - nl
	  jmax = i + nl
	  if( .not. bperiod ) then
	    jmin = max(2,jmin)
	    jmax = min(nl-1,jmax)
	  end if
c	  write(6,*) i,i,weight,rkern,tkern,val

	  rkern = 1.
	  dist = 0.
	  j = i
	  do while( rkern .gt. eps .and. j .lt. jmax )
	    j = j + 1
	    dist = dist + dxy(j-1)
	    rintv = 0.5 * ( dxy(j-1) + dxy(j) )
	    rkern = a * exp( b * dist * dist )
	    weight = rkern * rintv
	    tkern = tkern + weight
	    val = val + raux(j) * weight
c	    write(6,*) i,j,weight,rkern,tkern,val,dist
	  end do

c	  write(6,*) i,j-i,nl,dist,rkern

	  rkern = 1.
	  dist = 0.
	  j = i
	  do while( rkern .gt. eps .and. j .gt. jmin )
	    j = j - 1
	    dist = dist + dxy(j)
	    rintv = 0.5 * ( dxy(j-1) + dxy(j) )
	    rkern = a * exp( b * dist * dist )
	    weight = rkern * rintv
	    tkern = tkern + weight
	    val = val + raux(j) * weight
c	    write(6,*) i,j,weight,rkern,tkern,val,dist
	  end do

c	  write(6,*) i,i-j,nl,dist,rkern

c	  write(6,*) '*** ',i,val,rt(i),tkern

	  if( tkern .gt. 0. .and. ht(i) .ge. 0. ) then
	    rt(i) = val / tkern
	  else
	    rt(i) = raux(i)
	  end if

	end do

	end

c********************************************************

	subroutine reduce_line(reduct,xt,yt,ht,nl)

c shell for reduce points

	implicit none

	real reduct
	real xt(1), yt(1)
	real ht(1)
	integer nl

	integer nparts,nf,ipart,nn,i

	integer ndim
	parameter (ndim=10000)

	real xn(ndim), yn(ndim), hn(ndim)
	real xf(ndim), yf(ndim), hf(ndim)

	if( nl .gt. ndim ) stop 'error stop reduce_line: ndim'

	call prepare_parts(nl,ht,nparts)
	write(6,*) 'line consists of parts: ',nparts
	nf = 0

	do ipart=1,nparts
	  call extract_part(ipart,nl,xt,yt,ht,nn,xn,yn,hn)
	  write(6,*) 'extracted part: ',ipart,nl,nn,reduct
	  call reduce_part(nn,xn,yn,hn,reduct)
	  write(6,*) 'reduced part: ',ipart,nl,nn,reduct
	  call add_part(nf,xf,yf,hf,nn,xn,yn,hn)
	  write(6,*) 'added part: ',ipart
	  nf = nf - 1	!do not count last point
	end do

	nf = nf + 1	!consider last point
	if( nf .gt. ndim ) stop 'error stop reduce_line: ndim (2)'

	nl = nf
	do i=1,nf
	  xt(i) = xf(i)
	  yt(i) = yf(i)
	  ht(i) = hf(i)
	end do

	end

c********************************************************

	subroutine rotate_line(nl,xt,yt,ht,aux)

c rotates line to have fixed point as first point

	implicit none

	integer nl
	real xt(1), yt(1)
	real ht(1)
	real aux(1)

	integer i
	logical bperiod, is_periodic

	real h
	logical is_fixed
	is_fixed(h) = h .lt. 0. .and. h .gt. -990.

	bperiod = is_periodic(xt,yt,nl)
	if( .not. bperiod ) return

	do i=1,nl-1
	  if( is_fixed(ht(i)) ) goto 1
	end do
	return	!no fixed point found -> do not rotate

    1	continue

	call rotate(nl-1,i,xt,aux)
	call rotate(nl-1,i,yt,aux)
	call rotate(nl-1,i,ht,aux)

	xt(nl) = xt(1)
	yt(nl) = yt(1)
	ht(nl) = ht(1)

	end

c********************************************************

	subroutine rotate(n,ifix,a,aux)

	implicit none

	integer n,ifix
	real a(1), aux(1)

	integer i,ip

	ip = ifix
	do i=1,n
	  aux(i) = a(ip)
	  ip = mod(ip,n) + 1
	end do

	do i=1,n
	  a(i) = aux(i)
	end do

	end

c********************************************************

	subroutine prepare_parts(nl,ht,nparts)

c computes total number of parts of a line
c
c line must already be rotated

	implicit none

	integer nl
	real ht(1)
	integer nparts

	integer i

	real h
	logical is_fixed
	is_fixed(h) = h .lt. 0. .and. h .gt. -990.

	nparts = 1
	do i=2,nl-1
	  if( is_fixed(ht(i)) ) nparts = nparts + 1
	end do

	end

c********************************************************

	subroutine add_part(nf,xf,yf,hf,nn,xn,yn,hn)

	implicit none

	integer nf,nn
	real xf(1),yf(1),hf(1)
	real xn(1),yn(1),hn(1)

	integer i

	do i=1,nn
	  nf = nf + 1
	  xf(nf) = xn(i)
	  yf(nf) = yn(i)
	  hf(nf) = hn(i)
	end do

	end

c********************************************************

	subroutine extract_part(ipart,nl,xt,yt,ht,nn,xn,yn,hn)

	implicit none

	integer ipart
	integer nl,nn
	real xt(1),yt(1),ht(1)
	real xn(1),yn(1),hn(1)

	real h
	logical is_fixed
	is_fixed(h) = h .lt. 0. .and. h .gt. -990.

	integer i,is,ip

	if( ipart .eq. 1 ) then
	  is = 1
	else
	  ip = 1
	  do i=2,nl-1
	    if( is_fixed(ht(i)) ) ip = ip + 1
	    if( ip .eq. ipart ) goto 1
	  end do
	  write(6,*) 'ipart = ',ipart
	  stop 'error stop extract_part: no such part found'
    1	  continue
	  is = i
	end if

	call insert_in_part(is,xt,yt,ht,1,xn,yn,hn)
	do i=is+1,nl
	  nn = i - is + 1
	  call insert_in_part(i,xt,yt,ht,nn,xn,yn,hn)
	  if( is_fixed(ht(i)) ) return
	end do

	end

c********************************************************

	subroutine insert_in_part(ifrom,x1,y1,h1,ito,x2,y2,h2)

	implicit none

	integer ifrom,ito
	real x1(1),y1(1),h1(1)
	real x2(1),y2(1),h2(1)

	x2(ito) = x1(ifrom)
	y2(ito) = y1(ifrom)
	h2(ito) = h1(ifrom)

	end

c********************************************************

	subroutine reduce_part(nn,xn,yn,hn,reduct)

	implicit none

	integer nn
	real xn(nn),yn(nn),hn(nn)
	real reduct

	integer ndim
	parameter (ndim=10000)

	real xf(ndim), yf(ndim), hf(ndim)
	real dxy(ndim),dcum(ndim)
	real x2n(ndim), y2n(ndim), h2n(ndim)
	real aux(ndim)
	real dx,dy,dist,hh
	real dtot,rdist,rtot
	integer i,nf,np

	logical bspline
	parameter (bspline=.true.)
	!parameter (bspline=.false.)

	real h
	logical is_fixed
	is_fixed(h) = h .lt. 0. .and. h .gt. -990.

	if( nn .ge. 3 ) then
	  hn(1) = hn(2)
	  hn(nn) = hn(nn-1)
	else
	  if( is_fixed(hn(1)) ) hn(1) = 1.
	  if( is_fixed(hn(nn)) ) hn(nn) = 1.
	end if

	dxy(1) = 0.
	dcum(1) = 0.
	dtot = 0.
	do i=2,nn
	  hh = 0.5*(hn(i)+hn(i-1))
	  dx = xn(i) - xn(i-1)
	  dy = yn(i) - yn(i-1)
	  dist = sqrt( dx*dx + dy*dy )
	  dxy(i) = dist / hh
	  dtot = dtot + dxy(i)
	  dcum(i) = dtot
	end do

	do i=1,nn
	  write(6,*) 'dxy ',i,dcum(i),xn(i),yn(i),hn(i)
	end do

	if( bspline ) then
	  np = nint(dtot/reduct)
	  rdist = dtot/np		!here we adjust reduct
	  call spline_init(nn,dcum,xn,aux,x2n)
	  call spline_init(nn,dcum,yn,aux,y2n)
	  call spline_init(nn,dcum,hn,aux,h2n)
	else
	  rdist = reduct
	end if

	write(6,*) 'before insert 1'
	call insert_in_part(1,xn,yn,hn,1,xf,yf,hf)
	write(6,*) 'after insert 1'
	
	if( bspline ) then
	  rtot = 0.
	  do i=2,np
	    rtot = rtot + rdist
	    call spline_eval(nn,dcum,xn,x2n,rtot,xf(i))
	    call spline_eval(nn,dcum,yn,y2n,rtot,yf(i))
	    call spline_eval(nn,dcum,hn,h2n,rtot,hf(i))
	  end do
	  nf = np
	else
	  nf = 1
	  rtot = 0.
	  do i=2,nn-1
	    dist = dxy(i)
	    rtot = rtot + dist
	    if( rtot .gt. rdist ) then
	      nf = nf + 1
	      call insert_in_part(i,xn,yn,hn,nf,xf,yf,hf)
	      rtot = 0.
	      !rtot = rtot - rdist	!other possibility
	    end if
	  end do
	end if

	nf = nf + 1
	call insert_in_part(nn,xn,yn,hn,nf,xf,yf,hf)
	write(6,*) 'after insert last'

	nn = nf
	do i=1,nf
	  xn(i) = xf(i)
	  yn(i) = yf(i)
	  hn(i) = hf(i)
	end do

	end

c********************************************************

	subroutine reduce(reduct,xt,yt,ht,nl)

c reduces points in line

	implicit none

	real reduct
	real xt(1)
	real yt(1)
	real ht(1)
	integer nl

	integer i,nnew
	real rr,rtot
	real dist,dx,dy
	real h

c very simplicistic approach

	rr = reduct
	rtot = 1.

	if( nl .lt. 3 * rr ) then
	  rr = nl / 3.
	end if

	nnew = 0
	rtot = 0.
	if( ht(1) .lt. 0. .or. ht(1) .ge. 0. ) then	!bug fix 23.09.2004
	    nnew = nnew + 1
	    xt(nnew) = xt(1)
	    yt(nnew) = yt(1)
	    ht(nnew) = ht(1)
	end if
	do i=2,nl
	  dx = xt(i) - xt(i-1)
	  dy = yt(i) - yt(i-1)
	  h = ht(i)
	  dist = sqrt( dx*dx + dy*dy )
	  !write(6,*) dist,h,dist/h,rtot,reduct
	  if( h .gt. 0. ) dist = dist / h
	  rtot = rtot + dist
	  if( rtot .gt. reduct .or. h .lt. 0. ) then
	    nnew = nnew + 1
	    xt(nnew) = xt(i)
	    yt(nnew) = yt(i)
	    ht(nnew) = ht(i)
	    rtot = 0.
	  end if
	end do

c	if( nnew .lt. 3 ) stop 'error stop reduce: line too short...'
	if( nnew .lt. 3 ) then
	  write(6,*) 'line too short -> eliminated'
	  nnew = 0
	end if

	write(6,*) 'reduce: new nodes = ',nnew
	nl = nnew

	end

c********************************************************

	subroutine intpdep(ht,nl,bperiod)

c interpolates depth values for points in line
c
c negative values are not interpolated and are left alone

	implicit none

	real ht(1)
	integer nl
	logical bperiod

	integer i
	integer ifirst,ilast,inext,nonzero
	integer nval
	real value,vfirst,vlast
	real aux
	logical bdebug

	bdebug = .false.

c------------------------------------------------------------
c look for values greater than 0
c------------------------------------------------------------

	nonzero = 0
	ifirst = 0
	ilast = 0

	do i=1,nl
	  if( ht(i) .gt. 0. ) then
	    nonzero = nonzero + 1
	    if( nonzero .eq. 1 ) ifirst = i
	    ilast = i
	  end if
	end do

c------------------------------------------------------------
c if no value or 1 value found -> set everything constant
c------------------------------------------------------------

	if( nonzero .le. 1 ) then
	  if( nonzero .eq. 0 ) then
	    value = 1.
	  else
	    value = ht(ifirst)
	  end if
	  do i=1,nl
	    if( ht(i) .ge. 0. .or. ht(i) .lt. -990. ) ht(i) = value
	  end do
	  return
	end if

c------------------------------------------------------------
c interpolate extreme values
c------------------------------------------------------------

	vfirst = ht(ifirst)
	vlast = ht(ilast)
	nval = nl - ilast + ifirst
	aux = (vfirst-vlast)/nval
	inext = 0

	if( bdebug ) write(6,*) bperiod
	if( bdebug ) write(6,*) ifirst,ilast,nval,vfirst,vlast

	i = 1
	do while( ht(i) .lt. 0 )
	  i = i + 1
	end do
	if( bperiod ) then
	  inext = nl - ilast + i
	  value = vlast + inext * aux
	else
	  value = ht(ifirst)
	end if
	if( bdebug ) write(6,*) i,inext,value,ht(i)
	ht(i) = value
	ifirst = i

	i = nl
	do while( ht(i) .lt. 0 )
	  i = i - 1
	end do
	if( bperiod ) then
	  inext = i - ilast
	  value = vlast + inext * aux
	else
	  value = ht(ilast)
	end if
	if( bdebug ) write(6,*) i,inext,value,ht(i)
	ht(i) = value
	ilast = i

c------------------------------------------------------------
c interpolate central values
c------------------------------------------------------------

	do while( ifirst .lt. ilast )
	  do i=ifirst+1,nl
	    if( ht(i) .gt. 0. ) goto 1
	  end do
	  stop 'error stop intpdep: impossible branch (1)'
    1	  continue
	  inext = i
	  nval = inext - ifirst
	  if( bdebug ) write(6,*) ifirst,inext,ilast,ht(ifirst),ht(inext)
	  do i=ifirst+1,inext-1
	    value = ht(ifirst) + (i-ifirst) * (ht(inext)-ht(ifirst))/nval
	    if( ht(i) .ge. 0. ) ht(i) = value
	    if( bdebug ) write(6,*) i,value
	  end do
	  ifirst = inext
	end do

c------------------------------------------------------------
c do some sanity check
c------------------------------------------------------------

	do i=1,nl
	  if( ht(i) .eq. 0 ) then
	    write(6,*) i,nl,ht(i)
	    stop 'error stop intpdep: values not interpolated'
	  end if
	  if( bdebug ) write(6,*) i,ht(i)
	end do

c------------------------------------------------------------
c end of routine
c------------------------------------------------------------

	end

c********************************************************


