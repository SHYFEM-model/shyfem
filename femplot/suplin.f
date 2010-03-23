c
c $Id: suplin.f,v 1.5 2010-03-11 15:32:16 georg Exp $
c
c routines for section plot (scalars)
c
c revision log :
c
c 14.09.2009    ggu     routines written from scratch
c 09.10.2009    ggu     vertical plot nearly ready
c 14.10.2009    ggu     vertical plot finished (scalar and velocities)
c
c************************************************************************

	subroutine plot_sect(sv)

c plots section

	implicit none

        include 'param.h'

	real sv(nlvdim,nkndim)		!scalar to be plotted

	integer nldim
	parameter (nldim=200)

	real hev(1)
	common /hev/hev
	real xgv(1), ygv(1)
	common /xgv/xgv, /ygv/ygv
	real hlv(1)
	common /hlv/hlv
	integer ilhkv(1)
	common /ilhkv/ilhkv
        integer nlvdi,nlv
        common /level/ nlvdi,nlv

	real val(0:2*nlvdim,nldim)

	character*80 file
	character*80 string,line
	logical bhoriz
	integer it,i,k1,k2,ie,l,lbot,ltop,j
	integer ltot
	integer nr,nc,mode,ir
	real xmin,ymin,xmax,ymax
	real xrmin,yrmin,xrmax,yrmax
	real d,dd,dx,dy,htop,htot,hbot
	real x0,y0,x1,x2,y1,y2
	real vmin,vmax
	real rrl,rrd,x,y,xtick,ytick,ylast,rdist,rmax
	real xcm,ycm
	real fact
	integer ndec,nctick

	character*80 vtitle,xtitle,ytitle,rtitle,ltitle
	logical bgrid
	integer ivert
	integer n
	integer nodes(nldim)
	integer elems(nldim)
	integer lelems(nldim)
	real helems(nldim)
	integer ll,lvmax
	real rl,rd,hvmax
	real x0s,y0s,x1s,y1s
	save n,nodes,elems,helems,lelems,rl,rd,ll
	save bgrid,ivert,hvmax,lvmax
	save vtitle,xtitle,ytitle,ltitle,rtitle

	logical inboxdim_noabs
	integer gettime
	integer ialfa
	real getpar
	real hlog,divdist,roundm

	integer icall
	save icall
	data icall / 0 /

c----------------------------------------------------------------
c initialize
c----------------------------------------------------------------

	if( icall .eq. 0 ) then
	  call getfnm('vsect',file)
	  call line_read_nodes(file,nldim,n,nodes)
	  call line_find_elements(n,nodes,hev,hlv,elems,helems,lelems)
	  call line_find_min_max(n,nodes,helems,lelems,xgv,ygv,rl,rd,ll)

	  bgrid = nint(getpar('ivgrid')) .ne. 0
	  ivert = nint(getpar('ivert'))
	  lvmax = nint(getpar('lvmax'))
	  hvmax = getpar('hvmax')

	  call set_max_dep_lay(nlv,hlv,rd,ll,hvmax,lvmax) !vertical range

	  call getfnm('vtitle',vtitle)
	  call getfnm('xtitle',xtitle)
	  call getfnm('ytitle',ytitle)
	  call getfnm('ltitle',ltitle)
	  call getfnm('rtitle',rtitle)
	end if

	icall = icall + 1

	it = gettime()

        call qstart
        call annotes(it,'vertical plot')
	call annote

c----------------------------------------------------------------
c prepare data
c----------------------------------------------------------------

	write(6,*) 'plotting section: ',it,n,rl,rd

	call line_insert_scalars(n,nodes,ilhkv,nlvdim,sv,val,vmin,vmax)
	call colauto(vmin,vmax)
	write(6,*) 'min/max on line: ',vmin,vmax

c----------------------------------------------------------------
c set viewport
c----------------------------------------------------------------

	xrmin = 0.
	yrmin = 0.
	xrmax = 1.
	yrmax = 1.

	call qgetvp(xmin,ymin,xmax,ymax)
	write(6,*) 'plot_sect',xmin,ymin,xmax,ymax

	ymax = ymax / 2.			!empirical
	call qsetvp(xmin,ymin,xmax,ymax)
	call qworld(xrmin,yrmin,xrmax,yrmax)
	call pbox(xrmin,yrmin,xrmax,yrmax)	!plot outer box
	call bpixel

        bhoriz = inboxdim_noabs('sect',x0s,y0s,x1s,y1s)

        if( bhoriz ) then
	  xmin = xmin + 1.5
	  ymin = ymin + 1.5
	  xmax = xmax - 1.5
	  ymax = ymax - 1.5
	else
	  xmin = xmin + 1.5
	  ymin = ymin + 1.5
	  xmax = xmax - 2.5
	  ymax = ymax - 1.5
	end if

	call qsetvpnc(xmin,ymin,xmax,ymax)	!this is the plotting area

c-----------------------------------------------------------------
c set world coordinates
c-----------------------------------------------------------------

	xrmin = 0.
	xrmax = rl
	yrmin = -rd
	yrmax = 0.

	if( ivert .eq. 1 ) yrmin = -ll

	call qworld(xrmin,yrmin,xrmax,yrmax)

c--------------------------------------------------------------------
c plot scalar
c--------------------------------------------------------------------

	d = 0.
	do i=2,n
	  k1 = nodes(i-1)
	  k2 = nodes(i)
	  htot = helems(i)
	  ltot = lelems(i)
	  htot = min(htot,hvmax)
	  ltot = min(ltot,lvmax)
	  ie = elems(i)
	  dx = xgv(k2) - xgv(k1)
	  dy = ygv(k2) - ygv(k1)
	  dd = sqrt( dx**2 + dy**2 )
	  x1 = d
	  x2 = d + dd
	  y1 = yrmin
	  y2 = -htot
	  if( ivert .eq. 1 ) y2 = -ltot 
	  if( ivert .eq. 2 ) y2 = -hlog(htot,rd)
	  call qgray(0.5)
	  call qrfill(x1,y1,x2,y2)
	  htop = 0.
	  do l=1,ltot
	    hbot = hlv(l)
	    ltop = 2*l - 2
	    if( hbot .gt. htot ) hbot = htot
	    if( ivert .eq. 1 ) then
	      htop = l-1
	      hbot = l
	    else if( ivert .eq. 2 ) then
	      hbot = hlog(hbot,rd)
	    end if
	    call plot_rect(x1,-htop,x2,-hbot,val(ltop,i-1),val(ltop,i))
	    call qgray(0.5)
	    if( bgrid ) call pbox(x1,-htop,x2,-hbot)
	    htop = hbot
	  end do
	  d = d + dd
	end do

c--------------------------------------------------------------------
c labeling
c--------------------------------------------------------------------

        call qfont('Times-Roman')
        call qgray(0.0)
	call qtxts(9)
	call qcm(xcm,ycm)

	rrl = divdist(rl,5,0)
	nr = rl/rrl
	nc = -1		!do not write decimal point
	mode = -1	!left flushing
	y = yrmin
	ytick = 0.1 * ycm

	do i=0,nr			!x-axis
	  x = i * rrl
	  ir = ialfa(x,string,nc,mode)
	  call qline(x,y,x,y-ytick)
          call qtxtcr(0.,+2.5)
          call qtext(x,y,string(1:ir))
	end do

	ylast = -100.	!big negative number
	x = xrmin
	xtick = 0.1 * xcm
	rmax = -yrmin
	rrd = 0.
	rdist = max(1,nint(ll/7.))	!for ivert == 1

	do while( rrd .le. rmax )
	  y = rrd
	  if( ivert .eq. 2 ) y = hlog(y,rmax)
	  !write(6,*) rmax,rrd,y
          if( rrd .ne. 0.5 .and. y-ylast .gt. rmax/50. ) then
	    ir = ialfa(rrd,string,nc,mode)
	    call qline(x,-y,x-xtick,-y)
            call qtxtcr(+1.,0.)
            call qtext(x-2*xtick,-y,string(1:ir))
	    ylast = y
	  end if
	  if( ivert .eq. 1 ) then
	    rrd = rrd + rdist
	  else
	    rrd = roundm(rrd+0.5,1)
	  end if
	end do

c--------------------------------------------------------------------
c titles
c--------------------------------------------------------------------

	call qcomm('Titles')
        call qtxtcr(0.,0.)

	call qtxts(15)
        call qtext(rl/2.,-yrmin/15.,vtitle)	!title

	call qtxts(12)
        call qtext(rl/2.,yrmin*1.1,xtitle)	!x-axis

        call qtxtcr(+1.,0.)
        call qtxtcr(0.,0.)
        call qtext(0.,yrmin*1.07,ltitle)	!left point
        call qtxtcr(-1.,0.)
        call qtxtcr(0.,0.)
        call qtext(rl,yrmin*1.07,rtitle)	!right point

        call qtxtr(90.)
        call qtxtcr(0.,0.)
	if( ivert .eq. 1 .and. ytitle .eq. 'Depth [m]' ) then
		ytitle = 'Layers'
	end if
        call qtext(-rl/15.,yrmin*0.5,ytitle)	!y-axis

	call pbox(xrmin,yrmin,xrmax,yrmax)

c--------------------------------------------------------------------
c color bar
c--------------------------------------------------------------------

        if( bhoriz ) then
	  dx = xmax - xmin
	  dy = ymax - ymin
	  xmax = xmin + x1s * dx
	  ymax = ymin + y1s * dy
	  xmin = xmin + x0s * dx
	  ymin = ymin + y0s * dy
	else
	  dx = 0.2
	  dy = 1.
	  xmin = xmax + dx
	  ymin = ymin + dy
	  xmax = xmax + 2.5 - dx
	  ymax = ymax - dy
	end if

	xrmin = 0.
	yrmin = 0.
	xrmax = 1.
	yrmax = 1.

	call qsetvp(xmin,ymin,xmax,ymax)	!this is the plotting area
	call qworld(xrmin,yrmin,xrmax,yrmax)

	call qcomm('Color Bar')
	call pbox(x0,y0,x1,y1)

        fact = getpar('faccol')
        ndec = nint(getpar('ndccol'))
        nctick = nint(getpar('nctick'))
        call getfnm('legcol',line)

	call qgray(0.0)
	call colbar(line,bhoriz,nctick,ndec,fact,xrmin,yrmin,xrmax,yrmax)

c--------------------------------------------------------------------
c end of routine
c--------------------------------------------------------------------

        call qend

	end

c************************************************************************
c************************************************************************
c************************************************************************

	function hlog(x,rd)

	implicit none

	real hlog
	real x,rd

	hlog = rd * log(x+1.)/log(rd+1.)

	end

c************************************************************************

	subroutine plot_rect(x1,y1,x2,y2,v1,v2)

	implicit none

	real x1,y1,x2,y2
	real v1(3),v2(3)

	include 'color.h'

	logical bfirst
	real xm,ym,vm
	real x(3),y(3),f(3)

	bfirst = .false.

	xm = (x1+x2)/2.
	ym = (y1+y2)/2.
	vm = (v1(2)+v2(2))/2.

	if( bfirst ) then

	call setxyf(x1,x1,x2,y1,ym,y1,v1(1),v1(2),v2(1),x,y,f)
	call plcol(x,y,f,ciso,fiso,isoanz+1,fnull)
	call setxyf(x1,x2,x2,ym,ym,y1,v1(2),v2(2),v2(1),x,y,f)
	call plcol(x,y,f,ciso,fiso,isoanz+1,fnull)

	call setxyf(x1,x1,x2,ym,y2,ym,v1(2),v1(3),v2(2),x,y,f)
	call plcol(x,y,f,ciso,fiso,isoanz+1,fnull)
	call setxyf(x1,x2,x2,y2,y2,ym,v1(3),v2(3),v2(2),x,y,f)
	call plcol(x,y,f,ciso,fiso,isoanz+1,fnull)

	else

	call setxyf(x1,x1,xm,y1,ym,ym,v1(1),v1(2),vm,x,y,f)
	call plcol(x,y,f,ciso,fiso,isoanz+1,fnull)
	call setxyf(x1,xm,x2,y1,ym,y1,v1(1),vm,v2(1),x,y,f)
	call plcol(x,y,f,ciso,fiso,isoanz+1,fnull)
	call setxyf(xm,x2,x2,ym,ym,y1,vm,v2(2),v2(1),x,y,f)
	call plcol(x,y,f,ciso,fiso,isoanz+1,fnull)

	call setxyf(x1,x1,xm,ym,y2,ym,v1(2),v1(3),vm,x,y,f)
	call plcol(x,y,f,ciso,fiso,isoanz+1,fnull)
	call setxyf(x1,xm,x2,y2,ym,y2,v1(3),vm,v2(3),x,y,f)
	call plcol(x,y,f,ciso,fiso,isoanz+1,fnull)
	call setxyf(xm,x2,x2,ym,y2,ym,vm,v2(3),v2(2),x,y,f)
	call plcol(x,y,f,ciso,fiso,isoanz+1,fnull)

	end if

	end

c************************************************************************

	subroutine setxyf(x1,x2,x3,y1,y2,y3,f1,f2,f3,x,y,f)
	real x(3),y(3),f(3)
	x(1)=x1
	x(2)=x2
	x(3)=x3
	y(1)=y1
	y(2)=y2
	y(3)=y3
	f(1)=f1
	f(2)=f2
	f(3)=f3
	end

c************************************************************************

	subroutine line_insert_scalars(n,nodes,ilhkv,nlvdim,sv
     +					,values,vmin,vmax)

c inserts scalar values into matrix section

	implicit none

	integer n
	integer nodes(n)
	integer ilhkv(1)		!number of layers in node
	integer nlvdim
	real sv(nlvdim,1)
	real values(0:2*nlvdim,n)
	real vmin,vmax

	integer i,k,lmax,l,llayer
	integer ipext

	vmin = 1.e+30
	vmax = -1.e+30

	!write(6,*) 'line_insert_scalars'
	!write(6,*) n,nlvdim

	do i=1,n
	  k = nodes(i)
	  lmax = ilhkv(k)
	  values(0,i) = sv(1,k)
	  values(1,i) = sv(1,k)
	  vmin = min(vmin,sv(1,k))
	  vmax = max(vmax,sv(1,k))
	  do l=2,lmax
	    llayer = 2*l-1
	    values(llayer,i) = sv(l,k)
	    values(llayer-1,i) = 0.5 * ( sv(l,k) + sv(l-1,k) )
	    vmin = min(vmin,sv(l,k))
	    vmax = max(vmax,sv(l,k))
	  end do
	  values(2*lmax,i) = sv(lmax,k)
	  !write(78,*) k,ipext(k),lmax
	  !write(78,*) (sv(l,k),l=1,lmax)
	  !write(79,*) k,ipext(k),lmax,2*lmax
	  !write(79,*) (values(l,i),l=0,2*lmax)
	end do

	end

c************************************************************************

	subroutine line_find_min_max(n,nodes,helems,lelems
     +					,xgv,ygv,rl,rd,ll)

c finds length and max depth of line

	implicit none

	integer n
	integer nodes(n)
	real helems(n)		!depth in chosen elements
	integer lelems(n)
	real xgv(1), ygv(1)
	real rl,rd		!length and depth (return)
	integer ll

	integer i,k1,k2
	real dx,dy

	rl = 0.
	rd = 0.
	ll = 0

	do i=2,n
	  k1 = nodes(i-1)
	  k2 = nodes(i)
	  dx = xgv(k2) - xgv(k1)
	  dy = ygv(k2) - ygv(k1)
	  rl = rl + sqrt( dx**2 + dy**2 )
	  rd = max(rd,helems(i))
	  ll = max(ll,lelems(i))
	end do

	end

c************************************************************************

	subroutine line_find_elements(n,nodes,hev,hlv
     +					,elems,helems,lelems)

c finds elements along line given by nodes

	implicit none

	integer n
	integer nodes(n)
	real hev(1)		!depth in elements
	real hlv(1)		!layer structure
	integer elems(n)	!element number of chosen elements
	real helems(n)		!depth in chosen elements
	integer lelems(n)	!layers in element number of 

	integer i,k1,k2,ie1,ie2,l
	real h
	integer ipext,ieext

	do i=2,n
	  k1 = nodes(i-1)
	  k2 = nodes(i)
	  call find_elems_to_segment(k1,k2,ie1,ie2)
	  elems(i) = ie1
	  if( ie2 .gt. 0 ) then
	    if( hev(ie2) .gt. hev(ie1) ) elems(i) = ie2
	  end if
	  helems(i) = hev(elems(i))
	  !write(6,*) 'line_find_elements: ',k1,k2,ie1,ie2,hev(ie1)
	  !write(6,*) '                    ',ipext(k1),ipext(k2)
	  !write(6,*) '                    ',ieext(ie1),ieext(ie2)
	end do

	do i=2,n
	  h = helems(i)
	  l = 1
	  do while( hlv(l) .lt. h ) 
	    l = l + 1
	  end do
	  lelems(i) = l
	end do
	  
	end

c************************************************************************

	subroutine line_read_nodes(file,ndim,n,nodes)

c reads line given by nodes

	implicit none

	character*(*) file
	integer ndim
	integer n
	integer nodes(ndim)

	integer iunit,nn
	integer ifileo,ipint

	n = 0

	iunit = ifileo(0,file,'form','old')
	if( iunit .le. 0 ) goto 99 

    1	continue
	  read(iunit,*,end=2) nn
	  if( nn .le. 0 ) goto 2
	  n = n + 1
	  if( n .gt. ndim ) stop 'error stop line_read_nodes: ndim'
	  nodes(n) = ipint(nn)
	goto 1
    2	continue

	write(6,*) 'finished reading line section: ',n

	close(iunit)

	return
   99	continue
	write(6,*) 'Cannot open section file: '
	write(6,*) file
	stop 'error stop line_read_nodes: no such file'
	end

c************************************************************************

	subroutine set_max_dep_lay(nlv,hlv,rd,ll,hvmax,lvmax)

c sets hvmax and lvmax

	implicit none

	integer nlv
	real hlv(1)
	real rd
	integer ll
	real hvmax
	integer lvmax

	integer l

        if( hvmax .gt. 0. .and. lvmax .gt. 0 ) then
          write(6,*) 'hvmax,lvmax: ',hvmax,lvmax
          write(6,*) 'only one of hvmax - lvmax can be given'
          stop 'error stop set_max_dep_lay: hvmax/lvmax'
        end if

        if( hvmax .gt. 0. ) then
          do l=1,nlv
            if( hlv(l) .gt. hvmax ) goto 7
          end do
    7     continue
          lvmax = l             !could be one higher than max layer
          ll = lvmax
          rd = hvmax
        else if( lvmax .gt. 0 ) then
          if( lvmax .gt. nlv .and. nlv .gt. 1 ) then
            hvmax = hlv(nlv) + hlv(nlv) - hlv(nlv-1)
          else
            hvmax = hlv(lvmax)
          end if
          ll = lvmax
          rd = hvmax
        else
          hvmax = rd
          lvmax = ll
        end if

	end

c************************************************************************

