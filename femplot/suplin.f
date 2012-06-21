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
c 26.03.2010    ggu     vertical plot for velocity finished 
c 13.04.2010    ggu     adapted also to spherical coordinates
c 15.04.2010    ggu     fix bug where lower layer is plotted with value 0
c 29.09.2010    ggu     finished velocity plot with reference arrow
c 09.10.2010    ggu     better labeling of reference arrow
c 21.12.2010    ggu     plotting vertical vector of sigma (not finished)
c 18.08.2011    ggu     better error check of node list
c 24.08.2011    ggu     small changes to avoid run time error
c 24.08.2011    ggu     plot real depth for zeta layers
c 14.11.2011    ggu     hybrid levels introduced
c 23.11.2011    ggu     in line_find_elements() adjust depth for hybrid
c 27.01.2012    deb&ggu adjusted for hybrid coordinates
c 20.06.2012    ggu     plots bottom also for sigma layers (plot_bottom())
c
c************************************************************************

	subroutine plot_sect(bvel,sv)

c plots section

	implicit none

        include 'param.h'

	logical bvel			!plot velocities
	real sv(nlvdim,nkndim)		!scalar to be plotted

	integer nldim
	parameter (nldim=200)

	real hev(1)
	common /hev/hev
	real hm3v(3,1)
	common /hm3v/hm3v
	real xgv(1), ygv(1)
	common /xgv/xgv, /ygv/ygv
	real hlv(1)
	common /hlv/hlv
	integer nen3v(3,1)
	common /nen3v/nen3v
	integer ilhkv(1)
	common /ilhkv/ilhkv
        integer nlvdi,nlv
        common /level/ nlvdi,nlv

c elems(1) is not used, etc..

	real val(0:2*nlvdim,nldim)	!scalar value along line
	real vel(3,0:2*nlvdim,nldim)	!projected velocities along line
	integer nodes(nldim)		!nodes along line
	integer elems(nldim)		!elements along line
	integer lelems(nldim)		!levels in element
	integer lnodes(nldim)		!levels in nodes
	real helems(2,nldim)		!depth in elements (for both nodes)
	real dxy(2,nldim)		!direction of projection line
	real xy(nldim)			!linear distance

	real ya(2,0:nlvdim)

	character*80 file
	character*80 string,line
	logical bhoriz,barrow,btwo
	logical blayer,blog
	integer it,i,k1,k2,ie,l,lbot,ltop,j
	integer ltot
	integer nr,nc,mode,ir
	real xmin,ymin,xmax,ymax
	real xrmin,yrmin,xrmax,yrmax
	real xrrmin,yrrmin,xrrmax,yrrmax
	real xfmin,yfmin,xfmax,yfmax
	real d,dd,dx,dy
	real x0,y0,x1,x2,y1,y2
	real u,v
	real vmin,vmax
	real vhmin,vhmax
	real wscale,vv
	real rrlmax,rrdmax,x,y,xtick,ytick,ylast,rdist,rmax,xs
	real h1,h2,yb1,yb2,yt1,yt2,yt,yb
	real ytaux,ymid
	real xcm,ycm
	real fact,r
	integer ndec,nctick
	integer isphe

	real x0s,y0s,x1s,y1s
	real xmid,hmid,umid,wmid
	real scale

	integer n
	save n,nodes,elems,helems,lelems,lnodes,dxy,xy

	integer llmax,lvmax
	save llmax,lvmax
	real rlmax,rdmax,hvmax
	save rlmax,rdmax,hvmax
	logical bgrid
	integer ivert
	save bgrid,ivert

	character*80 vtitle,xtitle,ytitle,ltitle,rtitle
	save vtitle,xtitle,ytitle,ltitle,rtitle
	real ascale,rscale,stip,rwscal
	save ascale,rscale,stip,rwscal
	integer vmode
	save vmode
	real faccol
	save faccol
	real rxscal,ryscal
	save rxscal,ryscal

	logical inboxdim_noabs
	logical bdebug,bdebug_scalar
	integer gettime
	integer ialfa,ichanm
	real getpar
	real hlog,divdist,roundm

	integer icall
	save icall
	data icall / 0 /

c----------------------------------------------------------------
c initialize
c----------------------------------------------------------------

	bdebug = .true.
	bdebug = .false.

	if( icall .eq. 0 ) then
	  call getfnm('vsect',file)
	  isphe = nint(getpar('isphe'))
	  call line_read_nodes(file,nldim,n,nodes)
	  call line_find_elements(n,nodes,nlv,nen3v,hev,hm3v,hlv
     +			,elems,helems,lelems,lnodes)
	  call line_find_min_max(n,nodes,helems,lelems,xgv,ygv
     +			,isphe,rlmax,rdmax,llmax,xy)
	  call make_proj_dir(n,isphe,nodes,xgv,ygv,dxy)

	  bgrid = nint(getpar('ivgrid')) .ne. 0	!plot grid?
	  ivert = nint(getpar('ivert'))	!type of vertical coordinate
	  lvmax = nint(getpar('lvmax'))	!max layer to be plotted
	  hvmax = getpar('hvmax')	!max depth to be plotted

	  call set_max_dep_lay(nlv,hlv,rdmax,llmax,hvmax,lvmax) !vertical range

	  call getfnm('vtitle',vtitle)
	  call getfnm('xtitle',xtitle)
	  call getfnm('ytitle',ytitle)
	  call getfnm('ltitle',ltitle)
	  call getfnm('rtitle',rtitle)

	  ascale = getpar('avscal')	!absolute scale
	  rscale = getpar('rvscal')	!relative scale
	  rwscal = getpar('rwscal')	!extra scale for vertical
	  stip = getpar('svtip')	!arrow tip size
	  vmode = getpar('vmode')	!0=normal vel  1=tang vel as overlay

	  rxscal = getpar('rxscal')	!x scale for reference vector
	  ryscal = getpar('ryscal')	!y scale for reference vector

	  faccol = getpar('faccol')	!factor for velocity (for legend)
	end if

	icall = icall + 1

c----------------------------------------------------------------
c prepare logical variables
c----------------------------------------------------------------

c ivert:   0=depth  1=layers  2=log depth
c	   -1=layers scaled by sigma values

	btwo = .false.
	btwo = .true.				!plot two arrows for reference

	barrow = bvel .and. stip .ge. 0.	!plot arrow

	blayer = abs(ivert) .eq. 1
	blog = ivert .eq. 2

c----------------------------------------------------------------
c start plot
c----------------------------------------------------------------

	it = gettime()

        call qstart
        call annotes(it,'vertical plot')
	call annote

c----------------------------------------------------------------
c prepare data
c----------------------------------------------------------------

	write(6,*) 'plotting section: ',it,n,rlmax,rdmax

	if( bvel ) then
	  call proj_velox(vmode,n,nodes,lnodes,ilhkv,dxy,vel,val
     +					,vmin,vmax,vhmin,vhmax)
	else
	  call line_insert_scalars(n,nodes,lnodes,ilhkv,sv
     +					,val,vmin,vmax)
	  vhmin = vmin	!...
	  vhmax = vmax	!not needed, only to avoid run time error
	end if
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
	if( bdebug ) write(6,*) 'plot_sect: ',xmin,ymin,xmax,ymax

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

c	hvmax == rdmax

	xrmin = 0.
	xrmax = rlmax
	yrmin = -rdmax
	yrmax = 0.

	if( blayer ) yrmin = -llmax

	call qworld(xrmin,yrmin,xrmax,yrmax)

c--------------------------------------------------------------------
c plot scalar
c--------------------------------------------------------------------

	do i=2,n
	  ltot = lelems(i)
	  ltot = min(ltot,lvmax)
	  x1 = xy(i-1)
	  x2 = xy(i)

	  call make_segment_depth(ivert,ltot,helems(1,i),hvmax,hlv,ya)
c	write(77,*) i,ya(1,ltot),ya(2,ltot)

c	  call qgray(0.5)
c	  yb = yrmin
c	  yt = max(ya(1,ltot),ya(2,ltot))
c	  call qrfill(x1,yb,x2,yt)	!land (bottom)
c	write(6,*) 'line_ggu i',i,ltot,lvmax
c	write(6,*) 'line_ggu ya',ya
c	write(6,*) 'line_ggu hlv',nlv,(hlv(l),l=1,nlv)
c	write(6,*) 'line_ggu x',x1,yb,x2,yt

	  call plot_bottom(x1,x2,yrmin,ya(1,ltot))

	  do l=1,ltot
	    ltop = 2*l - 2
	    yt1 = ya(1,l-1)
	    yt2 = ya(2,l-1)
	    yb1 = ya(1,l)
	    yb2 = ya(2,l)
	    call plot_scal(x1,yt1,yb1,x2,yt2,yb2
     +				,val(ltop,i-1),val(ltop,i))
	  end do
	end do

	bdebug_scalar = .true.
	bdebug_scalar = .false.
	if( bdebug_scalar ) then
	  write(6,*) 'scalar: '
	  do i=1,n
	    ltot = lelems(i)
	    write(6,*) 'node ',i,ltot
	    do l=0,2*ltot
	      write(6,*) l,l/2.,val(l,i)
	    end do
	  end do
	end if

c--------------------------------------------------------------------
c plot vector
c--------------------------------------------------------------------

	vhmax = max(abs(vhmin),abs(vhmax))
	scale = rlmax/(2.*vhmax*(n-1))
	if( ascale .gt. 0. ) scale = ascale
	if( ascale .lt. 0. ) then		!scale in cm
	  call qcm(xcm,ycm)
	  scale = -ascale * xcm			!not yet documented
	end if
	scale = scale * rscale 			!adjust scale

	if( barrow ) then
	  write(6,*) 'arrow scale: ',vhmax,ascale,rscale,scale
	end if

	wscale = rwscal
	if( blayer ) wscale = wscale / llmax

	do i=2,n
	  ltot = lelems(i)
	  x1 = xy(i-1)
	  x2 = xy(i)
	  call make_segment_depth(ivert,ltot,helems(1,i),hvmax,hlv,ya)
	  do l=1,ltot
	    ltop = 2*l - 2
	    yt1 = ya(1,l-1)
	    yt2 = ya(2,l-1)
	    yb1 = ya(1,l)
	    yb2 = ya(2,l)

	    if( bdebug .and. l .eq. ltot ) then
	      write(6,*) 'last layer: ',yt1,yb1,yt2,yb2
	    end if

	    if( barrow ) then
	      ymid = 0.25*(yt1+yt2+yb1+yb2)
	      xmid = 0.5*(x1+x2)
	      umid = 0.5*(vel(2,ltop+1,i-1)+vel(2,ltop+1,i))
	      wmid = 0.5*(vel(3,ltop+1,i-1)+vel(3,ltop+1,i))
	      call qgray(0.0)
	      call plot_arrow(xmid,ymid,umid,wscale*wmid,scale,stip)
	    end if
	    if( bgrid ) then
	      call qgray(0.5)
	      call trapez(x1,yt1,yb1,x2,yt2,yb2)
	    end if
	  end do
	end do

c--------------------------------------------------------------------
c plot reference vector
c--------------------------------------------------------------------

	if( barrow ) then
	
	xrrmin = xrmax 
	xrrmax = xrmax + 2.5*(xrmax-xrmin)/(xmax-xmin)
	yrrmax = yrmin
	yrrmin = yrmin + 2.18*(yrmax-yrmin)/(ymax-ymin)

	xrrmin = xrrmin + 0.1*(xrrmax-xrrmin)
	xrrmax = xrrmax - 0.1*(xrrmax-xrrmin)
	yrrmin = yrrmin + 0.1*(yrrmax-yrrmin)
	!yrrmax = yrrmax - 0.1*(yrrmax-yrrmin)

        call qgray(0.0)
	!call pbox(xrrmin,yrrmin,xrrmax,yrrmax)		!debug

	dx = rxscal*(xrrmax-xrrmin)
	dy = ryscal*(yrrmax-yrrmin)
	u = dx/scale
	v = dy/scale
	v = v * wscale
	u = roundm(u,-1)
	v = roundm(v,-1)
	dx = u*scale
	dy = v*scale

	x = xrrmin + (xrrmax-xrrmin-dx)/2.
	y = yrrmax - (yrrmax-yrrmin-dy)/2.
	xfmin = xrrmin + 0.1*(xrrmax-xrrmin)
	yfmin = yrrmax - 0.1*(yrrmax-yrrmin)
	xfmax = xrrmin + 0.9*(xrrmax-xrrmin)
	yfmax = yrrmax - 0.9*(yrrmax-yrrmin)

	if( btwo ) then
	 call plot_arrow(xfmin,yfmin,u,0.,scale,stip)
	 call plot_arrow(xfmin,yfmin,0.,-v,scale,stip)
	else
	  call plot_arrow(xfmin,yfmin,u,-v,scale,stip)
	end if

	mode = -1	!left flushing
        call qfont('Times-Roman')
	call qtxts(10)

	u = u * faccol
	call find_nc(u,nc)
	ir = ialfa(u,string,nc,mode)
	write(6,*) 'label x: ',u,nc,ir,string
        call qtxtcr(1.,-2.5)
        !call qtext(x+dx,y,string(1:ir))
        !call qtext(x+dx/2.,y,string(1:ir))	!center in x
        call qtext(xfmax,yfmin,string(1:ir))	!center in x

	v = v * faccol
	call find_nc(-v,nc)
	ir = ialfa(-v,string,nc,mode)
	write(6,*) 'label x: ',-v,nc,ir,string
        call qtxtcr(-1.,-1.5)
        call qtext(xfmin,yfmin-dy,string(1:ir))

	x = xrrmin + (xrrmax-xrrmin)/2.
	y = yrrmax - (yrrmax-yrrmin)/2.
	call get_vel_unit(faccol,string)
	ir = ichanm(string)
        call qtxtcr(1.,-1.)
        call qtext(xfmax,yfmin-dy,string(1:ir))

	end if

c--------------------------------------------------------------------
c labeling
c--------------------------------------------------------------------

        call qfont('Times-Roman')
        call qgray(0.0)
	call qtxts(9)
	call qcm(xcm,ycm)

	rrlmax = divdist(rlmax,5,0)
	nr = rlmax/rrlmax
	nc = -1		!do not write decimal point
	mode = -1	!left flushing
	y = yrmin
	ytick = 0.1 * ycm

	do i=0,nr			!x-axis
	  xs = i * rrlmax		!number to be written
	  x = xs			!x-position
	  ir = ialfa(xs,string,nc,mode)
	  call qline(x,y,x,y-ytick)
          call qtxtcr(0.,+2.5)
          call qtext(x,y,string(1:ir))
	end do

	ylast = -100.	!big negative number
	x = xrmin
	xtick = 0.1 * xcm
	rmax = -yrmin
	rrdmax = 0.
	rdist = max(1,nint(llmax/7.))	!for blayer

	do while( rrdmax .le. rmax )
	  y = rrdmax
	  if( blog ) y = hlog(y,rmax)
	  !write(6,*) rmax,rrdmax,y
          if( rrdmax .ne. 0.5 .and. y-ylast .gt. rmax/50. ) then
	    ir = ialfa(rrdmax,string,nc,mode)
	    call qline(x,-y,x-xtick,-y)
            call qtxtcr(+1.,0.)
            call qtext(x-2*xtick,-y,string(1:ir))
	    ylast = y
	  end if
	  if( blayer ) then
	    rrdmax = rrdmax + rdist
	  else
	    rrdmax = roundm(rrdmax+0.5,1)
	  end if
	end do

c--------------------------------------------------------------------
c titles
c--------------------------------------------------------------------

	call qcomm('Titles')
        call qtxtcr(0.,0.)

	call qtxts(15)
        call qtext(rlmax/2.,-yrmin/15.,vtitle)	!title

	call qtxts(12)
        call qtext(rlmax/2.,yrmin*1.1,xtitle)	!x-axis

        call qtxtcr(+1.,0.)
        call qtxtcr(0.,0.)
        call qtext(0.,yrmin*1.07,ltitle)	!left point
        call qtxtcr(-1.,0.)
        call qtxtcr(0.,0.)
        call qtext(rlmax,yrmin*1.07,rtitle)	!right point

        call qtxtr(90.)
        call qtxtcr(0.,0.)
	if( blayer .and. ytitle .eq. 'Depth [m]' ) then
		ytitle = 'Layers'
	end if
        call qtext(-rlmax/15.,yrmin*0.5,ytitle)	!y-axis

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
	  if( barrow ) then	!make space for reference arrow)
	    ymin = ymin + 0.99
	    ymax = ymax + 0.99
	  end if
	end if

	xrmin = 0.
	yrmin = 0.
	xrmax = 1.
	yrmax = 1.

	call qsetvp(xmin,ymin,xmax,ymax)	!this is the plotting area
	call qworld(xrmin,yrmin,xrmax,yrmax)

	call qcomm('Color Bar')
	!call pbox(xrmin,yrmin,xrmax,yrmax)	!for debug

        fact = getpar('faccol')
        ndec = nint(getpar('ndccol'))
        nctick = nint(getpar('nctick'))
        call getfnm('legcol',line)

	call qgray(0.0)
	call colbar(line,bhoriz,nctick,ndec,fact,xrmin,yrmin,xrmax,yrmax)

c--------------------------------------------------------------------
c plot in/out indicator
c--------------------------------------------------------------------

	if( bvel .and. vmode .eq. 0 ) then
	  call qcomm('in/out indicator')
	  if( bhoriz ) then
	    dy = (ymax-ymin)/(xmax-xmin)
	    call qworld(xrmin,-dy/2.,xrmax,dy/2.)
	    x1 = 0.05
	    y1 = 0.
	    x2 = 0.95
	    y2 = 0.
	  else
	    dx = (xmax-xmin)/(ymax-ymin)
	    call qworld(-dx/2.,yrmin,dx/2.,yrmax)
	    x1 = 0.
	    y1 = 0.05
	    x2 = 0.
	    y2 = 0.95
	  end if
	  r = 0.03
	  call qgray(0.0)
	  call colminmax(vmin,vmax)
	  if( vmin .lt. 0. ) then
	    call pcircle(x1,y1,r)
	    call pcirclefill(x1,y1,r/10.)
	  else if( vmin .gt. 0. ) then
	    call pcircle(x1,y1,r)
	    call pcross(x1,y1,r)
	  end if
	  if( vmax .lt. 0. ) then
	    call pcircle(x2,y2,r)
	    call pcirclefill(x2,y2,r/10.)
	  else if( vmax .gt. 0. ) then
	    call pcircle(x2,y2,r)
	    call pcross(x2,y2,r)
	  end if
	end if

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

	subroutine plot_bottom(x1,x2,yrmin,ya_bot)

c plots bottom with grey color

	implicit none

	real x1,x2,yrmin,ya_bot(2)

	integer n
	parameter (n=4)
	real x(n),y(n)

	call qgray(0.5)

	x(1) = x1
	y(1) = yrmin
	x(2) = x2
	y(2) = yrmin
	x(3) = x2
	y(3) = ya_bot(2)
	x(4) = x1
	y(4) = ya_bot(1)
	
	call qafill(n,x,y)	!land (bottom)

	end

c************************************************************************

	subroutine plot_scal(x1,yt1,yb1,x2,yt2,yb2,v1,v2)

c x coords must be the same, but y coords may be different

	implicit none

	real x1,yt1,yb1,x2,yt2,yb2
	real v1(3),v2(3)

	include 'color.h'

	!integer icsave
	real xm,ym1,ym2,ym,vm
	real x(3),y(3),f(3)

	xm = (x1+x2)/2.
	ym1 = (yb1+yt1)/2.
	ym2 = (yb2+yt2)/2.
	ym = (ym1+ym2)/2.
	vm = (v1(2)+v2(2))/2.

	!write(6,*) '--------- plot_scal ---------'
	!write(6,*) x1,yt1,yb1
	!write(6,*) x2,yt2,yb2
	!write(6,*) xm,ym1,ym2,ym
	
	call set_auto_color_table

	!first plot upper triangles, than lower ones

	call setxyf(x1,x1,xm,yt1,ym1,ym,v1(1),v1(2),vm,x,y,f)
	call plcol(x,y,f,ciso,fiso,isoanz+1,fnull)
	call setxyf(x1,xm,x2,yt1,ym,yt2,v1(1),vm,v2(1),x,y,f)
	call plcol(x,y,f,ciso,fiso,isoanz+1,fnull)
	call setxyf(xm,x2,x2,ym,ym2,yt2,vm,v2(2),v2(1),x,y,f)
	call plcol(x,y,f,ciso,fiso,isoanz+1,fnull)

	call setxyf(x1,x1,xm,ym1,yb1,ym,v1(2),v1(3),vm,x,y,f)
	call plcol(x,y,f,ciso,fiso,isoanz+1,fnull)
	call setxyf(x1,x2,xm,yb1,yb2,ym,v1(3),v2(3),vm,x,y,f)
	call plcol(x,y,f,ciso,fiso,isoanz+1,fnull)
	call setxyf(xm,x2,x2,ym,yb2,ym2,vm,v2(3),v2(2),x,y,f)
	call plcol(x,y,f,ciso,fiso,isoanz+1,fnull)

	call reset_auto_color_table

	end

c************************************************************************

	subroutine plot_rect(x1,y1,x2,y2,v1,v2)

	implicit none

	real x1,y1,x2,y2
	real v1(3),v2(3)

	include 'color.h'

	logical bfirst
	!integer icsave
	real xm,ym,vm
	real x(3),y(3),f(3)

	bfirst = .false.	!what is this? -> old way: no middle point

	xm = (x1+x2)/2.
	ym = (y1+y2)/2.
	vm = (v1(2)+v2(2))/2.

	call set_auto_color_table

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

	!first plot upper triangles, than lower ones

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

	call reset_auto_color_table

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

	subroutine plot_arrow(x,y,u,w,scale,stip)

	implicit none

	real x,y,u,w,scale,stip

	real uu,ww,s

	uu = x + u*scale
	ww = y + w*scale
	s = stip

	call fcmpfeil(x,y,uu,ww,s)

	end

c************************************************************************

	subroutine make_proj_dir(n,isphe,nodes,xgv,ygv,dxy)

c computes projection line for nodes

	implicit none

	integer n
	integer isphe
	integer nodes(n)
	real xgv(1), ygv(1)
	real dxy(2,n)			!direction of line (for projection)

	integer i,i1,i2,k1,k2
	real dx,dy,dd
	real y,xfact,yfact

	do i=1,n
	  i1 = max(1,i-1)
	  i2 = min(n,i+1)
	  k1 = nodes(i1)
	  k2 = nodes(i2)
	  dx = xgv(k2) - xgv(k1)
	  dy = ygv(k2) - ygv(k1)
	  y = ( ygv(k2) + ygv(k1) ) / 2.
	  call make_dist_fact(isphe,y,xfact,yfact)
	  dx = xfact * dx
	  dy = yfact * dy
	  dd = sqrt( dx*dx + dy*dy )
	  dxy(1,i) = dx/dd
	  dxy(2,i) = dy/dd
	end do

	end

c************************************************************************

	subroutine proj_velox(mode,n,nodes,lnodes,ilhkv,dxy,vel,val
     +				,vmin,vmax,vhmin,vhmax)

c computes vel (vector) and val (scalar) for section
c in vel is normal/tangential/vertical velocity
c in val is scalar velocity used for overlay
c modes: 0=use normal vel   1=use tangent vel   as scalar velocity
 
	implicit none

	include 'param.h'

	integer mode			!what to use as scalar vel
	integer n			!total number of nodes
	integer nodes(n)		!node numbers
	integer lnodes(n)		!total nuber of layers
	integer ilhkv(1)		!number of layers in node
	real dxy(2,n)			!direction of line (for projection)
	real vel(3,0:2*nlvdim,n)	!projected velocities along line (ret)
	real val(0:2*nlvdim,n)		!scalar velocity for overlay (ret)
	real vmin,vmax			!min/max of scalar vel (ret)
	real vhmin,vhmax		!min/max of tangent vel (ret)

        real uprv(nlvdim,nkndim)
        common /uprv/uprv
        real vprv(nlvdim,nkndim)
        common /vprv/vprv
        real wprv(nlvdim,nkndim)
        common /wprv/wprv

	integer i,k,l,lmax,llayer
	real ut,un,w

	do i=1,n
	  k = nodes(i)
	  lmax = ilhkv(k)
	  do l=1,lmax
	    llayer = 2*l-1
	    ut =  dxy(1,i)*uprv(l,k) + dxy(2,i)*vprv(l,k)	!tangent
	    un = -dxy(2,i)*uprv(l,k) + dxy(1,i)*vprv(l,k)	!normal
	    w  =  wprv(l,k)
	    if( mode .eq. 0 ) then
	      val(llayer,i) = un		!use normal vel for overlay
	    else
	      val(llayer,i) = sqrt(ut**2+w**2)	!use tang vel for overlay
	    end if
	    vel(1,llayer,i) = un
	    vel(2,llayer,i) = ut
	    vel(3,llayer,i) = w
	  end do
	  call insert_between_layers(3,2*lmax,vel(1,0,i))
	  call insert_between_layers(1,2*lmax,val(0,i))

	  if( lnodes(i) .gt. lmax ) then	!more layers
	    do l=2*lmax+1,2*lnodes(i)
	      val(l,i) = val(2*lmax,i)
	      vel(1,l,i) = vel(1,2*lmax,i)
	      vel(2,l,i) = vel(2,2*lmax,i)
	      vel(3,l,i) = vel(3,2*lmax,i)
	    end do
	  end if

	end do

	vmin = val(0,1)
	vmax = val(0,1)
	vhmin = vel(2,0,1)
	vhmax = vel(2,0,1)
	do i=1,n
	  k = nodes(i)
	  lmax = ilhkv(k)
	  do l=0,2*lmax
	    vmin = min(vmin,val(l,i))
	    vmax = max(vmax,val(l,i))
	    vhmin = min(vhmin,vel(2,l,i))
	    vhmax = max(vhmax,vel(2,l,i))
	  end do
	end do

	end

c************************************************************************

	subroutine insert_between_layers(nd,nl,vals)

	implicit none

	integer nd,nl
	real vals(nd,0:2*nl)

	integer i,l

	do i=1,nd
	  vals(i,0) = vals(i,1)
	  do l=2,2*(nl-1),2
	    vals(i,l) = 0.5 * (vals(i,l-1)+vals(i,l+1))
	  end do
	  vals(i,2*nl) = vals(i,2*nl-1)
	end do

	end

c************************************************************************

	subroutine line_insert_scalars(n,nodes,lnodes,ilhkv,sv
     +					,values,vmin,vmax)

c inserts scalar values into matrix section

	implicit none

	include 'param.h'

	integer n
	integer nodes(n)
	integer lnodes(n)
	integer ilhkv(1)		!number of layers in node
	real sv(nlvdim,1)
	real values(0:2*nlvdim,n)
	real vmin,vmax

	integer i,k,lmax,l,llayer
	integer ipext

	vmin = 1.e+30
	vmax = -1.e+30

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
	  if( lnodes(i) .gt. lmax ) then	!more layers
	    do l=2*lmax+1,2*lnodes(i)
	      values(l,i) = values(2*lmax,i)
	    end do
	  end if
	  !write(78,*) k,ipext(k),lmax
	  !write(78,*) (sv(l,k),l=1,lmax)
	  !write(79,*) k,ipext(k),lmax,2*lmax
	  !write(79,*) (values(l,i),l=0,2*lmax)
	end do

	end

c************************************************************************

	subroutine line_find_min_max(n,nodes,helems,lelems,xgv,ygv
     +					,isphe,rlmax,rdmax,llmax,xy)

c finds length and max depth of line

	implicit none

	integer n
	integer nodes(n)
	real helems(2,n)	!depth in chosen elements
	integer lelems(n)	!maximum layer in element
	real xgv(1), ygv(1)	!coordinates
	integer isphe		!spherical coords?
	real rlmax,rdmax	!length and depth (return)
	integer llmax		!maximum layer
	real xy(n)		!distance of nodes from first node

	integer i,k1,k2
	real dx,dy
	real y,xfact,yfact

	rlmax = 0.
	rdmax = 0.
	llmax = 0
	xy(1) = 0.

	do i=2,n
	  k1 = nodes(i-1)
	  k2 = nodes(i)
	  dx = xgv(k2) - xgv(k1)
	  dy = ygv(k2) - ygv(k1)
	  y = (ygv(k1)+ygv(k2))/2.
	  call make_dist_fact(isphe,y,xfact,yfact)
	  dx = xfact * dx
	  dy = yfact * dy
	  rlmax = rlmax + sqrt( dx**2 + dy**2 )
	  xy(i) = rlmax
	  rdmax = max(rdmax,helems(1,i))
	  rdmax = max(rdmax,helems(2,i))
	  llmax = max(llmax,lelems(i))
	end do

	end

c************************************************************************

	subroutine make_dist_fact(isphe,y,xfact,yfact)

c computes factor for transformation from spherical to cartesian coordinates

	implicit none

	integer isphe			!spherical coordinates?
	real y				!average y coordinate
	real xfact,yfact		!factors for transformation

	real r,pi,rad
	parameter ( r = 6378206.4 , pi = 3.14159 , rad = pi/180. )

	if( isphe .ne. 0 ) then		!handle sperical coordinates
	  yfact = rad*r
	  xfact = yfact*cos(y*rad)
	else
	  xfact = 1.
	  yfact = 1.
	end if

	end

c************************************************************************

	subroutine line_find_elements(n,nodes,nlv,nen3v,hev,hm3v,hlv
     +					,elems,helems,lelems,lnodes)

c finds elements along line given by nodes
c
c deepest element is chosen

	implicit none

	integer n
	integer nodes(n)
	integer nlv		!number of layers
	integer nen3v(3,1)	!element index
	real hev(1)		!depth in elements
	real hm3v(3,1)		!depth in elements (on vertices)
	real hlv(1)		!layer structure
	integer elems(n)	!element number of chosen elements (return)
	real helems(2,n)	!depth in chosen elements (return)
	integer lelems(n)	!layers in element (return)
	integer lnodes(n)	!layers in node (return)

	logical bsigma,berror
	integer i,k1,k2,ie1,ie2,l
	integer ii,ie
	integer nsigma
	real hsigma
	real h
	integer ipext,ieext

	call get_sigma_info(nlv,nsigma,hsigma)
	bsigma = nsigma .gt. 0
	berror = .false.

c------------------------------------------------------------------
c set up depth structure in line (elements)
c------------------------------------------------------------------

	do i=2,n
	  k1 = nodes(i-1)
	  k2 = nodes(i)
	  call find_elems_to_segment(k1,k2,ie1,ie2)
	  if( ie2 .le. 0 ) then
	    if( ie1 .le. 0 ) then
	      write(6,*) '*** nodes: ',ipext(k1),ipext(k2)
	      berror =.true.
	    end if
	    elems(i) = ie1
	  else if( ie1 .le. 0 ) then
	    elems(i) = ie2
	  else
	    elems(i) = ie1
	    if( hev(ie2) .gt. hev(ie1) ) elems(i) = ie2
	  end if
	  ie = elems(i)
	  if( bsigma .and. hev(ie) .le. hsigma ) then
	    do ii=1,3
	      if( k1 .eq. nen3v(ii,ie) ) helems(1,i) = hm3v(ii,ie)
	      if( k2 .eq. nen3v(ii,ie) ) helems(2,i) = hm3v(ii,ie)
	    end do
	  else
	    helems(1,i) = hev(ie)
	    helems(2,i) = hev(ie)
	  end if
	  !write(6,*) 'line_find_elements: ',k1,k2,ie1,ie2,ie
	  !write(6,*) '                    ',hev(ie)
	  !write(6,*) '                    ',helems(1,i),helems(2,i)
	  !write(6,*) '                    ',ipext(k1),ipext(k2)
	  !write(6,*) '                    ',ieext(ie1),ieext(ie2)
	end do

	if( berror ) goto 99

c------------------------------------------------------------------
c compute total layers in line (elements)
c------------------------------------------------------------------

	do i=2,n
	  h = max(helems(1,i),helems(2,i))
	  l = max(1,nsigma)
	  if( bsigma .and. nlv .eq. nsigma ) then
	    lelems(i) = nsigma
	  else
	    do while( hlv(l) .lt. h ) 
	      l = l + 1
	    end do
	    lelems(i) = l
	  end if
	end do
	  
c------------------------------------------------------------------
c compute total layers in line (nodes)
c------------------------------------------------------------------

	lnodes(1) = lelems(2)
	do i=2,n-1
	  lnodes(i) = max(lelems(i),lelems(i+1))
	end do
	lnodes(n) = lelems(n)
	  
c------------------------------------------------------------------
c end of routine
c------------------------------------------------------------------

	return
   99	continue
	write(6,*) 'These nodes are not linked by any element'
	write(6,*) 'please revise your node list'
	stop 'error stop line_find_elements: cannot find elements'
	end

c************************************************************************

	subroutine line_read_nodes(file,ndim,n,nodes)

c reads line given by nodes

	implicit none

	character*(*) file
	integer ndim
	integer n
	integer nodes(ndim)

	logical berror
	integer iunit,nn,k
	integer ifileo,ipint

	n = 0
	berror = .false.

	iunit = ifileo(0,file,'form','old')
	if( iunit .le. 0 ) goto 99 

    1	continue
	  read(iunit,*,end=2) nn
	  if( nn .le. 0 ) goto 2
	  n = n + 1
	  if( n .gt. ndim ) stop 'error stop line_read_nodes: ndim'
	  k = ipint(nn)
	  if( k .le. 0 ) then
	    write(6,*) 'no such node: ',nn
	    berror = .true.
	  end if
	  nodes(n) = k
	goto 1
    2	continue

	write(6,*) 'finished reading line section: ',n

	close(iunit)

	if( berror ) goto 98

	return
   98	continue
	write(6,*) 'Some nodes are not existing'
	write(6,*) 'please revise your node list'
	stop 'error stop line_read_nodes: no such nodes'
   99	continue
	write(6,*) 'Cannot open section file: '
	write(6,*) file
	stop 'error stop line_read_nodes: no such file'
	end

c************************************************************************

	subroutine set_max_dep_lay(nlv,hlv,rdmax,llmax,hvmax,lvmax)

c sets hvmax and lvmax

	implicit none

	integer nlv
	real hlv(1)
	real rdmax
	integer llmax
	real hvmax
	integer lvmax

	integer l

        if( hvmax .gt. 0. .and. lvmax .gt. 0 ) then
          write(6,*) 'hvmax,lvmax: ',hvmax,lvmax
          write(6,*) 'only one of hvmax - lvmax can be given'
          stop 'error stop set_max_dep_lay: hvmax/lvmax'
        end if

        if( hvmax .gt. 0. ) then	!works also for sigma layers -> max
          do l=1,nlv
            if( hlv(l) .gt. hvmax ) goto 7
          end do
    7     continue
          lvmax = l             !could be one higher than max layer
          llmax = lvmax
          rdmax = hvmax
        else if( lvmax .gt. 0 ) then
          if( lvmax .gt. nlv .and. nlv .gt. 1 ) then	!more than available
            hvmax = hlv(nlv) + hlv(nlv) - hlv(nlv-1)	!add last layer
          else
            hvmax = hlv(lvmax)
          end if
          llmax = lvmax
          rdmax = hvmax
        else
          hvmax = rdmax
          lvmax = llmax
        end if

	end

c************************************************************************

	subroutine find_nc(x,nc)

c finds first significant position of x after decimal point - returns in nc

c 0.3 -> 2    0.03 -> 3       3 -> 0   30 -> 0

	implicit none

	real x
	integer nc

	real y

	nc = -1
	nc = +1
	y = abs(x)
	!write(6,*) 'find_nc: ',x,y,nc
	if( y .ge. 1. .or. y .eq. 0. ) return
	!write(6,*) 'find_nc: ',x,y,nc

	nc = 0
	do while( y .lt. 1. .or. nc .gt. 100 )
	  y = y * 10.
	  nc = nc + 1
	end do

	if( nc .gt. 100 ) stop 'error stop find_nc: too many iterations'

	end

c************************************************************************

	subroutine get_vel_unit(fact,string)

	implicit none

	real fact
	character*(*) string

	if( fact .eq. 1. ) then
	  string = 'm/s'
	else if( fact .eq. 100. ) then
	  string = 'cm/s'
	else if( fact .eq. 1000. ) then
	  string = 'mm/s'
	else
	  write(6,*) 'Cannot determine unit for fact = ',fact
	  string = ' '
	end if

	end

c************************************************************************

	subroutine make_segment_depth(ivert,nlv,hdep,hvmax,hlv,ya)

	implicit none

	integer ivert
	integer nlv
	real hdep(2)		!depth on two nodes
	real hvmax		!maximum depth
	real hlv(1)
	real ya(2,0:1)		!bottom depth of layers

	logical blayer,blog,bsigma
	integer i,l
	integer nlvaux,nsigma
	real hsigma
	real hd,h,hsig
	real hlog

	call get_sigma_info(nlvaux,nsigma,hsigma)
	!write(6,*) 'seg_depth sigma: ',nlvaux,nsigma,hsigma

	blayer = abs(ivert) .eq. 1
	blog = ivert .eq. 2

	do i=1,2
	  hd = hdep(i)
	  hsig = min(hsigma,hd)
	  ya(i,0) = 0.
	  do l=1,nlv
	    h = hlv(l)
	    bsigma = l .le. nsigma
	    if( bsigma .and. l .eq. nsigma ) h = -1.
	    if( blayer ) then
	      if( bsigma .and. ivert .eq. -1 ) then
		h = -h * l
	      else
		h = l
	      end if
	    else if( bsigma ) then
	      h = -h * hsig
	    else if( l .eq. nlv ) then
	      h = hd
	    end if
	    if( blog ) h = hlog(h,hvmax)
	    if( .not. blayer .and. h .gt. hvmax ) h = hvmax
	    ya(i,l) = -h
	  end do
	end do

	end

c************************************************************************

	subroutine trapez(x1,yt1,yb1,x2,yt2,yb2)

	implicit none

	real x1,yt1,yb1,x2,yt2,yb2

	call qline(x1,yt1,x1,yb1)
	call qline(x1,yb1,x2,yb2)
	call qline(x2,yb2,x2,yt2)
	call qline(x2,yt2,x1,yt1)

	end

c************************************************************************

