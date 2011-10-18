c
c $Id: supsim.f,v 1.23 2009-11-18 17:16:00 georg Exp $
c
c utility routines for plotsim
c
c contents :
c
c subroutine plonos(type)
c subroutine mkvarline(ivar,line)       makes variable description in nos file
c subroutine plozet
c subroutine plobar			plots barene
c subroutine plosim(bvel)
c subroutine ploval(nkn,parray,title)
c subroutine plozbar
c subroutine plovel(ivel)
c subroutine plo2vel(bvel)
c subroutine plo3vel(bvel)
c subroutine por2vel(n,up,vp,uv,vv,ht)
c subroutine plobas
c
c subroutine elkdep(nkn,hkv)
c subroutine eltype0(nel,iarv)
c subroutine eldepth(nel,hm3v,auxv,color,title)
c subroutine eltype(nel,iarv,color,title)
c
c subroutine plotel( ie , color )
c subroutine basplc( nel , color , line , descrr )
c
c subroutine bnd2val(a,val)			sets values on boundary to val
c subroutine moduv(u,v,uv,n,uvmax)		computes modulus and maximum
c subroutine normuv(uvnv,vvnv,vv,n)		normalize horizontal velocities
c subroutine intpn(vev,vnv,bwater)		compute nodal values vnv()
c subroutine apply_dry_mask(bmask,value,n,flag)	aplies mask to nodes
c
c revision log :
c
c 30.01.2002	ggu	better debug info, allow for vertical plotting
c 30.01.2002	ggu	new routines mkvarline, aspecial, traxy, rdspecial
c 05.10.2003	ggu	changes to allow for automatical plotting of vector
c 31.10.2003	ggu	subroutine plowind()
c 17.03.2003	ggu	use okvar to decide which var to plot (0->all)
c 05.10.2004	ggu	use znv instead xv
c 05.10.2004	ggu	inorm instead ivert, overlay over z/h, comp_scale()
c 16.12.2004	ggu	plot also regular net
c 04.03.2005	ggu	grey plot over bathy
c 14.03.2007	ggu	wave plotting with plowave
c 17.09.2008	ggu	plot last layer
c 06.12.2008	ggu	new routine get_minmax_flag(), bvel -> ivel for velsh
c 09.01.2009	ggu	deleted traref, tramin
c 27.01.2009	ggu	bug in plonos (input variable ivar has been changed)
c 07.05.2009	ggu	bug fix for grid and color plotting (plobas)
c 14.09.2009	ggu	section plot for scalars in plonos
c 09.10.2009	ggu	new routine plopres() for atmos. pressure
c 13.10.2009	ggu	section plot for velocities in plosim
c 26.03.2010	ggu	section plot for velocities in plosim adapted
c 17.12.2010    ggu     substituted hv with hkv
c 31.03.2011    ggu     no plotting in dry nodes implemented - read fvl file
c 17.05.2011    ggu     in plobas may plot node and element numbers
c 12.07.2011    ggu     eliminated all references to out routines
c 31.08.2011    ggu     new eos plotting (pleos,ploeval)
c 07.10.2011    ggu&dbf error calling extelev with nkn, and not nel
c
c**********************************************************
c**********************************************************
c**********************************************************

	subroutine plonos(type,ivar_in)

c 3D concentrations

	implicit none

	include 'param.h'

	character*(*) type
	integer ivar_in

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real p3(nlvdim,1)
	common /p3/p3
	real parray(1)
	common /parray/parray
	real fvlv(nlvdim,1)	!finite volume
	common /fvlv/fvlv

        character*80 line
	integer nrec,it
	integer ivaria
	integer level,ivar,isect
	logical nosnext,oktime,okvar
	integer getlev,getvar,getisec

	nrec = 0
	level = getlev()
        isect = getisec()
	ivar = ivar_in

	call nosopen(type)
	call fvlopen('.fvl')
	call timeask
	if( ivar .gt. 0 ) then
	  call setvar(ivar)
	else
	  call askvar
          ivar = getvar()
	end if

        call mkvarline(ivar,line)

	do while( nosnext(it,ivaria,nlvdim,p3) )
	  nrec = nrec + 1
	  write(6,*) nrec,it,ivaria,ivar
	  call fvlnext(it,nlvdim,fvlv)
	  if( oktime(it) .and. okvar(ivaria) ) then
            if( isect .eq. 0 ) then
	      write(6,*) '..........horizontal plotting nodes'
	      call extnlev(level,nlvdim,nkn,p3,parray)
	      call prepsim
	      call ploval(nkn,parray,line)
            else
	      write(6,*) '..........vertical plotting nodes'
              call plot_sect(.false.,p3)
            end if
	  end if
	end do

	call nosclose

	end

c**********************************************************

	subroutine ploeos(type,ivar_in)

c 3D concentrations (element values)

	implicit none

	include 'param.h'

	character*(*) type
	integer ivar_in

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real p3(nlvdim,1)
	common /p3/p3
	real parray(1)
	common /parray/parray
	real fvlv(nlvdim,1)	!finite volume
	common /fvlv/fvlv

        character*80 line
	integer nrec,it
	integer ivaria
	integer level,ivar,isect
	logical eosnext,oktime,okvar
	integer getlev,getvar,getisec

	nrec = 0
	level = getlev()
        isect = getisec()
	ivar = ivar_in

	call eosopen(type)
	call timeask
	if( ivar .gt. 0 ) then
	  call setvar(ivar)
	else
	  call askvar
          ivar = getvar()
	end if

        call mkvarline(ivar,line)

	do while( eosnext(it,ivaria,nlvdim,p3) )
	  nrec = nrec + 1
	  write(6,*) nrec,it,ivaria,ivar
	  if( oktime(it) .and. okvar(ivaria) ) then
            if( isect .eq. 0 ) then
	      write(6,*) '..........horizontal plotting elements'
	      !call extelev(level,nlvdim,nkn,p3,parray)
	      call extelev(level,nlvdim,nel,p3,parray)
	      call prepsim
	      call ploeval(nel,parray,line)
            else
	      write(6,*) '..........no vertical plotting for elements'
            end if
	  end if
	end do

	call eosclose

	end

c**********************************************************

        subroutine mkvarline(ivar,line)

c makes variable description in nos file

        implicit none

        integer ivar            !variable number
        character*(*) line      !description on return

        integer iaux
        integer ialfa

        line = 'variable = '
        iaux = ialfa(float(ivar),line(12:),-1,-1)

        if( ivar .eq. 10 ) then
          line = 'generic tracer'
        else if( ivar .eq. 11 ) then
          line = 'salinity'
        else if( ivar .eq. 12 ) then
          line = 'temperature'
        else if( ivar .eq. 15 ) then
          line = 'oxygen'
        else if( ivar .eq. 18 ) then
          line = 'rms velocity'
        end if

        end

c**********************************************************

	subroutine plozet

	implicit none

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real znv(1)
	common /znv/znv

	integer nrec,it,k
	logical ousnext,oktime

	nrec = 0

	call ousopen
	call timeask

	do while( ousnext(it) )
	  nrec = nrec + 1
	  write(6,*) nrec,it
	  if( oktime(it) ) then
	    write(6,*) '..........plotting'
	    call prepsim
	    call ploval(nkn,znv,'water levels')
	  end if
	end do

	call ousclose

	end

c**********************************************************

	subroutine plobar

c plots barene

	implicit none

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real parray(1)
	common /parray/parray

	integer nrec,it,k
	logical ousnext,oktime

	nrec = 0

	call ousopen
	call timeask

	do while( ousnext(it) )
	  nrec = nrec + 1
	  write(6,*) nrec,it
	  if( oktime(it) ) then
	    write(6,*) '..........plotting'
	    call plozbar
	  end if
	end do

	call ousclose

	end

c**********************************************************

	subroutine plopres

	implicit none

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real pres(1)
	common /pres/pres

	integer nrec,it
	logical windnext,oktime

	nrec = 0

	call windopen
	call timeask

	do while( windnext(it) )
	  nrec = nrec + 1
	  write(6,*) nrec,it
	  if( oktime(it) ) then
	    write(6,*) '..........plotting'
	    call resetsim
	    call ploval(nkn,pres,'atmospheric pressure')
	  end if
	end do

	call windclose

	end

c**********************************************************

	subroutine plowind

	implicit none

	integer nrec,it
        integer ivel
	logical windnext,oktime

	nrec = 0
        ivel = 3	!wind

	call windopen
	call timeask

	do while( windnext(it) )
	  nrec = nrec + 1
	  write(6,*) nrec,it
	  if( oktime(it) ) then
	    write(6,*) '..........plotting'
	    call resetsim
	    call plovel(ivel)
	  end if
	end do

	call windclose

	end

c**********************************************************

	subroutine plowave

	implicit none

	integer nrec,it
        integer ivel
	logical wavenext,oktime

	nrec = 0
        ivel = 4	!wave

	call waveopen
	call timeask

	do while( wavenext(it) )
	  nrec = nrec + 1
	  write(6,*) nrec,it
	  if( oktime(it) ) then
	    write(6,*) '..........plotting'
	    call resetsim
	    call plovel(ivel)
	  end if
	end do

	call waveclose

	end

c**********************************************************

	subroutine plosim(bvel)

	implicit none

	include 'param.h'

	logical bvel

	integer nrec,it
        integer ivel,isect

	real p3(nlvdim,1)
	common /p3/p3

	logical velnext,oktime
	integer getisec

        if( bvel ) then
          ivel = 1
        else
          ivel = 2
        end if
        isect = getisec()

	nrec = 0

	call velopen
	call timeask

	do while( velnext(it) )
	  nrec = nrec + 1
	  write(6,*) nrec,it
	  if( oktime(it) ) then
            if( isect .eq. 0 ) then
	      write(6,*) '..........horizontal plotting'
	      call prepsim
	      call plovel(ivel)
            else
	      write(6,*) '..........vertical plotting'
	      call prepare_vel(p3)
              call plot_sect(.true.,p3)
            end if
	  end if
	end do

	call velclose

	end

c**********************************************************

	subroutine ploval(nkn,parray,title)

c plots node values

	implicit none

	integer nkn
	real parray(1)
        character*(*) title

        logical bwater(1)		!mask for elements
        common /bwater/bwater
        logical bkwater(1)		!mask for nodes
        common /bkwater/bkwater

	real pmin,pmax,flag
	real getpar
        integer gettime

        call get_flag(flag)

	call qstart
        call annotes(gettime(),title)
	call bash(0)

	call get_minmax_flag(parray,nkn,pmin,pmax)
        write(6,*) 'min/max: ',nkn,pmin,pmax
	call apply_dry_mask(bkwater,parray,nkn,flag)	!flags to dry nodes
	call colauto(pmin,pmax)

        call qcomm('Plotting isolines')
        call isoline(parray,nkn,0.,2)
        call colsh

	!call bash(3)	! overlays grid (for debug)

	call bash(2)
	call qend

	end

c**********************************************************

	subroutine ploeval(nel,parray,title)

c plots element values

	implicit none

	integer nel
	real parray(1)
        character*(*) title

        logical bwater(1)		!mask for elements
        common /bwater/bwater
        logical bkwater(1)		!mask for nodes
        common /bkwater/bkwater

	real pmin,pmax,flag
	real getpar
        integer gettime

        call get_flag(flag)

	call qstart
        call annotes(gettime(),title)
	call bash(0)

	call get_minmax_flag(parray,nel,pmin,pmax)
        write(6,*) 'min/max: ',nel,pmin,pmax
	call apply_dry_mask(bwater,parray,nel,flag)	!flags to dry nodes
	call colauto(pmin,pmax)

        call qcomm('Plotting element values')
        call isoline(parray,nel,0.,3)			!plot on elements
        call colsh

	!call bash(3)	! overlays grid (for debug)

	call bash(2)
	call qend

	end

c**********************************************************

	subroutine plozbar

	implicit none

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real zenv(3,1)
        common /zenv/zenv
        real znv(1)
        common /znv/znv
        integer nen3v(3,1)
        common /nen3v/nen3v

	logical bdry
	integer ie,ii,k

	call qstart
	call bash(0)

        call qcomm('Plotting barene')
	do ie=1,nel
	  bdry = .false.
	  do ii=1,3
	    k = nen3v(ii,ie)
	    if( zenv(ii,ie) .ne. znv(k) ) bdry = .true.
	  end do
	  if( bdry ) then
	    call plotel(ie,0.1)
	  else
	    call plotel(ie,0.8)
	  end if
	end do

	call bash(2)
	call qend

	end

c**********************************************************

	subroutine plovel(ivel)

	implicit none

	integer ivel

	logical is2d

	if( is2d() ) then
	  write(6,*) 'plotting 2d...'
	  call plo2vel(ivel,'2D ')
	else
	  write(6,*) 'plotting 3d...'
	  call plo3vel(ivel)
	end if

	end

c**********************************************************

	subroutine plo2vel(ivel,title)

c plots entity with vectors
c
c ivel = 1	velocities
c ivel = 2	transports
c ivel = 3	wind
c ivel = 4	waves

	implicit none

	integer ivel			!what to plot
        character*(*) title

	integer nxdim,nydim
	parameter( nxdim = 300 , nydim = 300 )

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	real xgv(1), ygv(1)
	common /xgv/xgv, /ygv/ygv
        real usnv(1), vsnv(1)
        common /usnv/usnv, /vsnv/vsnv
        real wsnv(1)
        common /wsnv/wsnv
        real hetv(1)
        common /hetv/hetv
        real uvnv(1), vvnv(1)
        common /uvnv/uvnv, /vvnv/vvnv
	real vev(1)
	common /vev/vev
	real v1v(1)
	common /v1v/v1v
	real uv(1), vv(1)
	common /uv/uv, /vv/vv

        real hkv(1)
        common /hkv/hkv
        real znv(1)
        common /znv/znv

        logical bwater(1)		!mask for elements
        common /bwater/bwater
        logical bkwater(1)		!mask for nodes
        common /bkwater/bkwater

	real ureg(nxdim,nydim)
	real vreg(nxdim,nydim)

	logical boverl,bnode,bbound,bnorm
        logical bspecial
	logical bdebug
	logical bregular
	logical bvel,btrans,bwind,bwave
        character*80 line
        character*80 spcvel
	integer ie,k
	integer i,j
	integer level
	integer ioverl,inorm
	real u,v,xm,ym
	real x0,y0,dx,dy,flag
	real ut,vt
	real href,scale
	real velref,valref
	real velmin,valmin
	real typls,typlsf
	real color
	real xmin,ymin,xmax,ymax
	real pmin,pmax
	real uvmax,uvmed
	real typsca
	real fscale
	integer nx,ny

	real getpar,getcol
	integer gettime,getlev

	bdebug = .true.
	bdebug = .false.

c------------------------------------------------------------------
c set up parameters
c------------------------------------------------------------------

c	-----------------------------------------------------------
c	locally defined parameters (please change here)
c	-----------------------------------------------------------

	bnode = .false.		!plot arrows on nodes
	bnode = .true.		!plot arrows on nodes
	bbound = .true. 	!plot arrows on boundary nodes
	bbound = .false. 	!plot arrows on boundary nodes
	fscale = 3.		!extra empirical factor for computing arrows

c	-----------------------------------------------------------
c	global parameters from STR file
c	-----------------------------------------------------------

        bvel   = ivel .eq. 1
        btrans = ivel .eq. 2
        bwind  = ivel .eq. 3
        bwave  = ivel .eq. 4

	velref = getpar('velref')
	velmin = getpar('velmin')

	typlsf = getpar('typlsf')               !additional factor for length
	href = getpar('href')
        call getfnm('spcvel',spcvel)

	ioverl = nint(getpar('ioverl')) 	!mode of overlay plotting
	inorm  = nint(getpar('inorm'))	 	!normalization of arrows

	bspecial = spcvel .ne. ' '
	if( bspecial ) ioverl = 0

	boverl = ioverl .ne. 0			!overlay color
	bnorm  =  inorm .eq. 1			!normalize vectors

c	-----------------------------------------------------------
c	last initializations and start plotting
c	-----------------------------------------------------------

	call mkht(hetv,href)

	call qstart

        line = title
        if( bwind ) then
          line(4:) = 'wind velocity'
        else if( bwave ) then
          line(4:) = 'waves'
        else if( bvel ) then
          line(4:) = 'velocity'
        else if( btrans ) then
          line(4:) = 'transport'
        else
          stop 'error stop plo2vel: internal error (3)'
        end if
	call annotes(gettime(),line)
	call bash(0)

c------------------------------------------------------------------
c see if regular grid
c------------------------------------------------------------------

	call getgeo(x0,y0,dx,dy,flag)
	bregular = dx .gt. 0. .and. dy .gt. 0.
	if( bregular ) then
	  x0 = x0 - dx
	  y0 = y0 - dy
	  call getbas(xmin,ymin,xmax,ymax)
	  nx = nint((xmax-xmin)/dx)
	  ny = nint((ymax-ymin)/dy)
	  !write(6,*) 'ggu: x0,y0,dx,dy ',x0,y0,dx,dy,nx,ny
	  !write(6,*) 'ggu: ',xmin,ymin,xmax,ymax
	  if( nx .gt. nxdim .or. ny .gt. nydim ) then
	    write(6,*) '*** Warning: regular plot may not be complete'
	    write(6,*) ' nxdim, nydim : ',nxdim,nydim
	    write(6,*) ' nx   , ny    : ',nx   ,ny   
	  end if
	end if

c------------------------------------------------------------------
c prepare for velocity or transport
c------------------------------------------------------------------

	valref = velref
	valmin = velmin

	if( bvel ) then
	  call por2vel(nel,usnv,vsnv,uvnv,vvnv,hetv)
        else if( btrans ) then
	  do ie=1,nel
	    uvnv(ie) = usnv(ie)
	    vvnv(ie) = vsnv(ie)
	  end do
        else if( bwind .or. bwave ) then
	  !ok
	else
	  stop 'error stop plo2vel: internal error (77)'
	end if

c------------------------------------------------------------------
c compute values on nodes -> 0 for dry nodes
c------------------------------------------------------------------

        if( bwind .or. bwave ) then
	  call intpe(uv,uvnv,bwater)
	  call intpe(vv,vvnv,bwater)
        else
	  call intpn(uvnv,uv,bwater)
	  call intpn(vvnv,vv,bwater)
        end if

c------------------------------------------------------------------
c regular interpolation
c------------------------------------------------------------------

	if( bregular ) then
	  call av2amk(bwater,uv,ureg,nxdim,nydim)
	  call av2amk(bwater,vv,vreg,nxdim,nydim)
	end if

c------------------------------------------------------------------
c compute modulus and maximum of velocity
c------------------------------------------------------------------

	call moduv(uv,vv,v1v,nkn,uvmax,uvmed)	!compute modulus, max, aver

	if( .not. bbound ) then		!set boundary vectors to 0
	  call bnd2val(uv,0.)	
	  call bnd2val(vv,0.)	
	end if

c	v1v is used for overlay plot

c------------------------------------------------------------------
c underlying color 
c -> create vev/v1v (v1v is modulus of vel/tra) and plot
c------------------------------------------------------------------

	if( boverl ) then
	  if( ioverl .eq. 2 ) then		!vertical velocities
	    do k=1,nkn
	      v1v(k) = wsnv(k)
	    end do
          else if( ioverl .eq. 3 ) then		!water levels
	    do k=1,nkn
	      v1v(k) = znv(k)
	    end do
          else if( ioverl .eq. 4 ) then		!bathymetry
	    do k=1,nkn
	      v1v(k) = hkv(k)
	    end do
          else if( ioverl .ne. 1 ) then		!not horizontal velocity
            write(6,*) 'ioverl = ',ioverl
            stop 'error stop plo2vel: value not allowed for ioverl'
	  end if
	  call mima(v1v,nkn,pmin,pmax)
	  call apply_dry_mask(bkwater,v1v,nkn,flag)
	  call colauto(pmin,pmax)
	  write(6,*) 'plotting overlay color... ',pmin,pmax
          call qcomm('Plotting overlay')
          call isoline(v1v,nkn,0.,2)
	end if

c------------------------------------------------------------------
c get scale
c------------------------------------------------------------------

	typls  = getpar('typls')		!is available only now...
	if( valref .le. 0. ) then		!use maximum
	    if( bnorm ) then
	      valref = uvmax
	    else
	      valref = uvmax / fscale		!use empirical smaller value
	    end if
	end if
	typsca = typls * typlsf
	write(6,*) 'scale (type): ',typls,typlsf,ioverl
	write(6,*) 'scale (value): ',valref,valmin,uvmax,uvmed

c------------------------------------------------------------------
c plotting of arrow
c------------------------------------------------------------------

	call qgray(0.)

	if( bspecial ) then
	  scale = typsca / valref
          call aspecial(spcvel,xgv,ygv,uv,vv,scale)
        else if( bregular ) then
	  write(6,*) 'regular grid: ',dx,dy,nx,ny
	  do j=1,nydim				!FIXME -> maybe nx,ny is enough
	    do i=1,nxdim
	      if( ureg(i,j) .ne. flag ) then	!can plot
		xm = x0 + i * dx
		ym = y0 + j * dy
		u = ureg(i,j)
		v = vreg(i,j)
		call comp_scale(inorm,typsca,valmin,valref,u,v,scale)
	        call pfeil(xm,ym,ureg(i,j),vreg(i,j),scale)
	      end if
	    end do
	  end do
	else if( bnode ) then
	  do k=1,nkn
	    u = uv(k)
	    v = vv(k)
	    call comp_scale(inorm,typsca,valmin,valref,u,v,scale)
	    call pfeil(xgv(k),ygv(k),u,v,scale)
	  end do
	else
	  do ie=1,nel
	    call baric(ie,xm,ym)
	    u = uvnv(ie)
	    v = vvnv(ie)
	    call comp_scale(inorm,typsca,valmin,valref,u,v,scale)
	    call pfeil(xm,ym,u,v,scale)
	  end do
	end if

c------------------------------------------------------------------
c plotting of color bar or scale arrow
c------------------------------------------------------------------

	if( boverl ) then			!plot color scale
	  call colsh
	end if
	if( .not. bnorm ) then	                !vel arrows are still in scale
	  scale = typsca / valref		!arrows are in scale
	  call velsh(ivel,scale,valref)
	end if

c------------------------------------------------------------------
c end of plot
c------------------------------------------------------------------

	call bash(2)
	call qend

c------------------------------------------------------------------
c end of routine
c------------------------------------------------------------------

	end

c**********************************************************

	subroutine plo3vel(ivel)

	implicit none

	integer ivel

c	integer nlvdim
c	parameter ( nlvdim = 14 )

	include 'param.h'

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        integer nlvdi,nlv
        common /level/ nlvdi,nlv

        integer ilhv(1)
        common /ilhv/ilhv
        integer ilhkv(1)
        common /ilhkv/ilhkv

        real utlnv(nlvdim,1)
        common /utlnv/utlnv
        real vtlnv(nlvdim,1)
        common /vtlnv/vtlnv
        real wlnv(0:nlvdim,1)
        common /wlnv/wlnv
        real usnv(1), vsnv(1)
        common /usnv/usnv, /vsnv/vsnv
        real wsnv(1)
        common /wsnv/wsnv

	integer ie,ii,k,l
	integer level,lmax,lact
        real utot,vtot
	logical bdebug

	integer getlev

	bdebug = .true.
	bdebug = .false.

	if( nlvdi .ne. nlvdim ) stop 'error stop plo3vel: nlvdim'

c make vertical velocities

	call vertical

c handle level

	level = getlev()

	if( level .lt. -1 ) stop 'internal error (1)'

	if( level .eq. 0 ) then		!barotropic
          do ie=1,nel
            utot = 0.
            vtot = 0.
            lmax = ilhv(ie)
            do l=1,lmax
              utot = utot + utlnv(l,ie)
              vtot = vtot + vtlnv(l,ie)
            end do
            usnv(ie) = utot
            vsnv(ie) = vtot
          end do
	  do k=1,nkn
	    wsnv(k) = wlnv(level,k)
	  end do
	else
	  do ie=1,nel
	    lmax = ilhv(ie)
	    lact = level
	    if( lact .eq. -1 ) lact = lmax

	    if( level .gt. lmax ) then	!no such layer
	      usnv(ie) = 0.
	      vsnv(ie) = 0.
	    else
	      usnv(ie) = utlnv(lact,ie)
	      vsnv(ie) = vtlnv(lact,ie)
	    end if
	  end do
	  do k=1,nkn
	    lmax = ilhkv(k)
	    lact = level
	    if( lact .eq. -1 ) lact = lmax
	    wsnv(k) = wlnv(lact,k)
	  end do
	end if

c debug

	if( bdebug ) then
	  ie = 100
	  write(6,*) 'debugging 3d routine...'
	  write(6,*) ie,level,ilhv(ie),usnv(ie),vsnv(ie)
	  do l=1,ilhv(ie)
	    write(6,*) l,utlnv(l,ie),vtlnv(l,ie)
	  end do
	end if

c call 2d routine

	call plo2vel(ivel,'3D ')
	
	end

c**********************************************************

	subroutine por2vel(n,up,vp,uv,vv,ht)

	implicit none

	integer n
	real up(1),vp(1),uv(1),vv(1),ht(1)

	integer i
	real r

	do i=1,n
	  if( ht(i) .gt. 0. ) then
	    r = 1. / ht(i)
	    uv(i) = up(i) * r
	    vv(i) = vp(i) * r
	  else
	    uv(i) = 0.
	    vv(i) = 0.
	  end if
	end do

	end

c**********************************************************

	subroutine plobas

	implicit none

	character*80 descrr
	common /descrr/ descrr

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	real hkv(1)
	common /hkv/hkv
	integer iarv(1)
	common /iarv/iarv

	logical bnumber
        real pmin,pmax

	bnumber = .true.	! plot node and element numbers
	bnumber = .false.

c only boundary line

	call qstart
	call bash(0)
	call bash(2)
	call qend

c grid

	call qstart
	call bash(0)
	call bash(1)
	call qend

c grid with gray

	call qstart
	call bash(0)
	call bash(3)
	call qend

c bathymetry (gray or color)

	call resetsim
	call ploval(nkn,hkv,'basin')

c bathymetry with grid

	call qstart
	call bash(0)

	call mima(hkv,nkn,pmin,pmax)
	call colauto(pmin,pmax)
        call isoline(hkv,nkn,0.,2)

        call qcomm('before basin overlay');
        call qlwidth(0.001)
	call bash(1)    !black grid
        call colsh
        call qcomm('after basin overlay');
        call qlwidth(0.01)
	call qend

c bathymetry with grid (gray)

	call qstart
	call bash(0)

	call mima(hkv,nkn,pmin,pmax)
	call colauto(pmin,pmax)
        call isoline(hkv,nkn,0.,2)

        call qcomm('before basin overlay');
        call qlwidth(0.001)
	call bash(3)    !gray grid
        call colsh
        call qcomm('after basin overlay');
        call qlwidth(0.01)
	call qend

c special

	call eltype0(nel,iarv)

c boundary line with net

	call qstart
	call bash(0)
	call bash(2)
        call reggrid(10,0.,0.5)
	call qend

c here only debug (node and element numbers)

	if( bnumber ) then

	call qstart
	call bash(0)
	call bash(3)
	call basin_number(1)
	call qend

	call qstart
	call bash(0)
	call bash(3)
	call basin_number(2)
	call qend

	call qstart
	call bash(0)
	call bash(3)
	call basin_number(-1)
	call qend

	call qstart
	call bash(0)
	call bash(3)
	call basin_number(-2)
	call qend

	end if

c end of routine

	end

c**************************************************************

	subroutine elkdep(nkn,hkv)

	implicit none

	integer nkn
	real hkv(1)

	call qstart

	call bash(0)

        call qcomm('Plotting isolines')
	call isoline(hkv,nkn,0.,2)
	call colsh

        call bash(2)

	call qend

	end

c**************************************************************

	subroutine eltype0(nel,iarv)

	implicit none

	integer nel
	integer iarv(1)

	integer ie,ia

	call qstart

	call bash(0)

        call qcomm('Plotting element types')

	do ie=1,nel
	  ia = iarv(ie)
c	  if( ia .eq. 1 ) then					!channels
c	  if( ia .ge. 3 .and. ia .le. 5 .or. ia .eq. 9 ) then	!inlets
c	  if( ia .ne. 0 ) then					!special type
c	  if( ia .ge. 3 .and. ia .le. 5 .or. ia .eq. 9 ) then	!inlets
          if( ia .ge. 12 ) then                                 !valli pesca
            call plotel(ie,0.4)
          end if
	end do

        call bash(2)

	call qend

	end

c**************************************************************

	subroutine eldepth(nel,hm3v,auxv,color,title)

	implicit none

	integer nel
	real hm3v(3,1)
	real auxv(1)
	real color(1)
	character*(*) title

	integer ii,ie
	real hm

	do ie=1,nel
	  hm = 0.
	  do ii=1,3
	    hm = hm + hm3v(ii,ie)
	  end do
	  auxv(ie) = hm / 3.
	end do

	do ie=1,nel
	  hm = auxv(ie)
	  color(ie) = 1.
	  if( hm .le. 0. ) color(ie) = 0.7
	end do
	call basplc(nel,color,'depth = 0.',title)

	do ie=1,nel
	  hm = auxv(ie)
	  color(ie) = 1.
	  if( hm .le. 0.2 ) color(ie) = 0.7
	end do
	call basplc(nel,color,'depth = 0.2',title)

	do ie=1,nel
	  hm = auxv(ie)
	  color(ie) = 1.
	  if( hm .le. 0.4) color(ie) = 0.7
	end do
	call basplc(nel,color,'depth = 0.4',title)

	do ie=1,nel
	  hm = auxv(ie)
	  color(ie) = 1.
	  if( hm .le. 0.6 ) color(ie) = 0.7
	end do
	call basplc(nel,color,'depth = 0.6',title)

	end

c**************************************************************

	subroutine eltype(nel,iarv,color,title)

	implicit none

	integer nel,iarv(1)
	real color(1)
	character*(*) title

	integer ia,ie

	do ie=1,nel
	  ia = iarv(ie)
	  if( ia .eq. 0 ) ia = 1	!treat barene and bassifondi the same
	  color(ie) = 0.1 * ( 10 - ia )
	end do
	call basplc(nel,color,'type of elements (1=0)',title)

	do ie=1,nel
	  ia = iarv(ie)
	  if( ia .eq. 0 ) ia = 1	!treat barene and bassifondi the same
	  if( ia .ne. 1 ) ia = 5
	  color(ie) = 0.1 * ( 10 - ia )
	end do
	call basplc(nel,color,'channels and lagoon',title)

	do ie=1,nel
	  ia = iarv(ie)
	  if( ia .gt. 9 ) then
		ia = 1
	  else
		ia = 5
	  end if
	  color(ie) = 0.1 * ( 10 - ia )
	end do
	call basplc(nel,color,'strange types and islands',title)

	do ie=1,nel
	  ia = iarv(ie)
	  if( ia .eq. 1 ) then
		ia = 5
	  else
		ia = 1
	  end if
	  color(ie) = 0.1 * ( 10 - ia )
	end do
	call basplc(nel,color,'barene',title)

	end

c**************************************************************

	subroutine plotel( ie , color )

	implicit none

	integer ie
	real color

        real xgv(1), ygv(1)
        integer nen3v(3,1)
        common /xgv/xgv, /ygv/ygv
        common /nen3v/nen3v

	integer ii,k
	real xp(3),yp(3)

	if( color .lt. 0. ) return

	do ii=1,3
	  k = nen3v(ii,ie)
	  xp(ii) = xgv(k)
	  yp(ii) = ygv(k)
	end do

	call qsetc(color)
	call qafill(3,xp,yp)

	end

c**************************************************************

	subroutine basplc( nel , color , line , descrr )

	implicit none

	integer nel
	real color(nel)
	character*(*) line,descrr

	integer ie

	call qstart

        call basin(0)
	call frame(0)

        call qcomm('Plotting elements')
	do ie=1,nel
	  call plotel(ie,color(ie))
	end do

        call qcomm('Plotting basin')
        call basin(2)

	call qend

	end

c*****************************************************************

	subroutine bnd2val(a,val)	

c sets values on boundary to val

	real a(1)
	real val

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	integer kantv(2,1)
	common /kantv/kantv

	integer k

	logical isbound
	isbound(k) = kantv(1,k) .gt. 0

	do k=1,nkn
	  if( isbound(k) ) then
	    a(k) = val
	  end if
	end do

	end

c*****************************************************************

	subroutine comp_scale(mode,amax,vmin,vmax,u,v,scale)

c computes scale to apply for plot
c
c mode:
c       0       uses real length
c       1       normalizes vector (all vectors have same length)
c       2       uses linear scale
c       3       uses logarithmic scale

	implicit none

	integer mode		!type of scaling
	real amax		!velocity independent scale
	real vmin		!minimum velocity
	real vmax		!reference velocity -> gives amax length
	real u,v		!velocity values
	real scale		!computed scale (return)

	real vmod,aux

	vmod = sqrt(u*u+v*v)

	if( vmod .le. vmin ) then
	  scale = 0.
	  return
	end if

	if( mode .eq. 0 ) then
	  scale = amax / vmax
	else if( mode .eq. 1 ) then
	  scale = amax / vmod
	else if( mode .eq. 2 ) then
	  aux = (vmod-vmin) / (vmax-vmin)                       !aux in [0-1]
	  scale = ( amax / vmod ) * aux
	else if( mode .eq. 3 ) then
	  aux = log10( 1 + 9 * (vmod-vmin)/(vmax-vmin) )        !aux in [0-1]
	  scale = ( amax / vmod ) * aux
	else
	  write(6,*) 'mode = ',mode
	  stop 'error stop comp_scale: no such mode'
	end if

	end

c*****************************************************************

	subroutine moduv(u,v,uv,n,uvmax,uvmed)

c computes modulus for velocity vector and maximum and average

	implicit none

	integer n
	real u(1), v(1), uv(1)
	real uvmax
	real uvmed

	integer i
	real um,vm,rmod,rmax,rmed

	rmax = 0.
	rmed = 0.

	do i=1,n
	  um = u(i)
	  vm = v(i)
	  rmod = sqrt( um*um + vm*vm )
	  rmax = max(rmax,rmod)
	  rmed = rmed + rmod
	  uv(i) = rmod
	end do

	uvmax = rmax
	uvmed = rmed / n

	end

c*****************************************************************

	subroutine normuv(uvnv,vvnv,vv,n)

c normalize horizontal velocities

	implicit none

	integer n
	real uvnv(1), vvnv(1), vv(1)

	integer i
	real rmod

	do i=1,n
	  rmod = vv(i)
	  if( rmod .gt. 0. ) then
	    uvnv(i) = uvnv(i) / rmod
	    vvnv(i) = vvnv(i) / rmod
	  end if
	end do

	end

c*****************************************************************

	subroutine intpe(vnv,vev,bwater)

c compute elemental values vev()

	implicit none

        real vnv(1)
	real vev(1)
	logical bwater(1)

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	integer nen3v(3,1)
	common /nen3v/nen3v

	integer ie,ii,k
	real sum

	do ie=1,nel
	  if( bwater(ie) ) then
            sum = 0.
	    do ii=1,3
	      k = nen3v(ii,ie)
	      sum = sum + vnv(k)
	    end do
            vev(ie) = sum / 3.
          else
            vev(ie) = 0.
	  end if
	end do

	end

c*****************************************************************

	subroutine intpn(vev,vnv,bwater)

c compute nodal values vnv()

	implicit none

	real vev(1)
        real vnv(1)
	logical bwater(1)

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        real v2v(1)
        common /v2v/v2v
	integer nen3v(3,1)
	common /nen3v/nen3v

	integer ie,ii,k
	real r

	do k=1,nkn
	  vnv(k) = 0.
	  v2v(k) = 0.
	end do

	do ie=1,nel
	 if( bwater(ie) ) then
	  do ii=1,3
	    k = nen3v(ii,ie)
	    vnv(k) = vnv(k) + vev(ie)
	    v2v(k) = v2v(k) + 1.
	  end do
	 end if
	end do

	do k=1,nkn
	  r = v2v(k)
	  if( r .gt. 0. ) then
	    vnv(k) = vnv(k) / v2v(k)
	  else
	    vnv(k) = 0.
	  end if
	end do

	end

c*****************************************************************

	subroutine apply_dry_mask(bmask,value,n,flag)

c aplies mask to nodal value

	implicit none

	integer n
	real flag
	real value(1)
	logical bmask(1)

	integer i

	do i=1,n
	  if( .not. bmask(i) ) value(i) = flag
	end do

	end

c*****************************************************************

        subroutine aspecial(file,xgv,ygv,uv,vv,scale)

c plots arrow on special points

        implicit none

        character*80 file
        real scale
        real xgv(1)
        real ygv(1)
        real uv(1)
        real vv(1)

        integer ndim
        parameter (ndim=100)

        integer ns
        integer ks(ndim)
        real us(ndim),vs(ndim)
        real xs(ndim),ys(ndim)
        save ns,ks,us,vs,xs,ys

        logical berror
        integer i,k
        real xn,yn,dist

        integer icall
        save icall
        data icall / 0 /

        if( icall .lt. 0 ) return
        if( icall .eq. 0 ) then
          write(6,*) 'reading special velocity file...'
          write(6,*) file(1:79)
          ns = ndim
          call rdspecial(file,ns,ks,us,vs)
          call n2int(ns,ks,berror)
          if( berror ) stop 'error stop aspecial: no such node'
          do i=1,ns
            xs(i) = xgv(ks(i))
            ys(i) = ygv(ks(i))
          end do
          icall = 1
        end if

        !call qhue(1.)
        call qgray(0.5)
        dist = 50.
	do i=1,ns
          k = ks(i)
          call traxy(xgv(k),ygv(k),us(i),vs(i),dist,xn,yn)
	  call pfeil(xn,yn,us(i),vs(i),scale)
	end do

        call qgray(0.)
	do i=1,ns
          k = ks(i)
	  call pfeil(xgv(k),ygv(k),uv(k),vv(k),scale)
	end do

        call qgray(0.)

        end

c*****************************************************************

        subroutine traxy(xo,yo,u,v,dist,xn,yn)

c translates x/y to new coordinate to the left of direction (u,v) by dist

          sp = sqrt(u**2+v**2)
          un = - v/sp
          vn = + u/sp
          xn = xo + dist * un
          yn = yo + dist * vn

          end

c*****************************************************************

        subroutine rdspecial(file,n,ks,us,vs)

        implicit none

        character*80 file
        integer n
        integer ks(1)
        real us(1)
        real vs(1)

        integer ndim,i,k
        real pi,rad
        real sp,dir
        real f(10)
        character*80 line

        integer iscanf

        ndim = n
        n = 0
        pi = 4. * atan(1.)
        rad = pi / 180.

        open(9,file=file,form='formatted',status='old')

    1   continue
          !read(9,*,end=2) i,k,sp,dir
          read(9,'(a)',end=2) line
          if( line(1:1) .eq. '#' ) goto 1       !comment
          i = iscanf(line,f,10)
          if( i .eq. 0 ) goto 1                 !empty
          if( i .lt. 4 ) goto 99
          n = n + 1
          if( n .gt. ndim ) stop 'error stop rdspecial: ndim'
          k = nint(f(2))
          if( i .eq. 4 ) then
            sp = f(3)
            dir = f(4)
          else if( i .ge. 6 ) then
            sp = f(5)
            dir = f(6)
          end if
          ks(n) = k
          write(6,*) i,k,sp,dir
          sp = sp / 100.
          dir = 90. - dir
          if( dir .lt. 0 ) dir = dir + 360.
          dir = dir * rad
          us(n) = sp * cos(dir)
          vs(n) = sp * sin(dir)
          goto 1
    2   continue
        close(9)

        return
   99   continue
        write(6,*) 'i = ',i
        write(6,*) line
        stop 'error stop rdspecial: line'
        end

c*****************************************************************

	subroutine get_minmax_flag(parray,nkn,pmin,pmax)

        implicit none

        integer nkn
        real parray(nkn)
        real pmin,pmax

        integer k
        real flag,high,val

        call get_flag(flag)
        high = 1.e+30

        pmin = high
        pmax = -high
        
        do k=1,nkn
          val = parray(k)
          if( val .ne. flag ) then
            pmin = min(pmin,val)
            pmax = max(pmax,val)
          end if
        end do

        end

c*****************************************************************

