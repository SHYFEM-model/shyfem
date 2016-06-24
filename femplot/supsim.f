c
c $Id: supsim.f,v 1.23 2009-11-18 17:16:00 georg Exp $
c
c utility routines for plotsim
c
c contents :
c
c subroutine plofem(type,ivarin)
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
c subroutine normuv(u,v,vv,n)			normalize horizontal velocities
c subroutine apply_dry_mask(bmask,value,n,flag)	aplies mask to nodes
c
c revision log :
c
c 30.01.2002	ggu	better debug info, allow for vertical plotting
c 30.01.2002	ggu	new routines mkvarline, aspecial, traxy, rdspecial
c 05.10.2003	ggu	changes to allow for automatical plotting of vector
c 31.10.2003	ggu	subroutine plowind()
c 17.03.2004	ggu	use okvar to decide which var to plot (0->all)
c 05.10.2004	ggu	use znv instead xv
c 05.10.2004	ggu	inorm instead ivert, overlay over z/h, comp_scale()
c 16.12.2004	ggu	plot also regular net
c 04.03.2005	ggu	gray plot over bathy
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
c 10.02.2012    ggu	belem in plobas to plot bathymetry on elements
c 26.03.2012    ccf&ggu	call mkht only for bvel .or. btrans (plo2vel)
c 13.06.2013    ggu	new routine plofem()
c 05.09.2013    ggu	endtime() and nplot introduced
c 30.05.2014    ggu	flag no data points and do not plot
c 20.10.2014    ggu	new time management
c 05.06.2015    ggu	some plotting routines adjourned (flag)
c 14.09.2015    ggu	prepared for plotting velocities given in fem file
c 06.11.2015    ggu	set valref to 1 if vel field == 0
c
c notes :
c
c customize belem in plobas() to plot bathymetry on elements
c customize bnumber in plobas() to write node and element numbers
c
c**********************************************************
c**********************************************************
c**********************************************************

	subroutine plofem(type,ivarin)

c 3D concentrations

	use mod_plot2d
	use mod_plot3d
	use levels, only : nlvdi,nlv,ilhkv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	character*(*) type
	integer ivarin

        character*80 line
	integer nrec,ivel,nplot
	integer ivaria,ivar
	integer level,isect
	logical femnext,ptime_ok,ptime_end
	integer getlev,getisec
	double precision atime

	nrec = 0
	nplot = 0
	level = getlev()
        isect = getisec()

	call femopen(type)
	ivar = ivarin
	write(6,*) 'ivar: ',ivar
	call checkvar(ivar)

        call mkvarline(ivar,line)
	ivaria = ivar

	do while( femnext(atime,ivaria,nlvdi,nkn,p3) )
	  nrec = nrec + 1
	  write(6,*) 'new record: ',nrec,ivaria,ivar,atime
	  if( ivar .ne. ivaria ) goto 99
	  !call ptime_info
	  if( ptime_end() ) exit	!FIXME
	  if( ptime_ok() ) then
	    nplot = nplot + 1
            write(6,*) '==================================='
            write(6,*) 'creating new plot with plofem... ',nplot
            write(6,*) '==================================='
            if( isect .eq. 0 ) then
	      write(6,*) '..........horizontal plotting nodes '
	      if( ivar .eq. 21 .or. ivar .eq. 2 ) then	!wind or vel
		ivel = 3		!wind
		!if( ivar .eq. 2 ) ivel = 3
		if( ivar .eq. 2 ) ivel = 5	!plot 2d vel field
		call set_uvfem(nlvdi,nkn,level,ilhkv,p3)
	        call reset_dry_mask
	        call plovel(ivel)
	      else
	        call extnlev(level,nlvdi,nkn,p3,parray)
	        call prepare_dry_mask
	        call ploval(nkn,parray,line)
	      end if
            else
	      write(6,*) '..........vertical plotting nodes '
              call plot_sect(.false.,p3)
            end if
            write(6,*) '==================================='
            write(6,*) 'end of plot'
            write(6,*) '==================================='
	  end if
	end do

	call femclose

	write(6,*) 'Total number of plots: ',nplot

	return
   99	continue
	write(6,*) 'ivar,ivaria: ',ivar,ivaria
	stop 'error stop plofem: no such variable in file'
	end

c**********************************************************

	subroutine set_uvfem(nlvddi,nkn,level,ilhkv,p3)

	use mod_hydro_plot

	implicit none

	integer nlvddi
	integer nkn
	integer level
	integer ilhkv(nkn)
	real p3(nlvddi,nkn,2)

	integer k,lev

	lev = max(1,level)

	do k=1,nkn
	  if( level .le. ilhkv(k) ) then
	    uvnode(k) = p3(lev,k,1)
	    vvnode(k) = p3(lev,k,2)
	  else
	    uvnode(k) = 0.
	    vvnode(k) = 0.
	  end if
	end do

	end

c**********************************************************

	subroutine set_uv(nlvddi,nkn,p3)

	use mod_hydro_plot

	implicit none

	integer nlvddi
	integer nkn
	real p3(nlvddi,nkn,2)

	integer k

	do k=1,nkn
	  uvnode(k) = p3(1,k,1)
	  vvnode(k) = p3(1,k,2)
	end do

	end

c**********************************************************

	subroutine plonos(type,ivar_in)

c 3D concentrations

	use mod_plot2d
	use mod_plot3d
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	character*(*) type	!default extension for file
	integer ivar_in		!desired variable id


        character*80 line
	integer nrec,it,nplot
	integer ivaria
	integer level,ivar,isect
	logical nosnext,okvar,ptime_ok,ptime_end
	integer getlev,getisec

	nrec = 0
	nplot = 0
	level = getlev()
        isect = getisec()
	ivar = ivar_in

	call nosopen(type)
	call fvlopen('.fvl')
	call checkvar(ivar)

        call mkvarline(ivar,line)

	do while( nosnext(it,ivaria,nlvdi,p3) )
	  nrec = nrec + 1
	  call fvlnext(it,nlvdi,fvlv)
	  write(6,*) 'new record: ',nrec,it,ivaria,ivar
	  !call ptime_info()
	  if( ptime_end() ) exit
	  if( ptime_ok() .and. okvar(ivaria) ) then
	    nplot = nplot + 1
            write(6,*) '==================================='
            write(6,*) 'creating new plot with plonos... ',nplot
            write(6,*) '==================================='
            if( isect .eq. 0 ) then
	      write(6,*) '..........horizontal plotting nodes '
	      call extnlev(level,nlvdi,nkn,p3,parray)
	      call prepare_dry_mask
	      call ploval(nkn,parray,line)
            else
	      write(6,*) '..........vertical plotting nodes '
              call plot_sect(.false.,p3)
            end if
            write(6,*) '==================================='
            write(6,*) 'end of plot'
            write(6,*) '==================================='
	  end if
	end do

	call nosclose

	write(6,*) 'Total number of plots: ',nplot

	end

c**********************************************************

	subroutine ploeos(type,ivar_in)

c 3D concentrations (element values)

	use mod_plot2d
	use mod_plot3d
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	character*(*) type	!default extension for file
	integer ivar_in		!desired variable id


        character*80 line
	integer nrec,it,nplot
	integer ivaria
	integer level,ivar,isect
	logical eosnext,okvar,ptime_ok,ptime_end
	integer getlev,getisec

	nrec = 0
	nplot = 0
	level = getlev()
        isect = getisec()
	ivar = ivar_in

	call eosopen(type)
	call checkvar(ivar)

        call mkvarline(ivar,line)

	do while( eosnext(it,ivaria,nlvdi,p3) )
	  nrec = nrec + 1
	  write(6,*) nrec,it,ivaria,ivar
	  if( ptime_end() ) exit
	  if( ptime_ok() .and. okvar(ivaria) ) then
	    nplot = nplot + 1
            write(6,*) '==================================='
            write(6,*) 'creating new plot with ploeos... ',nplot
            write(6,*) '==================================='
            if( isect .eq. 0 ) then
	      write(6,*) '..........horizontal plotting elements '
	      !call extelev(level,nlvdi,nkn,p3,parray)
	      call extelev(level,nlvdi,nel,p3,parray)
	      call prepare_dry_mask
	      call ploeval(nel,parray,line)
            else
	      write(6,*) '..........no vertical plotting for elements'
	      nplot = nplot - 1
            end if
            write(6,*) '==================================='
            write(6,*) 'end of plot'
            write(6,*) '==================================='
	  end if
	end do

	call eosclose

	write(6,*) 'Total number of plots: ',nplot

	end

c**********************************************************

        subroutine mkvarline(ivar,line)

c makes variable description in nos file

        implicit none

        integer ivar            !variable number
        character*(*) line      !description on return

        integer iaux
        integer ialfa

	call ivar2string(ivar,line)

	if( line .eq. ' ' ) then
          line = 'variable = '
          iaux = ialfa(float(ivar),line(12:),-1,-1)
	end if

        end

c**********************************************************

	subroutine plozet

	use mod_hydro
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer nrec,it,k,nplot
	logical ousnext,ptime_ok,ptime_end

	nrec = 0
	nplot = 0

	call ousopen

	do while( ousnext(it) )
	  nrec = nrec + 1
	  write(6,*) nrec,it
	  if( ptime_end() ) exit
	  if( ptime_ok() ) then
	    nplot = nplot + 1
	    write(6,*) '..........plotting ',nplot
	    call prepare_dry_mask
	    call ploval(nkn,znv,'water levels')
	  end if
	end do

	call ousclose

	write(6,*) 'Total number of plots: ',nplot

	end

c**********************************************************

	subroutine plobar

c plots barene

	use mod_plot2d
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer nrec,it,k,nplot
	logical ousnext,ptime_ok,ptime_end

	nrec = 0
	nplot = 0

	call ousopen

	do while( ousnext(it) )
	  nrec = nrec + 1
	  write(6,*) nrec,it
	  if( ptime_end() ) exit
	  if( ptime_ok() ) then
	    nplot = nplot + 1
	    write(6,*) '..........plotting ',nplot
	    call plozbar
	  end if
	end do

	call ousclose

	write(6,*) 'Total number of plots: ',nplot

	end

c**********************************************************

	subroutine plowave

	implicit none

	integer nrec,it,nplot
        integer ivel
	logical wavenext,ptime_ok,ptime_end

	nrec = 0
	nplot = 0
        ivel = 4	!wave

	call waveopen

	do while( wavenext(it) )
	  nrec = nrec + 1
	  write(6,*) nrec,it
	  if( ptime_end() ) exit
	  if( ptime_ok() ) then
	    nplot = nplot + 1
	    write(6,*) '..........plotting ',nplot
	    call reset_dry_mask
	    call plovel(ivel)
	  end if
	end do

	call waveclose

	write(6,*) 'Total number of plots: ',nplot

	end

c**********************************************************

	subroutine plosim(bvel)

	use mod_plot2d
	use mod_plot3d

	implicit none

	logical bvel

	integer nrec,it,nplot
        integer ivel,isect


	logical velnext,ptime_ok,ptime_end
	integer getisec

        if( bvel ) then
          ivel = 1
        else
          ivel = 2
        end if
        isect = getisec()

	nrec = 0
	nplot = 0

	call velopen

	do while( velnext(it) )
	  nrec = nrec + 1
	  write(6,*) nrec,it
	  if( ptime_end() ) exit
	  if( ptime_ok() ) then
	    nplot = nplot + 1
            write(6,*) '==================================='
            write(6,*) 'creating new plot with plosim... ',nplot
            write(6,*) '==================================='
            if( isect .eq. 0 ) then
	      write(6,*) '..........horizontal plotting '
	      call prepare_dry_mask
	      call plovel(ivel)
            else
	      write(6,*) '..........vertical plotting '
	      call prepare_vel(p3)
              call plot_sect(.true.,p3)
            end if
            write(6,*) '==================================='
            write(6,*) 'end of plot'
            write(6,*) '==================================='
	  end if
	end do

	call velclose

	write(6,*) 'Total number of plots: ',nplot

	end

c**********************************************************

	subroutine ploreg(nreg,preg,regpar,title,bintp,bgrid)

c plots node values

	use basin
	use mod_plot2d

	implicit none

	integer nreg
	real preg(nreg)
	real regpar(7)
        character*(*) title
	logical bintp		!interpolate on fem grid
	logical bgrid		!overlay regular grid

	integer nx,ny
	integer ierr
	real x0,y0,dx,dy
	real pmin,pmax,flag
	integer np
	real, allocatable :: pa(:)

	if( bintp ) then
	  np = nkn
	  allocate(pa(np))
	  call getreg(regpar,nx,ny,x0,y0,dx,dy,flag)
	  call setregextend(.false.)
	  call setregextend(.true.)
	!write(6,*) nx,ny,x0,y0,dx,dy
	!write(6,*) preg
          call intp_reg_nodes(nx,ny,x0,y0,dx,dy,flag,preg
     +                          ,pa,ierr)
	  if( ierr /= 0 ) then
	    write(6,*) 'intp_reg_nodes ierr = ',ierr
	    !stop 'error stop ploreg: interpolation'
	  end if
	else
	  np = nreg
	  allocate(pa(np))
	!write(6,*) np
	!write(6,*) preg
	  pa = preg
	end if

	call setreg_grid(regpar)

	call qstart
        call annotes(title)
	call bash(0)

	call get_minmax_flag(pa,np,pmin,pmax)
        write(6,*) 'min/max: ',np,pmin,pmax
	call count_flag(pa,np,flag)
	call colauto(pmin,pmax)

	if( bintp ) then
          call qcomm('Plotting interpolated regular grid')
          call isoline(pa,np,0.,2)
	else
          call qcomm('Plotting regular grid')
          call isoreg(pa,np,regpar,0.,2)
	end if

        call colsh
	call bash(4)	! overlays grid
	call bash(2)

	call qend

	end

c**********************************************************

	subroutine ploval(nkn,pa,title)

c plots node values

	use mod_plot2d

	implicit none

	integer nkn
	real pa(nkn)
        character*(*) title

	real pmin,pmax,flag
	real getpar

        call get_flag(flag)

	call qstart
        call annotes(title)
	call bash(0)

	call get_minmax_flag(pa,nkn,pmin,pmax)
        write(6,*) 'min/max: ',nkn,pmin,pmax
	call apply_dry_mask(bkwater,pa,nkn,flag)	!flags to dry nodes
	call count_flag(pa,nkn,flag)
	call colauto(pmin,pmax)

        call qcomm('Plotting isolines')
        call isoline(pa,nkn,0.,2)
	call plot_dry_areas
        call colsh

	call bash(4)	! overlays grid

	call bash(2)
	call qend

	end

c**********************************************************

	subroutine ploeval(nel,pa,title)

c plots element values

	use mod_plot2d

	implicit none

	integer nel
	real pa(nel)
        character*(*) title

	real pmin,pmax,flag
	real getpar

        call get_flag(flag)

	call qstart
        call annotes(title)
	call bash(0)

	call get_minmax_flag(pa,nel,pmin,pmax)
        write(6,*) 'min/max: ',nel,pmin,pmax
	call apply_dry_mask(bwater,pa,nel,flag)	!flags to dry nodes
	call colauto(pmin,pmax)

        call qcomm('Plotting element values')
        call isoline(pa,nel,0.,3)			!plot on elements
	call plot_dry_areas
        call colsh

	call bash(4)	! overlays grid

	call bash(2)
	call qend

	end

c**********************************************************

	subroutine plozbar

	use mod_hydro
	use basin

	implicit none

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

	if( ivel == 3 ) then
	  write(6,*) 'plotting wind...'
	  call plo2vel(ivel,'2D ')
	else if( ivel == 4 ) then
	  write(6,*) 'plotting wave...'
	  call plo2vel(ivel,'2D ')
	else if( ivel == 5 ) then
	  write(6,*) 'plotting velocities...'
	  call plo2vel(ivel,'3D ')
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
c ivel = 5	velocities already prepared
c
c for ivel == 1,2:	utrans,vtrans must be set
c for other values:	uvnode,vvnode must be set

	use mod_hydro_plot
	use mod_plot2d
	use mod_depth
	use mod_hydro
	use basin

	implicit none

	integer ivel			!what to plot
        character*(*) title

	real uvmod(nkn)

	logical boverl,bnode,bbound,bnorm
        logical bspecial
	logical bdebug
	logical bregular
	logical bvel,btrans,bwind,bwave,bvelok
        character*80 anoline
        character*80 spcvel
	integer ie,k
	integer i,j,nnn
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

	character*13, save :: varname(5) = (/'velocity     '
     +				,'transport    ','wind velocity'
     +				,'waves        ','velocity     '/)

	real getpar,getcol
	integer getlev

c------------------------------------------------------------------
c set up parameters
c------------------------------------------------------------------

c	-----------------------------------------------------------
c	locally defined parameters (please change here)
c	-----------------------------------------------------------

	bdebug = .true.
	bdebug = .false.
	bnode = .true.		!plot arrows on nodes
	bnode = .false.		!plot arrows on nodes
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
        bvelok = ivel .eq. 5

	if( ivel < 1 .or. ivel > 5 ) then
          stop 'error stop plo2vel: internal error (3)'
        end if

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

        anoline = title//varname(ivel)

	if( bvel .or. btrans ) call mkht(hetv,href)

c	-----------------------------------------------------------
c	start plotting
c	-----------------------------------------------------------

	call qstart
	call annotes(anoline)
	call bash(0)

c------------------------------------------------------------------
c see if regular grid
c------------------------------------------------------------------

	call getgeoflag(flag)
	call make_regular(nx,ny,bregular)

c------------------------------------------------------------------
c prepare for velocity or transport
c------------------------------------------------------------------

	valref = velref
	valmin = velmin

	if( bvel ) then
	  call por2vel(nel,utrans,vtrans,uvelem,vvelem,hetv)
        else if( btrans ) then
	  uvelem = utrans
	  vvelem = vtrans
        else if( bwind .or. bwave .or. bvelok ) then
	  !ok
	else
	  stop 'error stop plo2vel: internal error (77)'
	end if

c------------------------------------------------------------------
c compute values on nodes -> 0 for dry nodes
c------------------------------------------------------------------

        if( bwind .or. bwave .or. bvelok ) then
	  call intp2elem(uvnode,uvelem,bwater)
	  call intp2elem(vvnode,vvelem,bwater)
        else
	  call intp2node(uvelem,uvnode,bwater)
	  call intp2node(vvelem,vvnode,bwater)
        end if

	nnn = 0
        if( nnn > 0 ) then
          write(112,*) 'writing uvelem ',nel
          write(112,*) (uvelem(i),i=1,nel,nnn)
          write(112,*) 'writing vvelem ',nel
          write(112,*) (vvelem(i),i=1,nel,nnn)
          write(112,*) 'writing utrans ',nel
          write(112,*) (utrans(i),i=1,nel,nnn)
          write(112,*) 'writing vtrans ',nel
          write(112,*) (vtrans(i),i=1,nel,nnn)
          write(112,*) 'writing uvnode ',nkn
          write(112,*) (uvnode(i),i=1,nkn,nnn)
          write(112,*) 'writing vvnode ',nkn
          write(112,*) (vvnode(i),i=1,nkn,nnn)
        end if

c------------------------------------------------------------------
c regular interpolation
c------------------------------------------------------------------

	if( bregular ) then
	  call av2amk(bwater,uvnode,ureg,nx,ny)
	  call av2amk(bwater,vvnode,vreg,nx,ny)
	end if

c------------------------------------------------------------------
c compute modulus and maximum of velocity
c------------------------------------------------------------------

	call moduv(uvnode,vvnode,uvmod,nkn,uvmax,uvmed) !compute mod/max/aver

	if( .not. bbound ) then		!set boundary vectors to 0
	  call bnd2val(uvnode,0.)	
	  call bnd2val(vvnode,0.)	
	end if

c------------------------------------------------------------------
c underlying color 
c------------------------------------------------------------------

	if( boverl ) then
	  if( ioverl .eq. 1 ) then		!horizontal velocities
	    uvover = uvmod
	  else if( ioverl .eq. 2 ) then		!vertical velocities
	    uvover = wsnv
          else if( ioverl .eq. 3 ) then		!water levels
	    uvover = znv
          else if( ioverl .eq. 4 ) then		!bathymetry
	    uvover = hkv
	  else
            write(6,*) 'ioverl = ',ioverl
            stop 'error stop plo2vel: value not allowed for ioverl'
	  end if
	  !call mima(uvover,nkn,pmin,pmax)
	  call get_minmax_flag(uvover,nkn,pmin,pmax)
	  call apply_dry_mask(bkwater,uvover,nkn,flag)
	  call colauto(pmin,pmax)
	  write(6,*) 'plotting overlay color... ',pmin,pmax
          call qcomm('Plotting overlay')
          call isoline(uvover,nkn,0.,2)
	end if

	call plot_dry_areas

c------------------------------------------------------------------
c get scale (if not given, typls is computed in bastlscale)
c------------------------------------------------------------------

	typls  = getpar('typls')		!is available only now...
	if( valref .le. 0. ) then		!use maximum
	    if( bnorm ) then			!normalize vectors
	      valref = uvmax
	    else
	      valref = uvmax / fscale		!use empirical smaller value
	    end if
	end if
	typsca = typls * typlsf
	if( valref <= 0. ) valref = 1		!if vel field == 0

	write(6,*) 'plo2vel: '
	write(6,*) 'scale (type): ',typls,typlsf,ioverl
	write(6,*) 'scale (value): ',valref,valmin,uvmax,uvmed

c------------------------------------------------------------------
c plotting of arrow
c------------------------------------------------------------------

	call qgray(0.)

	if( bspecial ) then
	  scale = typsca / valref
          call aspecial(spcvel,xgv,ygv,uvnode,vvnode,scale)
        else if( bregular ) then
	  call getgeo(x0,y0,dx,dy,flag)
	  write(6,*) 'regular grid: ',dx,dy,nx,ny
	  do j=1,ny				!FIXME -> maybe nx,ny is enough
	    do i=1,nx
	      u = ureg(i,j)
	      v = vreg(i,j)
	      if( u > flag .and. v > flag ) then
		xm = x0 + (i-1) * dx
		ym = y0 + (j-1) * dy
		call comp_scale(inorm,typsca,valmin,valref,u,v,scale)
	        call pfeil(xm,ym,ureg(i,j),vreg(i,j),scale)
	      end if
	    end do
	  end do
	else if( bnode ) then
	  do k=1,nkn
	    u = uvnode(k)
	    v = vvnode(k)
	    if( u > flag .and. v > flag ) then
	      call comp_scale(inorm,typsca,valmin,valref,u,v,scale)
	      call pfeil(xgv(k),ygv(k),u,v,scale)
	    end if
	  end do
	else
	  do ie=1,nel
	    call baric(ie,xm,ym)
	    u = uvelem(ie)
	    v = vvelem(ie)
	    if( u > flag .and. v > flag ) then
	      call comp_scale(inorm,typsca,valmin,valref,u,v,scale)
	      call pfeil(xm,ym,u,v,scale)
	    end if
	  end do
	end if

c------------------------------------------------------------------
c plotting of color bar or scale arrow
c------------------------------------------------------------------

	if( boverl ) then			!plot color scale
	  call colsh
	end if
	!if( .not. bnorm ) then	                !vel arrows are still in scale
	if( inorm == 0 ) then	                !vel arrows are still in scale
	  scale = typsca / valref		!arrows are in scale
	  call velsh(ivel,scale,valref)
	end if

c------------------------------------------------------------------
c end of plot
c------------------------------------------------------------------

	call bash(4)
	call bash(2)
	call qend

c------------------------------------------------------------------
c end of routine
c------------------------------------------------------------------

	end

c**********************************************************

	subroutine plo3vel(ivel)

	use mod_hydro_plot
	use mod_hydro_vel
	use mod_hydro
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer ivel

	integer ie,ii,k,l
	integer level,lmax,lact
        real utot,vtot
	logical bdebug

	integer getlev

	bdebug = .true.
	bdebug = .false.

c make vertical velocities

	call make_vertical_velocity

c handle level

	level = getlev()
	write(6,*) 'plo3vel: level = ',level,' ivel = ',ivel
	write(6,*) 'nlvdi: ',nlvdi,'  nlv: ',nlv

	wsnv = 0.

	if( ivel .eq. 3 .or. ivel .eq. 4 ) then	!wave or wind
	  call plo2vel(ivel,'3D ')
	  return
	end if

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
            utrans(ie) = utot
            vtrans(ie) = vtot
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
	      utrans(ie) = 0.
	      vtrans(ie) = 0.
	    else
	      utrans(ie) = utlnv(lact,ie)
	      vtrans(ie) = vtlnv(lact,ie)
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
	  write(6,*) ie,level,ilhv(ie),utrans(ie),vtrans(ie)
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
	real up(n),vp(n),uv(n),vv(n),ht(n)

	where( ht > 0. )
	  uv = up / ht
	  vv = vp / ht
	else where
	  uv = 0.
	  vv = 0.
	end where

	end

c**********************************************************

	subroutine plobas

	use mod_depth
	use basin

	implicit none

	logical bnumber,belem
        real pmin,pmax

	bnumber = .true.	! plot node and element numbers
	bnumber = .false.
	belem = .true.		! plot bathymetry on elements
	belem = .false.

c only boundary line

	call qstart
	call bash(0)
	call bash(2)
	call qend

c grid

	call qstart
	call bash(0)
	call bash(1)
	call bash(2)
	call qend

c grid with gray

	call qstart
	call bash(0)
	call bash(3)
	call bash(2)
	call qend

c bathymetry (gray or color)

	call reset_dry_mask

	if( belem ) then
	  call ploeval(nkn,hev,'basin')
	else
	  call ploval(nkn,hkv,'basin')
	end if

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
	call bash(2)
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
	call bash(2)
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
	call bash(2)
	call qend

	call qstart
	call bash(0)
	call bash(3)
	call basin_number(2)
	call bash(2)
	call qend

	call qstart
	call bash(0)
	call bash(3)
	call basin_number(-1)
	call bash(2)
	call qend

	call qstart
	call bash(0)
	call bash(3)
	call basin_number(-2)
	call bash(2)
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

	use basin

	implicit none

	integer ie
	real color

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

        call plot_basin(0)
	call frame(0)

        call qcomm('Plotting elements')
	do ie=1,nel
	  call plotel(ie,color(ie))
	end do

        call qcomm('Plotting basin')
        call plot_basin(2)

	call qend

	end

c*****************************************************************

	subroutine bnd2val(a,val)	

c sets values on boundary to val

	use mod_geom
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	real a(1)
	real val

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
	real u(n), v(n), uv(n)
	real uvmax
	real uvmed

	integer i,nn
	real um,vm,rmod,rmax,rmed
	real flag

        call get_flag(flag)
	rmax = 0.
	rmed = 0.
	nn = 0

	do i=1,n
	  um = u(i)
	  vm = v(i)
	  if( um > flag .and. vm > flag ) then
	    nn = nn + 1
	    rmod = sqrt( um*um + vm*vm )
	    rmax = max(rmax,rmod)
	    rmed = rmed + rmod
	    uv(i) = rmod
	  else
	    uv(i) = flag
	  end if
	end do

	uvmax = rmax
        uvmed = 0.
	if( nn > 0 ) uvmed = rmed / nn

	end

c*****************************************************************

	subroutine normuv(u,v,vv,n)

c normalize horizontal velocities

	implicit none

	integer n
	real u(n), v(n), vv(n)

	where( vv > 0. )
	  u = u / vv
	  v = v / vv
	end where

	end

c*****************************************************************

	subroutine intp2elem(vnv,vev,bwater)

c compute elemental values vev()

	use basin

	implicit none

        real vnv(nkn)
	real vev(nel)
	logical bwater(nel)

	integer ie,ii,k,iflag
	real sum,flag

        call get_flag(flag)

	do ie=1,nel
	  if( bwater(ie) ) then
            sum = 0.
	    iflag = 0
	    do ii=1,3
	      k = nen3v(ii,ie)
	      sum = sum + vnv(k)
	      if( vnv(k) > flag ) iflag = iflag + 1
	    end do
            vev(ie) = sum / 3.
	    if( iflag .ne. 3 ) vev(ie) = flag
          else
            vev(ie) = 0.
	  end if
	end do

	end

c*****************************************************************

	subroutine intp2node(vev,vnv,bwater)

c compute nodal values vnv()

	use basin

	implicit none

	real vev(nel)
        real vnv(nkn)
	logical bwater(nel)

	integer ie,ii,k
	real r,flag
	real v2v(nkn)

        call get_flag(flag)

	do k=1,nkn
	  vnv(k) = 0.
	  v2v(k) = 0.
	end do

	do ie=1,nel
	 if( bwater(ie) ) then
	  do ii=1,3
	    if( vev(ie) > flag ) then
	      k = nen3v(ii,ie)
	      vnv(k) = vnv(k) + vev(ie)
	      v2v(k) = v2v(k) + 1.
	    end if
	  end do
	 end if
	end do

	do k=1,nkn
	  r = v2v(k)
	  if( r .gt. 0. ) then
	    vnv(k) = vnv(k) / v2v(k)
	  else
	    vnv(k) = 0.
	    vnv(k) = flag
	  end if
	end do

	end

c*****************************************************************

	subroutine count_flag(value,n,flag)

c counts total number of flagged nodes

	implicit none

	integer n
	real flag
	real value(1)

	integer i,iflag

	iflag = 0
	do i=1,n
	  if( value(i) .eq. flag ) iflag = iflag + 1
	end do

	write(6,*) 'flagged values: (',flag,')',iflag,' of ',n

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

	implicit none

	real xo,yo,u,v,dist,xn,yn

	real sp,un,vn

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

	subroutine get_minmax_flag(pa,nkn,pmin,pmax)

        implicit none

        integer nkn
        real pa(nkn)
        real pmin,pmax

        integer k
        real flag,high,val

        call get_flag(flag)
        high = 1.e+30

        pmin = high
        pmax = -high
        
        do k=1,nkn
          val = pa(k)
          if( val .ne. flag ) then
            pmin = min(pmin,val)
            pmax = max(pmax,val)
          end if
        end do

        end

c*****************************************************************

	subroutine plot_dry_areas

c plots node values

	use mod_plot2d
	use basin

	implicit none

	integer ie
	real dgray,hgray,h
	real x(3),y(3)

	double precision dgetpar

	dgray = dgetpar('dgray')
	hgray = dgetpar('hgray')
	if( dgray < 0 ) return

	call qgray(dgray)

	do ie=1,nel
	  h = sum(hm3v(:,ie))/3.
	  if( bwater(ie) .and. h > hgray ) cycle
	  call set_xy(ie,x,y)
	  call qafill(3,x,y)
	end do

	call qgray(0.)

	end

c*****************************************************************

	subroutine make_regular(nx,ny,bregular)

	use mod_hydro_plot

	implicit none

	integer nx,ny
	logical bregular

	real x0,y0,dx,dy,flag
	real xmin,ymin,xmax,ymax

	nx = 0
	ny = 0
	call getgeo(x0,y0,dx,dy,flag)
	bregular = dx .gt. 0. .and. dy .gt. 0.
	if( .not. bregular ) return
	
	call getbas(xmin,ymin,xmax,ymax)
	nx = nint((xmax-xmin)/dx)
	ny = nint((ymax-ymin)/dy)

	call hydro_plot_regular(nx,ny)	!nx,ny may be changed here

	!write(6,*) 'ggu: x0,y0,dx,dy ',x0,y0,dx,dy,nx,ny
	!write(6,*) 'ggu: ',xmin,ymin,xmax,ymax

	end

c*****************************************************************

