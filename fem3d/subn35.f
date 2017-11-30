c
c $Id: subn35.f,v 1.17 2009-02-04 15:26:54 georg Exp $
c
c parameter changing area routines
c
c contents :
c
c function cdf(h,z0)		computes cd from h and z0
c
c subroutine rdarea		reads area section (chezy) from STR file
c subroutine ckarea		checks values for chezy parameters
c subroutine prarea		prints chezy values to log file
c subroutine tsarea		prints test message to terminal
c subroutine inarea		initializes chezy values
c
c revision log :
c
c revised on 31.08.88 by ggu  (writes real chezy on czv)
c revised on 29.11.88 by ggu  (new chezy, iarv array)
c revised on 12.04.90 by ggu  (href)
c revised on 03.06.90 by ggu  (austausch)
c revised on 26.06.97 by ggu  (implicit none, useless parts deleted)
c 25.05.1998	ggu	documentation started
c 21.08.1998	ggu	xv eliminated
c 25.05.1999	ggu	new routine bofric
c 20.01.2000    ggu     common block /dimdim/ eliminated
c 09.08.2003    ggu     bofric now returns array with friction
c 10.08.2003    ggu     completely restructured, counter from 0 to nczdum
c 04.09.2003    ggu     bug fix for missing return in get_chezy_values
c 11.01.2005    ggu     ausv eliminated (was not used anymore)
c 02.04.2007    ggu     in check -> warning only for cz=0 and Chezy/Strickler
c 10.12.2008    ggu     re-organized, deleted sp135r(), use bottom_friction()
c 29.01.2009    ggu     ausdef eliminated (chezy(5,.) is not used anymore)
c 16.02.2011    ggu     new routines to deal with nodal area code
c 21.06.2012    ggu&aar new friction for mud module
c 28.04.2015    ggu     czdef is default for all areas not given
c 12.05.2015    ggu     rewritten with modules and allocatable
c 10.04.2017    ggu     compute cd, normalized bottom stress and bottom stress
c 09.05.2017    ggu     bug fix for computing bottom stress
c 03.11.2017    mbj     new documentation for ireib
c
c***********************************************************
c***********************************************************
c***********************************************************

!==================================================================
        module chezy
!==================================================================

        implicit none

        integer, save :: nczdum = 0
        real, save, allocatable :: czdum(:,:)

        integer, save :: nz_lines = 0
	character*80, save, allocatable :: cz_lines(:)

!==================================================================
        contains
!==================================================================

	subroutine chezy_init(n)

	integer n

	nczdum = n
	allocate(czdum(6,0:n))

	czdum = -1.

	end subroutine chezy_init

!==================================================================
        end module chezy
!==================================================================

c
c-------------------------------------------------------------------
c
c DOCS  START   P_friction
c
c DOCS  FRICTION		Bottom friction
c
c The friction term in the momentum equations can be written as
c $Ru$ and $Rv$ where $R$ is the variable friction coefficient and
c $u,v$ are the velocities in $x,y$ direction respectively.
c The form of $R$ can be specified in various ways. The value of 
c |ireib| is choosing between the formulations. In the parameter
c input file a value $\lambda$ is specified that is used in 
c the formulas below. In a 2D simulation the Strickler (2) or the Chezy (3)
c formulation is the preferred option, while for a 3D simulation is is
c recommended to use the drag coefficient (5) or the roughness length
c formulation (6).
c
c |ireib|	Type of friction used (default 0):
c		\begin{description}
c		\item[0] No friction used
c		\item[1] $R=\lambda$ is constant
c		\item[2] $\lambda$ is the Strickler coefficient.
c			 In this formulation $R$ is written as
c			 $R = \frac{g}{C^2} \frac{\vert u \vert}{H}$
c			 with $C=k_s H^{1/6}$ and $\lambda=k_s$ is
c			 the Strickler coefficient. In the above
c			 formula $g$ is the gravitational acceleration,
c			 $\vert u \vert$ the modulus of the current velocity
c			 and $H$ the total water depth.
c		\item[3] $\lambda$ is the Chezy coefficient.
c			 In this formulation $R$ is written as
c			 $R = \frac{g}{C^2} \frac{\vert u \vert}{H}$
c			 and $\lambda=C$ is the Chezy coefficient.
c		\item[4] $R=\lambda/H$ with $H$ the total water depth
c		\item[5] $\lambda$ is a constant drag coefficient and $R$ is
c			 computed as $R=\lambda\frac{\vert u \vert}{H}$
c		\item[6] $\lambda$ is the bottom roughness length and $R$ is
c			 computed through the formula
c			 $R=C\frac{\vert u \vert}{H}$ with 
c			 $C=\big(\frac{0.4}{log(\frac{\lambda+0.5H}
c			 {\lambda})}\big)^2$
c		\item[7] If $\lambda \geq 1$ it specifies the Strickler 
c			coefficient (|ireib=2|), otherwise it specifies a 
c			constant drag coefficient (|ireib=5|).
c		\end{description}
c |czdef|	The default value for the friction parameter $\lambda$.
c		Depending on the value of |ireib| the coefficient $\lambda$
c		is representing linear friction, a constant drag coefficient,
c		the Chezy or Strickler parameter, or the roughness length
c		(default 0).
c |iczv|	Normally the bottom friction coefficient 
c		(such as Strickler, Chezy, etc.)
c		is evaluated at every time step (|iczv| = 1).
c		If for some reason this behavior is not desirable,
c		|iczv| = 0 evaluates this value only before the
c		first time step, keeping it constant for the
c		rest of the simulation. Please note that this is
c		only relevant if you have given more than one bottom
c		friction value (inflow/outflow) for an area. The
c		final value of $R$ is computed at every time step
c		anyway. (default 1)
c
c The value of $\lambda$ may be specified for the whole basin through
c the value of |czdef|. For more control over the friction parameter
c it can be also specified in section |area| where the friction
c parameter depending on the type of the element may be varied. Please
c see the paragraph on section |area| for more information.
c
c DOCS  END
c
c next are experimental settings
c
c		\item[8] As ireib = 6 but using a bottom roughness length 
c			computed by Sedtrans
c		\item[9] As ireib = 6 but using a bottom roughness length 
c			computed by the fluid mud module
c
c-------------------------------------------------------------------
c
c***********************************************************

	subroutine bottom_friction

c computes bottom friction
c
c rfric is given value (in czv)
c rcd is drag coefficient ( tau = rho * rcd * u**2 )
c
c for some formulations no rcd might exists

	use mod_fluidmud
	use mod_layer_thickness
	use mod_roughness
	use mod_diff_visc_fric
	use mod_hydro
	use levels
	use basin

	implicit none

	real, parameter :: drittl = 1./3.

	include 'pkonst.h'

	integer ie,ii,k,lmax
	integer ireib
	real hzg,alpha
	real hzoff
	real uso,vso,uv,uuvv
	real rfric,rcd,rr,ss,rho0

	real getpar,cdf

c-------------------------------------------------------------------
c get variables
c-------------------------------------------------------------------

	hzoff = getpar('hzoff')
	ireib = nint(getpar('ireib'))
	rho0 = rowass

c-------------------------------------------------------------------
c loop on elements
c-------------------------------------------------------------------

	do ie=1,nel

c         ----------------------------------------------------------
c	  get transport in layer
c         ----------------------------------------------------------

	  lmax = ilhv(ie)

          uso = utlov(lmax,ie)
          vso = vtlov(lmax,ie)
	  uv = sqrt(uso*uso+vso*vso)

c         ----------------------------------------------------------
c	  set total depth
c         ----------------------------------------------------------

	  hzg = hdeov(lmax,ie)
          if( hzg .lt. hzoff ) hzg = hzoff

c         ----------------------------------------------------------
c	  get friction parameter
c         ----------------------------------------------------------

	  rfric = czv(ie)

c         ----------------------------------------------------------
c	  compute friction
c         ----------------------------------------------------------

	  if(ireib.eq.0) then
		rr = 0.
		rcd = 0.
          else if(ireib.eq.1) then
                rr = rfric
		rcd = 0.
		if( uv > 0. ) rcd = rr*hzg*hzg/uv
	  else if(ireib.eq.2) then		! Strickler
		rcd = grav/((rfric**2)*(hzg**drittl))
		rr = rcd*uv/(hzg*hzg)
	  else if(ireib.eq.3) then		! Chezy
		rcd = grav/(rfric**2)
		rr = rcd*uv/(hzg*hzg)
          else if(ireib.eq.4) then
                rr = rfric/hzg
		rcd = 0.
		if( uv > 0. ) rcd = rr*hzg*hzg/uv
	  else if(ireib.eq.5) then		! constant drag coefficient
		rcd = rfric
		rr = rcd*uv/(hzg*hzg)
	  else if(ireib.eq.6) then		! rfric is z0
                rcd = cdf(hzg,rfric)
		rr = rcd*uv/(hzg*hzg)
          else if(ireib.eq.7) then		! mixed Strickler / drag
                if( rfric .ge. 1. ) then
		  rcd = grav/((rfric**2)*(hzg**drittl))
                else
		  rcd = rfric
                end if
		rr = rcd*uv/(hzg*hzg)
          else if(ireib.eq.8) then		! use z0 computed by sedtrans
                ss = 0.
                do ii=1,3
                  k = nen3v(ii,ie)
                  ss = ss + z0bk(k)
                end do
                ss = ss / 3.
                rcd = cdf(hzg,ss)
		rr = rcd*uv/(hzg*hzg)
          else if(ireib.eq.9) then		! function of fluid mud (AR:)
                ss = 0.
                do ii=1,3
                  k = nen3v(ii,ie)
                  lmax = ilhkv(k)
                  call set_mud_roughness(k,lmax,alpha) ! (ARON)
                  ss = ss + alpha * rfric ! rfric = ks for this parameterization
                end do
                ss = ss / 3.
                z0bk(k) = ss
                !z0bk(k) = max(z0bkmud(k),ss)
                !ss = rfric	!ARON: do you really need to compute ss above?
                rcd = cdf(hzg,ss)
                rr = rcd*uv/(hzg*hzg)
		!Well not really there are mainls two issues ...
		!1. Rougnes get reduced by mud this is taken into 
		!account by calling the routine above
		!2. We need to apply mixing length for the 1st grid-cell 
		!otherwise turbulence in gotm fully collapse since k-eps 
		!is only valid for isotropic turbulence. 
	  else
		write(6,*) 'unknown friction : ',ireib
		stop 'error stop bottom_friction'
	  end if

	  rfricv(ie) = rr
	  rcdv(ie) = rcd
	  uuvv = uv / hzg
	  bnstressv(ie) = rcd * uuvv * uuvv

	end do

c-------------------------------------------------------------------
c end of routine
c-------------------------------------------------------------------

	end

c***********************************************************

	subroutine current_bottom_stress(bstressv)

c computes bottom stress at nodes

	use basin
	use levels
	use evgeom
	use mod_diff_visc_fric
	use mod_ts

	implicit none

	include 'pkonst.h'

	real bstressv(nkn)

	integer ie,k,ii,lmax
	real area,rho,bnstress
	real aux(nkn)

	call bottom_friction	!bottom stress on elements (bnstressv)

	bstressv = 0.
	aux = 0.

	do ie=1,nel
	  lmax = ilhv(ie)
	  area = ev(10,ie)
	  bnstress = bnstressv(ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
            rho = rowass + rhov(lmax,k)
	    bstressv(k) = bstressv(k) + rho * bnstress * area
	    aux(k) = aux(k) + area
	  end do
	end do

	where( aux > 0. ) bstressv = bstressv / aux

	end

c***********************************************************

        function cdf(h,z0)

c computes cd from h and z0

        implicit none

        real cdf
        real h,z0

        real kappa,cds

        kappa = 0.4

        cds = kappa / log( (z0+0.5*h) / z0 )

        cdf = cds*cds

        end

c***********************************************************

	subroutine init_nodal_area_code

c interpolates area codes from elements to nodes (min or max)

	use basin
	use shympi

	implicit none

	integer init,mode
	integer k,ie,ii,ia

	mode = -1		! -1: use minimum   +1: use maximum

	init = 99999999
	if( mode .gt. 0 ) init = -init

	do k=1,nkn
	  iarnv(k) = init
	end do

	do ie=1,nel
	  ia = iarv(ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    if( mode .eq. -1 ) then
		iarnv(k) = min(iarnv(k),ia)
	    else
		iarnv(k) = max(iarnv(k),ia)
	    end if
	  end do
	end do

	!call shympi_comment('exchanging iarnv')
	call shympi_exchange_2d_node(iarnv)

        if( mode .eq. -1 ) then
          call shympi_exchange_2d_nodes_min(iarnv)
          !call shympi_comment('shympi_elem: exchange iarnv_min')
        else
          call shympi_exchange_2d_nodes_max(iarnv)
          !call shympi_comment('shympi_elem: exchange iarnv_max')
        end if

!       shympi_elem:   exchange min or max iarnv

	end

c***********************************************************

	subroutine get_nodal_area_code(k,ncode)

	use basin

	implicit none

	integer k	!node number
	integer ncode	!nodal area code (return)

	ncode = iarnv(k)

	end

c***********************************************************
c***********************************************************
c***********************************************************
c***********************************************************
c***********************************************************

	subroutine n_chezy_values(nareas)

	use chezy

	implicit none

	integer nareas

	nareas = nczdum

	end

c***********************************************************

	subroutine get_chezy_values(iar,valin,valout)

	use chezy

	implicit none

	integer iar
	real valin,valout

	if( iar .gt. nczdum ) goto 99

	valin = czdum(1,iar)
	valout = czdum(2,iar)

	return
   99	continue
	write(6,*) 'iar,nczdum: ',iar,nczdum
	stop 'error stop get_chezy_values'
	end

c***********************************************************

	subroutine set_chezy_values(iar,valin,valout)

	use chezy

	implicit none

	integer iar
	real valin,valout

	if( iar .gt. nczdum ) goto 99

	czdum(1,iar) = valin
	czdum(2,iar) = valout

	return
   99	continue
	write(6,*) 'iar,nczdum: ',iar,nczdum
	stop 'error stop set_chezy_values'
	end

c***********************************************************
c***********************************************************
c***********************************************************
c***********************************************************
c***********************************************************

	subroutine set_chezy

c initializes chezy arrays

	use mod_diff_visc_fric
	use basin
	use chezy

	implicit none

	integer ie,iar

	do ie=1,nel
	    iar=iarv(ie)
	    czv(ie)=czdum(6,iar)
	end do

	end

c***********************************************************

	subroutine init_chezy

c initializes chezy arrays

	use chezy

	implicit none

	logical bdebug
	integer i

	bdebug = .true.
	bdebug = .false.

	do i=0,nczdum
	  czdum(6,i)=czdum(1,i)
	end do

	if( bdebug ) call print_chezy

	call set_chezy

	end

c***********************************************************

	subroutine adjust_chezy

c adjusts chezy arrays

	use mod_hydro_print
	use basin
	use chezy

	implicit none

	logical bdebug
	integer i,k1,k2
	integer iczv
	real dx,dy,scal

	real getpar

	bdebug = .true.
	bdebug = .false.

	iczv=nint(getpar('iczv'))
	if( iczv .eq. 0 ) return	!chezy is not adjusted

	do i=0,nczdum
	    if(czdum(2,i).eq.0.) then
		czdum(6,i)=czdum(1,i)
	    else
		k1=nint(czdum(3,i))
		k2=nint(czdum(4,i))
		dx=xgv(k2)-xgv(k1)
		dy=ygv(k2)-ygv(k1)
		scal=dx*up0v(k1)+dy*vp0v(k1)
		if(scal.ge.0.) then
			czdum(6,i)=czdum(1,i)
		else
			czdum(6,i)=czdum(2,i)
		end if
	    end if
	end do

	if( bdebug ) call print_chezy

	call set_chezy

	end

c***********************************************************

	subroutine print_chezy

c prints chezy arrays

	use chezy

	implicit none

	integer i
	integer iunit

	iunit = 6

	write(iunit,*) 'Values for chezy (czv) :'
	do i=0,nczdum
	  write(iunit,*) i,czdum(6,i)
	end do

	end

c***********************************************************

	subroutine check_chezy

c checks chezy arrays

	use basin
	use chezy

	implicit none

	integer ie,iar
	integer i,j,k
	real cz

	do i=0,nczdum
	  cz = czdum(1,i)
	  if( cz .lt. 0. .or. cz .gt. 1.e+10 ) goto 99
	  cz = czdum(2,i)
	  if( cz .lt. 0. .or. cz .gt. 1.e+10 ) goto 99
	  if( cz .gt. 0. ) then
	    k = nint(czdum(3,i))
	    if( k .lt. 0. .or. k .gt. nkn ) goto 99
	    k = nint(czdum(4,i))
	    if( k .lt. 0. .or. k .gt. nkn ) goto 99
	  end if
	  cz = czdum(6,i)
	  if( cz .lt. 0. .or. cz .gt. 1.e+10 ) goto 99
	end do

	do ie=1,nel
	    iar=iarv(ie)
	    if( iar .gt. nczdum ) goto 98
	    cz=czdum(6,iar)
	    if( cz .lt. 0. .or. cz .gt. 1.e+10 ) goto 98
	end do

	return
   98	continue
	write(6,*) 'ie,iar,nczdum,cz: ',ie,iar,nczdum,cz
	write(6,*) (czdum(j,iar),j=1,6)
	stop 'error stop check_chezy: error in values (1)'
   99	continue
	write(6,*) 'i,iar,nczdum: ',i,i-1,nczdum
	write(6,*) (czdum(j,i),j=1,6)
	stop 'error stop check_chezy: error in values (2)'
	end

c***********************************************************
c***********************************************************
c***********************************************************
c***********************************************************
c***********************************************************

	subroutine rdarea

c reads area section (chezy) from STR file

	use nls
	use chezy

	implicit none

	integer n,i

	n = nls_read_table()

	if( n > 0 ) then
	  nz_lines = n
	  allocate(cz_lines(n))
	  call nls_copy_char_vect(n,cz_lines)
	end if

	end

c***********************************************************

	subroutine parse_area

c parses area section (chezy) from STR file

	use chezy

	implicit none

	character*80 line
	integer ianz,iar,i,n,il
	real f(10)

	integer iscanf

	do il=1,nz_lines
	  line = cz_lines(il)
	  !write(6,*) il,trim(line)
	  ianz = iscanf(line,f,10)
	  if( ianz .gt. 0 ) then
	    iar = nint(f(1))
            if(iar.lt.0) goto 88
	    if( ianz .gt. 7 ) goto 86
            if(iar.gt.nczdum) then
	      write(6,*) 'warning: no such area code... ignoring ',iar
	      cycle
	    end if
	    do i=2,ianz
	      czdum(i-1,iar) = f(i)
	    end do
	  else if( ianz .lt. 0 ) then
			goto 98
	  end if
	end do

	if( nz_lines > 0 ) then
	  nz_lines = 0
	  deallocate(cz_lines)	!we dont need it anymore
	end if

	return
   86   continue
        write(6,*) 'Too much data on line'
        write(6,*) line
        stop 'error stop : rdarea'
   88   continue
        write(6,*) 'Negative area code = ',iar,' not allowed'
        write(6,*) line
        stop 'error stop : rdarea'
   98   continue
        write(6,*) 'Read error in line :'
	write(6,*) line
        stop 'error stop : rdarea'
	end

c***********************************************************

	subroutine ckarea

c checks values for chezy parameters

	use basin
	use chezy

	implicit none

	integer i,knode,knodeh,ireib,nczmax
	logical bstop,bpos
	real czdef

	integer ipint
	real getpar

	bstop = .false.

c get default values

        ireib=nint(getpar('ireib'))
	!if( ireib .le. 0 ) return

	bpos = ireib .gt. 1 .and. ireib .ne. 5	!must be strictly positive

        czdef=getpar('czdef')

c compute maximum value of area code

	nczmax = 0
	do i=1,nel
          if(iarv(i).gt.nczmax) nczmax=iarv(i)
        end do

c allocate and parse arrays

	call chezy_init(nczmax)
	call parse_area

c check read in values

        do i=0,nczdum

         if(czdum(1,i).eq.-1.) czdum(1,i)=czdef

         if(czdum(1,i).lt.0.) then
                write(6,*) 'Friction value cannot be negative:'
                write(6,*) 'area = ',i,'  chezy = ',czdum(1,i)
                bstop=.true.
	 end if

	 if( bpos .and. czdum(1,i).eq.0.) then
                write(6,*) 'Friction value must be positive:'
                write(6,*) 'area = ',i,'  chezy = ',czdum(1,i)
                bstop=.true.
	 end if

         if(czdum(2,i).eq.-1.) czdum(2,i)=0.

         if(czdum(3,i).eq.-1. .or. czdum(3,i).eq.0.) then
           czdum(3,i)=0.
         else
           knodeh=nint(czdum(3,i))
           knode=ipint(knodeh)          !$$EXTINW
           if(knode.le.0) then
                write(6,*) 'section AREA : node not found ',knodeh
                bstop=.true.
           end if
           czdum(3,i)=knode
         end if

         if(czdum(4,i).eq.-1. .or. czdum(4,i).eq.0.) then
           czdum(4,i)=0.
         else
           knodeh=nint(czdum(4,i))
           knode=ipint(knodeh)          !$$EXTINW
           if(knode.le.0) then
                write(6,*) 'section AREA : node not found ',knodeh
                bstop=.true.
           end if
           czdum(4,i)=knode
         end if

         if(czdum(5,i).eq.-1.) czdum(5,i)=0.
         czdum(6,i)=0.

        end do

	if( bstop ) stop 'error stop ckarea'

	end

c***********************************************************

	subroutine prarea

c prints chezy values to log file

	use chezy

	implicit none

	integer ianf,i
	integer ipext

        ianf=0
        if(czdum(1,0).eq.0) ianf=1
        write(6,*)
        write(6,1007)

        do i=ianf,nczdum
            if(czdum(2,i).ne.0.) then			!with two chezy
                write(6,1008) i,czdum(1,i),czdum(2,i)
     +                          ,ipext(nint(czdum(3,i)))
     +                          ,ipext(nint(czdum(4,i)))
            else					!just one chezy
                write(6,1008) i,czdum(1,i)
            end if
        end do

	return
 1007   format(' area,cz1,cz2,k1,k2 : ')
 1008   format(i5,2e12.4,2i7,e12.4)
	end

c***********************************************************

	subroutine tsarea

c prints test message to terminal

	use chezy

	implicit none

	integer j,i

        write(6,*) '/chezy/'
        write(6,*) nczdum
        do j=0,nczdum
            write(6,'(1x,6e12.4)') (czdum(i,j),i=1,6)
        end do

	end

c***********************************************************

	subroutine inarea

c initializes chezy values

	use chezy

	implicit none

	end

c***********************************************************

