c
c $Id: lagrange_inout.f,v 1.2 2009-03-11 16:25:59 georg Exp $
c
c handles insertion and output of particles
c
c revision log :
c
c 05.02.2009    ggu	copied from lagrange_cont.f and integrated from others
c 16.12.2011    ggu	initialization with body id
c 23.01.2012    ggu	new call to insert_particle, new id_body
c 01.10.2012    ggu	output concentrations only for one station
c 28.03.2014    ggu	bug fix for insert with ie=0, area in concentrations
c 23.04.2014    ggu	new 3d insertion, new version of lgr file, new copy
c 06.05.2015    ccf	write relative total depth and type to output
c 07.05.2015    ccf	assign settling velocity to particles
c 07.10.2015    mic	seed 3d between surface and l_bot 
c 15.02.2016    ggu&fra	new release type ipvert=-1
c
c*******************************************************************

	subroutine insert_particle(ie,ity,lb,rtime,x,y,z)

c inserts particle at position x,y
c
c this is the ultimate routine that is called for insertion
c
c tin is insert time - compute with rtime
c z = 0.5 (corresponding to larva in water) 

	use mod_lagrange
        use basin, only : neldi,iarv
 
	implicit none

	include 'femtime.h'

	integer ie	!element of particle - if unknown use 0
	integer ity	!particle type
	integer lb	!layer [1-lmax]
	real rtime	!fraction of time step to be inserted (0 for all) [0-1]
	real x,y	!coordinates of particle to be inserted
	real z		!vertical (relative) coordinate [0-1]

	integer pt	!particle type
	real pc		!custom property
	double precision ps		!settling velocity of particle [m/s]
	logical bdebug
	logical lmytype			!use element type for particle type

c rtime is the relative time position when the particle is inserted
c 0 inserts at the beginning of the time step - advection for full time step
c 1 inserts at end of time step - no advection for this time step needed

	integer ip
	double precision xx,yy,xi(3)

	bdebug = .true.
	bdebug = .false.
	lmytype = .true.
	lmytype = .false.

	nbdy = nbdy + 1
	idbdy = idbdy + 1

	if( nbdymax > 0 .and. nbdy .gt. nbdymax ) goto 99
	call mod_lagrange_handle_alloc(nbdy)

        lgr_ar(nbdy)%tin = t_act - dt_act*(1.-rtime)
        !tin(nbdy) = it - real(idt)*(1.-rtime)
        !tin(nbdy) = it				!old

	if( ie .eq. 0 ) then
	  call find_element(x,y,ie)
	end if
	if( ie .eq. 0 ) then
	  write(6,*) 'inserting particle with ie = 0'
	end if

	!call lagr_connect_get_station(ie,ip,z)	!connectivity (z is color)

        lgr_ar(nbdy)%est = ie
        lgr_ar(nbdy)%xst = x
        lgr_ar(nbdy)%yst = y
        lgr_ar(nbdy)%zst = z
        lgr_ar(nbdy)%ie  = ie
        lgr_ar(nbdy)%l   = lb
        lgr_ar(nbdy)%x   = x
        lgr_ar(nbdy)%y   = y
        lgr_ar(nbdy)%z   = z

	xi = 1./3.
	if( ie > 0 ) then
	  xx = x
	  yy = y
	  call xy2xi(ie,xx,yy,xi)
	  call xi_move_inside(ie,xi)
	  lgr_ar(nbdy)%xi(:) = xi
	end if
	call track_xi_check('insert particle',idbdy,xi)
	if( bdebug ) then
	  write(6,*) 'debug insert particle: ',nbdy
	  write(6,*) ie,idbdy
	  write(6,*) x,y,z
	  write(6,*) ity,lb,rtime
	  write(6,*) lgr_ar(nbdy)%xst
	  write(6,*) lgr_ar(nbdy)%yst
	  write(6,*) lgr_ar(nbdy)%xi
	end if
	lgr_ar(nbdy)%id = idbdy

	call lagr_connect_bitmap_init(nbdy)

	pt = ity 	!by default use boundary number
	if( lmytype ) then
	  pt = iarv(ie) !use type of element ie
	end if 

	call lgr_set_properties(bsedim,blarvae,boilsim,pt,ps,pc)

	lgr_ar(nbdy)%ty = pt
	lgr_ar(nbdy)%sv = ps
	lgr_ar(nbdy)%c  = pc

! ..............

	return
   99	continue
	write(6,*) 'nbdymax = ',nbdymax
	write(6,*) 'Cannot insert more than nbdymax particles'
	write(6,*) 'Please change nbdymax in STR file'
	stop 'error stop insert_particle: nbdy'
	end

c*******************************************************************

	subroutine xi_move_inside(ie,xi)

c moves particle a small distance inside the element
c this is done to not fall onto the side or vertex

	integer ie
	double precision xi(3)

	integer is,it,i
	double precision dist,disttot
	double precision, parameter :: eps = 1.e-5

	is = 0
	in = 0
	do i=1,3
	  if( xi(i) == 0. ) then
	    is = is + i
	    in = in + 1
	  end if
	end do

	if( in == 0 ) return

	if( in == 1 ) then		!particle on side
	  disttot = 0.
	  do i=1,3
	    if( i /= is ) then
	      dist = eps*xi(i)
	      xi(i) =  xi(i) - dist
	      disttot = disttot + dist
	    end if
	  end do
	  xi(is) = xi(is) + disttot
	else if( in == 2 ) then		!particle on vertex
	  is = 6 - is
	  xi(is) = 1. - 2.*eps
	  do i=1,3
	    if( i /= is ) xi(i) = xi(i) + eps
	  end do
	else
	  stop 'error stop xi_move_inside: internal error'
	end if

	end

c*******************************************************************

	subroutine insert_particle_3d(ie,ity,rtime,x,y)

	use mod_lagrange
	use levels

	implicit none

	integer ::  ie	!element of particle - if unknown use 0
	integer :: ity	!type of particle to be inserted
	real :: rtime	!fraction of time step to be inserted (0 for all) [0-1]
	real :: x,y	!coordinates of particle to be inserted

c n = abs(itype)
c itype == 0	release one particle in surface layer
c itype > 0	release n particles regularly
c itype < 0	release n particles randomly
c itype == -1	release particles in every layer

	logical b2d,bdebug
	integer n,lmax,l,i
	real h,dh,hact,r
	integer itype	!type of vertical distribution
	integer lb	!layer [1-lmax]
	integer linf	!bottom layer [1-lmax]
	integer lsrf	!surface layer [1-lmax]
	real z		!vertical (relative) coordinate [0-1]
	real hl(nlv)
	real htot,htotz

	bdebug = .true.
	bdebug = .false.

	lmax = ilhv(ie)
	linf = linbot	 
	lsrf = lintop	 
	itype = ipvert
	b2d = nlv <= 1

	if( itype /= 0 ) then
	  if( linf .gt. lmax ) then 
	    linf = lmax 
	  else if ( linf .eq. 0 ) then 
	    linf = lmax 
	  end if	
	  if( lsrf .lt. 1 ) then 
	    lsrf = 1 
	  else if ( lsrf .eq. 0 ) then 
	    lsrf = 1 
	  end if	
	end if

	!-------------------------------------------------
	! handle 2d or surface situation
	!-------------------------------------------------

	if( blgrsurf .or. b2d .or. itype == 0 ) then
	  lb = 1
	  z = 0.5
	  call insert_particle(ie,ity,lb,rtime,x,y,z)
	  return
	end if

	if( itype == -1 ) then	! realease one particle in every layer
	  do l=lsrf,linf
	    z = 0.5
	    call insert_particle(ie,ity,l,rtime,x,y,z)
	  end do
	  return
	end if

	!-------------------------------------------------
	! handle 3d situation
	!-------------------------------------------------

	n = abs(itype)
	call lagr_layer_thickness(ie,lmax,hl,htot,htotz)
	h = 0.

	do l=1,linf 
	  h = h + hl(l)
	end do
	dh = h / n
	hact = -dh/2.

	if( bdebug ) write(6,*) '3d release for lagrange ',n,x,y

	do i=1,n
	  if( itype > 0 ) then
	    hact = hact + dh		!regular subdivision
	  else
	    call random_number(r)	!randomly seed
	    hact = r*h
	  end if
	  call find_vertical_position(lmax,hl,hact,lb,z)
	  call insert_particle(ie,ity,lb,rtime,x,y,z)
	  if( bdebug ) write(6,*) '   ',i,hact,lb,z
	end do
	
	!-------------------------------------------------
	! end of routine
	!-------------------------------------------------

	end

c*******************************************************************

	subroutine insert_particle_surface(ie,ity,rtime,x,y)

	implicit none

	integer ie	!element of particle - if unknown use 0
	integer ity	!type of particle to be inserted
	real rtime	!fraction of time step to be inserted (0 for all) [0-1]
	real x,y	!coordinates of particle to be inserted

	integer lb	!layer [1-lmax]
	real z		!vertical (relative) coordinate [0-1]

	lb = 1
	z = 0.5

	call insert_particle(ie,ity,lb,rtime,x,y,z)

	end

c*******************************************************************

	subroutine find_vertical_position(lmax,hl,hact,lb,z)

c finds vertical position for depth hact

	implicit none

	integer lmax		!total number of layers
	real hl(lmax)		!layer thickness
	real hact		!vertical depth of particle
	integer lb		!layer [1-lmax] (return)
	real z			!position in layer [0-1] (return, 1 bottom)

	integer l
	real hbot,htop

	hbot = 0.
	htop = 0.

	do l=1,lmax
	  hbot = hbot + hl(l)
	  if( hbot >= hact ) exit
	  htop = hbot
	end do
	if( l > lmax ) goto 99

	lb = l
	z = (hact-htop)/(hbot-htop) 	!relative position

	!write(6,*) 'lb,z,hl: ',lb,z,hl

	return
   99	continue
	write(6,*) 'lmax,hbot,hact: ',lmax,hbot,hact
	write(6,*) hl
	stop 'error stop find_vertical_position: hact > hbot'
	end

c*******************************************************************

	subroutine copy_particle(ifrom,ito)

c copies particle from ifrom to ito

	use mod_lagrange

	implicit none

	integer ifrom,ito

	if( ifrom .eq. ito ) return

        lgr_ar(ito) = lgr_ar(ifrom)

	return

        lgr_ar(ito)%est = lgr_ar(ifrom)%est
        lgr_ar(ito)%xst = lgr_ar(ifrom)%xst
        lgr_ar(ito)%yst = lgr_ar(ifrom)%yst
        lgr_ar(ito)%zst = lgr_ar(ifrom)%zst

        lgr_ar(ito)%ie  = lgr_ar(ifrom)%ie
        lgr_ar(ito)%l   = lgr_ar(ifrom)%l 
        lgr_ar(ito)%x   = lgr_ar(ifrom)%x 
        lgr_ar(ito)%y   = lgr_ar(ifrom)%y 
        lgr_ar(ito)%z   = lgr_ar(ifrom)%z 

	lgr_ar(ito)%tin = lgr_ar(ifrom)%tin
	lgr_ar(ito)%ty  = lgr_ar(ifrom)%ty
	lgr_ar(ito)%sv  = lgr_ar(ifrom)%sv
	lgr_ar(ito)%c   = lgr_ar(ifrom)%c
	lgr_ar(ito)%id  = lgr_ar(ifrom)%id
	
        lgr_ar(ito)%xi(:)  = lgr_ar(ifrom)%xi(:)

	call lagr_connect_bitmap_copy(ifrom,ito)

	end

c*******************************************************************

	subroutine delete_particle(ip)

c deletes particle ip

	use mod_lagrange

	implicit none

	integer ip

	if( lgr_ar(ip)%ie .gt. 0 ) lgr_ar(ip)%ie = 0

	end

c*******************************************************************

	subroutine ntot_particles(ntot)

c returns total number of particles

	use mod_lagrange

	implicit none

	integer ntot

	ntot = nbdy

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************
! set properties of particles
!   - settling velocity [m/s]
!   - particle type 
!   - curstom property

        subroutine lgr_set_properties(bsedim,blarvae,boilsim,pt,ps,pc)

        use lgr_sedim_module, only : lgr_set_sedim

        implicit none

	logical, intent(in)           :: bsedim	 !true for sediment lagrangian module
	logical, intent(in)           :: blarvae !true for larvae module
	logical, intent(in)           :: boilsim !true for oil module
        integer, intent(inout)        :: pt      !particle type
        double precision, intent(out) :: ps      !settling velocity [m/s]
        real, intent(out)             :: pc      !custom property

	ps = 0.
	pc = 0.
	
	if ( bsedim ) call lgr_set_sedim(pt,ps,pc)
	!if ( blarvae ) call lgr_set_larvae(pt,ps,pc) 	    !pc=length 
	!if ( boilsim ) call lgr_set_boilsim(pt,ps,pc) 	!TODO ccf

        end subroutine lgr_set_properties

c*******************************************************************

	subroutine lgr_output(iu,it)

c outputs particles to file
c
c if body has exited the element number is negativ (last element)
c once it has been written to output with negative ie, it is set to 0

c mic : 15/10/2015 write in v.5  hl = effective absolute depth 
c mic : need to add more start info: level depth custom 
c mic : suggests to write for each release a small file_ini.lgr

	use mod_lagrange
	use mod_layer_thickness
	use levels

	implicit none

	integer iu,it

	integer nvers,mtype
	parameter ( nvers = 5 , mtype = 367265 )

	integer nn,nout,ie,i,id,lb,l,ty,est
	real x,y,z,c,xst,yst,zst,tin
	real hl,ht,hr
	real getpar

	integer icall
	save icall
	data icall / 0 /

c----------------------------------------------------------------
c initialize params
c----------------------------------------------------------------

	if( icall .eq. 0 ) then
	  write(iu) mtype,nvers
	  write(iu) nlv			!from version 5 on
	  write(6,*) 'lgr file initialized...'
	end if

	icall = icall + 1

c----------------------------------------------------------------
c count active particles
c----------------------------------------------------------------

	nn = 0
	nout = 0
	do i=1,nbdy
	  ie = lgr_ar(i)%ie
	  if( ie .ne. 0 ) nn = nn + 1
	  if( ie .lt. 0 ) nout = nout + 1
	end do
	  
c----------------------------------------------------------------
c write to file
c----------------------------------------------------------------

	write(iu) it,nbdy,nn,nout

	do i=1,nbdy
          x   = lgr_ar(i)%x
          y   = lgr_ar(i)%y
          z   = lgr_ar(i)%z
	  ie  = lgr_ar(i)%ie
	  lb  = lgr_ar(i)%l
	  id  = lgr_ar(i)%id
	  ty  = lgr_ar(i)%ty
	  c   = lgr_ar(i)%c
          est = lgr_ar(i)%est
          xst = lgr_ar(i)%xst
          yst = lgr_ar(i)%yst
          zst = lgr_ar(i)%zst
	  tin = lgr_ar(i)%tin

	  call lgr_compute_depths(ie,lb,z,hl,ht,hr)

	  if( ie .ne. 0 ) then
	    if( nvers .eq. 3 ) then
	      write(iu) id,x,y,z,ie,xst,yst,zst,est
	    else if( nvers .eq. 4 ) then
	      write(iu) id,x,y,z,ie,xst,yst,zst,est,tin
	    else if( nvers .eq. 5 ) then
	      write(iu) id,x,y,z,ie,lb,hr,ty,c
     +			,xst,yst,zst,est,tin
	    else
	      write(6,*) 'nvers = ',nvers
	      stop 'error stop lgr_output: nvers unknown'
	    end if
	  end if
	  if( ie .lt. 0 ) lgr_ar(i)%ie = 0		!flag as out
	end do

c----------------------------------------------------------------
c end of routine
c----------------------------------------------------------------

	end

c*******************************************************************

	subroutine lgr_compute_depths(ie,lb,z,hl,ht,hr)

	use mod_layer_thickness
	use levels

	implicit none

	integer ie
	integer lb
	real z
	real hl,ht,hr

	integer ieh,l

	hl = 0. !absolute depth in water column
	ht = 0. !absolute depth to bottom of water column
	hr = 0. !relative depth in water column 

	ieh = abs(ie)
	if( ieh == 0 ) return

	do l = 1,lb-1
	  hl = hl + hdenv(l,ieh)
	  ht = ht + hdenv(l,ieh)
	end do

	hl = hl + z*hdenv(lb,ieh)

	do l = lb,ilhv(ieh)
	  ht = ht + hdenv(l,ieh)
	end do

	hr = hl / ht

	end

c*******************************************************************

        subroutine lgr_output_concentrations

c outputs particles as density (concentration) to NOS file

	use mod_lagrange
	use evgeom
	use basin

        implicit none

	include 'femtime.h'

	integer ie,ii,k
	integer ic,i
	integer ip,ip_station
	integer nvar,ivar,nlvdi
	real area_elem,area_node,z

        integer ecount(nel)
        integer kcount(nkn)
        real area(nkn)
        real density(nkn)

	integer ifem_open_file

        integer iunit
        save iunit
        data iunit / 0 /

c	ip_station = 5
	ip_station = 0	!if different from 0 -> output only this station

c---------------------------------------------------------
c initialize
c---------------------------------------------------------

	do ie=1,nel
	  ecount(ie) = 0
	end do
	do k=1,nkn
	  kcount(k) = 0
	  area(k) = 0.
	end do

c---------------------------------------------------------
c count particles in elements
c---------------------------------------------------------

	do i=1,nbdy
	  ie = lgr_ar(i)%ie
          z  = lgr_ar(i)%z
	  call lagr_connect_get_station(ie,ip,z)	! connectivity
	  if( ie .ne. 0 ) then
	    if( ip*ip_station .eq. 0 .or. ip .eq. ip_station ) then
	      ecount(ie) = ecount(ie) + 1
	    end if
	  end if
	end do

c---------------------------------------------------------
c copy to nodes
c---------------------------------------------------------

	do ie=1,nel
	  ic = ecount(ie)
	  area_node = 4*ev(10,ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    kcount(k) = kcount(k) + ic
	    area(k) = area(k) + area_node
	  end do
	end do

c---------------------------------------------------------
c compute density
c---------------------------------------------------------

	do k=1,nkn
	  density(k) = kcount(k) / area(k)	!in area is total area of node
	end do

c---------------------------------------------------------
c write to file
c---------------------------------------------------------

	nvar = 1
	ivar = 131
	nlvdi = 1
	call conwrite(iunit,'lgc',nvar,ivar,nlvdi,density)

c---------------------------------------------------------
c end of routine
c---------------------------------------------------------

        end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine print_elements(text)

c writes element numbers of particles to terminal

	use mod_lagrange

	implicit none

	character*(*) text

	integer i,ii

	write(6,*) text
	write(6,*) 'particles: ',nbdy
	do i=1,nbdy,10
	  write(6,*) (lgr_ar(ii)%ie,ii=i,min(i+9,nbdy))
	end do

	end

c*******************************************************************
	
	subroutine compress_particles

c deletes particles not in system and compresses array

	use mod_lagrange

	implicit none

	include 'femtime.h'

	integer ifree,icopy,i
	integer nbefore

	logical is_free,is_particle
	is_free(i) = lgr_ar(i)%ie .le. 0
	is_particle(i) = lgr_ar(i)%ie .gt. 0

	ifree = 0
	icopy = nbdy+1
	nbefore = nbdy

	do while( ifree .lt. icopy )

	  do i=ifree+1,icopy-1
	    if( is_free(i) ) goto 1
	  end do
    1	  continue
	  ifree = i
	  if( ifree .eq. icopy ) goto 3

	  do i=icopy-1,ifree+1,-1
	    if( is_particle(i) ) goto 2
	  end do
    2	  continue
	  icopy = i
	  if( ifree .eq. icopy ) goto 3

	  call copy_particle(icopy,ifree)
	  call delete_particle(icopy)

	end do
    3	continue

	nbdy = ifree - 1

c--------------------------------------------------------------
c check particle structure
c--------------------------------------------------------------

	write(lunit,*) 'compress_particles: ',it,nbefore,nbdy,nbefore-nbdy

	do i=1,nbdy
	  if( is_free(i) ) then
	    write(6,*) 'particle i = ',i,' is deleted'
	    stop 'error stop compress_particles: structure'
	  end if
	end do

	do i=nbdy+1,nbefore
	  if( is_particle(i) ) then
	    write(6,*) 'particle i = ',i,' is not deleted'
	    stop 'error stop compress_particles: structure'
	  end if
	end do

	call mod_lagrange_handle_alloc(nbdy)

	end

c*******************************************************************

	subroutine test_compress

c tests compress routine

	implicit none

	integer nadd,nloop,i,ndel
	integer iloop,ntot,ip,ie
	real ggrand

	nadd = 3
	nloop = 20
	ie = 0

	do iloop=1,nloop

	do i=1,nadd
	  ie = ie + 1
	  call insert_particle_surface(ie,0,0.,1.,1.)
	end do

	call ntot_particles(ntot)

	ndel = nint(ggrand(99)*ntot)

	do i=1,ndel
	  ip = nint(ggrand(99)*ntot)
	  call delete_particle(ip)
	end do

	call print_elements('before compress')
	call compress_particles
	call print_elements('after compress')

	end do

	end

c*******************************************************************
c	call test_compress
c	end
c*******************************************************************

