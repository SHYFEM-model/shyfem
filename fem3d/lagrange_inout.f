
!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2018  Georg Umgiesser
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

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
c 07.10.2015    mcg	seed 3d between surface and l_bot 
c 15.02.2016    ggu&fdp	new release type ipvert=-1
c 27.11.2017    ggu	inserted code to start from lgr files
c 20.07.2018    ccf	new output file format
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
	integer nc	!number of custom properties
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

        lgr_ar(nbdy)%init%time = t_act - dt_act*(1.-rtime)
        lgr_ar(nbdy)%actual%time = lgr_ar(nbdy)%init%time

	if( ie .eq. 0 ) then
	  call find_element(x,y,ie)
	end if
	if( ie .eq. 0 ) then
	  write(6,*) 'inserting particle with ie = 0'
	end if

	!call lagr_connect_get_station(ie,ip,z)	!connectivity (z is color)

        lgr_ar(nbdy)%init%ie   = -ie
        lgr_ar(nbdy)%init%l    = lb
        lgr_ar(nbdy)%init%z    = z
        lgr_ar(nbdy)%actual%ie = ie
        lgr_ar(nbdy)%actual%l  = lb
        lgr_ar(nbdy)%actual%z  = z

	xi = 1./3.
	if( ie > 0 ) then
	  call xi_init_particle(ie,x,y,xi)
	  lgr_ar(nbdy)%init%xi(:) = xi
	  lgr_ar(nbdy)%actual%xi(:) = xi
	end if
	call track_xi_check('insert particle',idbdy,xi)
	if( bdebug ) then
	  write(6,*) 'debug insert particle: ',nbdy
	  write(6,*) ie,idbdy
	  write(6,*) x,y,z
	  write(6,*) ity,lb,rtime
	  write(6,*) lgr_ar(nbdy)%init%xi
	end if
	lgr_ar(nbdy)%id = idbdy

	call lagr_connect_bitmap_init(nbdy)

	pt = ity 	!by default use boundary number
	if( lmytype ) then
	  pt = iarv(ie) !use type of element ie
	end if 

	call lgr_set_properties(bsedim,blarvae,boilsim,pt,ps,pc,nc)
        allocate( lgr_ar(nbdy)%custom(nc) )
	ncust = nc

	lgr_ar(nbdy)%type = pt
	lgr_ar(nbdy)%sinking = ps
	lgr_ar(nbdy)%custom(1)  = pc

! ..............

	return
   99	continue
	write(6,*) 'nbdymax = ',nbdymax
	write(6,*) 'Cannot insert more than nbdymax particles'
	write(6,*) 'Please change nbdymax in STR file'
	stop 'error stop insert_particle: nbdy'
	end

c*******************************************************************

	subroutine xi_init_particle(ie,x,y,xi)

	implicit none

	integer ie
	real x,y
	double precision xi(3)

	double precision xx,yy

	xx = x
	yy = y
	call xy2xi(ie,xx,yy,xi)
	call xi_move_inside(ie,xi)

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

        lgr_ar(ito)%init%ie      = lgr_ar(ifrom)%init%ie
        lgr_ar(ito)%init%z       = lgr_ar(ifrom)%init%z
        lgr_ar(ito)%init%l       = lgr_ar(ifrom)%init%l 
        lgr_ar(ito)%init%xi(:)   = lgr_ar(ifrom)%init%xi(:)
        lgr_ar(ito)%init%time    = lgr_ar(ifrom)%init%time

        lgr_ar(ito)%actual%ie    = lgr_ar(ifrom)%actual%ie
        lgr_ar(ito)%actual%z     = lgr_ar(ifrom)%actual%z
        lgr_ar(ito)%actual%l     = lgr_ar(ifrom)%actual%l 
        lgr_ar(ito)%actual%xi(:) = lgr_ar(ifrom)%actual%xi(:)
        lgr_ar(ito)%actual%time  = lgr_ar(ifrom)%actual%time

	lgr_ar(ito)%id           = lgr_ar(ifrom)%id
	lgr_ar(ito)%type         = lgr_ar(ifrom)%type
	lgr_ar(ito)%sinking      = lgr_ar(ifrom)%sinking
	lgr_ar(ito)%custom(:)    = lgr_ar(ifrom)%custom(:)

	call lagr_connect_bitmap_copy(ifrom,ito)

	end

c*******************************************************************

	subroutine delete_particle(ip)

c deletes particle ip

	use mod_lagrange

	implicit none

	integer ip

	if( lgr_ar(ip)%actual%ie .gt. 0 ) lgr_ar(ip)%actual%ie = 0

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

	subroutine lgr_insert_particle(i,id,z,ie,lb,ty,c,tin)

	use mod_lagrange

	implicit none

	integer ie,i,id,lb,ty
	real z,c,tin

        lgr_ar(nbdy)%init%ie     = -ie
        lgr_ar(nbdy)%init%l      = lb
        lgr_ar(nbdy)%init%z      = z
        lgr_ar(nbdy)%init%time   = tin

        lgr_ar(nbdy)%actual%ie   = ie
        lgr_ar(nbdy)%actual%l    = lb
        lgr_ar(nbdy)%actual%z    = z
        lgr_ar(nbdy)%actual%time = tin

        lgr_ar(nbdy)%type = ty
        lgr_ar(nbdy)%sinking = 0.
        lgr_ar(nbdy)%custom(1) = c

	end

c*******************************************************************
c TO BE DONE
	subroutine lgr_input_shell

	use mod_lagrange
        use shyfile

	implicit none

        !active particles
        integer, allocatable            :: ida(:)
        integer, allocatable            :: tya(:)
        double precision, allocatable   :: tta(:)
        real, allocatable               :: sa(:)
        integer, allocatable            :: iea(:)
        real, allocatable               :: xa(:),ya(:),za(:)
        integer, allocatable            :: lba(:)
        real, allocatable               :: hla(:)
        real, allocatable               :: ca(:,:)

	integer				:: i,ic
	integer				:: n_act,ntot,n_aux
	integer				:: id
	integer				:: iu
        integer                         :: iwhat
        integer                         :: ierr
	double precision 		:: xi(3)
	double precision 		:: dtime,ditime
	double precision 		:: dtlgin
        character*20                    :: aline
	character*80 			:: ifile

        INTERFACE
        subroutine lgr_alloc(nn,nc
     +          ,id,ty,tt,s,ie,x,y,z,lb,hl,c)
        integer nn,nc
        integer, allocatable            :: id(:)
        integer, allocatable            :: ty(:)
        double precision, allocatable   :: tt(:)
        real, allocatable               :: s(:)
        integer, allocatable            :: ie(:)
        real, allocatable               :: x(:),y(:),z(:)
        integer, allocatable            :: lb(:)
        real, allocatable               :: hl(:)
        real, allocatable               :: c(:,:)
        end subroutine
        END INTERFACE

        !--------------------------------------------------------------
        ! open input file
        !--------------------------------------------------------------

	call getfnm('lgrini',ifile)
	if( ifile == ' ' ) return
        call convert_date_d('itlgin',dtlgin)   !start of lagrangian sim
        call get_timeline(dtlgin,aline)

        id = shy_init(ifile)
        call shy_read_header(id,ierr)
        call shy_get_iunit(id,iu)
        read(iu) ncust

        !--------------------------------------------------------------
        ! loop over records
        !--------------------------------------------------------------
        do
          !----------------------------------------------------------------
          ! read lgr data block --> active particles (use only these particles)
          !----------------------------------------------------------------
          call lgr_peek_block_header(iu,ditime,n_act,iwhat,ierr)
          if( ierr > 0 ) write(6,*) 'error in reading file : ',ierr
          if( ierr /= 0 ) exit
          call lgr_alloc(n_act,ncust
     +          ,ida,tya,tta,sa,iea,xa,ya,za,lba,hla,ca)
          call lgr_get_block(iu,n_act,ncust,
     +                  ida,tya,tta,sa,iea,xa,ya,za,lba,hla,ca)
          if( ditime == dtlgin ) exit
          !----------------------------------------------------------------
          ! skip data block --> inserted particles
          !----------------------------------------------------------------
          call lgr_peek_block_header(iu,ditime,n_aux,iwhat,ierr)
          if( ierr /= 0 ) exit
          call lgr_skip_block(iu,n_aux,ncust)
          !----------------------------------------------------------------
          ! skip data block --> exited particles
          !----------------------------------------------------------------
          call lgr_peek_block_header(iu,ditime,n_aux,iwhat,ierr)
          if( ierr /= 0 ) exit
          call lgr_skip_block(iu,n_aux,ncust)
        end do

        if( ditime /= dtlgin ) goto 95

        call get_act_dtime(dtime)

	write(6,*) 'initializing lagrangian from file: '
	write(6,*) 'file = ',trim(ifile)
	write(6,*) 'time = ',aline

	ntot = nbdy + n_act
	if( nbdymax > 0 .and. ntot > nbdymax ) goto 94
	call mod_lagrange_handle_alloc(ntot)

	do i=1,n_act
          nbdy = nbdy + 1
          lgr_ar(nbdy)%init%ie     = -iea(i)
          lgr_ar(nbdy)%init%l      = lba(i)
          lgr_ar(nbdy)%init%z      = za(i)
          lgr_ar(nbdy)%init%time   = dtime
          lgr_ar(nbdy)%actual%ie   = iea(i)
          lgr_ar(nbdy)%actual%l    = lba(i)
          lgr_ar(nbdy)%actual%z    = za(i)
          lgr_ar(nbdy)%actual%time = dtime
          lgr_ar(nbdy)%id          = ida(i)
          lgr_ar(nbdy)%type        = tya(i)
          lgr_ar(nbdy)%sinking     = sa(i)

          allocate( lgr_ar(nbdy)%custom(ncust) )
          do ic = 1,ncust
            lgr_ar(nbdy)%custom(ic)= ca(i,ic)
          end do

          call xi_init_particle(iea(i),xa(i),ya(i),xi)
          lgr_ar(nbdy)%init%xi(:) = xi(:)
          lgr_ar(nbdy)%actual%xi(:) = xi(:)
	  call track_xi_check('insert particle',id,xi)

	  call lagr_connect_bitmap_init(nbdy)

	end do

        call shy_close(id)

	write(6,*) 'finished initializing lagrangian from file'

	return
   94	continue
	write(6,*) 'nbdymax = ',nbdymax
	write(6,*) 'Cannot insert more than nbdymax particles'
	write(6,*) 'Please change nbdymax in STR file'
	stop 'error stop lgr_input_shell: nbdy'
   95	continue
	write(6,*) 'cannot find time record: ',dtlgin
	stop 'error stop lgr_input_shell: no such time record'
   96	continue
	write(6,*) 'error reading file: ',iu
	stop 'error stop lgr_input_shell: read error'
   99	continue
	write(6,*) 'cannot open lgr file: ',trim(ifile)
	stop 'error stop lgr_input_shell: open file'
	end

c*******************************************************************

        subroutine lgr_write_header(iu,nlv,ncust)

        implicit none

	integer, intent(in) 		:: iu
	integer, intent(in) 		:: nlv
	integer, intent(in) 		:: ncust

        integer, parameter      	:: nvers = 6
        integer, parameter      	:: mtype = 367265


        write(iu) mtype,nvers
        write(iu) nlv,ncust                   !from version 5 on

        write(6,*) 'lgr file initialized...'

	end subroutine lgr_write_header

c*******************************************************************

	subroutine lgr_output(iu,dtime)

c outputs particles to file
c
c if body has exited the element number is negativ (last element)
c once it has been written to output with negative ie, it is set to 0

c mcg : 15/10/2015 write in v.5  hl = effective absolute depth 
c mcg : need to add more start info: level depth custom 
c mcg : suggests to write for each release a small file_ini.lgr
c ccf : new output format

	use mod_lagrange
	use mod_layer_thickness
	use levels

	implicit none

	integer, intent(in) 		:: iu
	double precision, intent(in)	:: dtime

	integer 		:: i,ip,id,lb,l,ty,ii,ie
	integer 		:: iwhat
	integer			:: n_new  	!particles inserted
	integer			:: n_act   	!particles inside the domain
	integer			:: n_ext     	!particles exited the domain
	real			:: x,y,z,s
	real 			:: hl,ht,hr
	real 			:: getpar
	double precision 	:: xx,yy
        double precision 	:: xi(3)
        double precision 	:: tt
	real, allocatable	:: c(:)
	integer, allocatable	:: list_new(:)
	integer, allocatable	:: list_ins(:)
	integer, allocatable	:: list_ext(:)
        double precision dgetpar
        integer, save 		:: icall = 0

c----------------------------------------------------------------
c count inserted, active and exited particles
c----------------------------------------------------------------

        if (nbdy == 0 ) return

	allocate(c(ncust))
	allocate(list_new(nbdy))
	allocate(list_ins(nbdy))
	allocate(list_ext(nbdy))

	n_new = 0
	n_act = 0
	n_ext = 0
	do i=1,nbdy
	  ie = lgr_ar(i)%init%ie
	  if( ie < 0 ) then
	    n_new = n_new + 1
	    list_new(n_new) = i
	  end if
	  ie = lgr_ar(i)%actual%ie
	  if( ie > 0 ) then
	    n_act = n_act + 1
	    list_ins(n_act) = i
	  end if
	  if( ie <= 0 ) then
	    n_ext = n_ext + 1
	    list_ext(n_ext) = i
	  end if
	end do

c----------------------------------------------------------------
c write to file - active particles, iwhat = 0
c----------------------------------------------------------------
	
	iwhat = 0
	write(iu) dtime,n_act,iwhat

	do i=1,n_act
	  ip  = list_ins(i)
	  ie  = lgr_ar(ip)%actual%ie
          do ii=1,3
            xi(ii) = lgr_ar(ip)%actual%xi(ii)
          end do
	  call xi2xy(abs(ie),xx,yy,xi)
          x   = xx
          y   = yy
          z   = lgr_ar(ip)%actual%z
	  lb  = lgr_ar(ip)%actual%l
	  tt  = lgr_ar(ip)%actual%time
	  id  = lgr_ar(ip)%id
	  ty  = lgr_ar(ip)%type
	  s   = lgr_ar(ip)%sinking
	  c   = lgr_ar(ip)%custom

	  call lgr_compute_depths(ie,lb,z,hl,ht,hr)

	  write(iu) id,ty,tt,s,ie
	  write(iu) x,y,z,lb,hl
	  write(iu) (c(ii), ii=1,ncust)
	end do

c----------------------------------------------------------------
c write to file - inserted particles, iwhat = -1
c----------------------------------------------------------------
	
	iwhat = -1
	write(iu) dtime,n_new,iwhat

	do i=1,n_new
	  ip  = list_new(i)
	  ie  = lgr_ar(ip)%init%ie
          do ii=1,3
            xi(ii) = lgr_ar(ip)%init%xi(ii)
          end do
	  call xi2xy(abs(ie),xx,yy,xi)
          x   = xx
          y   = yy
          z   = lgr_ar(ip)%init%z
	  lb  = lgr_ar(ip)%init%l
	  tt  = lgr_ar(ip)%init%time
	  id  = lgr_ar(ip)%id
	  ty  = lgr_ar(ip)%type
	  s   = lgr_ar(ip)%sinking
	  c   = lgr_ar(ip)%custom(1)

	  call lgr_compute_depths(ie,lb,z,hl,ht,hr)

	  write(iu) id,ty,tt,s,abs(ie)
	  write(iu) x,y,z,lb,hl
	  write(iu) (c(ii), ii=1,ncust)

	  lgr_ar(ip)%init%ie = -ie
	   !should we deleted the init information after writing it?
	end do

c----------------------------------------------------------------
c write to file - exited particles, iwhat = 1
c particles are then deleted by routine compress_particles 
c bcompress must be set to true
c----------------------------------------------------------------
	
	iwhat = 1
	write(iu) dtime,n_ext,iwhat

	do i=1,n_ext
	  ip  = list_ext(i)
	  ie  = lgr_ar(ip)%actual%ie
          do ii=1,3
            xi(ii) = lgr_ar(ip)%actual%xi(ii)
          end do
	  call xi2xy(abs(ie),xx,yy,xi)
          x   = xx
          y   = yy
          z   = lgr_ar(ip)%actual%z
	  lb  = lgr_ar(ip)%actual%l
          tt  = lgr_ar(ip)%actual%time
	  id  = lgr_ar(ip)%id
	  ty  = lgr_ar(ip)%type
	  s   = lgr_ar(ip)%sinking
	  c   = lgr_ar(ip)%custom(1)

	  call lgr_compute_depths(ie,lb,z,hl,ht,hr)

	  write(iu) id,ty,tt,s,abs(ie)
	  write(iu) x,y,z,lb,hl
	  write(iu) (c(ii), ii=1,ncust)

	  !particles are deletedexit particles information after writing it
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
	  write(6,*) (lgr_ar(ii)%actual%ie,ii=i,min(i+9,nbdy))
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
	is_free(i) = lgr_ar(i)%actual%ie .le. 0
	is_particle(i) = lgr_ar(i)%actual%ie .gt. 0

	ifree = 0
	icopy = nbdy+1
	nbefore = nbdy

	do while( ifree .lt. icopy )

	  do i=ifree+1,icopy-1
	    if( is_free(i) ) exit
	  end do
	  ifree = i
	  if( ifree .eq. icopy ) exit

	  do i=icopy-1,ifree+1,-1
	    if( is_particle(i) ) exit
	  end do
	  icopy = i
	  if( ifree .eq. icopy ) exit

	  call copy_particle(icopy,ifree)
	  call delete_particle(icopy)

	end do

	nbdy = ifree - 1

c--------------------------------------------------------------
c check particle structure
c--------------------------------------------------------------

	!write(lunit,*) 'compress_particles: ',it,nbefore,nbdy,nbefore-nbdy

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

