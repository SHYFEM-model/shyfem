c
c $Id: baswork.f,v 1.20 2010-03-22 15:29:31 georg Exp $
c
c general framework to work on files needing basin
c
c revision log :
c
c 06.04.1999	ggu	completely restructured
c 04.06.1999	ggu	new statistics are computed
c 28.11.2005	ggu	new call to makehkv
c 31.05.2007	ggu	added area and volume frequency curve
c 24.08.2007	ggu	added new routine write_grd_from_bas
c 06.04.2009    ggu     read param.h
c 12.06.2009    ggu     areatr in double precision - new algorithm
c 01.03.2010    ggu     new routine basqual() to compute grid quality
c 22.03.2010    ggu     write external element number in basqual()
c 17.05.2011    ggu     changes in freqdep()
c 12.07.2011    ggu     better treatment of freqdep()
c 16.11.2011    ggu     basin.h introduced
c 23.01.2012    ggu     new from basinf
c 25.09.2015    ggu     prepared for nudging velocities
c 29.09.2015    ggu     finished nudging velocities
c 04.11.2015    ggu     bug in velocitiy nudging fixed
c
c****************************************************************

!==================================================================
        module mod_nudge
!==================================================================

        implicit none

	integer, parameter :: ndgdim = 50		!max size of variables
	integer, parameter :: ndgdatdim	= 10*ndgdim	!max size for BC file

	integer, save, private :: nkn_nudge = 0

	integer, save :: idsurf = 0

	integer, save :: nvars = 0
	real, save :: tramp = 0.

	real, allocatable :: andg_data(:)		!data of observations
	integer, allocatable :: ndg_nodelist(:)		!nodes of obs
	integer, allocatable :: ndg_use(:)		!use observations
	character*40, allocatable :: ndg_names(:)	!name of stations

	real, allocatable :: andg_dist(:)	!distance to be used
	real, allocatable :: andg_weight(:)	!weight to be used
	real, allocatable :: andg_obs(:)	!observations to be used

	integer, allocatable :: ndg_nodes(:)	!nodes of influence
	integer, allocatable :: ndg_area(:)	!area of influence

!==================================================================
        contains
!==================================================================

	subroutine mod_nudge_init(nkn)

	integer nkn

	if( nkn == nkn_nudge ) return

	if( nkn_nudge > 0 ) then
	  deallocate(andg_data)
	  deallocate(ndg_nodelist)
	  deallocate(ndg_use)
	  deallocate(ndg_names)
	  deallocate(andg_dist)
	  deallocate(andg_weight)
	  deallocate(andg_obs)
	  deallocate(ndg_nodes)
	  deallocate(ndg_area)
	end if

	nkn_nudge = nkn

	if( nkn == 0 ) return

	allocate(andg_data(ndgdatdim))
	allocate(ndg_nodelist(ndgdim))
	allocate(ndg_use(ndgdim))
	allocate(ndg_names(ndgdim))
	allocate(andg_dist(nkn))
	allocate(andg_weight(nkn))
	allocate(andg_obs(nkn))
	allocate(ndg_nodes(nkn))
	allocate(ndg_area(nkn))
	
	end subroutine mod_nudge_init

!==================================================================
        end module mod_nudge
!==================================================================

	subroutine init_zeta_nudging

	use mod_nudge
	use mod_nudging
	use basin

	implicit none

	logical binfl,binsert
	integer nintp,nsize,ndim
	integer k,i,node,ivar
	real ttau,sigma
	character*40 file_obs,file_stations
	character*16 cname
	character*4 cdummy

	integer ipint

	nvars = 27
	nvars = 30
	nvars = 0
	file_obs='data_level_18-27_06_2008.dat'
	file_stations='mareo_input.txt'
	tramp = 43200
	ttau = 0
	ttau = 3600
	ttau = 600
	ttau = 100
	sigma = 0
	sigma = 2000
	ivar = 0

	binfl = sigma .gt. 0.

	andgzv = 0.

	if( nvars .gt. ndgdim ) stop 'error stop nudge_init: ndgdim'

	if( nvars .le. 0 ) return

	call mod_nudge_init(nkn)

	nintp = 4
	nsize = 1
	ndim = ndgdatdim

	call exffil(file_obs,nintp,nvars,nsize,ndim,andg_data)

	ivar = 0
	open(1,file='NUDGE',status='old',form='formatted',err=1)
	read(1,*) ivar
	close(1)
    1	continue

	do i=1,nvars
	  ndg_use(i) = 1			! use all data
	  !if( mod(i,2) .ne. 0 ) ndg_use(i) = 0	! only use even vars
	  !if( i .eq. 10 ) ndg_use(i) = 0	! no Saline
	  if( i .eq. ivar ) ndg_use(i) = 0
	end do

	open(1,file=file_stations,status='old',form='formatted')
	write(6,*) 'initializing zeta nudging... nvars = ',nvars
	do i=1,nvars
	  read(1,'(a16,a4,i6)') cname,cdummy,node
	  write(6,*) i,cname,node,ndg_use(i)
	  k = ipint(node)
	  if( k .le. 0 ) goto 99
	  ndg_names(i) = cname
	  ndg_nodelist(i) = k
	end do
	write(6,*) 'finished initializing zeta nudging...'
	close(1)

	do k=1,nkn
	  andg_dist(k) = 0.
	  andg_weight(k) = 0.
	  andg_obs(k) = 0.
	  ndg_nodes(k) = 0
	  ndg_area(k) = 0
	end do

	do i=1,nvars
	  binsert = ndg_use(i) .gt. 0
	  if( binsert ) then
	    k = ndg_nodelist(i)
	    ndg_nodes(k) = k
	    ndg_area(k) = i
	    andg_weight(k) = 0.
	    if( ttau .gt. 0. ) andg_weight(k) = 1./ttau
	  end if
	end do

	if( binfl ) then
	  call nudge_influence
	  call distance(ttau,sigma)
	end if

	return
   99	continue
	write(6,*) 'Cannot find internal node: ',node
	stop 'error stop nudge_init: no internal node'
	end

c****************************************************************

	subroutine set_zeta_nudging

	use mod_nudge
	use mod_nudging
	use mod_hydro
	use basin

	implicit none

	include 'femtime.h'

	integer i,k,kk,ia
	integer nnudge,iuse,iu
	real talpha,ttau,t,zobs,zcontrib,w
	real rint(ndgdim)
	real zeta(ndgdim)
	real cont(ndgdim)

	if( nvars .le. 0 ) return

	t = it

	call exfintp(andg_data,t,rint)

	talpha = 1.
	if( tramp .gt. 0. ) talpha = (it-itanf)/tramp
	if( talpha .gt. 1. ) talpha = 1.

	do i=1,nvars
	  k = ndg_nodelist(i)
	  zeta(i) = zov(k)
	end do

	write(88,'(i10,50f7.3)') it,(rint(i),i=1,nvars)
	write(89,'(i10,50f7.3)') it,(zeta(i),i=1,nvars)

	do i=1,nvars
	  k = ndg_nodelist(i)
	  iuse = ndg_use(i)
	  if( iuse .le. 0 ) then
	    iu = 100 + i
	    write(iu,'(i10,f7.3)') it,zov(k)
	  end if
	end do

	nnudge = 0
	do k=1,nkn
	  w = andg_weight(k)
	  !kk = ndg_nodes(k)
	  ia = ndg_area(k)
	  zcontrib = 0.
	  if( ia .gt. 0 ) then
	    zobs = rint(ia)
	    zcontrib = talpha * w * ( zobs - zov(k) )
	    nnudge = nnudge + 1
	  end if
	  !zcontrib = 0.
	  andgzv(k) = zcontrib
	end do

	write(6,*) 'nudging zeta: ',nnudge

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine distance(ttau,sigma)

c sets up andg_dist and andg_weight

	use mod_nudge
	use basin

	implicit none

	real ttau		!relaxation time
	real sigma		!sigma distance

	integer k,kk,iu
	real dx,dy,dist,dist2,w,s

	real d(nkn)
	real dw(nkn)

	s = 2.*sigma*sigma

	do k=1,nkn
	  kk = ndg_nodes(k)
	  dx = xgv(kk) - xgv(k)
	  dy = ygv(kk) - ygv(k)
	  dist2 = dx*dx+dy*dy
	  dist = sqrt(dist2)
	  andg_dist(k) = dist2
	  d(k) = dist			!only for output

	  w = exp(-dist2/s)
	  dw(k) = w

	  andg_weight(k) = w / ttau
	end do

	iu = 0
	call conwrite(iu,'.dst',1,988,1,d)
	call nos_close(iu)
	close(iu)

	iu = 0
	call conwrite(iu,'.www',1,989,1,dw)
	call nos_close(iu)
	close(iu)

	end

c*******************************************************************

	subroutine nudge_influence

c sets up icol (local) and ndg_nodes and ndg_area

	use mod_nudge
	use basin

	implicit none

	integer icol(nkn)
	integer icolaux(nkn)
	real d(nkn)

	logical binsert
	integer k,ie,ii,it,is,ic,isc,i
	integer iu
	integer ngood

	do k=1,nkn
	  icol(k) = 0
	  icolaux(k) = 0
	end do

	do i=1,nvars
	  binsert = ndg_use(i) .gt. 0
	  if( binsert ) then
	    k = ndg_nodelist(i)
	    icol(k) = i
	  end if
	end do

	ngood = 0
	do k=1,nkn
	  if( icol(k) .gt. 0 ) ngood = ngood + 1
	end do

	do while( ngood .lt. nkn )

	  write(6,*) 'influence: ',ngood,nkn

	  do ie=1,nel

	    it = 0
	    is = 0
	    do ii=1,3
	      k = nen3v(ii,ie)
	      if( icol(k) .gt. 0 ) then
	        it = it + 1
	        is = is + ii
	      end if
	    end do
	
	    if( it .eq. 1 ) then
	      k = nen3v(is,ie)
	      ic = icol(k)
	      do ii=1,3
	        k = nen3v(ii,ie)
	        icolaux(k) = ic
	      end do
	    else if( it .eq. 2 ) then
	      is = 6 - is
	      isc = mod(is,3) + 1
	      k = nen3v(isc,ie)
	      ic = icol(k)
	      k = nen3v(is,ie)
	      icolaux(k) = ic
	    end if

	    if( it .gt. 0 ) then	!check
	      do ii=1,3
	        k = nen3v(ii,ie)
	        if( icol(k) .le. 0 .and. icolaux(k) .le. 0 ) then
		  write(6,*) 'internal error...... '
		  write(6,*) ie,it,is
		  stop 'error stop: internal error'
	        end if
	      end do
	    end if

	  end do

	  ngood = 0
	  do k=1,nkn
	    if( icolaux(k) .gt. 0 ) then
		if( icol(k) .le. 0 ) icol(k) = icolaux(k)
	    end if
	    if( icol(k) .gt. 0 ) ngood = ngood + 1
	  end do

	end do

	do k=1,nkn
	  ic = icol(k)
	  ndg_area(k) = ic
	  ndg_nodes(k) = ndg_nodelist(ic)
	end do

	do k=1,nkn
	  d(k) = icol(k)
	end do

	iu = 0
	call conwrite(iu,'.col',1,987,1,d)
	call nos_close(iu)
	close(iu)

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine init_velocity_nudging

	use mod_nudge
	use mod_nudging
	use intp_fem_file
	use basin

	implicit none

	integer nvar,nintp,ibc
	real vconst(2)
	real uobs_surf(nel),vobs_surf(nel)
	double precision dtime0,dtime
	character*10 what
	character*80 surffile

	integer np,lmax
	integer nodes(1)

	include 'femtime.h'

	real getpar

	dtime0 = itanf

	nodes = 0
	nvar = 2
	nintp = 2
	np = nel
	lmax = 0
	ibc = 0				!no lateral boundary
	what = 'surfvel'
	vconst = (/0.,0./)

	call getfnm(what,surffile)
	if( surffile == ' ' ) return

	taudefvel = getpar('tauvel')

        call iff_init(dtime0,surffile,nvar,np,lmax,nintp
     +                          ,nodes,vconst,idsurf)
        call iff_set_description(idsurf,ibc,'surf vel')
        call iff_need_all_values(idsurf,.false.)

        call velocity_nudging_check_data(idsurf,nvar)

	lmax = 1
	dtime = t_act
        call iff_read_and_interpolate(idsurf,dtime0)
        call iff_time_interpolate(idsurf,dtime,1,np,lmax,uobs_surf)
        call iff_time_interpolate(idsurf,dtime,2,np,lmax,vobs_surf)

	uobs(1,:) = uobs_surf(:)
	vobs(1,:) = vobs_surf(:)

	end

c*******************************************************************

	subroutine set_velocity_nudging

	use mod_nudge
	use mod_nudging
	use basin
	use levels
	use mod_internal
	use mod_hydro
	use mod_layer_thickness
	use intp_fem_file

	implicit none

	include 'femtime.h'

	integer ie,l,lmax,iflag
	real h,tau,taudef
	real u,v,s,flag
	real smax
	real uobs_surf(nel),vobs_surf(nel)
	double precision dtime

	if( idsurf <= 0 ) return

	dtime = t_act
	lmax = 1

!------------------------------------------------------------------
! read and interpolate
!------------------------------------------------------------------

	if( .not. iff_is_constant(idsurf) ) then
          call iff_read_and_interpolate(idsurf,dtime)
          call iff_time_interpolate(idsurf,dtime,1,nel,lmax,uobs_surf)
          call iff_time_interpolate(idsurf,dtime,2,nel,lmax,vobs_surf)
	  uobs(1,:) = uobs_surf(:)
	  vobs(1,:) = vobs_surf(:)
	end if

!------------------------------------------------------------------
! set relaxation time
!------------------------------------------------------------------

	taudef = taudefvel

	tauvel(1,:) = 0.

	call iff_get_flag(idsurf,flag)

	iflag = 0
	smax = 0.
        do ie=1,nel
	  u = uobs(1,ie)
	  v = vobs(1,ie)
	  if( u == flag .or. v == flag ) then
	    iflag = iflag + 1
	  else
	    s = sqrt(u*u+v*v)
	    smax = max(smax,s)
	    tauvel(1,ie) = taudef	!good point - define tau
	  end if
	end do

	!write(6,*) 'flags found.... ',iflag,nel,smax

!------------------------------------------------------------------
! add contribution to explicit term
!------------------------------------------------------------------

        do ie=1,nel
          lmax = ilhv(ie)
          do l=1,lmax
	    tau = tauvel(l,ie)
	    if( tau > 0. ) then
	      h = hdenv(l,ie)
	      fxv(l,ie) = fxv(l,ie) - (h*uobs(l,ie)-utlov(l,ie))/tau
	      fyv(l,ie) = fyv(l,ie) - (h*vobs(l,ie)-vtlov(l,ie))/tau
	    end if
	  end do
	end do

!------------------------------------------------------------------
! end of routine
!------------------------------------------------------------------

	end 

c*******************************************************************

        subroutine velocity_nudging_check_data(id,nvar)
 
	use intp_fem_file

	implicit none

        integer id
        integer nvar

	integer ivtype
        character*60 string1,string2

!       ---------------------------------------------------------
!       check nvar and get parameters
!       ---------------------------------------------------------

        if( nvar /= 2 ) then
          write(6,*) 'no support for nvar = ',nvar
          stop 'error stop velocity_nudging_check_data'
        end if

!       ---------------------------------------------------------
!       handle velocities
!       ---------------------------------------------------------

        call iff_get_var_description(id,1,string1)
        call iff_get_var_description(id,2,string2)

	string1 = adjustl(string1)
	string2 = adjustl(string2)

        if( string1 == ' ' ) then        !TS file or constant
          ivtype = 0
          if( iff_has_file(id) ) ivtype = 1

          if( ivtype == 1 ) then
            call iff_set_var_description(id,1,'velocity x')
            call iff_set_var_description(id,2,'velocity y')
          end if
        else
          if( string1 == 'velocity x' 
     +			.and. string2 == 'velocity y' ) then
            ivtype = 1
          else
            write(6,*) 'description string for velocity not recognized:'
            write(6,*) trim(string1)
            write(6,*) trim(string2)
            stop 'error stop velocity_nudging_check_data: description'
          end if
        end if

!       ---------------------------------------------------------
!       remember values and write to monitor
!       ---------------------------------------------------------

        if( ivtype == 0 ) then
          write(6,*) 'no velocity file opened'
        else
          write(6,*) 'velocity file opened: ',ivtype
          call iff_get_var_description(id,1,string1)
          call iff_get_var_description(id,2,string2)
          write(6,*) 'content: '
          write(6,*) ' 1    ',trim(string1)
          write(6,*) ' 2    ',trim(string2)
        end if

!       ---------------------------------------------------------
!       end of routine
!       ---------------------------------------------------------

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine init_nudging

	implicit none

	call init_zeta_nudging
	call init_velocity_nudging

	end 

c*******************************************************************

	subroutine set_nudging

	implicit none

	call set_zeta_nudging
	call set_velocity_nudging

	end 

c*******************************************************************
