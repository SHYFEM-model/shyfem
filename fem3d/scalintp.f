c
c $Id: basbathy.f,v 1.7 2010-02-26 17:35:06 georg Exp $
c
c revision log :
c
c 06.04.1999	ggu	completely restructured
c 04.06.1999	ggu	new statistics are computed
c 08.09.2003	ggu	mode 5 -> write depth values from elements
c 23.09.2004    ggu     interpolq() changed for bathy interpolation
c 02.10.2004    ggu     interpole() for exponential interpolation
c 12.05.2005    ggu     pass hmin to interpolation functions
c 06.04.2009    ggu     read param.h
c 24.04.2009	ggu	new call to rdgrd()
c 21.05.2009	ggu	restructured to allow for nodal interpolation
c 16.12.2010	ggu	bug fix in transfer_depth()
c 02.12.2011	ggu	introduction of nminimum - hardcoded for now
c 16.03.2012	ggu	autoregression introduced (make_auto_corr,interpola)
c 16.03.2012	ggu	default value for umfact set to 3, new mode = 3
c 01.06.2012	ggu	some more changes
c 13.06.2013	ggu	copy_depth() renamed to transfer_depth()
c 13.02.2014	ggu	new data written, can read also bas file
c
c****************************************************************

        program scalintp

c performs scalar interpolation in basin
c
c takes care of lat/lon coordinates

	use mod_depth
	use evgeom
	use basin

	implicit none

	include 'param.h'
	include 'simul.h'
	include 'pkonst.h'

	integer ndim
	parameter(ndim=1200000)
	real xp(ndim)
	real yp(ndim)
	real dp(ndim)
	real ap(ndim)

        character*40 bfile,gfile,nfile
        character*60 line
	integer node,nit,ivar
	integer mode,np,n,i,k,nt
        integer ner,nco,nknh,nelh,nli
	integer nlidim,nlndim
	integer ike,idepth
	integer nminimum
	integer isphe
	real ufact,umfact
	real flag
	real f(5)
	logical bstop,bbasin
	integer iscanf

	real, allocatable :: xt(:)
	real, allocatable :: yt(:)
	real, allocatable :: at(:)
	real, allocatable :: ht(:)

c-----------------------------------------------------------------
c what to do
c-----------------------------------------------------------------

	nminimum = 5	!for dwejra
	nminimum = 1	!minimum number of points to be used for interpolation

        write(6,*)
        write(6,*) 'I need the name of the basin file '
        write(6,*)
	write(6,*) 'Enter file name: '
	read(5,'(a)') gfile
        if( gfile .eq. ' ' ) stop
	write(6,*) 'basin is read from file : ', gfile
        write(6,*)

        write(6,*)
        write(6,*) 'I need the name of the data file '
        write(6,*) '(the file must be in GRD or node/value format)'
        write(6,*)
	write(6,*) 'Enter file name: '
	read(5,'(a)') bfile
        if( bfile .eq. ' ' ) stop
	write(6,*) 'Data is read from file : ', bfile
        write(6,*)

        write(6,*)
        write(6,*) 'Different algorithms are available:'
        write(6,*) '  1   exponential interpolation (default)'
        write(6,*) '  2   uniform interpolation on squares'
        write(6,*) '  3   exponential interpolation (autocorrelation)'
        write(6,*)
	write(6,*) 'Enter choice: '
	read(5,'(i10)') mode
	if( mode .lt. 1 ) mode = 1
	if( mode .gt. 3 ) mode = 1
	write(6,*) 'Mode is : ', mode

	idepth = 2	!interpolate everywhere
	ike = 2		!interpolate on nodes

	ufact = -1.
	umfact = 2.	!old default
	umfact = 3.

	if( mode .eq. 1 .or. mode .eq. 3 ) then

        write(6,*)
	write(6,*) 'Enter parameters for expontential interpolation:'
        write(6,*)
	write(6,*) 'The std deviation is about the size of the elements'
	write(6,*) 'With ufact you can ultimately correct it (default=1)'
	write(6,*) 'The maximum radius is 3 times the standard deviation'
	write(6,*) 'With umfact you can correct it (default=3)'
        write(6,*)
	write(6,*) 'Enter params ufact and umfact (<CR> for default): '
        read(5,'(a)') line
        n = iscanf(line,f,2)
	if( n .lt. 0 .or. n .gt. 2 ) goto 95
        if( n .gt. 0 ) ufact = f(1)
        if( n .gt. 1 ) umfact = f(2)
        write(6,*) 'ufact,umfact :',ufact,umfact
        write(6,*)

	end if

c-----------------------------------------------------------------
c read in data file
c-----------------------------------------------------------------

	write(6,*) 'reading data file : ',bfile
	np = ndim
	call readgrd(bfile,np,xp,yp,dp)

c-----------------------------------------------------------------
c read in basin
c-----------------------------------------------------------------

	call check_basin_name(gfile,bbasin)

	write(6,*) 'reading basin as bas file...'
	call read_basin(gfile)

	allocate(xt(nel),yt(nel),at(nel),ht(nel))

c-----------------------------------------------------------------
c handling of depth and coordinates
c-----------------------------------------------------------------

	call check_spheric_ev			!sets lat/lon flag
	call get_coords_ev(isphe)
	call set_dist(isphe)
	call set_ev

c-----------------------------------------------------------------
c general info
c-----------------------------------------------------------------

        write(6,*)
        write(6,*) ' nkn  = ',nkn, '  nel  = ',nel
        write(6,*)

c-----------------------------------------------------------------
c prepare points where to interpolate
c-----------------------------------------------------------------

	flag = -999.
	do k=1,nkn
	  hkv(k) = flag
	end do

        if( ike .eq. 1 ) then           		!elementwise
          call prepare_on_elem(nt,xt,yt,at,ht,ufact)
        else						!nodewise
          call prepare_on_node(nt,xt,yt,at,ht,ufact)
        end if

c-----------------------------------------------------------------
c interpolate
c-----------------------------------------------------------------

	if( mode .eq. 1 ) then
	  call interpole(np,xp,yp,dp,nt,xt,yt,at,ht,umfact,nminimum)
        else if( mode .eq. 2 ) then
	  write(6,*) 'mode = 2 -> cannot use this interpolation mode'
	  stop 'error stop scalintp: mode 2 not possible'
	  !call interpolq(np,xp,yp,dp)
	else if( mode .eq. 3 ) then
	  call make_auto_corr(np,xp,yp,dp,ap,ufact)
	  call interpola(np,xp,yp,dp,ap,nt,xt,yt,at,ht)
        else
          write(6,*) 'wrong choice for mode : ',mode
          stop 'error stop'
	end if

c-----------------------------------------------------------------
c copy to depth array
c-----------------------------------------------------------------

	do k=1,nkn
	  hkv(k) = ht(k)
	end do

c-----------------------------------------------------------------
c write
c-----------------------------------------------------------------

	nfile = 'scalintp.grd'
	open(1,file=nfile,status='unknown',form='formatted')
	call wrgrd(1,ike)
	close(1)
        write(6,*) 'The data has been written to ',nfile

	ivar = 0
	call write_nos('scalintp.nos',nkn,nel,ivar,ht,hev)

	call write_data('scalintp.dat',nkn,hkv)
	!call write_xy('basbathy.xyz',nkn,ipv,xgv,ygv)

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	stop
   95	continue
	write(6,*) n,(f(i),i=1,n)
	write(6,*) line
	stop 'error stop basbathy: error in parameters'
	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

