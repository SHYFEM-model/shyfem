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
c 05.03.2014	ggu	subroutines copied to other routine
c 16.12.2015    ggu     depth ht is now passed in for square interpol
c
c****************************************************************

        subroutine basbathy

c performs bathymetry interpolation in basin
c
c takes care of lat/lon coordinates

	use mod_depth
	use evgeom
	use basin
	use grd
	use basutil

	implicit none

	integer np
	real, allocatable :: xp(:)
	real, allocatable :: yp(:)
	real, allocatable :: dp(:)
	real, allocatable :: ap(:)

	integer node,nit
	integer n,i,nn
	integer nk,ne,nl,nne,nnl
        integer ner,nco,nknh,nelh,nli
	integer isphe

	integer mode,ike,idepth
	integer nminimum
	real ufact,umfact
	real :: flag = -999

	integer nt
        real, allocatable :: xt(:)
        real, allocatable :: yt(:)
        real, allocatable :: at(:)
        real, allocatable :: ht(:)

c-----------------------------------------------------------------
c what to do
c-----------------------------------------------------------------

	!nminimum = 5	!for dwejra
	!nminimum = 1	!minimum number of points to be used for interpolation

        !write(6,*)
        !write(6,*) 'I need the name of the basin file '
        !write(6,*) '(the file can be in GRD or BAS format)'
        !write(6,*) '(please include extension - default is GRD)'
        !write(6,*)
	!write(6,*) 'Enter file name: '
	!read(5,'(a)') gfile
        !if( gfile .eq. ' ' ) stop
	!write(6,*) 'grid is read from file : ', gfile
        !write(6,*)

        !write(6,*)
        !write(6,*) 'I need the name of the bathymetry data file '
        !write(6,*) '(the file must be in GRD format)'
        !write(6,*)
	!write(6,*) 'Enter file name: '
	!read(5,'(a)') bfile
        !if( bfile .eq. ' ' ) stop
	!write(6,*) 'Bathymetry is read from file : ', bfile
        !write(6,*)

        !write(6,*)
        !write(6,*) 'Two different algorithms are available:'
        !write(6,*) '  1   exponential interpolation (default)'
        !write(6,*) '  2   uniform interpolation on squares'
        !write(6,*) '  3   exponential interpolation (autocorrelation)'
        !write(6,*)
	!write(6,*) 'Enter choice: '
	!read(5,'(i10)') mode
	!if( mode .lt. 1 ) mode = 1
	!if( mode .gt. 3 ) mode = 1
	!write(6,*) 'Mode is : ', mode

        !write(6,*)
        !write(6,*) 'If there are some values missing you can:'
        !write(6,*) '  1   interpolate on missing depth values (default)'
        !write(6,*) '  2   interpolate on all elements/nodes'
        !write(6,*)
	!write(6,*) 'Enter choice: '
	!read(5,'(i10)') idepth
	!if( idepth .ne. 2 ) idepth = 1
	!write(6,*) 'Choice is : ', idepth

	!ike = 1
	!ufact = 1.
	!umfact = 2.	!old default
	!umfact = 3.

	!if( mode .eq. 1 .or. mode .eq. 3 ) then

        !write(6,*)
        !write(6,*) 'For the exponential algorithm you can:'
        !write(6,*) '  1   interpolate on elements (default)'
        !write(6,*) '  2   interpolate on nodes'
        !write(6,*)
	!write(6,*) 'Enter choice: '
	!read(5,'(i10)') ike
	!if( ike .ne. 2 ) ike = 1
	!write(6,*) 'Choice is : ', ike

        !write(6,*)
	!write(6,*) 'Enter parameters for expontential interpolation:'
        !write(6,*)
	!write(6,*) 'The std deviation is about the size of the elements'
	!write(6,*) 'With ufact you can ultimately correct it (default=1)'
	!write(6,*) 'The maximum radius is 3 times the standard deviation'
	!write(6,*) 'With umfact you can correct it (default=3)'
        !write(6,*)
	!write(6,*) 'Enter params ufact and umfact (<CR> for default): '
        !read(5,'(a)') line
        !n = iscanf(line,f,2)
	!if( n .lt. 0 .or. n .gt. 2 ) goto 95
        !if( n .gt. 0 ) ufact = f(1)
        !if( n .gt. 1 ) umfact = f(2)
        !write(6,*) 'ufact,umfact :',ufact,umfact
        !write(6,*)

	!end if

	mode = bmode
	idepth = 1
	if( ball ) idepth = 2
	ike = 1
	if( bnode ) ike = 2
	ufact = usfact
	umfact = uxfact

	nminimum = 1	!minimum number of points to be used for interpolation

c-----------------------------------------------------------------
c read in bathymetry file
c-----------------------------------------------------------------

	write(6,*) 'reading bathymetry file : ',bfile

	call grd_read(bfile)
	call grd_get_params(nk,ne,nl,nne,nnl)

	np = nk
	allocate(xp(np),yp(np),dp(np),ap(np))

	call grd_get_nodes(np,xp,yp,dp)
	call grd_close

c-----------------------------------------------------------------
c allocate arrays for basin
c-----------------------------------------------------------------

	nn = max(nkn,nel)
	allocate(xt(nn),yt(nn),at(nn),ht(nn))

c-----------------------------------------------------------------
c handling of depth and coordinates
c-----------------------------------------------------------------

	call get_coords_ev(isphe)
	call set_dist(isphe)

	hkv = flag
	call set_depth_i(idepth,nknh,nelh)

c-----------------------------------------------------------------
c node_test
c-----------------------------------------------------------------

	call node_test

        if( ike .eq. 1 ) then                           !elementwise
          call prepare_on_elem(nt,xt,yt,at,ht,ufact)
        else                                            !nodewise
          call prepare_on_node(nt,xt,yt,at,ht,ufact)
        end if

c-----------------------------------------------------------------
c interpolate
c-----------------------------------------------------------------

	if( mode .eq. 1 ) then
	  call interpole(np,xp,yp,dp,nt,xt,yt,at,ht,umfact,nminimum)
        else if( mode .eq. 2 ) then
	  call interpolq(np,xp,yp,dp,ht)
	else if( mode .eq. 3 ) then
	  call make_auto_corr(np,xp,yp,dp,ap,ufact)
	  call interpola(np,xp,yp,dp,ap,nt,xt,yt,at,ht)
        else
          write(6,*) 'wrong choice for mode : ',mode
          stop 'error stop'
	end if

c-----------------------------------------------------------------
c transfer depth
c-----------------------------------------------------------------

        if( ike .eq. 1 ) then                           !elementwise
	  hev(1:nel) = ht(1:nel)
        else                                            !nodewise
	  hkv(1:nkn) = ht(1:nkn)
        end if

	call transfer_depth(ike)	!copy to nodes/elems (also sets hm3v)

c-----------------------------------------------------------------
c write
c-----------------------------------------------------------------

        call basin_to_grd

        call grd_write('basbathy.grd')
        write(6,*) 'The basin has been written to basbathy.grd'

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	return
	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

        subroutine set_depth_i(idepth,nknh,nelh)

c handles depth values

        use mod_depth
        use basin, only : nkn,nel,ngr,mbw

        implicit none

        integer idepth          !1: only at missing points 2: everywhere
        integer nknh,nelh       !return - depth values found

        integer k,ie
        real flag

        flag = -999.

        if( idepth .eq. 2 ) then
	  hkv = flag
	  hev = flag
        end if

        nknh = 0
        nelh = 0

        do k=1,nkn
          if( hkv(k) .le. flag ) nknh = nknh + 1
        end do
        do ie=1,nel
          if( hev(ie) .le. flag ) nelh = nelh + 1
        end do

        end

c*******************************************************************

