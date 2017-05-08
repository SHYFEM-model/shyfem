c
c $Id: plotsim.f,v 1.53 2010-03-11 15:33:09 georg Exp $
c
c revision log :
c
c 12.02.1999	ggu	adapted to auto mode
c 30.10.2003    ggu     added plowind()
c 22.03.2004	ggu	lagrang model output added (mode 12)
c 04.10.2004	ggu	temp and salt from nos file
c 02.12.2004    ggu     copyright statement added
c 02.03.2005    ggu     introduced flag for scalar plot
c 18.10.2005    ggu     data structures introduced to call set_geom
c 14.03.2007    ggu     wave plotting
c 16.04.2008    ggu     new Makefile structure
c 09.12.2008    ggu     changes from Malta integrated (annotation)
c 26.01.2009	ggu	new makedepend
c 06.04.2009	ggu	new param.h structure
c 20.04.2009	ggu	level.h eliminated
c 09.10.2009	ggu	plotting of atmos. pressure
c 13.10.2009	ggu	new arrays for velocity section plot
c 23.02.2010	ggu	new call to set_default_color_table()
c 17.12.2010	ggu	substituted hv with hkv
c 31.03.2011	ggu	new arrays fvlv, arfvlv for scalar plotting
c 31.08.2011	ggu	new copyright, eos plotting
c 14.11.2011	ggu	new array hl for depth structure
c 27.02.2013	ggu	deleted amat, better color handling
c 13.06.2013	ggu	new plotting for fem files
c 05.09.2013	ggu	better handling of variable
c 05.03.2014	ggu	bug fix: isphe was not stored
c 20.10.2014	ggu	handle absolute and plotting time
c
c*************************************************************

	program plotsim

c plots simulation

	use mod_hydro_plot
	use mod_geom
	use mod_depth
	use mod_hydro_print
	use mod_hydro_vel
	use mod_hydro
	use evgeom
	use levels
	use basin

	implicit none

c parameters

	include 'mkonst.h'
	include 'pkonst.h'

c local
	character*20 what
	integer mode,ivar
	integer ie,ii,k,l,i
	integer isphe
	integer date,time
	integer iapini
	real eps
	real sflag
	real getpar

c----------------------------------------------
c copyright
c----------------------------------------------

	call shyfem_copyright('plotsim - plotting maps on FE grids')

c----------------------------------------------
c initialize parameters
c----------------------------------------------

	eps=1.e-5
	eps1=1.e-5
	eps2=1.e-6
	pi=3.141592653
	flag = -9988765.0
	flag = -999.0
	sflag = -999.0
	high=1.e+35
	grav=9.81

	call set_flag(sflag)

c----------------------------------------------
c read basin
c----------------------------------------------

	if(iapini(7,0,0,0).eq.0) then
		stop 'error stop : iapini'
	end if
	write(6,*) 'Basin:   nkn = ',nkn,'  nel = ',nel

c----------------------------------------------
c set up elements
c----------------------------------------------

	call allocate_2d_arrays(nel)

	isphe = nint(getpar('isphe'))
	call set_coords_ev(isphe)
	call set_ev
	call set_geom
	call get_coords_ev(isphe)
	call putpar('isphe',float(isphe))

c----------------------------------------------
c make depth on nodes and elements
c----------------------------------------------

	call makehev(hev)
	call makehkv(hkv)

c----------------------------------------------
c time management
c----------------------------------------------

	call ptime_init
	call get_date_time(date,time)
	call ptime_set_date_time(date,time)
	call ptime_min_max
	call iff_init_global_2d(nkn,nel,hkv,hev,date,time)	!FIXME

c----------------------------------------------
c interactive set up
c----------------------------------------------

	call init_plot

	call ichoice(mode,ivar)

	call read_apn_file(ivar)
	call initialize_color

c----------------------------------------------
c open plot
c----------------------------------------------

	call qopen

c----------------------------------------------
c do plotting
c----------------------------------------------

	if( mode .eq. 1 )  call plobas
	if( mode .eq. 2 )  call plosim(.true.)
	if( mode .eq. 3 )  call plosim(.false.)
	if( mode .eq. 4 )  call plozet
c	if( mode .eq. 4 )  call plobar
	if( mode .eq. 5 )  call plonos('.con',10)
	if( mode .eq. 6 )  call plonos('.nos',12)
	if( mode .eq. 7 )  call plonos('.nos',11)
	if( mode .eq. 8 )  call plonos('.rms',18)
	if( mode .eq. 9 )  call plonos('.oxy',15)
	if( mode .eq. 10 ) call plonos('.nos',0)
	if( mode .eq. 11 ) call plofem('.fem',21)	!wind
	if( mode .eq. 12 ) call plolagr
	if( mode .eq. 13 ) call plowave
	if( mode .eq. 14 ) call plofem('.fem',20)	!pressure
	if( mode .eq. 15 ) call ploeos('.eos',0)
	if( mode .eq. 16 ) call plofem('.fem',0)

c----------------------------------------------
c close plot
c----------------------------------------------

	call qclose

c----------------------------------------------
c end of routine
c----------------------------------------------

	end

c*****************************************************************

	subroutine ichoice(mode,ivar)

	implicit none

	integer mode
	integer ivar

	integer ndim
	parameter (ndim=16)

	character*10 what
	character*20 string
	integer iwhat,iauto,isub
	integer ideflt
	real getpar

	character*10 whats(0:ndim)
	save whats
	data whats /' ','bath','vel','trans','zeta','conz'
     +			,'temp','salt','rms','oxygen','nos'
     +			,'wind','lgr','wave','pres','elem'
     +			,'fem'/
	
	iwhat = nint(getpar('iwhat'))
	iauto = nint(getpar('iauto'))

        if( iauto .eq. 0 .or. iwhat .eq. 0 ) then
	  write(6,*)
	  write(6,*) ' basin ...................  1'
	  write(6,*) ' velocity ................  2'
	  write(6,*) ' transport ...............  3'
	  write(6,*) ' water level .............  4'
	  write(6,*) ' concentration ...........  5'
	  write(6,*) ' temperature .............  6'
	  write(6,*) ' salinity ................  7'
	  write(6,*) ' rms .....................  8'
	  write(6,*) ' oxygen ..................  9'
	  write(6,*) ' generic scalar value .... 10'
	  write(6,*) ' wind field .............. 11'
	  write(6,*) ' lagrangian .............. 12'
	  write(6,*) ' wave .................... 13'
	  write(6,*) ' atmospheric pressure .... 14'
	  write(6,*) ' generic element values .. 15'
	  write(6,*) ' generic fem file ........ 16'
	  write(6,*)
          iwhat = ideflt(iwhat,'Enter choice : ')
        else
          write(6,*) 'Plotting with choice: ',iwhat
	  write(6,*)
        end if

	if( iwhat .lt. 1 .or. iwhat .gt. ndim ) then
	  write(6,*) 'iwhat = ',iwhat
	  stop 'error stop ichoice: iwhat'
	end if

	mode = iwhat
	what = whats(iwhat)
	call string2ivar(what,ivar)
	call checkvar(ivar)
	call ivar2string(ivar,string,isub)

	write(6,*) 'Plotting ',ivar,what,string

	end

c*****************************************************************

	subroutine initialize_color0

	implicit none

	integer icolor
	real getpar

	call colsetup
	icolor = nint(getpar('icolor'))
	call set_color_table( icolor )
	call set_default_color_table( icolor )

	end

c*****************************************************************

        subroutine initialize_color

        implicit none

        integer icolor
        real getpar
        logical has_color_table

        call colsetup
        call admin_color_table

        if( has_color_table() ) call putpar('icolor',8.)

        icolor = nint(getpar('icolor'))
        call set_color_table( icolor )
        call set_default_color_table( icolor )

        call write_color_table

        end

c*****************************************************************

        subroutine get_date_time(date,time)

        implicit none

        integer date,time

        double precision dgetpar

        date = nint(dgetpar('date'))
        time = nint(dgetpar('time'))

        end

c*****************************************************************

	subroutine allocate_2d_arrays(npd)

	use mod_hydro_plot
	use mod_geom
	use mod_depth
	use evgeom
	use basin, only : nkn,nel,ngr,mbw
	use levels

	implicit none

	integer npd

	integer np,nlvaux

	np = max(nel,npd)

	call ev_init(nel)
	call mod_geom_init(nkn,nel,ngr)

	call mod_depth_init(nkn,nel)

	nlvaux = 1
	call mod_hydro_plot_init(nkn,nel,nlvaux,np)

	write(6,*) 'allocate_2d_arrays: ',nkn,nel,ngr,np

	end

c*****************************************************************

	subroutine reallocate_2d_arrays(npd)

	use mod_hydro_plot
	use basin, only : nkn,nel,ngr,mbw
	use levels

	implicit none

	integer npd

	integer np

	np = max(nel,npd)
	call mod_hydro_plot_init(nkn,nel,nlv,np)

	write(6,*) 'reallocate_2d_arrays: ',nkn,nel,np

	end

c*****************************************************************

