
!--------------------------------------------------------------------------
!
!    Copyright (C) 1999,2003-2005,2007-2011,2013-2014  Georg Umgiesser
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

c revision log :
c
c 12.02.1999	ggu	adapted to auto mode
c 30.10.2003	ggu	added plowind()
c 22.03.2004	ggu	lagrang model output added (mode 12)
c 04.10.2004	ggu	temp and salt from nos file
c 02.12.2004	ggu	copyright statement added
c 02.03.2005	ggu	introduced flag for scalar plot
c 18.10.2005	ggu	data structures introduced to call set_geom
c 14.03.2007	ggu	wave plotting
c 16.04.2008	ggu	new Makefile structure
c 09.12.2008	ggu	changes from Malta integrated (annotation)
c 26.01.2009	ggu	new makedepend
c 06.04.2009	ggu	new param.h structure
c 20.04.2009	ggu	level.h eliminated
c 09.10.2009	ggu	plotting of atmos. pressure
c 13.10.2009	ggu	new arrays for velocity section plot
c 23.02.2010	ggu	new call to set_default_color_table()
c 23.03.2010	ggu	changed v6.1.1
c 14.04.2010	ggu	changed v6.1.4
c 08.10.2010	ggu	changed VERS_6_1_13
c 17.12.2010	ggu	substituted hv with hkv
c 31.03.2011	ggu	new arrays fvlv, arfvlv for scalar plotting
c 14.04.2011	ggu	changed VERS_6_1_22
c 31.08.2011	ggu	new copyright, eos plotting
c 01.09.2011	ggu	changed VERS_6_1_32
c 14.11.2011	ggu	new array hl for depth structure
c 22.11.2011	ggu	changed VERS_6_1_37
c 09.12.2011	ggu	changed VERS_6_1_38
c 24.01.2012	ggu	changed VERS_6_1_41
c 27.02.2013	ggu	deleted amat, better color handling
c 03.05.2013	ggu	changed VERS_6_1_63
c 13.06.2013	ggu	new plotting for fem files
c 05.09.2013	ggu	better handling of variable
c 12.09.2013	ggu	changed VERS_6_1_67
c 28.01.2014	ggu	changed VERS_6_1_71
c 05.03.2014	ggu	bug fix: isphe was not stored
c 30.05.2014	ggu	changed VERS_6_1_76
c 18.07.2014	ggu	changed VERS_7_0_1
c 20.10.2014	ggu	handle absolute and plotting time
c 23.12.2014	ggu	changed VERS_7_0_11
c 19.01.2015	ggu	changed VERS_7_1_2
c 19.01.2015	ggu	changed VERS_7_1_3
c 26.02.2015	ggu	changed VERS_7_1_5
c 10.07.2015	ggu	changed VERS_7_1_50
c 17.07.2015	ggu	changed VERS_7_1_52
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 24.07.2015	ggu	changed VERS_7_1_82
c 30.07.2015	ggu	changed VERS_7_1_83
c 29.09.2015	ggu	changed VERS_7_2_5
c 28.04.2016	ggu	changed VERS_7_5_9
c 25.05.2016	ggu	changed VERS_7_5_10
c 30.05.2016	ggu	changed VERS_7_5_11
c 17.06.2016	ggu	changed VERS_7_5_15
c 12.01.2017	ggu	changed VERS_7_5_21
c 31.03.2017	ggu	changed VERS_7_5_24
c 09.05.2017	ggu	changed VERS_7_5_26
c 14.11.2017	ggu	changed VERS_7_5_36
c 25.10.2018	ggu	changed VERS_7_5_51
c 18.12.2018	ggu	changed VERS_7_5_52
c 21.05.2019	ggu	changed VERS_7_5_62
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
!	if( mode .eq. 12 ) call plolagr
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
	  stop 'no interactivity anymore...'
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
          !iwhat = ideflt(iwhat,'Enter choice : ')
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

