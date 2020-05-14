
!--------------------------------------------------------------------------
!
!    Copyright (C) 2003,2005,2007,2009,2011,2013-2015  Georg Umgiesser
!    Copyright (C) 2018-2019  Georg Umgiesser
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

c writes time series to nc files
c
c revision log :
c
c 02.09.2003	ggu	adapted to new OUS format
c 24.01.2005	ggu	computes maximum velocities for 3D (only first level)
c 16.10.2007	ggu	new debug routine
c 27.10.2009	ggu	include evmain.h, compute volume
c 23.03.2011	ggu	compute real u/v-min/max of first level
c 25.09.2013	ggu	restructured
c 25.10.2013	ggu	changed VERS_6_1_68
c 23.12.2014	ggu	changed VERS_7_0_11
c 19.01.2015	ggu	changed VERS_7_1_3
c 17.07.2015	ggu	changed VERS_7_1_53
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 22.10.2015	ggu	changed VERS_7_3_7
c 03.04.2018	ggu	changed VERS_7_5_43
c 16.02.2019	ggu	changed VERS_7_5_60
c
c***************************************************************

	program ts2nc

c reads file and writes time series to nc file

	use mod_depth
	use mod_hydro
	use evgeom
	use levels
	use basin

	implicit none

        include 'param.h'
	include 'simul.h'

	integer ndsdim
	parameter(ndsdim=1)

        integer nvers,nin
	integer it,ie,i
        integer ierr,nread,ndry
	integer irec,maxrec
        integer nknous,nelous,nlvous
        real href,hzoff,hlvmin
	real volume
	real zmin,zmax
	real umin,umax
	real vmin,vmax

        integer ncid
        integer dimids_2d(2)
        integer dimids_1d(1)
        integer coord_varid(3)
        integer coordts_varid(2)
        integer rec_varid
        integer surge_id,tide_id
	integer date0,time0
	integer day,month,year,hour,minu,sec

	real msurge(ndsdim)
	real astro(ndsdim)
	real obs(ndsdim)

	integer nodes
	real lon(ndsdim)
	real lat(ndsdim)

	character*80 units
	character*160 std

c	integer rdous,rfous
	integer iapini,ideffi

c-----------------------------------------------------------------
c initialize basin and simulation
c-----------------------------------------------------------------

	call get_date_and_time(date0,time0)

	open(9,file='datats.dat',status='old',form='formatted')

	maxrec = 0		!max number of records to be written
	maxrec = 2		!max number of records to be written
	maxrec = 180		!max number of records to be written

	nread=0
	irec = 0

c-----------------------------------------------------------------
c prepare netcdf file
c-----------------------------------------------------------------

	nodes = 1
	lon(1) = 12.515
	lat(1) = 45.31333

        call nc_open_ts(ncid,nodes,date0,time0)
	call nc_global(ncid,descrp)

	std = 'sea_surface_height_correction_due_to_air_pressure' //
     +		'_and_wind_at_high_frequency'
	units = 'm'
	call nc_define_2d(ncid,'storm_surge',surge_id)
	call nc_define_attr(ncid,'units',units,surge_id)
	call nc_define_attr(ncid,'standard_name',std,surge_id)

	std = 'sea_surface_height_amplitude_due_to_non_equilibrium' //
     +		'_ocean_tide'
	units = 'm'
	call nc_define_2d(ncid,'astronomical_tide',tide_id)
	call nc_define_attr(ncid,'units',units,tide_id)
	call nc_define_attr(ncid,'standard_name',std,tide_id)

        call nc_end_define(ncid)
        call nc_write_coords_ts(ncid,lon,lat)

c-----------------------------------------------------------------
c loop on data of simulation
c-----------------------------------------------------------------

	it = -2*86400
	it = it - 3600

  300   continue

	read(9,*) day,month,year,hour,minu,sec,msurge(1),astro(1),obs(1)

	it= it+3600
	nread=nread+1

        irec = irec + 1
        call nc_write_time(ncid,irec,it)
        call nc_write_data_2d(ncid,surge_id,irec,1,msurge)
        call nc_write_data_2d(ncid,tide_id,irec,1,astro)

	if ( maxrec .gt. 0 .and. irec .ge. maxrec ) goto 100

	goto 300

  100	continue

c-----------------------------------------------------------------
c end of loop
c-----------------------------------------------------------------

	write(6,*)
	write(6,*) nread,' records read'
	write(6,*)

        call nc_close(ncid)

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	end

c******************************************************************

	subroutine get_date_and_time(date0,time0)

	implicit none

	integer date0,time0

	open(7,file='date0',status='old',form='formatted')
	read(7,*) date0
	close(7)
	open(8,file='time0',status='old',form='formatted')
	read(8,*) time0
	close(8)

	end

c******************************************************************

