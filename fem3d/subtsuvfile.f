
!--------------------------------------------------------------------------
!
!    Copyright (C) 2012-2016,2018-2020  Georg Umgiesser
!    Copyright (C) 2016  Erik Pascolo
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

c reading and interpolation of external files
c
c contents :
c
c revision log :
c
c 29.10.2012	ggu	created from scratch
c 05.11.2012	ggu	changed VERS_6_1_60
c 17.12.2012	ggu	changed VERS_6_1_61a
c 25.01.2013	ggu	changed VERS_6_1_62
c 03.05.2013	ggu	changed VERS_6_1_63
c 13.06.2013	ggu	changed VERS_6_1_65
c 17.06.2013	ggu	do not pass function into subroutine
c 18.06.2014	ggu	changed VERS_6_1_77
c 27.06.2014	ggu	changed VERS_6_1_78
c 02.07.2014	ggu	new framework finished
c 10.07.2014	ggu	only new file format allowed
c 18.07.2014	ggu	changed VERS_7_0_1
c 15.01.2015	ggu	changed VERS_7_1_1
c 26.02.2015	ggu	changed VERS_7_1_5
c 21.05.2015	ggu	changed VERS_7_1_11
c 13.07.2015	ggu	changed VERS_7_1_51
c 20.11.2015	ggu	changed VERS_7_3_15
c 22.02.2016	ggu&erp	new files for generic tracer (nvar>1)
c 06.06.2016	ggu	tracer_file routines changed
c 10.06.2016	ggu	changed VERS_7_5_13
c 09.09.2016	ggu	changed VERS_7_5_17
c 25.02.2018	ggu	file cleaned - time is now double
c 03.04.2018	ggu	changed VERS_7_5_43
c 16.02.2019	ggu	changed VERS_7_5_60
c 13.03.2019	ggu	changed VERS_7_5_61
c 14.02.2020	ggu	new routine ts_file_exists()
c 04.03.2020	ggu	iunit converted to id
c
c*******************************************************************	
c*******************************************************************	
c*******************************************************************	

	subroutine ts_file_descrp(id,name)
	use intp_fem_file
	implicit none
	integer id
	character*(*) name
	call iff_set_description(id,0,name)
	end

c*******************************************************************	
c*******************************************************************	
c*******************************************************************	

	subroutine ts_file_open(file,dtime,np,nlv,id)

c opens T/S file

	use intp_fem_file

	implicit none

	character*(*) file		!name of file
	double precision dtime		!initial time
	integer np			!number of points expected
	integer nlv			!vertical dimmension
	integer id			!unit number (return)

	integer nvar,nexp,lexp,nintp
	integer nodes(1)
	real vconst(1)

	nvar = 1
	nexp = np
	lexp = nlv
	nintp = 2
	nodes = 0
	vconst = 0.

!$OMP CRITICAL
	call iff_init(dtime,file,nvar,nexp,lexp,nintp
     +                                  ,nodes,vconst,id)
!$OMP END CRITICAL

	end

c*******************************************************************	

	subroutine ts_next_record(dtime,id,nlvddi,nkn,nlv,value)

	use intp_fem_file

	implicit none

	double precision dtime
	integer id
	integer nlvddi
	integer nkn
	integer nlv
	real value(nlvddi,nkn)

	integer ldim,ndim,ivar
        real vmin,vmax
	character*80 string

c--------------------------------------------------------------
c read new data
c--------------------------------------------------------------

	ivar = 1
	ndim = nkn
	ldim = nlvddi

	!write(6,*)'reading T/S values: ',dtime

	call iff_read_and_interpolate(id,dtime)
	call iff_time_interpolate(id,dtime,ivar,ndim,ldim,value)

c--------------------------------------------------------------
c some statistics
c--------------------------------------------------------------

        call conmima(nlvddi,value,vmin,vmax)

        !write(6,*) 'min/max: ',vmin,vmax

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	end

c*******************************************************************	

	subroutine ts_file_close(id)

c closes T/S file

	use intp_fem_file

	implicit none

	integer id

	call iff_forget_file(id)

	end

c*******************************************************************	

	subroutine ts_file_exists(file,bexist)

c checks if file exists (and no read error)

	use intp_fem_file

	implicit none

	character*(*) file		!name of file
	logical bexist

	call iff_file_exists(file,bexist)

	end

c*******************************************************************	
c*******************************************************************	
c*******************************************************************	
c****** new routines ***********************************************	
c*******************************************************************	
c*******************************************************************	
c*******************************************************************	

	subroutine tracer_file_open(file,dtime,nvar,np,nlv,val0,id)

c opens tracer file

	use intp_fem_file

	implicit none

	character*(*) file		!name of file
	double precision dtime		!time
	integer nvar			!number of (state) variables
	integer np			!number of points expected
	integer nlv			!number of vertical levels
	real val0(nvar)			!default initial condition
	integer id			!id of file (return)

	integer nexp,lexp,nintp
	integer nodes(nvar)

	nexp = np
	lexp = nlv
	nintp = 2
	nodes = 0

!$OMP CRITICAL
	call iff_init(dtime,file,nvar,nexp,lexp,nintp
     +                                  ,nodes,val0,id)
!$OMP END CRITICAL

	if( id <= 0 ) then
	  write(6,*) 'Cannot open file: ',file
	  stop 'error stop tracer_file_open: error file open'
	end if

	end

c*******************************************************************	

	subroutine tracer_file_next_record(dtime,id
     +					,nvar,nlvddi,nkn,nlv,value)

c reads next record of tracer

	use intp_fem_file

	implicit none

	double precision dtime
	integer id
	integer nvar
	integer nlvddi
	integer nkn
	integer nlv
	real value(nlvddi,nkn,nvar)

	integer ldim,ndim,ivar
	character*80 string

c--------------------------------------------------------------
c read new data
c--------------------------------------------------------------

	ndim = nkn
	ldim = nlvddi

	call iff_read_and_interpolate(id,dtime)
	do ivar=1,nvar
	  call iff_time_interpolate(id,dtime,ivar,ndim,ldim
     +					,value(:,:,ivar))
	end do

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	end

c*******************************************************************	

	subroutine tracer_file_descrp(id,text)

c sets description for file

	use intp_fem_file

	implicit none

	integer id
	character*(*) text

	call iff_set_description(id,0,text)

	end

c*******************************************************************	

	subroutine tracer_file_close(id)

c closes tracer file

	use intp_fem_file

	implicit none

	integer id

	call iff_forget_file(id)

	end

c*******************************************************************	
c*******************************************************************	
c*******************************************************************	

	subroutine tracer_file_init(what,file_init,dtime
     +				,nvar,nlvddi,nlv,nkn,val0,val)

c initialization of tracer from file

        implicit none

	character*(*) what
	character*(*) file_init
        double precision dtime
        integer nvar
        integer nlvddi
        integer nlv
        integer nkn
        real val0(nvar)			!default for vals if no file is given
        real val(nlvddi,nkn,nvar)

        integer id,iv
        character*80 file

        call getfnm(file_init,file)

	do iv=1,nvar
          val(:,:,iv) = val0(iv)
	end do

        if( file == ' ' ) return

        write(6,*) 'tracer_init: opening file for ',trim(what)
        write(6,*) '   file name: ',trim(file)
        write(6,*) '   variables: ',nvar

        call tracer_file_open(file,dtime,nvar,nkn,nlv,val0,id)
        call tracer_file_descrp(id,what)
        call tracer_file_next_record(dtime,id,nvar,nlvddi,nkn,nlv,val)
        call tracer_file_close(id)

	write(6,*) 'tracer_init: successful init for ',trim(what)

	end

c*******************************************************************	
c*******************************************************************	
c*******************************************************************	

