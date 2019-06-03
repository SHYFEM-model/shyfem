
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

c handles open boundary conditions for scalar variables
c
c contents :
c
c subroutine bnds_init(text,file,nintp,nvar,ndim,array,aconst)
c			initializes boundary condition
c subroutine bnds_set(text,t,ndim,array,aaux)
c			sets boundary condition
c subroutine bnds_trans(text,ndim,array,aaux,ivar,nlvddi,r3v)
c			transfers boundary condition to matrix
c subroutine bnds_set_def(text,ndim,array)
c			sets default value for boundaries
c subroutine bnds_print(text,ndim,array)
c			prints boundary condition
c
c revision log :
c
c 10.08.2003	ggu	exit if no file can be opened
c 03.03.2005	ggu	new call to exfini
c 01.02.2006	ggu	bugfix -> nextra was 3 -> not used anymore
c 03.02.2006	ggu	bugfix -> calling exfpres
c 17.02.2006	ggu	included text for debug, new call to get_bflux()
c 23.03.2006	ggu	changed time step to real
c 05.10.2007	ggu	definition of array(ndim,1) changed to array(ndim,0:1)
c 17.03.2008	ggu	routines re-arranged, new bnds_trans, bnds_set_def
c 17.04.2008	ggu	deleted bnds_set_global
c 23.04.2008	ggu	in bnds_set_def() eliminated aaux
c 23.03.2010	ggu	changed v6.1.1
c 14.02.2012	ggu	changed VERS_6_1_44
c 16.02.2012	ggu	new routine bnds_init0 (to force spatially const bound)
c 30.03.2012	ggu	changed VERS_6_1_51
c 13.06.2013	ggu	changed VERS_6_1_65
c 05.12.2013	ggu	changed VERS_6_1_70
c 25.06.2014	ggu	new routines bnds_init_new() and bnds_trans_new()
c 07.07.2014	ggu	changed VERS_6_1_79
c 10.07.2014	ggu	only new file format allowed
c 18.07.2014	ggu	changed VERS_7_0_1
c 05.11.2014	ggu	changed VERS_7_0_5
c 05.02.2015	ggu	check for number of variables read
c 10.02.2015	ggu	new routine bnds_read_new()
c 26.02.2015	ggu	changed VERS_7_1_5
c 10.07.2015	ggu	changed VERS_7_1_50
c 30.09.2015	ggu	new routine iff_flag_ok() for ambient value
c 23.10.2015	ggu	changed VERS_7_3_9
c 01.04.2016	ggu	changed VERS_7_5_7
c 19.04.2018	ggu	changed VERS_7_5_45
c 14.02.2019	ggu	changed VERS_7_5_56
c 16.02.2019	ggu	changed VERS_7_5_60
c
c******************************************************************

	subroutine bnds_init_new(what,dtime0,nintp,nvar,nkn,nlv
     +					,cdef,ids)

c initializes boundary condition for scalars

	use intp_fem_file

	implicit none

	character*(*) what	!what variable to initialize
	double precision dtime0	!time
	integer nintp		!degree of interpolation - same for all bounds
	integer nvar		!number of variables
	integer nkn		!number of points (max)
	integer nlv		!number of vertical levels
	real cdef(nvar)		!default values
	integer ids(*)		!boundary info id (return)

	character*80 file
	integer nbc,ibc
	integer nk,nsize
	integer iunit,n,i,id
	real val

	integer nodes(nkn)
	real aconst(nvar)

	integer nbnds,nkbnds,ifileo,kbnds
	integer nvar_orig
	logical exists_bnd_par
	logical bdebug

	bdebug = .true.
	bdebug = .false.

	if( bdebug ) then
	  write(6,*) '-------------------------'
	  write(6,*) 'Initialization for scalar:'
	  write(6,*) what,nvar
	  write(6,*) '-------------------------'
	end if

	nvar_orig = nvar
	nbc = nbnds()

	do ibc=1,nbc

	  ids(ibc) = 0
          nk = nkbnds(ibc)
	  if( nk .le. 0 ) cycle
          do i=1,nk
            nodes(i) = kbnds(ibc,i)
          end do

	  call get_boundary_file(ibc,what,file)

	  aconst = cdef
	  if( exists_bnd_par(what) ) then
	    call get_bnd_par(ibc,what,val)
	    aconst = val
	  end if

          call iff_init(dtime0,file,nvar,nk,nlv,nintp
     +                          ,nodes,aconst,id)
	  if( nvar /= nvar_orig ) goto 99
	  call iff_set_description(id,ibc,what)
	  call iff_flag_ok(id)		!can deal with ambient value

	  ids(ibc) = id

	  if( bdebug ) then
            write(6,*) 'boundary: ',ibc,id,what
            if( file .ne. ' ' ) then
	      write(6,*) '  file: ',file(1:60)
	    else
              write(6,*) '  def: ',aconst
	    end if
	  end if

        end do

	if( bdebug ) then
	  do ibc=1,nbc
	    call iff_print_info(ids(ibc))
	  end do
	end if

	return
   99	continue
	write(6,*) 'number of variables is wrong'
	write(6,*) 'expected: ',nvar_orig
	write(6,*) 'read from file: ',nvar
	write(6,*) 'type of boundary condition: ',what
	write(6,*) 'number of boundary: ',ibc
	write(6,*) 'file name: ',trim(file)
	call iff_print_info(id)
	stop 'error stop bnds_init_new: wrong number of variables'
	end

c******************************************************************

	subroutine bnds_read_new(text,ids,dtime)

c reads new boundary condition

	use intp_fem_file
	use iso8601

	implicit none

	character*(*), intent(in)	:: text		!for debug
	integer, intent(in)		:: ids(*)
	double precision, intent(in)	:: dtime

	integer nbc,ibc,id
	integer nbnds

	integer date,time,ierr
	double precision atime,atime0,dctime
	character*80 string

	nbc = nbnds()

	do ibc=1,nbc
	  id = ids(ibc)
	  if( id .le. 0 ) cycle
	  call iff_read_and_interpolate(id,dtime)
	end do

	return

! debug code

	string = '2018-12-31::22:40:00'
	call convert_to_dtime(string,dctime)
	if( dtime < dctime ) return

	do ibc=1,nbc
	  id = ids(ibc)
	  if( id .le. 0 ) cycle
	  call iff_info_on_data(id)
	end do

	end

c******************************************************************

	subroutine bnds_trans_new(text,ids,dtime,ivar,nkn,nlv,nlvddi,r3v)

c transfers boundary condition to matrix

	use intp_fem_file

	implicit none

	character*(*) text	!text for debug
	integer ids(*)
	double precision dtime
	integer ivar		!variable to use (can be 0 -> 1)
	integer nkn
	integer nlv
	integer nlvddi		!vertical dimension of levels
	real r3v(nlvddi,nkn)	!matrix to which BC values are transfered

	integer nbc,ibc
	integer nvar,nsize,ndata
	integer nk,iv,kn
	integer i,ip,id
	real t

	real vals(nlv,nkn)

	integer nbnds,nkbnds,kbnds

	call init_scal_bc(r3v)	!sets r3v to flag - dims nlv,nkn

	nbc = nbnds()
	iv = max(ivar,1)

	do ibc=1,nbc

	  nk = nkbnds(ibc)   !total number of nodes of this boundary
	  id = ids(ibc)
	  if( id .le. 0 ) cycle

	  call iff_time_interpolate(id,dtime,iv,nkn,nlv,vals)

	  do i=1,nk
	    kn = kbnds(ibc,i)
	    if( kn <= 0 ) cycle
	    call dist_3d(nlvddi,r3v,kn,nlv,vals(1,i))
	  end do

	end do

	!call iff_print_info(44,0,.true.)
	!call iff_print_boundary_info(9,0,.true.)

	end

c******************************************************************

