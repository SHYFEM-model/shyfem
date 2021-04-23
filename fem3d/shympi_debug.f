
!--------------------------------------------------------------------------
!
!    Copyright (C) 2015,2017-2019  Georg Umgiesser
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

! mpi debug routines
!
! revision log :
!
! 03.06.2020	ggu	written from scratch and debugged
!
!******************************************************************

!==================================================================
        module shympi_debug
!==================================================================

	implicit none

	public

	integer, save :: iu_debug = 0

	integer, parameter :: type_integer = 1
	integer, parameter :: type_real    = 2
	integer, parameter :: type_double  = 3

        INTERFACE shympi_write_debug_record
       	MODULE PROCEDURE shympi_write_debug_record_2d_i
     +			,shympi_write_debug_record_2d_r
        END INTERFACE

!==================================================================
        contains
!==================================================================

	subroutine shympi_write_debug_init

	use shympi

	double precision dtime

	iu_debug = 199

	end

!-----------------------------------------------------------

	subroutine shympi_write_debug_time(dtime)

	use shympi

	double precision dtime

	if( .not. shympi_is_master() ) return

	write(iu_debug) dtime

	end

!-----------------------------------------------------------

	subroutine shympi_write_debug_final

	use shympi

	double precision dtime

	if( .not. shympi_is_master() ) return

	write(iu_debug) 0,0,0

	end

!-----------------------------------------------------------

	subroutine shympi_write_debug_record_2d_i(text,array)

	use shympi

	character*(*) text
	integer array(:)

	integer nn,nn_global,lmax
	integer, parameter :: record_type = type_integer
	integer, allocatable :: garray(:)
	character*80 gtext

	lmax = 1
	nn = size(array,1)
	gtext = text

	if( nn == nkn_local ) then
	  nn_global = nkn_global
	else if( nn == nel_local ) then
	  nn_global = nel_global
	else
	  write(6,*) nn,nkn_local,nel_local
	  stop 'error stop write_debug_record: nn'
	end if

	allocate(garray(nn_global))
	call shympi_exchange_array(array,garray)

	if( .not. shympi_is_master() ) return

	write(6,*) 'writing record: ',trim(text)
     +			,nn_global,lmax,record_type

	write(iu_debug) nn_global,lmax,record_type
	write(iu_debug) gtext
	write(iu_debug) garray

	end

!-----------------------------------------------------------

	subroutine shympi_write_debug_record_2d_r(text,array)

	use shympi

	character*(*) text
	real array(:)

	integer nn,nn_global,lmax
	integer, parameter :: record_type = type_real
	real, allocatable :: garray(:)
	character*80 gtext

	lmax = 1
	nn = size(array,1)
	gtext = text

	if( nn == nkn_local ) then
	  nn_global = nkn_global
	else if( nn == nel_local ) then
	  nn_global = nel_global
	else
	  write(6,*) nn,nkn_local,nel_local
	  stop 'error stop write_debug_record: nn'
	end if

	allocate(garray(nn_global))
	call shympi_exchange_array(array,garray)

	if( .not. shympi_is_master() ) return

	write(6,*) 'writing record: ',trim(text)
     +			,nn_global,lmax,record_type

	write(iu_debug) nn_global,lmax,record_type
	write(iu_debug) gtext
	write(iu_debug) garray

	end

!==================================================================
        end module shympi_debug
!==================================================================

