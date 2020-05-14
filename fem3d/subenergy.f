
!--------------------------------------------------------------------------
!
!    Copyright (C) 2019  Georg Umgiesser
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

c routines to compute and write out nenergy (kinetic, potential, total)
c
c contents :
c
c revision log :
c
c 14.12.2019	ggu	started from scratch (still not finished)
c
c notes :
c
c********************************************************************
c********************************************************************
c********************************************************************
c********************************************************************
c********************************************************************
c
c********************************************************************

!====================================================================
	module energy
!====================================================================

	implicit none


	logical, parameter :: baver = .true.	!use avraged values
	integer, parameter :: nvar = 3		!number of variables

        integer, save :: iuinfo = 0

	real, save, allocatable :: ener(:,:,:)	!energy data
	double precision, save :: raccum
	double precision, save, allocatable :: aener(:,:,:)	!accumulated

	double precision, save :: da_out(4)

!====================================================================
	end module energy
!====================================================================

	subroutine handle_energy

c energy module

	use levels
	use basin
	use energy
	use shympi
 
	implicit none

        character*10 what

        real tpstot(npstate)              !for mass test
        real tsstot(nsstate)

        logical has_output_d,next_output_d

        what = 'energy'

c-------------------------------------------------------------------
c initialization
c-------------------------------------------------------------------

	if( icall .le. -1 ) return

	if( icall .eq. 0 ) then
	  iener = nint(getpar('iener'))
	  if( iener .le. 0 ) icall = -1
	  if( icall .le. -1 ) return
	  icall = 1

          call get_first_dtime(dtime0)

          if( iuinfo .eq. 0 ) then
            call getinfo(iuinfo)  !unit number of info file
          end if

c         --------------------------------------------------
c	  initialize output 
c         --------------------------------------------------

	  call energy_init_file_output
	  call energy_write_file_output(dtime0)

	  write(6,*) 'energy module initialized...'

	end if

c-------------------------------------------------------------------
c normal call
c-------------------------------------------------------------------

c	-------------------------------------------------------------------
c	compute bulk quantities and write to inf file
c	-------------------------------------------------------------------

        call energ3d(kenergy,penergy,ksurf,-1)

        kenergy = shympi_sum(kenergy)
        penergy = shympi_sum(penergy)
        ksurf = shympi_sum(ksurf)
        tenergy = kenergy + penergy

        if(shympi_is_master()) then
          call get_act_timeline(aline)
          write(iuinfo,1000) ' energy: ',aline
     +                          ,kenergy,penergy,tenergy,ksurf
 1000     format(a,a20,4e12.4)
        end if

	if( .not. has_output_d(da_out) return

c	-------------------------------------------------------------------
c	write of results
c	-------------------------------------------------------------------

	call energy_write_file_output(dtime)

c	-------------------------------------------------------------------
c	debug output
c	-------------------------------------------------------------------

c	-------------------------------------------------------------------
c	end of routine
c	-------------------------------------------------------------------

	end

c*************************************************************
c*************************************************************
c*************************************************************

        subroutine energy_init_file_output

        use basin
        use levels
        use energy

        implicit none

        integer id
        logical has_output_d

        call init_output_d('itmene','idtene',da_out)
        if( has_output_d(da_out) ) then
          call shyfem_init_scalar_file('energy',nvar,.false.,id)
          da_out(4) = id

	  allocate(ener(nlvdi,nkndi,nvar))
	  allocate(aener(nkndi,nkndi,nvar))
	  ener = 0.
	  aener = 0.
	  raccum = 0.
        end if

        write(6,*) 'energy init file output done'

        end 

c*************************************************************

	subroutine energy_write_file_output(dtime)

	use basin
	use levels
	use energy

	implicit none

	double precision dtime

	integer id,i
	logical next_output_d

        if( next_output_d(da_out) ) then
          id = nint(da_out(4))
          do i=1,nvar
            idc = 250 + i
            call shy_write_scalar_record(id,dtime,idc,nlvdi
     +                                          ,ener(1,1,i))
          end do
        end if

	end

c*************************************************************
c*************************************************************
c*************************************************************

