
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

c ecological dummy module
c
c revision log :
c
c 29.04.2008    ggu     bfm model integrated in main branch
c 18.02.2011    ggu     general framework for ecological model
c
c**************************************************************

        subroutine ecological_module

c general interface to ecological module

        implicit none

	integer ibfm,ibio
	real getpar

	integer, save :: icall = 0

	if( icall .eq. -1 ) return

        ibfm = nint(getpar('ibfm'))
	ibio = nint(getpar('ibio'))

	ibfm = 0
	if( ibfm .gt. 0 ) then
	  write(6,*) 'BFM module has not been linked'
	  write(6,*) 'ibfm = ',ibfm
	  write(6,*) 'You must enable this feature in Rules.make'
	  stop 'error stop ecological_module: ibfm'
	end if

	if( ibio .gt. 0 ) then
	  write(6,*) 'No ecological module has been linked'
	  write(6,*) 'ibio = ',ibio
	  write(6,*) 'You must enable this feature in Rules.make'
	  stop 'error stop ecological_module: ibio'
	end if

	icall = -1

        end

c**************************************************************

        subroutine write_restart_eco(iunit)
        implicit none
	integer iunit
	integer nstate,nkn,i
	nstate = 0
	nkn = 0
        write(iunit) nstate,nkn
	end
        subroutine skip_restart_eco(iunit)
        implicit none
	integer iunit
	integer nstate,nkn,i
        read(iunit) nstate,nkn
        do i=1,nstate
          read(iunit)
        end do
	end
        subroutine read_restart_eco(iunit)
        implicit none
	integer iunit
	call skip_restart_eco(iunit)
	end

c**************************************************************

