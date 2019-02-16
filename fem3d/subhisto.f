
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

c routines dealing with histogram
c
c contents :
c
c subroutine histo_init(nbin,bin0,dbin,rbin)
c subroutine histo_insert(value)
c subroutine histo_final(ic)
c 
c revision log :
c
c 16.03.2004    ggu     routines written from scratch
c
c****************************************************************

        subroutine histo_init(nbin,bin0,dbin,rbin)

        implicit none

        integer nbin            !total number of bins
        real bin0               !first bin (limit)
        real dbin               !regular bin size (0 => use rbin)
        real rbin(nbin)         !bin size limits (upper)

	include 'histo.h'

        integer i

        if( nbin .gt. ndim_histo ) stop 'error stop histo_init: nbin'

        if( dbin .gt. 0 ) then
          do i=1,nbin
            abin(i) = bin0 + (i-1) * dbin
          end do
        else
          do i=1,nbin
            abin(i) = rbin(i)
          end do
        end if

        ncbin = nbin
        do i=1,nbin+1
          icount(i) = 0
        end do

        end

c****************************************************************

        subroutine histo_insert(value)

        implicit none

        real value

	include 'histo.h'

        integer i

        do i=1,ncbin
          if( value .le. abin(i) ) then
            icount(i) = icount(i) + 1
            return
          end if
        end do

        i = ncbin+1
        icount(i) = icount(i) + 1

        end

c****************************************************************

        subroutine histo_final(ic)

        implicit none

        integer ic(1)

	include 'histo.h'

        integer i

        do i=1,ncbin+1
          ic(i) = icount(i)
        end do

        end

c****************************************************************

