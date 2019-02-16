
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

c**********************************************************

	program offinf

c shows content of offline data file

	implicit none

	integer it,nkn,nel,nrec,iu,i,type,irecs
	integer ios
	character*60 name
	!double precision buffer(5)
	real buffer(1)

        call off_init(name)

	nrec = 0
	iu = 1
	irecs = 9

	open(iu,file=name,status='old',form='unformatted',iostat=ios)
	if( ios /= 0 ) stop 'error stop offinf: opening file'

	do
          read(iu,iostat=ios) it,nkn,nel,type
	  if( ios /= 0 ) exit
	  nrec = nrec + 1
	  write(6,*) nrec,it,nkn,nel,type
	  if( type /= 3 ) stop 'error stop offinf: type /= 3'

	  do i=1,irecs
	    read(iu) buffer
	  end do

	  !write(6,'(5g14.6)') buffer
	end do

	if( ios > 0 ) stop 'error stop offinf: reading file'

	close(iu)

	end

c**********************************************************

        subroutine off_init(offfile)

        use clo

        implicit none

        character*(*) offfile

        call shyfem_copyright('offinf - info on offline file')

        call clo_init('offinf','offfile','1.0')

	call clo_add_info('returns info on records of offline file')

        call clo_parse_options

        call clo_check_files(1)
        call clo_get_file(1,offfile)

        end

c**********************************************************

