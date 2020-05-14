
!--------------------------------------------------------------------------
!
!    Copyright (C) 2017,2019  Georg Umgiesser
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

! revision log :
!
! 05.12.2017	ggu	changed VERS_7_5_39
! 14.02.2019	ggu	changed VERS_7_5_56

	program whdf

c writes hdf files

	parameter (irank=3)
	parameter (idimx=25,idimy=10,idimz=50)

	real u(idimx,idimy,idimz)
	real v(idimx,idimy,idimz)
	real w(idimx,idimy,idimz)
	real s(idimx,idimy+1,idimz)

	real a(idimx,idimy)
	real b(idimx,idimy+1)

	integer idim(irank)
	integer iaux(irank)
	data idim /idimx,idimy,0/
	data iaux /idimx,idimy,0/

	lu=0
	lv=0
	lw=0
	ls=0

    1	continue

	write(6,*) 'This program reads from stdin'

	read(5,end=99) it,ir,iwhat
	read(5) (idim(i),i=1,ir)
	write(6,*) it,ir,iwhat,(idim(i),i=1,ir)

	if(iwhat.eq.1) then

	  if(idim(1).ne.idimx.or.idim(2).ne.idimy) stop 'error dim'
	  read(5) a
	  lu=lu+1
	  if(lu.gt.idimz) stop 'error dim z'
	  do iy=1,idimy
	    do ix=1,idimx
	      u(ix,iy,lu) = a(ix,iy)
	    end do
	  end do

        else if(iwhat.eq.2) then

          if(idim(1).ne.idimx.or.idim(2).ne.idimy) stop 'error dim'
          read(5) a
          lv=lv+1
	  if(lv.gt.idimz) stop 'error dim z'
          do iy=1,idimy
            do ix=1,idimx
              v(ix,iy,lv) = a(ix,iy)
            end do
          end do

	else if(iwhat.eq.3) then

	  if(idim(1).ne.idimx.or.idim(2).ne.idimy) stop 'error dim'
	  read(5) a
	  lw=lw+1
	  if(lw.gt.idimz) stop 'error dim z'
	  do iy=1,idimy
	    do ix=1,idimx
	      w(ix,iy,lw) = a(ix,iy)
	    end do
	  end do
	  
	else if(iwhat.eq.5) then

	  if(idim(1).ne.idimx.or.idim(2).ne.idimy+1) stop 'error dim'
	  read(5) b
	  ls=ls+1
	  if(ls.gt.idimz) stop 'error dim z'
	  do iy=1,idimy+1
	    do ix=1,idimx
	      s(ix,iy,ls) = b(ix,iy)
	    end do
	  end do

	end if

	goto 1

   99	continue

	write(6,*) lu,lv,lw,ls

	iaux(3)=lu
	iret = dsadata('u.hdf',irank,iaux,u)

	iaux(3)=lv
	iret = dsadata('v.hdf',irank,iaux,v)

	iaux(3)=lw
	iret = dsadata('w.hdf',irank,iaux,w)

	iaux(2)=idimy+1
	iaux(3)=ls
	iret = dsadata('s.hdf',irank,iaux,s)

	stop
	end
