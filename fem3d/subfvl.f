
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

c routine to write volume
c
c revision log :
c
c 28.04.2010    ggu     written from scratch
c
c******************************************************************

	subroutine wrfvla

c write of finite volume data

	use mod_layer_thickness
	use mod_area
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer k,l,lmax,id,nvar,ishyff

	real getpar
	logical has_output,next_output
	real saux(nlvdi,nkn)

	integer ia_out(4)
	save ia_out

        integer, save :: icall = 0

c start of code

        if( icall .eq. -1 ) return

c initialization

        if( icall .eq. 0 ) then

	  ishyff = nint(getpar('ishyff'))
	  call init_output('itmcon','idtcon',ia_out)
	  if( ishyff == 1 ) icall = -1
	  if( .not. has_output(ia_out) ) icall = -1
	  if( icall .le. -1 ) return

	  nvar = 1
	  call open_scalar_file(ia_out,nlv,nvar,'fvl')

        end if

c normal call

        icall = icall + 1

	if( .not. next_output(ia_out) ) return

	do k=1,nkn
	  lmax = ilhkv(k)
	  do l=1,lmax
	    saux(l,k) = areakv(l,k) * hdknv(l,k)
	  end do
	end do

        id = 66       			!for finite volume
	call write_scalar_file(ia_out,id,nlvdi,saux)

	end

c******************************************************************

