
!--------------------------------------------------------------------------
!
!    Copyright (C) 2010-2011,2014-2016,2019  Georg Umgiesser
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
c 28.04.2010	ggu	written from scratch
c 03.05.2010	ggu	changed VERS_6_1_8
c 14.04.2011	ggu	changed VERS_6_1_22
c 22.11.2011	ggu	changed VERS_6_1_37
c 26.11.2014	ggu	changed VERS_7_0_7
c 23.12.2014	ggu	changed VERS_7_0_11
c 19.01.2015	ggu	changed VERS_7_1_3
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 24.07.2015	ggu	changed VERS_7_1_82
c 18.09.2015	ggu	changed VERS_7_2_3
c 07.06.2016	ggu	changed VERS_7_5_12
c 16.02.2019	ggu	changed VERS_7_5_60
c 01.07.2019	ggu	not writing file anymore
c 20.03.2022	ggu	upgraded to da_out
c
c******************************************************************

	subroutine wrfvla

c write of finite volume data

	use mod_layer_thickness
	use mod_area
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer k,l,lmax,id,nvar,idc

	real getpar
	logical has_output_d,next_output_d
	real saux(nlvdi,nkn)
	double precision dtime

	double precision, save :: da_out(4)

        integer, save :: icall = 0

c start of code

        if( icall .eq. -1 ) return

c initialization

        if( icall .eq. 0 ) then

	  call init_output_d('itmcon','idtcon',da_out)
	  if( .not. has_output_d(da_out) ) icall = -1
	  if( icall .le. -1 ) return

	  nvar = 1
	  call shyfem_init_scalar_file('fvl',nvar,.false.,id)
	  da_out(4) = id

        end if

c normal call

        icall = icall + 1

	if( .not. next_output_d(da_out) ) return

	do k=1,nkn
	  lmax = ilhkv(k)
	  do l=1,lmax
	    saux(l,k) = areakv(l,k) * hdknv(l,k)
	  end do
	end do

	id = nint(da_out(4))
        idc = 66       			!for finite volume
	call get_act_dtime(dtime)
	call shy_write_scalar_record(id,dtime,idc,nlvdi,saux)

	end

c******************************************************************

