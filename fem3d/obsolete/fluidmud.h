
!--------------------------------------------------------------------------
!
!    Copyright (C) 2015,2019  Georg Umgiesser
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
! 19.01.2015	ggu	changed VERS_7_1_3
! 16.02.2019	ggu	changed VERS_7_5_60

        double precision rhosed		!Mud primary particle density (kg/m3)
        common /rhosed/ rhosed
	save /rhosed/

 	double precision dm0,nf
 	common /dm0/ dm0
	common /nf/ nf
	save /dm0/,/nf/

        real z0bkmud(nkndim)       !bottom roughenss on nodes for mud
        common /z0bkmud/z0bkmud
        save /z0bkmud/

        real mudc(nlvdim,nkndim)        !Fluid mud concentrationarray (kg/m3)
        common /mudc/mudc
        double precision rhomud(nlvdim,nkndim) !Mud floc part. density (kg/m3)
        common /rhomud/rhomud
        save /mudc/,/rhomud/

        real visv_yield(0:nlvdim,nkndim) !viscosity (mud)
        common /visv_yield/visv_yield
        real difv_yield(0:nlvdim,nkndim) !diffusivity (mud)
        common /difv_yield/difv_yield
	save /visv_yield/,/difv_yield/

        real lambda(nlvdim,nkndim)    ! Structural parameter
        common /lambda/lambda
        real vts(0:nlvdim,nkndim)        ! Rheological Viscosity [m2/s]
        common /vts/vts
        real dmf_mud(nlvdim,nkndim)  ! Floc size array.
        common /dmf_mud/dmf_mud
	save /lambda/,/vts/,/dmf_mud/

        real wprvs(0:nlvdim,nkndim)    !water density
        common /wprvs/wprvs
	save /wprvs/

