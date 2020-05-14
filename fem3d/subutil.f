
!--------------------------------------------------------------------------
!
!    Copyright (C) 2020  Georg Umgiesser
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

! general utility routines (no other dependencies)
!
! contents :
!
! subroutine convert_uv_sd(u,v,s,d,bmeteo)	convert from uv to sd
! subroutine convert_sd_uv(s,d,u,v,bmeteo)	convert from sd to uv
!
! revision log :
!
! 17.04.2020	ggu	directional conversion routines copied here
!

!*********************************************************************

        subroutine convert_uv_sd(u,v,s,d,bmeteo)

! convert from uv to sd

        implicit none

        real u,v,s,d
        logical bmeteo          !use meteo convention

        real dir
        real pi,rad,rrad
        parameter(pi=3.14159,rad=pi/180.,rrad=1./rad)

        s = sqrt(u*u+v*v)
        d = 0.
        if( s > 0 ) d = atan2(v/s,u/s)

        d = 90. - d * rrad
        if( bmeteo ) d = d + 180.
        if( d < 0. ) d = d + 360.
        if( d > 360. ) d = d - 360.

        end

!*********************************************************************

        subroutine convert_sd_uv(s,d,u,v,bmeteo)

! convert from sd to uv

        implicit none

        real s,d,u,v
        logical bmeteo          !use meteo convention

        real dir
        !real pi,rad
        !parameter(pi=3.14159,rad=pi/180.)
	real, parameter :: pi = 4.*atan(1.)
	real, parameter :: rad = pi/180.

        dir = 90. - d
        if( bmeteo ) dir = dir + 180.
        do while( dir .lt. 0. )
          dir = dir + 360.
        end do
        dir = mod(dir,360.)

        u = s*cos(rad*dir)
        v = s*sin(rad*dir)

        end

!*********************************************************************

