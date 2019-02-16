
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

	real cpa,cp,emiss,bolz,Kelvin,const06,pi,deg2rad,rad2deg

        parameter     (             cpa=1008.)
        parameter     (             cp=3985.)
        parameter     (             emiss=0.97)
        parameter     (             bolz=5.67e-8)
        parameter     (             Kelvin=273.16)
        parameter     (             const06=0.62198)
        parameter     (             pi=3.14159265358979323846)
        parameter     (             deg2rad=pi/180.)
        parameter     (             rad2deg=180./pi)

	real a1,a2,a3,a4,a5,a6,a7,eps

        parameter     (             a1=6.107799961)
        parameter     (             a2=4.436518521e-1)
        parameter     (             a3=1.428945805e-2)
        parameter     (             a4=2.650648471e-4)
        parameter     (             a5=3.031240396e-6)
        parameter     (             a6=2.034080948e-8)
        parameter     (             a7=6.136820929e-11)
        parameter     (             eps=1.0e-12)

