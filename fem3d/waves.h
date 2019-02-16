
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

	integer iwave			!call for wave model
	integer iwwm			!call for coupling with wwm
	integer idcoup			!time step for syncronizing with wwm [s]
	common /wavcst/ iwave,iwwm,idcoup

        real waveh(nkndim)      	!sign. wave height [m]
        common /waveh/ waveh

        real wavep(nkndim)      	!mean wave period [s]
	common /wavep/ wavep

        real wavepp(nkndim)      	!peak wave period [s]
	common /wavepp/ wavepp

        real waved(nkndim)      	!mean wave direction
	common /waved/ waved

        real waveov(nkndim)     	!wave orbital velocity
	common /waveov/ waveov

        real wavefx(nlvdim,neldim)      !wave forcing terms
	common /wavefx/ wavefx

        real wavefy(nlvdim,neldim)
	common /wavefy/ wavefy

	save /wavcst/
	save /waveh/,/wavep/,/wavepp/,/waved/
	save /waveov/,/wavefx/,/wavefy/

