
!--------------------------------------------------------------------------
!
!    Copyright (C) 2018-2019  Georg Umgiesser
!    Copyright (C) 2018  Debora Bellafiore
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

! projection handling routines
!
! contents :
!
! subroutine shy_proj
!
! revision log :
!
! 20.07.2018	dbf	projection routines for shyelab
! 30.08.2018	ggu	bug fix if mode == 0
! 16.02.2019	ggu	changed VERS_7_5_60
!
!*****************************************************************

        subroutine shy_proj

        use clo
        use basin
        use shympi
        use projection
        use coordinates

        implicit none

        character*80 sproj
        integer iscand,jj,iproj1
        integer                         :: mode,i
        double precision		:: c_param(9)
        double precision		:: d(9)
        real                            :: xgeov1(nkn),ygeov1(nkn) 

	logical is_geographical
        
        !mode: +1: cart to geo  -1: geo to cart                
        !iproj1 1=GB 2=UTM 3=CPP 4=UTM non standard
        
        call clo_get_option('proj',sproj)

        jj=iscand(sproj,d,6)		!mode,iproj,c_param
        if( jj == 0 ) return

        mode = nint(d(1))
        if( mode == 0 ) return
        iproj1 = nint(d(2))
        c_param(1)=d(3)
        c_param(2)=d(4)
        c_param(3)=d(5)
        if( jj>=6 ) c_param(4)=d(6)

        !write(6,*) jj,mode,iproj1,c_param

        call init_coords(iproj1,c_param)

        if( mode == 1) then
                call convert_coords(1,nkn,xgv,ygv,xgeov1,ygeov1)
        else
		if( .not. is_geographical(nkn,xgv,ygv) ) goto 99
                call convert_coords(-1,nkn,xgv,ygv,xgeov1,ygeov1)
        endif

        xgv = xgeov1
        ygv = ygeov1

	return
   99	continue
        write(6,*)'coordinates are not lat,lon'
        write(6,*)'  xmin,xmax: ',minval(xgv),maxval(xgv)
        write(6,*)'  ymin,ymax: ',minval(ygv),maxval(ygv)
        stop' error stop shy_proj: not lat/lon'
        end subroutine shy_proj
       
!*****************************************************************

