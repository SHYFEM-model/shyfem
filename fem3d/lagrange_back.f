
!--------------------------------------------------------------------------
!
!    Copyright (C) 2005-2006,2009,2012  Andrea Cucco
!    Copyright (C) 2006,2010,2012,2014-2015,2018-2019  Georg Umgiesser
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

c subroutines for lagrangian backtracing
c
c revision log :
c
c 00.00.2005	aac	from scratch
c 24.03.2006	aac	use flxtype in lagr_vel to avoid backtracing at OB
c 12.06.2006	ggu	integrated in new version
c 29.11.2006	ggu	new version with different backtrace possibilities
c 27.01.2009	aac	inertial effects from high vel to low vel
c 04.02.2009	aac	cleaned up
c 23.03.2010	ggu	changed v6.1.1
c 23.03.2012	aac	new calling for track_body (id)
c 30.03.2012	ggu	changed VERS_6_1_51
c 23.12.2014	ggu	changed VERS_7_0_11
c 19.01.2015	ggu	changed VERS_7_1_3
c 01.04.2015	ggu	changed VERS_7_1_7
c 21.05.2015	ggu	changed VERS_7_1_11
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 16.11.2015	ggu	changed VERS_7_3_14
c 25.10.2018	ggu	changed VERS_7_5_51
c 16.02.2019	ggu	changed VERS_7_5_60
c
c**********************************************************************

	subroutine back_trace(uadv,vadv)

c handles lagrangian backtracing

	use mod_lagrange
	use basin, only : nkn,nel,ngr,mbw

	implicit none

        include 'param.h'

                
	real uadv(1), vadv(1)
        
	bback = .true.
	nback = nel

	call lagr_setup_timestep	!sets up length of sides and fluxes

	call bar_ie     	!initializes particles at centers of elements
	
 	call drogue_back

        call lagr_vel(uadv,vadv)
        
	end

c****************************************************************
        
        subroutine bar_ie

c initializes particles at baricenter

	use mod_lagrange

        implicit none
        
        include 'param.h'

        integer i
        real xb,yb
        
        do i=1,nback
          call baric(i,xb,yb)
          x_back(i)=xb
          y_back(i)=yb
          ie_back(i)=i                   
        end do

        end 

c**********************************************************************

	subroutine drogue_back

	use mod_lagrange

	implicit none
	
        include 'param.h'

	integer i,ie,id,lb
	real x,y,z
	real dt,ttime

	call get_timestep(dt)

        do i=1,nback
          x  = x_back(i)
          y  = y_back(i)
          z  = lgr_ar(i)%actual%z
          lb = lgr_ar(i)%actual%l

          ie = ie_back(i)
	  id = lgr_ar(i)%id

	  ttime=dt
          call track_body(i,id,x,y,z,lb,ie,ttime)

          x_back(i) = x
          y_back(i) = y
          lgr_ar(i)%actual%z = z
          lgr_ar(i)%actual%l = lb

          ie_back(i) = ie
        end do

	end	

c**********************************************************************

        subroutine lagr_vel(uadv,vadv)

c interpolation of velocities on the points that have been backtraced
        
	use mod_lagrange
	use mod_depth
	use mod_hydro_print
	use mod_hydro
	use basin

        implicit none

        include 'param.h'
        
	real uadv(1), vadv(1)

                
       
        

       
        real x,y
        real vold,vnow
        real up,vp,zp
        real u(3),v(3),z(3)
        integer i,ii,ie,k
	real uc,vc,zc
	integer ibtype,ib
	logical blimit,binertial
        
        integer flxtype
        
	ibtype = 4		!type of backtracing
	blimit = .false.	!use depth limiter from Andrea
	binertial = .false.	!use inertial limiter from Andrea

        do ie=1,nel

c	 --------------------------------------
c	 first compute values in central point of actual element
c	 --------------------------------------

	 uc = 0.
	 vc = 0.
	 zc = 0.
         do ii=1,3
            k=nen3v(ii,ie)
            uc=uc+uprv(1,k)
            vc=vc+vprv(1,k)
	    zc=zc+zeov(ii,ie)
         end do
	 uc = uc / 3.
	 vc = vc / 3.
	 zc = zc / 3.
	 zc = zc + hev(ie)

c	 --------------------------------------
c	 no do the same for point of back-advection
c	 --------------------------------------

         uadv(ie) = 0
         vadv(ie) = 0

         x=x_back(ie)
         y=y_back(ie)
         i=ie_back(ie)

         if(i.gt.0) then		! only if we found element
           do ii=1,3
            k=nen3v(ii,i)
            u(ii)=uprv(1,k)
            v(ii)=vprv(1,k)
           end do
           call femintp(i,u,x,y,up)
           call femintp(i,v,x,y,vp)
           u_lag(ie)=up
           v_lag(ie)=vp

           call femintp(i,zeov(1,i),x,y,zp)
	   zp = zp + hev(i)
 	   if( blimit .and. zp .ge. zc ) zp = zc

c	   --------------------------------------
c	   in u/vadv we now store the difference of transport
c	   --------------------------------------

	   ib = ibtype

           if( binertial ) then		!inertial effects from high to low vel
	     vold = sqrt((utlov(1,ie)**2 +vtlov(1,ie)**2))
             vnow = sqrt((u_lag(ie)**2 +v_lag(ie)**2))
             if( vold .lt. vnow ) ib = 0
	   end if

           if( ib .eq. 0 ) then
             uadv(ie) = 0
             vadv(ie) = 0
           else if( ib .eq. 1 ) then
             uadv(ie) = uc*zc - u_lag(ie) * zp
             vadv(ie) = vc*zc - v_lag(ie) * zp
           else if( ib .eq. 2 ) then
             uadv(ie) = ( uc - u_lag(ie) ) * zp
             vadv(ie) = ( vc - v_lag(ie) ) * zp
           else if( ib .eq. 3 ) then
             uadv(ie) = ( uc - u_lag(ie) ) * zc
             vadv(ie) = ( vc - v_lag(ie) ) * zc
           else if( ib .eq. 4 ) then               !original version
             uadv(ie) = utlov(1,ie) - u_lag(ie) * zp
             vadv(ie) = vtlov(1,ie) - v_lag(ie) * zp
           else
             stop 'error stop lagr_vel: ibtype'
           end if

         end if

c        check if we are on open boundary -> no backtracking

         do ii=1,3
          k=nen3v(ii,ie)                 
          if( flxtype(k) .ge. 3 ) then  ! we are close to open boundary
           uadv(ie) = 0
           vadv(ie) = 0
          endif
         end do

        end do
       
        end 

c**********************************************************************

