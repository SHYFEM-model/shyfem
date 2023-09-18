!
! $Id: lagrange_back.f,v 1.2 2009-03-11 16:25:59 georg Exp $
!
! subroutines for lagrangian backtracing
!
! revision log :
!
! 00.00.2005    aac	from scratch
! 24.03.2006    aac	use flxtype in lagr_vel to avoid backtracing at OB
! 12.06.2006    ggu	integrated in new version
! 29.11.2006    ggu	new version with different backtrace possibilities
! 27.01.2009    aac     inertial effects from high vel to low vel
! 04.02.2009    aac     cleaned up
! 23.03.2012    aac     new calling for track_body (id)
!
!**********************************************************************
!-----------------------------------------------------------------------
        module lagrange_back
!-----------------------------------------------------------------------
        contains
!-----------------------------------------------------------------------

	subroutine back_trace(uadv,vadv)

! handles lagrangian backtracing

	use lagrange_data
        use lagrange_flux
	use basin, only : nkn,nel,ngr,mbw

	implicit none

        include 'param.h'

                
	double precision uadv(1), vadv(1)
        
	bback = .true.
	nback = nel

	call lagr_setup_timestep	!sets up length of sides and fluxes

	call bar_ie     	!initializes particles at centers of elements
	
 	call drogue_back

        call lagr_vel(uadv,vadv)
        
	end

!****************************************************************
        
        subroutine bar_ie

! initializes particles at baricenter

	use lagrange_data
        use fem_util

        implicit none
        
        include 'param.h'

        integer i
        double precision xb,yb
        
        do i=1,nback
          call baric(i,xb,yb)
          x_back(i)=xb
          y_back(i)=yb
          ie_back(i)=i                   
        end do

        end 

!**********************************************************************

	subroutine drogue_back

	use lagrange_data
	use lagrange_main
        use time_util

	implicit none
	
        include 'param.h'

	integer i,ie,id,lb
	double precision x,y,z
	double precision dt,ttime

	call get_timestep(dt)

        do i=1,nback
          x  = x_back(i)
          y  = y_back(i)
          z  = lgr_ar(i)%z
          lb = lgr_ar(i)%l

          ie = ie_back(i)
	  id = lgr_ar(i)%id

	  ttime=dt
          call track_body(i,id,x,y,z,lb,ie,ttime)

          x_back(i) = x
          y_back(i) = y
          lgr_ar(i)%z = z
          lgr_ar(i)%l = lb

          ie_back(i) = ie
        end do

	end	

!**********************************************************************

        subroutine lagr_vel(uadv,vadv)

! interpolation of velocities on the points that have been backtraced
        
	use lagrange_data
	use depth
	use hydro_print
	use hydro_admin
	use basin
        use regular
        use discharge_flux

        implicit none

        include 'param.h'
        
	double precision uadv(1), vadv(1)
       
        double precision x,y
        double precision vold,vnow
        double precision up,vp,zp
        double precision u(3),v(3),z(3)
        integer i,ii,ie,k
	double precision uc,vc,zc
	integer ibtype,ib
	logical blimit,binertial
        
	ibtype = 4		!type of backtracing
	blimit = .false.	!use depth limiter from Andrea
	binertial = .false.	!use inertial limiter from Andrea

        do ie=1,nel

!	 --------------------------------------
!	 first compute values in central point of actual element
!	 --------------------------------------

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

!	 --------------------------------------
!	 no do the same for point of back-advection
!	 --------------------------------------

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

!	   --------------------------------------
!	   in u/vadv we now store the difference of transport
!	   --------------------------------------

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

!        check if we are on open boundary -> no backtracking

         do ii=1,3
          k=nen3v(ii,ie)                 
          if( flxtype(k) .ge. 3 ) then  ! we are close to open boundary
           uadv(ie) = 0
           vadv(ie) = 0
          endif
         end do

        end do
       
        end 

!**********************************************************************

!-----------------------------------------------------------------------
        end module lagrange_back
!-----------------------------------------------------------------------
