
!--------------------------------------------------------------------------
!
!    Copyright (C) 2020-2021  Micol Pucci
!    Copyright (C) 2021  Georg Umgiesser
!    Copyright (C) 2021  Debora Bellafiore
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

! notes :
!
! for 2D applications rfricv(ie) is changed
! for 3D applications ifricv(l,ie) is changed
!
! Pucci et al. https://doi.org/10.3390/jmse8121010
!
! revision log :
!
! 01.01.2020    mpc     started routines from scratch
! 01.07.2021    mpc     transfered routines from subcus.f
! 17.07.2021    ggu&dbf adapted to SHYFEM framework
!
!*****************************************************************

        subroutine turbine 

c	This routine allows to reproduce the presence of a tidal turbine (vertical axis) in 2D simulations.
c	In the first part, there is the process of identifying the grid elements belonging to the turbine.
C	Then, by using flow local condition, there is the evaluation of Lift and Drag coefficient (CL and CD) in static or dynamic conditions 
c	using the  SUBROUTINE FOR STATIC COEFFICIENT or SUBROUTINE FOR DYNAMIC COEFFICIENT respectively.
c 	The external text file called 'dati.txt' contains CL and CD tabulated as a function of the Reynolds number and the angle of attack.
c	The 'dati.txt' file must contain number of rows, numers of columns, Reynolds number values, angle of attack values, CL and CD value in the mentioned order.	
c	Pucci et al. https://doi.org/10.3390/jmse8121010
c	Before starting the routine check ALL the assigned parameters       

        use mod_hydro_vel
        use basin
        use mod_diff_visc_fric
        use evgeom
        
        use mod_nudging
        use mod_internal
        use mod_geom_dynamic
        use mod_depth
        use mod_bound_dynamic
        use mod_hydro_baro
        use mod_hydro
        use levels
      
        implicit none

c-----------------------------------------
c icall = 0	run the turbine model
c icall = -1	do not run the turbine model
c-----------------------------------------					

        integer, save :: icall = -1     ! number of time step

c-----------------------------------------
c turbine characteristics and fluid parameters
c-----------------------------------------					
        real x0, y0                     ! coordinates of the turbine centre [m]
        parameter(x0=72., y0=75.)
        real r, chord, nb, om           ! turbine characteristc radius [m], chord [m], number of blades, rotational speed [ras/s] 
        parameter(r=3., chord=0.4, nb=3.,om=1.575)
        real pi, ro, mu                 ! parameters, ro and mu are fluid density [kg/m^3] and viscosity [Pa*s] respectively 
        parameter(pi=3.141592654,  ro=998.2, mu=0.001003)
        integer segno_rot, n_ring, static_dyn_flag
        parameter(segno_rot=+1)         ! enter +1 for counterclockwise rotation or -1 for clockwise rotation
        parameter(n_ring=300)           ! enter the number of ring elements
        parameter(static_dyn_flag=0)    ! enter 0 for dynamic condition or 1 for static condition
        real U_inf                      ! undisturbed flow velocity [m/s]
        parameter(U_inf=1.75)
        real ring_thick                 ! thickness of the ring representing the 2D turbine [m]
        parameter(ring_thick=0.4)
c-----------------------------------------
c file 'dati.txt' variables
        real, save :: h(30),e(5),B(5,30),G(5,30)! file 'dati.txt' data, respectively angle of attack, Reynolds number, CL and CD, change array dimentions if needed 
        integer, save :: ic, ir         ! number of columns and rows in the 'dati.txt' file
c-----------------------------------------
        real x,y     ! coordinates of the element barycentre
        integer, save :: icount,disc_count   ! to count and check the elements found on the actuator ring and on the internal disc respectively
        integer, save, allocatable :: ring_ie(:)        ! ID of the ring elements belonging to the turbine 
        integer, save, allocatable :: disc_ie(:)        ! ID of the internal disc elements belonging to the turbine
        real, save, allocatable :: ring_theta(:)        ! azimuthal angle of the turbine ring elements 
        real, save, allocatable ::disc_theta(:)         ! azimuthal angle of the turbine disc elements 
        real theta,thetad   ! local value for the azimuthal angle in degrees of an element on the actuator ring and on the disc respectively 
        integer app_ie ! temporary value assigned to the element ID
        real app_ang   ! temporary value assigned to the theta angle
        real dist, aux3, aux4,aux5,aux6 ! dist is the distance between the element barycentre and the turbine centre
        integer, save :: ie0_75Dcount,ie1Dcount ! if necessary number of grid elements in wake at a distance of 1-3-5 diameters from the turbine (exception for 0.75D which is in front of the turbine)
        integer, save :: ie3Dcount,ie5Dcount
        integer, save :: ie_0_75D(1000),ie_1D(1000) ! ID of wake elements 
        integer, save :: ie_3D(1000),ie_5D(1000)
        real x1_0_75,x2_0_75,x1_1,x2_1
        real x1_3,x2_3, x1_5,x2_5
        parameter(x1_0_75=67.,x2_0_75=68.,x1_1=77.,x2_1=79.) ! enter the correct extreme to identify wake elements 
        parameter(x1_3=89.,x2_3=90., x1_5=101.,x2_5=103.)
        real yaux1,yaux2, deltay ! to virtually devide the turbine in horizontal stripes of thickness deltay [m]
        integer n_stripes
        parameter(n_stripes=18, deltay=0.36)    ! deltay=(2r+ring_thick)/n_stripes
        !DEB parameter(yaux1=71.8)                   ! yaux1=y0-r-ring_thick/2
        real, save :: theta_los(18)             ! azimuthal angle relative to each stripe
        real dist_2,dist_3
        real W,W_n,W_t,W_n1,W_n2,W_t1,W_t2      ! relative velocity to the blade [m/s]
        integer ie_prova
        real  f_Re, f_alpha, Cl_st2,Cd_st2      ! static subroutine variables 
        integer l,ie, i 
        integer j,el,los
        real qc,corda                           ! corda is the geometric cord of a circumference relative to an element with y distance abs(y-y0) from the centre 
        parameter(qc=1.2732)                    ! qc is the ratio between the area of a square with edge 2r and a disc with radius r 
        real Li, D, Re                          ! Lift and Drag force both [N/m^3] and Reynolds number
        real frict,frict_1,frict_2,frict_3! friction and auxiliary variables
        real area_ele,idt, foo, it, uv,u_abs   ! foo is the lift coefficient in dynamic condition
        real, allocatable :: alpha(:)  
        real al_min, al_max,al_in, al_in0,al_ds! dynamic subroutine variables
        real cl_fin(4001),Cl_st, Cd_st          ! dynamic subroutine variables
        real alpha_aux,area_frict
        real area_aux !DEB
           
c-----------------------------------------
c do we have to run turbine?
c
c if you wanto to run turbine, please set icall=0 above
c-----------------------------------------      

        if( icall < 0 ) return

c-----------------------------------------
c first time step
c-----------------------------------------      

         allocate (alpha(0:n_ring))    
        
        if(icall.eq.0) then

        allocate (ring_ie(n_ring)) 
        allocate (ring_theta(n_ring)) 
        allocate (disc_ie(n_ring))
        allocate (disc_theta(n_ring))
c        allocate (ie_1D(ie1Dcount)) 
c        allocate (ie_5D(ie5Dcount))
c-----------------------------------------
c read file 'dati.txt'
c-----------------------------------------        
        open(unit=55, file='dati.txt', status='old')  
        read(55,*) ir   
        read(55,*) ic   
        !allocate (e(ir))
        !allocate (h(ic))    
        !allocate (B(ir,ic))
        !allocate (G(ir,ic))
        read(55,*) (e(i), i=1, ir) !read Reynolds values
        read(55,*) (h(i), i=1, ic) !read alpha values
       
        do i=1, ir
          read(55,*) (B(i,j), j=1,ic)  !read Cl values
        end do
        do i=1, ir
          read(55,*) (G(i,j), j=1,ic)  !read Cd values
        end do
       
        close(55)
c-----------------------------------------
c turbine elements identification
c-----------------------------------------
        open(unit=82,file='ID_Celle_disco.txt',status='unknown')  
        icount=0
        ie1Dcount=0
        ie0_75Dcount=0
        disc_count=0
        
        do 10 ie=1,nel       
        call baric(ie,x,y)  
        dist = sqrt((x-x0)**2 + (y-y0)**2) 
        aux3=r+ring_thick/2. 
        aux4=r-ring_thick/2. 
       if(dist.lt.aux3.and.dist.gt.aux4)then  
          icount=icount+1               !a new element belonging to the actuator ring has been found
          ring_ie(icount)=ie
          theta=(asin((x0-x)/dist))*180/pi  !the function returns angles between -90 and +90 degrees
         if((segno_rot*(y-y0)/dist).lt.0)  then   !to achieve angles in the second half of upwind and in the first half of downwind
         theta=180-theta
         goto 99
         endif
         aux6=segno_rot*(y-y0)/dist
         if(aux6.ge.0.and.(asin((x0-x)/dist)*180/pi).le.0)then !to achieve angles in the second half of downwind
         theta=theta+360
         endif
        endif   
99      continue        
         ring_theta(icount)=theta
          if(dist.lt.aux4)then
          disc_count=disc_count+1       !a new element belonging to the internal disc has been found
          disc_ie(disc_count)=ie
         thetad=(asin((x0-x)/dist))*180/pi
         if((segno_rot*(y-y0)/dist).lt.0)  then   
        thetad=180-thetad
        goto 199
        endif
        aux6=segno_rot*(y-y0)/dist
        if(aux6.ge.0.and.(asin((x0-x)/dist)*180/pi).le.0)then 
        thetad=thetad+360
        endif
  199  continue  
        disc_theta(disc_count)=thetad  
          endif
10      continue
c-----------------------------------------
c if necessary to find the wake elements (exception for 0.75D which is in front of the turbine)
c-----------------------------------------
        do ie=1,nel
        call baric(ie,x,y)   
        open(unit=78,file='ID_Celle_1D.txt',status='unknown')
        open(unit=79,file='ID_Celle_075D.txt',status='unknown')
        open(unit=80,file='ID_Celle_3D.txt',status='unknown')
        open(unit=81,file='ID_Celle_5D.txt',status='unknown')
        if(x.gt.x1_1.and.x.lt.x2_1)then    !change extreme on the basis of the grid 
        ie1Dcount=ie1Dcount+1
        ie_1D(ie1Dcount)=ie
        write(78,*)ie_1D(ie1Dcount),ipev(ie_1D(ie1Dcount)),x,y
        end if
        if(x.gt.x1_0_75.and.x.lt.x2_0_75)then  !change extreme on the basis of the grid 
        ie0_75Dcount=ie0_75Dcount+1
        ie_0_75D(ie0_75Dcount)=ie
        write(79,*)ie_0_75D(ie0_75Dcount),
     +   ipev(ie_0_75D(ie0_75Dcount)),x,y
        end if
        if(x.gt.x1_3.and.x.lt.x2_3)then    !change extreme on the basis of the grid
        ie3Dcount=ie3Dcount+1
        ie_3D(ie3Dcount)=ie
        write(80,*)ie_3D(ie3Dcount),ipev(ie_3D(ie3Dcount)),x,y
        end if
        if(x.gt.x1_5.and.x.lt.x2_5)then    !change extreme on the basis of the grid
        ie5Dcount=ie5Dcount+1
        ie_5D(ie5Dcount)=ie
        write(81,*)ie_5D(ie5Dcount),ipev(ie_5D(ie5Dcount)),x,y
        end if
        
        end do

      
        if(icount.ne.n_ring) then
        write(6,*) 'error: element count is not the expected one'
        write(6,*) icount,n_ring
        endif

        open(unit=32,file='ID_celle_anello.txt',form='formatted',
     +  status='unknown')
        do 20 i=1,icount-1   !these two do-cycles ordinate both the arrays in ascending order of theta
           do  j=i+1,icount 
             if (ring_theta(i).gt.ring_theta(j)) then
                app_ang=ring_theta(i)
                app_ie=ring_ie(i)
                ring_theta(i)=ring_theta(j)
                ring_ie(i)=ring_ie(j)
                ring_theta(j)=app_ang
                ring_ie(j)=app_ie
             endif
             end do
        write(32,*) ring_ie(i),ipev(ring_ie(i)),ring_theta(i) 
20      continue
        close(32)

!the two arrays are now ordinated in ascending order of theta moreover their values have been memorized to be available for all the successive time steps
        open(unit=33, file='theta.txt', status='unknown')
        do i=1, icount 
        call baric(ring_ie(i),x,y)
        write(33,*) ring_theta(i), x, y
        end do  
        
        close(33)
c-----------------------------------------
c here the turbine is vitually devided in 18 horizontal stripes (change the number of stripes on the basis of n_ring) 
c to assign to the internal disc upwind elements (only those with theta<180 degrees) the same theta angle of the ring element belonging to the same stripe       
c-----------------------------------------      
        yaux1=71.8-0.36 !DEB
        yaux2=71.8 !DEB 
        !yaux1=yaux1-deltay
        !yaux2=yaux1
        do i=1,n_stripes
        yaux1=yaux1+0.36 !DEB
        yaux2=yaux2+0.36 !DEB
        !yaux1=yaux1+deltay
        !yaux2=yaux2+deltay
        area_aux=0
        do j=1,icount
        call baric(ring_ie(j),x,y)
            if(y.gt.yaux1.and.y.le.yaux2)then
            if(ring_theta(j).le.180)then
              theta_los(i)=ring_theta(j)
            end if
            end if
           
         end do 
        enddo !DEB        
    
         write(6,*)'theta_los',theta_los     
            
        end if   !if reltive to first time step icall=0

        icall=icall+1
        
        if(icall.gt.4000) then  !the fisrt 4000 time steps are used to converge the flow
        al_max=0
        al_min=0
        area_frict=0
c-----------------------------------------
c calculation of maximum and minimum alpha values (necessary for the dynamic subroutine)
c-----------------------------------------        
        do i=1,icount      
        l=1 
        W_t1=ulnv(l,ring_ie(i))*cos(ring_theta(i)*pi/180) 
        W_t2=(om*r+vlnv(l,ring_ie(i))*sin(ring_theta(i)*pi/180))        
        W_t=W_t1+segno_rot*W_t2                                         
        W_n1=ulnv(l,ring_ie(i))*sin(ring_theta(i)*pi/180) 
        W_n2=-segno_rot*vlnv(l,ring_ie(i))*cos(ring_theta(i)*pi/180)    
        W_n=W_n1+W_n2                                                   
        W=sqrt(W_t**2+W_n**2)
        alpha(i)=asin(W_n/W)*180/pi
        area_ele=12.*ev(10,ring_ie(i))
        area_frict=area_frict+area_ele
          if(alpha(i).gt.al_max) then       
          al_max=alpha(i)     
          endif
      
          if(alpha(i).lt.al_min) then         
          al_min=alpha(i)
          endif  
        end do
      
        if(al_max.gt.89.or.al_min.lt.-89) then
           write(6,*)'ERROR: alpha_max or alpha_min not acceptable'
           write(6,*) al_max, al_min 
        end if
c-----------------------------------------
c here starts the friction calculation for each actuator ring element
c-----------------------------------------     
        alpha(0)=0
        foo=0
        
        do i=1, icount      
        call baric(ring_ie(i),x,y)
        dist_2=abs(y0-y)
        l=1 !DEB
        W_t1=ulnv(l,ring_ie(i))*cos(ring_theta(i)*pi/180) 
        W_t2=(om*r+(vlnv(l,ring_ie(i)))*sin(ring_theta(i)*pi/180))      
        W_t=W_t1+segno_rot*W_t2
        W_n1=ulnv(l,ring_ie(i))*sin(ring_theta(i)*pi/180) 
        W_n2=-segno_rot*vlnv(l,ring_ie(i))*cos(ring_theta(i)*pi/180)    
        W_n=W_n1+W_n2 
        W=sqrt(W_t**2+W_n**2)
        alpha(i)=asin(W_n/W)*180/pi 
        if(abs(alpha(i)).lt.0.1)then 
        alpha(i)=0
        end if
        write(6,*)W_t,W_n,W,alpha(i)
        write(6,*)'u,v',ulnv(l,ring_ie(i)),vlnv(l,ring_ie(i))
        write(6,*)'theta',ring_theta(i)
        al_in=alpha(i)   
        al_in0=alpha(i-1)   
        alpha_aux=alpha(i)
        Re=ro*om*r*chord/mu 
        area_ele=12.*ev(10,ring_ie(i)) 
        u_abs=sqrt(ulnv(l,ring_ie(i))**2+vlnv(l,ring_ie(i))**2)
        Cd_st2=0
        uv=sqrt(utlnv(l,ring_ie(i))**2+vtlnv(l,ring_ie(i))**2)  

          call static(alpha(i), Re, Cl_st2, Cd_st2, ir, ic, h, e, B, G)

          if(abs(ring_theta(i)).lt.180)then
          call cl_dynamic(al_min, al_max,W,Re,al_in,al_in0,
     +  chord,om,Cl_st,ir,ic,h,e,B,G,
     +  foo,icall,alpha_aux,al_ds)
          write(6,*)'al_ds',al_ds,al_in
          write(6,*)'foo,ring_theta(i)',foo,ring_theta(i)
           Li=(0.5*foo*ro*chord*W**2)/area_ele     
           D=(0.5*Cd_st2*ro*chord*W**2)/area_ele
          else 
           Li=(0.5*Cl_st2*ro*chord*W**2)/area_ele
           D=(0.5*Cd_st2*ro*chord*W**2)/area_ele
          endif 

         frict_1=(Li**2+D**2)
         frict_3=(sqrt(frict_1))*6.4*nb/(qc*ro*uv) !6.4 is the "depth" change if necessary
         frict=frict_3/(disc_count+icount)
         if(abs(y-y0).gt.2.9.and.abs(y-y0).lt.3.2)then
         rfricv(ring_ie(i))=frict
         else
         corda=2*sqrt(r**2-dist_2**2)
         rfricv(ring_ie(i))=frict*2*r/corda
         end if

         
         open(unit=31, file='L_theta_partial.txt', status='unknown')    

        write(31,*) icall,ring_theta(i),alpha(i),foo,Cl_st2,Cd_st2,
     +  Li,D,W,area_ele,u_abs,rfricv(ring_ie(i)),utlnv(l,ring_ie(i)),
     +  vtlnv(l,ring_ie(i)),uv  

         end do !end of ring cycle
c-----------------------------------------
c here starts the friction calculation for each internal disc element
c-----------------------------------------         
         yaux1=yaux1-deltay
        yaux2=yaux1
        !do i=1,n_stripes
        do j=1,n_stripes !DEB
        yaux1=yaux1+deltay
        yaux2=yaux2+deltay
         do i=1,disc_count
         call baric(disc_ie(i),x,y)
         if(y.gt.yaux1.and.y.le.yaux2)then
         if(x.le.x0)then
         disc_theta(i)=theta_los(j)
         end if
        dist_2=abs(y0-y)
        l=1 
        W_t1=ulnv(l,disc_ie(i))*cos(disc_theta(i)*pi/180) 
        W_t2=(om*r+(vlnv(l,disc_ie(i)))*sin(disc_theta(i)*pi/180))      
        W_t=W_t1+segno_rot*W_t2                                         
        W_n1=ulnv(l,disc_ie(i))*sin(disc_theta(i)*pi/180) 
        W_n2=-segno_rot*vlnv(l,disc_ie(i))*cos(disc_theta(i)*pi/180)    
        W_n=W_n1+W_n2                                                   
        W=sqrt(W_t**2+W_n**2)
        alpha(i)=asin(W_n/W)*180/pi 
        write(6,*)'W_t,W_n,W,alpha(i)',W_t,W_n,W,alpha(i)
        if(abs(alpha(i)).lt.0.1)then 
        alpha(i)=0
        end if
        al_in=alpha(i)   
        al_in0=alpha(i-1)  
        alpha_aux=alpha(i)
        Re=ro*om*r*chord/mu 
        area_ele=12.*ev(10,disc_ie(i)) 
        u_abs=sqrt(ulnv(l,disc_ie(i))**2+vlnv(l,disc_ie(i))**2)
        Cd_st2=0
        uv=sqrt(utlnv(l,disc_ie(i))**2+vtlnv(l,disc_ie(i))**2) 
        write(6,*)'U,V,uv',utlnv(l,disc_ie(i)),vtlnv(l,disc_ie(i)),uv
        write(6,*)'u,v',ulnv(l,disc_ie(i)),vlnv(l,disc_ie(i))
        write(6,*)'dentro disco'
         write(6,*)'alpha_aux dynamic',alpha_aux
          write(6,*)'al_in,al_in0',al_in,al_in0
        call static(alpha(i), Re, Cl_st2, Cd_st2, ir, ic, h, e, B, G)
          if(abs(disc_theta(i)).lt.180)then
          call cl_dynamic(al_min, al_max,W,Re,al_in,al_in0,
     +  chord,om,Cl_st,ir,ic,h,e,B,G,
     +  foo,icall,alpha_aux,al_ds)
          write(6,*)'foo,disc_theta(i)',foo,disc_theta(i)
          Li=(0.5*foo*ro*chord*W**2)/area_ele    
          D=(0.5*Cd_st2*ro*chord*W**2)/area_ele
         else
          Li=(0.5*Cl_st2*ro*chord*W**2)/area_ele
          D=(0.5*Cd_st2*ro*chord*W**2)/area_ele
         end if

         
         frict_1=(Li**2+D**2)
         frict_3=(sqrt(frict_1))*6.4*nb/((icount+disc_count)*ro*uv)!6.4 is the "depth" change if necessary
         if(abs(y-y0).gt.2.9.and.abs(y-y0).lt.3.2)then
          rfricv(disc_ie(i))=frict_3/qc
         else
          corda=2*sqrt(r**2-dist_2**2)
          rfricv(disc_ie(i))=frict_3*2*r/(corda*qc)         
         endif 
       
        open(unit=41, file='L_theta_full.txt', status='unknown')
        write(41,*) icall,disc_theta(i),alpha(i),foo,Cl_st2,Cd_st2,
     +  Li,D,W,area_ele,u_abs,rfricv(disc_ie(i)),utlnv(l,disc_ie(i)),
     +  vtlnv(l,disc_ie(i)),uv                  
         
         else !if the elemente does not belong to the considered stripe
         goto 154
         end if
        
         
154      continue
         end do  !internal disc cycle
         end do  !stripes cycle  
              
  
163     continue 

        
        open(unit=68, file='velocità1D.txt', status='unknown')
        open(unit=69, file='velocità0.75D.txt', status='unknown')
        open(unit=70, file='velocità3D.txt', status='unknown')
        open(unit=71, file='velocità5D.txt', status='unknown')

        do el=1, ie1Dcount     
        write(68,*)icall,ulnv(l,ie_1D(el)),vlnv(l,ie_1D(el)),ie_1D(el),
     +  ipev(ie_1D(el))
        end do

        do el=1, ie0_75Dcount
        write(69,*)icall,ulnv(l,ie_0_75D(el)),vlnv(l,ie_0_75D(el)),
     +  ie_0_75D(el),ipev(ie_0_75D(el))
        end do

        do el=1, ie3Dcount
        write(70,*)icall,ulnv(l,ie_3D(el)),vlnv(l,ie_3D(el)),
     +  ie_3D(el),ipev(ie_3D(el))
        end do

        do el=1, ie5Dcount
        write(71,*)icall,ulnv(l,ie_5D(el)),vlnv(l,ie_5D(el)),
     +  ie_5D(el),ipev(ie_5D(el))
        end do
       close(68)
       close(69)
       close(70)
       close(71)
  
       open(unit=51, file='rfricv.txt', status='unknown')
        
         do j=1,disc_count
         write(51,*)icall,j,rfricv(disc_ie(j))!,utlnv(l,disc_ie(j)),
c     +   utlnv(l,disc_ie(j)),uv,frict_3    
         end do
  
102     continue
         
        else    !reltive to icall.lt.4000
        goto 114
        end if  !reltive to icall.lt.4000
       
114     continue          
        end   

c***************************************************************** 
c SUBROUTINE FOR STATIC COEFFICIENT 
c the  routine calculates the static lift and drag coefficient 
c using as input the angle of attack and the Reynolds number
c***************************************************************** 

        subroutine static(alpha, Re, Cl_st, Cd_st, ir, ic, h, e, B, G)  

        implicit none                                             
        real h(30), e(5), B(5,30), G(5,30), alpha, Re
        real, intent(out):: Cl_st, Cd_st
        real f_Re, f_alpha, Cl, Cd, angolo, Cl1, Cl2, Cl3
        real Cd1, Cd2, Cd3
        integer  j,  t, ir, ic, segno              
        angolo=abs(alpha)
 
        if (abs(alpha).eq.alpha) then
         segno=+1
        else
         segno=-1 
        endif
        
        if (angolo.gt.180)then
            angolo=180
        endif
        
        if(Re.gt.e(5))then         !this two "if" are usefull to avoid having values of Reynolds out of the range in the text file
           Re=e(5)
        endif
        if(Re.lt.e(1))then 
           Re=e(1)
        endif

        do j= 1, ic
                                                   
          if (angolo.ge.h(j-1).and.angolo.lt.h(j).or.    
     +  angolo.gt.h(j-1).and.angolo.le.h(j))then    
         f_alpha=(angolo-h(j-1))/(h(j)-h(j-1))  
         
        do t= 1, ir
 
          if (Re.ge.e(t-1).and.Re.lt.e(t).or.
     +  Re.gt.e(t-1).and.Re.le.e(t))then
        f_Re=(Re-e(t-1))/(e(t)-e(t-1))
                      
        Cl1=B(t-1,j-1)*(1-f_Re)*(1-f_alpha)           
        Cl2=B(t,j-1)*f_Re*(1-f_alpha)+B(t-1,j)*(1-f_Re)*f_alpha
        Cl3=B(t,j)*f_Re*f_alpha
        Cl=Cl1+Cl2+Cl3 
        Cl_st=segno*Cl              
                                  
        Cd1=G(t-1,j-1)*(1-f_Re)*(1-f_alpha)           
        Cd2=G(t,j-1)*f_Re*(1-f_alpha)+G(t-1,j)*(1-f_Re)*f_alpha
        Cd3=G(t,j)*f_Re*f_alpha
        Cd=Cd1+Cd2+Cd3              
        Cd_st=Cd
        
        goto 103                       
        end if                        
        end do

        end if
        end do

103     continue 

        end 

c***************************************************************** 
c    SUBROUTINE FOR DYNAMIC COEFFICIENT 
c This routine computes the whole Cl-alpha curve in dynamic conditions
c input values for each time step: alpha maximum and minimum, 
c W relative velocity, Omega, Reynolds number
c Rocchio et al. DOI:10.1002/we.2463
c*****************************************************************     

        subroutine cl_dynamic(al_min, al_max,W,Re,al_in,al_in0, 
     +   chord,om,Cl_st,ir,ic,h,e,B,G, 
     +   foo,icall,alpha,al_ds)

        implicit none 
        real tau, k, dt, t_max, U_lev, chord, om, pi, W, Re 
        real clmax, clmin, cl_min, Cl_st, Cd_st, al_min, al_max, Li, D
        real ro
        real al_dot(4000), cls(4000), cl0s(4000), clns
        real x_revers, x_lev(4000), al_in, al_in0
        integer phase(4000), ir, ic, icall 
        real al_ds, clsd(4000), cl0sd, f, fd, cl0d, al_m
        real der_cl(4000), cl_vort(4000), cl_shed(4000)
        real aux1,aux2 
        real aux6,aux7,aux8,aux9,aux10
        real  A(4000),cld(4000) 
        real segno, cl_alpha, ka, kf
        real, save :: al(4000),cl_fin(4000)
        integer imem, imem2
        real foo, foo1, A_amp, A_med
        real rapp_incr,delta_alpha
        real  tau_omega_6
        integer i_vortex_start,vortex_start  
        integer reattached
        real om_3, om_4, om_5, om_6
        real puls
        real cl_par0, cl_par, cl_stallo, alfa0, alfa, alpha, al_stallo
        integer Amp_init!, static_dyn_flag
        integer i,  nt, idt, nt_in
        parameter (pi=3.141592654, ro=998.2) 
        parameter (cl_alpha=0.109655805) !cl_alpha is 2pi*pi/180
        real h(30), e(5), B(5,30), G(5,30) 
        tau=chord/(2*W)
        k=(om*chord)/(2*W)
        dt=(2*pi)/(1000*abs(om))
        t_max=4*2*pi/abs(om)
        U_lev=W/(3*chord)
        nt=t_max/dt 
c-----------------------------------------
c Model parameters
c-----------------------------------------
        om_3=0.08
        om_5=0.12
        puls=2*pi*(0.235*W)/(chord) 
c----------------------------------------- 
c Local variables 
c Static stall angle of the cl-alpha curve
c-----------------------------------------
        cl_par0=0
        alfa=0   
        alfa0=0
        do i=1,26     
        alpha=alpha+1  
        call static(alpha, Re, Cl_st, Cd_st, ir, ic, h, e, B, G)
        cl_par=Cl_st 
        if(cl_par.gt.cl_par0) then
         cl_par0=cl_par
        alfa0=alpha
        endif
        end do
        
        al_stallo=alfa0
        cl_stallo=cl_par0
c----------------------------------------- 
c Initialization
c-----------------------------------------
                al(:)=0
                al_dot(:)=0   
                cls(:)=0
                cl0s(:)=0
                phase(:)=0
                x_lev(:)=0
                clsd(:)=0
                cld(:)=0
                der_cl(:)=0
                cl_vort(:)=0
                cl_shed(:)=0
                cl_fin(:)=0
                i_vortex_start=0
                vortex_start=0
                reattached=0
c-----------------------------------------
c calculation of the dynamic CL-alpha curve
c-----------------------------------------
         do i=1, nt
c        Amplitude and average value of the pitching motion 
          if(icall.eq.1) then 
               A_amp=36.0 ! deg 
               A_med=0.0  ! deg 
          else
                A_amp=0.5*(al_max-al_min)
                A_med=al_max-A_amp
          endif
c	 calcultion dynamic stall alpha value: 
          if(Re.lt.490000) then
                al_ds=320*(k*A_amp*pi/180)+12
          else if (Re.ge.490000.and.Re.lt.980000) then
                al_ds=320*(k*A_amp*pi/180)+14.5
          else if (Re.ge.980000.and.Re.lt.2500000) then
                al_ds=320*(k*A_amp*pi/180)+18
          else
                al_ds=320*(k*A_amp*pi/180)+19
          endif
        
c	calculation minimum alpha:
          if(Re.lt.490000) then
                al_m=-372*(k*A_amp*pi/180)+18
          else if (Re.ge.490000.and.Re.lt.980000)then
                al_m=-372*(k*A_amp*pi/180)+20
          else if (Re.ge.980000.and.Re.lt.2500000)then
                al_m=-372*(k*A_amp*pi/180)+22
          else
                al_m=-372*(k*A_amp*pi/180)+25
          end if
c	calculation of cl_minimo:
          if(Re.le.980000)then
                cl_min=-15*(k*A_amp*pi/180)+0.75
          else
                cl_min=-20*(k*A_amp*pi/180)+1.2
          endif
 
c	frequency omega_6: 
        tau_omega_6=2*pi*(A_amp+A_med-al_m)/(4*A_amp*3*om)
        
c	frequency omega_4: 
        om_4=(W/3)*tau 
  
c 	sinusoidal motion 
        al(i)=A_med+A_amp*sin(om*(i-1)*dt+0) 
  
c	alpha derivative in rad/s 
        al_dot(i)=+A_amp*om*cos(om*(i-1)*dt+0)*pi/180
c-----------------------------------------
c  calculation of static lift coefficient   
c-----------------------------------------              
        if (abs(al(i)).eq.al(i)) then
        segno=+1
        else
        segno=-1 
        endif
       
        call static(al(i), Re, Cl_st, Cd_st, ir, ic,  h, e, B, G)

        cls(i)=segno*Cl_st      
c	linear lift coefficient 
        cl0s(i)=cl_alpha*al(i)
        if(i.eq.0) then 
                cld(i)=cls(i)
                cl_vort(i)=0
                cl_shed(i)=0
                if(abs(cld(i)).gt.100)then
                write(6,*)'check A',cld(i),cld(i-1)
                endif                    
        else
                if(abs(al(i)).gt.abs(al(i-1)))then !increasing alpha
                 reattached=0        
                        if(vortex_start.eq.0)then
                        x_lev(i)=0
                        phase(i)=1
                        cld(i)=cl0s(i)
                        der_cl(i)=cl0s(i)-cls(i) 
                        aux10=(cl_vort(i-1)+(om_3/tau)*dt*der_cl(i))
                        cl_vort(i)=aux10/(1+(om_3/tau)*dt)         
        
                                if(abs(al(i)).gt.al_ds) then
                                vortex_start=1
                                i_vortex_start=i
                                end if
                        else
                         phase(i)=2
               cld(i)=(cld(i-1)+(om_4/tau)*dt*cls(i))/(1+(om_4/tau)*dt) 
               x_lev(i)=x_lev(i-1)+U_lev*dt
        
                  if(x_lev(i-1).lt.1) then
                  der_cl(i)=cld(i)-cls(i)
                  aux9=(cl_vort(i-1)+(om_3/tau)*dt*der_cl(i))
                  cl_vort(i)=aux9/(1+(om_3/tau)*dt)
                  cl_shed(i)=cl_vort(i)*sin(puls*(i-i_vortex_start)*dt)

                  else
        
                 der_cl(i)=0
                aux8=(cl_vort(i-1)+(om_3/tau)*dt*der_cl(i)) 
                cl_vort(i)=aux8/(1+(om_3/tau)*dt) 
              cl_shed(i)=cl_vort(i)*sin(puls*(i-i_vortex_start)*dt)

                  endif
               endif
        
        
            else !decreasing alpha

      if(abs(al(i)).gt.al_m.or.x_lev(i-1).lt.1.and.phase(i-1).eq.2)then
        phase(i)=3
        vortex_start=0
        reattached=0
        x_lev(i)=x_lev(i-1)+U_lev*dt
        tau_omega_6=2*pi*(A_amp+A_med-al_m)/(4*A_amp*3*om)
        om_6=1/(abs(tau_omega_6))
        cld(i)=(cld(i-1)+(om_6*dt*segno*cl_min))/(1+(om_6*dt))
        der_cl(i)=0
        aux1=abs(cl_shed(i-1))-abs(cl_shed(i-2)) 
            if(aux1.gt.0.and.x_lev(i).gt.1)then 
                aux7=(cl_vort(i-1)+(om_3/tau)*dt*der_cl(i))
                cl_vort(i)=aux7/(1+(om_3/tau)*dt)
            else
               aux2=cl_vort(i-1)+8*(om_3/tau)*dt*der_cl(i)
               cl_vort(i)=(aux2)/(1+8*(om_3/tau)*dt)
               cl_shed(i)=cl_vort(i)*sin(puls*(i-i_vortex_start)*dt)
           endif

        else !al lower than al_m

            if(reattached.eq.0) then
            cld(i)=(cld(i-1)+(om_5/tau)*dt*cls(i))/(1+(om_5/tau)*dt)
            else
            cld(i)=cls(i)
            x_lev(i)=x_lev(i-1)+U_lev*dt
            phase(i)=4
            der_cl(i)=0
               if(abs(cld(i)).ge.abs(cls(i))) then
                reattached=1
                cld(i)=cls(i)
                aux6=(cl_vort(i-1)+(om_3/tau)*dt*der_cl(i)) 
               cl_vort(i)=aux6/(1+(om_3/tau)*dt)   
               cl_shed(i)=cl_vort(i)*sin(puls*(i-i_vortex_start)*dt)
               endif
            endif

        endif      ! relative to alpha min
        endif      ! relative to check upstroke/downstroke 

        cl_fin(i)=cld(i)+cl_shed(i)
        
        endif      ! relative to i=0
        end do
c-----------------------------------------        
c this session uses the dynamic curve found and gives back the dynamic coefficient (called foo) relative to the alpha of interest (corresponding to al_in) 
c-----------------------------------------        
        nt_in=3*nt/4  
                       
        do i=nt_in, nt

                delta_alpha=al(i+1)-al(i-1)
                if(abs(delta_alpha).lt.0.001)then
                cl_fin(i)=0.5*(cl_fin(i+1)+cl_fin(i-1))
                else
                rapp_incr=(cl_fin(i+1)-cl_fin(i-1))/(delta_alpha)
                cl_fin(i)=cl_fin(i-1)+rapp_incr*(al(i)-al(i-1))
                end if
                if(((A_amp+A_med)-al(i)).lt.0.0001)then
                        clmax=cl_fin(i)  
                else if(abs((A_med-A_amp)-al(i)).lt.0.001)then
                        clmin=cl_fin(i)  
                end if
                
        end do
 
                if(abs(al_in).eq.al_in)then
                segno=+1
                else
                segno=-1
                end if
                 
           if(al_in.ge.(A_med+A_amp))then
                foo=clmax
                goto 108 
           else if (al_in.le.(A_med-A_amp))then
                foo=clmin
                goto 108 
           else if(al_in.ge.0)then
c-----------------------------------------        
c UPSTROKE: 
c-----------------------------------------
             if(al_in.ge.al_in0)then
                
               do i=nt_in, (nt_in+(nt-nt_in)/2) 
                    if(al(i-1).le.al_in.and.al(i).gt.al_in)then
        foo1=(cl_fin(i)-cl_fin(i-1))/(al(i)-al(i-1))*(al_in-al(i-1))
        foo=foo1+cl_fin(i-1)
                    goto 108
                    end if
                end do
                do i=(nt_in+(nt-nt_in)*3/4), nt 
                    if(al(i-1).le.al_in.and.al(i).gt.al_in)then
        foo1=(cl_fin(i)-cl_fin(i-1))/(al(i)-al(i-1))*(al_in-al(i-1))
        foo=foo1+cl_fin(i-1)
                    goto 108
                   end if
                end do
c-----------------------------------------        
c DOWNSTROKE: 
c-----------------------------------------
             else if(al_in.lt.al_in0)then

               do i=(nt_in+(nt-nt_in)/4), (nt_in+(nt-nt_in)*3/4)
                        
                  if(al(i-1).ge.al_in.and.al(i).lt.al_in)then
        foo1=(cl_fin(i)-cl_fin(i-1))/(al(i)-al(i-1))*(al_in-al(i-1))
        foo=foo1+cl_fin(i-1)
                  goto 108
                  end if
               end do
             end if !relative to upstroke or downstroke
          elseif(al_in.lt.0)then
c-----------------------------------------        
c UPSTROKE: 
c-----------------------------------------	        
               if(al_in.le.al_in0)then

                do i=(nt_in+(nt-nt_in)/2),nt
                    if(al(i).le.al_in.and.al(i-1).gt.al_in)then
        foo1=(cl_fin(i)-cl_fin(i-1))/(al(i)-al(i-1))*(al_in-al(i-1))
        foo=foo1+cl_fin(i-1)
                   goto 108
                   end if                          
                  end do       
c-----------------------------------------        
c DOWNSTROKE: 
c-----------------------------------------      
               elseif(al_in.gt.al_in0)then
                   do i=(nt_in+(nt-nt_in)/2), nt
                 
                    if(al(i).gt.al_in.and.al(i-1).le.al_in)then
        foo1=(cl_fin(i)-cl_fin(i-1))/(al(i)-al(i-1))*(al_in-al(i-1))
        foo=foo1+cl_fin(i-1)
                    goto 108
                    end if                              
                   end do                                 
                        
                   do i=nt_in, (nt_in+(nt-nt_in)/2)
                   
                    if(al(i).gt.al_in.and.al(i-1).le.al_in)then
        foo1=(cl_fin(i)-cl_fin(i-1))/(al(i)-al(i-1))*(al_in-al(i-1))
        foo=foo1+cl_fin(i-1)
                    goto 108
                    end if
                  end do                                
108     continue                             
               end if !relative to upstroke or downstroke   
          end if      !relative to al_in>0 or <0
        if(abs(foo).gt.100)then  
        write(6,*)'ERRORE abs(FOO)>100'
        foo=0.1
        endif
       write(6,*)'rissmooth', cl_fin(i),cl_fin(i-1),al(i),al(i-1),al_in
        end
