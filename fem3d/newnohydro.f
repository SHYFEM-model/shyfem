
!--------------------------------------------------------------------------
!
!    Copyright (C) 1991-1992,2001,2013-2016,2019  Georg Umgiesser
!    Copyright (C) 2013,2018  Debora Bellafiore
!    Copyright (C) 2016,2018  William McKiver
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

c routines for non hydrostatic terms
c
c revision log :
c
c 18.02.1991	ggu	(from scratch)
c 04.06.1991	ggu	(c=(1) : friction term has been corrected)
c 01.10.1992	ggu	(staggered FE - completely restructured)
c 12.01.2001	ggu	solve for znv and not level difference (ZNEW)
c 10.05.2013	dbf	written from scratch
c 31.05.2013	dbf	written from scratch
c 12.09.2013	ggu	changed VERS_6_1_67
c 18.06.2014	ggu	changed VERS_6_1_77
c 19.12.2014	ggu	changed VERS_7_0_10
c 23.12.2014	ggu	changed VERS_7_0_11
c 19.01.2015	ggu	changed VERS_7_1_3
c 30.04.2015	ggu	changed VERS_7_1_9
c 14.06.2016	ggu	changed VERS_7_5_14
c 17.06.2016	ggu&wmk	adapted to new version
c 18.12.2018	dbf&wmk	adapted to last version and inserted in develop
c 14.02.2019	ggu	changed VERS_7_5_56
c 16.02.2019	ggu	changed VERS_7_5_60
c 13.03.2019	ggu	changed VERS_7_5_61
c 21.05.2019	ggu	changed VERS_7_5_62
c 13.03.2021	clr&ggu	adapted for petsc solver
c 23.04.2021    clr     alternative implementation to replace pragma directives
c
c********************************************************************

	subroutine nonhydro_init


c initializes non hydrostatic pressure terms

        use mod_system !DWNH
        use mod_nohyd !DWNH
        use levels !DWNH
        use basin !DWNH

	implicit none

	integer inohyd
	real getpar

	inohyd = nint(getpar('inohyd'))

	!if( inohyd /= 0 ) then
	!if( inohyd == 0 ) then !DWNH
	!  write(6,*) 'inohyd = ',inohyd
	!  stop 'error stop nonhydro_init: cannot run non-hydrostatic'
	!end if

c       --------------------------------------------
c       initialize variables
c       --------------------------------------------

        bnohydro = ( inohyd /= 0 ) !DWNH

        bsys3d = bnohydro !DWNH

        qpov = 0. !DWNH
        qpnv = 0. !DWNH
	qdistv = 0.

	end

c********************************************************************

	subroutine nonhydro_get_flag(bnohyd)

	use mod_nohyd
	
	implicit none

	logical bnohyd

	bnohyd = bnohydro !DWNH

	end

c********************************************************************

	subroutine nonhydro_adjust

c integrates non hydrostatic adjustment to equations

        use mod_hydro
        use mod_nohyd
        use levels
        use basin
        use mod_bound_dynamic
        use mod_system

        implicit none

c local
        integer k,l,lmax
        integer iw,inhadj
        real aq
        real qvmax
        real getpar

c        non-hydrostatic semi-implicit weight
        aq = getpar('aqpar')
       
c        parameter for how to adjust u,v,eta after computation of NH pressure
        inhadj = nint(getpar('inhadj'))

c        solve for NH pressure q
        call system_init
        call nonhydro_prepare_matrix
	call system_solve_3d(nkn,nlvdi,nlv,qpnv)
	!call system_adjust_3d(nkn,nlvdi,nlv,qpnv)
	call system_get_3d(nkn,nlvdi,nlv,qpnv) !DWNH
	
c        Compute NH u,v,eta
        if (inhadj .eq. 1) then
c          recompute u,v,eta using updated NH pressure
	  do 
	    call hydro_transports
	    call setnod		
	    call set_link_info
	    call adjust_mass_flux	
            call system_init 
            call hydro_zeta(rqv) 
            !call system_solve_z(nkn,znv) 
            call system_solve(nkn,znv) !DWNH
            !call system_adjust_z(nkn,znv) 
            call system_get(nkn,znv) !DWNH
            call setweg(1,iw)
	    if( iw == 0 ) exit

	  end do
	  call hydro_transports_final
        else
c          make NH correction to u,v,eta using updated NH pressure
          call nonhydro_correct_uveta
        endif

        call setzev 
        call setuvd 
	call baro2l
  	call make_new_depth 
	call check_volume

c        Finally adjust NH pressure q
        qvmax=0.0
        do k = 1,nkn
          lmax = ilhkv(k)
          do l=1,lmax
            qpnv(l,k) = qpnv(l,k)-qpnv(1,k)
            qvmax=max(qvmax,abs(qpnv(l,k)))
          end do
        enddo
        write(6,*) 'qmax=',qvmax

	end


c********************************************************************

	subroutine nonhydro_copy

c copies new values of q to old time step

        use mod_hydro
        use mod_nohyd
        use levels
        use basin

        implicit none

        qpov = qpnv
        qpnv = 0.
        
        end

c********************************************************************

	subroutine sp256wnh 

c solves for w using momentum equation

        use mod_hydro
        use mod_nohyd
        use mod_depth
        use mod_area
        use mod_ts
        use mod_hydro_vel
        use mod_layer_thickness
        use mod_diff_visc_fric
        use evgeom
        use levels
        use basin
        use mod_internal

        implicit none

c parameters
        include 'pkonst.h'
        include 'femtime.h'

c local
	logical blast
        integer k,ie,ii,l
        integer ilevel,lmax
        integer kn(3)
        integer ivwadv,inhflx
	integer itot,isum

        real aj,aj4
        real b(3),c(3),f(3),fl(3)
	real wplus,wminus,wr
        real dt
        real wsum
        real dzo,dzn,dz
        real dzobb,dzott,dznbb,dzntt
        real wadvv,wadvh,wdifh,wq
        real wb,wc
        real wvmax,wmax,werr
        real ahpar,aq,aqt,az,azt
	real anu,rnu
	real us,vs
        real utn,ubn,uto,ubo
        real vtn,vbn,vto,vbo

        real wrhs(0:nlv,nkn),wdiag(0:nlv,nkn)
        real wdiagt(0:nlv,nkn),wdiagb(0:nlv,nkn)
        real wtop(nlv),wbot(nlv),wcen(nlv),wlrhs(nlv)
        real wgam(nlv),wsol(nlv)
        real wl(0:nlv,3)
        real hdn,hdpn,hdo,hdpo
        real getpar

c--------------------------------------------------------
c initialize parameters	
c--------------------------------------------------------

        call get_timestep(dt)

c        non-hydrostatic semi-implicit weight       
        aq = getpar('aqpar')
        aqt = 1. - aq
	az = getpar('azpar')
	azt = 1. - az

c        horizontal viscosity parameter
        ahpar = getpar('ahpar')

c        horizontal advection method (upwind flux scheme/averaging)
        inhflx = nint(getpar('inhflx'))

c        vertical advection method (upwind/finite diff)
        ivwadv = nint(getpar('ivwadv'))

c--------------------------------------------------------
c initialize arrays	
c--------------------------------------------------------

        wrhs   = 0. 
        wdiagt = 0.
        wdiag  = 0.
        wdiagb = 0.

c----------------------------------------------------------
c loop over elements to compute contributions 
c----------------------------------------------------------
        do ie=1,nel
          ilevel=ilhv(ie)
          do ii=1,3
            k=nen3v(ii,ie)
            kn(ii)=k
            b(ii)=ev(ii+3,ie)
            c(ii)=ev(ii+6,ie)
          end do
 
          aj=ev(10,ie) !area of triangle / 12
          aj4=4.*aj

c         --------------------------------------------------
c 	  set vertical velocities in element
c         --------------------------------------------------
	
	  do l=0,ilevel
	    do ii=1,3
	      k = kn(ii)
	      wl(l,ii) = wlov(l,k)
	    end do
	  end do

c         --------------------------------------------------
c 	  loop over levels
c         --------------------------------------------------
          do l=1,ilevel

	    blast = l .eq. ilevel

c            Set layer thicknesses, transport and viscosity values
	    hdn = hdenv(l,ie)
            hdo = hdeov(l,ie)
            utn=utlnv(l,ie)
            uto=utlov(l,ie)
            vtn=vtlnv(l,ie)
            vto=vtlov(l,ie)
            if (blast) then
              hdpn = 0.0
              hdpo = 0.0
              ubn=0.0
              ubo=0.0
              vbn=0.0
              vbo=0.0
	      anu = 0.5 * ahpar * difhv(l,ie)
            else
              hdpn = hdenv(l+1,ie)
              hdpo = hdeov(l+1,ie)
              ubn=utlnv(l+1,ie)
              ubo=utlov(l+1,ie)
              vbn=vtlnv(l+1,ie)
              vbo=vtlov(l+1,ie)
	      anu = 0.5 * ahpar * ( difhv(l,ie) + difhv(l+1,ie) )
            end if
            dzn = 0.5 * (hdn+hdpn)
            dzo = 0.5 * (hdo+hdpo)

c            semi-implicit (weighted) transport
            us=0.5*(az*(utn+ubn)/dzn+azt*(uto+ubo)/dzo)
            vs=0.5*(az*(vtn+vbn)/dzn+azt*(vto+vbo)/dzo)
            wb = 0.0
	    wc = 0.0
            wsum = 0.0
	    itot = 0
	    isum = 0
	    do ii=1,3
              f(ii)=us*b(ii)+vs*c(ii)
              if(f(ii).lt.0.) then  !flux out of node
                itot=itot+1
                isum=isum+ii
              end if
	      wb = wb + ( b(ii) * wl(l,ii) )
	      wc = wc + ( c(ii) * wl(l,ii) )
              wsum=wsum+wl(l,ii)
	    end do

            if (inhflx .eq. 1) then
c	      ----------------------------------------
c	      compute horizontal fluxes on top and bottom interface
c	      ----------------------------------------
              if(itot.eq.1) then
	        do ii=1,3
	          fl(ii) = f(ii) * wl(l,isum)
	        end do
              else if(itot.eq.2) then
                isum=6-isum
	        do ii=1,3
	          fl(ii) = f(ii) * wl(l,ii)
	        end do
                fl(isum) = 0.0
                fl(isum) = -(fl(1)+fl(2)+fl(3))
                isum=6-isum
	      else
	        fl = 0.0
	      end if
            end if

            do ii=1,3

              k=kn(ii)
c	      ----------------------------------------
c	      explicit non-hydrostatic pressure
c	      ----------------------------------------

	      if( blast ) then
		wq = 0.0
	      else
	        wq = - aqt * qdistv(k) * (qpov(l,k) - qpov(l+1,k))/dzo
	      end if

c             ----------------------------------------
c             vertical advection term
c             ----------------------------------------

	      wadvv = 0.
              if (ivwadv .ge. 1) then
                if (ivwadv .eq. 1) then
c                  vertical advection computed using upwind scheme
                  if (blast) then
                    wadvv = 0.0
                  else
	            wplus = 0.5 * (wlov(l,k)+wlov(l-1,k))/dzo
	            if( wplus .gt. 0. ) then
		      wadvv = wplus*wlov(l,k)
                    else
 		      wadvv = wplus*wlov(l-1,k)
	            end if
	            wminus = 0.5 * (wlov(l,k)+wlov(l+1,k))/dzo
	            if( wminus .gt. 0. ) then
	              wadvv = wadvv - wminus*wlov(l+1,k)
	            else
	  	      wadvv = wadvv - wminus*wlov(l,k)
	            end if
                  end if
                else         
c                  vertical advection computed using finite differencing
                  if (blast) then
                    wadvv = 0.0
                  else
                    wadvv= 0.5*wlov(l,k)*(wlov(l-1,k)-wlov(l+1,k))/dzo
                  end if 
                end if
              end if

c	      ----------------------------------------
c	      horizontal diffusion
c	      ----------------------------------------
	      wdifh = - 3. * anu * (  wb * b(ii) + wc * c(ii) )

c	      ----------------------------------------
c	      horizontal advection
c	      ----------------------------------------

              if (inhflx .eq. 1) then
	        wadvh = - 3. * fl(ii)
              else
	        wadvh = - (b(ii) * us + c(ii) * vs) * wsum 
              end if

c	      ----------------------------------------
c	      complete explicit contribution
c	      ----------------------------------------

	      wr = dzn * ( wlov(l,k) + dt * (wq + wdifh - wadvh - wadvv))

c             ----------------------------------------
c             define tri-diagonal coefficients
c             ----------------------------------------
	      wrhs(l,k)   = wrhs(l,k) + aj4  * wr

              if (blast) then
c                set w to zero on bottom layer 
       	        wdiag(l,k)  = 1.0
	        wdiagt(l,k) = 0.0
	        wdiagb(l,k) = 0.0
                wrhs(l,k) = 0.0
              else
                rnu=visv(l,k)
		wdiag(l,k)  = wdiag(l,k) + dzn * aj4 + dt*aj4*rnu
     +                       * ( 1./hdn + 1./hdpn )
	        wdiagt(l,k) = wdiagt(l,k) - dt*aj4*rnu/hdn 
	        wdiagb(l,k) = wdiagb(l,k) - dt*aj4*rnu/hdpn
              end if
c              output to check contributions 
c              if (k .eq. 1801) then
c                if (l.eq.1.or.l.eq.26.or.l.eq.51.or.l.eq.76) then
c                  write(688,*) t_act,l,wrhs(l,k),wq,wdifh,wadvh,wadvv
c                end if
c              end if
            end do
          end do
        end do

c-----------------------------------------------------------
c solve tri-diagonal system
c-----------------------------------------------------------
        wmax=0.0
        wvmax=0.0
	do k=1,nkn
	  lmax = ilhkv(k)
	  if (lmax .gt. 2) then
            do l=1,lmax-1
              wtop(l)=wdiagt(l,k)
              wbot(l)=wdiagb(l,k)
              wcen(l)=wdiag(l,k)
              wlrhs(l)=wrhs(l,k)
            enddo
c            call tridiag(wtop,wcen,wbot,wlrhs,wsol,lmax-1,werr)
	    call tridag(wtop,wcen,wbot,wlrhs,wsol,wgam,lmax-1)
            do l=1,lmax-1
              wlnv(l,k)=wsol(l)
              wvmax=max(wvmax,abs(wlnv(l,k)))
            enddo 
	  else
	    wlnv(1,k) = 0.0 
	    werr = 0.0
	  endif
c          wmax=max(werr,wmax)
c	  wlnv(0,k) = ( znv(k) - zov(k) ) / dt
   	  wlnv(0,k) = 0.0
          wlnv(lmax,k) = 0.0                    
        end do

        write(6,*) 'wmax=',wvmax

c--------------------------------------------------------
c end of routine
c--------------------------------------------------------

	end

c********************************************************************

	subroutine nonhydro_set_explicit 

c adds explicit part of non hydrostatic pressure to explict terms

        use mod_hydro
        use mod_nohyd
        use mod_layer_thickness
        use mod_internal
        use evgeom
        use levels
        use basin

        implicit none
        
        include 'pkonst.h'
        include 'femtime.h'

c local
        integer k,l,ie,ii,lmax
        real hhio,hhin
        real qxonh,qyonh ! non hydrostatic pressure gradient old time
        real qxnnh,qynnh ! non hydrostatic pressure gradient new time
        real qpo,qpn
        real b,c
        real aq,aqt
        real getpar

c        non-hydrostatic semi-implicit weight       
        aq = getpar('aqpar')
        aqt = 1. - aq

        do ie=1,nel
          lmax = ilhv(ie)
          do l=1,lmax
            qxonh = 0.0
            qyonh = 0.0
            qxnnh = 0.0
            qynnh = 0.0
            do ii=1,3
              k = nen3v(ii,ie)
              qpo = qdistv(k)*qpov(l,k)
              qpn = qdistv(k)*qpnv(l,k)
              b = ev(3+ii,ie)
              c = ev(6+ii,ie)
              qxonh = qxonh + b * qpo
              qyonh = qyonh + c * qpo
              qxnnh = qxnnh + b * qpn
              qynnh = qynnh + c * qpn
            end do
            hhio = hdeov(l,ie)  
            hhin = hdenv(l,ie)

            fxv(l,ie) = fxv(l,ie) + aqt * hhio * qxonh 
     +                            + aq  * hhin * qxnnh 
            fyv(l,ie) = fyv(l,ie) + aqt * hhio * qyonh
     +                            + aq  * hhin * qynnh 

          end do
        end do

	end 

c********************************************************************

	subroutine nonhydro_prepare_matrix

c assembles linear system matrix
c
c vqv		flux boundary condition vector
c
c semi-implicit scheme for 3d model

	use mod_internal
	use mod_depth
	use mod_area
	use evgeom
	use levels
	use basin
	use mod_layer_thickness
        use mod_hydro
        use mod_hydro_vel
        use mod_zeta_system, only : solver_type

	implicit none

	real pvar(nlvdi,nkn)

	include 'mkonst.h'
	include 'pkonst.h'
	include 'femtime.h'
 
	integer kn(3)
	integer ie,i,j,j1,j2,n,m,kk,l,k
	integer lmax
	real aj,aj4,aj12
	real ht
	real hik(3)
	real hia3d(-1:+1,3,3)
	real b(3),c(3),z(3)
        real dt
        real fhn,fho,ffh,fvn,fvo,ffv
        real aq,aqt,az,azt

	real hd,hdm,hdp,rhm,rhp,rhc
	real hldaux(0:nlv+1)
        real volo,voln,dvdt
	real volnode

	integer locsps,loclp,iround
	real getpar
	integer loccoo3d 

        if(trim(solver_type)=='PETSc')then
        write(6,*)'nonhydro_prepare_matrix uses function loccoo3d'
        write(6,*)'but no garanty is given of the return values'
        write(6,*)'given that PETSc solver is used instead of SPK'
        stop "Program ends, please check what happens with loccoo3d"    
        endif

	hldaux = 0.

c        non-hydrostatic semi-implicit weights
        aq = getpar('aqpar')
        aqt= 1. - aq

c        transport terms semi-implicit weights
	az = getpar('azpar')
	azt = 1. - az

        call get_timestep(dt)

c-------------------------------------------------------------
c loop over elements
c-------------------------------------------------------------

	do ie=1,nel

c	  ------------------------------------------------------
c	  initialize element values
c	  ------------------------------------------------------
	  hldaux(1:nlv)=hdenv(1:nlv,ie)
	  aj=ev(10,ie)
          aj4=4.*aj
          aj12=12.*aj
	  do i=1,3
	    kk=nen3v(i,ie)
	    kn(i)=kk
	    b(i)=ev(i+3,ie)
	    c(i)=ev(i+6,ie)
	  end do

c	  ------------------------------------------------------
c	  set element matrix and RHS
c	  ------------------------------------------------------

	  lmax = ilhv(ie)
	  do l=1,lmax
	    hia3d = 0.

	    hd = hldaux(l)
	    hdm = hldaux(l-1)
	    hdp = hldaux(l+1)

  	    rhm = 2. / ( hd + hdm )
	    if( l == 1 ) rhm = 0. 
	    rhp = 2. / ( hd + hdp )
	    if( l == lmax ) rhp = 0.
	    rhc = rhm + rhp

	    do n=1,3
	      do m=1,3
	        kk =loccoo3d(n,m,kn,l,ie)
	        hia3d(0,n,m) = -aj12*hd*dt*aq*az*(b(n)*b(m)+c(n)*c(m))

	        if (n.eq.m) then
  	          hia3d(0,n,m) = hia3d(0,n,m) - aj4*dt*aq*rhc
	          hia3d(-1,n,m) = aj4 * dt * aq * rhm
	          hia3d(+1,n,m) = aj4 * dt * aq * rhp
                endif
	      end do
              k=kn(n)

              fhn = utlnv(l,ie)*b(n) + vtlnv(l,ie)*c(n)
              fho = utlov(l,ie)*b(n) + vtlov(l,ie)*c(n)

              fvn = wlnv(l-1,k)-wlnv(l,k)
              fvo = wlov(l-1,k)-wlov(l,k)

              ffh = - (az * fhn + azt * fho)
              ffv = fvn

              hik(n) = aj4 * (3. * ffh + ffv) 

	    end do

c	    ------------------------------------------------------
c	    boundary conditions
c	    ------------------------------------------------------

            if (l == 1) then
  	      do i=1,3
                k=kn(i)
                hia3d(0,i,i)=hia3d(0,i,i)-aj4/(grav*dt)
                hik(i) = hik(i) + aj4*(znv(k)-zov(k))/dt
	      end do
	    end if

c	    ------------------------------------------------------
c	    in hia(i,j),hik(i),i,j=1,3 is system
c	    ------------------------------------------------------

	    call system_assemble_3d(ie,l,nlv,kn,hia3d,hik)

	  end do

	end do

c-------------------------------------------------------------
c end of loop over elements
c-------------------------------------------------------------

	call system_adjust_matrix_3d	

c-------------------------------------------------------------
c end of routine
c-------------------------------------------------------------

	end

c********************************************************************

        subroutine nonhydro_adjust_value

	end

c********************************************************************

	subroutine nonhydro_correct_uveta

        use mod_hydro
        use mod_nohyd
        use mod_layer_thickness
        use mod_area
        use mod_hydro_vel
        use evgeom
        use levels
        use basin

        implicit none

        include 'pkonst.h'

c local
        integer kn(3)
        integer ie,l,ilevel,ii,k,lmax
        real b(3),c(3)
        real dt
        real uqaux,vqaux,zqaux
        real getpar
        real aq

c        non-hydrostatic semi-implicit weight  
        aq = getpar('aqpar')

        call get_timestep(dt)

c        compute NH correction to znv
        do k = 1,nkn
          zqaux=qdistv(k)*qpnv(1,k) / grav
          znv(k) = znv(k) + zqaux
        enddo

c        compute NH correction to w
c        do  k = 1,nkn
c          lmax = ilhkv(k)
c          do l=1,lmax-1
c            dzcc = (hdknv(l+1,k)+hdknv(l,k))/2.
c	    wlnv(l,k)=wlnv(l,k) - aq * dt * (qpnv(l,k)-qpnv(l+1,k))/dzcc
c          enddo
c          wlnv(0,k)=(znv(k)-zov(k))/dt
c        enddo

c        compute NH correction to eta, u,v
        do ie=1,nel
          do ii=1,3
            k=nen3v(ii,ie)
            kn(ii)=k
            b(ii)=ev(ii+3,ie)
            c(ii)=ev(ii+6,ie)
          end do

          ilevel=ilhv(ie)
          do l=1,ilevel
            uqaux=0.0
            vqaux=0.0
            do ii=1,3
              uqaux = uqaux + (b(ii)*qdistv(kn(ii))*qpnv(l,kn(ii)))
              vqaux = vqaux + (c(ii)*qdistv(kn(ii))*qpnv(l,kn(ii)))
            enddo
            utlnv(l,ie) = utlnv(l,ie) - aq * dt * uqaux * hdenv(l,ie)
            vtlnv(l,ie) = vtlnv(l,ie) - aq * dt * vqaux * hdenv(l,ie) 
          enddo
        enddo
        
	end

c********************************************************************

	subroutine nh_handle_output(dtime)

	implicit none

	double precision dtime

	integer, save :: ia_out(4)
	double precision, save :: da_out(4)

	integer, save :: icall = 0
	integer, save :: iwvel,iqpnv

	real getpar

c       --------------------------------------------
c	initialize output files
c	--------------------------------------------

        if(icall.eq.0) then	!first time

          iwvel=nint(getpar('iwvel'))
          iqpnv=nint(getpar('iqpnv'))

	  call nh_open_output(ia_out,da_out,iwvel,iqpnv)

          icall=icall+1

	end if

c       ----------------------------------------------------------
c       write results to file
c       ----------------------------------------------------------

	call nh_write_output(dtime,ia_out,da_out,iwvel,iqpnv)

c       ----------------------------------------------------------
c       end of routine
c       ----------------------------------------------------------

	end

c********************************************************************

	subroutine nh_open_output(ia_out,da_out,iwvel,iqpnv)
	
c opens output of w/q

	use levels


	implicit none

	integer ia_out(4)
	double precision da_out(4)
	integer iwvel,iqpnv

	integer nvar,id
	logical has_output
	logical has_output_d
	real getpar

	nvar = 0
	if( iwvel .gt. 0 ) nvar = nvar + 1
	if( iqpnv .gt. 0 ) nvar = nvar + 1

	call init_output_d('itmcon','idtcon',da_out)

	if( has_output_d(da_out) ) then
	  call shyfem_init_scalar_file('nhyd',nvar,.false.,id)
	  da_out(4) = id
	end if

	end

c*******************************************************************	


	subroutine nh_write_output(dtime,ia_out,da_out,iwvel,iqpnv)

c writes output of wvel and qpnv

	use levels
	use mod_hydro_vel
	use mod_hydro_print
	use mod_nohyd

	implicit none

	double precision dtime
	integer ia_out(4)
	double precision da_out(4)
	integer iwvel,iqpnv
        integer iwrt,inhwrt
        save iwrt
        data iwrt /1/

	integer id,l
	logical next_output
	logical next_output_d
	real getpar

	inhwrt = nint(getpar('inhwrt'))

	if( iwvel .gt. 0 ) then
          do l=1,nlv
            wprv(l,:)=0.5*(wlnv(l,:)+wlnv(l-1,:))
          end do
	end if

        if (inhwrt .gt. 0) then
          if (iwrt .eq. inhwrt ) then
	    if( iwvel .gt. 0 ) then
	      call write_scalar_file(ia_out,14,nlvdi,wprv)
	    end if
	    if( iqpnv .gt. 0 ) then
	      call write_scalar_file(ia_out,15,nlvdi,qpnv)
	    end if
            iwrt=0
          end if
          iwrt=iwrt+1 
	else
	  if( next_output(ia_out) ) then
	    if( iwvel .gt. 0 ) then
	      call write_scalar_file(ia_out,14,nlvdi,wprv)
	    end if
	    if( iqpnv .gt. 0 ) then
	      call write_scalar_file(ia_out,15,nlvdi,qpnv)
	    end if
	  end if
	end if

	if( next_output_d(da_out) ) then
	  id = nint(da_out(4))
	  if( iwvel .gt. 0 ) then
	    call shy_write_scalar_record(id,dtime,14,nlvdi,wprv)
	  end if
	  if( iqpnv .gt. 0 ) then
	    call shy_write_scalar_record(id,dtime,15,nlvdi,qpnv)
	  end if
	end if

	end

c********************************************************************


        subroutine qhdist(qdist)

c makes distance array from open boundaries for NH pressure
c
c qdist is contained between 0 and 1
c if nqdist (=d) is not given qdist = 1 (default)
c otherwise the first d2=d/2 rows of nodes have qdist = 0
c and then the next ones have qdist = i/d with i = d2+1, d2+d
c the rest has again qdist = 1
c
c example: nqdist = d = 4,   d2 = d/2 = 2
c
c   row i:   1   2   3   4   5   6   7   8   ...
c   qdist:   0   0  1/4 2/4 3/4  1   1   1   ...

	use basin

	implicit none

        real qdist(nkn)

c local variables

        integer idist(nkn)

        integer i,k,kk
        integer nqdist,nad
        integer ibc,n,itype,nk
	integer nbc

	integer iapini,ipint
        integer nbnds,itybnd,nkbnds,kbnds
        real getpar

c-----------------------------------------------------------------
c get parameters
c-----------------------------------------------------------------

        do k=1,nkn
          qdist(k) = 1.
          idist(k) = 0
        end do

        nqdist = nint(getpar('nqdist'))		!global value

c-----------------------------------------------------------------
c gather open boundary nodes
c-----------------------------------------------------------------

        n = 0

	nbc = nbnds()

        do ibc=1,nbc
          itype = itybnd(ibc)
          if( itype .eq. 1 .or. itype .eq. 2 ) then
            nk = nkbnds(ibc)
	    call get_bnd_ipar(ibc,'nad',nad)	!local value
	    if( nad .lt. 0 ) nad = nqdist
	    if( nad .gt. 0 ) then
              do i=1,nk
                k = kbnds(ibc,i)
                idist(k) = 1			!old version
                idist(k) = nad
              end do
	    end if
          end if
        end do

c-----------------------------------------------------------------
c make distance
c-----------------------------------------------------------------

        write(6,*) 'Making distance qdist'
        call mkdist_new(nkn,idist,qdist)

c-----------------------------------------------------------------
c write dist (nos) file
c-----------------------------------------------------------------
 
	!old call... please adjourn
        !call wrnos2d('dist','distance from boundary nodes',qdist)

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

        end

c********************************************************************
