c
c routines for non hydrostatic terms
c
c revision log :
c
c 10.05.2013    dbf     written from scratch
c 31.05.2013    dbf     written from scratch
c
c********************************************************************

	subroutine nonhydro_init

c initializes non hydrostatic pressure terms

        implicit none

	include 'param.h'
	include 'ev.h'

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nlv,nlvdi,lmax
	common /level/ nlvdi,nlv
	integer k,l

	real*8 qpov(nlvdim,nkndim)
	common /qpov/qpov
	real*8 qpnv(nlvdim,nkndim)
	common /qpnv/qpnv
        save /qpov/,/qpnv/	!this should be enough for all appearances 

	do k = 1,nkn
	   do l=1,nlv
	    qpov(l,k) = 0. 
	    qpnv(l,k) = 0.
	   end do
	enddo

	end

c********************************************************************

	subroutine nonhydro_adjust

c integrates non hydrostatic adjustment to equations

	implicit none

	include 'param.h'
	include 'ev.h'

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nlv,nlvdi,lmax
	common /level/ nlvdi,nlv
	integer k,l

	real*8 qpov(nlvdim,nkndim)
	common /qpov/qpov
	real*8 qpnv(nlvdim,nkndim)
	common /qpnv/qpnv
	integer icall_conh
	save icall_conh
	data icall_conh /0/

	if (icall_conh.eq.0) then          ! only first time
	write(6,*)'system_init_nh3dmatrix'
	call system_init_nh3dmatrix
	write(6,*)'end system_init_nh3dmatrix'
	icall_conh = 1
	endif

	call nonhydro_prepare_matrix

	call nonhydro_solve_matrix
	!call nonhydro_solve_matrix_pard

	call nonhydro_adjust_value

	call nonhydro_correct_uveta

c BEGIN test per verificare che il resto sia corretto . Poi cancellare	

c	do k = 1,nkn
c	   do l=1,nlv
c	    qpnv(l,k) = qpov(l,k) 
c	   if(qpnv(l,k).ne.0.)write(6,*)'ups..',k,l,qpnv(l,k)
c	   end do
c	enddo

c END test per verificare che il resto sia corretto . Poi cancellare	

	end

c********************************************************************

	subroutine nonhydro_copy

c copies new values of q to old time step

	implicit none
	
	include 'param.h'
	include 'ev.h'

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nlv,nlvdi,lmax
	common /level/ nlvdi,nlv
	integer k,l

	real*8 qpov(nlvdim,nkndim)
	common /qpov/qpov
	real*8 qpnv(nlvdim,nkndim)
	common /qpnv/qpnv

	do k = 1,nkn
	   do l=1,nlv
	    qpov(l,k) = qpnv(l,k) 
	   end do
	enddo

	end

c********************************************************************

	subroutine sp256wnh 

c solves for w using prognostic (not diagnostic) equation

	implicit none

cc parameters
        include 'param.h'
	include 'ev.h'

c common 
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nlvdi,nlv
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /level/ nlvdi,nlv
	integer itanf,itend,idt,nits,niter,it
	common /femtim/itanf,itend,idt,nits,niter,it
        real grav,fcor,dcor,dirn,rowass,roluft
        common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft
        real rrho0	
	real difv(0:nlvdim,1)
	common /difv/difv
	real difhv(nlvdim,1)
	common /difhv/difhv
        real uprv(nlvdim,1), vprv(nlvdim,1)
	common /uprv/uprv, /vprv/vprv
        real wlov(0:nlvdim,1),wlnv(0:nlvdim,1)
	common /wlov/wlov, /wlnv/wlnv
        integer ilhv(1), ilhkv(1)
	common /ilhv/ilhv, /ilhkv/ilhkv
	real hdknv(nlvdim,nkndim)
	common /hdknv/hdknv
	real hdkov(nlvdim,1)
	common /hdkov/hdkov
	real hdenv(nlvdim,1)
	common /hdenv/hdenv
	real hdeov(nlvdim,1)
	common /hdeov/hdeov
        real areakv(nlvdim,1)
	common /areakv/areakv
	integer nen3v(3,1)
	common /nen3v/nen3v
	real zeov(3,1),zenv(3,1)
	common /zeov/zeov, /zenv/zenv
	real*8 qpov(nlvdim,nkndim)
	common /qpov/qpov
	real*8 qpnv(nlvdim,nkndim)
	common /qpnv/qpnv
        integer k,ie,ii,l,iii,ll
	integer ilevel,lmax
	integer lstart
	integer kn(3)
	double precision aj,rk3,rv,aj4,aj12
	double precision b(3),c(3),f(3)
	real dt,w,aux
	real dz, dzq
	real dzcc,dzbb,dztt
	real ddt,wdummy,wadvv,wadvh,wdifh,wq
	real wb,wc
	real wexpl(0:nlvdim,nkndim)
	real wdiag(0:nlvdim,nkndim)
	real wdiagt(0:nlvdim,nkndim),wdiagb(0:nlvdim,nkndim)

	call get_timestep(dt)

	rrho0 = 1./rowass
	write(6,*)'rhow',rrho0,rowass
c resolution of the explicit part
	
	do k = 1,nkn
	  wlnv(0,k) = 0.
	  lmax = ilhkv(k)
	  wlnv(lmax,k) = 0.
	  do l = 0,lmax
		wlnv(l,k) = 0.
	    	wexpl(l,k) = 0.
		wdiagt(l,k) = 0.
		wdiag(l,k) = 0.
		wdiagb(l,k) = 0.
	  enddo
	enddo

	do ie=1,nel
	  do ii=1,3
            k=nen3v(ii,ie)
	    kn(ii)=k
	    b(ii)=ev(ii+3,ie)
	    c(ii)=ev(ii+6,ie)
	  end do
          
	  aj=ev(10,ie) !area of triangle / 12
	  aj4=4.*aj
	  aj12=12.*aj
	  ilevel=ilhv(ie)
	  do l=1,ilevel
                wb = 0.
	        wc = 0.
	        do ii=1,3
	   		k=nen3v(ii,ie)
	   		wb = wb + ( b(ii) * wlov(l,k) )
	   		wc = wc + ( c(ii) * wlov(l,k) )
	   	enddo
	   
	   do ii=1,3
	    k=nen3v(ii,ie)

c start upwind vertical advection wadv=w*dw/dz

	    if (wlnv(l,k) .gt. 0. ) then
             if( l.ge.ilevel )then
	       wadvv = 0
             else
	       dz = hdenv(l+1,ie)
	       wadvv = wlov(l,k)* ( wlov(l+1,k) - wlov(l,k) ) / dz
!DEB NUOVO
             endif
            else
             if( l.eq.1 )then
		wadvv = 0.
	     else
	        dz = hdenv(l,ie)
	        wadvv = wlov(l,k) * ( wlov(l-1,k) - wlov(l,k) ) / dz
!DEB NUOVO
	     endif
            endif	  

c end upwind vertical advection wadv=w*dw/dz

c start computing wp = dq/dz  

	    if( l .eq. ilevel )then
	        dzq = 0.
	    else
                dzq = (hdenv(l,ie) + hdenv(l+1,ie)) / 2
	    endif
	       
            if(l .eq. ilevel ) then
		wq = 0.
	    else
	        wq = rrho0 * (qpov(l+1,k) - qpov(l,k)) / dzq 
	    endif

c end computing wp = dq/dz  

	      !wadvh = aj4 * wlov(l,k)* ( c(ii) + b(ii) ) 
	      wadvh = aj4 * wlov(l,k)*
     +			((uprv(l,k)* b(ii))+(vprv(l,k) * c(ii)))
!DEB NEW
	      wdifh = difhv(l,ie) * (  wb * b(ii) + wc * c(ii) ) *aj12
	      wdummy = wlov(l,k) + dt * ( -wadvh + wdifh  - wadvv + wq )
  	      wexpl(l,k) = wexpl(l,k) + wdummy	 
	      wlnv(l,k) = wexpl(l,k)
	      if(wlnv(l,k).gt.0.)then
	      write(688,*)l,k,wlnv(l,k),wexpl(l,k),wdummy,wdifh,wadvh,wq,
     +		wadvv
              endif
	   end do
	  end do

c----------------------------------------------------------------
c       set up vectors for use in assembling contributions
c----------------------------------------------------------------
          
	 do l=1,ilevel-1

	   dzcc = ( hdenv(l,ie) / 2 ) + ( hdenv(l+1,ie) / 2 ) 
           dztt = hdenv(l,ie)
	   dzbb = hdenv(l+1,ie)

	   do ii=1,3
              k=nen3v(ii,ie)

       	      wdiag(l,k) = wdiag(l,k) + aj4 +  dt * aj4 *(  
     +		  - difv(l,k) * (( 1. / ( dzbb * dzcc )) 
     +		  + ( 1. / ( dztt * dzcc ))))
	      wdiagt(l,k) = wdiagt(l,k) +  dt * aj4 * (
     +		  + difv(l,k) / ( dztt * dzcc ))
              wdiagb(l,k) = wdiagb(l,k) +  dt * aj4 * (
     +		  + difv(l,k) / ( dzbb * dzcc ))

	   enddo
	     
	 enddo

	end do

c a = diag b=diagb c=diagt
c  0        	    1 0 0 0 
c  1	  	    0 a c 0  
c                   . . . . 
c  ilevel-1	    0 b a 0  
c  ilevel	    0 0 b 1  

        do k=1,nkn
	  ilevel = ilhkv(k)
          wlnv(0,k) = 0.
          wlnv(ilevel,k) = 0.
	  wdiag(0,k) = 1.
	  wdiag(ilevel,k) = 1.
	  aux=1./wdiag(1,k)
	  wdiagt(0,k)=wdiagt(0,k)*aux
	  wlnv(0,k)=wlnv(0,k)*aux

	  do l=1,ilevel
	      aux=1./(wdiag(l,k)-wdiagb(l,k)*wdiagt(l-1,k))
	      wdiagt(l,k)=wdiagt(l,k)*aux
	      if(l.gt.1)then!DEB NUOVO
	      wlnv(l,k)=(wlnv(l,k)-wdiagb(l,k)*wlnv(l-1,k))*aux
              else !DEB NUOVO
	      wlnv(l,k)=(wlnv(l,k))*aux!DEB NUOVO
              endif !DEB NUOVO
	  end do
	   
	  lstart = ilevel-1
	  do l=lstart,0,-1
	      if(l.ne.-1)then !DEB NUOVO
	              wlnv(l,k)=wlnv(l,k)-wlnv(l+1,k)*wdiagt(l,k)
              else !DEB NUOVO
        	      wlnv(l,k)=wlnv(l,k)!DEB NUOVO
              endif !DEB NUOVO
c	  if(wlnv(l,k).gt.1.)write(6,*)'wlnv',l,k,wlnv(l,k)
	  end do
          wlnv(0,k) = 0. !DEB NUOVO
          wlnv(ilevel,k) = 0. !DEB NUOVO
	  wdiag(0,k) = 1. !DEB NUOVO
	  wdiag(ilevel,k) = 1. !DEB NUOVO
        end do

c end resolution of the implicit part

	end

c******************************************************************

	subroutine nonhydro_set_explicit 

c adds explicit part of non hydrostatic pressure to explict terms

	implicit none
        
	include 'param.h'
	include 'ev.h'
        
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real grav,fcor,dcor,dirn,rowass,roluft
        common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft
        integer nlv,nlvdi
        common /level/ nlvdi,nlv

        integer ilhv(1)
        common /ilhv/ilhv
        real fxv(nlvdim,1)      
        real fyv(nlvdim,1)
        common /fxv/fxv
        common /fyv/fyv
	integer nen3v(3,1)
	common /nen3v/nen3v
        real hdeov(nlvdim,1)
        common /hdeov/hdeov
	real*8 qpov(nlvdim,nkndim)
	common /qpov/qpov

        integer k,l,ie,ii,lmax	
	real hhi
        real rrho0	
	double precision aj,rk3,rv,aj4,aj12

	real qxnh,qynh ! non hydrostatic pressure gradient !DEB NH
	real qp
	real b,c
	real aq ,aqt

	real getpar

	aq = getpar('aqpar')
	aqt = 1. - aq

	rrho0=1./rowass

	do ie=1,nel
	  lmax = ilhv(ie)
          do l=1,lmax
	    qxnh = 0.
	    qynh = 0.
            do ii=1,3
                k = nen3v(ii,ie)
	        qp = qpov(l,k)
                b = ev(3+ii,ie)
                c = ev(6+ii,ie)
		qxnh = qxnh + b * qp
		qynh = qynh + c * qp
	    end do
	    qxnh = qxnh / 3.
	    qynh = qynh / 3.

	    hhi = hdeov(l,ie)  

            fxv(l,ie) = fxv(l,ie) + aqt * hhi * rrho0 * qxnh
            fyv(l,ie) = fyv(l,ie) + aqt * hhi * rrho0 * qynh
	  end do
	end do

	end 

c**********************************************************************


	subroutine nonhydro_prepare_matrix

	implicit none

cc parameters
        include 'param.h'
	include 'ev.h'

c common 
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nlvdi,nlv
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /level/ nlvdi,nlv
	integer itanf,itend,idt,nits,niter,it
	common /femtim/itanf,itend,idt,nits,niter,it
	real difv(0:nlvdim,1)
	common /difv/difv
	real difhv(nlvdim,1)
	common /difhv/difhv
        real uprv(nlvdim,1), vprv(nlvdim,1)
	common /uprv/uprv, /vprv/vprv
        real wlov(0:nlvdim,1),wlnv(0:nlvdim,1)
	common /wlov/wlov, /wlnv/wlnv
        integer ilhv(1), ilhkv(1)
	common /ilhv/ilhv, /ilhkv/ilhkv
	real hdknv(nlvdim,nkndim)
	common /hdknv/hdknv
	real hdkov(nlvdim,1)
	common /hdkov/hdkov
	real hdenv(nlvdim,1)
	common /hdenv/hdenv
	real hdeov(nlvdim,1)
	common /hdeov/hdeov
        real areakv(nlvdim,1)
	common /areakv/areakv
	integer nen3v(3,1)
	common /nen3v/nen3v
	real zeov(3,1),zenv(3,1)
	common /zeov/zeov, /zenv/zenv
	real*8 qpov(nlvdim,nkndim)
	common /qpov/qpov
	real*8 qpnv(nlvdim,nkndim)
	common /qpnv/qpnv

        integer k,ie,ii,l,iii,k2,k3,ii2,ii3,jj
	integer ilevel,lmax
	integer lstart
	integer kn(3)
	double precision aj,rk3,rv,aj4,aj12
	double precision b(3),c(3),f(3)
	real dt,w,aux
	real dz, dzq
	real dzcc,dzbb,dztt,dzcct,dztt1
	real ddt,wdummy,wadvv,wadvh,wdifh,wq
	real wb,wc
	real wexpl(0:nlvdim,nkndim)
	real wdiag(0:nlvdim,nkndim)
	real wdiagt(0:nlvdim,nkndim),wdiagb(0:nlvdim,nkndim)
	integer nl,n,m,ll,kk
	real*8 mat(3,3),matl(3,3)
	real*8 vnot(3)
	integer ia,ib,ic,iat,ibt,ict,iab,ibb,icb
	real aq
	 real getpar

	 aq = getpar('aqpar')
	
	call get_timestep(dt)


	do n=1,3
	  do m=1,3
	        mat(n,m)=0.
	  enddo
		vnot(n)=0. !DEB NUOVO
	enddo
	      

	do ie=1,nel
	  lmax = ilhv(ie)
	  aj=ev(10,ie) !area of triangle / 12
	  aj4=4.*aj
	  aj12=12.*aj
          do l=1,lmax
           dztt = hdenv(l,ie)
	   dzbb = hdenv(l+1,ie)
	   dzcc = ( hdenv(l,ie) / 2. ) + ( hdenv(l+1,ie) / 2. ) 
	   dzcct = ( hdenv(l-1,ie) / 2. ) + ( hdenv(l,ie) / 2. ) !NEW
	   dztt1 =  hdenv(l-1,ie) !NEW

           do ii=1,3
                kn(ii) = nen3v(ii,ie)
                b(ii)=ev(ii+3,ie)
		c(ii)=ev(ii+6,ie)
	   enddo
           do ii=1,3
                do jj=1,3
		wlnv(0,kn(jj))=0.
		wlnv(0,kn(ii))=0.
		wlnv(lmax,kn(jj))=0.
		wlnv(lmax,kn(ii))=0.

		mat(ii,jj)=+aj12*dt*dztt*(b(ii)*b(jj)+c(ii)*c(jj))*aq

		if(l.eq.1)then
	    		vnot(jj)= - aj4 * dztt * (uprv(l,kn(jj)) * b(jj) + 
     +			vprv(l,kn(jj)) * c(jj))
     +			- aj4 * 2* wlnv(l,kn(jj))  
		else
 	     		vnot(jj)= - aj4 * dztt * ( uprv(l,kn(jj)) * b(jj) + 
     +			vprv(l,kn(jj)) * c(jj))
     +			- aj4 * (wlnv(l,kn(jj))-wlnv(l-1,kn(jj)))    
		endif
		
		enddo          
		if(l.gt.1.and.l.lt.lmax)then
	         !layer center
     		 matl(ii,1) = aj4 * dt *(-(1./dzcct)+(1./dzcc)) *aq 
                 !layer top
	         matl(ii,2) = aj4 * dt * (1./dzcct)* aq 
                 !layer bottom
	         matl(ii,3) =  - aj4 * dt *(1./dzcc) * aq 
	       elseif(l.eq.1)then
     		 matl(ii,1) = aj4 * dt * (1./dzcc) * aq 
	         matl(ii,2) = 0. 
	         matl(ii,3) = - aj4 * dt * (1./dzcc) * aq 
	       else 
	         !layer center
     		 matl(ii,1) = aj4 * dt *(-1./dzcct) * aq 
                 !layer top
	         matl(ii,2) = aj4 * dt * (1./dzcct) * aq 	
                 !layer bottom
	         matl(ii,3) = 0.
		endif

	        if(mat(ii,ii).eq.0.)then
                    write(6,*)'mat',mat(ii,ii),ii,jj,aj12,dt
                    write(6,*) b(ii),c(ii)
		    stop
                endif
	  enddo
	  call system_assemble_nh_3d(ie,kn,l,mat,matl,vnot)
	 enddo
	enddo

	end

cc********************************************************************

         subroutine  nonhydro_adjust_value

         implicit none

	 include 'param.h'
	 include 'nohydlinks.h'
         integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	 common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
         integer nlvdi,nlv
	 common /level/ nlvdi,nlv
         integer k,l,lmax,nn
         integer ilhv(1), ilhkv(1)
         common /ilhv/ilhv, /ilhkv/ilhkv
         real*8 qpnv(nlvdim,nkndim)
         common /qpnv/qpnv

         nn=0
         do k=1,nkn
           lmax=ilhkv(k)
           do l=1,lmax
	     nn=nn+1
	     qpnv(l,k)=rvecnh(nn)
	     !if(l.eq.1.or.l.eq.lmax)qpnv(l,k)=0. !DEB NUOVO
	     !if(qpnv(l,k).ne.0.)then
	     !write(6,*)qpnv(l,k),'alloora'
             !endif
	   enddo
         enddo
	 
	 end

c********************************************************************
	 subroutine nonhydro_correct_uveta

	 implicit none

	 include 'param.h'
	 include 'ev.h'

         integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	 common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
         integer nlvdi,nlv
	 common /level/ nlvdi,nlv
         integer ilhv(1), ilhkv(1)
         common /ilhv/ilhv, /ilhkv/ilhkv
	 integer itanf,itend,idt,nits,niter,it
	 common /femtim/itanf,itend,idt,nits,niter,it
         integer nen3v(3,1)
	 common /nen3v/nen3v
         real*8 qpnv(nlvdim,nkndim)
         common /qpnv/qpnv
	 real hdknv(nlvdim,nkndim)
	 common /hdknv/hdknv
         real zov(1),znv(1),unv(1),vnv(1)
         real ulnv(nlvdim,1),vlnv(nlvdim,1)
         real wlov(0:nlvdim,1),wlnv(0:nlvdim,1)
	 common /wlov/wlov, /wlnv/wlnv
         real uprv(nlvdim,1), vprv(nlvdim,1)
	 common /uprv/uprv, /vprv/vprv
	 common /zov/zov, /znv/znv
	 common /ulnv/ulnv, /vlnv/vlnv
	 integer kn(3)
	 double precision aj,rk3,rv,aj4,aj12
	 double precision b(3),c(3),f(3)
	 real dt,dzcc
	 real*8 uqaux,vqaux
         real grav,fcor,dcor,dirn,rowass,roluft
         common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft
         real rrho0	
	 real getpar
	 real aq 
	 integer ie,l,ilevel,ii,k,lmax

	 aq = getpar('aqpar')

	 call get_timestep(dt)

	 rrho0=1./rowass

	 do ie=1,nel
	   do ii=1,3
            k=nen3v(ii,ie)
	    kn(ii)=k
	    b(ii)=ev(ii+3,ie)
	    c(ii)=ev(ii+6,ie)
	   end do
	  aj=ev(10,ie) !area of triangle / 12
	  aj4=4.*aj
	  aj12=12.*aj
	  ilevel=ilhv(ie)

          do l=1,ilevel
	  uqaux=0. !DEB NUOVO
	  vqaux=0. !DEB NUOVO
	   do ii=1,3
	      uqaux = uqaux + ( b(ii) *qpnv(l,kn(ii)))
	      vqaux = vqaux + (c(ii) *qpnv(l,kn(ii)))
	    !if(uqaux.ne.0.)write(6,*)'qpnv(l,kn(ii)))',qpnv(l,kn(ii))

	   enddo
	    !  if(uqaux.ne.0.)write(6,*)'uqaux',uqaux,'dt',dt
	   ulnv(l,ie) = ulnv(l,ie) - aq * dt * rrho0 * aj4 * uqaux
	   vlnv(l,ie) = vlnv(l,ie) - aq * dt * rrho0 * aj4 * vqaux
	   !write(6,*)'u,v',l,ie,ulnv(l,ie),vlnv(l,ie),aq,dt,uqaux,vqaux   
	  enddo
	 enddo

	 do k = 1,nkn
	      znv(k) = znv(k) + rrho0 * qpnv(1,k) / grav 
	      if(k.eq.1793)write(6,*)k,znv(k),dt
	    !  write(6,*)znv(k),rrho0 ,qpnv(1,k),grav
c new part
c	      lmax = ilhkv(k)
c	      do l=0,lmax
c              dzcc = (hdknv(l+1,k)/2.)+(hdknv(l,k)/2.)
c	      wlnv(l,k)=wlnv(l,k) - (aq * dt * rrho0 * 
c     +	 (qpnv(l+1,k)-qpnv(l,k))/dzcc)
c	      if(qpnv(l,k).ne.0.)then
c      !write(6,*)wlnv(l,k),aq,dt,qpnv(l+1,k),qpnv(l+1,k),dzcc
c              endif

c	      enddo
c end new part
	 enddo
	 end
c********************************************************************

	 subroutine system_assemble_nh_3d(ie,kn,l,mat,matl,vnot)
	
         implicit none

	 include 'param.h'
	 include 'links.h'
	 include 'nohydlinks.h'
         integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	 common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
         integer nlvdi,nlv
	 common /level/ nlvdi,nlv
         integer ilhv(1), ilhkv(1)
	 common /ilhv/ilhv, /ilhkv/ilhkv
         integer nen3v(3,1)
	 common /nen3v/nen3v

	 integer i,j,l,ipkk,ipkkt,ipkkb,ik,nk
	 integer kn(3)
	 real*8 mat(3,3),matl(3,3)
	 real*8 vnot(3) 
	 integer lmax
	 integer ie,nn

	 nn=0
	 do i=1,3
	  !lmax=ilhkv(kn(i))
	  lmax=ilhv(ie)
	  do j=1,3
	   nn=nn+1
	   ipkk=iopp(ie,nn,l)
	   ipkkt=ipkk - 1 !position for contribution on (k,l-1)
	   ipkkb=ipkk + 1 !position for contribution on (k,l+1)
	
	  if(ioii(ipkk).ne.kn(i)) then
	   write(6,*),'error stop ioii'
	   write(6,*)ioii(ipkk),kn(i),l,lmax,ipkk
	   stop
	  endif
	  if(iojj(ipkk).ne.kn(j)) then
	   write(6,*),'error stop iojj'
	   write(6,*)iojj(ipkk),kn(j),l,lmax,ipkk
	   stop
	  endif

	  if(ipkk.ne.0)then
	   if (l.ne.1.and.l.ne.lmax)then
	   if(ipkk.gt.0) conh(ipkk) = conh(ipkk) + mat(i,j) + matl(i,2)
	   if(ipkkt.gt.0) conh(ipkkt) = conh(ipkkt) + matl(i,1)
	   if(ipkkb.gt.0) conh(ipkkb) = conh(ipkkb) + matl(i,3)
	   elseif(l.eq.1)then
	   if(ipkk.gt.0) conh(ipkk) = conh(ipkk) + mat(i,j) + matl(i,2)
	   if(ipkkt.gt.0) conh(ipkkt) = conh(ipkkt) 
	   if(ipkkb.gt.0) conh(ipkkb) = conh(ipkkb) + matl(i,3)
	   elseif(l.eq.lmax)then
	   if(ipkk.gt.0) conh(ipkk) = conh(ipkk) + mat(i,j) + matl(i,2)
	   if(ipkkt.gt.0) conh(ipkkt) = conh(ipkkt) + matl(i,1)
	   if(ipkkb.gt.0) conh(ipkkb) = conh(ipkkb) 
	   else
	    write(6,*)'ERROR stop, ipkk pointer is zero'
            write(6,*) 'ie,kn(i),kn(j),l',ie,kn(i),kn(j),l
	    stop
	   endif
	  endif

	 if(ioii1(ipkk).eq.iojj1(ipkk).and.conh(ipkk).eq.0)then
	   write(6,*)'ERROR DIAG = 0',conh(ipkk),ipkk
           write(6,*)l,lmax,ioii1(ipkk),iojj1(ipkk)
           write(6,*)matl(i,1),matl(i,2),matl(i,3)
	   write(6,*)mat(i,j)
	   stop
	  endif

	  !vecnot(ipkk)=vecnot(ipkk) + vnot(j) !DEB NEW
	  !rvecnh(ipkk)=rvecnh(ipkk) + vnot(j) !DEB NEW
	  !rvecnh(iojj1(ipkk))=rvecnh(iojj1(ipkk)) + vnot(j) !DEB NEW
	  !rvecnh(ioii1(ipkk))=rvecnh(ioii1(ipkk)) + vnot(j) !DEB NEW
	  !write(6,*)'vecnot',ipkk,vecnot(ipkk),conh(ipkk)
	  !write(6,*)'vecnot',ipkk,rvecnh(iojj1(ipkk)),vnot(j),conh(ipkk)
	  enddo
	  rvecnh(iojj1(ipkk))=rvecnh(iojj1(ipkk)) + vnot(i) !DEB NEW
	  !vecnot(ipkk)=vecnot(ipkk) + vnot(i)
	 enddo
	 
	 end

c********************************************************************

	 subroutine system_init_nh3dmatrix

c construct pointers for matrix

         implicit none

	 include 'param.h'
	 include 'links.h'
	 include 'nohydlinks.h'
         integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	 common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
         integer nlvdi,nlv
	 common /level/ nlvdi,nlv
         integer ilhv(1), ilhkv(1)
         common /ilhv/ilhv, /ilhkv/ilhkv
         integer nen3v(3,1)
	 common /nen3v/nen3v
 	 integer kn(3)	 
	 integer kk,nk,i,k,l,ll,mm,ik,ip,n,nn,lmax
	 integer ii,jj,ie,nkmax(nlvdim,nkndim),ss,nnn
	 integer idiag(nlvdim,nkndim),lmax1,itest
	 logical debug

	 debug = .true.	

	 write(6,*) 'SOLVER: Sparskit for non hydro'
         n=0
         do k=1,nkn
          lmax=ilhkv(k)
          do l=1,lmax
	   n=n+1
           rvecnh(n) = 0.
	   rauxnh(n) = 0.
	  enddo
         end do

	 n=0
	 ipl(0) = 0
	 kpl(0) = 0 
	 lpl(0) = 0   
	 do k=1,nkn
	   lmax=ilhkv(k)
	   n=n+lmax
	   ipl(k)=n !aux sum of lmax for each node
	   iopup(0,k)=0
           do l = 1,nlv
	      idiag(l,k)=0
	      nkmax(l,k)=0
	      do nn=1,ngr+3
	       iop(nn,l,k)= 0
	       ikop(nn,l,k)= 0
	       ilop(nn,l,k)= 0
	      enddo
	   enddo
	 enddo

         nnn=0
	 do ie=1,nel
	   do nn=1,9
              do l = 1,nlv
	       nnn=nnn+1
	       iopp(ie,nn,l)=0 
	       ioii(nnn)=0
	       iojj(nnn)=0
	       ioii1(nnn)=0
	       iojj1(nnn)=0
	       conh(nnn)=0.
	       vecnot(nnn)=0.
	      enddo
	   enddo
	 enddo
	 
	 matdimmax=n
	 write(6,*)'Max dimension matrix non hydro = ',matdimmax
	 
	 mm = 0
	 do k=1,nkn

	  lmax=ilhkv(k)
	  call set_node_links(k,nk)

	  do l=1,lmax
	     do i=1,matdimmax
	        kpl(i) = 0 
	        lpl(i) = 0 
	        kpl1(i) = 0 
	        lpl1(i) = 0 
	        rownh(i) = 0
             enddo

             do kk=1,nk
	       ik = lnk_nodes(kk)
	       ip = ipl(ik-1) + l  !ip is index in array matdimmax
	       kpl(ip) = ik
	       lpl(ip) = l
	       rownh(ip)= 1
	     enddo

	     ip =  ipl(k-1) + l
	     kpl(ip) = k
	     lpl(ip) = l
	     rownh(ip)= 1 

	     if( lpl(ip) .eq. 1 )then
	     	  mm = mm + nk + 1 + 1
	     elseif( lpl(ip) .eq. lmax )then
	          mm = mm + nk + 1 + 1
	     else
	          mm = mm + nk + 1 + 2 ! max number of non zero values at the row (k,l)
	     endif

	     iopup(l,k)= mm

	   
	     nn=0
	     do i=1,matdimmax
	           if(rownh(i).eq.1)then
		     nn=nn+1
	  	     kpl1(nn) = kpl(i)
	             lpl1(nn) = lpl(i)
	           endif       
	     enddo

	     nkmax(l,k)=nn
	     call sorti(nn,kpl1)

	     do i=1,nkmax(l,k)
	          if(kpl1(i) .eq. k .and. lpl1(i) .eq.l)then
	               iop(i,l,k)=0
		       idiag(l,k) = i
	          endif
	     enddo

	     do i=1,nkmax(l,k)
	         if(i.lt.idiag(l,k))then
		    do ii=i,idiag(l,k)-1	
	              if(l.eq.1)then	    
	                iop(ii,l,k)= -idiag(l,k) + ii 
		      else
	                iop(ii,l,k)= -idiag(l,k) + ii - 1
		      endif
		    enddo
	         elseif(i.gt.idiag(l,k))then
	            do ii=i,nkmax(l,k)
	              if(l.eq.lmax)then	    
	                iop(ii,l,k)= -idiag(l,k) + ii 
	              else
	                iop(ii,l,k)= -idiag(l,k) + ii + 1
		      endif
		    enddo
	         endif

	     ikop(i,l,k) = kpl1(i)
	     ilop(i,l,k) = lpl1(i) 

	     enddo
          enddo
	 enddo

	 nnn=0
	 do ie = 1,nel
	   do ii = 1,3
               kn(ii) = nen3v(ii,ie)
	   enddo 
	   lmax = ilhv(ie)
	   do l=1,lmax
	      nn = 0
	      do ii=1,3
	        do jj=1,3
		    nn=nn+1
		  do ss=1,nkmax(l,kn(ii))

c to set position in the diagonal and on bands on the same layer (krow,l) (kcol,l)

	            if(kn(jj).eq.ikop(ss,l,kn(ii)))then
	             if(l.eq.1)then
		         if(kn(ii).eq.1)then
	                 iopp(ie,nn,l)= idiag(l,kn(ii)) + iop(ss,l,kn(ii)) 
		         else
	                 iopp(ie,nn,l)= iopup(ilhkv(kn(ii)-1),kn(ii)-1) + 
     +		         idiag(l,kn(ii)) + iop(ss,l,kn(ii))
		         endif 
	             else
	                 iopp(ie,nn,l)= iopup(l-1,kn(ii)) + 
     +		                  idiag(l,kn(ii)) + iop(ss,l,kn(ii)) +1
	             endif
		    endif
		  enddo 
	            if(iopp(ie,nn,l).ne.0)then
			   nnn=nnn+1
		    endif 
	            ioii(iopp(ie,nn,l))=kn(ii)	    
	            iojj(iopp(ie,nn,l))=kn(jj)	    
	            ioll(iopp(ie,nn,l))=l	    
	            ioii1(iopp(ie,nn,l))=ipl(kn(ii)-1) + l
	            iojj1(iopp(ie,nn,l))=ipl(kn(jj)-1) + l
	        enddo
	      enddo
	   if(iopp(ie,nn,l).eq.0)write(6,*)'ERRi iopp',ie,nn,l,lmax
	   enddo
	 enddo
	nn=0
        nnzeronh = nnn
	write(6,*)'Number of non zero non hydro =',nnzeronh


	if(debug)then
	do ie=1,nel
          do nn=1,9
	   do l=1,ilhv(ie)
	   if(iopp(ie,nn,l).eq.0)write(6,*)'ERRi iopp',ie,nn,l,ilhv(ie)
           enddo
	  enddo
	enddo
	endif

	 end

c******************************************************************

       subroutine sorti(n,ra)

c heapsort !DEB

      implicit none

      integer n
      integer ra(n)

      integer i,ir,j,l
      integer rra

      if(n.lt.2) return

      l=n/2+1
      ir=n
10    continue
        if(l.gt.1)then
          l=l-1
          rra=ra(l)
        else
          rra=ra(ir)
          ra(ir)=ra(1)
          ir=ir-1
          if(ir.eq.1)then
            ra(1)=rra
            return
          endif
        endif
        i=l
        j=l+l
20      continue
	if(j.le.ir)then
          if(j.lt.ir)then
            if(ra(j).lt.ra(j+1))j=j+1
          endif
          if(rra.lt.ra(j))then
            ra(i)=ra(j)
            i=j
            j=j+j
          else
            j=ir+1
          endif
        go to 20
        endif
        ra(i)=rra
      go to 10
      end

c********************************************************************
