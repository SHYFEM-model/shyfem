c
c routines for non hydrostatic terms
c
c revision log :
c
c 10.05.2013    dbf     written from scratch
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

	real qpov(nlvdim,nkndim)
	common /qpov/qpov
	real qpnv(nlvdim,nkndim)
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

	real qpov(nlvdim,nkndim)
	common /qpov/qpov
	real qpnv(nlvdim,nkndim)
	common /qpnv/qpnv

	!call prepare_matrix_q

	!call solve_matrix_q

	!call adjust_value_q

c BEGIN test per verificare che il resto sia corretto . Poi cancellare	

	do k = 1,nkn
	  do l=1,nlv
	    qpnv(l,k) = qpov(l,k)
	    if(qpnv(l,k).ne.0.) write(6,*) 'ups..',k,l,qpnv(l,k)
	  end do
	enddo

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

	real qpov(nlvdim,nkndim)
	common /qpov/qpov
	real qpnv(nlvdim,nkndim)
	common /qpnv/qpnv

	do k = 1,nkn
	  do l=1,nlv
	    qpov(l,k) = qpnv(l,k)
	  end do
	enddo

	end

c********************************************************************

	subroutine prepare_matrix_q

	implicit none

	end

c********************************************************************

	subroutine solve_matrix_q

	implicit none

	end

c********************************************************************

	subroutine  adjust_matrix_q

	implicit none

	end

c********************************************************************

	subroutine sp256wnh

c solves for w using prognostic (not diagnostic) equation

	implicit none

c parameters
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
	real qpov(nlvdim,1)
	common /qpov/qpov

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

c start upwind vertical avection wadv=w*dw/dz

	    if (wlnv(l,k) .gt. 0. ) then
             if( l.ge.ilevel )then
	       wadvv = 0
             else
	       dz = hdenv(l+1,ie)
	       wadvv = wlov(l,k)* ( wlov(l,k+1) - wlov(l,k) ) / dz
             endif
            else
             if( l.eq.1 )then
		wadvv = 0.
	     else
	        dz = hdenv(l-1,ie)
	        wadvv = wlov(l,k) * ( wlov(l,k-1) - wlov(l,k) ) / dz
	     endif
            endif	  

c end upwind vertical avection wadv=w*dw/dz

c start computing wp = dq/dz  

	    if( l .eq. ilevel )then
	        dzq = 0.
	    else
                dzq = (hdenv(l,ie) + hdenv(l+1,ie)) / 2
	    endif
	       
            if(l .eq. ilevel ) then
		wq = 0.
	    else
	        wq = (qpov(l,k) - qpov(l+1,k)) / dzq 
	    endif

c end computing wp = dq/dz  

	      wadvh = aj4 * wlov(l,k)* ( c(ii) + b(ii) ) 
	      wdifh = difhv(l,ie) * (  wb * b(ii) + wc * c(ii) ) *aj12
	      wdummy = wlov(l,k) + dt * ( wadvh - wdifh  + wadvv + wq )
  	      wexpl(l,k) = wexpl(l,k) + wdummy	 
	      wlnv(l,k) = wexpl(l,k)
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
     +		    difv(l,k) * (( 1. / ( dzbb * dzcc ))
     +		  + ( 1. / ( dztt * dzcc ))))
	      wdiagt(l,k) = wdiagt(l,k) +  dt * aj4 * (
     +		  - difv(l,k) / ( dztt * dzcc ))
              wdiagb(l,k) = wdiagb(l,k) +  dt * aj4 * (
     +		  - difv(l,k) / ( dzbb * dzcc ))

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
	      wlnv(l,k)=(wlnv(l,k)-wdiagb(l,k)*wlnv(l-1,k))*aux
	  end do
	   
	  lstart = ilevel-1
	  do l=lstart,0,-1
	      wlnv(l,k)=wlnv(l,k)-wlnv(l+1,k)*wdiagt(l,k)
	  end do
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
        integer nlv,nlvdi,lmax
        common /level/ nlvdi,nlv
        integer ilhv(1)
        common /ilhv/ilhv
        integer k,l,ie,ii		
        real fxv(nlvdim,1)
        real fyv(nlvdim,1)
        common /fxv/fxv
        common /fyv/fyv
	integer nen3v(3,1)
	common /nen3v/nen3v
	real qpov(nlvdim,nkndim)
	common /qpov/qpov

	real qxnh,qynh		! non hydrostatic pressure gradient
	real qp
	real b,c

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
            fxv(l,ie) = fxv(l,ie) + qxnh
            fyv(l,ie) = fyv(l,ie) + qynh
	  end do
	end do

	end 

c**********************************************************************

