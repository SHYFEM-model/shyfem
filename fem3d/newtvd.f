c
c $Id: newtvd.f,v 1.5 2009-04-07 10:43:57 georg Exp $
c
c tvd routines
c
c contents :
c
c subroutine tvd_init(itvd)				initializes tvd scheme
c subroutine tvd_grad_3d(cc,gx,gy,aux,nlvdi,nlv)	computes gradients 3D
c subroutine tvd_grad_2d(cc,gx,gy,aux)			computes gradients 2D
c subroutine tvd_get_upwind_c(ie,l,ic,id,cu,cv)		c of upwind node
c subroutine tvd_upwind_init				init x,y of upwind node
c subroutine tvd_fluxes(ie,l,itot,isum,dt,cl,cv,gxv,gyv,f,fl) tvd fluxes
c
c revision log :
c
c 02.02.2009	ggu&aac	all tvd routines into seperate file
c 24.03.2009	ggu	bug fix: isum -> 6; declaration of cl() was missing 0
c 30.03.2009	ggu	bug fix: ilhv was real in tvd_get_upwind()
c 31.03.2009	ggu	bug fix: do not use internal gradient (undershoot)
c 06.04.2009	ggu&ccf	bug fix: in tvd_fluxes() do not test for conc==cond
c
c*****************************************************************

        subroutine tvd_init(itvd)

c initializes tvd scheme

        implicit none

	integer itvd

	include 'param.h'
	include 'tvd.h'

	integer icall
	save icall
	data icall /0/

	if( icall .ne. 0 ) return

	icall = 1

	itvd_type = itvd

	if( itvd_type .eq. 2 ) call tvd_upwind_init

	write(6,*) 'TVD scheme: ',itvd

	end

c*****************************************************************

        subroutine tvd_grad_3d(cc,gx,gy,aux,nlvdi,nlv)

c computes gradients for scalar cc

        implicit none

	integer nlvdi,nlv
	real cc(nlvdi,1)
	real gx(nlvdi,1)
	real gy(nlvdi,1)
	real aux(nlvdi,1)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        
        integer nen3v(3,1)
        common /nen3v/nen3v        
	integer ilhv(1), ilhkv(1)
	common /ilhv/ilhv, /ilhkv/ilhkv
	include 'ev.h'
        
        integer k,l,ie,ii,lmax
	real b,c,area
	real ggx,ggy

	do k=1,nkn
	  lmax = ilhkv(k)
	  do l=1,lmax
	    gx(l,k) = 0.
	    gy(l,k) = 0.
	    aux(l,k) = 0.
	  end do
	end do

        do ie=1,nel
          area=ev(10,ie) 
	  lmax = ilhv(ie)
	  do l=1,lmax
            ggx=0
            ggy=0
            do ii=1,3
              k=nen3v(ii,ie)
              b=ev(ii+3,ie)
              c=ev(ii+6,ie)
              ggx=ggx+cc(l,k)*b
              ggy=ggy+cc(l,k)*c
              aux(l,k)=aux(l,k)+area
	    end do
            do ii=1,3
             k=nen3v(ii,ie)
             gx(l,k)=gx(l,k)+ggx*area
             gy(l,k)=gy(l,k)+ggy*area
            end do 
          end do
        end do

        do k=1,nkn
	  lmax = ilhkv(k)
	  do l=1,lmax
	    area = aux(l,k)
	    if( area .gt. 0. ) then
	      gx(l,k) = gx(l,k) / area
	      gy(l,k) = gy(l,k) / area
	    end if
	  end do
        end do

        end
        
c*****************************************************************

        subroutine tvd_grad_2d(cc,gx,gy,aux)

c computes gradients for scalar cc (only 2D - used in sedi3d)

        implicit none

	real cc(1)
	real gx(1)
	real gy(1)
	real aux(1)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        
        integer nen3v(3,1)
        common /nen3v/nen3v        
	include 'ev.h'
        
        integer k,ie,ii
	real b,c,area
	real ggx,ggy

	do k=1,nkn
	  gx(k) = 0.
	  gy(k) = 0.
	  aux(k) = 0.
	end do

        do ie=1,nel
          area=ev(10,ie) 
          ggx=0
          ggy=0
          do ii=1,3
              k=nen3v(ii,ie)
              b=ev(ii+3,ie)
              c=ev(ii+6,ie)
              ggx=ggx+cc(k)*b
              ggy=ggy+cc(k)*c
              aux(k)=aux(k)+area
	  end do
          do ii=1,3
             k=nen3v(ii,ie)
             gx(k)=gx(k)+ggx*area
             gy(k)=gy(k)+ggy*area
          end do 
        end do

        do k=1,nkn
	    area = aux(k)
	    if( area .gt. 0. ) then
	      gx(k) = gx(k) / area
	      gy(k) = gy(k) / area
	    end if
        end do

        end
        
c*****************************************************************

        subroutine tvd_get_upwind_c(ie,l,ic,id,cu,cv)

c computes concentration of upwind node

        implicit none

        include 'param.h'
        include 'tvd.h'

        integer ie,l
	integer ic,id
        real cu
        real cv(nlvdim,1)

        integer nen3v(3,1)
        common /nen3v/nen3v        
	integer ilhv(1)			!BUG fix
	common /ilhv/ilhv

        integer ienew
        integer ii,k
        real xu,yu
        real c(3)

        xu = tvdupx(id,ic,ie)
        yu = tvdupy(id,ic,ie)
        ienew = ietvdup(id,ic,ie)

        if( ienew .le. 0 ) return
	if( ilhv(ienew) .lt. l ) return		!TVD for 3D

        do ii=1,3
          k = nen3v(ii,ienew)
          c(ii) = cv(l,k)
        end do

        call femintp(ienew,c,xu,yu,cu)

        end

c*****************************************************************

        subroutine tvd_upwind_init

c initializes position of upwind node

        implicit none

        include 'param.h'
        include 'tvd.h'

        integer ie,ii,j,k
        integer ienew
        real xc,yc,xd,yd,xu,yu

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer nen3v(3,1)
        common /nen3v/nen3v
        real xgv(1)
        common /xgv/xgv
        real ygv(1)
        common /ygv/ygv

        write(6,*) 'setting up tvd upwind information...'

        do ie=1,nel
          do ii=1,3
            k = nen3v(ii,ie)
            xc = xgv(k)
            yc = ygv(k)

            j = mod(ii,3) + 1
            k = nen3v(j,ie)
            xd = xgv(k)
            yd = ygv(k)
            xu = 2*xc - xd
            yu = 2*yc - yd
            call find_elem_from_old(ie,xu,yu,ienew)
            tvdupx(j,ii,ie) = xu
            tvdupy(j,ii,ie) = yu
            ietvdup(j,ii,ie) = ienew

            j = mod(ii+1,3) + 1
            k = nen3v(j,ie)
            xd = xgv(k)
            yd = ygv(k)
            xu = 2*xc - xd
	    yu = 2*yc - yd

	    call find_elem_from_old(ie,xu,yu,ienew)
            tvdupx(j,ii,ie) = xu
            tvdupy(j,ii,ie) = yu
            ietvdup(j,ii,ie) = ienew

            tvdupx(ii,ii,ie) = 0.
            tvdupy(ii,ii,ie) = 0.
            ietvdup(ii,ii,ie) = 0

          end do
        end do

        write(6,*) '...tvd upwind setup done (itvd=2)'

        end

c*****************************************************************

	subroutine tvd_fluxes(ie,l,itot,isum,dt,cl,cv,gxv,gyv,f,fl)

c computes tvd fluxes for one element

	implicit none

	include 'param.h'
        include 'tvd.h'

	integer ie,l
	integer itot,isum
	double precision dt
	double precision cl(0:nlvdim+1,3)		!bug fix
	real cv(nlvdim,nkndim)
        real gxv(nlvdim,nkndim)
        real gyv(nlvdim,nkndim)
	double precision f(3)
	double precision fl(3)

	real eps
	parameter (eps=1.e-8)

        integer nen3v(3,1)
        common /nen3v/nen3v
        real xgv(1)
        common /xgv/xgv
        real ygv(1)
        common /ygv/ygv
        real ulnv(nlvdim,1), vlnv(nlvdim,1)
        common /ulnv/ulnv, /vlnv/vlnv

        logical bgradup
        logical bdebug
	integer ii,k
        integer ic,kc,id,kd,ip
        real term,fact,grad
        real conc,cond,conf,conu
        real gcx,gcy,dx,dy
        real u,v
        real rf,psi
        real alfa,dis
        real vel
        real gdx,gdy

	bgradup = .true.
	bgradup = itvd_type .eq. 2
	bdebug = .true.
	bdebug = .false.

	if( bdebug ) then
	  write(6,*) 'tvd: ',ie,l,itot,isum,dt
	  write(6,*) 'tvd: ',bgradup,itvd_type
	end if

	  do ii=1,3
	    fl(ii) = 0.
	  end do

	  if( itot .lt. 1 .or. itot .gt. 2 ) return

	  u = ulnv(l,ie)
          v = vlnv(l,ie)

            ip = isum
            if( itot .eq. 2 ) ip = 6 - ip		!bug fix

            do ii=1,3
              if( ii .ne. ip ) then
                if( itot .eq. 1 ) then
                  ic = ip
                  id = ii
                  fact = 1.
                else
                  id = ip
                  ic = ii
                  fact = -1.
                end if
                kc = nen3v(ic,ie)
                conc = cl(l,ic)
                kd = nen3v(id,ie)
                cond = cl(l,id)

                dx = xgv(kd) - xgv(kc)
                dy = ygv(kd) - ygv(kc)
                dis = sqrt(dx**2 +dy**2)
                !vel = sqrt(u**2 + v**2)                !total velocity
                vel = abs( u*dx + v*dy ) / dis          !projected velocity
                alfa = ( dt * vel  ) / dis

                if( bgradup ) then
                  conu = cond
                  !conu = 2.*conc - cond		!use internal gradient
                  call tvd_get_upwind_c(ie,l,ic,id,conu,cv)
                  grad = cond - conu
                else
                  gcx = gxv(l,kc)
                  gcy = gyv(l,kc)
                  grad = 2. * (gcx*dx + gcy*dy)
                end if

                if( abs(conc-cond) .lt. eps ) then	!BUG -> eps
                  rf = -1.
                else
                  rf = grad / (cond-conc) - 1.
                end if

                psi = max(0.,min(1.,2.*rf),min(2.,rf))  ! superbee
c               psi = ( rf + abs(rf)) / ( 1 + abs(rf))  ! muscl
c               psi = max(0.,min(2.,rf))                ! osher
c               psi = max(0.,min(1.,rf))                ! minmod

                conf = conc + 0.5*psi*(cond-conc)*(1.-alfa)
                term = fact * conf * f(ii)
                fl(ic) = fl(ic) - term
                fl(id) = fl(id) + term
              end if
            end do

	if( bdebug ) then
	  write(6,*) 'tvd: --------------'
	  write(6,*) 'tvd: ',vel,gcx,gcy,grad
	  write(6,*) 'tvd: ',rf,psi,alfa
	  write(6,*) 'tvd: ',conc,cond,conf
	  write(6,*) 'tvd: ',term,fact
	  write(6,*) 'tvd: ',f
	  write(6,*) 'tvd: ',fl
	  write(6,*) 'tvd: ',(cl(l,ii),ii=1,3)
	  write(6,*) 'tvd: ',ic,id,kc,kd
	  write(6,*) 'tvd: --------------'
	end if

	end

c*****************************************************************

