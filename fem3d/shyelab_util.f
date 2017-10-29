
!***************************************************************

        subroutine prepare_hydro(bvel,nndim,cv3all,znv,uprv,vprv)

        use basin
        use levels
        use mod_depth

        implicit none

        logical bvel
        integer nndim
        real cv3all(nlvdi,nndim,0:4)
        real znv(nkn)
        real uprv(nlvdi,nkn)
        real vprv(nlvdi,nkn)

        real, allocatable :: zenv(:)
        real, allocatable :: uv(:,:)
        real, allocatable :: vv(:,:)

        allocate(zenv(3*nel))
        allocate(uv(nlvdi,nel))
        allocate(vv(nlvdi,nel))

        znv(1:nkn)     = cv3all(1,1:nkn,1)
        zenv(1:3*nel)  = cv3all(1,1:3*nel,2)
        uv(:,1:nel)    = cv3all(:,1:nel,3)
        vv(:,1:nel)    = cv3all(:,1:nel,4)

        call shy_transp2vel(bvel,nel,nkn,nlv,nlvdi,hev,zenv,nen3v
     +                          ,ilhv,hlv,uv,vv
     +                          ,uprv,vprv)

        deallocate(zenv,uv,vv)

        end

!***************************************************************

        subroutine convert_to_speed(uprv,vprv,sv,dv)

        use basin
        use levels

        implicit none

        real uprv(nlvdi,nkn)
        real vprv(nlvdi,nkn)
        real sv(nlvdi,nkn)
        real dv(nlvdi,nkn)

        integer k,lmax,l
        real u,v,s,d

        do k=1,nkn
          lmax = ilhkv(k)
          do l=1,lmax
            u = uprv(l,k)
            v = vprv(l,k)
            call c2p_ocean(u,v,s,d)   !d is ocean convention
            sv(l,k) = s
            dv(l,k) = d
          end do
        end do

        end

!***************************************************************

        subroutine shy_transp2vel(bvel,nel,nkn,nlv,nlvddi
     +				,hev,zenv,nen3v
     +                          ,ilhv,hlv,utlnv,vtlnv
     +                          ,uprv,vprv)

c transforms transports at elements to velocities at nodes

        implicit none

	logical bvel			!if true compute velocities
        integer nel
        integer nkn
        integer nlv
        integer nlvddi
        real hev(nel)
        real zenv(3,nel)
        integer nen3v(3,nel)
        integer ilhv(nel)
        real hlv(nlvddi)
        real utlnv(nlvddi,nel)
        real vtlnv(nlvddi,nel)
        real uprv(nlvddi,nkn)
        real vprv(nlvddi,nkn)

        real weight(nlvddi,nkn)         !aux variable for weights
        real hl(nlvddi)                 !aux variable for real level thickness

        logical bsigma
        integer ie,ii,k,l,lmax,nsigma,nlvaux
        real hmed,u,v,area,zeta
        real hsigma

        real area_elem

        call get_sigma_info(nlvaux,nsigma,hsigma)
        if( nlvaux .gt. nlvddi ) stop 'error stop transp2vel: nlvddi'
        bsigma = nsigma .gt. 0

	weight = 0.
	uprv = 0.
	vprv = 0.
	hl = 1.		!in case of transports

        do ie=1,nel

          area = area_elem(ie)
          lmax = ilhv(ie)
	  if( bvel ) then
	    zeta = sum(zenv(:,ie)) / 3.	!average of zeta on element
	    call get_layer_thickness(lmax,nsigma,hsigma
     +				,zeta,hev(ie),hlv,hl)
	  end if

          do l=1,lmax
            hmed = hl(l)
            u = utlnv(l,ie) / hmed
            v = vtlnv(l,ie) / hmed
            do ii=1,3
              k = nen3v(ii,ie)
              uprv(l,k) = uprv(l,k) + area * u
              vprv(l,k) = vprv(l,k) + area * v
              weight(l,k) = weight(l,k) + area
            end do
          end do
        end do

	where( weight > 0. )
	  uprv = uprv / weight
	  vprv = vprv / weight
	end where

        end

c******************************************************************

