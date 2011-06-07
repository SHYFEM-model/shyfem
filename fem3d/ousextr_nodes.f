c
c $Id: ousextr.f,v 1.3 2009-09-14 08:20:58 georg Exp $
c
c extract nodes from OUS file
c
c revision log :
c
c 02.09.2003	ggu	adapted to new OUS format
c 24.01.2005	ggu	computes maximum velocities for 3D (only first level)
c 03.06.2011    ggu     routine adjourned
c
c***************************************************************

	program ousextr_nodes

c reads ous file and extracts nodes
c
c we would not even need to read basin
c
c still to be revised...

	implicit none

        include 'param.h'

	character*80 descrr,descrp
	common /descrr/ descrr
	common /descrp/ descrp
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	real xgv(nkndim), ygv(nkndim)
	real hm3v(3,neldim)
	integer nen3v(3,neldim)
	integer ipev(neldim), ipv(nkndim)
	integer iarv(neldim)
	common /xgv/xgv, /ygv/ygv
	common /hm3v/hm3v
	common /nen3v/nen3v
	common /ipev/ipev, /ipv/ipv
	common /iarv/iarv

	integer ilhv(neldim)
	integer ilhkv(nkndim)
	real hlv(nlvdim)
        real utlnv(nlvdim,neldim)
        real vtlnv(nlvdim,neldim)
	common /ilhv/ilhv
	common /ilhkv/ilhkv
	common /hlv/hlv
        common /utlnv/utlnv
        common /vtlnv/vtlnv

	real hev(neldim)

	real znv(nkndim)
	real zenv(3,neldim)

	real uprv(nlvdim,nkndim)
	real vprv(nlvdim,nkndim)
	real weight(nlvdim,nkndim)
	real ut2v(neldim)
	real vt2v(neldim)
	real u2v(neldim)
	real v2v(neldim)

	integer ilnv(nlvdim,nkndim)

        integer nvers,nin,nlv
        integer itanf,itend,idt,idtous
	integer it,ie,i
        integer ierr,nread,ndry
        integer nknous,nelous,nlvous
        real href,hzoff,hlvmin
	real zmin,zmax
	real umin,umax
	real vmin,vmax
	real xe,ye
	integer k,ke,ivar
	integer lmax,l

c	integer rdous,rfous
	integer iapini,ideffi
	logical berror

c---------------------------------------------------------------
c nodes for extraction
c---------------------------------------------------------------

        integer ndim
        integer nnodes
        parameter( ndim = nkndim )
        integer nodes(ndim)     !node numbers
        integer nodese(ndim)    !external node numbers

c-------------------------------------------------------------------
c initialize params
c-------------------------------------------------------------------

        nread=0

c---------------------------------------------------------------
c open simulation and basin
c---------------------------------------------------------------

	if(iapini(3,nkndim,neldim,0).eq.0) then
		stop 'error stop : iapini'
	end if

c--------------------------------------------------------------------
c open OUS file and read header
c--------------------------------------------------------------------

	nin=ideffi('datdir','runnam','.ous','unform','old')
	if(nin.le.0) goto 100

	nvers=1
        call rfous(nin
     +			,nvers
     +			,nknous,nelous,nlvous
     +			,href,hzoff
     +			,descrp
     +			,ierr)
	if(ierr.ne.0) goto 100

        write(6,*)
        write(6,*)   descrp
        write(6,*)
        write(6,*) ' nvers        : ',nvers
        write(6,*) ' href,hzoff   : ',href,hzoff
        write(6,*) ' nkn,nel      : ',nknous,nelous
        write(6,*) ' nlv          : ',nlvous
        write(6,*)

        if( nkn .ne. nknous .or. nel .ne. nelous ) goto 94

	nlv=nlvous
	call dimous(nin,nkndim,neldim,nlvdim)

	call rsous(nin,ilhv,hlv,hev,ierr)
        if(ierr.ne.0) goto 100

        call level_e2k(nkn,nel,nen3v,ilhv,ilhkv)

        write(6,*) 'Available levels: ',nlv
        write(6,*) (hlv(l),l=1,nlv)

c-------------------------------------------------------------------
c get nodes to extract from STDIN
c-------------------------------------------------------------------

        call get_nodes_from_stdin(ndim,nnodes,nodes,nodese)

        if( nnodes .le. 0 ) goto 100

c-------------------------------------------------------------------
c loop on input records
c-------------------------------------------------------------------

  300   continue

        call rdous(nin,it,nlvdim,ilhv,znv,zenv,utlnv,vtlnv,ierr)

        if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
        if(ierr.ne.0) goto 100

	nread=nread+1
	write(6,*) 'time : ',it,nread

	call comp_barotropic(nel,nlvdim,ilhv,utlnv,vtlnv,ut2v,vt2v)
        call comp_vel2d(nel,hev,zenv,ut2v,vt2v,u2v,v2v,umax,vmax)

        call transp2vel(nel,nkn,nlv,nlvdim,hev,zenv,nen3v
     +                          ,ilhv,hlv,utlnv,vtlnv
     +                          ,uprv,vprv,weight)

        do i=1,nnodes
          k = nodes(i)
          ke = nodese(i)
	  lmax = ilhkv(k)
          write(79,*) it,i,ke,k,lmax
          write(79,*) znv(k)
          write(79,*) (uprv(l,k),l=1,lmax)
          write(79,*) (vprv(l,k),l=1,lmax)
        end do

	write(80,*) it,(znv(nodes(k)),k=1,nnodes)

	goto 300

  100	continue

c-------------------------------------------------------------------
c end of loop
c-------------------------------------------------------------------

	write(6,*)
	write(6,*) nread,' records read'
	write(6,*)
	write(6,*) 'data written to file 79 and 80'
	write(6,*)

c-------------------------------------------------------------------
c end of routine
c-------------------------------------------------------------------

	stop
   94   continue
        write(6,*) 'incompatible simulation and basin'
        write(6,*) 'nkn: ',nkn,nknous
        write(6,*) 'nel: ',nel,nelous
        stop 'error stop ousextr_nodes: nkn,nel'
	end

c******************************************************************

        subroutine comp_vel2d(nel,hev,zenv,ut2v,vt2v,u2v,v2v,umax,vmax)

c computes velocity in elements for given level 
c
c returns result in uv,vv

        implicit none

        integer level		!level for which to compute velocity
        integer nel
        real hev(1)
        real zenv(3,1)
        real ut2v(1)
        real vt2v(1)
	real u2v(1), v2v(1)
        real umax,vmax

        integer ie,ii
        real zmed,hmed,u,v

        umax = 0.
        vmax = 0.

        do ie=1,nel
          zmed = 0.
          do ii=1,3
            zmed = zmed + zenv(ii,ie)
          end do
          zmed = zmed / 3.
          hmed = hev(ie) + zmed

          u = ut2v(ie) / hmed
          v = vt2v(ie) / hmed

	  u2v(ie) = u
	  v2v(ie) = v

          umax = max(umax,u)
          vmax = max(vmax,v)
        end do

        end

c******************************************************************

	subroutine comp_barotropic(nel,nlvdim,ilhv
     +			,utlnv,vtlnv,ut2v,vt2v)

c computes barotropic transport

	implicit none

	integer nel,nlvdim
	integer ilhv(1)
	real utlnv(nlvdim,1)
	real vtlnv(nlvdim,1)
	real ut2v(1)
	real vt2v(1)

	integer ie,l,lmax
	real utot,vtot

	do ie=1,nel
	  lmax = ilhv(ie)
	  utot = 0.
	  vtot = 0.
	  do l=1,lmax
	    utot = utot + utlnv(l,ie)
	    vtot = vtot + vtlnv(l,ie)
	  end do
	  ut2v(ie) = utot
	  vt2v(ie) = vtot
	end do

	end

c******************************************************************

