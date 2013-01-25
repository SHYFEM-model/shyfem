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
c 16.12.2011    ggu     bug: call to init_sigma_info and makehev (common hev)
c 09.03.2012    ggu     bug: zenv must be common
c 25.01.2013    ggu     code cleaned
c
c***************************************************************

	program ousextr_nodes

c extracts single nodes from OUS file -> creates time series

	implicit none

        include 'param.h'
        include 'basin.h'
        include 'evmain.h'

c--------------------------------------------------

	character*80 title

	integer ilhv(neldim)
	integer ilhkv(nkndim)
	real hlv(nlvdim)
	real hev(neldim)
	common /hev/hev, /hlv/hlv       ! need this for get_layer_thickness

	real zenv(3,neldim)
	common /zenv/zenv		! need this for get_layer_thickness

	real znv(nkndim)
        real utlnv(nlvdim,neldim)
        real vtlnv(nlvdim,neldim)
	real uprv(nlvdim,nkndim)
	real vprv(nlvdim,nkndim)
	real weight(nlvdim,nkndim)
	real hl(nlvdim)
	real ut2v(neldim)
	real vt2v(neldim)
	real u2v(neldim)
	real v2v(neldim)

        integer nread,ierr
        integer nvers,nin,nlv
        integer nknous,nelous

	integer it
	integer k,ke,i
	integer l,lmax

        real href,hzoff
	real umin,vmin,umax,vmax

	integer iapini,ideffi

c---------------------------------------------------------------
c nodes for extraction
c---------------------------------------------------------------

        integer ndim
        integer nnodes
        parameter( ndim = nkndim )
        integer nodes(ndim)     !node numbers
        integer nodese(ndim)    !external node numbers

c---------------------------------------------------------------
c initialize params
c---------------------------------------------------------------

        nread=0

c---------------------------------------------------------------
c open simulation and basin
c---------------------------------------------------------------

	if(iapini(3,nkndim,neldim,0).eq.0) then
		stop 'error stop : iapini'
	end if

	call set_ev

c---------------------------------------------------------------
c open OUS file and read header
c---------------------------------------------------------------

	nin=ideffi('datdir','runnam','.ous','unform','old')
	if(nin.le.0) goto 100

	nvers=1
        call rfous(nin
     +			,nvers
     +			,nknous,nelous,nlv
     +			,href,hzoff
     +			,title
     +			,ierr)
	if(ierr.ne.0) goto 100

        write(6,*)
        write(6,*) 'title        : ',title
        write(6,*)
        write(6,*) 'nvers        : ',nvers
        write(6,*) 'nkn,nel      : ',nknous,nelous
        write(6,*) 'nlv          : ',nlv
        write(6,*) 'href,hzoff   : ',href,hzoff
        write(6,*)

        if( nkn .ne. nknous .or. nel .ne. nelous ) goto 94

	call dimous(nin,nkndim,neldim,nlvdim)

	call rsous(nin,ilhv,hlv,hev,ierr)
        if(ierr.ne.0) goto 100

	call init_sigma_info(nlv,hlv)
        call level_e2k(nkn,nel,nen3v,ilhv,ilhkv)
	call makehev(hev)

        write(6,*) 'Available levels: ',nlv
        write(6,*) (hlv(l),l=1,nlv)

c---------------------------------------------------------------
c get nodes to extract from STDIN
c---------------------------------------------------------------

        call get_nodes_from_stdin(ndim,nnodes,nodes,nodese)

        if( nnodes .le. 0 ) goto 100

c---------------------------------------------------------------
c loop on input records
c---------------------------------------------------------------

  300   continue

        call rdous(nin,it,nlvdim,ilhv,znv,zenv,utlnv,vtlnv,ierr)

        if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
        if(ierr.ne.0) goto 100

	nread=nread+1
	write(6,*) 'time : ',it,nread

	call comp_barotropic(nel,nlvdim,ilhv,utlnv,vtlnv,ut2v,vt2v)
        call comp_vel2d(nel,hev,zenv,ut2v,vt2v,u2v,v2v
     +					,umin,vmin,umax,vmax)

        call transp2vel(nel,nkn,nlv,nlvdim,hev,zenv,nen3v
     +                          ,ilhv,hlv,utlnv,vtlnv
     +                          ,uprv,vprv,weight,hl)

c	---------------------------------------------------------
c	write to file
c	---------------------------------------------------------

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

c---------------------------------------------------------------
c end of loop
c---------------------------------------------------------------

	write(6,*)
	write(6,*) nread,' records read'
	write(6,*)
	write(6,*) 'data written to file 79 and 80'
	write(6,*)

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

	stop
   94   continue
        write(6,*) 'incompatible simulation and basin'
        write(6,*) 'nkn: ',nkn,nknous
        write(6,*) 'nel: ',nel,nelous
        stop 'error stop ousextr_nodes: nkn,nel'
	end

c***************************************************************

