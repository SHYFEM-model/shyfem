c
c $Id: nosextr_nodes.f,v 1.3 2009-04-07 10:43:57 georg Exp $
c
c extract nodes from NOS file
c
c revision log :
c
c 18.11.1998    ggu     check dimensions with dimnos
c 24.02.1999    ggu     use n2int for node number translation
c 03.12.2001    ggu     cleaned up, hakata bay
c 03.06.2011    ggu     routine adjourned
c 08.06.2011    ggu     routine rewritten for sediments
c
c****************************************************************

	program nosextr_fluxes

c computes fluxes of scalar (sediment) through section

	use basin !COMMON_GGU_SUBST

	implicit none

	include 'param.h'

c--------------------------------------------------
COMMON_GGU_DELETED	include 'basin.h'
	include 'simul.h'


c--------------------------------------------------

	character*80 title
	real cv(nkndim)
	real cv3(nlvdim,nkndim)

	integer ilhkv(nkndim)
	real hlv(nlvdim)
	real hev(neldim)

	integer ilhv(neldim)
        real utlnv(nlvdim,neldim)
        real vtlnv(nlvdim,neldim)
        real znv(nkndim)
        real zenv(3,neldim)
        real uprv(nlvdim,nkndim)
        real vprv(nlvdim,nkndim)
        real utprv(nlvdim,nkndim)
        real vtprv(nlvdim,nkndim)
        real weight(nlvdim,nkndim)

	logical berror
	integer i,n,k,ke,l,lmax
	integer nread
	integer nvers
	integer nlv,nvar,ivar,ierr
	integer nin,it
	integer ninnos,ninous
	integer nknous,nelous,nlvous
	real href,hzoff
	real cflux,cfold
	integer itnos,itous,itold
	double precision ctot
	double precision cft,cfp,cfm

	integer iapini,ideffi,ialfa

c---------------------------------------------------------------
c nodes for extraction
c---------------------------------------------------------------

	integer ndim
	integer nnodes
	parameter( ndim = nkndim )
	integer nodes(ndim)	!node numbers
	integer nodese(ndim)	!external node numbers
	real flx(2,ndim)

c---------------------------------------------------------------
c open simulation and basin
c---------------------------------------------------------------

	if(iapini(3,nkndim,neldim,0).eq.0) then
		stop 'error stop : iapini'
	end if

c---------------------------------------------------------------
c open NOS file and read header (only sediments)
c---------------------------------------------------------------

	ninnos=ideffi('datdir','runnam','.sco','unform','old')
	if(ninnos.le.0) goto 100

        nvers=3
	call rfnos(ninnos,nvers,nkn,nel,nlv,nvar,title,ierr)
        if(ierr.ne.0) goto 100

        write(6,*) 'nvers    : ',nvers
        write(6,*) 'nkn,nel  : ',nkn,nel
        write(6,*) 'nlv,nvar : ',nlv,nvar
        write(6,*) 'title    : ',title

        call dimnos(ninnos,nkndim,neldim,nlvdim)

	call rsnos(ninnos,ilhkv,hlv,hev,ierr)
        if(ierr.ne.0) goto 100

	write(6,*) 'Available levels: ',nlv
	write(6,*) (hlv(l),l=1,nlv)

c--------------------------------------------------------------------
c open OUS file and read header
c--------------------------------------------------------------------

        ninous=ideffi('datdir','runnam','.ous','unform','old')
        if(ninous.le.0) goto 100

        nvers=1
        call rfous(ninous
     +                  ,nvers
     +                  ,nknous,nelous,nlvous
     +                  ,href,hzoff
     +                  ,descrp
     +                  ,ierr)
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
        call dimous(ninous,nkndim,neldim,nlvdim)

        call rsous(ninous,ilhv,hlv,hev,ierr)
        if(ierr.ne.0) goto 100

        !call level_e2k(nkn,nel,nen3v,ilhv,ilhkv)

        write(6,*) 'Available levels: ',nlv
        write(6,*) (hlv(l),l=1,nlv)

c---------------------------------------------------------------
c initializing units and nodes to be extracted
c---------------------------------------------------------------

	nread=0

        call get_nodes_from_stdin(ndim,nnodes,nodes,nodese)
	call setup_fluxes(nnodes,nodes,flx)

	if( nnodes .le. 0 ) stop

c---------------------------------------------------------------
c loop on input records
c---------------------------------------------------------------

	itold = -1
	cfold = 0.
	ctot = 0.
	cft = 0.
	cfp = 0.
	cfm = 0.
	itnos = -87654321
	itous = -87654321

  300   continue

	do while( itnos .le. itous )
	  call rdnos(ninnos,itnos,ivar,nlvdim,ilhkv,cv3,ierr)
          if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
          if(ierr.ne.0) goto 100
	end do

	do while( itous .lt. itnos )
          call rdous(ninous,itous,nlvdim,ilhv,znv,zenv,utlnv,vtlnv,ierr)
          if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
          if(ierr.ne.0) goto 100
	end do

	if( itnos .ne. itous ) goto 300

	it = itnos
	nread=nread+1
	write(6,*) 'time : ',it,ivar

c	---------------------------------------------------------
c	write to file
c	---------------------------------------------------------

        call transp2nodes(nel,nkn,nlv,nlvdim,hev,zenv,nen3v
     +                          ,ilhv,hlv,utlnv,vtlnv
     +                          ,utprv,vtprv,weight)

	call compute_fluxes(nnodes,nodes,flx
     +			,nlvdim,ilhkv,utprv,vtprv,cv3,cflux)

	write(66,*) it,cflux

	if( itold .ne. -1 ) then
	  ctot = ctot + (it-itold)*0.5*(cfold+cflux)
	  call accum_fluxes(it-itold,cfold,cflux,cft,cfp,cfm)
	end if
	itold = it
	cfold = cflux

	!if( nread .gt. 100 ) goto 100
	goto 300

  100	continue

c---------------------------------------------------------------
c end of loop
c---------------------------------------------------------------

	write(6,*)
	write(6,*) nread,' records read'
	write(6,*)
	write(6,*) 'accumulated sediment flux: ',ctot
	write(6,*) cft,cfp,cfm
	write(99,*) cft,cfp,cfm
	write(6,*)
	write(6,*) 'data written to file 66 and 99'
	write(6,*)

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

	stop
   94   continue
        write(6,*) 'incompatible simulation and basin'
        write(6,*) 'nkn: ',nkn,nknous
        write(6,*) 'nel: ',nel,nelous
        stop 'error stop nosextr_fluxes: nkn,nel'
	end

c***************************************************************

	subroutine accum_fluxes(idt,cfoin,cfnin,cftout,cfpout,cfmout)

	implicit none

	integer idt
	real cfoin,cfnin
	double precision cftout,cfpout,cfmout

	double precision cfo,cfn
	double precision cft,cfp,cfm
	double precision frac,dt

	cfo = cfoin
	cfn = cfnin
	dt = 0.5*idt

	cft = dt*(cfo+cfn)
	cfp = 0.
	cfm = 0.
	frac = 0.

	if( cfo*cfn .gt. 0. ) then
	  if( cfo .gt. 0. ) then
	    cfp = dt*(cfo+cfn)
	  else
	    cfm = dt*(cfo+cfn)
	  end if
	else
	  if( cfo .gt. 0. ) then
	    frac = cfo/(cfo-cfn)
	    cfp = dt*frac*cfo
	    cfm = dt*(1.-frac)*cfn
	  else if( cfo .lt. 0. ) then
	    frac = cfo/(cfo-cfn)
	    cfm = dt*frac*cfo
	    cfp = dt*(1.-frac)*cfn
	  else	! if( cfo .eq. 0. ) then
	    if( cfn .gt. 0. ) then
	      cfp = dt*cfn
	    else
	      cfm = dt*cfn
	    end if
	  end if
	  if( frac .lt. 0. ) then
		write(6,*) 'frac: ',frac
		stop 'error stop accum_fluxes: frac'
	  end if
	end if

	if( abs(cft-(cfm+cfp)) .gt. 1.e-5 ) then
	  write(6,*) idt,cfo,cfn,frac
	  write(6,*) cft,cfp,cfm,cfp+cfm
	  stop 'error stop accum_fluxes: diff'
	end if

	cftout = cftout + cft
	cfpout = cfpout + cfp
	cfmout = cfmout + cfm

	end

c***************************************************************

	subroutine compute_fluxes(nnodes,nodes,flx
     +			,nlvdim,ilhkv,utprv,vtprv,cv3,cflux)

	implicit none

	integer nnodes
	integer nodes(nnodes)
	real flx(2,nnodes)
	integer nlvdim
	integer ilhkv(1)
	real utprv(nlvdim,1)
	real vtprv(nlvdim,1)
	real cv3(nlvdim,1)
	real cflux		!computed mass flux through section

	integer i,k,l,lmax
	real xn,yn,ut,vt,uvs
	double precision acum

	acum = 0.

	do i=1,nnodes
	  xn = flx(1,i)
	  yn = flx(2,i)
	  k = nodes(i)
	  lmax = ilhkv(k)
	  do l=1,lmax
	    ut = utprv(l,k)
	    vt = vtprv(l,k)
	    uvs = xn*ut + yn*vt
	    acum = acum + uvs * cv3(l,k)
	  end do
	end do

	cflux = acum

	end

c***************************************************************

	subroutine setup_fluxes(nnodes,nodes,flx)

	use basin !COMMON_GGU_SUBST

	implicit none

	include 'param.h'

	integer nnodes
	integer nodes(nnodes)
	real flx(2,nnodes)


COMMON_GGU_DELETED	include 'basin.h'

	integer i,ii,i1,i2,im,ip
	integer k,k1,k2
	integer ie,ief
	real x1,x2,y1,y2,xt,yt,xn,yn

	do i=2,nnodes
	  k1 = nodes(i-1)
	  k2 = nodes(i)
	  ief = 0
	  do ie=1,nel
	    do ii=1,3
	      i1 = mod(ii,3)+1
	      i2 = mod(i1,3)+1
	      if( nen3v(i1,ie) .eq. k1 .and. nen3v(i2,ie) .eq. k2 ) then
		ief = ie
	      end if
	    end do
	  end do
	  write(6,*) k1,k2,ief
	  if( ief .eq. 0 ) then
	    write(6,*) i,k1,k2
	    stop 'error stop: nodes'
	  end if
	end do

	do i=1,nnodes
	  im = max(1,i-1)
	  ip = min(nnodes,i+1)
	  k1 = nodes(im)
	  k2 = nodes(ip)
	  x1 = xgv(k1)
	  y1 = ygv(k1)
	  x2 = xgv(k2)
	  y2 = ygv(k2)

	  xt = x2 - x1
	  yt = y2 - y1
	  xn = -yt
	  yn =  xt

	  flx(1,i) = xn/2.
	  flx(2,i) = yn/2.
	end do

	end

c***************************************************************





