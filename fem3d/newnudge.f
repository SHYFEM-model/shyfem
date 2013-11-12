c
c $Id: baswork.f,v 1.20 2010-03-22 15:29:31 georg Exp $
c
c general framework to work on files needing basin
c
c revision log :
c
c 06.04.1999	ggu	completely restructured
c 04.06.1999	ggu	new statistics are computed
c 28.11.2005	ggu	new call to makehkv
c 31.05.2007	ggu	added area and volume frequency curve
c 24.08.2007	ggu	added new routine write_grd_from_bas
c 06.04.2009    ggu     read param.h
c 12.06.2009    ggu     areatr in double precision - new algorithm
c 01.03.2010    ggu     new routine basqual() to compute grid quality
c 22.03.2010    ggu     write external element number in basqual()
c 17.05.2011    ggu     changes in freqdep()
c 12.07.2011    ggu     better treatment of freqdep()
c 16.11.2011    ggu     basin.h introduced
c 23.01.2012    ggu     new from basinf
c
c****************************************************************

	subroutine nudge_init

	implicit none

	include 'param.h'
	include 'basin.h'
	include 'nudge.h'

	real andgzv(nkndim)		!contribution to computation
	common /andgzv/andgzv

	logical binfl,binsert
	integer nintp,nsize,ndim
	integer k,i,node,ivar
	real ttau,sigma
	character*40 file_obs,file_stations
	character*16 cname
	character*4 cdummy

	integer ipint

	nvar = 27
	nvar = 30
	nvar = 0
	file_obs='data_level_18-27_06_2008.dat'
	file_stations='mareo_input.txt'
	tramp = 43200
	ttau = 0
	ttau = 3600
	ttau = 600
	ttau = 100
	sigma = 0
	sigma = 2000
	ivar = 0

	binfl = sigma .gt. 0.

	do k=1,nkn
	  andgzv(k) = 0.
	end do

	if( nvar .gt. ndgdim ) stop 'error stop nudge_init: ndgdim'

	if( nvar .le. 0 ) return

	nintp = 4
	nsize = 1
	ndim = ndgdatdim

	call exffil(file_obs,nintp,nvar,nsize,ndim,andg_data)

	ivar = 0
	open(1,file='NUDGE',status='old',form='formatted',err=1)
	read(1,*) ivar
	close(1)
    1	continue

	do i=1,nvar
	  ndg_use(i) = 1			! use all data
	  !if( mod(i,2) .ne. 0 ) ndg_use(i) = 0	! only use even vars
	  !if( i .eq. 10 ) ndg_use(i) = 0	! no Saline
	  if( i .eq. ivar ) ndg_use(i) = 0
	end do

	open(1,file=file_stations,status='old',form='formatted')
	write(6,*) 'initializing zeta nudging... nvar = ',nvar
	do i=1,nvar
	  read(1,'(a16,a4,i6)') cname,cdummy,node
	  write(6,*) i,cname,node,ndg_use(i)
	  k = ipint(node)
	  if( k .le. 0 ) goto 99
	  ndg_names(i) = cname
	  ndg_nodelist(i) = k
	end do
	write(6,*) 'finished initializing zeta nudging...'
	close(1)

	do k=1,nkn
	  andg_dist(k) = 0.
	  andg_weight(k) = 0.
	  andg_obs(k) = 0.
	  ndg_nodes(k) = 0
	  ndg_area(k) = 0
	end do

	do i=1,nvar
	  binsert = ndg_use(i) .gt. 0
	  if( binsert ) then
	    k = ndg_nodelist(i)
	    ndg_nodes(k) = k
	    ndg_area(k) = i
	    andg_weight(k) = 0.
	    if( ttau .gt. 0. ) andg_weight(k) = 1./ttau
	  end if
	end do

	if( binfl ) then
	  call nudge_influence
	  call distance(ttau,sigma)
	end if

	return
   99	continue
	write(6,*) 'Cannot find internal node: ',node
	stop 'error stop nudge_init: no internal node'
	end

c****************************************************************

	subroutine nudge_zeta

	implicit none

	include 'param.h'
	include 'basin.h'
	include 'nudge.h'

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

	real andgzv(nkndim)		!contribution to z-computation
	common /andgzv/andgzv
	real zov(nkndim)
	common /zov/zov

	integer i,k,kk,ia
	integer nnudge,iuse,iu
	real talpha,ttau,t,zobs,zcontrib,w
	real rint(ndgdim)
	real zeta(ndgdim)
	real cont(ndgdim)

	if( nvar .le. 0 ) return

	t = it

	call exfintp(andg_data,t,rint)

	talpha = 1.
	if( tramp .gt. 0. ) talpha = (it-itanf)/tramp
	if( talpha .gt. 1. ) talpha = 1.

	do i=1,nvar
	  k = ndg_nodelist(i)
	  zeta(i) = zov(k)
	end do

	write(88,'(i10,50f7.3)') it,(rint(i),i=1,nvar)
	write(89,'(i10,50f7.3)') it,(zeta(i),i=1,nvar)

	do i=1,nvar
	  k = ndg_nodelist(i)
	  iuse = ndg_use(i)
	  if( iuse .le. 0 ) then
	    iu = 100 + i
	    write(iu,'(i10,f7.3)') it,zov(k)
	  end if
	end do

	nnudge = 0
	do k=1,nkn
	  w = andg_weight(k)
	  !kk = ndg_nodes(k)
	  ia = ndg_area(k)
	  zcontrib = 0.
	  if( ia .gt. 0 ) then
	    zobs = rint(ia)
	    zcontrib = talpha * w * ( zobs - zov(k) )
	    nnudge = nnudge + 1
	  end if
	  !zcontrib = 0.
	  andgzv(k) = zcontrib
	end do

	write(6,*) 'nudging zeta: ',nnudge

	end

c*******************************************************************

c        do ii=1,3
c          fnudge(ii) = robs * rtau(l,ii) * ( cob(l,ii) - cl(l,ii) )
c        end do
c          cdummy = aj4 * ( hold(l,ii)*cl(l,ii)
c     +                          + dt *  (
c     +                                      hold(l,ii)*fnudge(ii)
c          cn(l,k) = cn(l,k) + cdummy
c
c          co(l,k) = co(l,k) + dt * aj4 * hmed * robs * rtauv(l,k) !nudging

c*******************************************************************

	subroutine distance(ttau,sigma)

c sets up andg_dist and andg_weight

	implicit none

	include 'param.h'
	include 'basin.h'
	include 'nudge.h'

	real ttau		!relaxation time
	real sigma		!sigma distance

	integer k,kk,iu
	real dx,dy,dist,dist2,w,s

	real d(nkndim)
	real dw(nkndim)

	s = 2.*sigma*sigma

	do k=1,nkn
	  kk = ndg_nodes(k)
	  dx = xgv(kk) - xgv(k)
	  dy = ygv(kk) - ygv(k)
	  dist2 = dx*dx+dy*dy
	  dist = sqrt(dist2)
	  andg_dist(k) = dist2
	  d(k) = dist			!only for output

	  w = exp(-dist2/s)
	  dw(k) = w

	  andg_weight(k) = w / ttau
	end do

	iu = 0
	call conwrite(iu,'.dst',1,988,1,d)
	call nos_close(iu)
	close(iu)

	iu = 0
	call conwrite(iu,'.www',1,989,1,dw)
	call nos_close(iu)
	close(iu)

	end

c*******************************************************************

	subroutine nudge_influence

c sets up icol (local) and ndg_nodes and ndg_area

	implicit none

	include 'param.h'
	include 'basin.h'
	include 'nudge.h'

	integer icol(nkndim)
	integer icolaux(nkndim)
	real d(nkndim)

	logical binsert
	integer k,ie,ii,it,is,ic,isc,i
	integer iu
	integer ngood

	do k=1,nkn
	  icol(k) = 0
	  icolaux(k) = 0
	end do

	do i=1,nvar
	  binsert = ndg_use(i) .gt. 0
	  if( binsert ) then
	    k = ndg_nodelist(i)
	    icol(k) = i
	  end if
	end do

	ngood = 0
	do k=1,nkn
	  if( icol(k) .gt. 0 ) ngood = ngood + 1
	end do

	do while( ngood .lt. nkn )

	  write(6,*) 'influence: ',ngood,nkn

	  do ie=1,nel

	    it = 0
	    is = 0
	    do ii=1,3
	      k = nen3v(ii,ie)
	      if( icol(k) .gt. 0 ) then
	        it = it + 1
	        is = is + ii
	      end if
	    end do
	
	    if( it .eq. 1 ) then
	      k = nen3v(is,ie)
	      ic = icol(k)
	      do ii=1,3
	        k = nen3v(ii,ie)
	        icolaux(k) = ic
	      end do
	    else if( it .eq. 2 ) then
	      is = 6 - is
	      isc = mod(is,3) + 1
	      k = nen3v(isc,ie)
	      ic = icol(k)
	      k = nen3v(is,ie)
	      icolaux(k) = ic
	    end if

	    if( it .gt. 0 ) then	!check
	      do ii=1,3
	        k = nen3v(ii,ie)
	        if( icol(k) .le. 0 .and. icolaux(k) .le. 0 ) then
		  write(6,*) 'internal error...... '
		  write(6,*) ie,it,is
		  stop 'error stop: internal error'
	        end if
	      end do
	    end if

	  end do

	  ngood = 0
	  do k=1,nkn
	    if( icolaux(k) .gt. 0 ) then
		if( icol(k) .le. 0 ) icol(k) = icolaux(k)
	    end if
	    if( icol(k) .gt. 0 ) ngood = ngood + 1
	  end do

	end do

	do k=1,nkn
	  ic = icol(k)
	  ndg_area(k) = ic
	  ndg_nodes(k) = ndg_nodelist(ic)
	end do

	do k=1,nkn
	  d(k) = icol(k)
	end do

	iu = 0
	call conwrite(iu,'.col',1,987,1,d)
	call nos_close(iu)
	close(iu)

	end

c*******************************************************************

