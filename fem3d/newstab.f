c
c $Id: newstab.f,v 1.3 2010-03-08 17:46:45 georg Exp $
c
c routines for stability computations
c
c revision log :
c
c 19.02.2010    ggu     new file to contain stability computations
c 26.02.2010    ggu     internal_stability restructured (on element)
c 08.03.2010    ggu     run only down to avail layers in stability (bug fix)
c 26.01.2011    ggu     robs for nudging implemented
c 16.02.2011    ggu     pass robs to subroutines, write stb-ind to nos file
c
c*****************************************************************
c*****************************************************************
c*****************************************************************
c
c typical usage:
c
c call reset_stability at beginning of time loop
c call make_stability before advection of similar variables
c
c*****************************************************************

	subroutine compute_stability(robs,rkpar,azpar,rindex,saux)

c computes stability index

	implicit none

        include 'param.h'

	real robs
        real rkpar
        real azpar
        real rindex
        real saux(nlvdim,nkndim)

        integer nlvdi,nlv
        common /level/ nlvdi,nlv

        real cnv(nlvdim,1)
        common /cnv/cnv
        real difv(0:nlvdim,1)
        common /difv/difv
        real difhv(nlvdim,1)
	common /difhv/difhv

        real adpar,aapar
        real difmol
	real ddt
        integer isact,istot

        real getpar

c----------------------------------------------------------------
c set parameters
c----------------------------------------------------------------

	adpar=getpar('adpar')
	aapar=getpar('aapar')

        isact = 1
        istot = 1
        difmol = 0.
	ddt = 1.		!always for 1 sec
	rindex = 0.

c----------------------------------------------------------------
c call conzstab
c----------------------------------------------------------------

        call conzstab(cnv,saux
     +          ,ddt,robs,rkpar,difhv,difv
     +		,difmol,azpar,adpar,aapar
     +          ,rindex,istot,isact,nlvdi,nlv)

c----------------------------------------------------------------
c end of routine
c----------------------------------------------------------------

	end

c*****************************************************************

        subroutine make_stability(dt,robs,rkpar,rindex,istot,saux)

c gets stability index (if necessary computes it)

        implicit none

	include 'param.h'

	real dt
	real robs
        real rkpar
        real rindex
        integer istot
	real saux(nlvdim,nkndim)

	real azpar
	logical exist_stability

c----------------------------------------------------------------
c see if already computed
c----------------------------------------------------------------

	if( exist_stability(rkpar,rindex) ) goto 1
	  
c----------------------------------------------------------------
c compute stability index
c----------------------------------------------------------------

	call getaz(azpar)
	call compute_stability(robs,rkpar,azpar,rindex,saux)

c----------------------------------------------------------------
c insert stability index
c----------------------------------------------------------------

	call insert_stability(rkpar,rindex)

c----------------------------------------------------------------
c scale to real time step dt
c----------------------------------------------------------------

    1	continue
	rindex = dt * rindex
	istot = 1 + rindex

c----------------------------------------------------------------
c end of routine
c----------------------------------------------------------------

        end

c*****************************************************************

        subroutine info_stability(robs,dt,rkpar,rindex,istot,saux)

c gets stability index (if necessary computes it)

        implicit none

	include 'param.h'

	real robs
        real dt
        real rkpar
        real rindex
        integer istot
	real saux(nlvdim,nkndim)

        integer nlv,nlvdi
        common /level/ nlvdi,nlv
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	integer ia,iustab
	integer l,k
	real aindex
	real azpar
	logical exist_stability

c----------------------------------------------------------------
c compute stability index
c----------------------------------------------------------------

	call getaz(azpar)
	call compute_stability(robs,rkpar,azpar,rindex,saux)
	rindex = dt * rindex
	istot = 1 + rindex

c----------------------------------------------------------------
c write to terminal
c----------------------------------------------------------------

	ia = 1
	aindex = saux(1,1)

	do k=1,nkn
	  do l=1,nlv
	    if( saux(l,k) .gt. aindex ) then
	      aindex = saux(l,k)
	      ia = k
	    end if
	  end do
	end do

cggu protect
	write(6,*) 'info_stability'
	write(6,*) rkpar,azpar,rindex,istot
	write(6,*) ia,aindex,dt*aindex
	iustab = 0
	call conwrite(iustab,'.sta',1,778,nlvdi,saux)
cggu protect

c----------------------------------------------------------------
c end of routine
c----------------------------------------------------------------

	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

        subroutine reset_stability

c rests stability index

        implicit none

	integer ndim
	parameter (ndim = 50)
	integer nentry
	real rkind(2,ndim)
	common /istab/nentry, /rstab/rkind
	save /istab/, /rstab/

	nentry = 0

        end

c**********************************************************************

	subroutine insert_stability(rkpar,rindex)

c inserts stability index

	implicit none

	real rkpar,rindex

	integer ndim
	parameter (ndim = 50)
	integer nentry
	real rkind(2,ndim)
	common /istab/nentry, /rstab/rkind
	save /istab/, /rstab/

	integer i

	do i=1,nentry
	  if( rkind(1,i) .eq. rkpar ) return
	end do
	if( i .gt. ndim ) goto 99

	!return		!uncomment for debug -> never insert

cggu protect
	nentry = i
	rkind(1,i) = rkpar
	rkind(2,i) = rindex
cggu protect

	return
   99	continue
	write(6,*) 'nentry,ndim: ',i,ndim
	stop 'error stop insert_stability: ndim'
	end

c**********************************************************************

	function exist_stability(rkpar,rindex)

c tests if stability index has already been computed and returns it

	implicit none

	logical exist_stability
	real rkpar,rindex

	integer ndim
	parameter (ndim = 50)
	integer nentry
	real rkind(2,ndim)
	common /istab/nentry, /rstab/rkind
	save /istab/, /rstab/

	integer i

	do i=1,nentry
	  if( rkind(1,i) .eq. rkpar ) goto 1
	end do
    1	continue

	exist_stability = i .le. nentry

	if( exist_stability ) then
	  rindex = rkind(2,i)
	else
	  rindex = 0
	end if

	end

c**********************************************************************
c**********************************************************************
c**********************************************************************

        subroutine hydro_stability(dt,rindex)

c computes stability index for hydro timestep - no error output

        implicit none

        real dt
        real rindex

	call internal_stability(0,dt,rindex)

	end

c**********************************************************************

        subroutine error_stability(dt,rindex)

c computes stability index for hydro timestep - with error output

        implicit none

        real dt
        real rindex

	call internal_stability(1,dt,rindex)

	end

c**********************************************************************

        subroutine internal_stability(mode,dt,rindex)

c computes stability index for hydro timestep (internal)

        implicit none

	include 'param.h'

	integer mode		!0: normal call  1:error output
        real dt
        real rindex

        integer nlv,nlvdi
        common /level/ nlvdi,nlv
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        real sauxe1(nlvdim,neldim)
        common /sauxe1/sauxe1
        real sauxe2(nlvdim,neldim)
        common /sauxe2/sauxe2
	integer ilhv(1)
	common /ilhv/ilhv

	integer ie,l,lmax
        real rkpar,azpar,ahpar
	real dindex,aindex,tindex

	real getpar

        rkpar = 0.
	azpar = 1.
	ahpar = getpar('ahpar')

	do ie=1,nel
	  do l=1,nlv
	    sauxe1(l,ie) = 0.
	    sauxe2(l,ie) = 0.
	  end do
	end do

	!call compute_stability(robs,rkpar,azpar,rindex,saux1)
	call momentum_advective_stability(aindex,sauxe1)
	call momentum_viscous_stability(ahpar,dindex,sauxe2)

	aindex = aindex*dt
	dindex = dindex*dt

	tindex = 0.
	do ie=1,nel
	  lmax = ilhv(ie)
	  do l=1,lmax
	    tindex = max(tindex,sauxe1(l,ie)+sauxe2(l,ie))
	  end do
	end do
	tindex = tindex*dt

	if( mode .eq. 1 ) then		!error output
	  write(6,*) 'internal_stability: '
	  write(6,*) aindex,dindex,tindex
	  call output_errout_stability(dt,sauxe1,sauxe2)
	end if

	rindex = tindex

	!write(6,*) 'rindex = ',rindex,aindex,dindex

        end

c*****************************************************************

        subroutine output_errout_stability(dt,sauxe1,sauxe2)

c outputs stability index for hydro timestep (internal)

        implicit none

	include 'param.h'

        real dt
	real sauxe1(nlvdim,neldim)
	real sauxe2(nlvdim,neldim)
	real sauxn(nlvdim,nkndim)

        integer nlv,nlvdi
        common /level/ nlvdi,nlv
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer ilhv(1)
	common /ilhv/ilhv

	logical bnos
	integer ie,l,lmax
	integer ia,id,it
	real aindex,dindex,tindex

c set ifnos in order to have output to nos file

	integer icall,iustab,ifnos
	save icall,iustab,ifnos
	data icall,iustab,ifnos /0,0,0/

	icall = icall + 1

	ia = 0
	id = 0
	it = 0
	aindex = 0.
	dindex = 0.
	tindex = 0.

	do ie=1,nel
	  lmax = ilhv(ie)
	  do l=1,lmax
	    if( sauxe1(l,ie) .gt. aindex ) then
	      aindex = sauxe1(l,ie)
	      ia = ie
	    end if
	    if( sauxe2(l,ie) .gt. dindex ) then
	      dindex = sauxe2(l,ie)
	      id = ie
	    end if
	    if( sauxe1(l,ie)+sauxe2(l,ie) .gt. tindex ) then
	      tindex = sauxe1(l,ie)+sauxe2(l,ie)
	      it = ie
	    end if
	  end do
	end do

	write(6,*) 'errout_stability: int-node stab-index stab-index*dt'
	write(6,*) 'advective: ',ia,aindex,aindex*dt
	write(6,*) 'diffusive: ',id,dindex,dindex*dt
	write(6,*) 'total:     ',it,tindex,tindex*dt

	if( ia .ne. 0 ) then
	  call check_elem(ia)
	end if
	if( id .ne. 0 .and. id .ne. ia ) then
	  call check_elem(id)
	end if
	if( it .ne. 0 .and. it .ne. id .and. it .ne. ia ) then
	  call check_elem(it)
	end if

	if( ifnos .gt. 0 .and. mod(icall,ifnos) .eq. 0 ) then
	  call e2n3d_minmax(+1,nlvdim,sauxe1,sauxn)
	  call conwrite(iustab,'.sta',1,778,nlvdi,sauxn)
	  call e2n3d_minmax(+1,nlvdim,sauxe2,sauxn)
	  call conwrite(iustab,'.sta',1,778,nlvdi,sauxn)
	  do ie=1,nel
	    do l=1,nlv
	      sauxe1(l,ie) = sauxe1(l,ie) + sauxe2(l,ie)
	    end do
	  end do
	  call e2n3d_minmax(+1,nlvdim,sauxe1,sauxn)
	  call conwrite(iustab,'.sta',1,778,nlvdi,sauxn)
	end if

	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine parallel_test

c tests parallel implementation

	implicit none

	include 'param.h'

	real dt,rkpar,azpar,rindex
	real saux(nlvdim,nkndim)

	azpar = 0.
	rkpar = 0.

	write(6,*) 'parallel test...'
	call compute_stability(0.,rkpar,azpar,rindex,saux)
	write(6,*) 'parallel is ok.'

	end

c*****************************************************************

