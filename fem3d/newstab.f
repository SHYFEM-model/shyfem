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
c 20.05.2011    ggu     allow for elimination of elems due to high rstab
c 01.06.2011    ggu     wsink for stability integrated
c 12.07.2011    ggu     new rouitne output_stability()
c 14.07.2011    ggu     new rouitne output_stability_node()
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

	subroutine compute_stability(robs,wsink,rkpar,azpar,rindex,saux)

c computes stability index

	implicit none

        include 'param.h'

	real robs
	real wsink
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
     +          ,ddt,robs,wsink,rkpar,difhv,difv
     +		,difmol,azpar,adpar,aapar
     +          ,rindex,istot,isact,nlvdi,nlv)

c----------------------------------------------------------------
c end of routine
c----------------------------------------------------------------

	end

c*****************************************************************

        subroutine make_stability(dt,robs,wsink,rkpar,rindex,istot,saux)

c gets stability index (if necessary computes it)

        implicit none

	include 'param.h'

	real dt
	real robs
	real wsink
        real rkpar
        real rindex
        integer istot
	real saux(nlvdim,nkndim)

	real azpar
	logical exist_stability

c----------------------------------------------------------------
c see if already computed (only for wsink == 0)
c----------------------------------------------------------------

	if( wsink .eq. 0. .and. exist_stability(rkpar,rindex) ) goto 1
	  
c----------------------------------------------------------------
c compute stability index
c----------------------------------------------------------------

	call getaz(azpar)
	call compute_stability(robs,wsink,rkpar,azpar,rindex,saux)

c----------------------------------------------------------------
c insert stability index (only for wsink == 0)
c----------------------------------------------------------------

	if( wsink .eq. 0. ) call insert_stability(rkpar,rindex)

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

        subroutine info_stability(dt,robs,wsink,rkpar,rindex,istot,saux)

c gets stability index (if necessary computes it)

        implicit none

	include 'param.h'

        real dt
	real robs
	real wsink
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
	call compute_stability(robs,wsink,rkpar,azpar,rindex,saux)
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
c
c after this call the program should abort

        implicit none

        real dt
        real rindex

	call internal_stability(1,dt,rindex)

	end

c**********************************************************************

        subroutine eliminate_stability(rmax)

c eliminates elements with stability index higher than rmax

        implicit none

        real rmax

        real rindex,dt

	dt = 0.
	rindex = rmax
	call internal_stability(2,dt,rindex)

	end

c**********************************************************************

        subroutine internal_stability(mode,dt,rindex)

c computes stability index for hydro timestep (internal)
c
c mode = 0		normal call, compute stability
c mode = 1		error call, compute stability and write error message
c mode = 2		eliminate elements with r>rindex

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

	integer ie,l,lmax,iweg
        real rkpar,azpar,ahpar
	real dindex,aindex,tindex,sindex
	real rmax

	real getpar
	logical is_i_nan

        rkpar = 0.
	azpar = 1.
	ahpar = getpar('ahpar')

	do ie=1,nel
	  do l=1,nlv
	    sauxe1(l,ie) = 0.
	    sauxe2(l,ie) = 0.
	  end do
	end do

	rmax = 1.e+30
	if( mode .eq. 2 ) rmax = rindex
	if( mode .eq. 2 ) then
		write(6,*) 'eliminating rmax: ',rmax
	end if

	call momentum_advective_stability(aindex,sauxe1)
	call momentum_viscous_stability(ahpar,dindex,sauxe2)

        call output_stability(dt,sauxe1,sauxe2)	!in case write to file

	aindex = aindex*dt
	dindex = dindex*dt

	tindex = 0.
	do ie=1,nel
	  lmax = ilhv(ie)
	  iweg = 0
	  do l=1,lmax
	    sindex = sauxe1(l,ie)+sauxe2(l,ie)
	    if( sindex .ge. rmax ) iweg = 1
	    tindex = max(tindex,sindex)
	  end do
	  if( iweg .gt. 0 ) then
	    write(6,*) 'eliminating element for stability: ',ie
	    write(569,*) 'eliminating element for stability: ',ie
	    call check_set_unit(570)
	    call check_elem(ie)
	    call set_element_dry(ie)
	  end if
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

c outputs stability index for hydro timestep (internal) (error handling)

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

	call check_set_unit(6)
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

        subroutine output_stability_node(dt,cwrite)

c outputs stability index for hydro timestep (internal)

        implicit none

	include 'param.h'

        real dt
	real cwrite(nlvdim,nkndim)

        integer nlv,nlvdi
        common /level/ nlvdi,nlv
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nen3v(3,1)
	common /nen3v/nen3v
	integer ilhkv(1)
	common /ilhkv/ilhkv

	real smax(nkndim)
	save smax

	logical bnos
	integer ie,ii,k,l,lmax
	integer ia,id,it
	integer idtsti,itmsti
	real sindex,smin
	logical next_output

	integer icall,iustab,ia_out(4)
	save icall,iustab,ia_out
	data icall,iustab /0,0/

	real getpar

c	idtsti = 3600
c	idtsti = 0
c	itmsti = -1

	if( icall .lt. 0 ) return

	if( icall .eq. 0 ) then
	  idtsti = nint(getpar('idtsti'))
	  itmsti = nint(getpar('itmsti'))
	  if( idtsti .le. 0 ) icall = -1
	  if( icall .lt. 0 ) return
	  call set_output_frequency(itmsti,idtsti,ia_out)
	  call confop(iustab,itmsti,idtsti,1,1,'.stb')
	  !write(6,*) 'ia_out: ',ia_out
	  do k=1,nkn
	    smax(k) = 0.
	  end do
	end if

	do k=1,nkn
	  lmax = ilhkv(k)
	  do l=1,lmax
	    sindex = cwrite(l,k)
	    if( sindex .gt. smax(k) ) smax(k) = sindex
	  end do
	end do

	if( next_output(ia_out) ) then
	  smin = 1./dt
	  do k=1,nkn				!convert to time step
	    if( smax(k) .le. smin ) then
	      smax(k) = dt
	    else
	      smax(k) = 1./smax(k)
	    end if
	  end do
	  idtsti = ia_out(1)	!FIXME
	  itmsti = ia_out(2)
	  call confil(iustab,itmsti,idtsti,778,1,smax)
	  do k=1,nkn
	    smax(k) = 0.
	  end do
	end if

	icall = icall + 1

	end

c*****************************************************************

        subroutine output_stability(dt,sauxe1,sauxe2)

c outputs stability index for hydro timestep (internal)

        implicit none

	include 'param.h'

        real dt
	real sauxe1(nlvdim,neldim)
	real sauxe2(nlvdim,neldim)

        integer nlv,nlvdi
        common /level/ nlvdi,nlv
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nen3v(3,1)
	common /nen3v/nen3v
	integer ilhv(1)
	common /ilhv/ilhv

	real smax(nkndim)
	save smax

	logical bnos
	integer ie,ii,k,l,lmax
	integer ia,id,it
	integer idtsti,itmsti
	real sindex,smin
	logical next_output

	integer icall,iustab,ia_out(4)
	save icall,iustab,ia_out
	data icall,iustab /0,0/

	real getpar

c	idtsti = 3600
c	idtsti = 0
c	itmsti = -1

	if( icall .lt. 0 ) return

	if( icall .eq. 0 ) then
	  idtsti = nint(getpar('idtsti'))
	  itmsti = nint(getpar('itmsti'))
	  if( idtsti .le. 0 ) icall = -1
	  if( icall .lt. 0 ) return
	  call set_output_frequency(itmsti,idtsti,ia_out)
	  call confop(iustab,itmsti,idtsti,1,1,'.sti')
	  !write(6,*) 'ia_out: ',ia_out
	  do k=1,nkn
	    smax(k) = 0.
	  end do
	end if

	do ie=1,nel
	  lmax = ilhv(ie)
	  do l=1,lmax
	    sindex = sauxe1(l,ie)+sauxe2(l,ie)
	    do ii=1,3
	      k = nen3v(ii,ie)
	      if( sindex .gt. smax(k) ) smax(k) = sindex
	    end do
	  end do
	end do

	if( next_output(ia_out) ) then
	  smin = 1./dt
	  do k=1,nkn				!convert to time step
	    if( smax(k) .le. smin ) then
	      smax(k) = dt
	    else
	      smax(k) = 1./smax(k)
	    end if
	  end do
	  idtsti = ia_out(1)	!FIXME
	  itmsti = ia_out(2)
	  call confil(iustab,itmsti,idtsti,779,1,smax)
	  do k=1,nkn
	    smax(k) = 0.
	  end do
	end if

	icall = icall + 1

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
	call compute_stability(0.,0.,rkpar,azpar,rindex,saux)
	write(6,*) 'parallel is ok.'

	end

c*****************************************************************

