c
c $Id: bedload.f,v 1.4 2009-04-21 10:38:59 georg Exp $
c
c finite element model ht (version 3D)
c
c 24.11.2008    ggu	written exner to test bedload
c
c*****************************************************************

c----------------------------------------------------------------

	program bedload

	include 'param.h'

c----------------------------------------------------------------

	parameter (ibndim=100)
	parameter (mardim=nlvdim*10)

c variables and coefficients

	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /mkonst/ eps1,eps2,pi,flag,high,hihi
	common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft
	common /femtim/ itanf,itend,idt,nits,niter,it

	common /level/ nlvdi,nlv

c run and basin description

	character*80 descrp,descrr
	common /descrp/ descrp
	common /descrr/ descrr

c boundary file names			!$$ST	!$$DESCRP

	character*80 boundn(nbcdim)
	character*80 conzn(nbcdim)
	character*80 saltn(nbcdim)
	character*80 tempn(nbcdim)
        character*80 bio2dn(nbcdim)
        character*80 sed2dn(nbcdim)
        character*80 tox3dn(nbcdim)
	character*80 bfm1bc(nbcdim)
        character*80 bfm2bc(nbcdim)
        character*80 bfm3bc(nbcdim)

	common /boundn/ boundn
        common /conzn/ conzn
        common /saltn/ saltn
        common /tempn/ tempn
        common /bio2dn/ bio2dn
        common /sed2dn/ sed2dn
        common /tox3dn/ tox3dn
        common /bfm1bc/bfm1bc
        common /bfm2bc/bfm2bc
        common /bfm3bc/bfm3bc
  
c various arrays

	common /knausc/ knausm,knaus(nexdim)

        common /kfluxc/ nsect,kfluxm,kflux(nfxdim)
        common /iflux/ iflux(3,nfxdim)

        integer nvols,kvold,kvolm,kvol(nfxdim)
        common /kvolc/ nvols,kvold,kvolm,kvol
        integer ivolm,ivol(nfxdim)
        common /ivol/ivolm,ivol

c basin arrays

	common /nen3v/nen3v(3,neldim)
	common /iarv/iarv(neldim)
	common /ipv/ipv(nkndim), /ipev/ipev(neldim)
	common /xgv/xgv(nkndim), /ygv/ygv(nkndim)
	common /hm3v/hm3v(3,neldim)

c static geometry information

	include 'evmain.h'

        common /ilinkv/ilinkv(nkndim+1)
        common /lenkv/lenkv(nlkdim)
        common /linkv/linkv(nlkdim)

	common /ieltv/ieltv(3,neldim)
	common /kantv/kantv(2,nkndim)
	common /dxv/dxv(nkndim), /dyv/dyv(nkndim)

c dynamic geometry information

	common /iwegv/iwegv(neldim)
        common /inodv/inodv(nkndim)

c boundary arrays

	common /bnd/bnd(ibndim,nbcdim)
	common /ierv/ierv(2,nrbdim)
	common /rhv/rhv(nrbdim), /rlv/rlv(nrbdim)
	common /rrv/rrv(nrbdim), /irv/irv(nrbdim)
	common /rzv/rzv(nkndim), /rqv/rqv(nkndim)
	common /iopbnd/iopbnd(nkndim)

        real rqpsv(nkndim), rqdsv(nkndim)
        common /rqpsv/rqpsv, /rqdsv/rqdsv
        real evapv(nkndim)
        common /evapv/evapv
        real mfluxv(nlvdim,nkndim)
        common /mfluxv/mfluxv

c depth structure of levels

	common /ilhv/ilhv(neldim)
	common /ilhkv/ilhkv(nkndim)
	common /hlv/hlv(nlvdim), /hldv/hldv(nlvdim)

        integer ilmv(neldim)
        common /ilmv/ilmv
        integer ilmkv(nkndim)
        common /ilmkv/ilmkv

	common /hkv/hkv(nkndim), /hev/hev(neldim)

c new depth and area arrays

	real hdknv(nlvdim,nkndim)
	common /hdknv/hdknv
	real hdkov(nlvdim,nkndim)
	common /hdkov/hdkov

	real hdenv(nlvdim,neldim)
	common /hdenv/hdenv
	real hdeov(nlvdim,neldim)
	common /hdeov/hdeov

        real areakv(nlvdim,nkndim)
        common /areakv/areakv

c water level and velocity arrays

	common /zov/zov(nkndim), /znv/znv(nkndim)
	common /zeov/zeov(3,neldim), /zenv/zenv(3,neldim)	!$$ZEONV

	common /uov/uov(neldim), /vov/vov(neldim)
	common /unv/unv(neldim), /vnv/vnv(neldim)

	common /ulov/ulov(nlvdim,neldim)
	common /ulnv/ulnv(nlvdim,neldim)
	common /vlov/vlov(nlvdim,neldim)
	common /vlnv/vlnv(nlvdim,neldim)
	common /wlov/wlov(0:nlvdim,nkndim)
	common /wlnv/wlnv(0:nlvdim,nkndim)

	common /utlov/utlov(nlvdim,neldim)
	common /utlnv/utlnv(nlvdim,neldim)
	common /vtlov/vtlov(nlvdim,neldim)
	common /vtlnv/vtlnv(nlvdim,neldim)

	common /uprv/uprv(nlvdim,nkndim)
	common /vprv/vprv(nlvdim,nkndim)
	common /upro/upro(nlvdim,nkndim)
	common /vpro/vpro(nlvdim,nkndim)
	common /wprv/wprv(0:nlvdim,nkndim)

	common /up0v/up0v(nkndim)
	common /vp0v/vp0v(nkndim)
        
        common /fxv/fxv(nlvdim,neldim)		!new HYDRO debora
        common /fyv/fyv(nlvdim,neldim)

	common /xv/xv(3,nkndim)

c concentration, salinity and temperature

	common /saltv/saltv(nlvdim,nkndim)
	common /tempv/tempv(nlvdim,nkndim)

	common /rhov/rhov(nlvdim,nkndim)
	common /bpresv/bpresv(nlvdim,nkndim)
        common /bpresxv/bpresxv(nlvdim,neldim)	!deb110407	debora
        common /bpresyv/bpresyv(nlvdim,neldim)	!deb110407

c coriolis parameter

	common /fcorv/fcorv(neldim)

c friction and diffusion

	common /czv/czv(neldim)
	common /austv/austv(neldim)			!$$AUST

	common /wdifhv/wdifhv(3,3,neldim)   	!weights for horizontal diff.

	common /difhv/difhv(nlvdim,neldim)   	!horizontal diffusion - 3D

        real rdistv(nkndim)
        common /rdistv/rdistv

	common /visv/visv(0:nlvdim,nkndim)	!viscosity (momentum)
	common /difv/difv(0:nlvdim,nkndim)	!diffusivity (scalars)

c special boundary arrays

	common /ruv/ruv(nkndim), /rvv/rvv(nkndim)	!momentum input (2D)
	common /crad/crad(neldim)			!$$GWI (radiation)

c meteo (wind and pressure)

	common /tauxnv/tauxnv(nkndim), /tauynv/tauynv(nkndim)
	common /wxov/wxov(nkndim), /wyov/wyov(nkndim)
	common /wxnv/wxnv(nkndim), /wynv/wynv(nkndim)

	common /ppv/ppv(nkndim)
	common /pov/pov(nkndim), /pnv/pnv(nkndim)

c tidal potential

        real xgeov(nkndim), ygeov(nkndim)
        common /xgeov/xgeov, /ygeov/ygeov
        real zeqv(nkndim)
        common /zeqv/zeqv

c wave sub-module

        real waveh(nkndim)      !wave height [m]
        real wavep(nkndim)      !wave period [s]
        real waved(nkndim)      !wave direction (same as wind direction)
        real waveov(nkndim)     !orbital velocity
        real stokesx(nkndim)    !stokes velocity x
        real stokesy(nkndim)    !stokes velocity y

        common /waveh/waveh, /wavep/wavep, /waved/waved, /waveov/waveov
        common /stokesx/stokesx, /stokesy/stokesy
	save /waveh/,/wavep/,/waved/,/waveov/,/stokesx/,/stokesy/

        real z0bk(nkndim)                   !bottom roughenss on nodes
        common /z0bk/z0bk
	save /z0bk/

c variables for pipe

	integer ipipe,idcoup
	save ipipe, idcoup

c radiation stress

        real radx(neldim),rady(neldim)
        common /radx/radx,/rady/rady
	save /radx/,/rady/

c global matrix

	common /rmat/rmat(mardim)
	common /rvec/rvec(2*nlvdim)

c auxiliary arrays

	common /v1v/v1v(nkndim), /v2v/v2v(nkndim)
	common /v3v/v3v(nkndim), /v4v/v4v(nkndim)
	common /ve1v/ve1v(neldim)
	common /saux1/saux1(nlvdim,nkndim)
	common /saux2/saux2(nlvdim,nkndim)
	common /saux3/saux3(nlvdim,nkndim)
	common /saux4/saux4(nlvdim,nkndim)

c deleted arrays

c	real*8 urv,vrv,zrv
c	common /urv/urv(neldim), /vrv/vrv(neldim), /zrv/zrv(neldim)

	common /tken/tken(0:nlvdim,nkndim)	!turbulent kinetic energy
	common /eps/eps(0:nlvdim,nkndim)	!dissipation rate
	common /rls/rls(0:nlvdim,nkndim)	!length scale

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%%%%%%% code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	call vers3d	!set version of model

c-----------------------------------------------------------
c dimensions
c-----------------------------------------------------------

	nlvdi=nlvdim

	call cstdim(nkndim,neldim,nrbdim,nbcdim
     +			,mbwdim,ngrdim,nardim,nexdim
     +			,nfxdim,nlkdim)

c-----------------------------------------------------------
c read STR file
c-----------------------------------------------------------

	call cstinit
	call cstfile(nkndim,neldim)

c-----------------------------------------------------------
c check dimensions
c-----------------------------------------------------------

	call sp131a(nkndim,neldim)
	call sp131b(nrbdim,nbcdim)
        call sp131m(mbwdim)
	call sp131g(mardim,2,2*nlvdim,2,0)

	call cstcheck	!this also sets time and coriolis

c	call pritst(1)

c-----------------------------------------------------------
c initialize triangles
c-----------------------------------------------------------

	call set_ev
	call set_geom

c-----------------------------------------------------------
c initialize barene data structures
c-----------------------------------------------------------

	call setweg(-1,n)
	call setnod
	call update_geom	!update ieltv - needs inodv

c-----------------------------------------------------------
c inititialize time independent vertical arrays
c-----------------------------------------------------------

	call init_vertical

c-----------------------------------------------------------
c initialize boundary conditions
c-----------------------------------------------------------

	call sp111(1)           !here zenv, utlnv, vtlnv are initialized

	call handle_sigma_init

c-----------------------------------------------------------
c initialize depth arrays and barene data structure
c-----------------------------------------------------------

	call setweg(0,n)
	call setznv		! -> change znv since zenv has changed

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%%%%% time loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	do while( it .lt. itend )

	  it = it + idt

	  call exner

	end do

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%% end of time loop %%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	call endtime

	end

c*****************************************************************

	subroutine exner

	implicit none

	include 'param.h'
	include 'evmain.h'

	integer itanf,itend,idt,nits,niter,it
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /femtim/ itanf,itend,idt,nits,niter,it
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	real v1v(nkndim),v2v(nkndim),v3v(nkndim),v4v(nkndim)
	common /v1v/v1v, /v2v/v2v, /v3v/v3v, /v4v/v4v
	real hkv(nkndim)
	common /hkv/hkv
	real hev(neldim)
	common /hev/hev
	real xgv(nkndim)
	common /xgv/xgv
	real ygv(nkndim)
	common /ygv/ygv
	integer nen3v(3,neldim)
	common /nen3v/nen3v

	logical bnode,bslope
	integer ie,k,ii,i,n,itout
	real ymin,ymax,y
	real h,dt,hmed,hmin,hmax
	real ut,vt,u,v,uv,z,a,ao
	real b,c,con
	real qbx,qby
	real pi,high
	real t,dx

	real zn(3),f(3),fl(3)
	integer kn(3)

	integer icall
	save icall
	data icall / 0 /

	pi = 4. * atan(1.)
	high = 1.e30
	dx = 100.
	bnode = .true.
	bslope = .true.

	if( icall .eq. 0 ) then

	  ymin = high
	  ymax = -high
	  do k=1,nkn
	    y = ygv(k)
	    ymin = min(ymin,y)
	    ymax = max(ymax,y)
	  end do

	  write(6,*) 'ymin,ymax: ',ymin,ymax

	  do k=1,nkn
	    y = ygv(k)
	    h = 2. + cos((y-ymin)*2*pi/(ymax-ymin))
	    hkv(k) = h
	    v1v(k) = 3. - h
	  end do

	  do ie=1,nel
	    h = 0.
	    do ii=1,3
	      k = nen3v(ii,ie)
	      h = h + hkv(k)
	    end do
	    hev(ie) = h/3.
	  end do

	  t = 0.
	  call write_zeta(t,dx,nkn,v1v)

	end if

	icall = icall + 1

	t = it
	dt = idt
	ut = 0.
	vt = 1.
	a = 1.
	a = 0.05
	itout = 10000

	do k=1,nkn
	  v1v(k) = 0.
	  v2v(k) = 0.
	end do

	do ie=1,nel

	  ao = ev(10,ie)
	  u = ut / hev(ie)
	  v = vt / hev(ie)
	  uv = sqrt(u*u+v*v)
	  qbx = a * u
	  qby = a * v

	  z = 3. - hev(ie)

	  if( bnode ) then

	  con = 0.
	  do ii=1,3
	    k = nen3v(ii,ie)
	    b = ev(3+ii,ie)
	    c = ev(6+ii,ie)
	    h = hkv(k)
	    qbx = a * ut / h
	    qby = a * vt / h
	    con = con + ( b * qbx + c * qby )
	  end do

	  z = z - dt * con
	  hev(ie) = 3. - z

	  else

	  do ii=1,3
	    k = nen3v(ii,ie)
	    b = ev(3+ii,ie)
	    c = ev(6+ii,ie)
	    z = 3. - hkv(k)
	    con = 12. * ao * (qbx*b+qby*c)
	    v1v(k) = v1v(k) + 4.*ao*z + dt * con
	    v2v(k) = v2v(k) + 4. * ao
	  end do

	  end if
	end do

	if( bnode ) then

	do k=1,nkn
	  v1v(k) = 0.
	  hkv(k) = 0.
	end do
	do ie=1,nel
	  ao = ev(10,ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    hkv(k) = hkv(k) + ao * hev(ie)
	    v1v(k) = v1v(k) + ao 
	  end do
	end do
	do k=1,nkn
	  hkv(k) = hkv(k) / v1v(k)
	end do

	else 

	do k=1,nkn
	  v1v(k) = v1v(k) / v2v(k)
	  hkv(k) = 3. - v1v(k)
	end do
	do k=1,5
	  i = nkn+1-k
	  v1v(k) = 0.
	  hkv(k) = 3. - v1v(k)
	  v1v(i) = 0.
	  hkv(i) = 3. - v1v(i)
	end do

	do ie=1,nel
	  h = 0.
	  do ii=1,3
	    k = nen3v(ii,ie)
	    h = h + hkv(k)
	  end do
	  hev(ie) = h/3.
	end do

	end if

c slope limiter

	if( bslope ) then

	if( bnode ) then

	do ie=1,nel
	  hmin = high
	  hmax = -high
	  h = 0.
	  do ii=1,3
	    k = nen3v(ii,ie)
	    hmin = min(hmin,hkv(k))
	    hmax = max(hmax,hkv(k))
	    h = h + hkv(k)
	  end do
	  if( hev(ie) .lt. hmin .or. hev(ie) .gt. hmax ) hev(ie) = h/3.
	end do

	else

	do k=1,nkn
	  v1v(k) = high
	  v2v(k) = -high
	  v3v(k) = 0.
	  v4v(k) = 0.
	end do

	do ie=1,nel
	  ao = ev(10,ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    v1v(k) = min(v1v(k),hev(ie))
	    v2v(k) = max(v2v(k),hev(ie))
	    v3v(k) = v3v(k) + ao * hev(ie)
	    v4v(k) = v4v(k) + ao
	  end do
	end do

	do k=1,nkn
	  hmed = v3v(k) / v4v(k)
	  !if( hkv(k) .lt. v1v(k) ) hkv(k) = hmed
	  !if( hkv(k) .gt. v2v(k) ) hkv(k) = hmed
	  if( hkv(k) .lt. v1v(k) .or. hkv(k) .gt. v2v(k) ) then
	    write(6,*) 'slope: ',k,hkv(k),v1v(k),v2v(k),hmed
	    hkv(k) = hmed
	  end if
	end do

	do ie=1,nel
	  h = 0.
	  do ii=1,3
	    k = nen3v(ii,ie)
	    h = h + hkv(k)
	  end do
	  hev(ie) = h/3.
	end do

	end if

	end if

	do k=1,nkn
	  v1v(k) = 3. - hkv(k)
	end do

	if( mod(it,itout) .eq. 0 ) then
	  call write_zeta(t,dx,nkn,v1v)
	end if

	end

c*****************************************************************

	subroutine write_zeta(t,dx,n,zeta)

	implicit none

	include 'param.h'

	real t
	real dx
	integer n
	real zeta(0:n)

	real ygv(nkndim)
	common /ygv/ygv

	integer i,it
	real x
	character*8 time
	character*50 file

	it = nint(t)
	write(time,'(i8)') it

	do i=1,8
	  if( time(i:i) .eq. ' ' ) time(i:i) = '_'
	end do

	file = 'bedload' // time // '.dat'
	write(6,*) 'writing file: ',file
	write(6,*) n,t,dx

	open(1,file=file)

	do i=3,n,5
	  x = ygv(i)
	  write(1,*) x,zeta(i)
	end do

	close(1)

	end

c********************************************************************



