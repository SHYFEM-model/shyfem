c
c $Id: subres.f,v 1.24 2008-10-10 09:29:54 georg Exp $
c
c residual currents
c
c contents :
c
c subroutine resid		computes residual currents
c subroutine rmsvel		computes rms currents
c subroutine dischar		computes discharge (old)
c function resi(zov,znv,n)	computes residuum
c subroutine tsmed 		computes average of scalar values
c
c revision log :
c
c 21.12.1993	ggu	written residual currents
c 15.05.1997	ggu	$$IDTRMS - new parameters for RMS velocity
c 21.05.1997	ggu	$$AVERZR - average zr for dry/wet algorithm
c 21.05.1997	ggu	$$HREFBUG - bug in getting href, hzoff
c 16.06.1997	ggu	$$RESCU - not used subroutine rescu deleted
c 28.04.1998	ggu	dischar deleted
c 20.05.1998	ggu	use ifileo to open file
c 26.01.1999	ggu	use 3D NOS routines for rms velocity
c 18.10.1999	ggu	function resi copied from newres
c 08.10.2002	ggu	subroutine tsmed to compute average of T/S
c 10.10.2002	ggu	subroutine tsmed to compute min/max of T/S
c 19.08.2003	ggu	some cleanup in subroutine rmsvel
c 19.08.2003	ggu	new routines cmed_init, cmed_accum, ts_shell
c 04.03.2004	ggu	bug fix in cmed_accum() for array with more variables
c 10.08.2004	ggu	cmed_init, cmed_accum adjusted also for 2D
c 25.11.2004	ggu	resid converted to 3D (not tested)
c 09.10.2008    ggu     new call to confop
c 20.01.2014    ggu     new calls to ous routines
c
c********************************************************************
c
	subroutine resid
c
c computes residual currents
c
	implicit none
c
c parameter
        include 'param.h'
c common
	character*80 descrp
	common /descrp/ descrp
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer itanf,itend,idt,nits,niter,it
	common /femtim/ itanf,itend,idt,nits,niter,it
        integer nlvdi,nlv
        common /level/ nlvdi,nlv

        real utlnv(nlvdim,1),vtlnv(nlvdim,1)
        common /utlnv/utlnv, /vtlnv/vtlnv
        integer ilhv(1)
        common /ilhv/ilhv
        real znv(1)
        common /znv/znv
        real zenv(3,1)
        common /zenv/zenv
        real hlv(1)
        common /hlv/hlv
        real hev(1)
        common /hev/hev

c local
	character*80 nam,dir,file
	integer ierr,ii,ie,k
        integer nvers,lmax,l
	integer date,time
	integer nout
	real href,hzoff,rr,hm
	character*80 title,femver
c function
	integer iround
	integer wfout,wrout,ifileo
	real getpar
	double precision dgetpar
	integer ifemop
	logical has_output,is_over_output,next_output
c save
	real ur(nlvdim,neldim),vr(nlvdim,neldim)
        real znr(nkndim),zer(3,neldim)
	save ur,vr,znr,zer
	integer icall,nr
	save icall,nr
	data icall,nr /0,0/
	integer ia_out(4)
	save ia_out
	logical bfirst,bdebug
	save bfirst,bdebug
	data bfirst /.true./

	if(icall.eq.-1) return

c---------------------------------------------------------------------
c initialize routine and variables
c---------------------------------------------------------------------

	if(icall.eq.0) then
          call init_output('itmres','idtres',ia_out)
	  call increase_output(ia_out)	!itres=itmres+idtres
          if( .not. has_output(ia_out) ) icall = -1

	  if(icall.eq.-1) return

	  if(nkn.gt.nkndim.or.nel.gt.neldim) goto 98

	  bdebug = iround(getpar('levdbg')) .ge. 1

	  nout = ifemop('.res','unformatted','new')
	  if( nout .le. 0 ) goto 97
	  ia_out(4) = nout

          nvers = 2
	  href=getpar('href')		!$$HREFBUG
	  hzoff=getpar('hzoff')		!$$HREFBUG
          date = nint(dgetpar('date'))
          time = nint(dgetpar('time'))
          title = descrp
          call get_shyfem_version(femver)

          call ous_init(nout,nvers)
          call ous_set_title(nout,title)
          call ous_set_date(nout,date,time)
          call ous_set_femver(nout,femver)
	  call ous_set_hparams(nout,href,hzoff)
          call ous_write_header(nout,nkn,nel,nlv,ierr)
          if(ierr.gt.0) goto 78
          call ous_write_header2(nout,ilhv,hlv,hev,ierr)
          if(ierr.gt.0) goto 75

	  if( bdebug ) write(6,*) 'resid : res-file opened ',it

	  nr=0
	  do ie=1,nel
            lmax = ilhv(ie)
            do l=1,lmax
	      ur(l,ie)=0.
	      vr(l,ie)=0.
            end do
	    do ii=1,3
	      zer(ii,ie)=0.
	    end do
	  end do
	  do k=1,nkn
            znr(k) = 0.
	  end do
	end if

c---------------------------------------------------------------------
c already ready for adding?
c---------------------------------------------------------------------

	icall=icall+1

	if( .not. is_over_output(ia_out) ) return	!before start of accum

	if( bfirst ) then
	   write(6,*) 'resid: starting summing ',it
	   bfirst = .false.
	end if

c---------------------------------------------------------------------
c sum transports
c---------------------------------------------------------------------

	nr=nr+1
	do ie=1,nel
          lmax = ilhv(ie)
          do l=1,lmax
	    ur(l,ie)=ur(l,ie)+utlnv(l,ie)
	    vr(l,ie)=vr(l,ie)+vtlnv(l,ie)
          end do
	  do ii=1,3
	    zer(ii,ie)=zer(ii,ie)+zenv(ii,ie)
	  end do
	end do
        do k=1,nkn
	  znr(k)=znr(k)+znv(k)
	end do

c---------------------------------------------------------------------
c is it time to write file ?
c---------------------------------------------------------------------

	if( .not. next_output(ia_out) ) return

c---------------------------------------------------------------------
c write results into file
c---------------------------------------------------------------------

	  if( bdebug ) write(6,*) 'resid : res-file written ',it,nr

	  rr=1./nr

	  do ie=1,nel
            lmax = ilhv(ie)
            do l=1,lmax
	      ur(l,ie)=ur(l,ie)*rr
	      vr(l,ie)=vr(l,ie)*rr
            end do
	    do ii=1,3
	      zer(ii,ie)=zer(ii,ie)*rr
	    end do
	  end do
          do k=1,nkn
	    znr(k)=znr(k)*rr
	  end do

	  nout = ia_out(4)
          call ous_write_record(nout,it,nlvdim,ilhv,znr,zer
     +                                  ,ur,vr,ierr)
          if(ierr.ne.0.) goto 79

c---------------------------------------------------------------------
c reset variables
c---------------------------------------------------------------------

	  nr=0
	  do ie=1,nel
            lmax = ilhv(ie)
            do l=1,lmax
	      ur(l,ie)=0.
	      vr(l,ie)=0.
            end do
	    do ii=1,3
	      zer(ii,ie)=0.
	    end do
	  end do
	  do k=1,nkn
            znr(k) = 0.
	  end do

c---------------------------------------------------------------------
c end of routine
c---------------------------------------------------------------------

	return
   97   continue
	write(6,*) 'Error opening residual file : '
	write(6,*) file
	stop 'error stop : resid'
   98   continue
	write(6,*) 'Error in dimension :'
	write(6,*) 'nkn,nkndim : ',nkn,nkndim
	write(6,*) 'nel,neldim : ',nel,neldim
	write(6,*) 'Please change parameters NKNDIM and NELDIM'
	write(6,*) 'in routine RESID (file SUBRES)'
	stop 'error stop : resid'
   78   continue
	write(6,*) 'Error writing record 1 of residual file: ',ierr
	stop 'error stop : resid'
   75   continue
	write(6,*) 'Error writing record 2 of residual file: ',ierr
	stop 'error stop : resid'
   79   continue
	write(6,*) 'Error writing data record of residual file: ',ierr
	stop 'error stop : resid'
	end
c
c********************************************************************
c
	subroutine rmsvel
c
c computes rms currents
c
	implicit none
c
c parameter
	include 'param.h'
c common
	character*80 descrp
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer itanf,itend,idt,nits,niter,it
	real zenv(3,1),hm3v(3,1)
	real unv(1),vnv(1)
	integer nen3v(3,1)
	real v1v(1),v2v(1)
	common /descrp/ descrp
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /femtim/ itanf,itend,idt,nits,niter,it
	common /zenv/zenv, /hm3v/hm3v
	common /unv/unv, /vnv/vnv
	common /nen3v/nen3v
	common /v1v/v1v, /v2v/v2v
c local
	double precision rr
	integer ii,ie,k,id
	real hm,u,v
c function
	integer iround
c	integer wfnov,wrnov,ifileo
	real getpar
	logical has_output,is_over_output,next_output
c save
	double precision rms(neldim)
	save rms
	logical bdebug
	save bdebug
	integer ia_out(4)
	save ia_out
	integer icall,nr
	save icall,nr
	data icall,nr /0,0/

	if(icall.eq.-1) then
	  return
	else if(icall.eq.0) then
          call init_output('itmrms','idtrms',ia_out)
	  call increase_output(ia_out)	!itres=itmres+idtres
          if( .not. has_output(ia_out) ) icall = -1

	  if(icall.eq.-1) return

	  if(nkn.gt.nkndim.or.nel.gt.neldim) goto 98

	  bdebug = iround(getpar('levdbg')) .ge. 1

	  call open_scalar_file(ia_out,1,1,'rms')

	  if( bdebug ) write(6,*) 'rmsvel : rms-file opened ',it

	  nr=0
	  do ie=1,nel
	    rms(ie)=0.
	  end do
	end if

	icall=icall+1

	if( .not. is_over_output(ia_out) ) return

c       sum transports

	nr=nr+1
	do ie=1,nel
	  hm=0.
	  do ii=1,3
	    hm=hm+hm3v(ii,ie)+zenv(ii,ie)
	  end do
	  hm=3./hm
	  u=unv(ie)*hm
	  v=vnv(ie)*hm
	  rms(ie) = rms(ie) + u*u + v*v
	end do

c       if time write output to rms file

	if( next_output(ia_out) ) then

	  if( bdebug ) write(6,*) 'rmsvel : rms-file written ',it,nr

	  rr=1./nr

	  do ie=1,nel
	    rms(ie)=sqrt(rms(ie)*rr)
	  end do

	  do k=1,nkn
	    v1v(k)=0.
	    v2v(k)=0.
	  end do

	  do ie=1,nel
	    do ii=1,3
	      k=nen3v(ii,ie)
	      v1v(k)=v1v(k)+rms(ie)
	      v2v(k)=v2v(k)+1.
	    end do
	  end do

	  do k=1,nkn
	    v1v(k)=v1v(k)/v2v(k)
	  end do

c rms velocity is in v1v

	  id = 18
	  call write_scalar_file(ia_out,id,1,v1v)

	  nr=0
	  do ie=1,nel
	    rms(ie)=0.
	  end do
	end if

	return
   98   continue
	write(6,*) 'Error in dimension :'
	write(6,*) 'nkn,nkndim : ',nkn,nkndim
	write(6,*) 'nel,neldim : ',nel,neldim
	write(6,*) 'Please change parameters NKNDIM and NELDIM'
	write(6,*) 'in routine RMSVEL (file SUBRES)'
	stop 'error stop : rmsvel'
	end

c********************************************************************
c
c	subroutine dischar
c
c computes discharge -> old, use section flux instead
c
c parameter
c	integer npdim,nedim
c	parameter (npdim=83,nedim=1000)
c+++++++++++++++++++++++++++++++++++++++++++++++
c        data node /
c     +           -149,150,173,172,151,-152,0                      !chioggia
c     +          ,-4310,4311,4315,4314,4309,-4308,0                !mala
c     +          ,-2749,2750,2757,2758,2751,-2752,0                !lido
c     +          ,-4358,3293,3371,3370,3269,-3270,0                !treporti
c     +          ,-3742,3741,3785,3229,3247,3244,3176,3177
c     +          ,3178,3256,3259,3020,-3301,0                     !lido-nord
c     +          ,-2173,2172,2332,2333,2409,2261,2260,2418,1890
c     +          ,2014,2010,2008,1932,1933,1934,2053,2024
c     +          ,-2023,0                                         !mala-lido
c     +          ,-4330,4329,921,1138,1139,1133,1125,1123,1122
c     +          ,1121,1426,1403,878,880,886,893,890,772,773,618
c     +          ,-617,0                                          !mala-chio
c     +          /
c+++++++++++++++++++++++++++++++++++++++++++++++
c
c********************************************************************

        function resi(zov,znv,n)

c computes residuum

        implicit none

	real resi
        integer n
        real zov(1),znv(1)

        integer i
        real res,var,epsr
        data epsr /1.e-5/

        res=0.

        do i=1,n
           var=zov(i)
           if(abs(var).gt.epsr) then
                res=res+abs((var-znv(i))/var)
           else if(abs(znv(i)).gt.epsr) then
                res=res+abs((var-znv(i))/znv(i))
           end if
        end do

        resi=res

        end

c********************************************************************

	subroutine tsmed

c computes average of scalar values

	implicit none

c parameter

	include 'param.h'

c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer itanf,itend,idt,nits,niter,it
	common /femtim/ itanf,itend,idt,nits,niter,it
        integer nlvdi,nlv
        common /level/ nlvdi,nlv

	real saux1(nlvdim,nkndim)
	real saux2(nlvdim,nkndim)
	common /saux1/saux1
	common /saux2/saux2
	real tempv(nlvdim,nkndim)
	real saltv(nlvdim,nkndim)
	common /tempv/tempv
	common /saltv/saltv
	integer ilhkv(1)
	common /ilhkv/ilhkv
c local
	double precision rr
	real tact,sact
	integer k,l,nvar,nlev
	integer itsmed
c function
	real getpar
c save
	double precision tacu(nlvdim,nkndim)
	double precision sacu(nlvdim,nkndim)
	real tmin(nlvdim,nkndim), tmax(nlvdim,nkndim)
	real smin(nlvdim,nkndim), smax(nlvdim,nkndim)
	integer icall,nout,nr
	integer idtsca,itmsca,itsca
	real high
	logical bdebug
	save tacu,sacu
	save tmin,tmax
	save smin,smax
	save icall,nout,nr
	save idtsca,itmsca,itsca
	save high
	save bdebug

	data icall / 0 /
	data high  / 1.e+30 /

c-------------------------------------------------------------
c must compute ?
c-------------------------------------------------------------

	if( icall .eq. -1 ) return

c-------------------------------------------------------------
c first call
c-------------------------------------------------------------

	if(icall.eq.0) then

	  bdebug = .false.
	  bdebug = .true.

	  itsmed=nint(getpar('itsmed'))
	  idtsca=nint(getpar('idtcon'))
	  itmsca=nint(getpar('itmcon'))

	  if(itmsca.lt.itanf) itmsca=itanf
	  if( itsmed .le. 0 ) icall=-1
	  if(idtsca.le.0) icall=-1
	  if(itmsca+idtsca.gt.itend) icall=-1

	  if(icall.eq.-1) return

	  nout = 0
	  nvar = 6
          call confop(nout,itmsca,idtsca,nlv,nvar,'tsa')

	  if( bdebug ) write(6,*) 'tsmed : tsa file opened ',it

	  itsca=itmsca+idtsca

	  nr=0
	  do k=1,nkn
	    nlev = ilhkv(k)
	    do l=1,nlev
	      tacu(l,k) = 0.
	      sacu(l,k) = 0.
	      tmin(l,k) = high
	      smin(l,k) = high
	      tmax(l,k) = -high
	      smax(l,k) = -high
	    end do
	  end do

	end if

c-------------------------------------------------------------
c normal call
c-------------------------------------------------------------

	icall = icall + 1

	if( it .le. itmsca ) return

c-------------------------------------------------------------
c accumulate results
c-------------------------------------------------------------

	nr=nr+1
	do k=1,nkn
	  nlev = ilhkv(k)
	  do l=1,nlev
	    tact = tempv(l,k)
	    sact = saltv(l,k)
	    tacu(l,k) = tacu(l,k) + tact
	    sacu(l,k) = sacu(l,k) + sact
	    if( tact .lt. tmin(l,k) ) tmin(l,k) = tact
	    if( sact .lt. smin(l,k) ) smin(l,k) = sact
	    if( tact .gt. tmax(l,k) ) tmax(l,k) = tact
	    if( sact .gt. smax(l,k) ) smax(l,k) = sact
	  end do
	end do

	if( bdebug ) then
	  l = 1
	  k = 1000
	  write(98,*) 'tsmed : ',saltv(l,k),sacu(l,k)/nr
	end if

	if( it .lt. itsca ) return

c-------------------------------------------------------------
c write output to file
c-------------------------------------------------------------

	if( bdebug ) write(6,*) 'tsmed : tsa file written ',it,nr

	if( bdebug ) then
	  l = 1
	  k = 1000
	  write(98,*) 'tsmed : ',sacu(l,k)/nr
	end if

	itsca=itsca+idtsca

	rr=1./nr

	do k=1,nkn
	  nlev = ilhkv(k)
	  do l=1,nlev
	    saux1(l,k) = tacu(l,k) * rr
	    saux2(l,k) = sacu(l,k) * rr
	  end do
	end do

	call confil(nout,itmsca,idtsca,25,nlvdim,saux1)
	call confil(nout,itmsca,idtsca,26,nlvdim,saux2)

	call confil(nout,itmsca,idtsca,31,nlvdim,tmin)
	call confil(nout,itmsca,idtsca,32,nlvdim,tmax)
	call confil(nout,itmsca,idtsca,35,nlvdim,smin)
	call confil(nout,itmsca,idtsca,36,nlvdim,smax)

	nr=0
	do k=1,nkn
	  nlev = ilhkv(k)
	  do l=1,nlev
	    tacu(l,k) = 0.
	    sacu(l,k) = 0.
	    tmin(l,k) = high
	    smin(l,k) = high
	    tmax(l,k) = -high
	    smax(l,k) = -high
	  end do
	end do

c-------------------------------------------------------------
c end of routine
c-------------------------------------------------------------

	return
	end

c********************************************************************

	subroutine cmed_init(ext,id,nvar,nlvdii,idtc,itmc
     +                          ,cmed,cmin,cmax,ivect)

c computes average of scalar values - initialization
c
c for 2D arrays call with nlvdii = 1
c for 3D arrays call with nlvdii = nlvdim

	implicit none

c parameter

	include 'param.h'

	character*(*) ext	!extension of file
	integer id		!id number for variables to be written
	integer nvar		!number of variables to be handled
	integer nlvdii		!number of layers (either nlvdim or 1)
	integer idtc		!frequency of file to be written
	integer itmc		!start time for accumulation
	double precision cmed(nlvdii,nkndim,1)	!average
	real cmin(nlvdii,nkndim,1)		!minimum
	real cmax(nlvdii,nkndim,1)		!maximum
	integer ivect(8)	!info array that is set up

c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer itanf,itend,idt,nits,niter,it
	common /femtim/ itanf,itend,idt,nits,niter,it
        integer nlvdi,nlv
        common /level/ nlvdi,nlv

	integer ilhkv(1)
	common /ilhkv/ilhkv
c local
	logical bdebug
	integer i,k,l,nlev,nlvuse
	integer nout,itc,nr
	real high

	high = 1.e+30

	bdebug = .false.
	bdebug = .true.

c-------------------------------------------------------------
c initialize parameter array
c-------------------------------------------------------------

	if( bdebug ) then
	  write(6,*) 'cmed_init called... ',ext
	  write(6,*) id,nvar,nlvdi,idtc,itmc,itanf,itend,idt
	end if

	do i=1,8
	  ivect(i) = 0
	end do

c-------------------------------------------------------------
c check levels
c-------------------------------------------------------------

        if( nlvdi .ne. 1 .and. nlvdi .ne. nlvdim ) then
          write(6,*) 'nlvdi,nlvdim: ',nlvdi,nlvdim
          stop 'error stop cmed_init: invalid nlvdi'
        end if

c-------------------------------------------------------------
c initialize parameter array
c-------------------------------------------------------------

	ivect(2) = id
	ivect(3) = nvar
	ivect(5) = idtc
	ivect(6) = itmc
	ivect(8) = nlvdi

	if(itmc.lt.itanf) itmc=itanf
	if(idtc.le.0) return
	if(itmc+idtc.gt.itend) return

	nout = 0
	nlvuse = min(nlvdii,nlv)
        call confop(nout,itmc,idtc,nlvuse,3*nvar,ext)

	if( nout .le. 0 ) then
	  write(6,*) ext,nout,id,nvar,idtc,itmc
	  stop 'error stop cmed_init: error opening file'
	end if

	if( bdebug ) write(6,*) 'cmed_init : file opened ',ext,id,it

	itc=itmc+idtc

c-------------------------------------------------------------
c initialize arrays array
c-------------------------------------------------------------

	nr=0
	do i=1,nvar
	  do k=1,nkn
	    nlev = min(nlvdi,ilhkv(k))
	    do l=1,nlev
	      cmed(l,k,i) = 0.
	      cmin(l,k,i) = high
	      cmax(l,k,i) = -high
	    end do
	  end do
	end do

c-------------------------------------------------------------
c set parameter array
c-------------------------------------------------------------

	ivect(1) = nout
	ivect(2) = id
	ivect(3) = nvar
	ivect(4) = nr
	ivect(5) = idtc
	ivect(6) = itmc
	ivect(7) = itc
	ivect(8) = nlvdi

c-------------------------------------------------------------
c end of routine
c-------------------------------------------------------------

	end

c********************************************************************

	subroutine cmed_accum(nlvdi,cvec,cmed,cmin,cmax,ivect)

c computes average of scalar values - accumulation and writing
c
c for 2D arrays call with nlvdi = 1
c for 3D arrays call with nlvdi = nlvdim

	implicit none

c parameter

	include 'param.h'

	integer nlvdi		                !number of layers (nlvdim or 1)
	real cvec(nlvdi,nkndim,1)		!array with concentration
	double precision cmed(nlvdi,nkndim,1)	!average
	real cmin(nlvdi,nkndim,1)		!minimum
	real cmax(nlvdi,nkndim,1)		!maximum
	integer ivect(8)                	!info array that is set up

c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer itanf,itend,idt,nits,niter,it
	common /femtim/ itanf,itend,idt,nits,niter,it

	real saux1(nlvdim,nkndim)
	common /saux1/saux1
	real v3v(nkndim)
	common /v3v/v3v

	integer ilhkv(1)
	common /ilhkv/ilhkv
c local
	logical bdebug
	integer nout,id
	integer nvar,nr
	integer idtc,itmc,itc
	integer i,k,l,nlev,nlv
	real high
	real c
	double precision rr

	high = 1.e+30

	bdebug = .false.
	bdebug = .true.

c-------------------------------------------------------------
c get parameters and see what to do
c-------------------------------------------------------------

	nout = ivect(1)
	id   = ivect(2)
	nvar = ivect(3)
	nr   = ivect(4)
	idtc = ivect(5)
	itmc = ivect(6)
	itc  = ivect(7)
	nlv  = ivect(8)

	if( nout .le. 0 ) return

        if( nlvdi .ne. nlv ) then
          write(6,*) 'nlvdi,nlv: ',nlvdi,nlv
          stop 'error stop cmed_accum: invalid nlvdi or nlv'
        end if

c	if( bdebug ) write(6,*) it,nout,id,nvar,nr,idtc,itmc,itc,nlv

	if( it .le. itmc ) return

c-------------------------------------------------------------
c accumulate results
c-------------------------------------------------------------

	nr=nr+1
	do i=1,nvar
	  do k=1,nkn
	    nlev = min(nlvdi,ilhkv(k))
	    do l=1,nlev
	      c = cvec(l,k,i)
	      cmed(l,k,i) = cmed(l,k,i) + c
	      if( c .lt. cmin(l,k,i) ) cmin(l,k,i) = c
	      if( c .gt. cmax(l,k,i) ) cmax(l,k,i) = c
	    end do
	  end do
	end do

	ivect(4) = nr

	if( it .lt. itc ) return

c-------------------------------------------------------------
c write output to file
c-------------------------------------------------------------

	if( bdebug ) write(6,*) 'cmed_accum : file written ',id,it

	itc=itc+idtc

	rr=1./nr

	do i=1,nvar
	  do k=1,nkn
	    nlev = min(nlvdi,ilhkv(k))
	    do l=1,nlev
	      saux1(l,k) = cmed(l,k,i) * rr     !needed because cmed is real*8
              if( l .eq. 1 ) v3v(k) = saux1(l,k)
	    end do
	  end do
          !if( nlvdi .eq. 1 ) then
	  !  call confil(nout,itmc,idtc,id+1,nlvdi,v3v)
          !else
	  !  call confil(nout,itmc,idtc,id+1,nlvdi,saux1)
          !end if
	  call confil(nout,itmc,idtc,id+1,nlvdi,saux1)
	  call confil(nout,itmc,idtc,id+2,nlvdi,cmin(1,1,i))
	  call confil(nout,itmc,idtc,id+3,nlvdi,cmax(1,1,i))
	  id = id + 3
	end do

c	-------------------------------------------------------------
c 	re-initialize
c	-------------------------------------------------------------

	nr=0
	do i=1,nvar
	  do k=1,nkn
	    nlev = min(nlvdi,ilhkv(k))
	    do l=1,nlev
	      cmed(l,k,i) = 0.
	      cmin(l,k,i) = high
	      cmax(l,k,i) = -high
	    end do
	  end do
	end do

	ivect(4) = nr
	ivect(7) = itc

c-------------------------------------------------------------
c end of routine
c-------------------------------------------------------------

	end

c********************************************************************

	subroutine ts_shell

	implicit none

c parameter

	include 'param.h'

	real tempv(nlvdim,nkndim)
	real saltv(nlvdim,nkndim)
	common /tempv/tempv
	common /saltv/saltv
c local
	integer idtc,itmc,itsmed
	integer id,nvar
c function
	real getpar
c save
	double precision tacu(nlvdim,nkndim)
	double precision sacu(nlvdim,nkndim)
	real tmin(nlvdim,nkndim), tmax(nlvdim,nkndim)
	real smin(nlvdim,nkndim), smax(nlvdim,nkndim)
	integer itvect(8)
	integer isvect(8)

	save tacu,sacu
	save tmin,tmax
	save smin,smax
	save itvect,isvect

	integer icall
	save icall

	data icall / 0 /

	if( icall .lt. 0 ) return

	if( icall .eq. 0 ) then

	  itsmed=nint(getpar('itsmed'))
	  if( itsmed .le. 0 ) then
	    icall = -1
	    return
	  end if

	  idtc=nint(getpar('idtcon'))
	  itmc=nint(getpar('itmcon'))

	  nvar = 1

	  id = 160
	  call cmed_init('tav',id,nvar,nlvdim,idtc,itmc
     +                          ,tacu,tmin,tmax,itvect)
	  id = 170
	  call cmed_init('sav',id,nvar,nlvdim,idtc,itmc
     +                          ,sacu,smin,smax,isvect)

	  icall = 1
	end if

	call cmed_accum(nlvdim,tempv,tacu,tmin,tmax,itvect)
	call cmed_accum(nlvdim,saltv,sacu,smin,smax,isvect)

	end

c********************************************************************

