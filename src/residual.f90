!
! $Id: subres.f,v 1.24 2008-10-10 09:29:54 georg Exp $
!
! residual currents
!
! contents :
!
! subroutine resid		computes residual currents
! subroutine rmsvel		computes rms currents
! subroutine dischar		computes discharge (old)
! function resi(zov,znv,n)	computes residuum
! subroutine tsmed 		computes average of scalar values
!
! revision log :
!
! 21.12.1993	ggu	written residual currents
! 15.05.1997	ggu	$$IDTRMS - new parameters for RMS velocity
! 21.05.1997	ggu	$$AVERZR - average zr for dry/wet algorithm
! 21.05.1997	ggu	$$HREFBUG - bug in getting href, hzoff
! 16.06.1997	ggu	$$RESCU - not used subroutine rescu deleted
! 28.04.1998	ggu	dischar deleted
! 20.05.1998	ggu	use ifileo to open file
! 26.01.1999	ggu	use 3D NOS routines for rms velocity
! 18.10.1999	ggu	function resi copied from newres
! 08.10.2002	ggu	subroutine tsmed to compute average of T/S
! 10.10.2002	ggu	subroutine tsmed to compute min/max of T/S
! 19.08.2003	ggu	some cleanup in subroutine rmsvel
! 19.08.2003	ggu	new routines cmed_init, cmed_accum, ts_shell
! 04.03.2004	ggu	bug fix in cmed_accum() for array with more variables
! 10.08.2004	ggu	cmed_init, cmed_accum adjusted also for 2D
! 25.11.2004	ggu	resid converted to 3D (not tested)
! 09.10.2008    ggu     new call to confop
! 20.01.2014    ggu     new calls to ous routines
!
!********************************************************************
!--------------------------------------------------------------------
        module residual
!--------------------------------------------------------------------
        contains
!--------------------------------------------------------------------

	subroutine resid
!
! computes residual currents
!
	use depth
	use hydro_admin
	use levels
        use defnames
        use utility
        use para
	use basin, only : nkn,nel,nkndi,neldi,ngr,mbw
        use output
        use version
        use ioous
        use shympi
        use mpi_io_admin

	implicit none
!
! parameter
        include 'param.h'
! common
	include 'simul.h'
	include 'femtime.h'


! local
	character*80 nam,dir,file
	integer ierr,ii,ie,k
        integer nvers,lmax,l
	integer date,time
	integer nout
	double precision href,hzoff,rr,hm
	character*80 title,femver
! function
	integer wfout,wrout,ifileo
! save
	double precision, save, allocatable :: ur(:,:)
	double precision, save, allocatable :: vr(:,:)
	double precision, save, allocatable :: znr(:)
	double precision, save, allocatable :: zer(:,:)

	integer icall,nr
	save icall,nr
	data icall,nr /0,0/
	integer ia_out(4)
	save ia_out
	logical bfirst,bdebug
	save bfirst,bdebug
	data bfirst /.true./

	if(icall.eq.-1) return

!---------------------------------------------------------------------
! initialize routine and variables
!---------------------------------------------------------------------

	if(icall.eq.0) then
          call init_output('itmres','idtres',ia_out)
	  call increase_output(ia_out)	!itres=itmres+idtres
          if( .not. has_output(ia_out) ) icall = -1

	  if(icall.eq.-1) return

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
          call ous_write_header(nout,nkndi,neldi,nlvdi,ierr)
          if(ierr.gt.0) goto 78
          call ous_write_header2(nout,ilhv,hlv,hev,ierr)
          if(ierr.gt.0) goto 75

	  if( bdebug ) write(6,*) 'resid : res-file opened ',it

	  allocate(ur(nlvdi,nel))
	  allocate(vr(nlvdi,nel))
	  allocate(znr(nkn))
	  allocate(zer(3,nel))
	  nr=0
	  ur = 0.
	  vr = 0.
	  znr = 0.
	  zer = 0.
	end if

!---------------------------------------------------------------------
! already ready for adding?
!---------------------------------------------------------------------

	icall=icall+1

	if( .not. is_over_output(ia_out) ) return	!before start of accum

	if( bfirst ) then
	   write(6,*) 'resid: starting summing ',it
	   bfirst = .false.
	end if

!---------------------------------------------------------------------
! sum transports
!---------------------------------------------------------------------

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

!---------------------------------------------------------------------
! is it time to write file ?
!---------------------------------------------------------------------

	if( .not. next_output(ia_out) ) return

!---------------------------------------------------------------------
! write results into file
!---------------------------------------------------------------------

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


        ierr=0
        if(shympi_partition_on_elements()) then
          call rebuild_structures(ilhv,znr,zer,ur,vr)
          if(shympi_is_master()) then
             call ous_write_record(nout,it,nlvdi,outIlhv,outZnv,outZenv,outUtlnv,outVtlnv,ierr)
             if(ierr.ne.0.) goto 79
          end if
        else
          call ous_write_record(nout,it,nlvdi,ilhv,znr,zer,ur,vr,ierr)
          if(ierr.ne.0.) goto 79
        end if


!---------------------------------------------------------------------
! reset variables
!---------------------------------------------------------------------

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

!---------------------------------------------------------------------
! end of routine
!---------------------------------------------------------------------

	return
   97   continue
	write(6,*) 'Error opening residual file : '
	write(6,*) file
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
!
!********************************************************************
!
	subroutine rmsvel
!
! computes rms currents
!
	use hydro_baro
	use hydro_admin
	use basin
        use utility
        use para
        use output
        use conz_util

	implicit none
!
! parameter
	include 'param.h'
! common
	include 'simul.h'
	include 'femtime.h'
! local
	double precision rr
	integer ii,ie,k,id
	double precision hm,u,v
	double precision v1v(nkn),v2v(nkn)
! save
	double precision, save, allocatable :: rms(:)
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

	  bdebug = iround(getpar('levdbg')) .ge. 1

	  call open_scalar_file(ia_out,1,1,'rms')

	  if( bdebug ) write(6,*) 'rmsvel : rms-file opened ',it

	  allocate(rms(nel))
	  nr=0
	  rms = 0.
	end if

	icall=icall+1

	if( .not. is_over_output(ia_out) ) return

!       sum transports

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

!       if time write output to rms file

	if( next_output(ia_out) ) then

	  if( bdebug ) write(6,*) 'rmsvel : rms-file written ',it,nr

	  rr=1./nr

	  do ie=1,nel
	    rms(ie)=sqrt(rms(ie)*rr)
	  end do

	  v1v=0.
	  v2v=0.

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

! rms velocity is in v1v

	  id = 18
	  call write_scalar_file(ia_out,id,1,v1v)

	  nr=0
	  do ie=1,nel
	    rms(ie)=0.
	  end do
	end if

	return
	end

!********************************************************************
!
!	subroutine dischar
!
! computes discharge -> old, use section flux instead
!
! parameter
!	integer npdim,nedim
!	parameter (npdim=83,nedim=1000)
!+++++++++++++++++++++++++++++++++++++++++++++++
!        data node /
!     +           -149,150,173,172,151,-152,0                      !chioggia
!     +          ,-4310,4311,4315,4314,4309,-4308,0                !mala
!     +          ,-2749,2750,2757,2758,2751,-2752,0                !lido
!     +          ,-4358,3293,3371,3370,3269,-3270,0                !treporti
!     +          ,-3742,3741,3785,3229,3247,3244,3176,3177
!     +          ,3178,3256,3259,3020,-3301,0                     !lido-nord
!     +          ,-2173,2172,2332,2333,2409,2261,2260,2418,1890
!     +          ,2014,2010,2008,1932,1933,1934,2053,2024
!     +          ,-2023,0                                         !mala-lido
!     +          ,-4330,4329,921,1138,1139,1133,1125,1123,1122
!     +          ,1121,1426,1403,878,880,886,893,890,772,773,618
!     +          ,-617,0                                          !mala-chio
!     +          /
!+++++++++++++++++++++++++++++++++++++++++++++++
!
!********************************************************************

        function resi(zov,znv,n)

! computes residuum

        implicit none

	double precision resi
        integer n
        double precision zov(n),znv(n)

        integer i
        double precision res,var,epsr
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

!********************************************************************

	subroutine tsmed

! computes average of scalar values

	use ts
	use levels
        use para
        use conz_util
	use basin, only : nkn,nel,ngr,mbw

	implicit none

! parameter

	include 'param.h'

! common
	include 'femtime.h'

! local
	double precision rr
	double precision tact,sact
	integer k,l,nvar,nlev
	integer itsmed
! save
	double precision, save, allocatable :: tacu(:,:)
	double precision, save, allocatable :: sacu(:,:)
	double precision, save, allocatable :: tmin(:,:)
	double precision, save, allocatable :: tmax(:,:)
	double precision, save, allocatable :: smin(:,:)
	double precision, save, allocatable :: smax(:,:)
	double precision, save, allocatable :: saux(:,:)

	integer icall,nout,nr
	integer idtsca,itmsca,itsca
	double precision high
	logical bdebug
	save icall,nout,nr
	save idtsca,itmsca,itsca
	save high
	save bdebug

	data icall / 0 /
	data high  / 1.e+30 /

!-------------------------------------------------------------
! must compute ?
!-------------------------------------------------------------

	if( icall .eq. -1 ) return

!-------------------------------------------------------------
! first call
!-------------------------------------------------------------

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

	  allocate(tacu(nlvdi,nkn))
	  allocate(sacu(nlvdi,nkn))
	  allocate(tmin(nlvdi,nkn))
	  allocate(tmax(nlvdi,nkn))
	  allocate(smin(nlvdi,nkn))
	  allocate(smax(nlvdi,nkn))
	  allocate(saux(nlvdi,nkn))

	  nr=0
	  tacu = 0.
	  sacu = 0.
	  tmin = high
	  tmax = high
	  smin = -high
	  smax = -high

	end if

!-------------------------------------------------------------
! normal call
!-------------------------------------------------------------

	icall = icall + 1

	if( it .le. itmsca ) return

!-------------------------------------------------------------
! accumulate results
!-------------------------------------------------------------

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

	if( it .lt. itsca ) return

!-------------------------------------------------------------
! write output to file
!-------------------------------------------------------------

	if( bdebug ) write(6,*) 'tsmed : tsa file written ',it,nr

	itsca=itsca+idtsca

	rr=1./nr

	saux = tacu * rr
	call confil(nout,itmsca,idtsca,25,nlvdi,saux)
	saux = sacu * rr
	call confil(nout,itmsca,idtsca,26,nlvdi,saux)

	call confil(nout,itmsca,idtsca,31,nlvdi,tmin)
	call confil(nout,itmsca,idtsca,32,nlvdi,tmax)
	call confil(nout,itmsca,idtsca,35,nlvdi,smin)
	call confil(nout,itmsca,idtsca,36,nlvdi,smax)

	nr=0
	tacu = 0.
	sacu = 0.
	tmin = high
	tmax = high
	smin = -high
	smax = -high

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------

	end

!********************************************************************

	subroutine cmed_init(ext,id,nvar,nlvddi,idtc,itmc,cmed,cmin,cmax,ivect)

! computes average of scalar values - initialization
!
! for 2D arrays call with nlvddi = 1
! for 3D arrays call with nlvddi = nlvdi

	use levels
	use basin
        use conz_util

	implicit none

! parameter

	include 'param.h'
	include 'femtime.h'

	character*(*) ext	!extension of file
	integer id		!id number for variables to be written
	integer nvar		!number of variables to be handled
	integer nlvddi		!number of layers (either nlvdi or 1)
	integer idtc		!frequency of file to be written
	integer itmc		!start time for accumulation
	double precision cmed(nlvddi,nkndi,nvar)	!average
	double precision cmin(nlvddi,nkndi,nvar)		!minimum
	double precision cmax(nlvddi,nkndi,nvar)		!maximum
	integer ivect(8)	!info array that is set up

! local
	logical bdebug
	integer i,k,l,nlev,nlvuse
	integer nout,itc,nr
	double precision high

	high = 1.e+30

	bdebug = .false.
	bdebug = .true.

!-------------------------------------------------------------
! initialize parameter array
!-------------------------------------------------------------

	if( bdebug ) then
	  write(6,*) 'cmed_init called... ',ext
	  write(6,*) id,nvar,nlvddi,idtc,itmc,itanf,itend,idt
	end if

	do i=1,8
	  ivect(i) = 0
	end do

!-------------------------------------------------------------
! check levels
!-------------------------------------------------------------

        if( nlvddi .ne. 1 .and. nlvddi .ne. nlvdi ) then
          write(6,*) 'nlvddi,nlvdi: ',nlvddi,nlvdi
          stop 'error stop cmed_init: invalid nlvddi'
        end if

!-------------------------------------------------------------
! initialize parameter array
!-------------------------------------------------------------

	ivect(2) = id
	ivect(3) = nvar
	ivect(5) = idtc
	ivect(6) = itmc
	ivect(8) = nlvddi

	if(itmc.lt.itanf) itmc=itanf
	if(idtc.le.0) return
	if(itmc+idtc.gt.itend) return

	nout = 0
	nlvuse = min(nlvddi,nlv)
        call confop(nout,itmc,idtc,nlvuse,3*nvar,ext)

	if( nout .le. 0 ) then
	  write(6,*) ext,nout,id,nvar,idtc,itmc
	  stop 'error stop cmed_init: error opening file'
	end if

	if( bdebug ) write(6,*) 'cmed_init : file opened ',ext,id,it

	itc=itmc+idtc

!-------------------------------------------------------------
! initialize arrays array
!-------------------------------------------------------------

	nr=0
	do i=1,nvar
	  do k=1,nkn
	    nlev = min(nlvddi,ilhkv(k))
	    do l=1,nlev
	      cmed(l,k,i) = 0.
	      cmin(l,k,i) = high
	      cmax(l,k,i) = -high
	    end do
	  end do
	end do

!-------------------------------------------------------------
! set parameter array
!-------------------------------------------------------------

	ivect(1) = nout
	ivect(2) = id
	ivect(3) = nvar
	ivect(4) = nr
	ivect(5) = idtc
	ivect(6) = itmc
	ivect(7) = itc
	ivect(8) = nlvddi

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------

	end

!********************************************************************

	subroutine cmed_accum(nlvddi,cvec,cmed,cmin,cmax,ivect)

! computes average of scalar values - accumulation and writing
!
! for 2D arrays call with nlvddi = 1
! for 3D arrays call with nlvddi = nlvdi

	use levels
        use conz_util
	use basin, only : nkn,nel,ngr,mbw

	implicit none

! parameter

	include 'param.h'
	include 'femtime.h'

	integer nlvddi		                !number of layers (nlvdi or 1)
	double precision cvec(nlvddi,nkn,1)			!array with concentration
	double precision cmed(nlvddi,nkn,1)	!average
	double precision cmin(nlvddi,nkn,1)			!minimum
	double precision cmax(nlvddi,nkn,1)			!maximum
	integer ivect(8)                	!info array that is set up

! local
	logical bdebug
	integer nout,id
	integer nvar,nr
	integer idtc,itmc,itc
	integer i,k,l,nlev
	double precision high
	double precision c
	double precision rr
	double precision, allocatable :: saux(:,:)

	high = 1.e+30

	bdebug = .false.
	bdebug = .true.

!-------------------------------------------------------------
! get parameters and see what to do
!-------------------------------------------------------------

	nout = ivect(1)
	id   = ivect(2)
	nvar = ivect(3)
	nr   = ivect(4)
	idtc = ivect(5)
	itmc = ivect(6)
	itc  = ivect(7)
	nlev = ivect(8)

	if( nout .le. 0 ) return

        if( nlvddi .ne. nlev ) then
          write(6,*) 'nlvddi,nlev: ',nlvddi,nlev
          stop 'error stop cmed_accum: invalid nlvddi or nlv'
        end if

!	if( bdebug ) write(6,*) it,nout,id,nvar,nr,idtc,itmc,itc,nlv

	if( it .le. itmc ) return

!-------------------------------------------------------------
! accumulate results
!-------------------------------------------------------------

	nr=nr+1
	do i=1,nvar
	  do k=1,nkn
	    nlev = min(nlvddi,ilhkv(k))
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

!-------------------------------------------------------------
! write output to file
!-------------------------------------------------------------

	if( bdebug ) write(6,*) 'cmed_accum : file written ',id,it

	itc=itc+idtc

	rr=1./nr
	allocate(saux(nlvdi,nkn))

	do i=1,nvar
	  saux = cmed(:,:,i) * rr
	  call confil(nout,itmc,idtc,id+1,nlvddi,saux)
	  call confil(nout,itmc,idtc,id+2,nlvddi,cmin(1,1,i))
	  call confil(nout,itmc,idtc,id+3,nlvddi,cmax(1,1,i))
	  id = id + 3
	end do

	deallocate(saux)

!	-------------------------------------------------------------
! 	re-initialize
!	-------------------------------------------------------------

	nr=0
	cmed = 0.
	cmin = high
	cmax = -high

	ivect(4) = nr
	ivect(7) = itc

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------

	end

!********************************************************************

	subroutine ts_shell

	use ts
	use levels, only : nlvdi,nlv
	use basin
        use para

	implicit none

! parameter

	include 'param.h'
! local
	integer idtc,itmc,itsmed
	integer id,nvar
! save
	double precision, save, allocatable :: tacu(:,:)
	double precision, save, allocatable :: sacu(:,:)
	double precision, save, allocatable :: tmin(:,:)
	double precision, save, allocatable :: tmax(:,:)
	double precision, save, allocatable :: smin(:,:)
	double precision, save, allocatable :: smax(:,:)
	integer itvect(8)
	integer isvect(8)

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

	  allocate(tacu(nlvdi,nkn))
	  allocate(sacu(nlvdi,nkn))
	  allocate(tmin(nlvdi,nkn))
	  allocate(tmax(nlvdi,nkn))
	  allocate(smin(nlvdi,nkn))
	  allocate(smax(nlvdi,nkn))

	  tacu = 0.
	  sacu = 0.
	  tmin = 0.
	  tmax = 0.
	  smin = 0.
	  smax = 0.

	  id = 160
	  call cmed_init('tav',id,nvar,nlvdi,idtc,itmc,tacu,tmin,tmax,itvect)
	  id = 170
	  call cmed_init('sav',id,nvar,nlvdi,idtc,itmc,sacu,smin,smax,isvect)

	  icall = 1
	end if

	call cmed_accum(nlvdi,tempv,tacu,tmin,tmax,itvect)
	call cmed_accum(nlvdi,saltv,sacu,smin,smax,isvect)

	end

!********************************************************************

!--------------------------------------------------------------------
        end module residual
!--------------------------------------------------------------------
