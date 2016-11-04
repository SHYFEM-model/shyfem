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
c 17.05.2016    ggu     in ts_shell() save isvect/itvect
c 31.10.2016    ggu     new output format in rms files
c
c********************************************************************
c
	subroutine resid
c
c computes residual currents
c
	use mod_depth
	use mod_hydro
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

c common
	include 'simul.h'
	include 'femtime.h'


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
	real, save, allocatable :: ur(:,:)
	real, save, allocatable :: vr(:,:)
	real, save, allocatable :: znr(:)
	real, save, allocatable :: zer(:,:)

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
          call ous_write_record(nout,it,nlvdi,ilhv,znr,zer
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
	use mod_hydro_baro
	use mod_hydro
	use levels
	use basin
	use mod_layer_thickness

	implicit none

c common
	include 'simul.h'
	include 'femtime.h'
c local
	double precision rr,dtime
	logical bout
	integer ii,ie,k,id,idc
	integer l,lmax
	integer ishyff
	real hm,u,v
	real v1v(nkn),v2v(nkn)
	real rmse(nlvdi,nel)
	real rmsn(nlvdi,nkn)
	real aux(nlvdi,nkn)
c function
	real getpar
	logical has_output,is_over_output,next_output
	logical has_output_d,is_over_output_d,next_output_d
c save
	double precision, save, allocatable :: rms(:,:)
	logical, save :: bdebug
	integer, save :: ia_out(4)
	double precision, save :: da_out(4)
	integer, save :: icall = 0
	integer, save :: nr = 0

	if(icall.eq.-1) return

c-----------------------------------------------------------
c initialization
c-----------------------------------------------------------

	if(icall.eq.0) then
	  ishyff = nint(getpar('ishyff'))

          call init_output('itmrms','idtrms',ia_out)
	  if( ishyff == 1 ) ia_out = 0
	  if( has_output(ia_out) ) then
	    call open_scalar_file(ia_out,nlv,1,'rms')
	    call increase_output(ia_out)	!itres=itmres+idtres
	  end if

          call init_output_d('itmrms','idtrms',da_out)
	  if( ishyff == 0 ) da_out = 0
	  if( has_output_d(da_out) ) then
	    call shyfem_init_scalar_file('rms',1,.false.,id)
	    call increase_output_d(da_out)	!itres=itmres+idtres
	    da_out(4) = id
	  end if

	  bout = has_output(ia_out) .or. has_output_d(da_out)
	  if( .not. bout ) icall = -1
	  if(icall.eq.-1) return

	  bdebug = nint(getpar('levdbg')) .ge. 1
	  bdebug = .true.
	  if( bdebug ) write(6,*) 'rmsvel : rms-file opened ',it

	  allocate(rms(nlvdi,nel))
	  nr=0
	  rms = 0.
	end if

	icall=icall+1

	bout = is_over_output(ia_out) .or. is_over_output_d(da_out)
	if( .not. bout ) return

c-----------------------------------------------------------
c sum transports
c-----------------------------------------------------------

	nr=nr+1
	do ie=1,nel
	  lmax = ilhv(ie)
	  do l=1,lmax
	    hm=hdenv(l,ie)
	    u=utlnv(l,ie)/hm
	    v=vtlnv(l,ie)/hm
	    rms(l,ie) = rms(l,ie) + u*u + v*v
	  end do
	end do

c-----------------------------------------------------------
c if time write output to rms file
c-----------------------------------------------------------

	bout = next_output(ia_out) .or. next_output_d(da_out)
	if( bout ) then

	  if( bdebug ) write(6,*) 'rmsvel : rms-file written ',it,nr

	  idc = 18
	  dtime = t_act

	  rr=1./nr
	  rmse = sqrt(rms*rr)

	  call e2n3d(nlvdi,rmse,rmsn,aux)

c-----------------------------------------------------------
c rms velocity at nodes is in rmsn
c-----------------------------------------------------------

	  if( has_output(ia_out) ) then
	    call write_scalar_file(ia_out,idc,nlvdi,rmsn)
	  end if

          if( has_output_d(da_out) ) then
            id = nint(da_out(4))
	    call shy_write_scalar_record(id,dtime,idc,nlvdi,rmsn)
	  end if

	  nr=0
	  rms = 0.
	end if

c-----------------------------------------------------------
c end of routine
c-----------------------------------------------------------

	end

c********************************************************************

        function resi(zov,znv,n)

c computes residuum

        implicit none

	real resi
        integer n
        real zov(n),znv(n)

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
c********************************************************************
c********************************************************************

	subroutine cmed_init(ext,id,nvar,nlvddi,idtc,itmc
     +                          ,cmed,cmin,cmax,ivect)

c computes average of scalar values - initialization
c
c for 2D arrays call with nlvddi = 1
c for 3D arrays call with nlvddi = nlvdi

	use levels
	use basin

	implicit none

c parameter

	include 'femtime.h'

	character*(*) ext	!extension of file
	integer id		!id number for variables to be written
	integer nvar		!number of variables to be handled
	integer nlvddi		!number of layers (either nlvdi or 1)
	integer idtc		!frequency of file to be written
	integer itmc		!start time for accumulation
	double precision cmed(nlvddi,nkn,nvar)	!average
	real cmin(nlvddi,nkn,nvar)		!minimum
	real cmax(nlvddi,nkn,nvar)		!maximum
	integer ivect(8)	!info array that is set up

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
	  write(6,*) id,nvar,nlvddi,idtc,itmc,itanf,itend,idt
	end if

	ivect = 0

c-------------------------------------------------------------
c check levels
c-------------------------------------------------------------

        if( nlvddi .ne. 1 .and. nlvddi .ne. nlvdi ) then
          write(6,*) 'nlvddi,nlvdi: ',nlvddi,nlvdi
          stop 'error stop cmed_init: invalid nlvddi'
        end if

c-------------------------------------------------------------
c initialize parameter array
c-------------------------------------------------------------

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

c-------------------------------------------------------------
c initialize arrays array
c-------------------------------------------------------------

	nr=0
	cmed = 0.
	cmin = high
	cmax = -high

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
	ivect(8) = nlvddi

c-------------------------------------------------------------
c end of routine
c-------------------------------------------------------------

	end

c********************************************************************

	subroutine cmed_accum(nlvddi,cvec,cmed,cmin,cmax,ivect)

c computes average of scalar values - accumulation and writing
c
c for 2D arrays call with nlvddi = 1
c for 3D arrays call with nlvddi = nlvdi

	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

c parameter

	include 'femtime.h'

	integer nlvddi		                !number of layers (nlvdi or 1)
	real cvec(nlvddi,nkn,*)			!array with concentration
	double precision cmed(nlvddi,nkn,*)	!average
	real cmin(nlvddi,nkn,*)			!minimum
	real cmax(nlvddi,nkn,*)			!maximum
	integer ivect(8)                	!info array that is set up

c local
	logical bdebug
	integer nout,id
	integer nvar,nr
	integer idtc,itmc,itc
	integer i,k,l,nlev
	real high
	real c
	double precision rr
	real, allocatable :: saux(:,:)

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
	nlev = ivect(8)

	if( nout .le. 0 ) return

        if( nlvddi .ne. nlev ) then
          write(6,*) 'nlvddi,nlev: ',nlvddi,nlev
          stop 'error stop cmed_accum: invalid nlvddi or nlv'
        end if

c	if( bdebug ) write(6,*) it,nout,id,nvar,nr,idtc,itmc,itc,nlv

	if( it .le. itmc ) return

c-------------------------------------------------------------
c accumulate results
c-------------------------------------------------------------

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

c-------------------------------------------------------------
c write output to file
c-------------------------------------------------------------

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

c	-------------------------------------------------------------
c 	re-initialize
c	-------------------------------------------------------------

	nr=0
	do i=1,nvar
	  cmed(:,:,i) = 0.
	  cmin(:,:,i) = high
	  cmax(:,:,i) = -high
	end do

	ivect(4) = nr
	ivect(7) = itc

c-------------------------------------------------------------
c end of routine
c-------------------------------------------------------------

	end

c********************************************************************

	subroutine ts_shell

	use mod_ts
	use levels, only : nlvdi,nlv
	use basin

	implicit none

c local
	integer idtc,itmc,itsmed
	integer id,nvar
c function
	real getpar
c save
	double precision, save, allocatable :: tacu(:,:)
	double precision, save, allocatable :: sacu(:,:)
	real, save, allocatable :: tmin(:,:)
	real, save, allocatable :: tmax(:,:)
	real, save, allocatable :: smin(:,:)
	real, save, allocatable :: smax(:,:)
	integer, save :: itvect(8)
	integer, save :: isvect(8)

	integer, save :: icall = 0

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
	  call cmed_init('tav',id,nvar,nlvdi,idtc,itmc
     +                          ,tacu,tmin,tmax,itvect)
	  id = 170
	  call cmed_init('sav',id,nvar,nlvdi,idtc,itmc
     +                          ,sacu,smin,smax,isvect)

	  icall = 1
	end if

	call cmed_accum(nlvdi,tempv,tacu,tmin,tmax,itvect)
	call cmed_accum(nlvdi,saltv,sacu,smin,smax,isvect)

	end

c********************************************************************

