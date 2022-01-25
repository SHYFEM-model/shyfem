
!--------------------------------------------------------------------------
!
!    Copyright (C) 1993,1997-1999,2002-2004,2008,2010  Georg Umgiesser
!    Copyright (C) 2014-2019  Georg Umgiesser
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

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
c 09.10.2008	ggu	new call to confop
c 23.03.2010	ggu	changed v6.1.1
c 20.01.2014	ggu	new calls to ous routines
c 28.01.2014	ggu	changed VERS_6_1_71
c 26.11.2014	ggu	changed VERS_7_0_7
c 19.12.2014	ggu	changed VERS_7_0_10
c 23.12.2014	ggu	changed VERS_7_0_11
c 19.01.2015	ggu	changed VERS_7_1_3
c 05.06.2015	ggu	changed VERS_7_1_12
c 10.07.2015	ggu	changed VERS_7_1_50
c 13.07.2015	ggu	changed VERS_7_1_51
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 18.09.2015	ggu	changed VERS_7_2_3
c 23.09.2015	ggu	changed VERS_7_2_4
c 28.04.2016	ggu	changed VERS_7_5_9
c 17.05.2016	ggu	in ts_shell() save isvect/itvect
c 25.05.2016	ggu	changed VERS_7_5_10
c 31.10.2016	ggu	new output format in rms files
c 12.01.2017	ggu	changed VERS_7_5_21
c 03.04.2018	ggu	changed VERS_7_5_43
c 11.05.2018	ggu	changed VERS_7_5_47
c 14.02.2019	ggu	changed VERS_7_5_56
c 16.02.2019	ggu	changed VERS_7_5_60
c 15.07.2021	ggu	do not use -> it must be converted to dtime
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
c local
	character*80 nam,dir,file
	integer ierr,ii,ie,k,id
        integer nvers,lmax,l
	integer date,time
	integer ftype,nvar
	real href,hzoff,rr,hm
	double precision dtime
	character*80 title,femver
	character*20 aline
c function
	integer iround
	integer wfout,wrout,ifileo
	real getpar
	double precision dgetpar
	integer ifemop
	logical has_output_d,is_over_output_d,next_output_d
c save
	real, save, allocatable :: ur(:,:)
	real, save, allocatable :: vr(:,:)
	real, save, allocatable :: znr(:)
	real, save, allocatable :: zer(:,:)

	integer, save :: icall = 0
	integer, save :: nr = 0
	logical, save :: bfirst = .true.
	logical, save :: bdebug = .false.
	double precision, save :: da_out(4)

	if(icall.eq.-1) return

c---------------------------------------------------------------------
c initialize routine and variables
c---------------------------------------------------------------------

	if(icall.eq.0) then
          call init_output_d('itmres','idtres',da_out)
	  call increase_output_d(da_out)	!itres=itmres+idtres
          if( .not. has_output_d(da_out) ) icall = -1
	  if(icall.eq.-1) return

	  call shyfem_init_hydro_file('resid',.false.,id)
          da_out(4) = id

	  bdebug = iround(getpar('levdbg')) .ge. 1

	  if( bdebug ) then
	    call get_act_timeline(aline)
	    write(6,*) 'rmsvel : res-file opened ',aline
	  end if

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

	if( .not. is_over_output_d(da_out) ) return	!before start of accum

	if( bfirst ) then
	   call get_act_timeline(aline)
	   write(6,*) 'resid: starting summing ',aline
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

	if( .not. next_output_d(da_out) ) return

c---------------------------------------------------------------------
c write results into file
c---------------------------------------------------------------------

	if( bfirst ) then
	   call get_act_timeline(aline)
	   write(6,*) 'resid: res-file written ',aline,nr
	end if

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

          id = nint(da_out(4))
          call get_act_dtime(dtime)
          call shy_write_hydro_records(id,dtime,nlvdi,znr,zer
     +                                  ,ur,vr)

c---------------------------------------------------------------------
c reset variables
c---------------------------------------------------------------------

	  nr=0
	  ur = 0.
	  vr = 0.
	  znr = 0.
	  zer = 0.

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
c local
	double precision rr,dtime
	logical bout
	integer ii,ie,k,id,idc
	integer l,lmax
	real hm,u,v
	real v1v(nkn),v2v(nkn)
	real rmse(nlvdi,nel)
	real rmsn(nlvdi,nkn)
	real aux(nlvdi,nkn)
	character*20 aline
c function
	real getpar
	logical has_output_d,is_over_output_d,next_output_d
c save
	double precision, save, allocatable :: rms(:,:)
	logical, save :: bdebug
	double precision, save :: da_out(4)
	integer, save :: icall = 0
	integer, save :: nr = 0

	if(icall.eq.-1) return

c-----------------------------------------------------------
c initialization
c-----------------------------------------------------------

	if(icall.eq.0) then
          call init_output_d('itmrms','idtrms',da_out)
	  call increase_output_d(da_out)	!itres=itmres+idtres

	  if( .not. has_output_d(da_out) ) icall = -1
	  if(icall.eq.-1) return

	  call shyfem_init_scalar_file('rms',1,.false.,id)
	  da_out(4) = id

	  bdebug = nint(getpar('levdbg')) .ge. 1
	  if( bdebug ) then
	    call get_act_timeline(aline)
	    write(6,*) 'rmsvel : rms-file opened ',aline
	  end if

	  allocate(rms(nlvdi,nel))
	  nr=0
	  rms = 0.
	end if

c-----------------------------------------------------------
c see if we have to start summing
c-----------------------------------------------------------

	icall=icall+1

	bout = is_over_output_d(da_out)
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

	bout = next_output_d(da_out)
	if( bout ) then
	  if( bdebug ) then
	    call get_act_timeline(aline)
	    write(6,*) 'rmsvel : rms-file written ',aline,nr
	  end if

	  rr=1./nr
	  rmse = sqrt(rms*rr)

	  call e2n3d(nlvdi,rmse,rmsn,aux)

	  call get_act_dtime(dtime)
	  idc = 18
          id = nint(da_out(4))
	  call shy_write_scalar_record(id,dtime,idc,nlvdi,rmsn)
	  call shy_sync(id)

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
	integer itanf,itend
	double precision dtime
	real high

	high = 1.e+30

	bdebug = .false.
	bdebug = .true.

c-------------------------------------------------------------
c initialize parameter array
c-------------------------------------------------------------

	ivect = 0

c-------------------------------------------------------------
c check levels
c-------------------------------------------------------------

	stop 'error stop cmed_init: not ready for dtime'

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

	call get_first_dtime(dtime)
	itanf = nint(dtime)
	call get_last_dtime(dtime)
	itend = nint(dtime)

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
	integer idtc,itmc,itc,it
	integer i,k,l,nlev
	real high
	real c
	double precision rr
	double precision dtime
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

	call get_act_dtime(dtime)
	it = nint(dtime)

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

