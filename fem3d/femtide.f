
!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2018  Georg Umgiesser
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

! This program performs tidal analysis of water levels in fem format
! The considered constituents are: 
!  - semi-diurnal: M2, S2, N2, K2, NU2, MU2, L2, T2
!  - diurnal:      K1, O1, P1, Q1, J1, OO1, S1
!  - long period:  MF, MM, SSA, MSM, MSM, SA
!
! It produces two files:
!  - out.tid.fem  containing the tidal signal
!  - out.res.fem  containing the residual signal
!
! revision log :
!
! 12.03.2019	ccf	written from scratch

!*****************************************************************

	program femtide

	use clo
	use elabutil
	use elabtime
	use tide
        use fem_util

	implicit none

	character*80 name,string
	integer np,iunit,iout1,iout2
	integer nvers,lmax,nvar,ntype,nlvdi
	integer nvar0,lmax0,np0
	integer idt
	double precision dtime
	double precision atime,atold,atfirst,atlast,atnew
	integer ierr
	integer nrec,iv,ich,isk,nrecs,l,i,ivar
	integer itype(2)
	integer iformat,iformout
	integer date,time
	integer datetime(2)
	integer np_out
	real regpar(7)
	logical bdtok
	logical breg
	integer, allocatable :: ivars(:)
	character*80, allocatable :: strings(:)
	character*80, allocatable :: strings_out(:)
	character*20 dline,fline
	real,allocatable :: data(:,:,:)
	real,allocatable :: hd(:)
	real,allocatable :: hlv(:)
	integer,allocatable :: ilhkv(:)
	integer,allocatable :: ius(:)

	integer nfile
        integer nequ
        double precision	        :: lat   	  !representative latitude
        double precision, allocatable   :: x(:,:)
	double precision, allocatable   :: acov(:,:,:)
	double precision, allocatable   :: bcov(:,:,:,:)
        double precision, allocatable   :: tvar(:,:,:)
        real, allocatable       	:: tideh(:,:,:,:)
        real, allocatable 	        :: tideg(:,:,:,:)
        real, allocatable    	        :: datat(:,:,:) 	!tidal level
        real, allocatable     	        :: datas(:,:,:) 	!residual

	real				:: dtday	!time series lenght [day]

        nequ=ntd*2+1
        datetime = 0
	dtime = 0.
        nrec = 0
	regpar = 0.

!--------------------------------------------------------------
! set command line options
!--------------------------------------------------------------

        call clo_init('femtide','fem-file','1.0')

        call clo_add_info('performs tidal analysis on fem-file')

        call clo_add_sep('options in/output')
        call clo_add_option('verb',.false.,'be more verbose')
        call clo_add_option('silent',.false.,'be silent')
        call clo_add_option('quiet',.false.,'do not write time records')

        call clo_add_sep('additional options')

        call clo_add_option('lat',45.
     +             ,'latitude for nodal correction (default 45)')

!--------------------------------------------------------------
! parse command line options
!--------------------------------------------------------------

        call clo_parse_options(1)  !expecting (at least) 1 file after options

!--------------------------------------------------------------
! get command line options
!--------------------------------------------------------------

        call clo_get_option('verb',bverb)
        call clo_get_option('quiet',bquiet)
        call clo_get_option('silent',bsilent)
        call clo_get_option('lat',lat)

        if( bsilent ) bquiet = .true.

!--------------------------------------------------------------
! initialize tidal parameters
!--------------------------------------------------------------

	call tidepar_init

!--------------------------------------------------------------
! open file I
!--------------------------------------------------------------

        nfile = clo_number_of_files()
        if( nfile > 1 ) then
          write(6,*) 'Only one file allowed.. exiting'
          stop 'error stop femtide'
        end if

	call clo_get_file(1,infile)

	np = 0
	call fem_file_read_open(infile,np,iformat,iunit)
	if( iunit .le. 0 ) stop

	if( .not. bquiet ) then
	  write(6,*) 'file name: ',infile(1:len_trim(infile))
	  call fem_file_get_format_description(iformat,fline)
	  write(6,*) 'format: ',iformat,"  (",trim(fline),")"
	end if

!--------------------------------------------------------------
! prepare for output
!--------------------------------------------------------------

	iformout = iformat
        iout1 = iunit + 1
        iout2 = iunit + 2

        if( iformout .eq. 1 ) then
          open(iout1,file='out.tid.fem',status='unknown',
     +          form='formatted')
          open(iout2,file='out.res.fem',status='unknown',
     +          form='formatted')
        else
          open(iout1,file='out.tid.fem',status='unknown',
     +          form='unformatted')
          open(iout2,file='out.res.fem',status='unknown',
     +          form='unformatted')
        end if

        date = 0
        time = 0
        call elabtime_date_and_time(date,time)  !we work with absolute time
	call elabtime_set_minmax(stmin,stmax)

!--------------------------------------------------------------
! read first record
!--------------------------------------------------------------

        call fem_file_read_params(iformat,iunit,dtime
     +                        ,nvers,np,lmax,nvar,ntype,datetime,ierr)

	if( ierr .ne. 0 ) goto 99

	if( .not. bquiet ) then
	  write(6,*) 'nvers:  ',nvers
	  write(6,*) 'np:     ',np
	  write(6,*) 'lmax:   ',lmax
	  write(6,*) 'nvar:   ',nvar
	  write(6,*) 'ntype:  ',ntype
	end if

	allocate(hlv(lmax))
	call fem_file_make_type(ntype,2,itype)

	call fem_file_read_2header(iformat,iunit,ntype,lmax
     +			,hlv,regpar,ierr)
	if( ierr .ne. 0 ) goto 98

	if( .not. bquiet ) then
	 if( bverb .and. lmax > 1 ) then
	  write(6,*) 'levels: ',lmax
	  write(6,'(5f12.4)') hlv
	 end if
	 if( itype(1) .gt. 0 ) then
	  write(6,*) 'date and time: ',datetime
	 end if
	 breg = .false.
	 if( itype(2) .gt. 0 ) then
	  breg = .true.
	  write(6,*) 'regpar: '
	  call printreg(regpar)
	 end if
	end if

	call dts_convert_to_atime(datetime,dtime,atime)

	atfirst = atime
	atlast = atime
	atold = atime

	nvar0 = nvar
	lmax0 = lmax
	nlvdi = lmax
	np0 = np
	allocate(ivars(nvar))
	allocate(strings(nvar))
	allocate(strings_out(nvar))
	allocate(data(nlvdi,np,nvar))
	allocate(hd(np))
	allocate(ilhkv(np))

        allocate(x(nequ,np))
        allocate(acov(nequ,nequ,np))
        allocate(bcov(nequ,nvar,nlvdi,np))
        allocate(tvar(nvar,nlvdi,np))
        allocate(tideh(nvar,nlvdi,np,ntd))
        allocate(tideg(nvar,nlvdi,np,ntd))
        allocate(datat(nlvdi,np,nvar))
        allocate(datas(nlvdi,np,nvar))

	datat = 0.
	datas = 0.
	x = 0.
	acov = 0.
	tvar = 0.
	bcov = 0.

!--------------------------------------------------------------
! write info to terminal
!--------------------------------------------------------------

	if( .not. bquiet ) then
          write(6,*) 'available variables contained in file: '
          write(6,*) 'total number of variables: ',nvar
          write(6,*) '   varnum     varid    varname'
	end if

	do iv=1,nvar
	  call fem_file_skip_data(iformat,iunit
     +                          ,nvers,np,lmax,string,ierr)
	  if( ierr .ne. 0 ) goto 97
	  call string2ivar(string,ivar)
	  if( .not. bquiet ) then
	    write(6,'(2i10,4x,a)') iv,ivar,trim(string)
	  end if
	  ivars(iv) = ivar
	  strings(iv) = string
	end do

	strings_out = strings

!--------------------------------------------------------------
! close and re-open file II
!--------------------------------------------------------------

	close(iunit)

	np = 0
	call fem_file_read_open(infile,np,iformat,iunit)
	if( iunit .le. 0 ) stop

!--------------------------------------------------------------
! loop on all records
!--------------------------------------------------------------

	atime = atold
	nrec = 0
	idt = 0
	ich = 0
	isk = 0
        dtday = 0.d0

	do 
	  atold = atime
          call fem_file_read_params(iformat,iunit,dtime
     +                          ,nvers,np,lmax,nvar,ntype,datetime,ierr)
	  if( ierr .lt. 0 ) exit
	  if( ierr .gt. 0 ) goto 99
	  if( nvar .ne. nvar0 ) goto 96
	  if( lmax .ne. lmax0 ) goto 96
	  if( np .ne. np0 ) goto 96
	  nrec = nrec + 1

	  call dts_convert_to_atime(datetime,dtime,atime)
	  call dts_format_abs_time(atime,dline)
	  atlast = atime
	  atnew = atime

	  call fem_file_read_2header(iformat,iunit,ntype,lmax
     +			,hlv,regpar,ierr)
	  if( ierr .ne. 0 ) goto 98

	  do iv=1,nvar
            call fem_file_read_data(iformat,iunit
     +                          ,nvers,np,lmax
     +                          ,string
     +                          ,ilhkv,hd
     +                          ,nlvdi,data(1,1,iv)
     +                          ,ierr)
	    if( ierr .ne. 0 ) goto 97
	    if( string .ne. strings(iv) ) goto 95
	  end do

	  call fem_get_new_atime(iformat,iunit,atnew,ierr)	!peeking
	  if( ierr < 0 ) atnew = atime
	  if( ierr > 0 ) goto 94

          if( elabtime_over_time(atime,atnew,atold) ) exit
          if( .not. elabtime_in_time(atime,atnew,atold) ) cycle

	  bdtok = atime == atfirst .or. atime > atold
	  call check_dt(atime,atold,bcheckdt,nrec,idt,ich,isk)
	  if( .not. bdtok ) cycle
	  atlast = atime

          call set_acov(atime,np,nequ,lat,x,acov)
	  call accumul_tide(nequ,lmax,nvar,np,data,x,tvar,bcov)

	end do

!--------------------------------------------------------------
! finish loop - check Rayleigh criterion and compute tidal constants
!--------------------------------------------------------------

        dtday = (atlast - atfirst) / 86400.

	call rayleigh_crit(dtday)

	call fem_tidana(np,nvar,lmax,nequ,nrec,acov,tvar,bcov,
     +		      tideh,tideg)

        write(*,*)'Constant_name   Amplitude   Phase'
	do i = 1,ntd
          write(*,555)const_ar(i)%name,tideh(1,1,1,i)*100.,
     +       tideg(1,1,1,i)
        end do
 555	format(a5,f12.3,f12.2)

!--------------------------------------------------------------
! close and re-open file III
!--------------------------------------------------------------

        close(iunit)

        call fem_file_read_open(infile,np,iformat,iunit)
        if( iunit .le. 0 ) stop

!--------------------------------------------------------------
! loop on all records
!--------------------------------------------------------------

        do
          atold = atime
          call fem_file_read_params(iformat,iunit,dtime
     +                          ,nvers,np,lmax,nvar,ntype,datetime,ierr)
	  if( ierr .lt. 0 ) exit
	  if( ierr .gt. 0 ) goto 99
	  if( nvar .ne. nvar0 ) goto 96
	  if( lmax .ne. lmax0 ) goto 96
	  if( np .ne. np0 ) goto 96
	  nrec = nrec + 1

	  call dts_convert_to_atime(datetime,dtime,atime)
	  call dts_format_abs_time(atime,dline)
	  atlast = atime
	  atnew = atime

	  call fem_file_read_2header(iformat,iunit,ntype,lmax
     +			,hlv,regpar,ierr)
	  if( ierr .ne. 0 ) goto 98

	  call fem_file_make_type(ntype,2,itype)
	  breg = ( itype(2) .gt. 0 )
	  flag = -999.
	  if( breg ) flag = regpar(7)

          np_out = np
          call fem_file_write_header(iformout,iout1,dtime
     +                        ,0,np_out,lmax,nvar,ntype,lmax
     +                        ,hlv,datetime,regpar)
          call fem_file_write_header(iformout,iout2,dtime
     +                        ,0,np_out,lmax,nvar,ntype,lmax
     +                        ,hlv,datetime,regpar)

          do iv=1,nvar
            call fem_file_read_data(iformat,iunit
     +                      ,nvers,np,lmax
     +                      ,string
     +                      ,ilhkv,hd
     +                      ,nlvdi,data(1,1,iv)
     +                      ,ierr)

            call tidal_pred(np,lmax,
     +          atime,lat,tideh(iv,:,:,:),tideg(iv,:,:,:),
     +          data(1,1,iv),datat,datas)

            call fem_file_write_data(iformout,iout1
     +                          ,0,np_out,lmax
     +                          ,string
     +                          ,ilhkv,hd
     +                          ,nlvdi,datat)

            call fem_file_write_data(iformout,iout2
     +                          ,0,np_out,lmax
     +                          ,string
     +                          ,ilhkv,hd
     +                          ,nlvdi,datas)

          end do

          call fem_get_new_atime(iformat,iunit,atnew,ierr)      !peeking
          if( ierr < 0 ) atnew = atime
          if( ierr > 0 ) goto 94

          if( elabtime_over_time(atime,atnew,atold) ) exit
          if( .not. elabtime_in_time(atime,atnew,atold) ) cycle

          atlast = atime

        end do

!--------------------------------------------------------------
! finish loop - info on time records
!--------------------------------------------------------------

        close(iunit)

        if( .not. bsilent ) then
          write(6,*)
          call dts_format_abs_time(atfirst,dline)
          write(6,*) 'first time record: ',dline
          call dts_format_abs_time(atlast,dline)
          write(6,*) 'last time record:  ',dline

          write(6,*)
          write(6,*) nrec ,' time records read'
          !write(6,*) nelab,' time records elaborated'
          !write(6,*) nout ,' time records written to file'
          write(6,*)
        end if

	if( isk .gt. 0 ) then
	  write(6,*) '*** warning: records eliminated: ',isk
	end if

	write(6,*) 'tidal values written to file:    out.tid.fem'
	write(6,*) 'residual values written to file: out.res.fem'

!--------------------------------------------------------------
! end of routine
!--------------------------------------------------------------

	stop

   94	continue
	write(6,*) 'record: ',nrec+1
	write(6,*) 'cannot peek into header of file'
	stop 'error stop femtide'
   95	continue
	write(6,*) 'strings not in same sequence: ',iv
        write(6,*) string
        write(6,*) strings(iv)
	stop 'error stop femtide: strings'
   96	continue
	write(6,*) 'nvar,nvar0: ',nvar,nvar0
	write(6,*) 'lmax,lmax0: ',lmax,lmax0	!this might be relaxed
	write(6,*) 'np,np0:     ',np,np0	!this might be relaxed
	write(6,*) 'cannot change number of variables'
	stop 'error stop femtide'
   97	continue
	write(6,*) 'record: ',nrec+1
	write(6,*) 'cannot read data record of file'
	stop 'error stop femtide'
   98	continue
	write(6,*) 'record: ',nrec+1
	write(6,*) 'cannot read second header of file'
	stop 'error stop femtide'
   99	continue
	write(6,*) 'record: ',nrec+1
	write(6,*) 'cannot read header of file'
	stop 'error stop femtide'
	end

!*****************************************************************
! Get astronomical arguments from current year

	subroutine set_acov(atime,np,nequ,lat,x,acov)

        use tide

	implicit none

!  input variables
	double precision, intent(in)   :: atime
	integer, intent(in)            :: np
	integer, intent(in)            :: nequ
        double precision, intent(in)   :: lat   !latitude

!  output variables
	double precision, intent(out)   :: x(nequ,np)
	double precision, intent(inout) :: acov(nequ,nequ,np)

!  local variables
	double precision, allocatable   :: v(:)	!ASTRONOMICAL ARGUMENT ADJUSTMENT FOR PHASE
	double precision, allocatable   :: u(:)	!NODAL MODULATION FACTOR FOR PHASE
	double precision, allocatable   :: f(:)	!NODAL MODULATION FACTOR FOR AMPLITUDE
        double precision, allocatable   :: xx(:)
        double precision, allocatable   :: acovl(:,:)	!
	double precision                :: arg
        double precision                :: d    !Universal time (days since 1989-12-31::12)
	double precision                :: ux,vx,fx
	integer		                :: icon,ic,is,i1,i2,k
        double precision, parameter     :: d2r = 4.0*atan(1.0)/180.0  !Degrees to radians
        double precision, parameter     :: twopi=8.d0*atan(1.d0)

	allocate(v(ntd))
	allocate(u(ntd))
	allocate(f(ntd))
	allocate(xx(nequ))
        allocate(acovl(nequ,nequ))

	call tide_time(atime,d)
        call astro(d)
        call setvuf(lat,v,u,f)

        call get_overarray(nequ,v,u,f,xx)

	forall (i1=1:nequ,i2=1:nequ) acovl(i1,i2) = xx(i1)*xx(i2)

	do k = 1,np
          acov(:,:,k) = acov(:,:,k) + acovl(:,:)
          x(:,k) = xx
	end do
 
	end subroutine set_acov

!**********************************************************************
! Acumulate values

	subroutine accumul_tide(nequ,lmax,nvar,np,data,x,
     +			        tvar,bcov)

        implicit none

        integer, intent(in)             :: nequ   !ntd*2+1
        integer, intent(in)             :: lmax   !number of layers
        integer, intent(in)             :: nvar   !number of variables
        integer, intent(in)             :: np     !number of points
        real, intent(in)                :: data(lmax,np,nvar)
	double precision, intent(in)    :: x(nequ,np)

        real, parameter			:: flag = -999.
        double precision, intent(inout) :: tvar(nvar,lmax,np)
        double precision, intent(inout) :: bcov(nequ,nvar,lmax,np)

	integer 	     :: iv,l,i, k
	real		     :: dat

        do iv=1,nvar
          do k = 1,np
            if (data(1,k,iv) == flag ) then
               tvar(iv,:,k) = -999.d0
               bcov(:,iv,:,k) = -999.d0
            else
	      do l = 1,lmax
                dat = data(l,k,iv)
                tvar(iv,l,k) = tvar(iv,l,k) + dat*dat
                do i=1,nequ
                  bcov(i,iv,l,k) = bcov(i,iv,l,k) + x(i,k)*dat
                end do
              end do
            end if
          end do
        end do

	end subroutine accumul_tide

!**********************************************************************
! Compute amplitude and phase for all vars 

	subroutine fem_tidana(np,nvar,lmax,nequ,nrec,acov,tvar,
     +			  bcov,tideh,tideg)

        use tide

        implicit none

	integer, intent(in)  		:: np		!number of points
	integer, intent(in)  		:: nvar		!number of variables
	integer, intent(in) 		:: lmax		!number of layers
	integer, intent(in)  		:: nequ		!ntd*2+1
	integer, intent(in)  		:: nrec		!number of records
        double precision, intent(in)    :: acov(nequ,nequ,np)
        double precision, intent(in)    :: tvar(nvar,lmax,np)
        double precision, intent(in)    :: bcov(nequ,nvar,lmax,np)

        real, intent(out)    		:: tideh(nvar,lmax,np,ntd) !amplitude
        real, intent(out)    		:: tideg(nvar,lmax,np,ntd) !phase

        integer                         :: iv,l,k
        double precision                :: tvarl	!
        double precision, allocatable   :: bcovl(:)	!
        double precision, allocatable   :: acovl(:,:)	!
        real, allocatable               :: h(:)		!amplitude
        real, allocatable               :: g(:)		!phase

	allocate(acovl(nequ,nequ))
	allocate(bcovl(nequ))
	allocate(h(ntd))
	allocate(g(ntd))

        do iv=1,nvar
	  do k = 1,np
            if ( tvar(iv,1,k) == -999.d0 ) then
              tideh(iv,:,k,:) = -999.
              tideg(iv,:,k,:) = -999.
            else
	      do l = 1,lmax
	        tvarl = tvar(iv,l,k)
	        acovl = acov(:,:,k)
	        bcovl = bcov(:,iv,l,k)

	        call tda_solve(nequ,nrec,tvarl,acovl,bcovl,h,g)

	        tideh(iv,l,k,:) = h(:)
	        tideg(iv,l,k,:) = g(:)
              end do
            end if
          end do
        end do

	end subroutine fem_tidana

!**********************************************************************
! General tide prediction subroutine

        subroutine tidal_pred(np,lmax,atime,lat,tideh,tideg,
     +		   data,datat,datas)

	use tide

        implicit none

	integer, intent(in)           :: np	   	    !number of points
	integer, intent(in)           :: lmax	   	    !number of layers
	double precision, intent(in)  :: atime		    !time record
        double precision, intent(in)  :: lat   		    !representative latitude
        real, intent(in)              :: tideh(lmax,np,ntd) !amplitude
        real, intent(in)              :: tideg(lmax,np,ntd) !phase
        real, intent(in)              :: data(lmax,np)      !original data

        real, intent(out)    	      :: datat(lmax,np)     !tidal levels
        real, intent(out)             :: datas(lmax,np)     !residuals

	integer		     :: k,l,ic
	double precision, allocatable :: v(:)	!ASTRONOMICAL ARGUMENT ADJUSTMENT FOR PHASE
	double precision, allocatable :: u(:)	!NODAL MODULATION FACTOR FOR PHASE
	double precision, allocatable :: f(:)	!NODAL MODULATION FACTOR FOR AMPLITUDE
        double precision     :: d      	        !Universal time (days since 1989-12-31::12)
	double precision     :: arg
	real		     :: tp,h,g
        real, parameter	  	     :: flag = -999.
        double precision, parameter  :: d2r = 4.0*atan(1.0)/180.0  !Degrees to radians
        double precision, parameter  :: twopi=8.d0*atan(1.d0)

	allocate(v(ntd))
	allocate(u(ntd))
	allocate(f(ntd))

  	call tide_time(atime,d)
        call astro(d)
        call setvuf(lat,v,u,f)

	do k = 1,np
          if ( tideh(1,k,1) == flag ) then
            datat(:,k) = flag
            datas(:,k) = flag
          else 
            do l = 1,lmax
 	      tp = 0.
              do ic=1,ntd
	        h = tideh(l,k,ic)
	        g = tideg(l,k,ic)*d2r
                arg = (v(ic)+u(ic))*twopi - g
                tp = tp + h*f(ic)*cos(arg)*const_ar(ic)%mask
              end do
	      datat(l,k) = tp
	      datas(l,k) = data(l,k) - tp
            end do
          end if
	end do

	end subroutine tidal_pred

!**********************************************************************
