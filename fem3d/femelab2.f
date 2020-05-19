
!--------------------------------------------------------------------------
!
!    Copyright (C) 2015-2020  Georg Umgiesser
!    Copyright (C) 2017  Marco Bajo
!    Copyright (C) 2019  Christian Ferrarin
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

! elaborates fem files
!
! revision log :
!
! 14.01.2015	ggu	adapted from feminf
! 20.05.2015	ggu	use bhuman to convert to human readable time
! 05.06.2015	ggu	iextract to extract nodal value
! 05.11.2015	ggu	new option chform to change format
! 04.10.2016	ggu	output flags now similar to shyelab
! 05.10.2016	ggu	allow for expansion of regular grid
! 11.10.2016	ggu	introduced flag for min/max/med computation
! 31.10.2016	ggu	new flag condense (bcondense)
! 16.05.2017	ggu&mbj	better handling of points to extract
! 31.08.2017	ggu	new flag -grd to write grd from fem file
! 14.11.2017	ggu	changed VERS_7_5_36
! 17.11.2017	ggu	changed VERS_7_5_37
! 05.12.2017	ggu	changed VERS_7_5_39
! 22.02.2018	ggu	changed VERS_7_5_42
! 03.04.2018	ggu	changed VERS_7_5_43
! 18.12.2018	ggu	changed VERS_7_5_52
! 14.02.2019	ggu	changed VERS_7_5_56
! 16.02.2019	ggu	changed VERS_7_5_60
! 01.03.2019	ccf	lmax in function of considered var (needed for split)
! 13.03.2019	ggu	changed VERS_7_5_61
! 24.04.2019	ggu	write also fem file for extracted node
! 24.04.2019	ggu	correct regpar if lon > 180
! 21.05.2019	ggu	changed VERS_7_5_62
! 17.07.2019	ggu	changes to custom_elab()
! 22.07.2019	ggu	new routines for handling time step check
! 13.12.2019	ggu	new option checkrain and routine rain_elab()
! 03.03.2020	ggu	do not open out.fem if bextract is true
! 29.03.2020	ggu	with bcondense write also out.txt if possible
! 17.04.2020	ggu	with bsplit also write speed/dir if directional
! 18.05.2020	ggu	completely restructured using module fem_util
!
!******************************************************************

c*****************************************************************
c*****************************************************************
c*****************************************************************

	program femelab2_main
	call femelab2
	end

	subroutine femelab2

c writes info on fem file

	use clo
	use elabutil
	use elabtime
	use fem_util

	implicit none

	character*80 name,string
	integer np,iunit,iout
	integer nvers,lmax,nvar,ntype
	integer nvar0,np0
	integer idt,idtact
	double precision dtime
	double precision atime,atold,atfirst,atlast,atnew
	real dmin,dmax,dmed
	integer ierr
	integer nfile
	integer nrec,iv,nrecs,l,i,ivar,lextr,nout
	integer iformat,iformout
	integer date,time
	integer datetime(2),dateanf(2),dateend(2)
	integer iextract,it
	integer ie,nx,ny,ix,iy
	integer np_out,ntype_out
	real x0,y0,dx,dy,x1,y1
	real regpar(7)
	real xp,yp
	real ffact
	real depth
	logical bfirst,bskip,btskip,bnewstring
	logical bhuman,blayer
	logical bdtok,bextract,bexpand
	logical bread,breg,bterm
	integer, allocatable :: ivars(:)
	character*80, allocatable :: strings(:)
	character*80, allocatable :: strings_out(:)
	character*20 aline
	character*80 line
	!character*80 stmin,stmax
	real,allocatable :: facts(:)
	real,allocatable :: d3dext(:,:)

	type(fem_type) :: ffem,ffem_old,ffem_out
	type(femrec_type) :: frec,frec_out

	integer ifileo

!--------------------------------------------------------------------

	bhuman = .true.		!convert time in written fem file to dtime=0
	blayer = .false.
	blayer = .true.		!write layer structure - should be given by CLO

	iextract = 0

        datetime = 0
	dtime = 0.
        nrec = 0
	regpar = 0.

c--------------------------------------------------------------
c initialization
c--------------------------------------------------------------

	call elabutil_init('FEM','femelab')

	bexpand = regexpand > -1
	bextract = ( snode /= ' ' .or. scoord /= ' ' )
        bnewstring = ( newstring /= ' ' )

c--------------------------------------------------------------
c open file
c--------------------------------------------------------------

        call clo_reset_files
        call clo_get_next_file(infile)
	if( infile .eq. ' ' ) stop

	np = 0
	call femutil_open_for_read(infile,np,ffem,ierr)
	if( ierr .ne. 0 ) stop

c--------------------------------------------------------------
c prepare for output if needed
c--------------------------------------------------------------

        iout = 0
	iformat = ffem%femfile%iformat
	iformout = iformat
	if( bchform ) iformout = 1 - iformat
	if( iformout < 0 ) iformout = iformat

!	boutput is true if we have to write out.fem

        boutput = bout
	boutput = boutput .or. bchform
	boutput = boutput .or. bnewstring
	boutput = boutput .or. bexpand
	boutput = boutput .or. bcondense
	boutput = boutput .or. bextract

        if( boutput ) then
	  call femutil_open_for_write('out.fem',iformout,ffem_out)
          call handle_open_output_file(999,'out.fem')
        end if

        date = 0
        time = 0
        call elabtime_date_and_time(date,time)  !we work with absolute time
	call elabtime_set_minmax(stmin,stmax)

	call elabtime_set_inclusive(binclusive)

c--------------------------------------------------------------
c read first record
c--------------------------------------------------------------

	call femutil_read_record_fem(ffem,ierr,.true.)
	if( ierr .ne. 0 ) goto 99
	frec = ffem%femrec

	if( .not. bquiet ) then
	  write(6,*) 'file name: ',trim(infile)
	  call femutil_info(ffem)
	end if

	np0   = frec%np
	nvar0 = frec%nvar
	nvar = nvar0

	allocate(ivars(nvar))
	allocate(strings(nvar))
	allocate(strings_out(nvar))
	allocate(facts(nvar))

	breg = femutil_is_regular(frec)
	if( bexpand .and. .not. breg ) then
	  stop 'error stop femelab: for expand need regular domain'
	end if
	if( bgrd .and. .not. breg ) then
	  stop 'error stop: for bgrd grid must be regular'
	end if

	call femutil_get_time(ffem,atime)

	atfirst = atime
	atlast = atime
	atold = atime

	if( bextract ) then
	  call handle_extract(breg,bquiet,np0,frec%regpar,iextract)
	end if

c--------------------------------------------------------------
c write info to terminal
c--------------------------------------------------------------

	do iv=1,nvar
	  string = ffem%femrec%strings(iv)
	  call string2ivar(string,ivar)
	  ivars(iv) = ivar
	  strings(iv) = string
	end do

	if( .not. bquiet ) then
          write(6,*) 'available variables contained in file: '
          write(6,*) 'total number of variables: ',nvar
          write(6,*) '   varnum     varid    varname'
	  do iv=1,nvar
	    write(6,'(2i10,4x,a)') iv,ivars(iv),trim(strings(iv))
	  end do
	end if

	strings_out = strings
        if( bnewstring ) then
	  call change_strings(nvar,strings_out,newstring)
        end if
	call set_facts(nvar,facts,factstring)

	if( binfo ) return

c--------------------------------------------------------------
c close and re-open file
c--------------------------------------------------------------

	call femutil_close(ffem)

	np = 0
	call femutil_open_for_read(infile,np,ffem,ierr)
	if( ierr .ne. 0 ) stop

c--------------------------------------------------------------
c loop on all records
c--------------------------------------------------------------

	bread = bwrite .or. bextract .or. boutput
	bread = bread .or. bsplit .or. bgrd .or. bcheck
	bread = bread .or. bcheckrain
	bread = bread .or. bcondense
	bskip = .not. bread

	atime = atfirst
	nrec = 0
	idt = 0

	do 
	  atold = atime
	  ffem_old = ffem
	  call femutil_read_record_fem(ffem,ierr,bskip)
	  if( ierr .lt. 0 ) exit
	  if( ierr .gt. 0 ) goto 99

	  nrec = nrec + 1

	  frec = ffem%femrec
	  np   = frec%np
	  lmax = frec%lmax
	  flag = frec%flag

	  call femutil_get_time(frec,atime)
	  call dts_format_abs_time(atime,aline)
	  atlast = atime
	  atnew = atime

	  call femutil_peek_time(ffem,atnew,ierr)
	  if( ierr < 0 ) atnew = atime
	  if( ierr > 0 ) goto 94

          if( elabtime_over_time(atime,atnew,atold) ) exit
          if( .not. elabtime_in_time(atime,atnew,atold) ) cycle

	  call handle_timestep(atime,bcheckdt,btskip)
	  if( btskip ) cycle
	  atlast = atime

	  if( bverb ) then
            write(6,'(a,i8,f20.2,3x,a20)') 'time : ',nrec,atime,aline
	  end if

	  if( ffem%femrec%bchanged ) then
	    if( breg .neqv. femutil_is_regular(frec) ) goto 93
	    if( bextract .and. np .ne. np0 ) goto 89
	    if( nvar .ne. nvar0 ) goto 96
	    bterm = nrec > 1 .and. .not. bquiet
	    call handle_changed(ffem,ffem_old,bterm)
	    if( allocated(d3dext) ) deallocate(d3dext)
	    allocate(d3dext(lmax,nvar))
	  end if

	  do iv=1,nvar
	    !call custom_elab(frec)
	    if( bcheckrain ) call rain_elab(ffem,flag)
	    string = frec%strings(iv)
	    if( .not. bnewstring .and. string .ne. strings(iv) ) goto 95
	    ffact = facts(iv)
	    if( ffact /= 1. ) then
	      where( frec%data(:,:,iv) /= flag )
	        frec%data(:,:,iv) = frec%data(:,:,iv) * ffact
	      end where
	    end if
	  end do

          if( boutput ) then
	    if( bnewstring ) then
	      frec%strings = strings_out
	    end if
	    if( bexpand ) then
	      call reg_expand_shell(frec,regexpand)
	    end if
	    if( bcondense ) then
	      call fem_condense(lmax,nvar,frec,d3dext)
	      call copy_one_point(0,lmax,nvar,frec,d3dext)
	    end if
	    if( bextract ) then
	      call extract_point(iextract,lmax,nvar,frec,d3dext)
	      call copy_one_point(iextract,lmax,nvar,frec,d3dext)
	    end if
	    ffem_out%femrec = frec
	    call femutil_write_record(ffem_out)
          end if

	  do iv=1,nvar
	    if( bwrite ) then
	      write(6,*) nrec,iv,ivars(iv),trim(strings(iv))
	      do l=1,lmax
        	call minmax_data(l,iv,flag,frec,dmin,dmax,dmed)
	        write(6,1000) 'l,min,aver,max : ',l,dmin,dmed,dmax
 1000	        format(a,i5,3g16.6)
	      end do
	    end if
	  end do

	  if( bsplit ) then
	    call femsplit(ffem)
	  end if

	  if( bextract .or. bcondense ) then
            call write_ts(atime,nvar,lmax,strings,d3dext)
	  end if

	  if( bcheck ) then
	    call fem_check_shell(atime,frec,scheck,bquiet)
	  end if

	  if( bgrd ) then
	    do iv=1,nvar
	      write(6,*) 'writing grd file: ',iv,nrec
	      call make_grd_from_data(iv,nrec,lmax,np
     +			,regpar,frec%data(:,:,iv))
	    end do
	  end if

	end do

c--------------------------------------------------------------
c finish loop - last writes and info on time records
c--------------------------------------------------------------

	if( bcheck ) then	!write final data
	  atime = -1.
	  call fem_check_shell(atime,frec,scheck,bquiet)
	end if

        if( .not. bsilent ) then
          write(6,*)
          call dts_format_abs_time(atfirst,aline)
          write(6,*) 'first time record: ',aline
          call dts_format_abs_time(atlast,aline)
          write(6,*) 'last time record:  ',aline

	  call handle_timestep_last(bcheckdt)

          write(6,*)
          write(6,*) nrec ,' time records read'
          write(6,*)
        end if

	call femutil_close(ffem)
	if( boutput ) call femutil_close(ffem_out)

	if( bcheck ) then	!write final message
	  atime = -2.
	  call fem_check_shell(atime,frec,scheck,bquiet)
	end if

	if( .not. bquiet ) then
          if( boutput .or. bsplit ) then
            write(6,*) 'the following files have been written:'
            call handle_open_output_file(0,' ')
	  end if
	end if

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	return
   89	continue
	write(6,*) 'iectract,np,np0: ',iextract,np,np0
	stop 'error stop femelab: cannot change np'
   91	continue
	write(6,*) 'iectract,np: ',iextract,np
	stop 'error stop femelab: no such node'
   93	continue
	write(6,*) 'record: ',nrec+1
	write(6,*) 'old and new records are not equivalent with breg'
	write(6,*) 'old breg: ',breg
	write(6,*) 'new breg: ',femutil_is_regular(frec)
	stop 'error stop femelab'
   94	continue
	write(6,*) 'record: ',nrec+1
	write(6,*) 'cannot peek into header of file'
	stop 'error stop femelab'
   95	continue
	write(6,*) 'strings have changed: ',iv
        write(6,*) 'old: ',trim(strings(iv))
        write(6,*) 'new: ',trim(string)
	stop 'error stop femelab: strings'
   96	continue
	write(6,*) 'nvar,nvar0: ',nvar,nvar0
	!write(6,*) 'np,np0:     ',np,np0	!this might be relaxed
	write(6,*) 'cannot change number of variables'
	stop 'error stop femelab'
   97	continue
	write(6,*) 'record: ',nrec+1
	write(6,*) 'cannot read data record of file'
	stop 'error stop femelab'
   98	continue
	write(6,*) 'record: ',nrec+1
	write(6,*) 'cannot read second header of file'
	stop 'error stop femelab'
   99	continue
	write(6,*) 'record: ',nrec+1
	write(6,*) 'cannot read header of file'
	stop 'error stop femelab'
	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

        subroutine minmax_data(level,iv,flag,frec
     +				,vmin,vmax,vmed)

	use fem_util

        implicit none

	integer level		!level for which minmax to compute (0 for all)
	integer iv
	real flag
	type(femrec_type) :: frec
	real vmin,vmax,vmed

        integer nlvdi,np
        integer k,l,lmin,lmax,lm,ntot
        real v,high
	double precision vtot

	np = frec%np
	nlvdi = frec%lmax

	lmin = max(1,level)
	lmax = level
	if( level == 0 ) lmax = nlvdi

	ntot = 0
	vtot = 0.
	high = 1.e+30
        vmin = high
        vmax = -high

        do k=1,np
          lm = min(frec%ilhkv(k),lmax)
          do l=lmin,lm
            v = frec%data(l,k,iv)
	    if( v == flag ) cycle
	    ntot = ntot + 1
	    vtot = vtot + v
            vmax = max(vmax,v)
            vmin = min(vmin,v)
          end do
        end do

	if( ntot > 0 ) then
	  vmed = vtot / ntot
	else
	  vmin = 0.
	  vmax = 0.
	  vmed = 0.
	end if

        !write(86,*) 'min/max: ',it,vmin,vmax,vmed

        end

c*****************************************************************

        subroutine write_ts(atime,nvar,lmax,strings,data)

	!use iso8601

	implicit none

	double precision atime
	integer nvar,lmax
	character*(*) strings(nvar)
	real data(lmax,nvar)

	integer, save :: iu = 0

	integer it,ierr,n

	character*20 aline
	character*80, save :: eformat
	character*80 varline,file

	integer ifileo

!	------------------------------------------------------------
!	initialize output
!	------------------------------------------------------------

	if( iu == 0 ) then
	  file = 'out.txt'
	  iu = ifileo(88,file,'form','new')
          n = max(nvar,lmax)
	  write(eformat,'(a,i3,a)') '(a20,',n,'g14.6)'
	  !write(6,*) 'using format: ',trim(eformat)
	  call make_varline(nvar,strings,varline)
	  if( varline /= ' ' ) write(iu,'(a)') trim(varline)
          call handle_open_output_file(iu,file)
	end if

	call dts_format_abs_time(atime,aline)

        if( nvar == 1 ) then
          write(iu,eformat) aline,data(1:lmax,1)
	else
          write(iu,eformat) aline,data(1,1:nvar)	!only surface
        end if
          
        end

c*****************************************************************

	subroutine make_varline(nvar,strings,varline)

	use shyfem_strings

	implicit none

	integer nvar
	character*(*) strings(nvar)
	character*(*) varline

	integer iv,iend
	character*80 name
	character*14 short
        character*2, save :: xy(2)=(/'-x','-y'/)
        integer, save :: idir = 0

        logical has_direction

	varline = ' '
	iend = 0

	do iv=1,nvar
	  name = strings(iv)
	  call strings_get_short_name(name,short)
	  varline = varline(1:iend) // short
	  iend = iend + 14
          if( has_direction(name) ) then
            idir = mod(idir,2) + 1
	    varline = trim(varline) // xy(idir)
          end if
	end do

	varline = '#vars: date_and_time   ' // varline

	end

c*****************************************************************

	subroutine get_ivars(nvar,strings,ivars,ivarorig,ivx,ivy)

! looks for directional variable, changes ivar in speed and dir
! returns the position of x,y in ivx,ivy

	implicit none

	integer nvar
	character*(*) strings(nvar)
	integer ivars(nvar)
	integer ivarorig
	integer ivx,ivy

	integer iv,idir
	integer ivar,ivarspeed,ivardir
	character*1 dir
	character*10 unit
	character*80 string

	idir = 0
	ivarorig = 0

	do iv=1,nvar
	  string = strings(iv)
	  call string_direction_and_unit(string,dir,unit)
	  call string2ivar(string,ivar)
	  call get_direction_ivars(ivar,ivarspeed,ivardir)
	  if( ivarspeed /= 0 ) then			!directional variable
	    idir = idir + 1
	    if( idir > 2 ) goto 98
	    if( dir == 'x' .and. .not. idir == 1 ) goto 99
	    if( dir == 'y' .and. .not. idir == 2 ) goto 99
	    if( dir /= ' ' .and. dir /= 'x' .and. dir /= 'y' ) goto 99
	    if( idir == 1 ) then
	      ivx = iv
	      ivars(iv) = ivarspeed
	    else if( idir == 2 ) then
	      ivy = iv
	      ivars(iv) = ivardir
	    end if
	    ivarorig = ivar
	  else
	    ivars(iv) = ivar
	  end if
	end do

	iv = 0
	if( idir > 0 .and. idir /= 2 ) goto 99

	return
   98	continue
	write(6,*) 'only one directional variable allowed in fem file'
	stop 'error stop get_ivars: directional var'
   99	continue
	write(6,*) 'error in naming of variable: '
	write(6,*) 'idir,iv,dir: ',idir,iv,'  ',dir
	stop 'error stop get_ivars: structure of fem file'
	end

c*****************************************************************

	subroutine femsplit(ffem)

! splits fem file into single variables - directional vars are kept in one file

	use fem_util

	implicit none

	type(fem_type) :: ffem

	integer iformat
	integer np,lmax
	integer ivar,nvar,iv,isub
	double precision atime
	character*80 file,extra,filename,string
	character*1 dir
	character*10 unit
	type(femrec_type) :: frec

	integer, save :: icall = 0
	integer, save :: ivarorig,ivx,ivy
	integer, save, allocatable :: ivars(:)
	character*80, save, allocatable :: strings(:)
	type(fem_type), save, allocatable :: ffemout(:)
	type(fem_type), save :: ffemdirout

	frec = ffem%femrec

	iformat = ffem%femfile%iformat

	nvar = frec%nvar
	lmax = frec%nvar
	np   = frec%np

	if( icall == 0 ) then
	  icall = 1
	  allocate(ffemout(nvar))
	  allocate(ivars(nvar))
	  allocate(strings(nvar))
	  strings = frec%strings
	  call get_ivars(nvar,strings,ivars,ivarorig,ivx,ivy)
	  !write(6,*) 'ivars: ',ivars
	  do iv=1,nvar
	    ivar = ivars(iv)
            call ivar2filename(ivar,filename)
            call strings_pop_direction(filename)
	    call ivar2string(ivar,strings(iv),isub)
	    file = 'out.' // trim(filename) // '.fem'
	    !write(6,*) 'opening file: ',trim(file)
	    call femutil_open_for_write(file,iformat,ffemout(iv))
            call handle_open_output_file(999,file)
	    !ffemout(iv)%femrec = frec
	  end do
	end if

	call femsplit_sd(ffem,ivarorig,ivx,ivy)

	call femutil_get_time(frec,atime)

	do iv=1,nvar
	  if( frec%bchanged ) then
	    ffemout(iv)%femrec = ffem%femrec
	    ffemout(iv)%femrec%nvar = 1
	    ffemout(iv)%femrec%strings(1) = strings(iv)
	  end if
	  call femutil_set_time(ffemout(iv),atime)
	  ffemout(iv)%femrec%data(:,:,1) = ffem%femrec%data(:,:,iv)
	  call femutil_write_record_fem(ffemout(iv))
	end do

	end

c*****************************************************************

	subroutine femsplit_sd(ffem,ivar,ivx,ivy)

! splits directional variables in fem file into speed and direction

	use fem_util

	implicit none

	type(fem_type) :: ffem
	integer ivar,ivx,ivy

	integer iformat
	double precision atime
	integer np,lmax
	integer isub
	integer k,l
	character*80 file,filename
	real u,v,s,d
	real flag

	logical, save :: bmeteo
	character*80, save :: string
	type(fem_type), save :: ffemuv
	integer, save :: icall = 0

	if( icall < 0 ) return

	if( icall == 0 ) then
	  if( ivar == 0 ) icall = -1
	  if( icall < 0 ) return
	  icall = 1
	  iformat = ffem%femfile%iformat
	  call strings_meteo_convention(ivar,bmeteo)
          call ivar2filename(ivar,filename)
	  call ivar2string(ivar,string,isub)
	  file = 'out.' // trim(filename) // '.fem'
	  !write(6,*) 'opening file: ',trim(file)
	  call femutil_open_for_write(file,iformat,ffemuv)
          call handle_open_output_file(999,file)
	end if

	if( ffem%femrec%bchanged ) then
	  ffemuv%femrec = ffem%femrec
	  ffemuv%femrec%nvar = 2
	  ffemuv%femrec%strings(1) = trim(string)//' - x'
	  ffemuv%femrec%strings(2) = trim(string)//' - y'
	end if

	call femutil_get_time(ffem,atime)
	call femutil_set_time(ffemuv,atime)

	np = ffem%femrec%np
	lmax = ffem%femrec%lmax
	flag = ffem%femrec%flag

	do k=1,np
	  do l=1,lmax
	    u = ffem%femrec%data(l,k,ivx)
	    v = ffem%femrec%data(l,k,ivy)
	    if( u == flag ) then
	      s = flag
	      d = flag
	    else
	      call convert_uv_sd(u,v,s,d,bmeteo)
	    end if
	    ffem%femrec%data(l,k,ivx) = s
	    ffem%femrec%data(l,k,ivy) = d
	    ffemuv%femrec%data(l,k,1) = u
	    ffemuv%femrec%data(l,k,2) = v
	  end do
	end do

	call femutil_write_record_fem(ffemuv)

	return
	end

c*****************************************************************

	subroutine custom_elab(frec)

	use fem_util

	implicit none

	type(femrec_type) :: frec

	integer nvar,iv
	real fact,offset
	real flag
	character*80 descr,string

	flag = frec%flag
	nvar = frec%nvar

	descr = 'water level [m]'
	descr = 'temperature [C]'
	descr = 'salinity [psu]'

	do iv=1,nvar
          string = frec%strings(iv)
	  if( string /= descr ) then
	    write(6,*) 'changing description... : ',trim(string)
	    frec%strings(iv) = descr
	  end if
	end do

	fact = 1.
	offset = 0.20

	write(6,*) 'custom_elab: ',offset,fact

	where( frec%data /= flag ) 
	  frec%data = fact * frec%data + offset
	end where

	end

c*****************************************************************

	subroutine reg_expand_shell(frec,regexpand)

c shell to call expansion routine

	use fem_util

	implicit none

	type(femrec_type) :: frec
	integer regexpand

	integer np,lmax,nvar
	integer iv
	integer nx,ny
	real flag

	call reg_set_flag(frec,iv)

	np = frec%np
	lmax = frec%lmax
	nvar = frec%nvar

        nx = nint(frec%regpar(1))
        ny = nint(frec%regpar(2))
        flag = frec%regpar(7)

	if( nx*ny /= np ) then
	  write(6,*) 'np,nx,ny,nx*ny: ',np,nx,ny,nx*ny
	  stop 'error stop reg_expand_shell: incompatible params'
	end if

	write(6,*) 'expanding regular grid: ',nx,ny,regexpand

	do iv=1,nvar
	  call reg_expand_3d(lmax,nx,ny,lmax,regexpand,flag
     +				,frec%data(1,1,iv))
	  call adjust_reg_vertical(lmax,nx,ny,flag
     +				,frec%data(1,1,iv),frec%ilhkv)
	end do

	end

c*****************************************************************

	subroutine reg_set_flag(frec,iv)

	use fem_util

	implicit none

	type(femrec_type) :: frec
	integer iv

	integer k,lmax,lm,np
	real flag

	np = frec%np
	lmax = frec%lmax
	flag = frec%regpar(7)

	do k=1,np
	  lm = frec%ilhkv(k)
	  frec%data(lm+1:lmax,k,iv) = flag
	end do

	end

c*****************************************************************

	subroutine copy_one_point(ip,lmax,nvar,frec,data3d)

	use fem_util

	implicit none

	integer ip,lmax,nvar
	type(femrec_type) :: frec
	real data3d(lmax,nvar)

	integer iv,ilmax
	real depth

	if( ip == 0 ) then			!condense - all nodes
	  ilmax = maxval(frec%ilhkv)
	  depth = maxval(frec%hd)
	else
	  ilmax = frec%ilhkv(ip)
	  depth = frec%hd(ip)
	end if

	call femutil_alloc_record_np(frec,1)
        call fem_file_set_ntype(frec%ntype,2,0)
	frec%ilhkv(1) = ilmax
	frec%hd(1) = depth

	do iv=1,nvar
	  frec%data(:,1,iv) = data3d(:,iv)
	end do

	end

c*****************************************************************

	subroutine fem_condense(lmax,nvar,frec,data3d)

	use fem_util

	implicit none

	integer lmax,nvar
	type(femrec_type) :: frec
	real data3d(lmax,nvar)

	integer l,i,iv,np
	integer nacu
	real flag
	double precision val,acu

	np = frec%np
	flag = frec%flag

	do iv=1,nvar
	 do l=1,lmax
	  nacu = 0
	  acu = 0.
	  do i=1,np
	    val = frec%data(l,i,iv)
	    if( val /= flag ) then
	      nacu = nacu + 1
	      acu = acu + val
	    end if
	  end do
	  if( nacu == 0 ) then
	    val = flag
	  else
	    val = acu / nacu
	  end if
	  data3d(l,iv) = val
	 end do
	end do

	end

c*****************************************************************

	subroutine handle_extract(breg,bquiet,np,regpar,iextract)

	use clo

	implicit none

	logical breg
	logical bquiet
	integer np
	real regpar(7)
	integer iextract

	logical berror
	integer inode,icoord,ie
	integer nx,ny,ix,iy,ixx,iyy
	real x0,y0,dx,dy,xp,yp,x,y
	real cx,cy
	real dist,xydist
	real f(3)
	character*80 snode,scoord

	integer iscanf

	x0 = 0.
	y0 = 0.
	dx = 0.
	dy = 0.

        call clo_get_option('nodei',snode)
        call clo_get_option('coord',scoord)

	inode = iscanf(snode,f,3)
	icoord = iscanf(scoord,f,3)

	iextract = 0

	if( inode == 0 .and. icoord == 0 ) return

	if( inode /= 0 .and. icoord /= 0 ) then
	  write(6,*) 'only one of -node and -coord can be given'
	  stop 'error stop handle_extract: cannot extract'
	end if
	if( inode > 0 .and. inode > 2 ) then
	  write(6,*) 'parse error for -node: need one or two numbers'
	  write(6,*) '-node:  ',trim(snode)
	  stop 'error stop handle_extract: need 1 or 2 values'
	end if
	if( icoord > 0 .and. icoord /= 2 ) then
	  write(6,*) 'parse error for -coord: need two numbers'
	  write(6,*) '-coord:  ',trim(scoord)
	  stop 'error stop handle_extract: need 2 values'
	end if

	if( inode == 1 ) then		!continuous numbering
	  iextract = nint(f(1))
	  if( iextract > np .or. iextract < 1 ) then
	    write(6,*) 'cannot extract this node: ',iextract
	    write(6,*) 'min,max possible: ',1,np
	    stop 'error stop handle_extract: iextract out of range'
	  end if
	end if

	if( ( inode == 2 .or. icoord == 2 ) .and. .not. breg ) then
	  write(6,*) 'need regular file to extract with these values:'
	  write(6,*) '-node:  ',trim(snode)
	  write(6,*) '-coord: ',trim(scoord)
	  stop 'error stop handle_extract: need regular file'
	end if

	if( breg ) then
	  nx = nint(regpar(1))
	  ny = nint(regpar(2))
	  x0 = regpar(3)
	  y0 = regpar(4)
	  dx = regpar(5)
	  dy = regpar(6)
	end if

	if( inode == 2 ) then
	  ix = nint(f(1))
	  iy = nint(f(2))
	  berror = .false.
	  if( ix < 1 .or. ix > nx ) berror = .true.
	  if( iy < 1 .or. iy > ny ) berror = .true.
	  if( berror ) then
	    write(6,*) 'ix,iy out of range'
	    write(6,*) 'ix,iy,nx,ny: ',ix,iy,nx,ny
	    stop 'error stop handle_extract: out of range'
	  end if
	  iextract = (iy-1)*nx + ix
	end if

	if( icoord == 2 ) then
	  xp = f(1)
	  yp = f(2)
	  cx = xp
	  cy = yp
	  ixx = 1
	  iyy = 1
	  xydist = (x0-xp)**2 + (y0-yp)**2
	  do iy=1,ny
	    y = y0 + (iy-1)*dy
	    do ix=1,nx
	      x = x0 + (ix-1)*dx
	      dist = (x-xp)**2 + (y-yp)**2
	      if( dist < xydist ) then
		xydist = dist
	        ixx = ix
	        iyy = iy
	      end if
	    end do
	  end do
	  iextract = (iyy-1)*nx + ixx
	end if

	if( breg ) then
	  ie = iextract
	  iy =  1 + (ie-1)/nx
	  ix = ie - (iy-1)*nx
	  xp = x0 + (ix - 1) * dx
	  yp = y0 + (iy - 1) * dy
	  if( .not. bquiet ) then
	    write(6,*)   'regular grid:          ',nx,ny
	    if( icoord > 0 ) then
	      write(6,*) 'requested coords:      ',cx,cy
	    end if
	    write(6,*)   'extracting point:      ',ix,iy
	    write(6,*)   'extracting coords:     ',xp,yp
	  end if
	end if

	if( .not. bquiet ) then
	  write(6,*) 'extracting point (1d): ',iextract
	end if

	end

c*****************************************************************

	subroutine make_grd_from_data(iv,nrec,nlvdi,np
     +			,regpar,data)

	implicit none

	integer iv,nrec
	integer nlvdi,np
	real regpar(7)
	real data(nlvdi,np)
	
	integer nx,ny
	real x0,y0,dx,dy,flag
	real x,y,val
	integer ix,iy,ip,it
	character*80 name

	nx = nint(regpar(1))
	ny = nint(regpar(2))
	x0 = regpar(3)
	y0 = regpar(4)
	dx = regpar(5)
	dy = regpar(6)
	flag = regpar(7)

	call make_name(iv,nrec,name)
	open(101,file=name,status='unknown',form='formatted')

	ip = 0
	do iy=1,ny
	  y = y0 + (iy-1)*dy
	  do ix=1,nx
	    x = x0 + (ix-1)*dx
	    ip = ip + 1
	    val = data(1,ip)
	    it = 1
	    if( val == flag ) it = 3
	    write(101,1000) 1,ip,it,x,y,val
	  end do
	end do

 1000	format(i1,i10,i5,3f14.6)
	close(101)

	end

c*****************************************************************

	subroutine make_name(iv,nrec,name)

	implicit none

	integer iv,nrec
	character*(*) name

	integer i

	write(name,'(a,i5,i7,a)') 'data',iv,nrec,'.grd'

	do i=1,len(trim(name))
	  if( name(i:i) == ' ' ) name(i:i) = '_'
	end do

	end

c*****************************************************************

	subroutine fem_get_new_atime(iformat,iunit,atime,ierr)

	implicit none

	integer iformat,iunit
	double precision atime
	integer ierr

	double precision dtime
	integer nvers,np,lmax,nvar,ntype
	integer datetime(2)

	atime = -1

	call fem_file_peek_params(iformat,iunit,dtime
     +                          ,nvers,np,lmax,nvar,ntype,datetime,ierr)

	if( ierr /= 0 ) return

	call dts_convert_to_atime(datetime,dtime,atime)

	end

c*****************************************************************

	subroutine set_facts(nvar,facts,factstring)

	implicit none

	integer nvar
	real facts(nvar)
	character*(*) factstring

        integer ianz
        integer iscanf

	facts = 1.

	if( factstring == ' ' ) return

        ianz = iscanf(factstring,facts,4)

        if( ianz /= nvar ) then
          write(6,*) 'looking for ',nvar,' factors'
          write(6,*) 'factors given: ',ianz
          write(6,*) 'string given: ',trim(factstring)
          stop 'error stop set_facts: nvar'
        end if


	end

c*****************************************************************

	subroutine change_strings(nvar,strings_out,newstring)

	implicit none

	integer nvar
	character*(*) strings_out(nvar)
	character*(*) newstring

	integer ia,ic,ics,i
	character*80 string

	if( newstring == ' ' ) return

	write(6,*) 'using string: ',trim(newstring)
	ia = 1
	ics = 0
	do
	  ic = scan(newstring(ia:),',')
	  if( ic == 0 ) then
	    string = newstring(ia:)
	  else
	    ic = ic + ia - 1
	    string = newstring(ia:ic-1)
	  end if
	  ics = ics + 1
	  !write(6,*) ia,ic,ics,trim(string)
	  if( string /= ' ' ) strings_out(ics) = string
	  if( ic == 0 .or. ics == nvar ) exit
	  ia = ic + 1
	end do

	return

	write(6,*) 'new strings description: '
	do i=1,nvar
	  write(6,*) i,'  ',trim(strings_out(i))
	end do

	end

c*****************************************************************

	subroutine correct_regpar(regpar)

	implicit none

	real regpar(7)

	if( regpar(3) > 180. ) then
	  regpar(3) = regpar(3) - 360.
	end if

	end

c*****************************************************************

	subroutine rain_elab(ffem,flag)

	use fem_util

	implicit none

	type(fem_type) :: ffem
	real flag

	integer ys(8)
	integer nr,i,np
	integer, save :: year = 0
	integer, save :: day = 0
	integer, save :: ny = 0
	double precision rtot,tperiod,tyear,totyear,atime
	double precision, save :: tot = 0.
	double precision, save :: atstart = 0.
	double precision, save :: atold = 0.
	integer, allocatable, save :: nralloc(:)
	real, allocatable, save :: rain(:)
	double precision, allocatable, save :: ralloc(:)
	type(fem_type), save :: ffem_out

	logical, save :: brewrite = .false.
	integer, save :: icall = 0
	integer, save :: npalloc = 0
	integer lmax,nvar,ntype

        if( brewrite .and. icall == 0 ) then
	  call femutil_open_for_write_fem('outrain.fem'
     +				,ffem%femfile%iformat,ffem_out)
	  call femutil_copy(ffem,ffem_out)
	  np = ffem%femrec%np
	  nvar = ffem%femrec%nvar
	  lmax = ffem%femrec%lmax
	  if( nvar > 1 .or. lmax > 1 ) then
	    write(6,*) 'nvar,lmax: ',nvar,lmax
	    stop 'error stop rain_elab: nvar /= 1 /= lmax'
	  end if
	  allocate(nralloc(np),ralloc(np),rain(np))
	  npalloc = np
	  ralloc = 0.
	  nralloc = 0
	  call femutil_get_time(ffem,atime)
	  atold = atime
	  atstart = atime
        end if

	call femutil_get_time(ffem,atime)
	call dts_from_abs_time_to_ys(atime,ys)

	rain = ffem%femrec%data(1,:,1)

	nr = 0
	rtot = 0
	np = ffem%femrec%np
	do i=1,np
	  if( rain(i) == flag ) cycle
	  nr = nr + 1
	  rtot = rtot + rain(i)
	end do
	rtot = rtot / nr

	!write(6,*) year,ny,nr,rtot
        !write(6,*) i,atime,ys(1)

        ny = ny + 1
	rtot = rtot / 86400.			!convert to mm/s
	tot = tot + rtot * (atime - atold)	!integrate
	atold = atime

        if( year /= ys(1) ) then
	  if( icall > 0 ) then
	    tperiod = (atime - atstart)
	    tyear = 365 * 86400.
	    if( 4*(year/4) == year ) tyear = tyear + 86400.
	    totyear = tot * tyear / tperiod	!for total year
            write(6,*) year,ny,tot,totyear
	  else
            write(6,*) '       year       nrecs   accumulated' //
     +			'             yearly'
	  end if
          ny = 0
	  atstart = atime
          tot = 0.
          year = ys(1)
        end if

	icall = 1

	if( .not. brewrite ) return

! from here on rewriting rain file because of errors in original

	if( np /= npalloc ) stop 'error stop: np/=npalloc'

	where ( rain /= flag )
	  ralloc = ralloc + rain
	  nralloc = nralloc + 1
	end where

	if( day /= ys(3) ) then			!average and write daily
	  where ( nralloc > 0 )
	    !ralloc = ralloc / nralloc		!no averaging
	  else where
	    ralloc = flag
	  end where
	  ffem_out%femrec%data(1,:,1) = ralloc(:)
	  call femutil_set_time(ffem_out,atime)
	  call femutil_write_record(ffem_out)
	  day = ys(3)
	  ralloc = 0.
	  nralloc = 0
	end if

	end

c*****************************************************************

        subroutine handle_open_output_file(ius,file)

! remebers opened files to write them to terminal at end of run

        implicit none

        integer ius			!0: write to terminal >0: insert
        character*(*) file

        integer, parameter :: nmax = 100
        integer, save :: nfill = 0
        character*80, save :: openfiles(nmax)

        integer i

        if( ius == 0 ) then
          do i=1,nfill
            write(6,*) '  ',trim(openfiles(i))
          end do
        else if( ius > 0 ) then
          nfill = nfill + 1
          if( nfill > nmax ) then
            stop 'error stop handle_open_output_file: nfill>nmax'
          end if
          openfiles(nfill) = file
        end if

        end

c*****************************************************************

	subroutine handle_changed(ffem,ffem_old,bterm)

	use fem_util

	type(fem_type) :: ffem
	type(fem_type) :: ffem_old
	logical bterm

	double precision atime
	character*20 aline
	integer, save :: np0 = 0
	integer, save :: lmax0 = 0
	integer, save :: icall = 0

	if( icall == 0 ) then
	  icall = 1
	  np0   = ffem%femrec%np
	  lmax0 = ffem%femrec%lmax
	end if

	np = ffem%femrec%np
	lmax = ffem%femrec%lmax

	call femutil_get_time(ffem,atime)

	if( lmax .ne. lmax0 .or. np .ne. np0 ) then
	    if( bterm ) then
	      call dts_format_abs_time(atime,aline)
	      write(6,*) 
	      write(6,*) '*** warning: parameters have changed'
	      write(6,*) 'time: ',trim(aline)
	      write(6,'(a)') ' parameter     old       new'
	      write(6,1300) ' np     ',np0,np
	      write(6,1300) ' lmax   ',lmax0,lmax
 1300	      format(a,2i10)
	      if( lmax > 1 ) then
	        write(6,*) 'levels: ',lmax
	        write(6,'(5f12.4)') hlv
	      end if
	    end if
	end if

	lmax0 = lmax
	np0 = np

        end

c*****************************************************************

	subroutine extract_point(ip,lmax,nvar,frec,array)

! ip > 0 => copy to array

	use fem_util

	implicit none

	integer ip,lmax,nvar
	type(femrec_type) :: frec
	real array(lmax,nvar)

	integer iv

	do iv=1,nvar
	  array(:,iv) = frec%data(:,ip,iv)
	end do

	end

c*****************************************************************

	subroutine fem_check_shell(atime,frec,scheck,bquiet)

	use fem_util

	implicit none

	double precision atime
	type(femrec_type) :: frec
        character*(*) scheck            !period on which to average
        logical bquiet

	integer np,lmax,nvar
	real flag

	np = frec%np
	lmax = frec%lmax
	nvar = frec%nvar
	flag = frec%flag

	call fem_check(atime,np,lmax,nvar,frec%data,flag
     +				,frec%strings,scheck,bquiet)

	end

c*****************************************************************

