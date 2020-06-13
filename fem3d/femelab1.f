
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
! 30.04.2020	ggu	bugfix for breg and iextract
! 13.06.2020	ggu	handle case with no data (NODATA)
!
!******************************************************************

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine femelab

c writes info on fem file

	use clo
	use elabutil
	use elabtime

	implicit none

	character*80 name,string
	integer np,iunit,iout
	integer nvers,lmax,nvar,ntype,nlvdi
	integer nvar0,lmax0,np0
	integer idt,idtact
	double precision dtime,atime0
	double precision atime,atold,atfirst,atlast,atnew
	real dmin,dmax,dmed
	integer ierr
	integer nfile
	integer nrec,iv,nrecs,l,i,ivar,lextr,nout
	integer itype(2)
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
	logical bdtok,bextract,bexpand,bcondense_txt
	logical bread,breg
	integer, allocatable :: ivars(:)
	character*80, allocatable :: strings(:)
	character*80, allocatable :: strings_out(:)
	character*20 aline,fline
	character*80 line
	!character*80 stmin,stmax
	real,allocatable :: facts(:)
	real,allocatable :: data(:,:,:)
	real,allocatable :: data_profile(:)
	real,allocatable :: dext(:)
	real,allocatable :: d3dext(:,:)
	real,allocatable :: hd(:)
	real,allocatable :: hlv(:)
	integer,allocatable :: ilhkv(:)
	integer,allocatable :: ius(:)
	integer,allocatable :: ius_sd(:)
	integer,allocatable :: llmax(:)

	integer ifileo

!--------------------------------------------------------------------
	INTERFACE
	subroutine allocate_vars(nvar,np,lmax,hlv,hd,ilhkv
     +			,data_profile,d3dext,data)
	integer nvar,np,lmax
	real, allocatable :: hlv(:)
	real, allocatable :: hd(:)
	integer, allocatable :: ilhkv(:)
	real, allocatable :: data_profile(:)
	real, allocatable :: d3dext(:,:)
	real, allocatable :: data(:,:,:)
	end
	END INTERFACE
!--------------------------------------------------------------------

	bhuman = .true.		!convert time in written fem file to dtime=0
	blayer = .false.
	blayer = .true.		!write layer structure - should be given by CLO
        bnewstring = .false.

	iextract = 0

        datetime = 0
	dtime = 0.
        nrec = 0
	regpar = 0.

c--------------------------------------------------------------
c initialization
c--------------------------------------------------------------

	call elabutil_init('FEM','femelab')

	if( bcondense ) bout = .true.
        bcondense_txt = .false.
	bexpand = regexpand > -1
	bskip = .not. bwrite
	if( bout ) bskip = .false.
	bextract = ( snode /= ' ' .or. scoord /= ' ' )

c--------------------------------------------------------------
c open file
c--------------------------------------------------------------

        call clo_reset_files
        call clo_get_next_file(infile)
	if( infile .eq. ' ' ) stop

	np = 0
	call fem_file_read_open(infile,np,iformat,iunit)
	if( iunit .le. 0 ) stop

	if( .not. bquiet ) then
	  write(6,*) 'file name: ',infile(1:len_trim(infile))
	  call fem_file_get_format_description(iformat,fline)
	  write(6,*) 'format: ',iformat,"  (",trim(fline),")"
	end if

c--------------------------------------------------------------
c prepare for output if needed
c--------------------------------------------------------------

        iout = 0
	iformout = iformat
	if( bchform ) iformout = 1 - iformat
	if( iformout < 0 ) iformout = iformat

        boutput = bout
	boutput = boutput .or. bchform
	boutput = boutput .or. newstring /= ' '
	boutput = boutput .or. bexpand
	if( bextract ) boutput = .false.

        if( boutput ) then
          iout = iunit + 1
          if( iformout .eq. 1 ) then
	    open(iout,file='out.fem',status='unknown',form='formatted')
          else
	    open(iout,file='out.fem',status='unknown',form='unformatted')
          end if
        end if

        date = 0
        time = 0
        call elabtime_date_and_time(date,time)  !we work with absolute time
	call elabtime_set_minmax(stmin,stmax)

	call elabtime_set_inclusive(binclusive)

c--------------------------------------------------------------
c read first record
c--------------------------------------------------------------

        call fem_file_read_params(iformat,iunit,dtime
     +                          ,nvers,np,lmax,nvar,ntype,datetime,ierr)

	if( ierr .ne. 0 ) goto 99

	if( .not. bquiet ) then
	  write(6,*) 'nvers:  ',nvers
	  write(6,*) 'np:     ',np
	  write(6,*) 'lmax:   ',lmax
	  write(6,*) 'nvar:   ',nvar
	  write(6,*) 'ntype:  ',ntype
	end if

	nvar0 = nvar
	lmax0 = lmax
	nlvdi = lmax
	np0 = np

	call allocate_vars(nvar,np,lmax,hlv,hd,ilhkv
     +			,data_profile,d3dext,data)

	allocate(ius(nvar))
	allocate(ius_sd(nvar))
	allocate(ivars(nvar))
	allocate(strings(nvar))
	allocate(strings_out(nvar))
	allocate(facts(nvar))
	allocate(dext(nvar))
	allocate(llmax(nvar))

	call fem_file_make_type(ntype,2,itype)

	call fem_file_read_2header(iformat,iunit,ntype,lmax
     +			,hlv,regpar,ierr)
	if( ierr .ne. 0 ) goto 98
	call correct_regpar(regpar)

	breg = itype(2) .gt. 0

	if( .not. bquiet ) then
	 if( bverb .and. lmax > 1 ) then
	  write(6,*) 'levels: ',lmax
	  write(6,'(5f12.4)') hlv
	 end if
	 if( itype(1) .gt. 0 ) then
	  write(6,*) 'date and time: ',datetime
	 end if
	 if( breg ) then
	  write(6,*) 'regpar: '
	  call printreg(regpar)
	 end if
	end if

	call dts_convert_to_atime(datetime,dtime,atime)

	atime0 = atime		!absolute time of first record
	atfirst = atime
	atlast = atime
	atold = atime

	if( bextract ) then
	  bskip = .false.
	  call handle_extract(breg,bquiet,np,regpar,iextract)
	end if

	ius = 0
	ius_sd = 0

c--------------------------------------------------------------
c write info to terminal
c--------------------------------------------------------------

	if( .not. bquiet ) then
          write(6,*) 'available variables contained in file: '
          write(6,*) 'total number of variables: ',nvar
          write(6,*) '   varnum     varid    varname'
	end if

	do iv=1,nvar
	  call fem_file_skip_data(iformat,iunit
     +                          ,nvers,np,lmax,string,ierr)
	  if( ierr .ne. 0 ) goto 97
	  !string = adjustl(string)
	  call string2ivar(string,ivar)
	  if( .not. bquiet ) then
	    write(6,'(2i10,4x,a)') iv,ivar,trim(string)
	  end if
	  ivars(iv) = ivar
	  strings(iv) = string
	end do

	strings_out = strings
        if( newstring /= ' ' ) then
          bnewstring = .true.
	  call change_strings(nvar,strings_out,newstring)
        end if
	call set_facts(nvar,facts,factstring)

	if( binfo ) return

c--------------------------------------------------------------
c close and re-open file
c--------------------------------------------------------------

	close(iunit)

	np = 0
	call fem_file_read_open(infile,np,iformat,iunit)
	if( iunit .le. 0 ) stop

c--------------------------------------------------------------
c loop on all records
c--------------------------------------------------------------

	bread = bwrite .or. bextract .or. boutput
	bread = bread .or. bsplit .or. bgrd .or. bcheck
	bread = bread .or. bcheckrain
	bskip = .not. bread

	atime = atold
	nrec = 0
	idt = 0

	do 
	  atold = atime
          call fem_file_read_params(iformat,iunit,dtime
     +                          ,nvers,np,lmax,nvar,ntype,datetime,ierr)
	  if( ierr .lt. 0 ) exit
	  if( ierr .gt. 0 ) goto 99
	  if( nvar .ne. nvar0 ) goto 96

	  nrec = nrec + 1

	  call dts_convert_to_atime(datetime,dtime,atime)
	  call dts_format_abs_time(atime,aline)
	  atlast = atime
	  atnew = atime

	  if( lmax .ne. lmax0 .or. np .ne. np0 ) then
	    if( .not. bquiet ) then
	      write(6,*) 
	      write(6,*) '*** warning: parameters have changed'
	      write(6,*) 'time: ',trim(aline)
	      write(6,'(a)') ' parameter     old       new'
	      write(6,1300) ' np     ',np0,np
	      write(6,1300) ' lmax   ',lmax0,lmax
	      write(6,1300) ' nvar   ',nvar0,nvar
 1300	      format(a,2i10)
	      if( bverb .and. lmax > 1 ) then
	       write(6,*) 'levels: ',lmax
	       write(6,'(5f12.4)') hlv
	      end if
	    end if
	    lmax0 = lmax
	    nlvdi = lmax
	    np0 = np
	    call allocate_vars(nvar,np,lmax,hlv,hd,ilhkv
     +			,data_profile,d3dext,data)
	  end if

	  call fem_file_read_2header(iformat,iunit,ntype,lmax
     +			,hlv,regpar,ierr)
	  if( ierr .ne. 0 ) goto 98
	  call correct_regpar(regpar)

	  call fem_file_make_type(ntype,2,itype)
	  breg = ( itype(2) .gt. 0 )
	  flag = -999.
	  if( breg ) flag = regpar(7)

	  do iv=1,nvar
	    if( bskip ) then
	      call fem_file_skip_data(iformat,iunit
     +                          ,nvers,np,lmax,string,ierr)
	    else
              call fem_file_read_data(iformat,iunit
     +                          ,nvers,np,llmax(iv)
     +                          ,string
     +                          ,ilhkv,hd
     +                          ,nlvdi,data(1,1,iv)
     +                          ,ierr)
	    end if
	    if( ierr .ne. 0 ) goto 97
	    !call custom_elab(nlvdi,np,string,iv,flag,data(1,1,iv))
	    if( bcheckrain ) then
	      call rain_elab(np,atime,string,regpar,flag,data(1,1,iv))
	    end if
	    if( .not. bnewstring .and. string .ne. strings(iv) ) goto 95
	    ffact = facts(iv)
	    if( ffact /= 1. ) then
	      where( data(:,:,iv) /= flag )
	        data(:,:,iv) = data(:,:,iv) * ffact
	      end where
	    end if
	  end do

	  call fem_get_new_atime(iformat,iunit,atnew,ierr)	!peeking
	  if( ierr < 0 ) atnew = atime
	  if( ierr > 0 ) goto 94

          if( elabtime_over_time(atime,atnew,atold) ) exit
          if( .not. elabtime_in_time(atime,atnew,atold) ) cycle

	  if( bverb ) then
            write(6,'(a,i8,f20.2,3x,a20)') 'time : ',nrec,atime,aline
	  end if

	  call handle_timestep(atime,bcheckdt,btskip)
	  if( btskip ) cycle

	  atlast = atime

          if( boutput ) then
	    if( bhuman ) then
	      call dts_convert_from_atime(datetime,dtime,atime)
	    end if
	    np_out = np
            ntype_out = ntype
	    if( bcondense ) then
              call fem_file_set_ntype(ntype_out,2,0)
              np_out = 1
              bcondense_txt = ( nvar == 1 .or. lmax == 1 )
            end if
            call fem_file_write_header(iformout,iout,dtime
     +                          ,0,np_out,lmax,nvar,ntype_out,lmax
     +                          ,hlv,datetime,regpar)
          end if

	  do iv=1,nvar

	    string = strings_out(iv)
	    !write(6,*) iv,'  ',trim(string)
            if( boutput ) then
	      !call custom_elab(nlvdi,np,string,iv,flag,data(1,1,iv))
	      if( breg .and. bexpand ) then
		call reg_set_flag(nlvdi,np,ilhkv,regpar,data(1,1,iv))
		call reg_expand_shell(nlvdi,np,llmax(iv),regexpand
     +					,regpar,ilhkv,data(1,1,iv))
	      end if
	      if( bcondense ) then
		call fem_condense(np,llmax(iv),data(1,1,iv),flag,
     +				data_profile)
                call fem_file_write_data(iformout,iout
     +                          ,0,np_out,llmax(iv)
     +                          ,string
     +                          ,ilhkv,hd
     +                          ,nlvdi,data_profile)
                d3dext(:,iv) = data_profile
	      else
                call fem_file_write_data(iformout,iout
     +                          ,0,np_out,llmax(iv)
     +                          ,string
     +                          ,ilhkv,hd
     +                          ,nlvdi,data(1,1,iv))
	      end if
            end if
	    if( bwrite ) then
	      write(6,*) nrec,iv,ivars(iv),trim(strings(iv))
	      do l=1,llmax(iv)
                call minmax_data(l,llmax(iv),np,flag,ilhkv,data(1,1,iv)
     +					,dmin,dmax,dmed)
	        write(6,1000) 'l,min,aver,max : ',l,dmin,dmed,dmax
 1000	        format(a,i5,3g16.6)
	      end do
	    end if
	    if( bextract ) then
	      dext(iv) = data(1,iextract,iv)
	      d3dext(:,iv) = data(:,iextract,iv)
	    end if
	    if( bsplit ) then
	      call femsplit(iformout,ius(iv),dtime,nvers,np
     +			,llmax(iv),nlvdi,ntype
     +			,hlv,datetime,regpar,string
     +			,ilhkv,hd,data(:,:,iv))
	      call femsplit_sd(iformout,ius_sd(iv),dtime,nvers,np
     +			,llmax(iv),nlvdi,ntype
     +			,hlv,datetime,regpar,string
     +			,ilhkv,hd,data(:,:,iv))
	    end if
	  end do

	  if( bextract ) then
	    lextr = ilhkv(iextract)
	    depth = hd(iextract)
	    call write_extract(atime,nvar,lmax,strings
     +			,lextr,hlv,depth,dext,d3dext)
	  end if

	  if( bcondense_txt ) then
            call write_ts(atime,nvar,lmax,strings,d3dext)
          end if

	  if( bcheck ) then
	    call fem_check(atime,np,lmax,nvar,data,flag
     +				,strings,scheck,bquiet)
	  end if

	  if( bgrd ) then
	    if( .not. breg ) then
	      stop 'error stop: for bgrd grid must be regular'
	    end if
	    do iv=1,nvar
	      write(6,*) 'writing grd file: ',iv,nrec
	      call make_grd_from_data(iv,nrec,nlvdi,np
     +			,regpar,data(:,:,iv))
	    end do
	  end if

	end do

	if( bcheck ) then	!write final data
	  atime = -1.
	  call fem_check(atime,np,lmax,nvar,data,flag
     +				,strings,scheck,bquiet)
	end if

c--------------------------------------------------------------
c finish loop - info on time records
c--------------------------------------------------------------

        if( .not. bsilent ) then
          write(6,*)
          call dts_format_abs_time(atfirst,aline)
          write(6,*) 'first time record: ',aline
          call dts_format_abs_time(atlast,aline)
          write(6,*) 'last time record:  ',aline

	  call handle_timestep_last(bcheckdt)

          write(6,*)
          write(6,*) nrec ,' time records read'
          !write(6,*) nelab,' time records elaborated'
          !write(6,*) nout ,' time records written to file'
          write(6,*)
        end if

	close(iunit)
	if( iout > 0 ) close(iout)

	if( bcheck ) then	!write final message
	  atime = -2.
	  call fem_check(atime,np,lmax,nvar,data,flag
     +				,strings,scheck,bsilent)
	end if

	if( boutput .and. .not. bquiet ) then
          line = 'output written to out.fem'
          if( bcondense_txt ) then
            line = trim(line) // ' and out.txt'
          end if
	  write(6,*) trim(line)
        else if( bsplit ) then
          write(6,*) 'the following files have been written:'
          call handle_open_output_file(0,' ')
	end if

	if( bextract .and. .not. bquiet ) then
	  write(6,*) 'iextract = ',iextract
	  write(6,*) 'data written to out.fem and out.txt'
	end if

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	return
   91	continue
	write(6,*) 'iectract,np: ',iextract,np
	stop 'error stop femelab: no such node'
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
	!write(6,*) 'lmax,lmax0: ',lmax,lmax0	!this might be relaxed
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

        subroutine minmax_data(level,nlvddi,np,flag,ilhkv,data
     +				,vmin,vmax,vmed)

        implicit none

	integer level		!level for which minmax to compute (0 for all)
        integer nlvddi,np
	real flag
        integer ilhkv(1)
        real data(nlvddi,1)
	real vmin,vmax,vmed

        integer k,l,lmin,lmax,lm,ntot
        real v,high
	double precision vtot

	lmin = max(1,level)
	lmax = level
	if( level == 0 ) lmax = nlvddi

	ntot = 0
	vtot = 0.
	high = 1.e+30
        vmin = high
        vmax = -high

        do k=1,np
          lm = min(ilhkv(k),lmax)
          do l=lmin,lm
            v = data(l,k)
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
	character*80 varline

	integer ifileo

!	------------------------------------------------------------
!	initialize output
!	------------------------------------------------------------

	if( iu == 0 ) then
	  iu = ifileo(88,'out.txt','form','new')
          n = max(nvar,lmax)
	  write(eformat,'(a,i3,a)') '(a20,',n,'g14.6)'
	  !write(6,*) 'using format: ',trim(eformat)
	  call make_varline(nvar,strings,varline)
	  if( varline /= ' ' ) write(iu,'(a)') trim(varline)
	end if

	call dts_format_abs_time(atime,aline)

        if( nvar == 1 ) then
          write(iu,eformat) aline,data(1:lmax,1)
        else if( lmax == 1 ) then
          write(iu,eformat) aline,data(1,1:nvar)
        else
          write(6,*) nvar,lmax
          stop 'error stop write_ts: nvar and lmax incompatible'
        end if
          
        end

c*****************************************************************

        subroutine write_extract(atime,nvar,nlvddi,strings
     +                  ,lextr_in,hlv,depth,dext,d3dext)

	use iso8601

	implicit none

	double precision atime
	integer nvar,nlvddi
	character*(*) strings(nvar)
	integer lextr_in
	real hlv(nlvddi)
	real depth
	real dext(nvar)
	real d3dext(nlvddi,nvar)

	integer, save :: iu2d = 0
	integer, save :: iu3d = 0

	integer it,ierr
	integer iformat,nvers,np,ntype,ivar,lmax,lextr
	integer ilhkv(1)
	double precision dtime
	real, save :: regpar(7) = 0.
	integer datetime(2)
	real flag
	real hd(1)
	character*20 aline
	character*80, save :: eformat
	character*80 varline

	integer ifileo

!	------------------------------------------------------------
!	initialize output
!	------------------------------------------------------------

	if( iu2d == 0 ) then
	  iu2d = ifileo(88,'out.txt','form','new')
	  write(eformat,'(a,i3,a)') '(a20,',nvar,'g14.6)'
	  !write(6,*) 'using format: ',trim(eformat)
	  call make_varline(nvar,strings,varline)
	  if( varline /= ' ' ) write(iu2d,'(a80)') varline
	end if
	if( iu3d == 0 ) then
	  iu3d = ifileo(89,'out.fem','form','new')
	end if
	if( iu2d < 0 .or. iu3d < 0 ) then
	  write(6,*) 'iu2d = ',iu2d
	  write(6,*) 'iu3d = ',iu3d
	  stop 'error stop write_extract: error opening file'
	end if

!	------------------------------------------------------------
!	handle no data (NODATA)
!	------------------------------------------------------------

	lextr = lextr_in
	if( lextr < 1 ) then
	  lextr = 1
	  flag = -999.
	  d3dext(1,:) = flag
	end if
	regpar(7) = flag

!	------------------------------------------------------------
!	write 2d output
!	------------------------------------------------------------

	dext(:) = d3dext(1,:)
	call dts_format_abs_time(atime,aline)
	write(iu2d,eformat) aline,dext

!	------------------------------------------------------------
!	write 3d output
!	------------------------------------------------------------

	nvers = 0
	iformat = 1
	dtime = 0.
	np = 1
	ntype = 1
	lmax = lextr
	ilhkv(1) = lextr
	hd(1) = depth

	call string2date(aline,datetime,ierr)
	if( ierr /= 0 ) stop 'error stop write_extract: converting time'

        call fem_file_write_header(iformat,iu3d,dtime
     +                          ,nvers,np,lextr
     +                          ,nvar,ntype
     +                          ,lmax,hlv,datetime,regpar)
	do ivar=1,nvar
          call fem_file_write_data(iformat,iu3d
     +                          ,nvers,np,lmax
     +                          ,strings(ivar)
     +                          ,ilhkv,hd
     +                          ,lmax,d3dext(1:lmax,ivar))
	end do

!	------------------------------------------------------------
!	end of routine
!	------------------------------------------------------------

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

	subroutine femsplit(iformout,ius,dtime,nvers,np
     +			,lmax,nlvddi,ntype
     +			,hlv,datetime,regpar,string
     +			,ilhkv,hd,data)

! splits fem file into single variables - directional vars are kept in one file

	implicit none

	integer iformout,ius
	double precision dtime
	integer nvers,np,lmax,nlvddi,ntype
	real hlv(lmax)
	integer datetime(2)
	real regpar(7)
	character(*) string
	integer ilhkv(np)
	real hd(np)
	real data(nlvddi,np)

	integer ivar,nvar
	character*80 file,extra,filename
	character*1 dir
	character*10 unit
	integer, save :: iusold

	call string_direction_and_unit(string,dir,unit)

	if( ius == 0 ) then
	  if( dir == 'y' ) then		!is second part of vector
	    ius = iusold
	  else
	    call string2ivar(string,ivar)
            call ivar2filename(ivar,filename)
            call strings_pop_direction(filename)
	    file = 'out.' // trim(filename) // '.fem'
	    call fem_file_write_open(file,iformout,ius)
	    if( ius <= 0 ) goto 99
            call handle_open_output_file(ius,file)
	    iusold = ius
	  end if
	end if

	nvar = 1
	if( dir /= ' ' ) nvar = 2

	if( dir /= 'y' ) then
          call fem_file_write_header(iformout,ius,dtime
     +                          ,nvers,np,lmax
     +                          ,nvar,ntype
     +                          ,nlvddi,hlv,datetime,regpar)
	end if

        call fem_file_write_data(iformout,ius
     +                          ,nvers,np,lmax
     +                          ,string
     +                          ,ilhkv,hd
     +                          ,nlvddi,data)

	return
   99	continue
	write(6,*) 'cannot open file ',trim(file)
	stop 'error stop femsplit: cannot open output file'
	end

c*****************************************************************

	subroutine femsplit_sd(iformout,ius,dtime,nvers,np
     +			,lmax,nlvddi,ntype
     +			,hlv,datetime,regpar,string
     +			,ilhkv,hd,data)

! splits directional variables in fem file into speed and direction

	implicit none

	integer iformout,ius
	double precision dtime
	integer nvers,np,lmax,nlvddi,ntype
	real hlv(lmax)
	integer datetime(2)
	real regpar(7)
	character(*) string
	integer ilhkv(np)
	real hd(np)
	real data(nlvddi,np)

	integer ivar,nvar,ivars,ivard,isub
	integer k,l
	character*80 file,extra,filename
	character*1 dir
	character*10 unit
	real u,v,s,d
	real, parameter :: flag = -999.
        real, allocatable, save :: datax(:,:),datas(:,:),datad(:,:)
	integer, save :: iuss,iusd,iusx
	logical, save :: bmeteo
        character*80, save :: strings,stringd

	if( ius < 0 ) return

	if( ius == 0 ) then
	  call string2ivar(string,ivar)
	  call get_direction_ivars(ivar,ivars,ivard)
	  if( ivars == 0 ) then
	    ius = -1
	    return
	  end if
	  call string_direction_and_unit(string,dir,unit)
	  if( dir == 'x' ) then
            call ivar2filename(ivars,filename)
	    call ivar2string(ivars,strings,isub)
	  else if( dir == 'y' ) then
            call ivar2filename(ivard,filename)
	    call ivar2string(ivard,stringd,isub)
	  else
            write(6,*) 'unknown direction: ',trim(string),'  ',dir
            stop 'error stop femsplit_sd: unknown direction'
	  end if
	  call strings_meteo_convention(ivar,bmeteo)
	  file = 'out.' // trim(filename) // '.fem'
	  call fem_file_write_open(file,iformout,ius)
	  if( ius <= 0 ) goto 99
          call handle_open_output_file(ius,file)
	  if( dir == 'x' ) iusx = ius
	end if

	nvar = 1
	if( .not. allocated(datax) ) then
	  allocate(datax(nlvddi,np))
	  allocate(datas(nlvddi,np))
	  allocate(datad(nlvddi,np))
	end if

	if( ius == iusx ) then
	  datax = data
	  iuss = ius
	  return
	else
	  iusd = ius
	end if

	do k=1,np
	  do l=1,lmax
	    u = datax(l,k)
	    v = data(l,k)
	    if( u == flag ) then
	      s = flag
	      d = flag
	    else
	      call convert_uv_sd(u,v,s,d,bmeteo)
	    end if
	    datas(l,k) = s
	    datad(l,k) = d
	  end do
	end do

	ius = iuss
	string = strings

        call fem_file_write_header(iformout,ius,dtime
     +                          ,nvers,np,lmax
     +                          ,nvar,ntype
     +                          ,nlvddi,hlv,datetime,regpar)

        call fem_file_write_data(iformout,ius
     +                          ,nvers,np,lmax
     +                          ,string
     +                          ,ilhkv,hd
     +                          ,nlvddi,datas)

	ius = iusd
	string = stringd

        call fem_file_write_header(iformout,ius,dtime
     +                          ,nvers,np,lmax
     +                          ,nvar,ntype
     +                          ,nlvddi,hlv,datetime,regpar)

        call fem_file_write_data(iformout,ius
     +                          ,nvers,np,lmax
     +                          ,string
     +                          ,ilhkv,hd
     +                          ,nlvddi,datad)

	return
   99	continue
	write(6,*) 'cannot open file ',trim(file)
	stop 'error stop femsplit_sd: cannot open output file'
	end

c*****************************************************************

	subroutine custom_elab(nlvdi,np,string,iv,flag,data)

	implicit none

	integer nlvdi,np,iv
	character*(*) string
	real flag
	real data(nlvdi,np)

	real fact,offset
	character*80 descr

	!return

	!if( string(1:13) /= 'wind velocity' ) return
	!if( iv < 1 .or. iv > 2 ) return
	!if( string(1:11) /= 'water level' ) return

	descr = 'water level [m]'
	descr = 'temperature [C]'
	descr = 'salinity [psu]'
	if( string /= descr ) then
	  write(6,*) 'changing description... : ',trim(string)
	  string = descr
	end if
	return

	fact = 1.
	offset = 0.20

	!write(6,*) iv,'  ',trim(string)
	!write(6,*) 'attention: wind speed changed by a factor of ',fact
	write(6,*) 'custom_elab: ',offset,fact

	where( data /= flag ) 
	  data = fact * data + offset
	end where

	end

c*****************************************************************

	subroutine reg_expand_shell(nlvddi,np,lmax,regexpand
     +					,regpar,il,data)

c shell to call expansion routine

	implicit none

	integer nlvddi,np,lmax,regexpand
	real regpar(7)
	integer il(np)
	real data(nlvddi,np)

	integer nx,ny
	real flag

        nx = nint(regpar(1))
        ny = nint(regpar(2))
        flag = regpar(7)

	if( nx*ny /= np ) then
	  write(6,*) 'np,nx,ny,nx*ny: ',np,nx,ny,nx*ny
	  stop 'error stop reg_expand_shell: incompatible params'
	end if

	write(6,*) 'expanding regular grid: ',nx,ny,regexpand

	call reg_expand_3d(nlvddi,nx,ny,lmax,regexpand,flag,data)

	call adjust_reg_vertical(nlvddi,nx,ny,flag,data,il)

	end

c*****************************************************************

	subroutine reg_set_flag(nlvddi,np,il,regpar,data)

	implicit none

	integer nlvddi,np
	integer il(np)
	real regpar(7)
	real data(nlvddi,np)

	integer k,lmax
	real flag

	flag = regpar(7)

	do k=1,np
	  lmax = il(k)
	  data(lmax+1:nlvddi,k) = flag
	end do

	end

c*****************************************************************

	subroutine fem_condense(np,lmax,data,flag,data_profile)

	implicit none

	integer np,lmax
	real data(lmax,np)
	real flag
	real data_profile(lmax)

	integer l,i,nacu
	double precision val,acu

	do l=1,lmax
	  nacu = 0
	  acu = 0.
	  do i=1,np
	    val = data(l,i)
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
	  data_profile(l) = val
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
	   write(6,*) 'regular grid:          ',nx,ny
	   if( icoord > 0 ) write(6,*) 'requested coords:      ',cx,cy
	   write(6,*) 'extracting point:      ',ix,iy
	   write(6,*) 'extracting coords:     ',xp,yp
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

	subroutine allocate_vars(nvar,np,lmax,hlv,hd,ilhkv
     +			,data_profile,d3dext,data)

	implicit none

	integer nvar,np,lmax
	real, allocatable :: hlv(:)
	real, allocatable :: hd(:)
	integer, allocatable :: ilhkv(:)
	real, allocatable :: data_profile(:)
	real, allocatable :: d3dext(:,:)
	real, allocatable :: data(:,:,:)

	if( allocated(hlv) ) then
	  deallocate(hlv,hd,ilhkv,data_profile,d3dext,data)
	end if

	allocate(hlv(lmax))
	allocate(hd(np))
	allocate(ilhkv(np))
	allocate(data_profile(lmax))
	allocate(d3dext(lmax,nvar))
	allocate(data(lmax,np,nvar))

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

	subroutine rain_elab(np,atime,string,regpar,flag,data)

	implicit none

	integer np
	double precision atime
	character*(*) string
	real regpar(7)
	real flag
	real data(1,np)

	real rain(np)
	integer ys(8)
	integer nr,i
	integer, save :: year = 0
	integer, save :: day = 0
	integer, save :: ny = 0
	double precision, save :: tot = 0.
	double precision, save :: tstart = 0.
	double precision, save :: aold = 0.
	double precision rtot,tperiod,tyear,totyear,dtime
	integer, allocatable, save :: nralloc(:)
	double precision, allocatable, save :: ralloc(:)

	logical, save :: brewrite = .false.
	integer, save :: iformat = 0
	integer, save :: iout = 999
	integer, save :: icall = 0
	integer, save :: npalloc = 0
	integer lmax,nvar,ntype
	integer datetime(2)
	real hlv(1)
	real hd(1)
	integer ilhkv(1)

        if( brewrite .and. icall == 0 ) then
          if( iformat .eq. 1 ) then
	    open(iout,file='out1.fem',status='unknown',form='formatted')
          else
	    open(iout,file='out1.fem',status='unknown',form='unformatted')
          end if
	  allocate(nralloc(np),ralloc(np))
	  npalloc = np
	  ralloc = 0.
	  nralloc = 0
        end if

	call dts_from_abs_time_to_ys(atime,ys)

	rain = data(1,:)

	nr = 0
	rtot = 0
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
	tot = tot + rtot * (atime - aold)	!integrate
	aold = atime

        if( year /= ys(1) ) then
	  if( icall > 0 ) then
	    tperiod = (atime - tstart)
	    tyear = 365 * 86400.
	    if( 4*(year/4) == year ) tyear = tyear + 86400.
	    totyear = tot * tyear / tperiod	!for total year
            write(6,*) year,ny,tot,totyear
	  else
            write(6,*) '       year       nrecs   accumulated' //
     +			'             yearly'
	  end if
          ny = 0
	  tstart = atime
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

	if( day /= ys(3) ) then			!average daily
	  where ( nralloc > 0 )
	    !ralloc = ralloc / nralloc		!no averaging
	  else where
	    ralloc = flag
	  end where
	  rain = ralloc
	  dtime = 0.
	  lmax = 1
	  nvar = 1
	  ntype = 1
	  if( regpar(1) > 0 ) ntype = 11
	  hlv(1) = 10000.
	  hd = 10000.
	  ilhkv = 1
	  call dts_from_abs_time(datetime(1),datetime(2),atime)
          call fem_file_write_header(iformat,iout,dtime
     +                          ,0,np,lmax,nvar,ntype,lmax
     +                          ,hlv,datetime,regpar)
          call fem_file_write_data(iformat,iout
     +                          ,0,np,lmax
     +                          ,string
     +                          ,ilhkv,hd
     +                          ,lmax,rain)
	  day = ys(3)
	  ralloc = 0.
	  nralloc = 0
	end if

	end

c*****************************************************************

        subroutine handle_open_output_file(ius,file)

        implicit none

        integer ius
        character*(*) file

        integer, parameter :: nmax = 10
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

