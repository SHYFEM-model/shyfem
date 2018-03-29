!
! elaborates fem files
!
! revision log :
!
! 14.01.2015    ggu     adapted from feminf
! 20.05.2015    ggu     use bhuman to convert to human readable time
! 05.06.2015    ggu     iextract to extract nodal value
! 05.11.2015    ggu     new option chform to change format
! 04.10.2016    ggu     output flags now similar to shyelab
! 05.10.2016    ggu     allow for expansion of regular grid
! 11.10.2016    ggu     introduced flag for min/max/med computation
! 31.10.2016    ggu     new flag condense (bcondense)
! 16.05.2017    ggu&mbj better handling of points to extract
! 31.08.2017    ggu     new flag -grd to write grd from fem file
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
	integer nrec,iv,ich,isk,nrecs,l,i,ivar
	integer itype(2)
	integer iformat,iformout
	integer date,time
	integer datetime(2),dateanf(2),dateend(2)
	integer iextract,it
	integer ie,nx,ny,ix,iy
	integer np_out
	real x0,y0,dx,dy,x1,y1
	real regpar(7)
	real xp,yp
	real ffact
	logical bfirst,bskip
	logical bhuman,blayer
	logical bdtok,bextract,bexpand
	logical bread,breg
	integer, allocatable :: ivars(:)
	character*80, allocatable :: strings(:)
	character*80, allocatable :: strings_out(:)
	character*20 dline,fline
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

	integer ifileo

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

	if( bcondense ) bout = .true.
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

	atime0 = atime		!absolute time of first record
	atfirst = atime
	atlast = atime
	atold = atime

	if( bextract ) then
	  bskip = .false.
	  call handle_extract(breg,np,regpar,iextract)
	end if

	nvar0 = nvar
	lmax0 = lmax
	nlvdi = lmax
	np0 = np
	allocate(ivars(nvar))
	allocate(strings(nvar))
	allocate(strings_out(nvar))
	allocate(facts(nvar))
	allocate(dext(nvar))
	allocate(d3dext(nlvdi,nvar))
	allocate(data(nlvdi,np,nvar))
	allocate(data_profile(nlvdi))
	allocate(hd(np))
	allocate(ilhkv(np))
	allocate(ius(nvar))
	ius = 0

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
	call change_strings(nvar,strings_out,newstring)
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
	bskip = .not. bread

	atime = atold
	nrec = 0
	idt = 0
	ich = 0
	isk = 0

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

	  do iv=1,nvar
	    if( bskip ) then
	      call fem_file_skip_data(iformat,iunit
     +                          ,nvers,np,lmax,string,ierr)
	    else
              call fem_file_read_data(iformat,iunit
     +                          ,nvers,np,lmax
     +                          ,string
     +                          ,ilhkv,hd
     +                          ,nlvdi,data(1,1,iv)
     +                          ,ierr)
	    end if
	    if( ierr .ne. 0 ) goto 97
	    if( string .ne. strings(iv) ) goto 95
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
            write(6,'(a,i8,f20.2,3x,a20)') 'time : ',nrec,atime,dline
	  end if

	  bdtok = atime == atfirst .or. atime > atold
	  call check_dt(atime,atold,bcheckdt,nrec,idt,ich,isk)
	  if( .not. bdtok ) cycle
	  atlast = atime

          if( boutput ) then
	    if( bhuman ) then
	      call dts_convert_from_atime(datetime,dtime,atime)
	    end if
	    np_out = np
	    if( bcondense ) np_out = 1
            call fem_file_write_header(iformout,iout,dtime
     +                          ,0,np_out,lmax,nvar,ntype,lmax
     +                          ,hlv,datetime,regpar)
          end if

	  do iv=1,nvar

	    string = strings_out(iv)
	    !write(6,*) iv,'  ',trim(string)
            if( boutput ) then
	      !call custom_elab(nlvdi,np,string,iv,data(1,1,iv))
	      if( breg .and. bexpand ) then
		call reg_set_flag(nlvdi,np,ilhkv,regpar,data(1,1,iv))
		call reg_expand_shell(nlvdi,np,lmax,regexpand
     +					,regpar,ilhkv,data(1,1,iv))
	      end if
	      if( bcondense ) then
		call fem_condense(np,lmax,data(1,1,iv),flag,data_profile)
                call fem_file_write_data(iformout,iout
     +                          ,0,np_out,lmax
     +                          ,string
     +                          ,ilhkv,hd
     +                          ,nlvdi,data_profile)
	      else
                call fem_file_write_data(iformout,iout
     +                          ,0,np_out,lmax
     +                          ,string
     +                          ,ilhkv,hd
     +                          ,nlvdi,data(1,1,iv))
	      end if
            end if
	    if( bwrite ) then
	      write(6,*) nrec,iv,ivars(iv),trim(strings(iv))
	      do l=1,lmax
                call minmax_data(l,lmax,np,flag,ilhkv,data(1,1,iv)
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
     +			,lmax,nlvdi,ntype
     +			,hlv,datetime,regpar,string
     +			,ilhkv,hd,data(:,:,iv))
	    end if
	  end do

	  if( bextract ) then
	    call write_extract(atime,nvar,lmax,strings,dext,d3dext)
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

	if( bcheckdt ) then
         if( ich > 0 ) then
          write(6,*) 'idt:     irregular ',ich,isk
         else if( .not. bquiet ) then
          write(6,*) 'idt:    ',idt
	 end if
        end if

	if( isk .gt. 0 ) then
	  write(6,*) '*** warning: records eliminated: ',isk
	end if

	close(iunit)
	if( iout > 0 ) close(iout)

	if( bcheck ) then	!write final message
	  atime = -2.
	  call fem_check(atime,np,lmax,nvar,data,flag
     +				,strings,scheck,bsilent)
	end if

	if( boutput .and. .not. bquiet ) then
	  write(6,*) 'output written to file out.fem'
	end if

	if( bextract .and. .not. bquiet ) then
	  write(6,*) 'iextract = ',iextract
	  write(6,*) 'data written to out.txt'
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
	write(6,*) 'strings not in same sequence: ',iv
        write(6,*) string
        write(6,*) strings(iv)
	stop 'error stop femelab: strings'
   96	continue
	write(6,*) 'nvar,nvar0: ',nvar,nvar0
	write(6,*) 'lmax,lmax0: ',lmax,lmax0	!this might be relaxed
	write(6,*) 'np,np0:     ',np,np0	!this might be relaxed
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

	subroutine check_dt(atime,atold,bcheckdt,nrec,idt,ich,isk)

	implicit none

	double precision atime,atold
	logical bcheckdt
	integer nrec,idt,ich,isk

	integer idtact
	character*20 dline

          if( nrec > 1 ) then
            if( nrec == 2 ) idt = nint(atime-atold)
            idtact = nint(atime-atold)
            if( idtact .ne. idt ) then
              ich = ich + 1
              if( bcheckdt ) then
		call dts_format_abs_time(atime,dline)
                write(6,'(a,3i10,a,a)') '* change in time step: '
     +                          ,nrec,idt,idtact,'  ',dline
              end if
              idt = idtact
            end if
            if( idt <= 0 ) then
	      isk = isk + 1
	      call dts_format_abs_time(atime,dline)
              write(6,*) '*** zero or negative time step: ',nrec,idt
     +                          ,atime,atold,'  ',dline
            end if
          end if

	end

c*****************************************************************

	subroutine write_extract(atime,nvar,lmax,strings,dext,d3dext)

	implicit none

	double precision atime
	integer nvar,lmax
	character*(*) strings(nvar)
	real dext(nvar)
	real d3dext(lmax,nvar)

	integer, save :: iu2d = 0
	integer, save :: iu3d = 0

	integer it
	double precision dtime
	character*20 dline
	character*80, save :: eformat
	character*80 varline

	integer ifileo

	if( iu2d == 0 ) then
	  iu2d = ifileo(88,'out.txt','form','new')
	  write(eformat,'(a,i3,a)') '(a20,',nvar,'g14.6)'
	  !write(6,*) 'using format: ',trim(eformat)
	  call make_varline(nvar,strings,varline)
	  if( varline /= ' ' ) write(iu2d,'(a80)') varline
	end if
	if( iu3d == 0 ) then
	  !iu3d = ifileo(89,'out.fem','form','new')
	end if

	dext(:) = d3dext(1,:)
	call dts_format_abs_time(atime,dline)
	write(iu2d,eformat) dline,dext

! need more data here for 3d
!	nvers
!        call fem_file_write_header(iformat,iu3d,dtime
!     +                          ,nvers,np,lmax
!     +                          ,nvar,ntype
!     +                          ,nlvddi,hlv,datetime,regpar)
!        call fem_file_write_data(iformat,iu3d
!     +                          ,nvers,np,lmax
!     +                          ,string
!     +                          ,ilhkv,hd
!     +                          ,nlvddi,temp1)

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

	varline = ' '
	iend = 0

	do iv=1,nvar
	  name = strings(iv)
	  call strings_get_short_name(name,short)
	  varline = varline(1:iend) // short
	  iend = iend + 14
	end do

	varline = '#vars: date_and_time   ' // varline

	end

c*****************************************************************

	subroutine femsplit(iformout,ius,dtime,nvers,np
     +			,lmax,nlvddi,ntype
     +			,hlv,datetime,regpar,string
     +			,ilhkv,hd,data)

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
	character*80 file,extra
	character*1 dir
	character*10 unit
	integer, save :: iusold

	if( ius == 0 ) then
	  call string2ivar(string,ivar)
	  call string_direction_and_unit(string,dir,unit)
	  if( dir == 'y' ) then		!is second part of vector
	    ius = iusold
	  else
	    call alpha(ivar,extra)
	    file = 'out.' // trim(extra) // '.fem'
	    call fem_file_write_open(file,iformout,ius)
	    if( ius <= 0 ) goto 99
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

	subroutine custom_elab(nlvdi,np,string,iv,data)

	implicit none

	integer nlvdi,np,iv
	character*(*) string
	real data(nlvdi,np)

	real fact

	return

	if( string(1:13) /= 'wind velocity' ) return
	if( iv < 1 .or. iv > 2 ) return

	fact = 2.

	!write(6,*) iv,'  ',trim(string)
	write(6,*) 'attention: wind speed changed by a factor of ',fact

	data = fact * data

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

	subroutine handle_extract(breg,np,regpar,iextract)

	use clo

	implicit none

	logical breg
	integer np
	real regpar(7)
	integer iextract

	logical berror
	integer inode,icoord,ie
	integer nx,ny,ix,iy,ixx,iyy
	real x0,y0,dx,dy,xp,yp,x,y
	real dist,xydist
	real f(3)
	character*80 snode,scoord

	integer iscanf

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
	  write(6,*) 'regular grid:          ',nx,ny
	  write(6,*) 'extracting point:      ',ix,iy
	  write(6,*) 'extracting coords:     ',xp,yp
	end if

	write(6,*) 'extracting point (1d): ',iextract

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

