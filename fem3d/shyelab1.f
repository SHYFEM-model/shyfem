
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

! revision log :
!
! 18.11.1998    ggu     check dimensions with dimnos
! 06.04.1999    ggu     some cosmetic changes
! 03.12.2001    ggu     some extra output -> place of min/max
! 09.12.2003    ggu     check for NaN introduced
! 07.03.2007    ggu     easier call
! 08.11.2008    ggu     do not compute min/max in non-existing layers
! 07.12.2010    ggu     write statistics on depth distribution (depth_stats)
! 06.05.2015    ggu     noselab started
! 05.06.2015    ggu     many more features added
! 10.09.2015    ggu     std and rms for averaging implemented
! 11.09.2015    ggu     write in gis format
! 23.09.2015    ggu     handle more than one file (look for itstart)
! 16.10.2015    ggu     started shyelab
! 10.06.2016    ggu     shydiff included
! 08.09.2016    ggu     custom dates, map_influence
! 11.05.2017    ggu     use catmode to concatenate files
! 05.10.2017    ggu     implement silent option
! 07.10.2017    ggu     new names for -split option of hydro file
! 11.05.2018    ggu     call shympi_init later (after basin)
!
!**************************************************************

	subroutine shyelab1

	use clo
	use elabutil
	use elabtime
	use shyfile
	use shyutil
	use custom_dates

        use basin
        use mod_depth
        use evgeom
        use levels
        use shympi

	implicit none

	real, allocatable :: cv2(:)
	real, allocatable :: cv3(:,:)
	real, allocatable :: cv3all(:,:,:)
	real, allocatable :: cv3diff(:,:,:)

	integer, allocatable :: idims(:,:)
	integer, allocatable :: ivars(:)
	character*80, allocatable :: strings(:)
	integer, allocatable :: il(:)

	real, allocatable :: znv(:)
	real, allocatable :: uprv(:,:)
	real, allocatable :: vprv(:,:)
	real, allocatable :: sv(:,:)
	real, allocatable :: dv(:,:)

	logical bhydro,bscalar
	logical blastrecord,bforce
	integer nwrite,nwtime,nread,nelab,nrec,nin,nold,ndiff
	integer nvers
	integer nvar,npr
	integer ierr
	integer date,time
	integer it
	integer ivar,iaux
	integer iv,j,l,k,lmax,node
	integer ip
	integer ifile,ftype
	integer id,idout,iddiff
	integer n,m,nndim,nn
	integer naccum
	character*80 title,name,file
	character*80 basnam,simnam
	character*20 aline
	real rnull
	real cmin,cmax,cmed,cstd,atot,vtot
	double precision dtime,dtstart,dtnew,ddtime
	double precision atfirst,atlast
	double precision atime,atstart,atnew,atold

 	!logical, parameter :: bmap = .false.
 	real, parameter :: pthresh = 30.
 	real, parameter :: cthresh = 20.
 	!real, parameter :: cthresh = 0.

	integer iapini,i
	integer ifem_open_file
	logical concat_cycle_a

!--------------------------------------------------------------
! initialize everything
!--------------------------------------------------------------

	nread=0
	nelab=0
	nrec=0
	ndiff = 0
	rnull=0.
	rnull=-1.
	bopen = .false.
	bzeta = .false.		!file has zeta information
	ifile = 0
	id = 0
	bforce = .false.

	!--------------------------------------------------------------
	! set command line parameters
	!--------------------------------------------------------------

	call elabutil_init('SHY','shyelab')

	!call shympi_init(.false.)

	!--------------------------------------------------------------
	! open input files
	!--------------------------------------------------------------

	call open_new_file(ifile,id,atstart)	!atstart=-1 if no new file
	if( bverb ) call shy_write_filename(id)
	if( atstart /= -1 ) then
	  call dts_format_abs_time(atstart,aline)
	  if( .not. bsilent ) then
	    write(6,*) 'initial date for next file: ',aline
	  end if
	end if

	if( bdiff ) then
	  if( .not. clo_exist_file(ifile+1) ) goto 66
	  if( clo_exist_file(ifile+2) ) goto 66
	  call open_next_file(ifile+1,id,iddiff)
	  atstart = -1
	end if

	!--------------------------------------------------------------
	! set up params and modules
	!--------------------------------------------------------------

	call shy_get_params(id,nkn,nel,npr,nlv,nvar)
	call shy_get_ftype(id,ftype)

	if( .not. bquiet ) call shy_info(id)

        call basin_init(nkn,nel)
        call levels_init(nkn,nel,nlv)
        call mod_depth_init(nkn,nel)
	call basin_set_read_basin(.true.)
	call shy_copy_basin_from_shy(id)
	call shy_copy_levels_from_shy(id)

	call shy_proj

	!call test_internal_numbering(id)

	call shympi_init(.false.)		!call after basin has been read
	call shympi_set_hlv(nlv,hlv)

	call ev_set_verbose(.not.bquiet)
        call ev_init(nel)
	call set_ev

	!if( bverb ) write(6,*) 'hlv: ',nlv,hlv

	!--------------------------------------------------------------
	! set dimensions and allocate arrays
	!--------------------------------------------------------------

	bhydro = ftype == 1
	bscalar = ftype == 2

	if( bhydro ) then		!OUS
	  if( nvar /= 4 ) goto 71
	  nndim = 3*nel
	  allocate(il(nel))
	  il = ilhv
	else if( bscalar ) then		!NOS
	  nndim = nkn
	  allocate(il(nkn))
	  il = ilhkv
	else
	  goto 76	!relax later
	end if

        call shy_check_nvar(id,nvar)

	allocate(ieflag(nel))
	allocate(ikflag(nkn))
	allocate(cv2(nndim))
	allocate(cv3(nlv,nndim))
	allocate(cv3all(nlv,nndim,0:nvar))
	allocate(idims(4,nvar))
        allocate(ivars(nvar),strings(nvar))
	allocate(znv(nkn),uprv(nlv,nkn),vprv(nlv,nkn))
	allocate(sv(nlv,nkn),dv(nlv,nkn))
	if( bdiff ) then
	  allocate(cv3diff(nlv,nndim,0:nvar))
	else
	  allocate(cv3diff(1,1,0:1))
	end if

	!--------------------------------------------------------------
	! set up aux arrays, sigma info and depth values
	!--------------------------------------------------------------

	call shyutil_init(nkn,nel,nlv)

	call init_sigma_info(nlv,hlv)

	call shy_make_area
	!call shy_check_area
	call outfile_make_depth(nkn,nel,nen3v,hm3v,hev,hkv)

	!--------------------------------------------------------------
	! write info to terminal
	!--------------------------------------------------------------

	call shy_peek_record(id,dtime,iaux,iaux,iaux,iaux,ierr)
	if( ierr > 0 ) goto 99
	if( ierr < 0 ) goto 98
        call shy_get_string_descriptions(id,nvar,ivars,strings)

	if( bverb ) call depth_stats(nkn,nlvdi,ilhkv)

	if( .not. bquiet ) then
	  call shy_print_descriptions(nvar,ivars,strings)
	end if

	if( binfo ) return

	!--------------------------------------------------------------
	! setup node handling
	!--------------------------------------------------------------

	call initialize_nodes	!single node output

	!--------------------------------------------------------------
	! time averaging
	!--------------------------------------------------------------

	call elabutil_set_averaging(nvar)	!sets btrans and avermode
	call elabutil_check_options		!see if b2d and btrans

	if( btrans ) then
	  call shyutil_init_accum(avermode,nlvdi,nndim,nvar,istep)
	else
	  call shyutil_init_accum(avermode,1,1,1,1)
	end if

	!--------------------------------------------------------------
	! initialize volume
	!--------------------------------------------------------------

	shy_zeta = 0.
	call shy_make_volume

	!--------------------------------------------------------------
	! time handling
	!--------------------------------------------------------------

	call shy_get_date(id,date,time)
	call elabtime_date_and_time(date,time)
	call elabtime_set_minmax(stmin,stmax)
	call elabtime_set_inclusive(binclusive)
	
	!--------------------------------------------------------------
	! read single areas if defined
	!--------------------------------------------------------------

        call handle_area

	!--------------------------------------------------------------
	! open output file
	!--------------------------------------------------------------

	boutput = boutput .or. btrans

	if( bsumvar ) then
	  call shyelab_init_output(id,idout,ftype,1,(/10/))
	else if( bmap ) then
	  call shyelab_init_output(id,idout,ftype,1,(/75/))
	else
	  call shyelab_init_output(id,idout,ftype,nvar,ivars)
	end if

!--------------------------------------------------------------
! loop on data
!--------------------------------------------------------------

	dtime = 0.
	call shy_peek_record(id,dtime,iaux,iaux,iaux,iaux,ierr)
	call dts_convert_to_atime(datetime_elab,dtime,atime)
	it = dtime
	atfirst = atime
	atlast = atime - 1	!do as if atlast has been read
	call custom_dates_init(atime,datefile)

	cv3 = 0.
	cv3all = 0.

	do

	 atold = atime

	 !--------------------------------------------------------------
	 ! read new data set
	 !--------------------------------------------------------------

	 call read_records(id,dtime,nvar,nndim,nlvdi,idims
     +				,cv3,cv3all,ierr)

         if(ierr.ne.0) then	!EOF - see if we have to read another file
	   if( ierr > 0 .or. atstart == -1. ) exit
	   call open_new_file(ifile,id,atstart)
	   if( .not. bsilent ) call shy_write_filename(id)
	   cycle
	 end if

	 call dts_convert_to_atime(datetime_elab,dtime,atime)
	 call dts_format_abs_time(atime,aline)
	 if( concat_cycle_a(atime,atlast,atstart) ) cycle
	 atlast = atime

	 !--------------------------------------------------------------
	 ! handle diffs
	 !--------------------------------------------------------------

	 if( bdiff ) then
	   call read_records(iddiff,ddtime,nvar,nndim,nlvdi,idims
     +				,cv3,cv3diff,ierr)
	   if( ierr /= 0 ) goto 62
	   !if( dtime /= ddtime ) goto 61
	   cv3all = cv3all - cv3diff
	   call check_diff(nlv,nndim,nvar,cv3all,deps,ierr)
	   ndiff = ndiff + ierr
	   if( ierr /= 0 .and. .not. boutput ) goto 60
	 end if

	 nread = nread + 1
	 nrec = nrec + nvar

	 !--------------------------------------------------------------
	 ! look for new record and see if we are in time window
	 !--------------------------------------------------------------

	 call shy_peek_record(id,dtnew,iaux,iaux,iaux,iaux,ierr)
	 if( ierr .ne. 0 ) dtnew = dtime
	 blastrecord = ierr < 0 .and. atstart == -1
	 call dts_convert_to_atime(datetime_elab,dtnew,atnew)

	 if( elabtime_over_time(atime,atnew,atold) ) exit
	 if( .not. elabtime_in_time(atime,atnew,atold) ) cycle

	 call shy_make_zeta(ftype)
	 call shy_make_volume		!comment for constant volume

	 !--------------------------------------------------------------
	 ! initialize record header for output
	 !--------------------------------------------------------------

	 call shyelab_header_output(idout,ftype,dtime,nvar)

	 it = dtime
	 call custom_dates_over(atime,bforce)

	 !--------------------------------------------------------------
	 ! loop over single variables
	 !--------------------------------------------------------------

	 do iv=1,nvar

	  n = idims(1,iv)
	  m = idims(2,iv)
	  lmax = idims(3,iv)
	  ivar = idims(4,iv)
	  nn = n * m

	  cv3(:,:) = cv3all(:,:,iv)


	  if( iv == 1 ) nelab = nelab + 1

	  if( bverb .and. iv == 1 ) then
	    call shy_write_time(.true.,dtime,atime,0)
	  end if

	  if( bwrite ) then
	    call shy_write_min_max(nlvdi,nn,il,lmax,ivar,cv3)
	  end if

	  if( btrans ) then
	    call shy_time_aver(bforce,avermode,iv,nread,ifreq,istep,nndim
     +			,idims,threshold,cv3,boutput,bverb)
	  end if

	  if( b2d ) then
	    call shy_make_vert_aver(idims(:,iv),nndim,cv3,cv2)
	    call shyelab_record_output(id,idout,dtime,ivar,iv,n,m
     +						,1,1,cv2)
	  else if( bsumvar .or. bmap ) then
	    ! only write at end of loop over variables
	  else
	    call shyelab_record_output(id,idout,dtime,ivar,iv,n,m
     +						,lmax,nlvdi,cv3)
	  end if

	  if( baverbas .and. bscalar ) then
	    call shy_assert(nndim==nkn,'shyelab internal error (123)')
	    call shy_make_basin_aver(idims(:,iv),nlv,nndim,cv3,ikflag
     +                          ,cmin,cmax,cmed,cstd,atot,vtot)
	    call shy_write_aver(aline,nvar,iv,ivar
     +				,cmin,cmax,cmed,cstd,atot,vtot)
	  end if

	 end do		!loop on ivar

	 !--------------------------------------------------------------
	 ! finished loop over single variables - handle hydro file
	 !--------------------------------------------------------------

	 if( baverbas .and. bhydro ) then
           call shy_make_hydro_aver(aline,nndim,cv3all,ikflag
     +                  ,znv,uprv,vprv,sv,dv)
	 end if

	 if( bmap ) then
           ivar = 75
	   iv = 1
           call comp_map(nlvdi,nkn,nvar,pthresh,cthresh,cv3all,cv3)
	   call shyelab_record_output(id,idout,dtime,ivar,iv,n,m
     +						,lmax,nlvdi,cv3)
	 end if

	 !if( bsumvar ) then
	 !  ivar = 10
	 !  iv = 1
	 !  cv3 = sum(cv3all,3)
	!write(6,*) 'writing sum output: ',ivar,iv
	 !  call shyelab_record_output(id,idout,dtime,ivar,iv,n,m
   !  +						,nlv,nlvdi,cv3)
	 !end if

	 if( bnodes ) then	!nodal output
           call write_nodes(atime,ftype,nndim,nvar,ivars,cv3all)
	 end if
 
	 ! bsumvar is also handled in here
	 call shyelab_post_output(id,idout,dtime,nvar,n,m,nndim
     +                                  ,lmax,nlvdi,cv3all)

	end do		!time do loop

!--------------------------------------------------------------
! end of loop on data
!--------------------------------------------------------------

!--------------------------------------------------------------
! final write of variables
!--------------------------------------------------------------

	if( btrans ) then
	  do ip=1,istep
	   do iv=1,nvar
	    naccum = naccu(iv,ip)
	    if( naccum > 0 ) then
	      !call shyelab_increase_nwrite	!done in *_output
	      if( bverb ) write(6,*) 'final aver: ',ip,iv,naccum
	      call shy_time_aver(bforce,-avermode,iv,ip,0,istep,nndim
     +			,idims,threshold,cv3,boutput,bverb)
	      n = idims(1,iv)
	      m = idims(2,iv)
	      lmax = idims(3,iv)
	      ivar = idims(4,iv)
	      call shyelab_record_output(id,idout,dtime,ivar,iv,n,m
     +						,lmax,nlvdi,cv3)
	    end if
	   end do
	  end do
	end if

!--------------------------------------------------------------
! write final message
!--------------------------------------------------------------

	if( .not. bsilent ) then

	write(6,*)
	call dts_format_abs_time(atfirst,aline)
	write(6,*) 'first time record: ',aline
	call dts_format_abs_time(atlast,aline)
	write(6,*) 'last time record:  ',aline

	call shyelab_get_nwrite(nwrite,nwtime)

	write(6,*)
	write(6,*) ifile, ' file(s) read'
	write(6,*) nrec,  ' data records read'
	write(6,*) nread, ' time records read'
	write(6,*) nelab, ' time records elaborated'
	write(6,*) nwrite,' data records written'
	write(6,*) nwtime,' time records written'
	write(6,*)

	end if

	call shyelab_final_output(id,idout,nvar,ivars)

	if( bnodes .and. .not. bquiet ) then
	  call write_nodes_final(ftype,nvar,ivars)
	end if

	if( bdiff ) then
	  if( ndiff > 0 ) goto 60
	  write(6,*) 'no difference found > ',deps
	  write(6,*) 'files are identical'
	  call exit(99)	!99 means no difference
	end if

!--------------------------------------------------------------
! end of routine
!--------------------------------------------------------------

	stop
   60	continue
	write(6,*) 'difference found > ',deps
	write(6,*) '*** files are different'
	call exit(1)
	stop 'stop shyelab: difference'
   61	continue
	write(6,*) 'difference of time between files'
	write(6,*) 'time1,time2: ',dtime,ddtime
	stop 'error stop shyelab: time record in diffing'
   62	continue
	write(6,*) 'cannot read record in diff file'
	stop 'error stop shyelab: no record in diffing'
   66	continue
	write(6,*) 'for computing difference need exactly 2 files'
	stop 'error stop shyelab: need 2 files'
   71	continue
	write(6,*) 'ftype = ',ftype,'  nvar = ',nvar
	write(6,*) 'nvar should be 4'
	stop 'error stop shyelab: ftype,nvar'
   74	continue
	stop 'error stop shyelab: general error...'
   75	continue
	write(6,*) 'error writing header, ierr = ',ierr
	write(6,*) 'file = ',trim(file)
	stop 'error stop shyelab: writing header'
   76	continue
	write(6,*) 'ftype = ',ftype,'  expecting 1 or 2'
	stop 'error stop shyelab: ftype'
   77	continue
	write(6,*) 'error reading header, ierr = ',ierr
	write(6,*) 'file = ',trim(file)
	stop 'error stop shyelab: reading header'
   98	continue
	write(6,*) 'error reading file, ierr = ',ierr
	write(6,*) 'file = ',trim(file)
	stop 'error stop shyelab: file contains no data'
   99	continue
	write(6,*) 'error reading file, ierr = ',ierr
	write(6,*) 'file = ',trim(file)
	stop 'error stop shyelab: reading first record'
	end

!***************************************************************
!***************************************************************
!***************************************************************

        subroutine shy_make_hydro_aver(aline,nndim,cv3all,ikflag
     +                  ,znv,uprv,vprv,sv,dv)

        use basin
        use levels
        use mod_depth
	use shyutil

        implicit none

        integer, parameter :: nvar = 4
	character*20 aline
	integer iv
        integer nndim
        integer idims(4,nvar)
        real cv3all(nlvdi,nndim,0:nvar)
	integer ikflag(nkn)
        real znv(nkn)
        real uprv(nlvdi,nkn)
        real vprv(nlvdi,nkn)
        real sv(nlvdi,nkn)
        real dv(nlvdi,nkn)

        integer ivar,idim(4)
        real cmin,cmax,cmed,cstd,atot,vtot

        call prepare_hydro(.true.,nndim,cv3all,znv,uprv,vprv)
        call convert_to_speed(uprv,vprv,sv,dv)

	iv = 1
        ivar = 1
        idim = (/nkn,1,1,ivar/)
        call shy_make_basin_aver(idim,1,nkn,znv,ikflag
     +                          ,cmin,cmax,cmed,cstd,atot,vtot)
        call shy_write_aver(aline,nvar,iv,ivar
     +				,cmin,cmax,cmed,cstd,atot,vtot)

	iv = 2
        ivar = 2
        idim = (/nkn,1,nlv,ivar/)
        call shy_make_basin_aver(idim,nlv,nkn,uprv,ikflag
     +                          ,cmin,cmax,cmed,cstd,atot,vtot)
        call shy_write_aver(aline,nvar,iv,ivar
     +				,cmin,cmax,cmed,cstd,atot,vtot)

	iv = 3
        call shy_make_basin_aver(idim,nlv,nkn,vprv,ikflag
     +                          ,cmin,cmax,cmed,cstd,atot,vtot)
        call shy_write_aver(aline,nvar,iv,ivar
     +				,cmin,cmax,cmed,cstd,atot,vtot)

	iv = 4
        ivar = 6
        idim = (/nkn,1,nlv,ivar/)
        call shy_make_basin_aver(idim,nlv,nkn,sv,ikflag
     +                          ,cmin,cmax,cmed,cstd,atot,vtot)
	!vtot = 0.
        call shy_write_aver(aline,nvar,iv,ivar
     +				,cmin,cmax,cmed,cstd,atot,vtot)

        end

!***************************************************************
!***************************************************************
!***************************************************************

	subroutine shy_split_hydro(id,dtime,znv,uprv,vprv,sv,dv)

	use basin
	use levels

	implicit none

	integer id
	double precision dtime
	real znv(nkn)
	real uprv(nlvdi,nkn)
	real vprv(nlvdi,nkn)
	real sv(nlvdi,nkn)
	real dv(nlvdi,nkn)

	integer, save :: icall = 0
	integer, save :: idz,idu,idv,ids,idd

	if( icall == 0 ) then
	  call shy_split_internal(id,'zeta.shy',.true.,idz)
	  call shy_split_internal(id,'velx.shy',.false.,idu)
	  call shy_split_internal(id,'vely.shy',.false.,idv)
	  call shy_split_internal(id,'speed.shy',.false.,ids)
	  call shy_split_internal(id,'dir.shy',.false.,idd)
	  icall = 1
	end if

	call shy_write_output_record(idz,dtime,1,nkn,1,1,1,znv)
	call shy_write_output_record(idu,dtime,2,nkn,1,nlv,nlv,uprv)
	call shy_write_output_record(idv,dtime,2,nkn,1,nlv,nlv,vprv)
	call shy_write_output_record(ids,dtime,6,nkn,1,nlv,nlv,sv)
	call shy_write_output_record(idd,dtime,7,nkn,1,nlv,nlv,dv)

	end

!***************************************************************

	subroutine shy_split_internal(id,file,b2d,id_out)

	use shyfile

	implicit none

	integer id
	character*(*) file
	logical b2d
	integer id_out		!id of output file (return)

	integer ierr
	integer nk,ne,np,nl,nv

	id_out = shy_init(file)
	if( id_out <= 0 ) goto 99

        call shy_clone(id,id_out)
	call shy_get_params(id_out,nk,ne,np,nl,nv)
	call shy_set_params(id_out,nk,ne,1,nl,1)
	if( b2d ) call shy_convert_2d(id_out)
	call shy_set_ftype(id_out,2)

        call shy_write_header(id_out,ierr)
	if( ierr /= 0 ) goto 98


	return
   98	continue
	stop 'error stop shy_split_internal: cannot write header'
   99	continue
	stop 'error stop shy_split_internal: cannot open file'
	end

!***************************************************************

	subroutine shy_split_id(ivar,id_in,id_out)

	use shyfile

	implicit none

	integer ivar
	integer id_in
	integer id_out

	integer, save :: nsplit = 0
	integer, save, allocatable :: iusplit(:)
	integer, allocatable :: iuaux(:)
	integer ierr,ip
	character*80 name

        if( nsplit == 0 ) then
	  nsplit = max(ivar,100)
	  allocate(iusplit(nsplit))
	  iusplit = 0
	end if

        if( ivar > nsplit ) then
          allocate(iuaux(2*ivar))
	  iuaux = 0
          iuaux(1:nsplit) = iusplit(1:nsplit)
          call move_alloc(iuaux,iusplit)
	  nsplit = 2*ivar
        end if

	if( ivar < 0 ) then	!just check - do not open
	  if( -ivar > nsplit ) then
	    id_out = -1
	  else
	    id_out = iusplit(-ivar)
	  end if
	  return
	end if

	id_out = iusplit(ivar)

	!write(6,*) 'shy_split_id: ',ivar,id_out

	if( id_out == 0 ) then
	  call ivar2filename(ivar,name)
	  name = trim(name) // '.shy'
	  id_out = shy_init(name)
	  if( id_out <= 0 ) goto 99
          call shy_clone(id_in,id_out)
       	  call shy_convert_1var(id_out)
          call shy_write_header(id_out,ierr)
	  if( ierr /= 0 ) goto 98
	  iusplit(ivar) = id_out
	end if

	return
   98	continue
	stop 'error stop shy_split_id: cannot write header'
   99	continue
	stop 'error stop shy_split_id: cannot open file'
	end

!***************************************************************
!***************************************************************
!***************************************************************

	subroutine check_diff(nlv,nn,nvar,cv3all,deps,ndiff)

        implicit none

        integer nlv,nn,nvar
        real cv3all(nlv,nn,0:nvar)
	real deps
        integer ndiff

        real diff
        integer iloc(3)

        diff = maxval( abs(cv3all) )

        if( diff > deps ) then
          ndiff = ndiff + 1
          iloc = maxloc( abs(cv3all) )
          write(6,*) '*** record differing: ',diff,iloc
        end if

	end

!***************************************************************

        subroutine comp_map(nlvddi,nkn,nvar,pt,ct,cvv,valri)

c compute dominant discharge and put index in valri (custom routine)

	use levels

        implicit none

        integer nlvddi,nkn,nvar
        real pt,ct
        real cvv(nlvddi,nkn,0:nvar)
        real valri(nlvddi,nkn)

        integer iv,k,ismax,l,lmax
        real conz, pconz
        real sum,rmax
        real cthresh,pthresh

        pthresh = pt     !threshold on percentage
        cthresh = ct     !threshold on concentration - 0 for everywhere

        !pthresh = 30.
        !cthresh = 20.

	valri = 0.

        do k=1,nkn
	  lmax = ilhkv(k)
          do l=1,lmax
                sum = 0.
                rmax = 0.
                ismax = 0
                do iv=1,nvar
                   conz = cvv(l,k,iv)
                   sum = sum + conz
                   if( conz .gt. rmax ) then
                        rmax = conz
                        ismax = iv
                   end if
                end do

                conz = 0.
                if( ismax .gt. 0 ) conz = cvv(l,k,ismax)
                pconz = 0.
                if( sum .gt. 0. ) pconz = (conz/sum)*100

                valri(l,k) = 0.
                if( conz .gt. cthresh ) then
                  if( pconz .gt. pthresh ) then
                    valri(l,k) = ismax
                  end if
                end if
          end do
        end do

        end

!***************************************************************

        subroutine shy_assert(bval,text)

        logical bval
        character*(*) text

        if( .not. bval ) then
          write(6,*) 'assertion violated'
          write(6,*) text
          stop 'error stop shy_assert: assertion violated'
        end if

        end subroutine shy_assert

!***************************************************************

	subroutine shy_write_filename(id)

	use shyfile

	integer id

        character*80 file

	call shy_get_filename(id,file)

        write(6,*) '================================'
        write(6,*) 'new file: ',trim(file)
        write(6,*) '================================'

	end subroutine shy_write_filename

!***************************************************************

	subroutine test_internal_numbering(id)

	use basin

	implicit none

	integer id
	integer i

	write(6,*) 'test_internal_numbering: '

	do i=1,nel,nel/10
	  write(6,*) i,ipev(i)
	end do

	do i=1,nkn,nkn/10
	  write(6,*) i,ipv(i)
	end do

	stop

	end

!***************************************************************

	subroutine shy_check_area

	use shyutil

	implicit none

	real area_k,area_e

	area_k = sum(areak)
	area_e = sum(areae)
	write(6,*) 'areas: ',area_k,area_e

	end

!***************************************************************

