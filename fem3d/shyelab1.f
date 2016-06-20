!
! $Id: noselab.f,v 1.8 2008-11-20 10:51:34 georg Exp $
!
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
!
!**************************************************************

	subroutine shyelab1

	use clo
	use elabutil
	use elabtime
	use shyfile
	use shyutil

        use basin
        use mod_depth
        use evgeom
        use levels

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
	logical blastrecord
	integer nwrite,nread,nelab,nrec,nin,nold,ndiff
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
	real rnull
	real cmin,cmax,cmed,vtot
	double precision dtime,dtstart,dtnew,ddtime
	double precision atime,atstart,atnew,atold

	integer iapini
	integer ifem_open_file

!--------------------------------------------------------------
! initialize everything
!--------------------------------------------------------------

	nread=0
	nwrite=0
	nelab=0
	nrec=0
	ndiff = 0
	rnull=0.
	rnull=-1.
	bopen = .false.
	bzeta = .false.		!file has zeta information
	ifile = 0
	id = 0

	!--------------------------------------------------------------
	! set command line parameters
	!--------------------------------------------------------------

	call elabutil_init('SHY')

	!--------------------------------------------------------------
	! open input files
	!--------------------------------------------------------------

	call open_new_file(ifile,id,atstart)	!atstart=-1 if no new file

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

	call shy_info(id)
	!if( bdiff ) call shy_info(iddiff)

        call basin_init(nkn,nel)
        call levels_init(nkn,nel,nlv)
        call mod_depth_init(nkn,nel)
	call shy_copy_basin_from_shy(id)
	call shy_copy_levels_from_shy(id)
        call ev_init(nel)
	call set_ev

	if( bverb ) write(6,*) 'hlv: ',nlv,hlv

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
	call outfile_make_depth(nkn,nel,nen3v,hm3v,hev,hkv)

	if( bverb ) call depth_stats(nkn,nlvdi,ilhkv)

	!--------------------------------------------------------------
	! setup node handling
	!--------------------------------------------------------------

	call handle_nodes	!single node output

	!--------------------------------------------------------------
	! time averaging
	!--------------------------------------------------------------

	call elabutil_set_averaging(nvar)	!sets btrans
	call elabutil_check_options		!see if b2d and btrans

	if( btrans ) then
	  call shyutil_init_accum(nlvdi,nndim,nvar,istep)
	else
	  call shyutil_init_accum(1,1,1,1)
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
	call elabtime_minmax(stmin,stmax)
	call elabtime_set_inclusive(binclusive)
	
	!--------------------------------------------------------------
	! open output file
	!--------------------------------------------------------------

	boutput = boutput .or. btrans .or. bsplit .or. bsumvar

	call shy_peek_record(id,dtime,iaux,iaux,iaux,iaux,ierr)
	if( ierr > 0 ) goto 99
	if( ierr < 0 ) goto 98
        call shy_get_string_descriptions(id,nvar,ivars,strings)
	call shyelab_init_output(id,idout,nvar,ivars)

	!--------------------------------------------------------------
	! write info to terminal
	!--------------------------------------------------------------

        write(6,*) 'available variables contained in file: '
        write(6,*) 'total number of variables: ',nvar
        write(6,*) '   varnum     varid    varname'
        do iv=1,nvar
          ivar = ivars(iv)
          write(6,'(2i10,4x,a)') iv,ivar,trim(strings(iv))
        end do

!--------------------------------------------------------------
! loop on data
!--------------------------------------------------------------

	dtime = 0.
	call shy_peek_record(id,dtime,iaux,iaux,iaux,iaux,ierr)
	call dts_convert_to_atime(datetime_elab,dtime,atime)

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
	   cycle
	 end if

	 call dts_convert_to_atime(datetime_elab,dtime,atime)

	 !--------------------------------------------------------------
	 ! handle diffs
	 !--------------------------------------------------------------

	 if( bdiff ) then
	   call read_records(iddiff,ddtime,nvar,nndim,nlvdi,idims
     +				,cv3,cv3diff,ierr)
	   if( ierr /= 0 ) goto 62
	   if( dtime /= ddtime ) goto 61
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
	 !call shy_make_volume		!comment for constant volume

	 !--------------------------------------------------------------
	 ! initialize record header for output
	 !--------------------------------------------------------------

	 call shyelab_header_output(id,idout,dtime,nvar)

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

	  nelab = nelab + 1

	  if( .not. bquiet ) then
	    call shy_write_time(.true.,dtime,atime,ivar)
	  end if

	  if( bwrite ) then
	    call shy_write_min_max(nlvdi,nn,lmax,cv3)
	  end if

	  if( btrans ) then
	    call shy_time_aver(mode,iv,nread,ifreq,istep,nndim
     +			,idims,threshold,cv3,boutput)
	  end if

	  if( b2d ) then
	    call shy_make_vert_aver(idims(:,iv),nndim,cv3,cv2)
	    call shyelab_record_output(id,idout,dtime,ivar,iv,n,m
     +						,1,1,cv2)
	  else
	    call shyelab_record_output(id,idout,dtime,ivar,iv,n,m
     +						,nlv,nlvdi,cv3)
	  end if

	  if( baverbas .and. bscalar ) then
	    call shy_make_basin_aver(idims(:,iv),nndim,cv3
     +                          ,cmin,cmax,cmed,vtot)
	    call shy_write_aver(dtime,ivar,cmin,cmax,cmed,vtot)
	  end if

	  if( bnodes .and. bscalar ) then	!scalar output
	    call write_nodes(dtime,ivar,cv3)
	  end if

	 end do		!loop on ivar

	 !--------------------------------------------------------------
	 ! finished loop over single variables - handle hydro file
	 !--------------------------------------------------------------

	 if( baverbas .and. bhydro ) then
           call shy_make_hydro_aver(dtime,nndim,cv3all
     +                  ,znv,uprv,vprv,sv,dv)
	 end if

	 if( bnodes ) then	!nodal output
	   if( bhydro ) then	!hydro output
	     call prepare_hydro(.true.,nndim,cv3all,znv,uprv,vprv)
	     call write_nodes_vel(dtime,znv,uprv,vprv)
	     call write_nodes_hydro(dtime,znv,uprv,vprv)
	   else if( bscalar ) then
	     call write_nodes_scalar(dtime,nvar,nndim,strings,cv3all)
	   end if
	 end if
 
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
	      nwrite = nwrite + 1
	      write(6,*) 'final aver: ',ip,iv,naccum
	      call shy_time_aver(-mode,iv,ip,0,istep,nndim
     +			,idims,threshold,cv3,boutput)
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

	write(6,*)
	write(6,*) nrec,  ' records read'
	write(6,*) nread, ' unique time records read'
	write(6,*) nelab, ' records elaborated'
	write(6,*) ifile, ' file(s) read'
	write(6,*) nwrite,' records written'
	write(6,*)

	call shyelab_final_output(id,idout,nvar)

	if( bnodes ) then
	  write(6,*) 'output of node(s) has been written to out.fem'
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

	subroutine convert_to_speed(uprv,vprv,sv,dv)

	use basin
	use levels
	
	implicit none

	real uprv(nlvdi,nkn)
	real vprv(nlvdi,nkn)
	real sv(nlvdi,nkn)
	real dv(nlvdi,nkn)

	integer k,lmax,l
	real u,v,s,d

        do k=1,nkn
          lmax = ilhkv(k)
          do l=1,lmax
            u = uprv(l,k)
            v = vprv(l,k)
            call c2p(u,v,s,d)   !d is meteo convention
            d = d + 180.
            if( d > 360. ) d = d - 360.
            sv(l,k) = s
            dv(l,k) = d
          end do
        end do

	end

!***************************************************************

	subroutine prepare_hydro(bvel,nndim,cv3all,znv,uprv,vprv)

	use basin
	use levels
	use mod_depth
	
	implicit none

	logical bvel
	integer nndim
	real cv3all(nlvdi,nndim,0:4)
	real znv(nkn)
	real uprv(nlvdi,nkn)
	real vprv(nlvdi,nkn)

	real, allocatable :: zenv(:)
	real, allocatable :: uv(:,:)
	real, allocatable :: vv(:,:)

	allocate(zenv(3*nel))
	allocate(uv(nlvdi,nel))
	allocate(vv(nlvdi,nel))

        znv(1:nkn)     = cv3all(1,1:nkn,1)
        zenv(1:3*nel)  = cv3all(1,1:3*nel,2)
        uv(:,1:nel)    = cv3all(:,1:nel,3)
        vv(:,1:nel)    = cv3all(:,1:nel,4)

	call shy_transp2vel(bvel,nel,nkn,nlv,nlvdi,hev,zenv,nen3v
     +                          ,ilhv,hlv,uv,vv
     +                          ,uprv,vprv)

	deallocate(zenv,uv,vv)

	end

!***************************************************************

        subroutine shy_make_hydro_aver(dtime,nndim,cv3all
     +                  ,znv,uprv,vprv,sv,dv)

        use basin
        use levels
        use mod_depth

        implicit none

        integer, parameter :: nvar = 4
        double precision dtime
        integer nndim
        integer idims(4,nvar)
        real cv3all(nlvdi,nndim,0:nvar)
        real znv(nkn)
        real uprv(nlvdi,nkn)
        real vprv(nlvdi,nkn)
        real sv(nlvdi,nkn)
        real dv(nlvdi,nkn)

        integer ivar,idim(4)
        real cmin,cmax,cmed,vtot

        call prepare_hydro(.true.,nndim,cv3all,znv,uprv,vprv)
        call convert_to_speed(uprv,vprv,sv,dv)

        ivar = 1
        idim = (/nkn,1,1,ivar/)
        call shy_make_basin_aver(idim,nkn,znv
     +                          ,cmin,cmax,cmed,vtot)
	vtot = 0.
        call shy_write_aver(dtime,ivar,cmin,cmax,cmed,vtot)

        ivar = 2
        idim = (/nkn,1,nlv,ivar/)
        call shy_make_basin_aver(idim,nkn,uprv
     +                          ,cmin,cmax,cmed,vtot)
	vtot = 0.
        call shy_write_aver(dtime,ivar,cmin,cmax,cmed,vtot)

        call shy_make_basin_aver(idim,nkn,vprv
     +                          ,cmin,cmax,cmed,vtot)
	vtot = 0.
        call shy_write_aver(dtime,ivar,cmin,cmax,cmed,vtot)

        ivar = 6
        idim = (/nkn,1,nlv,ivar/)
        call shy_make_basin_aver(idim,nkn,sv
     +                          ,cmin,cmax,cmed,vtot)
	vtot = 0.
        call shy_write_aver(dtime,ivar,cmin,cmax,cmed,vtot)

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
	  call shy_split_internal(id,'z.shy',.true.,idz)
	  call shy_split_internal(id,'u.shy',.false.,idu)
	  call shy_split_internal(id,'v.shy',.false.,idv)
	  call shy_split_internal(id,'s.shy',.false.,ids)
	  call shy_split_internal(id,'d.shy',.false.,idd)
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
	  write(name,'(i4,a)') ivar,'.shy'
	  id_out = shy_init(adjustl(name))
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

