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
!
!**************************************************************

	program shyplot

	use plotutil

	implicit none

	call plotutil_init('SHY')
	call classify_files

	if( shyfilename /= ' ' ) then
	  call plot_shy_file
	else if( basfilename /= ' ' ) then
	  call plot_bas_file
	end if

	end

!**************************************************************

	subroutine plot_bas_file

        use mod_depth
        use mod_geom
        use evgeom
        use basin
        use plotutil

	implicit none

	integer ivar

	write(6,*) 'not yet ready...'

        call read_command_line_file(basfilename)

        call ev_init(nel)
        call set_ev

        call mod_geom_init(nkn,nel,ngr)
        call set_geom

        call mod_depth_init(nkn,nel)

        call makehev(hev)
        call makehkv(hkv)
        call allocate_2d_arrays(nel)

        call init_plot

	ivar = 5			!bathymetry
	call init_nls_fnm
	call read_str_files(-1)
	call read_str_files(ivar)
        call initialize_color

        call qopen
	call plobas
        call qclose

	end

!**************************************************************

	subroutine plot_shy_file

	use clo
	!use elabutil
	use plotutil
	use elabtime
	use shyfile
	use shyutil

        use basin
        use levels
        use evgeom
        use mod_depth
	use mod_hydro_plot
        use mod_hydro
        use mod_hydro_vel
        use mod_hydro_print

	implicit none

	real, allocatable :: cv2(:)
	real, allocatable :: cv3(:,:)
	real, allocatable :: cv3all(:,:,:)

	integer, allocatable :: idims(:,:)
	integer, allocatable :: il(:)
	integer, allocatable :: ivars(:)
	character*80, allocatable :: strings(:)

	!real, allocatable :: znv(:)
	!real, allocatable :: uprv(:,:)
	!real, allocatable :: vprv(:,:)
	real, allocatable :: sv(:,:)
	real, allocatable :: dv(:,:)

	logical bhydro,bscalar
	logical blastrecord
	integer nwrite,nread,nelab,nrec,nin,nold
	integer nvers
	integer nvar,npr
	integer ierr,ivel,iarrow
	!integer it,itvar,itnew,itold,itstart
	integer it
	integer ivar,iaux,nv
	integer iv,j,l,k,lmax,node
	integer ip
	integer ifile,ftype
	integer id,idout,idold
	integer n,m,nndim,nn
	integer naccum
	integer isphe
	integer date,time
	character*80 title,name,file
	character*80 basnam,simnam,varline
	real rnull
	real cmin,cmax,cmed,vtot
	double precision dtime,dtstart,dtnew
	double precision atime,atstart,atnew,atold

	integer iapini
	integer ifem_open_file
	real getpar

!--------------------------------------------------------------
! initialize everything
!--------------------------------------------------------------

	nread=0
	nwrite=0
	nelab=0
	nrec=0
	rnull=0.
	rnull=-1.
	bopen = .false.
	bzeta = .false.		!file has zeta information
	ifile = 0
	id = 0
	idold = 0

	!--------------------------------------------------------------
	! set command line parameters
	!--------------------------------------------------------------

	call init_nls_fnm
	call read_str_files(-1)
	call read_str_files(ivar3)
	write(6,*) 'icolor: ',nint(getpar('icolor'))

	!--------------------------------------------------------------
	! open input files
	!--------------------------------------------------------------

	call open_next_file_by_name(shyfilename,idold,id)
	if( id == 0 ) stop

	!--------------------------------------------------------------
	! set up params and modules
	!--------------------------------------------------------------

	call shy_get_params(id,nkn,nel,npr,nlv,nvar)
	call shy_get_ftype(id,ftype)

	call shy_info(id)

        call basin_init(nkn,nel)
        call levels_init(nkn,nel,nlv)

	call shy_copy_basin_from_shy(id)
	call shy_copy_levels_from_shy(id)

        call mod_depth_init(nkn,nel)
	call allocate_2d_arrays(nel)
	call allocate_simulation(0)
	call mod_hydro_plot_init(nkn,nel)

        isphe = nint(getpar('isphe'))
        call set_coords_ev(isphe)
        call set_ev
        call set_geom
        call get_coords_ev(isphe)
        call putpar('isphe',float(isphe))

	!--------------------------------------------------------------
	! set time
	!--------------------------------------------------------------

        call ptime_init
	call shy_get_date(id,date,time)
        call ptime_set_date_time(date,time)

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
	!allocate(znv(nkn),uprv(nlv,nkn),vprv(nlv,nkn))
	allocate(sv(nlv,nkn),dv(nlv,nkn))

	!--------------------------------------------------------------
	! set up aux arrays, sigma info and depth values
	!--------------------------------------------------------------

	call shyutil_init(nkn,nel,nlv)

	call init_sigma_info(nlv,hlv)

	call shy_make_area
	call outfile_make_depth(nkn,nel,nen3v,hm3v,hev,hkv)

	if( bverb ) call depth_stats(nkn,nlvdi,ilhkv)

	!--------------------------------------------------------------
	! initialize plot
	!--------------------------------------------------------------

	call init_plot

	allocate(ivars(nvar),strings(nvar))
	call shy_get_string_descriptions(id,nvar,ivars,strings)
	write(6,*) 'available variables: ',nvar
	nv = 0
	do iv=1,nvar
	  ivar = ivars(iv)
	  write(6,*) iv,ivar,trim(strings(iv))
	  if( ivar3 == ivars(iv) ) nv = nv + 1
	  if( ivar3 == 2 .and. ivars(iv) == 3 ) nv = nv + 1
	end do

        if( ivar3 > 0 ) then
	  if( nv == 0 ) then
	    call ivar2string(ivar3,varline)
            write(6,*) 'no such variable in file: ',ivar3,varline
            stop 'error stop shyplot'
	  end if
	else if( nvar == 1 ) then
	  ivar3 = ivars(1)
	  call read_str_files(ivar3)
        else
          write(6,*) 'no variable given to be plotted: ',ivar3
          stop 'error stop shyplot'
        end if
	if( layer > nlv ) then
          write(6,*) 'no such layer: ',layer
          write(6,*) 'maximum layer available: ',nlv
          stop 'error stop shyplot'
	end if

	call mkvarline(ivar3,varline)
	write(6,*) 'varline: ',trim(varline)
	write(6,*) 'ivar3: ',ivar3
	write(6,*) 'layer: ',layer
	call setlev(layer)

	ivel = 0
	if( ivar3 == 2 ) ivel = 1
	if( ivar3 == 3 ) ivel = 2
	if( ivel > 0 .and. nv /= 2 ) then
	  write(6,*) 'could not read all data for arrow: ',ivel,nv
          stop 'error stop shyplot'
	end if

	!--------------------------------------------------------------
	! initialize volume
	!--------------------------------------------------------------

	shy_zeta = 0.
	call shy_make_volume

	!--------------------------------------------------------------
	! initialize plot
	!--------------------------------------------------------------

        call initialize_color

	call qopen

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

         if(ierr.ne.0) exit

	 call dts_convert_to_atime(datetime_elab,dtime,atime)

	 iarrow = 0
	 nread = nread + 1
	 nrec = nrec + nvar

	 !--------------------------------------------------------------
	 ! look for new record and see if we are in time window
	 !--------------------------------------------------------------

	 call shy_peek_record(id,dtnew,iaux,iaux,iaux,iaux,ierr)
	 if( ierr .ne. 0 ) dtnew = dtime
	 blastrecord = ierr < 0 .and. atstart == -1
	 call dts_convert_to_atime(datetime_elab,dtnew,atnew)

	 if( elabtime_over_time_a(atime,atnew,atold) ) exit
	 if( .not. elabtime_check_time_a(atime,atnew,atold) ) cycle

	 call shy_make_zeta(ftype)
	 !call shy_make_volume		!comment for constant volume

	 if( ifreq > 0 .and. mod(nread,ifreq) == 0 ) cycle

	 call ptime_set_dtime(dtime)

	 !--------------------------------------------------------------
	 ! loop over single variables
	 !--------------------------------------------------------------

	 do iv=1,nvar

	  n = idims(1,iv)
	  m = idims(2,iv)
	  lmax = idims(3,iv)
	  ivar = idims(4,iv)
	  nn = n * m

	  if( ivar /= ivar3 .and. ivel == 0 ) cycle
	  if( ivar == 1 .and. m == 3 ) cycle	!water level in element
	  if( ivar == 1 .and. ivel > 0 ) cycle	!want to plot vel/trans

	  cv3(:,:) = cv3all(:,:,iv)

	  nelab = nelab + 1

	  if( .not. bquiet ) then
	    call shy_write_time(.true.,dtime,atime,ivar)
	  end if

	  if( b2d ) then
	    call shy_make_vert_aver(idims(:,iv),nndim,cv3,cv2)
	  else
	    if( m /= 1 ) stop 'error stop: m/= 1'
	    if( n == nkn .and. ivar == 1 ) then
	      cv2(:) = cv3(1,:)
	    else if( n == nkn ) then
	      call extnlev(layer,nlvdi,nkn,cv3,cv2)
	    else if( n == nel ) then
	      call extelev(layer,nlvdi,nel,cv3,cv2)
	    else
	      write(6,*) 'n,nkn,nel: ',n,nkn,nel
	      stop 'error stop: n'
	    end if
	  end if

	  if( ivel > 0 .and. ivar == 3 ) then
	    iarrow = iarrow + 1
	    if( iarrow == 1 ) utrans(1:nel) = cv2(1:nel)
	    if( iarrow == 2 ) vtrans(1:nel) = cv2(1:nel)
	    !if( iarrow == 2 ) call make_vertical_velocity
	    if( iarrow == 2 ) wsnv = 0.	!FIXME - no vertical velocity
	  end if
	  if( ivel > 0 .and. iarrow /= 2 ) cycle

	  write(6,*) 'plotting: ',ivar,layer

          !call prepare_dry_mask
	  !call reset_dry_mask
	  call make_mask(layer)
	  if( m /= 1 ) stop 'error stop: m/= 1'
	  if( ivel > 0 ) then
	    call plo2vel(ivel,'3D ')
	  else if( n == nkn ) then
            call ploval(nkn,cv2,varline)
	  else if( n == nel ) then
            call ploeval(nel,cv2,varline)
	  else
	    write(6,*) 'n,nkn,nel: ',n,nkn,nel
	    stop 'error stop: n'
	  end if

	 end do		!loop on ivar

	 !--------------------------------------------------------------
	 ! finished loop over single variables - handle hydro file
	 !--------------------------------------------------------------

	end do		!time do loop

!--------------------------------------------------------------
! end of loop on data
!--------------------------------------------------------------

	call qclose
	write(6,*) 'total number of plots: ',nelab

!--------------------------------------------------------------
! final write of variables
!--------------------------------------------------------------

!--------------------------------------------------------------
! end of routine
!--------------------------------------------------------------

	stop
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
	end

!***************************************************************
!***************************************************************
!***************************************************************


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

	subroutine make_mask(level)

	use levels
        use mod_plot2d
        !use mod_hydro

	implicit none

	integer level

        call reset_dry_mask

        call set_level_mask(bwater,ilhv,level)        !element has this level
        call make_dry_node_mask(bwater,bkwater)       !copy elem to node mask

        call adjust_no_plot_area
        call make_dry_node_mask(bwater,bkwater) !copy elem to node mask
        call info_dry_mask(bwater,bkwater)

	end

!***************************************************************
!***************************************************************
!***************************************************************

        subroutine initialize_color

        implicit none

        integer icolor
        real getpar

        call colsetup
        icolor = nint(getpar('icolor'))
        call set_color_table( icolor )
        call set_default_color_table( icolor )

        end

!***************************************************************

c*****************************************************************

        subroutine allocate_2d_arrays(npd)

        use mod_hydro_plot
        use mod_plot2d
        use mod_geom
        use mod_depth
        use evgeom
        use basin, only : nkn,nel,ngr,mbw

        implicit none

        integer npd

        integer np

        np = max(nel,npd)

        call ev_init(nel)
        call mod_geom_init(nkn,nel,ngr)

        call mod_depth_init(nkn,nel)

        call mod_plot2d_init(nkn,nel,np)
        call mod_hydro_plot_init(nkn,nel)

        write(6,*) 'allocate_2d_arrays: ',nkn,nel,ngr,np

        end

c*****************************************************************

        subroutine read_command_line_file(file)

        use basin
        !use basutil

        implicit none

        character*(*) file
        logical is_grd_file

        if( basin_is_basin(file) ) then
          write(6,*) 'reading BAS file ',trim(file)
          call basin_read(file)
          !breadbas = .true.
        else if( is_grd_file(file) ) then
          write(6,*) 'reading GRD file ',trim(file)
          call grd_read(file)
          call grd_to_basin
          call estimate_ngr(ngr)
          !breadbas = .false.
        else
          write(6,*) 'Cannot read this file: ',trim(file)
          stop 'error stop read_given_file: format not recognized'
        end if

        end

c*****************************************************************

