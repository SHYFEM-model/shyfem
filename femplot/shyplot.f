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
! 10.06.2016    ggu     shyplot now plots fem files
! 13.06.2016    ggu     shyplot now plots barotropic vars (layer==0)
! 31.10.2016    ggu     shyplot restructured... directional plot still broken
! 14.02.2017    ggu     bug fix in plotting regular fem files - introduced il
!
!**************************************************************

	program shyplot

	use plotutil

	implicit none

	call plotutil_init('SHY')
	call classify_files

	if( shyfilename /= ' ' ) then
	  call plot_shy_file
	else if( femfilename /= ' ' ) then
	  call plot_fem_file
	else if( basfilename /= ' ' ) then
	  call plot_bas_file
	end if

	end

!**************************************************************

	subroutine plot_bas_file

        use mod_depth
        use mod_hydro_plot
        use mod_geom
        use evgeom
        use basin
        use plotutil

	implicit none

	integer ivar

        call read_command_line_file(basfilename)

        call ev_init(nel)
        call set_ev

        call mod_geom_init(nkn,nel,ngr)
        call set_geom

        call mod_depth_init(nkn,nel)
        call mod_hydro_plot_init(nkn,nel,1,nel)

        call makehev(hev)
        call makehkv(hkv)
        call allocate_2d_arrays

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

	logical bhydro,bscalar,bsect,bvect,bvel
	integer irec,nplot,nread,nin,nold
	integer nvers
	integer nvar,npr
	integer ierr,ivel,iarrow
	integer it
	integer ivar,iaux,nv
	integer ivarplot(2),ivs(2)
	integer iv,i,j,l,k,lmax,node
	integer ip,np
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
	double precision dtime
	double precision atime

	integer iapini
	integer ifem_open_file
	integer getisec
	real getpar

!--------------------------------------------------------------
! initialize everything
!--------------------------------------------------------------

	irec=0
	nplot=0
	nread=0
	rnull=0.
	rnull=-1.
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
	call allocate_2d_arrays
	np = nel
	call mod_hydro_plot_init(nkn,nel,nlv,np)

	call mod_hydro_init(nkn,nel,nlv)
	call mod_hydro_vel_init(nkn,nel,nlv)
	call mod_hydro_print_init(nkn,nlv)
	znv = 0.
	zenv = 0.

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
        call elabtime_date_and_time(date,time)
        call elabtime_minmax(stmin,stmax)
	call elabtime_set_inclusive(.false.)

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

	call iff_init_global_2d(nkn,nel,hkv,hev,date,time)

	call init_plot

	call shy_get_string_descriptions(id,nvar,ivars,strings)
	call choose_var(nvar,ivars,strings,varline,ivarplot,bvect)

	bsect = getisec() /= 0
	call setlev(layer)

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
! loop on data (time loop)
!--------------------------------------------------------------

	dtime = 0.
	cv3 = 0.
	cv3all = 0.

	do

	 !--------------------------------------------------------------
	 ! read new data set
	 !--------------------------------------------------------------

	 call read_records(id,dtime,nvar,nndim,nlvdi,idims
     +				,cv3,cv3all,ierr)

         if(ierr.ne.0) exit

	 call dts_convert_to_atime(datetime_elab,dtime,atime)

	 if( .not. bquiet ) call shy_write_time(.true.,dtime,atime,0)

	 irec = irec + 1
	 nread = nread + nvar

	 !--------------------------------------------------------------
	 ! see if we are in time window
	 !--------------------------------------------------------------

	 if( elabtime_over_time(atime) ) exit
	 if( .not. elabtime_in_time(atime) ) cycle

	 call shy_make_zeta(ftype)
	 !call shy_make_volume		!comment for constant volume

	 if( ifreq > 0 .and. mod(irec,ifreq) /= 0 ) cycle

	 call ptime_set_dtime(dtime)

	 !--------------------------------------------------------------
	 ! find out whay to plot
	 !--------------------------------------------------------------

	 ivars(:) = idims(4,:)
	 call get_vars_to_plot(nvar,ivars,ivar3,ivarplot,bvect,ivnum,ivs)

	 ivar = ivarplot(1)
	 iv = ivs(1)
	 if( iv == 0 ) then
	   write(6,*) 'no such variable in time record: ',ivar,iv
	   cycle
	 end if

	 !--------------------------------------------------------------
	 ! extract variables
	 !--------------------------------------------------------------

	 if( bhydro ) then
	   znv(1:nkn) = cv3all(1,1:nkn,1)
	   zenv = reshape(cv3all(1,:,2),(/3,nel/))
	 end if

	 do i=1,2

	 iv = ivs(i)
	 if( iv == 0 ) cycle		!this happens for scalar
	 n = idims(1,iv)
	 m = idims(2,iv)
	 lmax = idims(3,iv)
	 ivar = idims(4,iv)
	 nn = n * m
	 if( m /= 1 ) stop 'error stop: m/= 1'
	 if( n /= nkn .and. n /= nel ) then
	      write(6,*) 'n,nkn,nel: ',n,nkn,nel
	      stop 'error stop: n'
	 end if

	 cv3(:,:) = cv3all(:,:,iv)
	 if( b2d ) then
	    call shy_make_vert_aver(idims(:,iv),nndim,cv3,cv2)
	 else if( lmax == 1 ) then
	   cv2(:) = cv3(1,:)
	 else
	   call extlev(layer,nlvdi,n,il,cv3,cv2)
	 end if

	 if( bvect ) then
	   call directional_insert(bvect,ivar,ivar3,ivarplot,n,cv2,ivel)
	   if( bsect ) then
	      if( iv == 3 ) utlnv(:,1:nel) = cv3(:,1:nel)
	      if( iv == 4 ) vtlnv(:,1:nel) = cv3(:,1:nel)
	   end if
	 end if

	 end do

	 !------------------------------------------
	 ! from here plotting
	 !------------------------------------------

	 if( .not. bvect ) ivel = 0
	 bvel = ivel > 0
	 if( bvect .and. .not. bvel ) stop 'error stop: ivel == 0'

	 nplot = nplot + 1

	 if( .not. bquiet ) then
	   write(6,*) 'plotting: ',ivar,layer,n,ivel
	   write(6,*) 'plotting: ',ivar3,ivarplot
	   call shy_write_time(.true.,dtime,atime,ivar)
	 end if

	 call make_mask(layer)

	 if( bsect ) then
	   call plot_sect(bvel,cv3)
	 else if( bvel ) then
	   call plo2vel(ivel,'3D ')
	 else
           call plo_scal_val(n,cv2,varline)
	 end if

	 !------------------------------------------
	 ! end of plotting
	 !------------------------------------------

	end do		!time do loop

!--------------------------------------------------------------
! end of loop on data (time loop)
!--------------------------------------------------------------

	call qclose
	write(6,*) 'total number of plots: ',nplot

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

	subroutine plot_fem_file

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
        use mod_geom
        use mod_hydro_plot

	implicit none

	logical bhasbasin,breg,bvect,bvel
	logical bsect,bskip
	logical bintp,bplotreg
	integer i,ierr,iformat,irec,l,nplot,ivel
	integer isphe,iunit,lmax,lmax0,np,np0,npaux
	integer ntype,nvar,nvar0,nvers
	integer date,time
	integer datetime(2)
	integer itype(2)
	integer ivarplot(2),ivs(2)
	real regpar(7)
	real flag
	double precision dtime,atime,atime0
	character*80 line,string,varline
        real,allocatable :: data2d(:)
        real,allocatable :: data3d(:,:)
        real,allocatable :: data2ddir(:)
        real,allocatable :: data3ddir(:,:)
        real,allocatable :: data(:,:,:)
        real,allocatable :: dext(:)
        real,allocatable :: hd(:)
        integer,allocatable :: il(:)
        !real,allocatable :: hlv(:)
        !integer,allocatable :: ilhkv(:)
        integer,allocatable :: ivars(:)
        character*80, allocatable :: strings(:)

	integer getisec
	real getpar

        !--------------------------------------------------------------
        ! set command line parameters
        !--------------------------------------------------------------

        call init_nls_fnm
        call read_str_files(-1)
        call read_str_files(ivar3)

        !--------------------------------------------------------------
        ! open input files
        !--------------------------------------------------------------

	infile = femfilename
        if( infile .eq. ' ' ) stop

        np = 0
        call fem_file_read_open(infile,np,iformat,iunit)
        if( iunit .le. 0 ) stop

        write(6,*) 'file name: ',infile(1:len_trim(infile))
        call fem_file_get_format_description(iformat,line)
        write(6,*) 'format: ',iformat,"  (",line(1:len_trim(line)),")"

        !--------------------------------------------------------------
        ! set up params and modules
        !--------------------------------------------------------------



	!--------------------------------------------------------------
	! read first record
	!--------------------------------------------------------------

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

	nlv = lmax
        call levels_init(np,np,nlv)	!first call - will be changed later
        call fem_file_make_type(ntype,2,itype)

        call fem_file_read_2header(iformat,iunit,ntype,lmax
     +                  ,hlv,regpar,ierr)
        if( ierr .ne. 0 ) goto 98

        if( lmax > 1 .and. .not. bquiet ) then
          write(6,*) 'vertical layers: ',lmax
          write(6,*) hlv
        end if
        if( itype(1) .gt. 0 .and. .not. bquiet ) then
          write(6,*) 'date and time: ',datetime
        end if
	breg = ( itype(2) > 0 )
	if( breg ) dflag = regpar(7)
        if( breg .and. .not. bquiet ) then
          write(6,*) 'regpar: ',regpar
        end if

        !--------------------------------------------------------------
        ! configure basin
        !--------------------------------------------------------------

	bhasbasin = basfilename /= ' '

	if( bhasbasin ) then
          call read_command_line_file(basfilename)
	else if( breg ) then
	  call bas_insert_regular(regpar)
	else	!should not be possible
	  write(6,*) 'internal error: ',bhasbasin,breg
	  stop 'error stop plot_fem_file: internal error (7)'
	end if

	if( bhasbasin .or. breg ) then
          call ev_init(nel)
          isphe = nint(getpar('isphe'))
          call set_coords_ev(isphe)
          call set_ev
          call get_coords_ev(isphe)
          call putpar('isphe',float(isphe))

          call mod_geom_init(nkn,nel,ngr)
          call set_geom

          call levels_init(nkn,nel,nlv)
          call mod_depth_init(nkn,nel)

	  npaux = nel
	  call mod_hydro_plot_init(nkn,nel,nlv,npaux)

          call makehev(hev)
          call makehkv(hkv)
          call allocate_2d_arrays
	end if

	if( breg ) then		!data is in regular format
	  bplotreg = .true.
	  bintp = .false.
	  if( bhasbasin ) bintp = .true.
	  if( bregall ) bintp = .false.
	else			!data is 1 value, BC or on nodes
	  bplotreg = .false.
	  bintp = .true.
	  if( bhasbasin ) then
	    if( nkn /= np ) goto 93
	    if( basintype /= 'bas' ) goto 92
	  else
	    !goto 94
	  end if
	end if

        nvar0 = nvar
        lmax0 = lmax
        np0 = np
        allocate(strings(nvar))
        allocate(ivars(nvar))
        allocate(dext(nvar))
        allocate(data2d(np))
        allocate(data3d(lmax,np))
        allocate(data2ddir(np))
        allocate(data3ddir(lmax,np))
        allocate(data(lmax,np,nvar))
        allocate(hd(np))
        allocate(il(np))

	!--------------------------------------------------------------
	! choose variable to plot
	!--------------------------------------------------------------

        do i=1,nvar
          call fem_file_skip_data(iformat,iunit
     +                          ,nvers,np,lmax,string,ierr)
          if( ierr .ne. 0 ) goto 97
          strings(i) = string
        end do

	call get_vars_from_string(nvar,strings,ivars)
	call choose_var(nvar,ivars,strings,varline,ivarplot,bvect)
	call get_vars_to_plot(nvar,ivars,ivar3,ivarplot,bvect,ivnum,ivs)
	call set_ivel(ivar3,ivel)
	bvel = ivel > 0
	if( bvect .and. .not. bvel ) stop 'error stop: ivel == 0'
        if( ivs(1) == 0 ) then
          write(6,*) 'no such variable in file: ',ivar3,ivarplot,ivs
	  stop 'error stop plot_fem_file'
	else if( ivs(2) == 0 .and. bvect ) then
          write(6,*) 'need two variables for directional plot ',ivs
	  stop 'error stop plot_fem_file'
	end if
	write(6,*) 'what to plot: ',ivar3,ivarplot,ivs

	if( .not. breg .and. .not. bhasbasin ) goto 94

	!--------------------------------------------------------------
	! close and re-open file
	!--------------------------------------------------------------

        close(iunit)

        np = 0
        call fem_file_read_open(infile,np,iformat,iunit)
        if( iunit .le. 0 ) stop

        !--------------------------------------------------------------
        ! set time
        !--------------------------------------------------------------

        call dts_convert_to_atime(datetime,dtime,atime)
        atime0 = atime          !absolute time of first record

	date = datetime(1)
	time = datetime(2)
        call ptime_init
        call ptime_set_date_time(date,time)
        call elabtime_date_and_time(date,time)
        call elabtime_minmax(stmin,stmax)
        call elabtime_set_inclusive(.false.)

	!--------------------------------------------------------------
	! initialize plot
	!--------------------------------------------------------------

        call init_plot

	bsect = getisec() /= 0
	call setlev(layer)
	b2d = layer == 0

        call initialize_color

	call qopen

        !--------------------------------------------------------------
        ! loop on records
        !--------------------------------------------------------------

        irec = 0
	nplot = 0

        do
          call fem_file_read_params(iformat,iunit,dtime
     +                          ,nvers,np,lmax,nvar,ntype,datetime,ierr)

          if( ierr .lt. 0 ) exit
          irec = irec + 1
          if( ierr .gt. 0 ) goto 99
          if( nvar .ne. nvar0 ) goto 96
          if( lmax .ne. lmax0 ) goto 96
          if( np .ne. np0 ) goto 96

          call dts_convert_to_atime(datetime,dtime,atime)
          call dts_format_abs_time(atime,line)
	  call ptime_set_atime(atime)

          if( bdebug ) write(6,*) irec,atime,trim(line)

          call fem_file_read_2header(iformat,iunit,ntype,lmax
     +                  ,hlv,regpar,ierr)
          if( ierr .ne. 0 ) goto 98
	  call init_sigma_info(lmax,hlv)

	  bskip = .false.
	  if( elabtime_over_time(atime) ) exit
	  if( .not. elabtime_in_time(atime) ) bskip = .true.
	  if( ifreq > 0 .and. mod(irec,ifreq) /= 0 ) bskip = .true.

	  write(6,*) irec,atime,trim(line)

          do i=1,nvar
            if( bskip ) then
              call fem_file_skip_data(iformat,iunit
     +                          ,nvers,np,lmax,string,ierr)
            else
              call fem_file_read_data(iformat,iunit
     +                          ,nvers,np,lmax
     +                          ,string
     +                          ,il,hd
     +                          ,lmax,data(1,1,i)
     +                          ,ierr)
            end if
            if( ierr .ne. 0 ) goto 97
            if( string .ne. strings(i) ) goto 95
          end do

	  if( bskip ) cycle

	  nplot = nplot + 1
	  data3d = data(:,:,ivs(1))
	  if( bvect ) data3ddir = data(:,:,ivs(2))

	  flag = dflag
	  call adjust_levels_with_flag(nlvdi,np,il,flag,data3d)

	  if( lmax == 1 ) then
	    data2d(:) = data3d(1,:)
	    if( bvect ) data2ddir(:) = data3ddir(1,:)
	  else if( b2d ) then
	    call fem_average_vertical(nlvdi,np,lmax,il,hlv,hd
     +					,data3d,data2d)
	    if( bvect ) then
	      call fem_average_vertical(nlvdi,np,lmax,il,hlv,hd
     +					,data3ddir,data2ddir)
	    end if
	  else
	    call extlev(layer,nlvdi,np,il,data3d,data2d)
	    if( bvect ) then
	      call extlev(layer,nlvdi,np,il,data3ddir,data2ddir)
	    end if
	  end if

	  !write(6,*) 'data: ',lmax,b2d,layer,(data2d(i),i=1,np,np/10)

	  if( bplotreg ) then
	    call ploreg(np,data2d,regpar,varline,bintp,.true.)
	  else
            !call outfile_make_hkv(nkn,nel,nen3v,hm3v,hev,hkv)
	    if( np /= nkn ) stop 'cannot handle yet: internal error (9)'
	    ilhkv = il
            call ilhk2e(nkn,nel,nen3v,ilhkv,ilhv)
            call adjust_layer_index(nel,nlv,hev,hlv,ilhv)
	    call make_mask(layer)
	    if( bvect ) then
	      call plo2vel(ivel,'3D ')
	    else
	      call ploval(np,data2d,varline)
	    end if
	  end if

	end do

        !--------------------------------------------------------------
        ! end of routine
        !--------------------------------------------------------------

	return
   92   continue
        write(6,*) 'for non regular file we need bas, not grd file'
        stop 'error stop plot_fem_file: basin'
   93   continue
        write(6,*) 'incompatible node numbers: ',nkn,np
        stop 'error stop plot_fem_file: basin'
   94   continue
        write(6,*) 'fem file with non regular data needs basin'
        write(6,*) 'please specify basin on command line'
        stop 'error stop plot_fem_file: basin'
   95   continue
        write(6,*) 'strings not in same sequence: ',i
        write(6,*) string
        write(6,*) strings(i)
        stop 'error stop plot_fem_file: strings'
   96   continue
        write(6,*) 'nvar,nvar0: ',nvar,nvar0
        write(6,*) 'lmax,lmax0: ',lmax,lmax0    !this might be relaxed
        write(6,*) 'np,np0:     ',np,np0        !this might be relaxed
        write(6,*) 'cannot change number of variables'
        stop 'error stop plot_fem_file'
   97   continue
        write(6,*) 'record: ',irec
        write(6,*) 'cannot read data record of file'
        stop 'error stop plot_fem_file'
   98   continue
        write(6,*) 'record: ',irec
        write(6,*) 'cannot read second header of file'
        stop 'error stop plot_fem_file'
   99   continue
        write(6,*) 'record: ',irec
        write(6,*) 'cannot read header of file'
        stop 'error stop plot_fem_file'
	end

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

	subroutine make_mask(layer)

	use levels
        use mod_hydro_plot
        !use mod_hydro

	implicit none

	integer layer

	integer level

	level = layer
	if( level == 0 ) level = 1

        call reset_dry_mask

        call set_level_mask(bwater,ilhv,level)        !element has this level
        call make_dry_node_mask(bwater,bkwater)       !copy elem to node mask

        call adjust_no_plot_area
        call make_dry_node_mask(bwater,bkwater)	      !copy elem to node mask
        call info_dry_mask(bwater,bkwater)

	end

!***************************************************************
!***************************************************************
!***************************************************************

        subroutine initialize_color

        implicit none

        integer icolor
        real getpar
	logical has_color_table

        call colsetup
	call admin_color_table

	if( has_color_table() ) call putpar('icolor',8.)

        icolor = nint(getpar('icolor'))
        call set_color_table( icolor )
        call set_default_color_table( icolor )

	call write_color_table

        end

!***************************************************************

c*****************************************************************

        subroutine allocate_2d_arrays

        use mod_hydro_plot
        use mod_geom
        use mod_depth
        use evgeom
        use basin, only : nkn,nel,ngr,mbw

        implicit none

        call ev_init(nel)
        call mod_geom_init(nkn,nel,ngr)

        call mod_depth_init(nkn,nel)

        write(6,*) 'allocate_2d_arrays: ',nkn,nel,ngr

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
	  call basin_set_read_basin(.true.)
          !breadbas = .false.
        else
          write(6,*) 'Cannot read this file: ',trim(file)
          stop 'error stop read_given_file: format not recognized'
        end if

        end

c*****************************************************************
c*****************************************************************
c*****************************************************************
c routines for directional plots
c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine set_ivel(ivar3,ivel)

! ivel == 3	wind
! ivel == 4	wave
! ivel == 5	velocities already prepared

	implicit none

	integer ivar3
	integer ivel

	if( ivar3 == 2 .or. ivar3 == 3 ) then
	  ivel = 5
	else if( ivar3 > 230 .and. ivar3 < 240 ) then
	  ivel = 4
	else if( ivar3 == 21 ) then
	  ivel = 3
	else
	  ivel = 0
	end if

	end

c*****************************************************************

	subroutine directional_init(nvar,ivars,ivar3,bvect,ivarplot)

! from ivar3 sets up bvect and ivarplot

	implicit none

	integer nvar			!total number of variables in ivars
	integer ivars(nvar)		!id of variables in file
	integer ivar3			!id of what to plot
	logical bvect			!it is a vector variable (out)
	integer ivarplot(2)		!what variable id to use (out)

	logical bvel,bwave,bwind
	integer ivar,iv,nv,ivel

	ivarplot = ivar3

	bvel  = ivar3 == 2 .or. ivar3 == 3		!velocity/transport
	bwave = ivar3 > 230 .and. ivar3 < 240		!wave plot
	bwind = ivar3 == 21				!wind

	bvect = .false.
	bvect = bvect .or. bvel
	bvect = bvect .or. bwave
	bvect = bvect .or. bwind

	if( .not. bvect ) return

	if( bvel ) then
	  ivarplot = 2
	  if( any( ivars == 3 ) ) ivarplot = 3 	!transports in shy file
	end if
	if( bwave ) ivarplot = (/ivar3,233/)
	if( bwind ) ivarplot = 21

	nv = 0
	do iv=1,nvar
	  ivar = ivars(iv)
	  if( any( ivarplot == ivar ) ) then
	    nv = nv + 1
	  end if
	end do

	if( nv == 2 ) then
	  write(6,*) 'can plot vector: ',ivarplot
	else if( nv == 1 ) then
	  if( bvect ) then
	    bvect = .false.
	    write(6,*) 'can plot vector only as scalar: ',ivarplot
	  end if
	else if( nv == 0 ) then
	  write(6,*) '*** file does not contain needed varid: ',ivar3
          !stop 'error stop shyplot'
	else
	  write(6,*) '*** error in variables: ',ivar3,nv,nvar
	  write(6,*) ivars(:)
          stop 'error stop shyplot'
	end if

	end

c*****************************************************************

	subroutine directional_insert(bvect,ivar,ivar3,ivarplot
     +					,n,cv2,ivel)

	use basin
	use mod_hydro_plot

	implicit none

	logical bvect		!can we plot a vector variable?
	integer ivar		!variable id read
	integer ivar3		!variable id requested
	integer ivarplot(2)	!variable ids needed
	integer n		!total number of data points
	real cv2(n)		!data (2d)
	integer ivel		!on return indicates if and what to plot

! ivel == 1	velocities (from transports)
! ivel == 2	transports
! ivel == 3	wind
! ivel == 4	wave
! ivel == 5	velocities already prepared
!
! ivar == 2	velocities
! ivar == 3	transports

	logical bwave,bvel,bwind
	integer, save :: iarrow = 0

	ivel = 0
	bonelem = .false.
	bisreg = .false.
	bistrans = .false.

	if( .not. bvect ) return

	bvel  = ivar3 == 2 .or. ivar3 == 3		!velocity/transport
	bwave = ivar3 > 230 .and. ivar3 < 240		!wave plot
	bwind = ivar3 == 21				!wind

	if( bvel ) then
	  if( ivar == 3 ) then				!we read transports
	    if( n /= nel ) goto 99			!only for shy files
	    iarrow = iarrow + 1
	    bonelem = .true.
	    bistrans = .true.
	    if( iarrow == 1 ) utrans(1:nel) = cv2(1:nel)
	    if( iarrow == 2 ) vtrans(1:nel) = cv2(1:nel)
	  end if
	  if( ivar == 2 ) then				!we read velocities
	    iarrow = iarrow + 1
	    if( n == nel ) then
	      bonelem = .true.
	      if( iarrow == 1 ) uvelem(1:n) = cv2(1:n)
	      if( iarrow == 2 ) vvelem(1:n) = cv2(1:n)
	    else
	      bisreg = ( n /= nkn )
	      if( iarrow == 1 ) uvnode(1:n) = cv2(1:n)
	      if( iarrow == 2 ) vvnode(1:n) = cv2(1:n)
	    end if
	  end if
	  if( iarrow == 2 ) then
	    ivel = ivar3 - 1
	    !call make_vertical_velocity
	    wsnv = 0.		!FIXME
	  end if
	end if

	if( bwave ) then
	  !if( n /= nkn ) goto 97
	  bisreg = ( n /= nkn )
	  if( ivarplot(1) == ivar ) then
	    iarrow = iarrow + 1
	    uvspeed(1:n) = cv2(1:n)
	  else if( ivarplot(2) == ivar ) then
	    iarrow = iarrow + 1
	    uvdir(1:n) = cv2(1:n)
	  end if
	  if( iarrow == 2 ) then
	    ivel = 4
	    call polar2xy(n,uvspeed,uvdir,uvnode,vvnode)
	  end if
	end if

	if( bwind ) then
	  !if( n /= nkn ) goto 96
	  bisreg = ( n /= nkn )
	  if( ivar == 21 ) then
	    iarrow = iarrow + 1
	    if( iarrow == 1 ) uvnode(1:n) = cv2(1:n)
	    if( iarrow == 2 ) vvnode(1:n) = cv2(1:n)
	  end if
	  if( iarrow == 2 ) then
	    ivel = 3
	    uvspeed = sqrt( uvnode**2 + vvnode**2 )
	  end if
	end if

	if( ivel > 0 ) iarrow = 0	!reset for next records

	return
   96	continue
	write(6,*) 'can read wind data only on nodes: ',n,nkn
	stop 'error stop directional_insert: wind data not on nodes'
   97	continue
	write(6,*) 'can read wave data only on nodes: ',n,nkn
	stop 'error stop directional_insert: wave data not on nodes'
   98	continue
	write(6,*) 'can read velocities only on nodes: ',n,nkn
	stop 'error stop directional_insert: velocities not on nodes'
   99	continue
	write(6,*) 'can read transports only on elements: ',n,nel
	stop 'error stop directional_insert: transports not on elems'
	end

c*****************************************************************

	subroutine choose_var(nvar,ivars,strings,varline,ivarplot,bvect)

	use plotutil
	use levels

! choses variable to be plotted
!
! either ivnum (sequential number) or ivar3 must have been set
! both of them are set on return
!
! cwarguments set on return:  varline, ivarplot
! global values set: ivar3, ivnum

	implicit none

	integer nvar
	integer ivars(nvar)
	character*80 strings(nvar)
	character*80 varline
	integer ivarplot(2)

	logical bvect
	integer nv,iv,ivar,isub

!	---------------------------------------------------
!	write contents of file to terminal
!	---------------------------------------------------

	call shy_print_descriptions(nvar,ivars,strings)

!	---------------------------------------------------
!	if sequential number is given, use this one
!	---------------------------------------------------

	if( ivnum > 0 ) then
	  if( ivnum > nvar ) then
	    write(6,*) 'ivnum too big for nvar'
	    write(6,*) 'ivnum,nvar: ',ivnum,nvar
            stop 'error stop shyplot'
	  end if
	  ivar3 = ivars(ivnum)
	end if

!	---------------------------------------------------
!	only one variable in file - plot this one
!	---------------------------------------------------

	if( nvar == 1 .and. ivar3 == 0 ) then
	  ivnum = 1
	  ivar3 = ivars(ivnum)
	end if
        if( ivar3 == 0 ) then
          write(6,*) 'no variable given to be plotted: ',ivar3
          stop 'error stop shyplot'
        end if

!	---------------------------------------------------
!	prepare for directional plot and read STR file
!	---------------------------------------------------

	call directional_init(nvar,ivars,ivar3,bvect,ivarplot)
	call read_str_files(ivar3)

!	---------------------------------------------------
!	see if ivar3 is in file
!	---------------------------------------------------

	write(6,*) 
	write(6,*) 'varid to be plotted:       ',ivar3
	nv = 0
	do iv=1,nvar
	  ivar = ivars(iv)
	  !write(6,'(2i10,4x,a)') iv,ivar,trim(strings(iv))
	  !if( ivar == ivar3 ) nv = nv + 1
	  if( ivar3 == ivar .or. any( ivarplot == ivar ) ) then
	    nv = nv + 1
	    if( ivnum == 0 ) ivnum = iv
	  end if
	end do

!	---------------------------------------------------
!	error check
!	---------------------------------------------------

	!if( nv == 0 .and. .not. bvect ) then
	if( nv == 0 ) then
	  call ivar2string(ivar3,varline,isub)
          write(6,*) 'no such variable in file: ',ivar3,varline
          stop 'error stop shyplot'
        end if

	if( layer > nlv ) then
          write(6,*) 'no such layer: ',layer
          write(6,*) 'maximum layer available: ',nlv
          stop 'error stop shyplot'
	end if

!	---------------------------------------------------
!	make varline and write to terminal
!	---------------------------------------------------

	call mkvarline(ivar3,varline)
	write(6,*) 
	write(6,*) 'information for plotting:'
	write(6,*) 'varline: ',trim(varline)
	write(6,*) 'ivnum: ',ivnum
	write(6,*) 'ivar3: ',ivar3
	write(6,*) 'layer: ',layer
	write(6,*) 

!	---------------------------------------------------
!	end of routine
!	---------------------------------------------------

	end

c*****************************************************************

	subroutine fem_average_vertical(nlvddi,np,lmax,ilhkv,hlv,hd
     +					,data,data2d)

	implicit none

	integer nlvddi
	integer np,lmax
	integer ilhkv(np)
	real hlv(nlvddi)
	real hd(np)
	real data(nlvddi,np)
	real data2d(np)

	real hl(nlvddi)
	integer nsigma,k,lm,l,nlv
	real z,hsigma,h,hh,zeps
	double precision vacu,dacu

        call get_sigma_info(nlv,nsigma,hsigma)
	z = 0.
	zeps=0.01

	do k=1,np

	  lm = min(lmax,ilhkv(k))
	  if( lm <= 1 ) then
	    data2d(k) = data(1,k)
	    cycle
	  end if
	  h = hd(k)
	  if( h < -990. ) h = hlv(lm)
	  !if( h == -1. ) h = 1.
	  if( h+z<zeps ) z = zeps-h
          call get_layer_thickness(lm,nsigma,hsigma
     +                          ,z,h,hlv,hl)

	  vacu = 0.
	  dacu = 0.
	  do l=1,lm
	    hh = hl(l)
	    vacu = vacu + data(l,k)*hh
	    dacu = dacu + hh
	  end do
	  data2d(k) = vacu / dacu

	end do

	end

c*****************************************************************

	subroutine adjust_levels_with_flag(nlvddi,np,ilhkv,flag,data3d)

	implicit none

	integer nlvddi,np
	integer ilhkv(np)
	real flag
	real data3d(nlvddi,np)

	integer i

	do i=1,np
	  if( ilhkv(i) <= 0 ) then
	    ilhkv(i) = 1
	    data3d(1,i) = flag
	  end if
	end do

	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine get_var_to_plot(nvar,ivars,ivar,iv)

	implicit none

	integer nvar
	integer ivars(nvar)
	integer ivar		!desired variable
	integer iv		!return

	do iv=1,nvar
	  if( ivars(iv) == ivar ) return
	end do

	iv = 0

	end

c*****************************************************************

	subroutine get_vars_to_plot(nvar,ivars,ivar,ivarplot
     +					,bvect,ivnum,ivs)

	implicit none

	integer nvar
	integer ivars(nvar)
	integer ivar		!desired variable
	integer ivarplot(2)	!desired variable for vector plot
	logical bvect
	integer ivnum
	integer ivs(2)		!return

	integer iv,ia

	ia = 1
	ivs = 0

	if( bvect ) then
	  do iv=1,nvar
	    if( ivars(iv) == ivarplot(ia) ) then
	      if( ia <= 2 ) ivs(ia) = iv
	      ia = ia + 1
	    end if
	  end do
	  if( ia /= 3 ) then
	    write(6,*) 'for vector plot need 2 variables to read: ',ia-1
	    ivs = 0
	  end if
	else
	  !do iv=1,nvar
	  !  if( ivars(iv) == ivar ) then
	  !    ivs(1) = iv
	  !    return
	  !  end if
	  !end do
	  ivs(1) = ivnum
	end if

	end

c*****************************************************************

