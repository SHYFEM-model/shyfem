!
! $Id: nos2nc.f,v 1.15 2009-11-18 16:50:37 georg Exp $
!
! convert NOS to NC file
!
! revision log :
!
! 02.09.2003	ggu	adapted to new OUS format
! 24.01.2005	ggu	computes maximum velocities for 3D (only first level)
! 16.10.2007	ggu	new debug routine
! 27.10.2009    ggu     include evmain.h, compute volume
! 23.03.2011    ggu     compute double precision u/v-min/max of first level
! 21.01.2013    ggu     restructured
! 25.01.2013    ggu     regular and fem outout in one routine
! 20.02.2013    ggu     choose period implemented
! 15.04.2016    ggu     ndim eliminated (allocatable)
! 27.04.2016    mbj     adapted to new framework
!
!***************************************************************

	program nos2nc

! reads nos file and writes NetCDF file

	use depth
	use hydro_admin
	!use evgeom_2nc
	use levels
	use basin
        use version
        use apnfile
        use depth_util_2nc
        use netcdf_util
        use dts,        only:   dtsini
        use nosutil
        use ionos_2nc
        use regular_2nc
        use output_util
        use netcdf


	implicit none

!-------------------------------------------------

	character*80 title

	integer, allocatable :: ivars(:)
	integer, allocatable :: var_ids(:)

	double precision, allocatable :: uprv(:,:)
	double precision, allocatable :: vprv(:,:)
	double precision, allocatable :: ut2v(:)
	double precision, allocatable :: vt2v(:)
	double precision, allocatable :: u2v(:)
	double precision, allocatable :: v2v(:)
	double precision, allocatable :: cv3(:,:)

        integer nvers,nin,lmax,l
        integer itanf,itend,idt
	integer it,ie,i
        integer ierr,nread,ndry
	integer irec,maxrec,iwrite
        integer nknnos,nelnos
	integer iztype
	integer nvar,ivar
        double precision href,hzoff,hlvmin
	double precision volume
	double precision zmin,zmax
	double precision umin,umax
	double precision vmin,vmax
#ifdef SINGLEP        
	real flag
	real, allocatable :: xlon(:)
	real, allocatable :: ylat(:)
	real, allocatable :: fm(:,:,:)
	real, allocatable :: ndepth(:,:)
	real, allocatable :: value3d(:,:,:)
	real, allocatable :: vnc3d(:,:,:)
	real, allocatable :: var3d(:)
	real x0,y0,dx,dy
#else
	double precision flag
	double precision, allocatable :: xlon(:)
	double precision, allocatable :: ylat(:)
	double precision, allocatable :: fm(:,:,:)
	double precision, allocatable :: ndepth(:,:)
	double precision, allocatable :: value3d(:,:,:)
	double precision, allocatable :: vnc3d(:,:,:)
	double precision, allocatable :: var3d(:)
	double precision x0,y0,dx,dy
#endif

	integer nx,ny,nxymax
	double precision, allocatable :: value2d(:,:)

	logical breg
	logical bdate
        integer ncid
        integer var_id
	integer date00,time00
	integer date,time
	integer it0
	integer iperiod,its,ite,nfreq
	logical bwrite
	logical bfirst

	character*80 units,std

	call shyfem_copyright('nos2nc - netcdf output')

!-----------------------------------------------------------------
! initialize basin and simulation
!-----------------------------------------------------------------

	maxrec = 2		!max number of records to be written
	maxrec = 0		!max number of records to be written (0 -> all)
	nxymax = 0		!max size of regular grid (0 -> any)
	!nxymax = 400		!max size of regular grid (0 -> any)

	it0 = 0			!subtract from it (hack)

!-----------------------------------------------------------------
! do not change anything below here
!-----------------------------------------------------------------

        time00 = 0
        date00 = 2000
        bdate = .true.          !if true read from terminal

#ifdef SINGLEP        
        flag = -999.
#else
        flag = -999.d0
#endif

	nread = 0
	irec = 0
	iwrite = 0

	call ap_init(.false.,3,0,0)

!-----------------------------------------------------------------
! Init modules
!-----------------------------------------------------------------

	call ev_init(nel)
	call mod_depth_init(nkn,nel)

	call set_ev

	call makehkv_minmax(hkv,1)
	call makehev(hev)

!-----------------------------------------------------------------
! get input from terminal about what to do
!-----------------------------------------------------------------

	if( bdate ) call read_date_and_time(date00,time00)
	call dtsini(date00,time00)

        call get_dimensions(nx,ny,x0,y0,dx,dy)
        breg = nx .gt. 0 .and. ny .gt. 0          !regular output
        if( breg ) then
	  write(6,*) 'NETCDF output: regular ',dx,dy,nx,ny
	  if( nxymax > 0 .and. max(nx,ny) > nxymax ) goto 95
	  allocate(xlon(nx),ylat(ny))
	  allocate(ndepth(nx,ny),value2d(nx,ny))
	  allocate(fm(4,nx,ny))
	  call set_reg_xy(nx,ny,x0,y0,dx,dy,xlon,ylat)
          call setgeo(x0,y0,dx,dy,flag)
          call av2fm(fm,nx,ny)
	else
	  write(6,*) 'NETCDF output: unstructured (fem)'
        end if

	call get_period(iperiod,its,ite,nfreq)

!-----------------------------------------------------------------
! first read of nos file to get the dimensions
!-----------------------------------------------------------------

	call open_nos_type('.nos','old',nin)

	call nos_is_nos_file(nin,nvers)
	if( nvers .le. 0 ) then
          write(6,*) 'nvers: ',nvers
          stop 'error stop noselab: not a valid nos file'
        end if

	call peek_nos_header(nin,nknnos,nelnos,nlv,nvar)

        if( nkn /= nknnos .or. nel /= nelnos ) goto 94

!-----------------------------------------------------------------
! allocate arrays
!-----------------------------------------------------------------

	call levels_init(nkn,nel,nlv)

	allocate(uprv(nlv,nkn),vprv(nlv,nkn))
	allocate(ut2v(nel),vt2v(nel),u2v(nel),v2v(nel))
	allocate(var3d(nlv*nkn))
	allocate(cv3(nlv,nkn))

	if( breg ) allocate(value3d(nlv,nx,ny),vnc3d(nx,ny,nlv))

!-----------------------------------------------------------------
! read header of simulation
!-----------------------------------------------------------------

        call read_nos_header(nin,nkn,nel,nlv,ilhkv,hlv,hev)
        call nos_get_params(nin,nknnos,nelnos,nlv,nvar)
	call nos_get_title(nin,title)
	call nos_get_date(nin,date,time)
	if( date .gt. 0 ) then
	  date00 = date
	  time00 = time
	  call dtsini(date00,time00)
	end if

	call init_sigma_info(nlv,hlv)
	call level_k2e(nkn,nel,nen3v,ilhkv,ilhv)
	call compute_iztype(iztype)

	allocate(ivars(nvar),var_ids(nvar))
	call nos_get_vars(nin,nvar,ivars)

        write(6,*) 'Available variables: ',nvar
        write(6,*) (ivars(i),i=1,nvar)

	if( breg ) call get_lmax_reg(nx,ny,fm,ilhv,lmax)

!-----------------------------------------------------------------
! prepare netcdf file
!-----------------------------------------------------------------

        if( breg ) then
          call nc_open_reg(ncid,nx,ny,lmax,flag,date00,time00,iztype)
        else
          call nc_open(ncid,nkn,nel,nlv,date00,time00,iztype)
        end if
        call nc_global(ncid,title)

	do i=1,nvar
	  call nc_init_variable(ncid,breg,3,ivars(i),flag,var_ids(i))
	end do

        call nc_end_define(ncid)
        if( breg ) then
          call fm2am2d(hkv,nx,ny,fm,ndepth)
          call nc_write_coords_reg(ncid,nx,ny,xlon,ylat,ndepth)
        else
          call nc_write_coords(ncid)
        end if

!-----------------------------------------------------------------
! loop on data of simulation
!-----------------------------------------------------------------

  300   continue

	call rdnos(nin,it,ivar,nlv,ilhkv,cv3,ierr)

	it = it - it0

        if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
        if(ierr.ne.0) goto 100

	nread=nread+1

	i = mod(nread,nvar)
	if( i .eq. 0 ) i = nvar
	bfirst = i .eq. 1

	if( bfirst ) then	!first variable
          irec = irec + 1
          call write_time(it)
	  call check_period(it,iperiod,its,ite,nfreq,bwrite)
	end if

	if( .not. bwrite ) goto 300

	if( ivars(i) .ne. ivar ) goto 92
	var_id = var_ids(i)

        if( bfirst ) then
	  iwrite = iwrite + 1
	  if ( maxrec .gt. 0 .and. iwrite .gt. maxrec ) goto 100
	  call nc_write_time(ncid,iwrite,it)
	end if

	write(6,*) '   writing: ',ivar,i,irec,nread,iwrite

	if( breg ) then
	  call fm2am3d(nlv,ilhv,cv3,lmax,nx,ny,fm,value3d)
	  call nc_rewrite_3d_reg(lmax,nx,ny,value3d,vnc3d)
	  call nc_write_data_3d_reg(ncid,var_id,iwrite,lmax,nx,ny,vnc3d)
	else
	  call nc_compact_3d(nlv,nlv,nkn,cv3,var3d)
          call nc_write_data_3d(ncid,var_id,iwrite,nlv,nkn,var3d)
	end if

	goto 300

  100	continue

!-----------------------------------------------------------------
! end of loop
!-----------------------------------------------------------------

	write(6,*)
	write(6,*) nread,' records read'
	write(6,*)

	if( breg ) call write_dimensions(nx,ny,x0,y0,dx,dy)

        call nc_close(ncid)

!-----------------------------------------------------------------
! end of routine
!-----------------------------------------------------------------

        stop
   92   continue
        write(6,*) 'wrong order of variables'
        write(6,*) 'read: ',ivar,'   expected: ',ivars(i)
        stop 'error stop nos2nc: variables'
   94   continue
        write(6,*) 'incompatible basin and simulation'
        write(6,*) 'nkn,nknnos: ',nkn,nknnos
        write(6,*) 'nel,nelnos: ',nel,nelnos
        stop 'error stop nos2nc: parameter mismatch'
   95   continue
        write(6,*) 'regular grid too big'
        write(6,*) 'nx,ny: ',nx,ny,'   nxymax: ',nxymax
        write(6,*) 'please increase nxymax or set to zero for any size'
        stop 'error stop nos2nc: size regular grid'
        end

!******************************************************************

