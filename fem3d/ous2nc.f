c
c $Id: ous2nc.f,v 1.15 2009-11-18 16:50:37 georg Exp $
c
c convert OUS to NC file
c
c revision log :
c
c 02.09.2003	ggu	adapted to new OUS format
c 24.01.2005	ggu	computes maximum velocities for 3D (only first level)
c 16.10.2007	ggu	new debug routine
c 27.10.2009    ggu     include evmain.h, compute volume
c 23.03.2011    ggu     compute real u/v-min/max of first level
c 21.01.2013    ggu     restructured
c 25.01.2013    ggu     regular and fem outout in one routine
c 20.02.2013    ggu     choose period implemented
c 27.04.2016    mbj     adapted to new framework
c
c***************************************************************

	program ous2nc

c reads ous file and writes NetCDF file

	use mod_depth
	use mod_hydro
	use evgeom
	use levels
	use basin

	implicit none

c-------------------------------------------------

	character*80 title

	real, allocatable :: uprv(:,:)
	real, allocatable :: vprv(:,:)
	real, allocatable :: ut2v(:)
	real, allocatable :: vt2v(:)
	real, allocatable :: u2v(:)
	real, allocatable :: v2v(:)

	real, allocatable :: var3d(:)

        integer nvers,nin,lmax,l
        integer itanf,itend,idt
	integer it,ie,i
        integer ierr,nread,ndry
	integer irec,maxrec,iwrite
        integer nknous,nelous
	integer iztype
        real href,hzoff,hlvmin
	real volume
	real zmin,zmax
	real umin,umax
	real vmin,vmax
	real flag

	integer nx,ny,nxymax
	real, allocatable :: xlon(:)
	real, allocatable :: ylat(:)
	real, allocatable :: depth(:,:)
	real, allocatable :: value2d(:,:)
	real, allocatable :: fm(:,:,:)
	real, allocatable :: value3d(:,:,:)
	real, allocatable :: vnc3d(:,:,:)
	real x0,y0,dx,dy

	logical breg
	logical bdate
        integer ncid
        integer dimids_2d(2)
        integer coord_varid(3)
        integer rec_varid
        integer z_id,u_id,v_id
	integer date0,time0
	integer date,time
	integer it0
	integer iperiod,its,ite,nfreq
	logical bwrite

	character*80 units,std

	integer iapini,ideffi

	call shyfem_copyright('ous2nc - netcdf output')

c-----------------------------------------------------------------
c initialize basin and simulation
c-----------------------------------------------------------------

	maxrec = 2		!max number of records to be written
	maxrec = 0		!max number of records to be written (0 -> all)
	nxymax = 0		!max size of regular grid (0 -> any)
	nxymax = 400		!max size of regular grid (0 -> any)

	it0 = 0			!subtract from it (hack)

c-----------------------------------------------------------------
c do not change anything below here
c-----------------------------------------------------------------

        time0 = 0
        date0 = 2000
        bdate = .true.          !if true read from terminal
        flag = -999.

	nread = 0
	irec = 0
	iwrite = 0

	call ap_init(.false.,3,0,0)

c-----------------------------------------------------------------
c Init modules
c-----------------------------------------------------------------

	call ev_init(nel)
	call mod_depth_init(nkn,nel)

	call set_ev

	call makehkv_minmax(hkv,1)
	call makehev(hev)

c-----------------------------------------------------------------
c get input from terminal about what to do
c-----------------------------------------------------------------

	if( bdate ) call read_date_and_time(date0,time0)
	call dtsini(date0,time0)

	call get_dimensions(nx,ny,x0,y0,dx,dy)
        breg = nx .gt. 0 .and. ny .gt. 0          !regular output
        if( breg ) then
          write(6,*) 'NETCDF output: regular ',dx,dy,nx,ny
	  if( nxymax > 0 .and. max(nx,ny) > nxymax ) goto 95
	  allocate(xlon(nx),ylat(ny))
	  allocate(depth(nx,ny),value2d(nx,ny))
	  allocate(fm(4,nx,ny))
	  call set_reg_xy(nx,ny,x0,y0,dx,dy,xlon,ylat)
          call setgeo(x0,y0,dx,dy,flag)
          call av2fm(fm,nx,ny)
        else
          write(6,*) 'NETCDF output: unstructured (fem)'
        end if

	call get_period(iperiod,its,ite,nfreq)

c-----------------------------------------------------------------
c first read of ous file to get the dimensions
c-----------------------------------------------------------------

	call open_ous_type('.ous','old',nin)

	call ous_is_ous_file(nin,nvers)
	if( nvers .le. 0 ) then
          write(6,*) 'nvers: ',nvers
          stop 'error stop ouselab: not a valid ous file'
        end if

	call peek_ous_header(nin,nknous,nelous,nlv)

        if( nkn /= nknous .or. nel /= nelous ) goto 94

c-----------------------------------------------------------------
c allocate arrays
c-----------------------------------------------------------------

	call levels_init(nkn,nel,nlv)
	call mod_hydro_init(nkn,nel,nlv)

	allocate(uprv(nlv,nkn),vprv(nlv,nkn))
	allocate(ut2v(nel),vt2v(nel),u2v(nel),v2v(nel))
	allocate(var3d(nlv*nkn))

	if( breg ) allocate(value3d(nlv,nx,ny),vnc3d(nx,ny,nlv))

c-----------------------------------------------------------------
c read header of simulation
c-----------------------------------------------------------------

        call read_ous_header(nin,nkn,nel,nlv,ilhv,hlv,hev)
        call ous_get_params(nin,nknous,nelous,nlv)
	call ous_get_title(nin,title)
	call ous_get_date(nin,date,time)
	if( date .gt. 0 ) then
	  date0 = date
	  time0 = time
	  call dtsini(date0,time0)
	end if

	call init_sigma_info(nlv,hlv)
	call level_e2k(nkn,nel,nen3v,ilhv,ilhkv)
	call compute_iztype(iztype)

	if( breg ) call get_lmax_reg(nx,ny,fm,ilhv,lmax)

c-----------------------------------------------------------------
c prepare netcdf file
c-----------------------------------------------------------------

	if( breg ) then
	  call nc_open_reg(ncid,nx,ny,lmax,flag,date0,time0,iztype)
	else
          call nc_open(ncid,nkn,nel,nlv,date0,time0,iztype)
	end if
	call nc_global(ncid,title)

	call nc_init_variable(ncid,breg,2,1,flag,z_id)
	call nc_init_variable(ncid,breg,3,2,flag,u_id)
	call nc_init_variable(ncid,breg,3,3,flag,v_id)

        call nc_end_define(ncid)
	if( breg ) then
	  call fm2am2d(hkv,nx,ny,fm,depth)
	  call nc_write_coords_reg(ncid,nx,ny,xlon,ylat,depth)
	else
          call nc_write_coords(ncid)
	end if

c-----------------------------------------------------------------
c loop on data of simulation
c-----------------------------------------------------------------

  300   continue

	call ous_read_record(nin,it,nlv,ilhv,znv,zenv,utlnv,vtlnv,ierr)

	it = it - it0

        if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
        if(ierr.ne.0) goto 100

	nread=nread+1

        irec = irec + 1
        call write_time(it)
	call check_period(it,iperiod,its,ite,nfreq,bwrite)

	if( .not. bwrite ) goto 300

	call mima(znv,nknous,zmin,zmax)
        call comp_barotropic(nel,nlv,ilhv,utlnv,vtlnv,ut2v,vt2v)
	call comp_vel2d(nel,hev,zenv,ut2v,vt2v,u2v,v2v
     +				,umin,vmin,umax,vmax)
	call compute_volume(nel,zenv,hev,volume)

c        call debug_write_node(0,it,nread,nkn,nel,nlv,nkn,nel,nlv
c     +          ,nen3v,zenv,znv,utlnv,vtlnv)

        call transp2vel(nel,nkn,nlv,nlv,hev,zenv,nen3v
     +                          ,ilhv,hlv,utlnv,vtlnv
     +                          ,uprv,vprv)

	iwrite = iwrite + 1
	if ( maxrec .gt. 0 .and. iwrite .gt. maxrec ) goto 100
        call nc_write_time(ncid,iwrite,it)

	write(6,*) '   writing: ',irec,nread,iwrite

	if( breg ) then
	  call fm2am2d(znv,nx,ny,fm,value2d)
          call nc_write_data_2d_reg(ncid,z_id,iwrite,nx,ny,value2d)

	  call fm2am3d(nlv,ilhv,uprv,lmax,nx,ny,fm,value3d)
	  call nc_rewrite_3d_reg(lmax,nx,ny,value3d,vnc3d)
	  call nc_write_data_3d_reg(ncid,u_id,iwrite,lmax,nx,ny,vnc3d)

	  call fm2am3d(nlv,ilhv,vprv,lmax,nx,ny,fm,value3d)
	  call nc_rewrite_3d_reg(lmax,nx,ny,value3d,vnc3d)
          call nc_write_data_3d_reg(ncid,v_id,iwrite,lmax,nx,ny,vnc3d)
	else
          call nc_write_data_2d(ncid,z_id,iwrite,nkn,znv)

	  call nc_compact_3d(nlv,nlv,nkn,uprv,var3d)
          call nc_write_data_3d(ncid,u_id,iwrite,nlv,nkn,var3d)

	  call nc_compact_3d(nlv,nlv,nkn,vprv,var3d)
          call nc_write_data_3d(ncid,v_id,iwrite,nlv,nkn,var3d)
	end if

	goto 300

  100	continue

c-----------------------------------------------------------------
c end of loop
c-----------------------------------------------------------------

	write(6,*)
	write(6,*) nread,' records read'
	write(6,*)

        if( breg ) call write_dimensions(nx,ny,x0,y0,dx,dy)

        call nc_close(ncid)

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

        stop
   94   continue
        write(6,*) 'incompatible basin and simulation'
        write(6,*) 'nkn,nknous: ',nkn,nknous
        write(6,*) 'nel,nelous: ',nel,nelous
        stop 'error stop ous2nc: parameter mismatch'
   95   continue
        write(6,*) 'regular grid too big'
        write(6,*) 'nx,ny: ',nx,ny,'   nxymax: ',nxymax
        write(6,*) 'please increase nxymax or set to zero for any size'
        stop 'error stop nos2nc: size regular grid'
        end

c******************************************************************

