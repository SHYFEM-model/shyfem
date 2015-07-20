c
c $Id: nos2nc.f,v 1.15 2009-11-18 16:50:37 georg Exp $
c
c convert NOS to NC file
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
c
c***************************************************************

	program nos2nc

c reads nos file and writes NetCDF file

	use mod_depth
	use mod_hydro
	use evgeom
	use levels
	use basin

	implicit none

        include 'param.h'

c-------------------------------------------------

	integer ndim
	parameter (ndim=100)

	integer nxdim,nydim
	parameter (nxdim=400,nydim=400)

	character*80 title

	real cv3(nlvdim,nkndim)
	integer ivars(ndim)
	integer var_ids(ndim)



	real uprv(nlvdim,nkndim)
	real vprv(nlvdim,nkndim)
	real ut2v(neldim)
	real vt2v(neldim)
	real u2v(neldim)
	real v2v(neldim)

	real haux(nkndim)
	real weight(nlvdim,nkndim)
	real hl(nlvdim)

	real var3d(nlvdim*nkndim)

        integer nvers,nin,lmax,l
        integer itanf,itend,idt
	integer it,ie,i
        integer ierr,nread,ndry
	integer irec,maxrec,iwrite
        integer nknnos,nelnos
	integer iztype
	integer nvar,ivar
        real href,hzoff,hlvmin
	real volume
	real zmin,zmax
	real umin,umax
	real vmin,vmax
	real flag

	integer nx,ny
	real xlon(nxdim)
	real ylat(nydim)
	real depth(nxdim,nydim)
	real value2d(nxdim,nydim)
	real value3d(nlvdim,nxdim,nydim)
	real vnc3d(nxdim,nydim,nlvdim)
	real fm(4,nxdim,nydim)
	real x0,y0,dx,dy

	logical breg
	logical bdate
        integer ncid
        integer dimids_2d(2)
        integer coord_varid(3)
        integer rec_varid
        integer var_id
	integer date0,time0
	integer date,time
	integer it0
	integer iperiod,its,ite,nfreq
	logical bwrite
	logical bfirst

	character*80 units,std

	integer iapini,ideffi

	call shyfem_copyright('nos2nc - netcdf output')

c-----------------------------------------------------------------
c initialize basin and simulation
c-----------------------------------------------------------------

	maxrec = 2		!max number of records to be written
	maxrec = 0		!max number of records to be written (0 -> all)

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

	if(iapini(3,nkndim,neldim,0).eq.0) then
		stop 'error stop : iapini'
	end if

	call set_ev

	call makehkv_minmax(hkv,haux,1)
	call makehev(hev)

c-----------------------------------------------------------------
c get input from terminal about what to do
c-----------------------------------------------------------------

	if( bdate ) call read_date_and_time(date0,time0)
	call dtsini(date0,time0)

        call get_dimensions(nxdim,nydim,nx,ny,x0,y0,dx,dy,xlon,ylat)
        breg = nx .gt. 0 .and. ny .gt. 0          !regular output
        write(6,*) 'breg: ',breg
        if( breg ) then
	  write(6,*) 'NETCDF output: regular ',dx,dy,nx,ny
          call setgeo(x0,y0,dx,dy,flag)
          call av2fm(fm,nx,ny)
	else
	  write(6,*) 'NETCDF output: unstructured (fem)'
        end if

	call get_period(iperiod,its,ite,nfreq)

c-----------------------------------------------------------------
c read header of simulation
c-----------------------------------------------------------------

        call open_nos_type('.nos','old',nin)

        call read_nos_header(nin,nkndim,neldim,nlvdim,ilhkv,hlv,hev)
        call nos_get_params(nin,nknnos,nelnos,nlv,nvar)
	call nos_get_date(nin,date,time)
	if( date .gt. 0 ) then
	  date0 = date
	  time0 = time
	  call dtsini(date0,time0)
	end if

	if( nkn .ne. nknnos .or. nel .ne. nelnos ) goto 94
	if( nvar .gt. ndim ) goto 95

	call init_sigma_info(nlv,hlv)
	call level_k2e(nkn,nel,nen3v,ilhkv,ilhv)
	call compute_iztype(iztype)

	call nos_get_vars(nin,nvar,ivars)

        write(6,*) 'Available variables: ',nvar
        write(6,*) (ivars(i),i=1,nvar)

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

	do i=1,nvar
	  call nc_init_variable(ncid,breg,3,ivars(i),flag,var_ids(i))
	end do

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

	call rdnos(nin,it,ivar,nlvdim,ilhkv,cv3,ierr)

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
	  call fm2am3d(nlvdim,ilhv,cv3,lmax,nx,ny,fm,value3d)
	  call nc_rewrite_3d_reg(lmax,nx,ny,value3d,vnc3d)
	  call nc_write_data_3d_reg(ncid,var_id,iwrite,lmax,nx,ny,vnc3d)
	else
	  call nc_compact_3d(nlvdim,nlv,nkn,cv3,var3d)
          call nc_write_data_3d(ncid,var_id,iwrite,nlv,nkn,var3d)
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
   92   continue
        write(6,*) 'wrong order of variables'
        write(6,*) 'read: ',ivar,'   expected: ',ivars(i)
        stop 'error stop nos2nc: variables'
   94   continue
        write(6,*) 'incompatible simulation and basin'
        write(6,*) 'nkn: ',nkn,nknnos
        write(6,*) 'nel: ',nel,nelnos
        stop 'error stop nos2nc: nkn,nel'
   95   continue
        write(6,*) 'nvar,ndim: ',nvar,ndim
        stop 'error stop nos2nc: ndim'
        end

c******************************************************************

