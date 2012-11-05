c
c $Id: ous2nc.f,v 1.15 2009-11-18 16:50:37 georg Exp $
c
c info on OUS files
c
c revision log :
c
c 02.09.2003	ggu	adapted to new OUS format
c 24.01.2005	ggu	computes maximum velocities for 3D (only first level)
c 16.10.2007	ggu	new debug routine
c 27.10.2009    ggu     include evmain.h, compute volume
c 23.03.2011    ggu     compute real u/v-min/max of first level
c
c***************************************************************

	program ous2ncreg

c reads ous file and writes NetCDF file (regular)

	implicit none

        include 'param.h'
	include 'evmain.h'

	integer nxdim,nydim
	parameter (nxdim=400,nydim=400)

	character*80 descrr,descrp
	common /descrr/ descrr
	common /descrp/ descrp
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	real xgv(nkndim), ygv(nkndim)
	real hm3v(3,neldim)
	integer nen3v(3,neldim)
	integer ipev(neldim), ipv(nkndim)
	integer iarv(neldim)
	common /xgv/xgv, /ygv/ygv
	common /hm3v/hm3v
	common /nen3v/nen3v
	common /ipev/ipev, /ipv/ipv
	common /iarv/iarv

	integer ilhv(neldim)
	real hlv(nlvdim)
        real utlnv(nlvdim,neldim)
        real vtlnv(nlvdim,neldim)
	common /ilhv/ilhv
	common /hlv/hlv
        common /utlnv/utlnv
        common /vtlnv/vtlnv
	real zenv(3,neldim)
	common /zenv/zenv
	real znv(nkndim)
	common /znv/znv

	real hkv(nkndim)
	common /hkv/hkv
	real hev(neldim)
	common /hev/hev

	real haux(nkndim)
	real uprv(nlvdim,nkndim)
	real vprv(nlvdim,nkndim)

	real var3d(nlvdim*nkndim)

	integer nx,ny
	real xlon(nxdim)
	real ylat(nydim)
	real depth(nxdim,nydim)
	real value2d(nxdim,nydim)
	real value3d(nlvdim,nxdim,nydim)
	real vnc3d(nxdim,nydim,nlvdim)
	real fm(4,nxdim,nydim)

        integer nvers,nin,nlv,lmax
        integer itanf,itend,idt,idtous
	integer it,ie,i
        integer ierr,nread,ndry
	integer irec,maxrec
        integer nknous,nelous,nlvous
        real href,hzoff,hlvmin
	real volume
	real zmin,zmax
	real umin,umax
	real vmin,vmax

	real x0,y0,dx,dy
	real flag

	logical bdate
        integer ncid
        integer dimids_2d(2)
        integer coord_varid(3)
        integer rec_varid
        integer level_id,u_id,v_id
	integer date0,time0
	integer it0

	character*80 units,std

c	integer rdous,rfous
	integer iapini,ideffi

	call shyfem_copyright('ous2cn_reg - regular netcdf output')

c-----------------------------------------------------------------
c initialize basin and simulation
c-----------------------------------------------------------------

	date0 = 2009
	time0 = 0
	maxrec = 2		!max number of records to be written
	maxrec = 0		!max number of records to be written (0 -> all)

	it0 = 0			!subtract from it (hack)
	bdate = .false.		!if true read from file

	flag = -999.

c-----------------------------------------------------------------
c do not change anything below here
c-----------------------------------------------------------------

	nread = 0
	irec = 0

	if(iapini(3,nkndim,neldim,0).eq.0) then
		stop 'error stop : iapini'
	end if

	call get_dimensions(nxdim,nydim,nx,ny,x0,y0,dx,dy,xlon,ylat)
	call setgeo(x0,y0,dx,dy,flag)

	call set_ev

	call makehkv_minmax(hkv,haux,1)
	call makehev(hev)

	call av2fm(fm,nx,ny)

	nin=ideffi('datdir','runnam','.ous','unform','old')
	if(nin.le.0) goto 100

	call dtsini(date0,time0)

c-----------------------------------------------------------------
c read header of simulation
c-----------------------------------------------------------------

	nvers=1
        call rfous(nin
     +			,nvers
     +			,nknous,nelous,nlvous
     +			,href,hzoff
     +			,descrp
     +			,ierr)

	nlv=nlvous
	call dimous(nin,nkndim,neldim,nlvdim)

        write(6,*)
        write(6,*)   descrp
        write(6,*)
        write(6,*) ' nvers        : ',nvers
        write(6,*) ' href,hzoff   : ',href,hzoff
        write(6,*) ' nkn,nel      : ',nknous,nelous
        write(6,*) ' nlv          : ',nlvous
        write(6,*)

	call rsous(nin,ilhv,hlv,hev,ierr)

	call init_sigma_info(nlv,hlv)

	call get_lmax_reg(nx,ny,fm,ilhv,lmax)

c-----------------------------------------------------------------
c prepare netcdf file
c-----------------------------------------------------------------

        call nc_open_reg(ncid,nx,ny,lmax,flag,date0,time0)
	call nc_global(ncid, descrp)

	std = 'water_surface_height_above_reference_datum'
	units = 'm'
	call nc_define_2d_reg(ncid,'water_level',level_id)
	call nc_define_attr(ncid,'units',units,level_id)
	call nc_define_attr(ncid,'standard_name',std,level_id)
	call nc_define_range(ncid,-10.0,+10.0,flag,level_id)

	std = 'eastward_sea_water_velocity_assuming_no_tide'
	units = 'm s-1'
	call nc_define_3d_reg(ncid,'u_velocity',u_id)
	call nc_define_attr(ncid,'units',units,u_id)
	call nc_define_attr(ncid,'standard_name',std,u_id)
	call nc_define_range(ncid,-10.0,+10.0,flag,u_id)

	std = 'northward_sea_water_velocity_assuming_no_tide'
	units = 'm s-1'
	call nc_define_3d_reg(ncid,'v_velocity',v_id)
	call nc_define_attr(ncid,'units',units,v_id)
	call nc_define_attr(ncid,'standard_name',std,v_id)
	call nc_define_range(ncid,-10.0,+10.0,flag,v_id)

        call nc_end_define(ncid)
	call fm2am2d(hkv,nx,ny,fm,depth)
        call nc_write_coords_reg(ncid,nx,ny,xlon,ylat,depth)

c-----------------------------------------------------------------
c loop on data of simulation
c-----------------------------------------------------------------

  300   continue

        call rdous(nin,it,nlvdim,ilhv,znv,zenv,utlnv,vtlnv,ierr)

	it = it - it0

        if(ierr.gt.0) then
		write(6,*) 'error in reading file : ',ierr
		goto 100
        else if(ierr.lt.0) then
		goto 100
	end if

	nread=nread+1

	call mima(znv,nknous,zmin,zmax)
        call comp_vel(1,nel,hev,zenv,nlvdim,utlnv,vtlnv
     +			,umin,vmin,umax,vmax)
	call compute_volume(nel,zenv,hev,volume)

c        call debug_write_node(it,nread,nkndim,neldim,nlvdim,nkn,nel,nlv
c     +          ,nen3v,zenv,znv,utlnv,vtlnv)

	write(6,*) 
	!write(6,*) 'time : ',it
	call write_time(it)
	write(6,*) 
	write(6,*) 'zmin/zmax : ',zmin,zmax
	write(6,*) 'umin/umax : ',umin,umax
	write(6,*) 'vmin/vmax : ',vmin,vmax
	write(6,*) 'volume    : ',volume

        call make_vel_on_nodes(nlv,nkn,nel
     +			,ilhv,nen3v
     +			,utlnv,vtlnv
     +			,uprv,vprv)

        irec = irec + 1
        call nc_write_time(ncid,irec,it)

	call fm2am2d(znv,nx,ny,fm,value2d)
        call nc_write_data_2d_reg(ncid,level_id,irec,nx,ny,value2d)

	call fm2am3d(nlvdim,ilhv,uprv,lmax,nx,ny,fm,value3d)
	call nc_rewrite_3d_reg(lmax,nx,ny,value3d,vnc3d)
	call nc_write_data_3d_reg(ncid,u_id,irec,lmax,nx,ny,vnc3d)

	call fm2am3d(nlvdim,ilhv,vprv,lmax,nx,ny,fm,value3d)
	call nc_rewrite_3d_reg(lmax,nx,ny,value3d,vnc3d)
        call nc_write_data_3d_reg(ncid,v_id,irec,lmax,nx,ny,vnc3d)

	if ( maxrec .gt. 0 .and. irec .ge. maxrec ) goto 100

	goto 300

  100	continue

c-----------------------------------------------------------------
c end of loop
c-----------------------------------------------------------------

	write(6,*)
	write(6,*) nread,' records read'
	write(6,*)

        call nc_close(ncid)

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	stop
	end

c******************************************************************

	subroutine compute_volume(nel,zenv,hev,volume)

	implicit none

	include 'param.h'
	include 'evmain.h'

	integer nel
	real zenv(3,neldim)
	real hev(neldim)
	real volume

	integer ie,ii
	real zav,area
	double precision vol,voltot,areatot

	voltot = 0.
	areatot = 0.

	do ie=1,nel
	  zav = 0.
	  do ii=1,3
	    zav = zav + zenv(ii,ie)
	  end do
	  area = 12. * ev(10,ie)
	  vol = area * (hev(ie) + zav/3.)
	  voltot = voltot + vol
	  !areatot = areatot + area
	end do

	volume = voltot

	end

c******************************************************************

        subroutine make_vel_on_nodes(nlv,nkn,nel
     +			,ilhv,nen3v
     +			,utlnv,vtlnv
     +			,uprv,vprv)

        implicit none

	include 'param.h'
	include 'evmain.h'

        integer nlv,nkn,nel
	integer ilhv(neldim)
	integer nen3v(3,neldim)
        real utlnv(nlvdim,neldim)
        real vtlnv(nlvdim,neldim)
        real uprv(nlvdim,nkndim)
        real vprv(nlvdim,nkndim)

        real auxv(nlvdim,nkndim)
	real hl(nlvdim)

	logical bzeta,bsigma
        integer ie,ii,k,l,lmax
	integer nlvaux,nsigma
	real hsigma
	real area
        real u,v

	bzeta = .true.		!use water levels

	call get_sigma_info(nlvaux,nsigma,hsigma)
	bsigma = nsigma .gt. 0

	do k=1,nkn
	  do l=1,nlv
	    uprv(l,k) = 0.
	    vprv(l,k) = 0.
	    auxv(l,k) = 0.
	  end do
	end do

        do ie=1,nel

	  area = 12. * ev(10,ie)
	  lmax = ilhv(ie)
	  call get_layer_thickness(ie,lmax,bzeta,nsigma,hsigma,hl)

	  do l=1,lmax
            u = utlnv(l,ie) / hl(l)
            v = vtlnv(l,ie) / hl(l)
	    do ii=1,3
	      k = nen3v(ii,ie)
	      uprv(l,k) = uprv(l,k) + area * u
	      vprv(l,k) = vprv(l,k) + area * v
	      auxv(l,k) = auxv(l,k) + area
	    end do
	  end do

        end do

	do k=1,nkn
	  do l=1,nlv
	    area = auxv(l,k)
	    if( area .gt. 0. ) then
	      uprv(l,k) = uprv(l,k) / area
	      vprv(l,k) = vprv(l,k) / area
	    end if
	  end do
	end do

        end

c******************************************************************

        subroutine comp_vel(level,nel,hev,zenv,nlvdim,utlnv,vtlnv
     +			,umin,vmin,umax,vmax)

        implicit none

        integer level
        integer nel
        real hev(1)
        real zenv(3,1)
        integer nlvdim
        real utlnv(nlvdim,1)
        real vtlnv(nlvdim,1)
        real umin,vmin
        real umax,vmax

        integer ie,ii
        real zmed,hmed,u,v

	umin =  1.e+30
	vmin =  1.e+30
        umax = -1.e+30
        vmax = -1.e+30

        do ie=1,nel
          zmed = 0.
          do ii=1,3
            zmed = zmed + zenv(ii,ie)
          end do
          zmed = zmed / 3.
          hmed = hev(ie) + zmed
          if( hmed .le. 0. ) stop 'error stop hmed...'

          u = utlnv(level,ie) / hmed
          v = vtlnv(level,ie) / hmed

          umin = min(umin,u)
          vmin = min(vmin,v)
          umax = max(umax,u)
          vmax = max(vmax,v)
        end do

        end

c******************************************************************

        subroutine debug_write_node(it,nrec
     +		,nkndim,neldim,nlvdim,nkn,nel,nlv
     +          ,nen3v,zenv,znv,utlnv,vtlnv)

c debug write

        implicit none

        integer it,nrec
        integer nkndim,neldim,nlvdim,nkn,nel,nlv
        integer nen3v(3,neldim)
        real znv(nkndim)
        real zenv(3,neldim)
        real utlnv(nlvdim,neldim)
        real vtlnv(nlvdim,neldim)

        integer ie,ii,k,l,ks
        logical bk

        ks = 6068

        write(66,*) 'time: ',it,nrec
        write(66,*) 'kkk: ',znv(ks)

        do ie=1,nel
          bk = .false.
          do ii=1,3
            k = nen3v(ii,ie)
            if( k .eq. ks ) then
              write(66,*) 'ii: ',ii,ie,zenv(ii,ie)
              bk = .true.
            end if
          end do
          if( bk ) then
          do l=1,nlv
            write(66,*) 'ie: ',ie,l,utlnv(l,ie),vtlnv(l,ie)
          end do
          end if
        end do

        end

c******************************************************************

	subroutine get_lmax_reg(nx,ny,fm,ilhv,lmax)

c computes max lmax for regular domain

	implicit none

	integer nx,ny
	real fm(4,nx,ny)
	integer ilhv(1)
	integer lmax		!max level (return)

	integer i,j,ie

	lmax = 0

	do j=1,ny
	  do i=1,nx
	    ie = nint(fm(4,i,j))
	    lmax = max(lmax,ilhv(ie))
	  end do
	end do

	end

c******************************************************************

	subroutine get_dimensions(nxdim,nydim,nx,ny,x0,y0,dx,dy,xlon,ylat)

c gets dimensions for reguar grid

	implicit none

	include 'param.h'

	integer nxdim,nydim,nx,ny
	real x0,y0,dx,dy
	real xlon(nxdim)
	real ylat(nydim)

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	real xgv(nkndim), ygv(nkndim)
	common /xgv/xgv, /ygv/ygv

	integer ichoose,i
	real x1,y1
	real xmin,ymin,xmax,ymax

	call mima(xgv,nkn,xmin,xmax)
	call mima(ygv,nkn,ymin,ymax)

	write(6,*) 'xmin/xmax: ',xmin,xmax
	write(6,*) 'ymin/ymax: ',ymin,ymax

	write(6,*) 'enter dx,dy: '
	read(5,*) dx,dy
	write(6,*) 'Want to choose domain? yes/no -> 1/0'
	read(5,'(i10)') ichoose

	if( ichoose .eq. 1 ) then
	  write(6,*) 'enter x0,y0,x1,y1 (min,max)'
	  read(5,*) x0,y0,x1,y1
	else
	  x0 = dx * (int(xmin/dx))
	  y0 = dy * (int(ymin/dy))
	  x1 = dx * (int(xmax/dx)+1)
	  y1 = dy * (int(ymax/dy)+1)
	end if

	nx = 1 + nint((x1-x0)/dx)
	ny = 1 + nint((y1-y0)/dy)

	write(6,*) 'dx,dy: ',dx,dy
	write(6,*) 'x0,y0,x1,y1: ',x0,y0,x1,y1
	write(6,*) 'nx,ny: ',nx,ny

	if( nx .gt. nxdim ) goto 99
	if( ny .gt. nydim ) goto 99

	do i=1,nx
	  xlon(i) = x0 + (i-1)*dx
	end do

	do i=1,ny
	  ylat(i) = y0 + (i-1)*dy
	end do

	return
   99	continue
	write(6,*) 'nx,nxdim: ',nx,nxdim
	write(6,*) 'ny,nydim: ',ny,nydim
	stop 'error stop get_dimensions: nx/ydim'
	end

c******************************************************************

	subroutine write_time(it)

	implicit none

	integer it

	character*40 line

	call dtsgf(it,line)
	write(6,*) 'time: ',it,'   ',line

	end

c******************************************************************

        subroutine read_date_and_time(date,time)

        implicit none

        integer date,time

        open(1,file='date0',status='old',form='formatted')
        read(1,*) date
        close(1)

        open(1,file='time0',status='old',form='formatted')
        read(1,*) time
        close(1)

        end

c******************************************************************

