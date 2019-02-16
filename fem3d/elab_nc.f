
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

! netcdf routines for shyelab: elab_nc
!
! revision log :
!
! 25.10.2018    ccf     write field already on regular grid
!
!************************************************************

	subroutine nc_output_init(ncid,title,nvar,ivars)

	use basin
	use levels
	use mod_depth
	use shyelab_out
        use elabutil

	implicit none

	integer ncid			!id of file (return)
	character*(*) title		!name of simulation
	integer nvar			!total number of variables to be written
	integer ivars(nvar)		!variable id of SHYFEM

	logical bhydro
	integer date0,time0
	integer lmax,i,iztype,idim,ivar,ncnlv
	real, save :: ncflag = -999.

	title = 'adriatic tiresias'
	call dts_get_date(date0,time0)
	call compute_iztype(iztype)

        ncnlv = nlv
	if ( b2d ) ncnlv = 1

	if( breg ) then
	  allocate(value2d(nxreg,nyreg))
	  allocate(value3d(ncnlv,nxreg,nyreg))
	  allocate(vnc3d(nxreg,nyreg,ncnlv))
	  call get_lmax_reg(nxreg,nyreg,fmreg,ilhv,lmax)
	  lmax = max(1,lmax)	!at least one layer
	  if ( b2d ) lmax = 1
	  lmaxreg = lmax
	  call nc_open_reg(ncid,nxreg,nyreg,lmax
     +				,ncflag,date0,time0,iztype)
	else
	  allocate(var3d(ncnlv*nkn))
	  call nc_open_fem(ncid,nkn,nel,ncnlv,date0,time0,iztype)
	end if

	call nc_global(ncid,title)

	bhydro = ivars(1) == 1		!this is a hydro file
	allocate(var_ids(nvar))
	var_ids = 0

        do i=1,nvar
	  ivar = ivars(i)
	  idim = 3
	  if( bhydro ) then
	    if( i == 1 ) idim = 2
	    if( i == 2 ) cycle	!do not write second level
	    if( i == 3 ) ivar = 2
	  end if
          call nc_init_variable(ncid,breg,idim,ivar,ncflag,var_ids(i))
        end do

        call nc_end_define(ncid)

        if( breg ) then
          call nc_write_coords_reg(ncid,nxreg,nyreg,ncnlv
     +					,xlon,ylat,hcoord,hlv)
        else
          call nc_write_coords_fem(ncid,nkn,nel,ncnlv,xgv,ygv
     +					,hkv,nen3v,hlv)
        end if

	end

!********************************************************************

	subroutine nc_output_record(ncid,var_id,cv3)

	use basin
	use levels
	use shyelab_out

	implicit none

	integer ncid
	integer var_id
	real cv3(nlvdi,nkn)

	integer lmax,nx,ny,iwrite

	iwrite = iwrite_nc

        if( breg ) then
	  lmax = lmaxreg
	  nx = nxreg
	  ny = nyreg
          call fm2am3d(nlv,ilhv,cv3,lmax,nx,ny,fmreg,value3d)
          call nc_rewrite_3d_reg(lmax,nx,ny,value3d,vnc3d)
          call nc_write_data_3d_reg(ncid,var_id,iwrite,lmax,nx,ny,vnc3d)
        else
          call nc_compact_3d(nlv,nlv,nkn,cv3,var3d)
          call nc_write_data_3d(ncid,var_id,iwrite,nlv,nkn,var3d)
        end if

	end

!********************************************************************
! nc output in which cv3 if breg is already on regular grid and does 
! not need to be reinterpolated

	subroutine nc_output_record_reg(ncid,var_id,nlvd,np,cv3)

	use shyelab_out

	implicit none

	integer ncid
	integer var_id
	integer nlvd
	integer np
	real cv3(nlvd,np)

	integer lmax,nx,ny,iwrite
        integer i,j,l,ii

	iwrite = iwrite_nc

        if( breg ) then
          ii = 0
	  lmax = lmaxreg
	  nx = nxreg
	  ny = nyreg
          do j=1,ny
            do i=1,nx
              ii = ii + 1
              do l=1,lmax
                vnc3d(i,j,l) = cv3(l,ii)
              end do
            end do
          end do
          call nc_write_data_3d_reg(ncid,var_id,iwrite,lmax,nx,ny,vnc3d)
        else
          call nc_compact_3d(nlvd,nlvd,np,cv3,var3d)
          call nc_write_data_3d(ncid,var_id,iwrite,nlvd,np,var3d)
        end if

	end

!********************************************************************

	subroutine nc_output_hydro(ncid,znv,uprv,vprv)

	use basin
	use levels
	use shyelab_out

	implicit none

	integer ncid
	real znv(nkn)
	real uprv(nlvdi,nkn)
	real vprv(nlvdi,nkn)

	integer lmax,nx,ny,var_id,iwrite

	iwrite = iwrite_nc

        if( breg ) then
	  lmax = lmaxreg
	  nx = nxreg
	  ny = nyreg

	  var_id = var_ids(1)
          call fm2am2d(znv,nx,ny,fmreg,value2d)
          call nc_write_data_2d_reg(ncid,var_id,iwrite,nx,ny,value2d)

	  var_id = var_ids(3)
          call fm2am3d(nlv,ilhv,uprv,lmax,nx,ny,fmreg,value3d)
          call nc_rewrite_3d_reg(lmax,nx,ny,value3d,vnc3d)
          call nc_write_data_3d_reg(ncid,var_id,iwrite,lmax,nx,ny,vnc3d)

	  var_id = var_ids(4)
          call fm2am3d(nlv,ilhv,vprv,lmax,nx,ny,fmreg,value3d)
          call nc_rewrite_3d_reg(lmax,nx,ny,value3d,vnc3d)
          call nc_write_data_3d_reg(ncid,var_id,iwrite,lmax,nx,ny,vnc3d)
        else
	  var_id = var_ids(1)
          call nc_write_data_2d(ncid,var_id,iwrite,nkn,znv)

	  var_id = var_ids(3)
          call nc_compact_3d(nlv,nlv,nkn,uprv,var3d)
          call nc_write_data_3d(ncid,var_id,iwrite,nlv,nkn,var3d)

	  var_id = var_ids(4)
          call nc_compact_3d(nlv,nlv,nkn,vprv,var3d)
          call nc_write_data_3d(ncid,var_id,iwrite,nlv,nkn,var3d)
        end if

	end

!********************************************************************

	subroutine nc_output_time(ncid,dtime)

	use shyelab_out

	implicit none

	integer ncid
	double precision dtime

	integer, save :: icall = 0
	double precision, save :: dtime_old = 0

	if( icall == 0 ) dtime_old = dtime - 1.
	icall = icall + 1
	if( dtime_old == dtime ) return			!already written

	iwrite_nc = iwrite_nc + 1

	call nc_write_dtime(ncid,iwrite_nc,dtime)

	end

!********************************************************************

	subroutine nc_output_final(ncid)

	implicit none

	integer ncid

	call nc_close(ncid)

	end

!********************************************************************

