
!--------------------------------------------------------------------------
!
!    Copyright (C) 2005-2018  Petras Zemlys, Georg Umgiesser
!
!    This file is part of SHYFEM. (m)
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

! AQUABC utilities mostly related to interface with SHYFEM
! Contains:
! subroutine open_scalar_file_bio (not used more)
! subroutine write_scalar_file_bio(not used more)
! subroutine dump_aquabc
! subroutine fluxes_aquabc

c*************************************************************
c*************************************************************

	subroutine open_scalar_file_bio(ia_out,nlv,ilhkv_,hlv_,
     + hev_,nvar,type)

c opens (NOS) file

! Modified by Petras for bottom sediments kinetics

c on return iu = -1 means that no file has been opened and is not written

	!use mod_depth
	!use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer ia_out(4)	!time information		       (in/out)
	integer nlv		!vertical dimension of scalar          (in)
	integer nvar		!total number of variables to write    (in)
	character*(*) type	!extension of file		       (in)

	!include 'param.h' !COMMON_GGU_SUBST
	include 'femtime.h'

	include 'simul.h'
! 	include 'nbasin.h'
! 	include 'levels.h'
! 	include 'depth.h'

	integer nvers
	integer date,time
	integer iu,ierr
	character*80 title,femver

	integer ifemop
	double precision dgetpar

	!Variables to run for bottom sediment kinetics	
	!integer ilhkv_(nkn)
	integer ilhkv_(1309)
	!real hlv_(nlv)
	real hlv_(1)
	!real hev_(nel)
	real hev_(2024)
    
c-----------------------------------------------------
c open file
c-----------------------------------------------------

	iu = ifemop(type,'unformatted','new')
	if( iu .le. 0 ) goto 98
	ia_out(4) = iu

c-----------------------------------------------------
c initialize parameters
c-----------------------------------------------------

	nvers = 5
	date = nint(dgetpar('date'))
	time = nint(dgetpar('time'))
	title = descrp
	call get_shyfem_version(femver)

c-----------------------------------------------------
c write header of file
c-----------------------------------------------------

	call nos_init(iu,nvers)
	call nos_set_title(iu,title)
	call nos_set_date(iu,date,time)
	call nos_set_femver(iu,femver)
	call nos_write_header(iu,nkn,nel,nlv,nvar,ierr)
        if(ierr.gt.0) goto 99
    
    ! ilhkv(nkndim) from levels  - number of layers for each node
    ! hlv(nlvdim)   from levels  - absolute(from surface) depth of layer bottom
    ! hev(neldim)   from depth.h - depth of element
       
	call nos_write_header2(iu,ilhkv_,hlv_,hev_,ierr)
        if(ierr.gt.0) goto 99

c-----------------------------------------------------
c write informational message to terminal
c-----------------------------------------------------

        write(6,*) 'open_scalar_file_bio: ',type,' file opened ',it

c-----------------------------------------------------
c end of routine
c-----------------------------------------------------

	return
   98	continue
	write(6,*) 'error opening file with type ',type
	stop 'error stop open_scalar_file_bio'
   99	continue
	write(6,*) 'error ',ierr,' writing file with type ',type
	stop 'error stop open_scalar_file_bio'
	end

c*************************************************************
c*************************************************************
c*************************************************************

	subroutine write_scalar_file_bio(ia_out,ivar,nlvdi,ilhkv_,c)

c writes NOS file
! Modified by Petras to 
c
c the file must be open, the file will be written unconditionally

	use basin, only: nkn

	implicit none

	integer ia_out(4)	!time information
	integer ivar		!id of variable to be written
	integer nlvdi		!vertical dimension of c
	real c(nlvdi,1)		!scalar to write

	!include 'param.h' !COMMON_GGU_SUBST
	include 'femtime.h'

	!include 'levels.h'

	logical binfo
	integer iu,ierr
	
	!Variables to run for bottom sediment kinetics		
	integer ilhkv_(nkn)
	
	binfo = .false.



c-----------------------------------------------------
c check if files has to be written
c-----------------------------------------------------

	iu = ia_out(4)
	if( iu .le. 0 ) return

c-----------------------------------------------------
c write file
c-----------------------------------------------------

	call nos_write_record(iu,it,ivar,nlvdi,ilhkv_,c,ierr)
	if(ierr.gt.0) goto 99

c-----------------------------------------------------
c write informational message to terminal
c-----------------------------------------------------

	if( binfo ) then
      write(6,*) 'write_scalar_file_bio: ivar = ',
     + ivar,' written at ',it
	end if

c-----------------------------------------------------
c end of routine
c-----------------------------------------------------

	return
   99	continue
	write(6,*) 'error ',ierr,' writing file at unit ',iu
	stop 'error stop write_scalar_file_bio: error in writing record'
	end

c*************************************************************
c*************************************************************



c********************************************************************
c********************************************************************
c   Dumps state variables for repeated runs
c
       subroutine dump_aquabc(file_name,state,
     +                 nkndim,nkn,nlvdim,nvar)

      implicit none
      character*(*) file_name
      integer nkndim,nkn,nlvdim,nvar
      integer un
      real state(nlvdim,nkndim,nvar)

      integer NUM_COLS, FILE_TYPE
      character * 5  NUM_STRING
      character * 30 FORMAT_STRING
      integer ifileo
      integer i,j,k

      NUM_COLS = nlvdim
      write(NUM_STRING,100) NUM_COLS
      FORMAT_STRING = '(' // NUM_STRING // 'F20.8)'

      un = ifileo(55,file_name,'form','unknown')

      write(un,*) '***************************************'
      write(un,*) '* DUMP FILE CREATED BY SHYFEM-AQUABC  *'
      write(un,*) '* !!! PLEASE DO NOT EDIT MANUALLY !!! *'
      write(un,*) '***************************************'

      FILE_TYPE = 4

      write(un,101) nkn, NUM_COLS, nvar, FILE_TYPE


      do i = 1, nvar
          do j = 1, nkn
              !write(un,FORMAT_STRING) (state(k,j,i),k=1,NUM_COLS)
              write(un,*) (state(k,j,i),k=1,NUM_COLS)
          end do
      end do


  100 format(i5)
  101 format(4i5)
      close(un)

      return
      end

c**********************************************************************
c**********************************************************************


c**********************************************************************
c**********************************************************************
c**********************************************************************
c**********************************************************************

	subroutine fluxes_aquabc(ext,ivbase,nscal,scal)
	
! created from template located in subflxa.f

! administers writing of flux data
!
! serves as a template for new variables
! please adapt to your needs
!
! this routine must be called after the other fluxes have already been set up
!
! here the number of scalars and the scalar values are passed into the routine
! you can also import them through other means (modules, etc..)
!
	use levels, only : nlvdi,nlv
	use basin, only : nkn
	use flux

	implicit none

	include 'simul.h'

	character(len=4) :: ext		!extension of new file (e.g., '.fxw')
	integer ivbase			!base for variable numbers
	integer nscal			!how many tracers to compute/write
	real scal(nlvdi,nkn,nscal)	!tracer

	integer i,nlmax,ivar,nvers
	integer idtflx
	integer nvar,ierr
	integer kext(kfluxm)
	real az,dt
	double precision atime,atime0
	character*80 title,femver

        double precision, save :: da_out(4)
        integer, save :: nbflx = 0

	real, save, allocatable :: trs(:)
	real, save, allocatable :: scalt(:,:,:,:)	!accumulator array

	integer ifemop,ipext
	logical has_output_d,next_output_d,is_over_output_d
	double precision dgetpar

!-----------------------------------------------------------------
! start of code
!-----------------------------------------------------------------

        if( nbflx .eq. -1 ) return

!-----------------------------------------------------------------
! initialization
!-----------------------------------------------------------------

        if( nbflx .eq. 0 ) then

                call init_output_d('itmflx','idtflx',da_out)
                call increase_output_d(da_out)
                if( .not. has_output_d(da_out) ) nbflx = -1

                if( kfluxm .le. 0 ) nbflx = -1
                if( nsect .le. 0 ) nbflx = -1
                if( nbflx .eq. -1 ) return

		if( .not. bflxinit ) goto 94

        	allocate(trs(nscal))
        	allocate(scalt(0:nlvdi,3,nsect,nscal))

        	call flux_alloc_arrays(nlvdi,nsect)
		call get_nlayers(kfluxm,kflux,nlayers,nlmax)

		do i=1,nscal
		  call fluxes_init(nlvdi,nsect,nlayers,trs(i)
     +				,scalt(0,1,1,i))
		end do

                nbflx = ifemop(ext,'unform','new')
                if( nbflx .le. 0 ) goto 99
		write(6,*) 'flux file opened: ',nbflx,' ',ext
		da_out(4) = nbflx

	        nvers = 0
		nvar = nscal
		idtflx = nint(da_out(1))
                call flx_write_header(nbflx,0,nsect,kfluxm,idtflx
     +                                  ,nlmax,nvar,ierr)
                if( ierr /= 0 ) goto 98

                title = descrp
                call get_shyfem_version_and_commit(femver)
                call get_absolute_ref_time(atime0)

                do i=1,kfluxm
                  kext(i) = ipext(kflux(i))
                end do

                call flx_write_header2(nbflx,0,nsect,kfluxm
     +                          ,kext,nlayers
     +                          ,atime0,title,femver,chflx,ierr)
                if( ierr /= 0 ) goto 98

        end if

!-----------------------------------------------------------------
! normal call
!-----------------------------------------------------------------

        if( .not. is_over_output_d(da_out) ) return

	call get_timestep(dt)
	call getaz(az)

!	-------------------------------------------------------
!	accumulate results
!	-------------------------------------------------------

	do i=1,nscal
	  ivar = ivbase + i
	  call flxscs(kfluxm,kflux,iflux,az,fluxes,ivar,scal(1,1,i))
	  call fluxes_accum(nlvdi,nsect,nlayers,dt,trs(i)
     +			,scalt(0,1,1,i),fluxes)
	end do

!	-------------------------------------------------------
!	time for output?
!	-------------------------------------------------------

        if( .not. next_output_d(da_out) ) return

!	-------------------------------------------------------
!	average and write results
!	-------------------------------------------------------

        call get_absolute_act_time(atime)

	do i=1,nscal
	  ivar = ivbase + i
	  call fluxes_aver(nlvdi,nsect,nlayers,trs(i)
     +			,scalt(0,1,1,i),fluxes)
          call flx_write_record(nbflx,nvers,atime,nlvdi,nsect,ivar
     +                          ,nlayers,fluxes,ierr)
          if( ierr /= 0 ) goto 97
	end do

!	-------------------------------------------------------
!	reset variables
!	-------------------------------------------------------

	do i=1,nscal
	  call fluxes_init(nlvdi,nsect,nlayers,trs(i)
     +			,scalt(0,1,1,i))
	end do

!-----------------------------------------------------------------
! end of routine
!-----------------------------------------------------------------

	return
   94   continue
        write(6,*) 'Flux section has not been initialized'
        stop 'error stop fluxes_template: no initialization'
   97   continue
        write(6,*) 'Error writing data record of FLX file'
        write(6,*) 'unit,ierr :',nbflx,ierr
        stop 'error stop fluxes_template: writing flx record'
   98   continue
        write(6,*) 'Error writing headers of FLX file'
        write(6,*) 'unit,ierr :',nbflx,ierr
        stop 'error stop fluxes_template: writing flx header'
   99	continue
	write(6,*) 'extension: ',ext
        stop 'error stop fluxes_template: cannot open flx file'
	end

!******************************************************************
!******************************************************************




c**********************************************************************
! This is commented old version
! 	subroutine fluxes_aquabc(it,nscal,scal)
! 
! c  writing aquabc variables flux data to *.csc
! c
! c created from template located in subflxa.f
! c
! c please copy to extra file and adapt to your needs
! c
! c this routine must be called after the other fluxes have already been set up
! c
! c here the number of scalars and the scalar values are passed into the routine
! c you can also import them through other means (modules, etc..)
! c
! c to change for adaptation:
! c
! c ext		extension for file
! c ivar_base	base of variable numbering
! 
! 	use levels, only : nlvdi,nlv
! 	use basin, only : nkn
! 	use flux
! 
! 	implicit none
! 
! 	integer it			!time
! 	integer nscal			!how many tracers to compute/write
! 	real scal(nlvdi,nkn,nscal)	!tracer
! 
! 	integer itend
! 	integer i,nlmax,ivar,nvers
! 	integer idtflx
! 	real az,azpar
! 
!         integer, save :: ia_out(4)
!         integer, save :: nbflx = 0
! 
! 	integer, parameter :: ivar_base = 200	!base of variable numbering
! 	character*4, parameter :: ext = '.csc'	!extension for file
! 
! 	integer, save, allocatable :: nrs(:)
! 	real, save, allocatable :: scalt(:,:,:,:)	!accumulator array
! 
! 	integer ifemop
! 	logical has_output,next_output,is_over_output
! 	double precision dgetpar
! 
! c-----------------------------------------------------------------
! c start of code
! c-----------------------------------------------------------------
! 
!         if( nbflx .eq. -1 ) return
! 
! c-----------------------------------------------------------------
! c initialization
! c-----------------------------------------------------------------
! 
!         if( nbflx .eq. 0 ) then
! 
!                 call init_output('itmflx','idtflx',ia_out)
!                 call increase_output(ia_out)
!                 if( .not. has_output(ia_out) ) nbflx = -1
! 
!                 if( kfluxm .le. 0 ) nbflx = -1
!                 if( nsect .le. 0 ) nbflx = -1
!                 if( nbflx .eq. -1 ) return
! 
!         	allocate(nrs(nscal))
!         	allocate(scalt(0:nlvdi,3,nsect,nscal))
! 
!         	call flux_alloc_arrays(nlvdi,nsect)
! 		call get_nlayers(kfluxm,kflux,nlayers,nlmax)
! 
! 		do i=1,nscal
! 		  call fluxes_init(nlvdi,nsect,nlayers
!      +				,nrs(i),scalt(0,1,1,i))
! 		end do
! 
!                 nbflx = ifemop(ext,'unform','new')
!                 if( nbflx .le. 0 ) goto 99
! 		write(6,*) 'fluxes_aquabc: flux file opened: ',nbflx,' ',ext
! 
! 	        nvers = 5
! 		idtflx = ia_out(1)
!                 call wfflx      (nbflx,nvers
!      +                          ,nsect,kfluxm,idtflx,nlmax
!      +                          ,kflux
!      +                          ,nlayers
!      +                          )
! 
!         end if
! 
! c-----------------------------------------------------------------
! c normal call
! c-----------------------------------------------------------------
! 
!         if( .not. is_over_output(ia_out) ) return
! 
! 	call getaz(azpar)
! 	az = azpar
! 
! c	-------------------------------------------------------
! c	accumulate results
! c	-------------------------------------------------------
! 
! 	do i=1,nscal
! 	  ivar = ivar_base + i
! 	  call flxscs(kfluxm,kflux,iflux,az,fluxes,ivar,scal(1,1,i))
! 	  call fluxes_accum(nlvdi,nsect,nlayers
!      +			,nrs(i),scalt(0,1,1,i),fluxes)
! 	end do
! 
! c	-------------------------------------------------------
! c	time for output?
! c	-------------------------------------------------------
! 
!         if( .not. next_output(ia_out) ) return
! 
! c	-------------------------------------------------------
! c	average and write results
! c	-------------------------------------------------------
! 
! 	do i=1,nscal
! 	  ivar = ivar_base + i
! 	  call fluxes_aver(nlvdi,nsect,nlayers
!      +			,nrs(i),scalt(0,1,1,i),fluxes)
! 	  call wrflx(nbflx,it,nlvdi,nsect,ivar,nlayers,fluxes)
! 	end do
! 
! c	-------------------------------------------------------
! c	reset variables
! c	-------------------------------------------------------
! 
! 	do i=1,nscal
! 	  call fluxes_init(nlvdi,nsect,nlayers
!      +			,nrs(i),scalt(0,1,1,i))
! 	end do
! 
! c-----------------------------------------------------------------
! c end of routine
! c-----------------------------------------------------------------
! 
! 	return
!    99	continue
! 	write(6,*) 'extension: ',ext
!         stop 'fluxes_aquabc: Cannot open fluxes file'
! 	end
! 
! c******************************************************************
! c******************************************************************
