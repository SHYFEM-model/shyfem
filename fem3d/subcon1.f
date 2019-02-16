
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

c routines for concentration (utilities) (old newcon1.f)
c
c contents :
c
c subroutine conini0(nlvddi,c,cref)                      sets initial conditions
c subroutine conini(nlvddi,c,cref,cstrat)		sets initial conditions
c
c subroutine conbnd(nlvddi,c,rbc)			boundary condition
c subroutine con3bnd(nlvddi,c,nlvbnd,rbc)		boundary condition (3D)
c
c subroutine confop(iu,itmcon,idtcon,nlv,nvar,type)	opens (NOS) file
c subroutine confil(iu,itmcon,idtcon,ivar,nlvddi,c)	writes NOS file
c subroutine conwrite(iu,type,nvar,ivar,nlvddi,c)        shell for writing file
c
c subroutine conmima(nlvddi,c,cmin,cmax)                 computes min/max
c subroutine conmimas(nlvddi,c,cmin,cmax)                computes scalar min/max
c
c notes :
c
c conbnd and con3bnd are not used and can be deleted
c
c revision log :
c
c 19.08.1998	ggu	call to conzfi changed
c 20.08.1998	ggu	makew removed (routine used is sp256w)
c 24.08.1998    ggu     levdbg used for debug
c 26.08.1998    ggu     subroutine convol, tstvol transferred to newchk
c 26.08.1998    ggu     all subroutines re-written more generally
c 26.01.1999    ggu     can be used also with 2D routines
c 16.11.2001    ggu     subroutine conmima and diffstab
c 05.12.2001    ggu     new routines diffstab,diffstab1,difflimit
c 11.10.2002    ggu     commented diffset
c 09.09.2003    ggu     new routine con3bnd
c 10.03.2004    ggu     new routine conwrite()
c 13.03.2004    ggu     new routines set_c_bound, distribute_vertically
c 13.03.2004    ggu     exec routine con3bnd() only for level BC (LEVELBC)
c 14.03.2004    ggu     new routines open_b_flux
c 05.01.2005    ggu     routine to write 2d nos file into subnosa.f
c 07.01.2005    ggu     routine diffwrite deleted
c 14.01.2005    ggu     new file for diffusion routines (copied to subdif.f)
c 23.03.2006    ggu     changed time step to real
c 31.05.2007    ggu     reset BC of flux type to old way (DEBHELP)
c 07.04.2008    ggu     deleted set_c_bound
c 08.04.2008    ggu     cleaned, deleted distribute_vertically, open_b_flux
c 09.10.2008    ggu&ccf call to confop changed -> nlv
c 20.11.2009    ggu	in conwrite only write needed (nlv) layers
c 20.01.2014    ggu	new writing format for nos files in confop, confil
c 03.11.2017	ggu	new routines to write shy files scalar_output_*()
c
c*****************************************************************

	subroutine conini(nlvddi,c,cref,cstrat,hdko)

c sets initial conditions (with stratification)

	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer nlvddi		!vertical dimension of c
	real c(nlvddi,nkn)	!variable to initialize
	real cref		!reference value
	real cstrat		!stratification [conc/km]
	real hdko(nlvddi,nkn)	!layer thickness
c local
	integer k,l
	real depth,hlayer

	do k=1,nkn
	  depth=0.
	  do l=1,nlvddi
	    hlayer = 0.5 * hdko(l,k)
	    depth = depth + hlayer
	    c(l,k) = cref + cstrat*depth/1000.
	    depth = depth + hlayer
	  end do
	end do

	end

c*************************************************************
c*************************************************************
c*************************************************************

	subroutine scalar_output_init(da_out,nl,nvar,type,ierr)

! opens scalar file for write - da_out already set, return ierr
!
! ierr: 0=ok  -1=no output  1=error

	implicit none

	double precision da_out(4)	!info on output frequency (already set)
	integer nl			!vertical dimension of scalar        
	integer nvar			!total number of variables to write 
	character*(*) type		!type of file (extension)
	integer ierr			!error code (return)

	logical b2d
	integer id
	logical has_output_d

	ierr = -1			!no output

        if( has_output_d(da_out) ) then
	  b2d = ( nl <= 1 )
          call shyfem_init_scalar_file(type,nvar,b2d,id)
          da_out(4) = id
	  ierr = 0			!success
	  if( id <= 0 ) ierr = 1	!error
        end if

	end

c*************************************************************

	subroutine scalar_output_open(dtanf,ddt,nl,nvar,type,da_out,ierr)

! opens scalar file for write - returns da_out and ierr
!
! ierr: 0=ok  -1=no output  1=error

	implicit none

	double precision dtanf		!time of first write		       
	double precision ddt		!time intervall of writes	      
	integer nl			!vertical dimension of scalar        
	integer nvar			!total number of variables to write 
	character*(*) type		!type of file (extension)
	double precision da_out(4)	!info on output frequency (return)
	integer ierr			!error code (return)

	logical b2d
	integer id
	logical has_output_d

	call set_output_frequency_d(dtanf,ddt,da_out)
	call scalar_output_init(da_out,nl,nvar,type,ierr)

	end

c*************************************************************

	subroutine scalar_output_write(dtime,da_out,ivar,nlvddi,val)

! writes scalar file

	implicit none

	double precision dtime
	double precision da_out(4)
	integer ivar
	integer nlvddi
	real val(nlvddi,*)

	integer id
	logical next_output_d

        if( next_output_d(da_out) ) then
          id = nint(da_out(4))
          call shy_write_scalar_record(id,dtime,ivar,nlvddi,val)
        end if

	end

c*************************************************************

	subroutine scalar_output_file(da_out,type,nvar,ivar,nlvddi,val)

c shell for writing file unconditionally to disk

	use levels, only : nlvdi,nlv

        implicit none

	double precision da_out(4)
        character*(*) type      !type of file
	integer nvar		!total number of variables
	integer ivar		!id of variable to be written
	integer nlvddi		!vertical dimension of c
	real val(nlvddi,*)	!concentration to write

        integer id,ierr
	double precision dtime,dtime0,ddt

	ddt = 1.		!this fakes write at this time step
	call get_first_dtime(dtime0)
	call get_act_dtime(dtime)

	id = nint(da_out(4))
        if( id .eq. 0 ) then
	  call scalar_output_open(dtime0,ddt,nlvddi,nvar,type,da_out,ierr)
	  if( ierr > 0 ) then
	    write(6,*) 'error opening file: ',trim(type)
	    stop 'error stop : error opening file'
	  end if
          id = nint(da_out(4))
        end if

        call shy_write_scalar_record(id,dtime,ivar,nlvddi,val)

	end

c*************************************************************

	subroutine scalar_output_once_2d(type,ivar,val)

c shell for writing file unconditionally to disk

	use levels, only : nlvdi,nlv

        implicit none

        character*(*) type      !type of file
	integer ivar		!id of variable to be written
	real val(*)		!concentration to write

	integer nvar		!total number of variables
	integer nlvddi		!vertical dimension of c
	double precision da_out(4)

        integer id,ierr

	nvar = 1
	nlvddi = 1
	da_out = 0

	call scalar_output_file(da_out,type,nvar,ivar,nlvddi,val)
        id = nint(da_out(4))
	call shy_close_output_file(id)

	end

c*************************************************************
c*************************************************************
c*************************************************************

	subroutine confop(iu,itmcon,idtcon,nl,nvar,type)

c opens (NOS) file

c on return iu = -1 means that no file has been opened and is not written

	use mod_depth
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer iu		!unit				       (in/out)
	integer itmcon		!time of first write		       (in/out)
	integer idtcon		!time intervall of writes	       (in/out)
	integer nl		!vertical dimension of scalar          (in)
	integer nvar		!total number of variables to write    (in)
	character*(*) type	!extension of file		       (in)

	include 'simul.h'

	integer nvers
	integer date,time
	integer ierr
	integer itcon
	double precision dtime
	!character*80 dir,nam,file
	character*80 title,femver

	integer ifemop
	real getpar
	double precision dgetpar

c-----------------------------------------------------
c check idtcon and itmcon and adjust
c-----------------------------------------------------

	call adjust_itmidt(itmcon,idtcon,itcon)

	iu = -1
        if( idtcon .le. 0 ) return

c-----------------------------------------------------
c open file
c-----------------------------------------------------

	iu = ifemop(type,'unformatted','new')
	if( iu .le. 0 ) goto 98

c-----------------------------------------------------
c initialize parameters
c-----------------------------------------------------

	nvers = 5
	date = nint(dgetpar('date'))
	time = nint(dgetpar('time'))
	title = descrp
	call get_shyfem_version_and_commit(femver)

c-----------------------------------------------------
c write header of file
c-----------------------------------------------------

	call nos_init(iu,nvers)
	call nos_set_title(iu,title)
	call nos_set_date(iu,date,time)
	call nos_set_femver(iu,femver)
	call nos_write_header(iu,nkn,nel,nl,nvar,ierr)
        if(ierr.gt.0) goto 99
	call nos_write_header2(iu,ilhkv,hlv,hev,ierr)
        if(ierr.gt.0) goto 99

c-----------------------------------------------------
c write informational message to terminal
c-----------------------------------------------------

	call get_act_dtime(dtime)
        write(6,*) 'confop: ',type,' file opened ',dtime

c-----------------------------------------------------
c end of routine
c-----------------------------------------------------

	return
   98	continue
	write(6,*) 'error opening file with type ',type
	stop 'error stop confop'
   99	continue
	write(6,*) 'error ',ierr,' writing file with type ',type
	stop 'error stop confop'
	end

c*************************************************************

	subroutine confil(iu,itmcon,idtcon,ivar,nlvddi,c)

c writes NOS file

	use levels

	implicit none

	integer iu		!unit
	integer itmcon		!time of first write
	integer idtcon		!time intervall of writes
	integer ivar		!id of variable to be written
	integer nlvddi		!vertical dimension of c
	real c(nlvddi,*)		!scalar to write

	include 'femtime.h'

	logical binfo
	integer ierr

	binfo = .false.

c-----------------------------------------------------
c check if files has to be written
c-----------------------------------------------------

	if( iu .le. 0 ) return
	if( it .lt. itmcon ) return
	if( mod(it-itmcon,idtcon) .ne. 0 ) return

c-----------------------------------------------------
c write file
c-----------------------------------------------------

	call nos_write_record(iu,it,ivar,nlvddi,ilhkv,c,ierr)
	if(ierr.gt.0) goto 99

c-----------------------------------------------------
c write informational message to terminal
c-----------------------------------------------------

	if( binfo ) then
          write(6,*) 'confil: variable ',ivar,' written at ',it
	end if

c-----------------------------------------------------
c end of routine
c-----------------------------------------------------

	return
   99	continue
	write(6,*) 'error ',ierr,' writing file at unit ',iu
	stop 'error stop confil'
	end

c*************************************************************

	subroutine conwrite(iu,type,nvar,ivar,nlvddi,c)

c shell for writing file unconditionally to disk

	use levels, only : nlvdi,nlv

        implicit none

	integer iu		!unit (0 for first call, set on return)
        character*(*) type      !type of file
	integer nvar		!total number of variables
	integer ivar		!id of variable to be written
	integer nlvddi		!vertical dimension of c
	real c(nlvddi,*)	!concentration to write

	include 'femtime.h'

        integer itmcon,idtcon,lmax

        itmcon = itanf
	idtcon = idt
	idtcon = 1
	lmax = min(nlvddi,nlv)

        if( iu .eq. 0 ) then
	  call confop(iu,itmcon,idtcon,lmax,nvar,type)
        end if

	call confil(iu,itmcon,idtcon,ivar,nlvddi,c)

        end

c*************************************************************
c*************************************************************
c*************************************************************

	subroutine open_scalar_file(ia_out,nl,nvar,type)
	implicit none
	integer ia_out(4)	!time information		       (in/out)
	integer nl		!vertical dimension of scalar          (in)
	integer nvar		!total number of variables to write    (in)
	character*(*) type	!extension of file		       (in)
	double precision da_out(4)
	da_out = ia_out
	call open_scalar_file_d(da_out,nl,nvar,type)
	ia_out = nint(da_out)
	end

c*************************************************************

	subroutine open_scalar_file_d(da_out,nl,nvar,type)

c opens (NOS) file

c on return iu = -1 means that no file has been opened and is not written

	use mod_depth
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	double precision da_out(4)	!time information	       (in/out)
	integer nl		!vertical dimension of scalar          (in)
	integer nvar		!total number of variables to write    (in)
	character*(*) type	!extension of file		       (in)

	include 'femtime.h'

	include 'simul.h'

	integer nvers
	integer date,time
	integer iu,ierr
	character*80 title,femver

	integer ifemop
	double precision dgetpar

c-----------------------------------------------------
c open file
c-----------------------------------------------------

	iu = ifemop(type,'unformatted','new')
	if( iu .le. 0 ) goto 98
	da_out(4) = iu

c-----------------------------------------------------
c initialize parameters
c-----------------------------------------------------

	nvers = 5
	date = nint(dgetpar('date'))
	time = nint(dgetpar('time'))
	title = descrp
	call get_shyfem_version_and_commit(femver)

c-----------------------------------------------------
c write header of file
c-----------------------------------------------------

	call nos_init(iu,nvers)
	call nos_set_title(iu,title)
	call nos_set_date(iu,date,time)
	call nos_set_femver(iu,femver)
	call nos_write_header(iu,nkn,nel,nl,nvar,ierr)
        if(ierr.gt.0) goto 99
	call nos_write_header2(iu,ilhkv,hlv,hev,ierr)
        if(ierr.gt.0) goto 99

c-----------------------------------------------------
c write informational message to terminal
c-----------------------------------------------------

        write(6,*) 'open_scalar_file: ',type,' file opened ',it

c-----------------------------------------------------
c end of routine
c-----------------------------------------------------

	return
   98	continue
	write(6,*) 'error opening file with type ',type
	stop 'error stop open_scalar_file'
   99	continue
	write(6,*) 'error ',ierr,' writing file with type ',type
	stop 'error stop open_scalar_file'
	end

c*************************************************************

	subroutine write_scalar_file(ia_out,ivar,nlvddi,c)
	implicit none
	integer ia_out(4)	!time information
	integer ivar		!id of variable to be written
	integer nlvddi		!vertical dimension of c
	real c(nlvddi,*)		!scalar to write
	double precision da_out(4)
	da_out = ia_out
	call write_scalar_file_d(da_out,ivar,nlvddi,c)
	ia_out = nint(da_out)
	end

c*************************************************************

	subroutine write_scalar_file_d(da_out,ivar,nlvddi,c)

c writes NOS file
c
c the file must be open, the file will be written unconditionally

	use levels

	implicit none

	double precision da_out(4)	!time information
	integer ivar			!id of variable to be written
	integer nlvddi			!vertical dimension of c
	real c(nlvddi,*)		!scalar to write

	include 'femtime.h'

	logical binfo
	integer iu,ierr

	binfo = .false.

c-----------------------------------------------------
c check if files has to be written
c-----------------------------------------------------

	iu = nint(da_out(4))
	if( iu .le. 0 ) return

c-----------------------------------------------------
c write file
c-----------------------------------------------------

	call nos_write_record(iu,it,ivar,nlvddi,ilhkv,c,ierr)
	if(ierr.gt.0) goto 99

c-----------------------------------------------------
c write informational message to terminal
c-----------------------------------------------------

	if( binfo ) then
          write(6,*) 'write_scalar_file: ivar = ',ivar,' written at ',it
	end if

c-----------------------------------------------------
c end of routine
c-----------------------------------------------------

	return
   99	continue
	write(6,*) 'error ',ierr,' writing file at unit ',iu
	stop 'error stop write_scalar_file: error in writing record'
	end

c*************************************************************
c*************************************************************
c*************************************************************

        subroutine conmima(nlvddi,c,cmin,cmax)

c computes min/max for scalar field

	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

c arguments
	integer nlvddi		!vertical dimension of c
	real c(nlvddi,nkn)		!concentration (cconz,salt,temp,...)
        real cmin,cmax
c local
	integer k,l,lmax
	real cc
        logical debug
        integer kcmin,lcmin,kcmax,lcmax

        debug = .false.
        cmin = c(1,1)
        cmax = c(1,1)

	do k=1,nkn
	  lmax=ilhkv(k)
	  do l=1,lmax
	    cc = c(l,k)
            if( debug ) then
              if( cc .lt. cmin ) then
                    kcmin = k
                    lcmin = l
              end if
              if( cc .gt. cmax ) then
                    kcmax = k
                    lcmax = l
              end if
            end if
            cmin = min(cmin,cc)
            cmax = max(cmax,cc)
	  end do
	end do

        if( debug ) then
          write(6,*) 'conmima: ',kcmin,lcmin,cmin
          write(6,*) 'conmima: ',kcmax,lcmax,cmax
        end if

        end

c*************************************************************

        subroutine conmimas(nlvddi,c,cmin,cmax)

c computes min/max for scalar field -> writes some info

	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

c arguments
	integer nlvddi		!vertical dimension of c
	real c(nlvddi,nkn)		!concentration (cconz,salt,temp,...)
        real cmin,cmax
c common
	include 'femtime.h'
c local
	integer k,l,lmax
        integer ntot
	real cc
        logical debug
        integer kcmin,lcmin,kcmax,lcmax

        debug = .false.
        cmin = c(1,1)
        cmax = c(1,1)

        ntot = 0
	do k=1,nkn
	  lmax=ilhkv(k)
	  do l=1,lmax
	    cc = c(l,k)
            if( debug ) then
              if( cc .lt. cmin ) then
                    kcmin = k
                    lcmin = l
              end if
              if( cc .gt. cmax ) then
                    kcmax = k
                    lcmax = l
              end if
            end if
            cmin = min(cmin,cc)
            cmax = max(cmax,cc)
            if( cc .le. 0. ) then
                    ntot = ntot + 1
                    write(96,*) it,l,k,cc,ntot
            end if
	  end do
	end do

        if( ntot .gt. 0 ) then
                write(96,*) 'ntot: ',it,ntot
        end if

        if( debug ) then
          write(6,*) 'conmima: ',kcmin,lcmin,cmin
          write(6,*) 'conmima: ',kcmax,lcmax,cmax
        end if

        end

c*************************************************************

