
!--------------------------------------------------------------------------
!
!    Copyright (C) 2018-2019  Georg Umgiesser
!    Copyright (C) 2018  Christian Ferrarin
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

! handles reading of particles and routines used also by shyplot
!
! revision log :
!
! 20.07.2018	ccf	written from scratch
! 25.10.2018	ggu	changed VERS_7_5_51
! 16.02.2019	ggu	changed VERS_7_5_60
! 25.06.2021	ggu	minor changes
! 21.03.2022	ggu	new routine lgr_rewind_block()
!
!**********************************************************
! check if file is lagrangian
!**********************************************************

        function lgr_is_lgr_file(file)

	implicit none

        logical lgr_is_lgr_file
        character*(*) file

        integer iu,ios
        integer ntype,nvers,ftype

        lgr_is_lgr_file = .false.

!--------------------------------------------------------------
! open file, read header and close
!--------------------------------------------------------------
        iu = 0
        open(iu,file=file,status='old',form='unformatted')

        read(iu,iostat=ios) ntype,nvers
        read(iu,iostat=ios) ftype

        if( ios /= 0 ) return

!--------------------------------------------------------------
! check for type and version
!--------------------------------------------------------------
        if( ntype /= 1617 .or. ftype /= 3 ) return

        lgr_is_lgr_file = .true.

        rewind(iu)

        end function lgr_is_lgr_file

!***********************************************************
! read lgr version header
!**********************************************************

        subroutine lgr_get_header(iu,nvers,lmax,ncust)

        implicit none

        integer, intent(in) 	:: iu
	integer, intent(out)	:: nvers
	integer, intent(out)	:: lmax
	integer, intent(out)	:: ncust

        integer mtype

        lmax = 1

        read(iu,end=98,err=99) mtype,nvers

        if( mtype /= 367265 ) goto 97
        if( nvers < 6 ) goto 95
        read(iu) lmax,ncust

        return
   95   continue
        write(6,*) mtype,nvers
        stop 'error stop lgr_get_header: cannot read this version'
   97   continue
        write(6,*) mtype,nvers
        stop 'error stop lgr_get_header: unknown mtype'
   98   continue
        stop 'error stop lgr_get_header: file is empty'
   99   continue
        stop 'error stop lgr_get_header: error reading header'

        end subroutine lgr_get_header

!**********************************************************
! peek header of the data block
!**********************************************************

        subroutine lgr_peek_block_header(iu,atime,nn,iwhat,ierr)

	implicit none

	integer			:: iu

	double precision	:: atime
	integer			:: nn
	integer			:: iwhat
	integer			:: ierr

        ierr = 1

        read(iu,iostat=ierr) atime,nn,iwhat
	if ( ierr > 0 ) goto 99
	if ( ierr < 0 ) return
        backspace(iu)

        return
   99   continue
        stop 'error stop lgr_peek_block_header: error reading record'

        end subroutine lgr_peek_block_header

!**********************************************************
! allocate variables 
!**********************************************************

        subroutine lgr_alloc(nn,nc
     +          ,id,ty,tt,s,ie,x,y,z,lb,hl,c)

        implicit none

        integer 			:: nn      !number of particles
	integer 			:: nc

	integer, allocatable		:: id(:)
	integer, allocatable		:: ty(:)
	double precision, allocatable	:: tt(:)
        real, allocatable 		:: s(:)
	integer, allocatable		:: ie(:)
        real, allocatable 		:: x(:),y(:),z(:)
	integer, allocatable		:: lb(:)
        real, allocatable 		:: hl(:)
        real, allocatable 		:: c(:,:)

	write(6,*) 'allocating lgr data structure: ',nn,nc

        if( allocated(id) ) then
          deallocate(id)
          deallocate(ty)
          deallocate(tt)
          deallocate(s)
          deallocate(ie)
          deallocate(x)
          deallocate(y)
          deallocate(z)
          deallocate(lb)
          deallocate(hl)
          deallocate(c)
        end if

        allocate(id(nn))
        allocate(ty(nn))
        allocate(tt(nn))
        allocate(s(nn))
        allocate(ie(nn))
        allocate(x(nn))
        allocate(y(nn))
        allocate(z(nn))
        allocate(lb(nn))
        allocate(hl(nn))
        allocate(c(nn,nc))

        end subroutine lgr_alloc

!**********************************************************
! return data block
!**********************************************************

        subroutine lgr_get_block(iu,nn,nc,
     +			idb,tyb,ttb,sb,ieb,xb,yb,zb,lbb,hlb,cb)

        implicit none

        integer 		:: iu
        integer 		:: nn      !number of particles
	integer 		:: nc	   !number of custom properties

	integer			:: idb(nn)
	integer			:: tyb(nn)
	double precision	:: ttb(nn)
        real 			:: sb(nn)
	integer			:: ieb(nn)
        real 			:: xb(nn),yb(nn),zb(nn)
	integer			:: lbb(nn)
        real 			:: hlb(nn)
        real 			:: cb(nn,nc)

	integer			:: id,ty,lb,ie
	real			:: x,y,z,hl,s
	double precision	:: tt,atime
	integer			:: iwhat
	real, allocatable	:: c(:)
	integer 		:: i,ii

!-----------------------------------------------------
! read in block header
!-----------------------------------------------------

        read(iu,end=100,err=99) atime,nn,iwhat

!-----------------------------------------------------
! loop over particles
!-----------------------------------------------------

	allocate(c(nc))

        do i=1,nn
	  read(iu) id,ty,tt,s,ie
          read(iu) x,y,z,lb,hl
          read(iu) (c(ii), ii=1,nc)
	  idb(i) = id
	  tyb(i) = ty
	  ttb(i) = tt
	  sb(i)  = s
	  ieb(i) = ie
	  xb(i)  = x
	  yb(i)  = y
	  zb(i)  = z
	  lbb(i) = lb
	  hlb(i) = hl
	  cb(i,:) = c(:)
        end do

        return
  100   continue
        return

   98   continue
        stop 'error stop lgr_get_xy_new: error reading data record'
   99   continue
        stop 'error stop lgr_get_xy_new: error reading time record'

        end subroutine lgr_get_block

!**********************************************************
! skip data block
!**********************************************************

        subroutine lgr_skip_block(iu,nn,nc)

        implicit none

        integer 		:: iu
        integer 		:: nn      !number of particles
	integer 		:: nc	   !number of custom properties

	integer			:: id,ty,lb,ie
	real			:: x,y,z,hl,s
	double precision	:: tt,atime
	integer			:: iwhat
	real, allocatable	:: c(:)
	integer 		:: i,ii

!-----------------------------------------------------
! read in block header
!-----------------------------------------------------

        read(iu,end=100,err=99) atime,nn,iwhat

!-----------------------------------------------------
! loop over particles
!-----------------------------------------------------

	allocate(c(nc))
        do i=1,nn
	  read(iu) id,ty,tt,s,ie
          read(iu) x,y,z,lb,hl
          read(iu) (c(ii), ii=1,nc)
        end do

        return
  100   continue
        return

   98   continue
        stop 'error stop lgr_get_xy_new: error reading data record'
   99   continue
        stop 'error stop lgr_get_xy_new: error reading time record'

        end subroutine lgr_skip_block

!**********************************************************
! rewind data block
!**********************************************************

        subroutine lgr_rewind_block(iu,nn)

        implicit none

        integer 		:: iu
        integer 		:: nn      !number of particles

	integer 		:: i

!-----------------------------------------------------
! rewind particles
!-----------------------------------------------------

        do i=1,nn
	  backspace(iu)
	  backspace(iu)
	  backspace(iu)
        end do

!-----------------------------------------------------
! rewind header
!-----------------------------------------------------

	backspace(iu)

        end subroutine lgr_rewind_block

!**********************************************************
! Get mean position per subroup of particles (type, assuming
! that type starts from 1).
!**********************************************************

        subroutine lgr_mean_posit(n,nc,tya,agea,sa,xa,ya,hla,ca,
     +				 nm,tym,agem,sm,xm,ym,hlm,cm)

        implicit none

	integer, intent(in)           :: n		!number of active parts
	integer, intent(in)           :: nc		!number of custom prop

        integer, intent(in)           :: tya(n)		!type
        real, intent(in)  	      :: agea(n)	!age [day]
        real, intent(in)              :: sa(n)		!settling
        real, intent(in)              :: xa(n),ya(n)	!x,y
        real, intent(in)              :: hla(n)		!absolute depth [m]
        real, intent(in)              :: ca(n,nc)	!custom props
 
	integer, intent(in)           :: nm		!number of types
        integer, intent(out)          :: tym(nm)
        real, intent(out) 	      :: agem(nm)
        real, intent(out)             :: sm(nm)
        real, intent(out)             :: xm(nm),ym(nm)
        real, intent(out)             :: hlm(nm)
        real, intent(out)             :: cm(nm,nc)

	real, allocatable	      :: npt(:)

        integer			      :: i,ty

	allocate(npt(nm))

	npt  = 0
        tym  = 0
        agem = 0.
        sm   = 0.
        xm   = 0.
        ym   = 0.
        hlm  = 0.
        cm   = 0.
 
        do i = 1,n
          ty = tya(i)
	  if ( ty <= 0 ) cycle
!          if ( hla(i) <= 4 .or. hla(i) >= 7 ) cycle
          npt(ty)  = npt(ty) + 1
          tym(ty)  = ty
          agem(ty) = agem(ty) + agea(i)
          sm(ty)   = sm(ty)   + sa(i)
          xm(ty)   = xm(ty)   + xa(i)
          ym(ty)   = ym(ty)   + ya(i)
          hlm(ty)  = hlm(ty)  + hla(i)
          cm(ty,:) = cm(ty,:) + ca(i,:)
	end do

        agem = agem / npt
        sm   = sm   / npt
        xm   = xm   / npt
        ym   = ym   / npt
        hlm  = hlm  / npt
	do i = 1,nc
          cm(:,i)  = cm(:,i)  / npt
	end do

	end subroutine lgr_mean_posit 

!**********************************************************
