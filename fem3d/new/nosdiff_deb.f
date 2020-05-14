
!--------------------------------------------------------------------------
!
!    Copyright (C) 2010,2019  Georg Umgiesser
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

c 18.11.1998    ggu     check dimensions with dimnos
c 06.04.1999    ggu     some cosmetic changes
c 03.12.2001    ggu     some extra output -> place of min/max
c 09.12.2003    ggu     check for NaN introduced
c
c**************************************************************

! revision log :
!
! 23.03.2010	ggu	changed v6.1.1
! 16.02.2019	ggu	changed VERS_7_5_60

c**************************************************************

	program nosdiff

c computes difference between two NOS files - same grid
c
c still to be finished

	implicit none

	include 'param.h'

	real cv(nkndim)
	real cv3(nlvdim,nkndim)
	real cc(nlvdim,nkndim)
	real cv3_new(nlvdim,nkndim)
	real cv3diff(nlvdim,nkndim)
        integer iu,iunos
	integer ilhkv(nkndim)
	real hlv(nlvdim)
	real hev(neldim)

	integer nread,nin
	integer nvers
	integer nkn,nel,nvar
        integer nlvdi,nlv
        common /level/ nlvdi,nlv

	integer ierr
	integer it,ivar
	integer l,k
	character*80 title
	real rnull
	real cmin,cmax

        integer iapini,ifileo
	integer ifem_open_file
        character*80 file1,fileout, descrp
c--------------------------------------------------------------

	nread=0
	rnull=0.
	rnull=-1.

           do l=1,nlvdim
              do k=1,nkndim
                cv3(l,k) = 0. 
              enddo
           enddo
c--------------------------------------------------------------
c open basin and simulation
c--------------------------------------------------------------

	if(iapini(2,nkndim,neldim,0).eq.0) then
		stop 'error stop : iapini'
	end if

        nvers=3
        write(6,*) 'Enter name of first set file: '
        read (5,'(a)') file1
        nin=ifileo(55,file1,'unform','old')
        print*,nin
c       call rhnos(nin,nvers,nkndim,neldim,nlvdim,nkn,nel,nlv,nvar
c     +				,ilhkv,hlv,hev,title)
        call rfnos(nin,nvers,nkn,nel,nlv,nvar,title,ierr)
        call rsnos(nin,ilhkv,hlv,hev,ierr)
c        print*,nin,title,nvers,nvar,ierr

c open file
        iunos = 56
        fileout='ciccio.nos'
        
        iunos = ifileo(iunos,fileout,'unform','new')
	if( iunos .le. 0 ) goto 98
c write header of file


	call wfnos(iunos,3,nkn,nel,nlv,nvar,descrp,ierr)
        print*,iunos,nkn,nkndim,nel,neldim,nvar,descrp
        if(ierr.gt.0) goto 99
	call wsnos(iunos,ilhkv,hlv,hev,ierr)
        if(ierr.gt.0) goto 99

c--------------------------------------------------------------
c loop on data
c--------------------------------------------------------------

	do while(.true.)
           do l=1,nlv
              do k=1,nkn
                cv3_new(l,k) = cv3(l,k)
              enddo
           enddo

	   call rdnos(nin,it,ivar,nlvdim,ilhkv,cv3,ierr)
           if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
           if(ierr.ne.0) goto 100

           
	   !it = it1
	   !ivar = ivar1

	   nread=nread+1
	   write(6,*) 'time : ',it,'   ivar : ',ivar

	   do l=1,nlv
	     do k=1,nkn
	       cv3diff(l,k) =  cv3(l,k)- cv3_new(l,k)
           !if(cv3diff(l,k) .gt.0.01) print*,k,l,cv3diff(1,k),cv3(l,k)
             end do
	   end do

           call wrnos(iunos,it,ivar,nlv,ilhkv,cv3diff,ierr)

	end do	!do while

c--------------------------------------------------------------
c end of loop on data
c--------------------------------------------------------------

  100	continue

	write(6,*)
	write(6,*) nread,' records read'
	write(6,*)

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------
        return
   98	continue
	write(6,*) 'error opening file 1',file1
	stop 'error stop confop'
   99	continue
	write(6,*) 'error ',ierr,' writing file ',fileout
	stop 'error stop confop'
	end

c***************************************************************

        subroutine mimar(xx,n,xmin,xmax,rnull)

c computes min/max of vector
c
c xx            vector
c n             dimension of vector
c xmin,xmax     min/max value in vector
c rnull		invalid value

        implicit none

        integer n,i,nmin
        real xx(n)
        real xmin,xmax,x,rnull

	do i=1,n
	  if(xx(i).ne.rnull) goto 1
	end do
    1	continue

	if(i.le.n) then
	  xmax=xx(i)
	  xmin=xx(i)
	else
	  xmax=rnull
	  xmin=rnull
	end if

	nmin=i+1

        do i=nmin,n
          x=xx(i)
	  if(x.ne.rnull) then
            if(x.gt.xmax) xmax=x
            if(x.lt.xmin) xmin=x
	  end if
        end do

        end

c***************************************************************

        subroutine mimar_s(xx,nlvdim,n,xmin,xmax,rnull)

c computes min/max of vector
c
c xx            vector
c n             dimension of vector
c xmin,xmax     min/max value in vector
c rnull		invalid value

        implicit none

        integer n
        integer nlvdim
        real xx(nlvdim,n)
        real xmin,xmax,rnull

        integer k,l
        real x

        do k=1,n
          do l=1,nlvdim
            x=xx(l,k)
	    if(x.ne.rnull) then
              if( x .lt. xmin .or. x .gt. xmax ) then
                write(6,*) l,k,x
              end if
	    end if
          end do
        end do

        end

c***************************************************************

