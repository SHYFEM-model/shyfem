
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

c utility routines for operations on vectors
c
c contents :
c
c subroutine mima(xx,n,xmin,xmax)			min/max of vector
c subroutine mimar(xx,n,xmin,xmax,rnull)		min/max of vector
c subroutine mimari(xx,n,xmin,xmax,imin,imax,rnull)	min/max and index
c subroutine aver(xx,n,xaver,rnull)			aver of vector
c
c revision log :
c
c 26.08.1998	ggu	routines mimari transferred from newbcl0
c 31.05.1999	ggu	new comodity routine mima2i
c 23.03.2010	ggu	changed v6.1.1
c 19.06.2016	ggu	new routine aver, some other routines deleted
c 27.06.2016	ggu	changed VERS_7_5_16
c 09.09.2016	ggu	changed VERS_7_5_17
c 16.02.2019	ggu	changed VERS_7_5_60
c
c*******************************************
c
	subroutine mima(xx,n,xmin,xmax)
c
c computes min/max of vector
c
c xx		vector
c n		dimension of vector
c xmin,xmax	min/max value in vector
c
        implicit none
c
        integer n,i
        real xx(n)
        real xmin,xmax,x
c
	xmax=xx(1)
	xmin=xmax
c
	do i=1,n
          x=xx(i)
          if(x.gt.xmax) xmax=x
          if(x.lt.xmin) xmin=x
	end do
c
	return
	end
c
c*******************************************

        subroutine mimar(xx,n,xmin,xmax,rnull)

c computes min/max of vector (2d)
c
c xx            vector
c n             dimension of vector
c xmin,xmax     min/max value in vector
c rnull         invalid value

        implicit none

        integer n,i,nmin
        real xx(n)
        real xmin,xmax,x,rnull

        do i=1,n
          if(xx(i).ne.rnull) exit
        end do

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

c*******************************************

	subroutine mima3d(n,nlvddi,il,val,vmin,vmax)

! computes min/max of 3d vector

	implicit none

	integer n,nlvddi
	integer il(n)
	real val(nlvddi,n)
	real vmin,vmax

	integer k,l,lmax
	real v

	vmin = val(1,1)
	vmax = val(1,1)

	do k=1,n
	  lmax = il(k)
	  do l=1,lmax
	    v = val(l,k)
	    vmin = min(vmin,v)
	    vmax = max(vmax,v)
	  end do
	end do

	end

c*******************************************
c
        subroutine mimari(xx,n,xmin,xmax,imin,imax,rnull)
c
c computes min/max of vector
c
c xx            vector
c n             dimension of vector
c xmin,xmax     min/max value in vector
c imin,imax     pointer to min/max value in vector
c rnull         invalid value
c
        implicit none
c
        integer n,i,nmin
        integer imin,imax
        real xx(n)
        real xmin,xmax,x,rnull

        do i=1,n
          if(xx(i).ne.rnull) goto 1
        end do
    1   continue

        if(i.le.n) then
          xmax=xx(i)
          xmin=xx(i)
          imin=i
          imax=i
        else
          xmax=rnull
          xmin=rnull
          imin=0
          imax=0
        end if

        nmin=i+1

        do i=nmin,n
          x=xx(i)
          if(x.ne.rnull) then
            if(x.gt.xmax) then
                xmax=x
                imax=i
            end if
            if(x.lt.xmin) then
                xmin=x
                imin=i
            end if
          end if
        end do
c
        return
        end

c*******************************************
c*******************************************
c*******************************************

        subroutine aver(xx,n,xaver,rnull)

c computes aver of vector
c
c xx            vector
c n             dimension of vector
c xmin,xmax     min/max value in vector
c rnull         invalid value

        implicit none

        integer n
        real xx(n)
        real xaver,rnull

        integer i,nacu
        double precision acu

        nacu = 0
        acu = 0.
        xaver = rnull

        do i=1,n
          if(xx(i).ne.rnull) then
            acu = acu + xx(i)
            nacu = nacu + 1
          end if
        end do

        if( nacu .gt. 0 ) xaver = acu / nacu

        end

c*******************************************
c*******************************************
c*******************************************

