
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

! revision log :
!
! 23.03.2010	ggu	changed v6.1.1
! 17.07.2015	ggu	changed VERS_7_1_53
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62

!**************************************************************************

	subroutine anneal(nkn,ngrddi,kphv,ng,iknot,maux,iaux)
c
c simulated annealing
c
	implicit none

	integer nkn,ngrddi
	integer kphv(nkn)
	integer ng(nkn),iknot(ngrddi,nkn)
	integer maux(nkn),iaux(nkn)

	integer nover,nlimit,itmax
	integer iseed,mbw,m,i
	integer itemp,iover,nsucc
	integer i1,i2,id,m1,m2,md
	integer j,k,kk,mh
	real t,tfactr

	integer mbandk,mband
	real ran
	logical metrop

	nover=100*nkn
	nlimit=10*nkn
	itmax=100
	tfactr=0.9
	t=0.2
	iseed=111

	do i=1,nkn
	  maux(i)=0
	end do

	mbw=0
	do j=1,nkn
	  k=kphv(j)
	  do i=1,ng(j)
	    kk=kphv(iknot(i,j))
	    mh=k-kk
	    if(mh.gt.0) maux(mh)=maux(mh)+1
	    if(mh.gt.mbw) mbw=mh
	  end do
	end do

	write(6,*) 'initial mbw : ',mbw

	mbw=mband(nkn,ngrddi,ng,iknot,kphv)

	write(6,*) 'initial mbw : ',mbw

	do itemp=1,itmax
	  nsucc=0
	  do iover=1,nover
	    i1=1+int(nkn*ran(iseed))
	    id=1+int(ng(i1)*ran(iseed))
	    i2=iknot(id,i1)
	    if(i1.lt.1.or.i1.gt.nkn) stop '1'
	    if(id.lt.1.or.id.gt.ng(i1)) stop '2'
	    m1=mbandk(nkn,ngrddi,i1,kphv(i2),ng,iknot,kphv)
	    m2=mbandk(nkn,ngrddi,i2,kphv(i1),ng,iknot,kphv)
	    m=max(m1,m2,abs(kphv(i2)-kphv(i1)))
	    md=m-mbw
	    if(metrop(float(md),t)) then
		nsucc=nsucc+1
		call exchi(nkn,ngrddi,i1,i2,kphv,ng,iknot,maux,mbw)
	    end if
c	    call mcontr(nkn,ngrddi,ng,iknot,kphv,maux,iaux,mbw)
c	    write(6,*) iover,i1,i2,md,mbw,nsucc
c	if(mod(iover,100).eq.0) read(5,'(i10)') i
	    if(nsucc.gt.nlimit) goto 1
	  end do
    1	  continue
	  write(6,*) '*** ',t,exp(-1./t),mbw,iover,nsucc
c	read(5,'(i10)') i
	  t=t*tfactr
	  if(nsucc.eq.0) return
	end do

	return
	end

c*****************************************************************

	function metrop(md,t)

c issues verdict for simulated annealing

	logical metrop
	real md,t

	real ran
	real aux,raux
	integer iseed
	save iseed
	data iseed / 97 /

c	metrop=.true.
c	if(md.lt.0.) return
c	aux=exp(-md/t)
c	raux=ran(iseed)
c	write(6,*) md,raux,aux
c	metrop = (md.lt.0) .or. (raux.lt.aux)

	metrop = .true.
	if (md.lt.0) return
	if (ran(iseed).lt.exp(-md/t)) return
	metrop = .false.

c	metrop = (md.lt.0) .or. (ran(iseed).lt.exp(-md/t))

	return
	end

c*****************************************************************

	subroutine exchi(nkn,ngrddi,i1,i2,kphv,ng,iknot,maux,mbw)

c exchanges indices

	implicit none

	integer nkn,ngrddi,mbw
	integer i1,i2
	integer kphv(nkn)
	integer ng(nkn)
	integer iknot(ngrddi,nkn)
	integer maux(nkn)

	integer k1,k2,k,i
	integer mold,mnew,mbwaux

	k1=kphv(i1)
	k2=kphv(i2)

	mbwaux=mbw

	do i=1,ng(i1)
	  k=kphv(iknot(i,i1))
	  mold=abs(k1-k)
	  mnew=abs(k2-k)
	  if(mold.ne.mnew.and.mnew.gt.0) then
	    maux(mold)=maux(mold)-1
	    maux(mnew)=maux(mnew)+1
	    if(mnew.gt.mbwaux) mbwaux=mnew
	  end if
	end do

	do i=1,ng(i2)
	  k=kphv(iknot(i,i2))
	  mold=abs(k2-k)
	  mnew=abs(k1-k)
	  if(mold.ne.mnew.and.mnew.gt.0) then
	    maux(mold)=maux(mold)-1
	    maux(mnew)=maux(mnew)+1
	    if(mnew.gt.mbwaux) mbwaux=mnew
	  end if
	end do

	if(mbwaux.gt.mbw) then
	  mbw=mbwaux
	else
	  do while(maux(mbw).eq.0)
	    mbw=mbw-1
	  end do
	end if

	kphv(i1)=k2
	kphv(i2)=k1

	return
	end

c*****************************************************************

	function mband(nkn,ngrddi,ng,iknot,kphv)

c computes bandwidth

	implicit none

	integer mband
	integer nkn,ngrddi
	integer ng(nkn)
	integer iknot(ngrddi,nkn)
	integer kphv(nkn)

	integer mbw,i,m
	integer mbandk

	mbw=0
	do i=1,nkn
	  m=mbandk(nkn,ngrddi,i,kphv(i),ng,iknot,kphv)
	  if(m.gt.mbw) mbw=m
	end do

	mband=mbw

	return
	end

c*****************************************************************

	function mbandk(nkn,ngrddi,j,k,ng,iknot,kphv)

c computes bandwidth for node at index j and node number k

c if real bandwidth is desired, call function with k=kphv(j)

	implicit none

	integer mbandk
	integer nkn,ngrddi
	integer j,k
	integer ng(nkn)
	integer iknot(ngrddi,nkn)
	integer kphv(nkn)

	integer m,mh,i,kk

	m=0
	do i=1,ng(j)
	  kk=kphv(iknot(i,j))
	  mh=abs(k-kk)
	  if(mh.gt.m) m=mh
	end do

	mbandk=m

	return
	end

c****************************************************************

	subroutine mcontr(nkn,ngrddi,ng,iknot,kphv,maux,iaux,mbw)

c controlls array maux if it is up to date

	implicit none

	integer nkn,ngrddi,mbw
	integer ng(nkn)
	integer iknot(ngrddi,nkn)
	integer kphv(nkn)
	integer maux(nkn),iaux(nkn)

	integer m,i,mbwaux
	integer j,k,kk,mh
	integer mbandk

	do i=1,nkn
	  iaux(i)=0
	end do

	mbwaux=0
	do j=1,nkn
	  k=kphv(j)
	  do i=1,ng(j)
	    kk=kphv(iknot(i,j))
	    mh=k-kk
	    if(mh.gt.0) iaux(mh)=iaux(mh)+1
	    if(mh.gt.mbwaux) mbwaux=mh
	  end do
	end do

	if(mbw.ne.mbwaux) then
	  write(6,*) 'mbw : ',mbw,mbwaux
	  mbw=mbwaux
	end if

	do i=1,nkn
	  if(maux(i).ne.iaux(i)) then
		write(6,*) 'mcontr : ',i,maux(i),iaux(i)
	  end if
	end do

	return
	end

c*********************************************************************

      FUNCTION RAN(ISEED)
      PARAMETER(IA=7141,IC=54773,IM=259200)
      ISEED=MOD(ISEED*IA+IC,IM)
      RAN=FLOAT(ISEED)/FLOAT(IM)
      RETURN
      END
