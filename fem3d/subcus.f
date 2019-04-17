
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

c custom routines
c
c contents :
c
c subroutine custom( it )	custom routines
c
c subroutine totvol( it )	computes total volume in basin
c subroutine georg( it )	test some nodes and elements
c subroutine temp( it )		q-flux test
c subroutine bocche		adjust depth at inlets
c
c revision log :
c
c 22.01.1998	ggu	custom routine called at end of time step
c 20.03.1998	ggu	custom routine in own file to avoid compiler warnings
c 08.05.1998	ggu	new custom routine to check for mass balance
c 25.06.1998	ggu	test for q-flux (zerlina)
c 21.08.1998    ggu     xv eliminated
c 03.09.1998    ggu     subroutine bocche to adjust depth at Venice inlets
c 03.10.1998    ggu     subroutine salt to check salt
c 19.04.1999    ggu     subroutine impli to change weighting factor
c 27.05.1999    ggu     use icust to call custom routines
c 05.12.2001    ggu     fixed compiler error with -Wall -pedantic
c 10.08.2003    ggu     use accessor routines for chezy values in anpa()
c 14.08.2003    ggu     new routines test_hakata (icust=26) and node3d
c 14.08.2003    ggu     only 3D in concmass1
c 02.09.2003    ggu     new routine lago (76)
c 04.09.2003    ggu     some fixes in routine anpa
c 05.10.2004    ggu     new routine aldo to set conz in area
c 22.02.2005    ggu     subroutines deleted: salt
c 14.03.2005    ggu     subroutine traccia
c 30.06.2005    ggu     traccia changed, new routines jamal, sedimt
c 01.12.2005    ggu     more changes in traccia
c 23.03.2006    ggu     changed time step to real
c 18.10.2006    ggu     jamal has been updated and comented
c 22.02.2007    ggu     traccie routines updated
c 23.05.2007    ggu     new routines oscillation, kreis, debora
c 03.08.2007    ggu     in jamal reset introduced...
c 12.03.2008    ggu     jamal restructured -> computes projected res time
c 17.03.2008    ggu     new routines zinit, cprint
c 26.06.2008    ggu     routien for testing diffusion
c 09.10.2008    ggu     new call to confop
c 11.10.2008	ggu	pass zfranco into subroutine
c 12.11.2008	ggu	new routine joel
c 19.11.2008	ggu	new routine viscos
c 06.12.2008	ggu	new routine vdiffus()
c 28.01.2009	aac	new routine andreac()
c 24.03.2009	ggu	new zinit, tvd_test
c 23.05.2009	ggu	bug fix in jamal - reset also outer area
c 16.10.2009	ggu	some changes and documentation for traccia
c 18.11.2009	ggu	residence time computation in jamal with correction
c 12.02.2010	ggu	new routines for horizontal diffusion tests (diffus2d)
c 01.03.2010	ggu	new routines jamal_fra() to reset conz
c 15.12.2010	ggu	traccia routines moved to subtrace.f
c 26.01.2011	ggu	set rtauv for black sea (black_sea_nudge)
c 17.05.2011	ggu	new routines skadar_debug() and wet_dry()
c 12.07.2011	ggu	new routines init_ts()
c 09.03.2012	ggu	no call to jamal anymore
c 25.03.2014	ggu	new routine conz_decay_curonian()
c 21.10.2016	ccf	set ttauv and stauv for med-mar-bla
c 16.04.2019	ggu	small changes in close_inlets1() (bfill)
c
c******************************************************************

	subroutine custom( dtime )

c custom routines

	implicit none

	double precision dtime

	integer it
	real getpar

	integer, save :: icall = 0

	if( icall .lt. 0 ) return

	if( icall .eq. 0 ) then
	  icall = nint(getpar('icust'))
	  if( icall .le. 0 ) icall = -1
	  if( icall .lt. 0 ) return
	end if

	it = nint(dtime)

	if( icall .eq.  6 ) call impli(it)
	if( icall .eq.  7 ) call anpa(it)
	if( icall .eq.  9 ) call channel(it)
	if( icall .eq. 10 ) call sedcust(it)
	if( icall .eq. 11 ) call velchan(it)
	if( icall .eq. 15 ) call concmass(it)
	if( icall .eq. 16 ) call concmass1(it)
	if( icall .eq. 21 ) call hakata(it)
	!if( icall .eq. 25 ) call close_inlets
	if( icall .eq. 25 ) call close_inlets1
	if( icall .eq. 26 ) call test_hakata(it)
	if( icall .eq. 27 ) call traccia
	if( icall .eq. 28 ) call oscillation
	if( icall .eq. 29 ) call kreis
	if( icall .eq. 31 ) call zinit
	if( icall .eq. 32 ) call cprint(it)
	if( icall .eq. 76 ) call lago(it)
	if( icall .eq. 80 ) call aldo(it)
	if( icall .eq. 81 ) call jamal
	!if( icall .eq. 811 ) call jamal_fra
	if( icall .eq. 82 ) call sedimt
	if( icall .eq. 83 ) call joel
	if( icall .eq. 90 ) call diffus
	if( icall .eq. 91 ) call viscos
	if( icall .eq. 92 ) call vdiffus(1)
	if( icall .eq. 93 ) call vdiffus(2)
	if( icall .eq. 94 ) call diffus2d
        if( icall .eq. 95 ) call ggu_ginevra
	if( icall .eq. 101 ) call black_sea_nudge
        if( icall .eq. 102 ) call marmara_sea_nudge
	if( icall .eq. 110 ) call mpi_test_basin(1)
	if( icall .eq. 111 ) call mpi_test_basin(2)
	if( icall .eq. 112 ) call mpi_test_basin(3)
	if( icall .eq. 113 ) call mpi_test_basin(4)
	if( icall .eq. 201 ) call conz_decay_curonian
	if( icall .eq. 301 ) call cyano_diana
        if( icall .eq. 883 ) call debora(it)
        if( icall .eq. 884 ) call tsinitdebora(it)
        if( icall .eq. 888 ) call uv_bottom
        if( icall .eq. 601 ) call andreac
        if( icall .eq. 602 ) call tvd_test(it)
        if( icall .eq. 603 ) call wet_dry
        if( icall .eq. 604 ) call skadar_debug
        if( icall .eq. 700 ) call init_ts
	if( icall .eq. 710 ) call fluid_mud
	if( icall .eq. 901 ) call test_par
	if( icall .eq. 905 ) call lock_exchange

c	call totvol(it)
c	call georg(it)
c	call temp(it)

	end

c*****************************************************************

	subroutine totvol( it )

c computes total volume in basin

	use mod_geom_dynamic
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer it

	include 'param.h'

	integer iwet,ie
	real area,voldry,volwet,vol,a,v

	real areaele,volele

	iwet = 0
	area = 0.
	voldry  = 0.
	volwet = 0.

	do ie=1,nel

	  a = areaele(ie)
	  v = volele(ie,+1)

	  if( iwegv(ie) .eq. 0 ) then
	    iwet = iwet + 1
	    area = area + a
	    volwet = volwet + v
	  else
	    voldry = voldry + v
	  end if

	end do

	vol = voldry + volwet

	write(99,'(i10,i5,4e16.8)') 'totvol: '
     +			,it,iwet,area,volwet,voldry,vol

	end

c*****************************************************************

	subroutine georg( it )

c test some nodes and elements

	use mod_geom_dynamic
	use mod_hydro

	implicit none

	integer it

	integer ii
	integer ipint,ieint
	integer n,ibase

	include 'param.h'

	integer icall,k1,k2,ie1,ie2,ie3
	save icall,k1,k2,ie1,ie2,ie3
c	data icall,k1,k2,ie1,ie2,ie3 /0,556,555,777,788,788/
	data icall,k1,k2,ie1,ie2,ie3 /0,337,368,513,514,511/

	if( icall .eq. 0 ) then
	  write(97,*) k1,k2,ie1,ie2,ie3
	  k1 = ipint(k1)
	  k2 = ipint(k2)
	  ie1 = ieint(ie1)
	  ie2 = ieint(ie2)
	  ie3 = ieint(ie3)
	  write(97,*) k1,k2,ie1,ie2,ie3
	  icall = 1
	end if

	write(97,*) it,iwegv(ie1),iwegv(ie2),iwegv(ie3)
	write(97,*) znv(k1),znv(k2)
	write(97,*) (zenv(ii,ie1),ii=1,3)
	write(97,*) (zenv(ii,ie2),ii=1,3)
	write(97,*) (zenv(ii,ie3),ii=1,3)

	end

c*****************************************************************

	subroutine temp( it )

c q-flux test

	use mod_ts

	implicit none

	integer it

	include 'param.h'

	integer i
	real rit

	integer ndim
	parameter(ndim=5)


	integer nodes(ndim)
	save nodes
	data nodes /5,59,113,167,221/

	rit = it/3600.
	write(97,'(6f12.5)') rit,(tempv(1,nodes(i)),i=1,ndim)

	end

c*****************************************************************

	subroutine bocche

c adjust depth at inlets -> must be called from cst routines

	use mod_depth
	use basin

	implicit none

	integer ie,ii,itype
	real hm,h,hlido,hmala,hchio
	real getpar


	include 'param.h'

	write(6,*) 'adjusting depth at inlets...'

	hlido = getpar('hlido')
	hmala = getpar('hmala')
	hchio = getpar('hchio')

	do ie = 1,nel

	  hm = hev(ie)
	  itype = iarv(ie)

	  if( itype .eq. 3 ) then	!chioggia
	    h = min(hm,hchio)
	  else if( itype .eq. 4 ) then	!malamocco
	    h = min(hm,hmala)
	  else if( itype .eq. 5 .or. itype .eq. 9 ) then	!lido
	    h = min(hm,hlido)
	  else
	    h = hm
	  end if

	  if( h .ne. hm ) then
	    write(6,*) 'depth changed: ',ie,itype,hm,h
	  end if

	  call setdepele(ie,h)

	end do

	end

c*****************************************************************

	subroutine impli(it)

c changes implicity

	implicit none

	integer it

	integer it1,it2
	real r
	double precision d

	save d
	data d /1./

	it1 = 30000000
	it2 = 60000000

	if( it .lt. it1 ) then
	  if( d .ne. 1. ) then
	    d = 1.
	    r = d
	    call putpar('azpar',r)
	    call putpar('ampar',r)
	    write(6,*) 'weight changed : ',r
	  end if
	else if( it .gt. it2 ) then
	  if( d .ne. 0.5 ) then
	    d = 0.5
	    r = d
	    call putpar('azpar',r)
	    call putpar('ampar',r)
	    write(6,*) 'weight changed : ',r
	  end if
	else
	  d = it-it1
	  d = d / (it2-it1)
	  d = 1. - 0.5 * d
	  r = d
	  call putpar('azpar',r)
	  call putpar('ampar',r)
	  write(6,*) 'weight changed : ',r
	end if

	end

c*****************************************************************

	subroutine anpa(it)

c changes chezy

	implicit none

	integer it

	integer itmax,itdt
	integer i
	real czmin,czact,czmax
	real czin,czout
	real czp(0:100),czm(0:100)
	save czp,czm
	integer nareas
	save nareas

	integer icall
	save icall
	data icall / 0 /

	itmax = 600000
	itdt = 7200
	czmin = 10.
	czmin = 5.
	czmin = 4.
	czmin = 3.

	if( it .lt. 0 ) return
	if( it .gt. itmax ) return

	if( icall .eq. 0 ) then
	  icall = 1

	  call n_chezy_values(nareas)
	  
	  do i=0,nareas
	    call get_chezy_values(i,czp(i),czm(i))
	  end do

	  write(6,*) 'anpa initialized...  ',nareas
	  write(6,*) (czp(i),i=0,20)
	  write(6,*) (czm(i),i=0,20)
	end if

	if( mod(it,itdt) .ne. 0 ) return

	do i=0,nareas
	  call get_chezy_values(i,czp(i),czm(i))
	end do

	write(15,'(i10,12i5)') it,(nint(czp(i)),i=0,nareas)
	write(15,'(i10,12i5)') it,(nint(czm(i)),i=0,nareas)

	do i=0,nareas
	  if( ( i .ge. 3 .and. i .le. 5 ) .or. i .eq. 9 ) then
	    call get_chezy_values(i,czin,czout)
	    if( czin  .gt. czmin ) czin  = czin  - 1.
	    if( czout .gt. czmin ) czout = czout - 1.
	    call set_chezy_values(i,czin,czout)
	  end if
	end do

	end

c*****************************************************************

	subroutine channel(it)

c channel output

	use mod_ts
	use mod_hydro_vel
	use levels

	implicit none

	include 'param.h'

	integer it

 
        real num(0:nlvdi)
        real nuh(0:nlvdi)
        real tk(0:nlvdi)
        real ep(0:nlvdi)
        real rl(0:nlvdi)

	logical berror
	integer i,k,l,nlev,nlev1
	integer ipext

	integer ndim
	parameter( ndim = 3 )
	integer nodes(ndim)
	save nodes
	integer icall
	save icall

	data nodes / 25 , 28 , 31 /
	data icall / 0 /

	if( icall .eq. 0 ) then
	  icall = 1
	  call n2int(ndim,nodes,berror)
	  if( berror ) stop 'error stop channel: n2int'
	  write(6,*) 'total nodes in channel: ',ndim
	  do i=1,ndim
	    write(6,*) i,nodes(i),ipext(nodes(i))
	  end do
	  open(81,file='gvel.dat',status='unknown',form='formatted')
	  open(82,file='grho.dat',status='unknown',form='formatted')
	  open(89,file='gnum.dat',status='unknown',form='formatted')
	  open(84,file='gnuh.dat',status='unknown',form='formatted')
	  open(85,file='gtke.dat',status='unknown',form='formatted')
	  open(86,file='geps.dat',status='unknown',form='formatted')
	  open(87,file='glen.dat',status='unknown',form='formatted')
	end if
	
c	write(86,*) it,ndim
c	write(88,*) it,ndim
c	do i=1,ndim
c	  k = nodes(i)
c	  nlev = ilhkv(k)
c
c	  write(86,*) i,k,nlev
c	  do l=1,nlev
c	    write(86,*) l,ulnv(l,k),vlnv(l,k)
c	  end do
c
c	  write(88,*) i,k,nlev
c	  do l=1,nlev-1
c	    write(88,*) l,visv(l,k),difv(l,k),tken(l,k),eps(l,k),rls(l,k)
c	  end do
c
c	end do

	k = nodes(2)
	nlev = ilhkv(k)
	nlev1 = nlev - 1

	call gotm_get(k,nlev,num,nuh,tk,ep,rl)

	write(81,*) it,nlev,1.
	write(81,*) (ulnv(l,k),l=1,nlev)

	write(82,*) it,nlev,1.
	write(82,*) (rhov(l,k),l=1,nlev)

	write(89,*) it,nlev1,1.
	write(89,*) (num(l),l=1,nlev1)

	write(84,*) it,nlev1,1.
	write(84,*) (nuh(l),l=1,nlev1)

	write(85,*) it,nlev1,1.
	write(85,*) (tk(l),l=1,nlev1)

	write(86,*) it,nlev1,1.
	write(86,*) (ep(l),l=1,nlev1)

	write(87,*) it,nlev1,1.
	write(87,*) (rl(l),l=1,nlev1)

	end

c*****************************************************************

	subroutine sedcust(it)

c channel output

	use mod_hydro_baro

	implicit none

	include 'param.h'

	integer it


	integer ie,i
	integer iunit
	real u,v,h
	integer ieint
	real depele

	integer ndim
	parameter(ndim=5)
	integer iesp,iisp
	save iesp,iisp
	integer iespv(ndim),iispv(ndim)
	save iespv,iispv
	integer icall
	save icall

	data iesp / 7106 /
	data iespv / 100 , 6941 , 1989 , 4284 , 6237 /
	data icall / 0 /

	if( icall .eq. 0 ) then
	  iisp = ieint(iesp)
c	  write(78,*) 'element: ',iisp,iesp
	  do i=1,ndim
	    iispv(i) = ieint(iespv(i))
	  end do
	  icall = 1
	end if

	ie = iisp
	u = unv(ie)
	v = vnv(ie)
	h = depele(ie,+1)
	write(78,*) it,u,v,h

	do i=1,ndim
	  ie = iispv(i)
	  u = unv(ie)
	  v = vnv(ie)
	  h = depele(ie,+1)
	  iunit = 80 + i
	  write(iunit,*) it,u,v,h
	end do

	end

c*****************************************************************

	subroutine velchan(it)

c channel velocity output

	use mod_ts
	use mod_hydro_print
	use mod_hydro
	use levels
	use basin

	implicit none

	include 'param.h'

	integer it



	logical berror
	logical b1,b2,b3,b4
	integer k,i,l,nlev
	integer ivert

	integer ndim,ndim1
	parameter(ndim=5,ndim1=1)

	integer k1(ndim),k2(ndim),k3(ndim1),k4(ndim1)
	save k1,k2,k3,k4

	integer icall
	save icall

	data icall / 0 /
	data k1 / 45,34,23,12,1 /
	!data k2 / 11,22,33,44,55 /
	data k2 / 203 , 228 , 253 , 278 , 303 /
	!data k3 / 28 /
	data k3 / 253 /
	data k4 / 558 /

	b1 = .false.
	b2 = .true.
	b3 = .true.
	b4 = .false.

	if( icall .eq. 0 ) then
	  call n2int(ndim,k1,berror)
	  if( berror ) stop 'error stop velchan: k1'
	  call n2int(ndim,k2,berror)
	  if( berror ) stop 'error stop velchan: k2'
	  call n2int(ndim1,k3,berror)
	  if( berror ) stop 'error stop velchan: k3'
	  icall = 1

	  if( b3 ) then
	    k = k3(1)
	    nlev = ilhkv(k)
	    ivert = 1
	    call animhead(94,nlev,ivert,0.,float(nlev),-2.,2.)	!vel
	    call animhead(96,nlev,ivert,0.,float(nlev),-5.,35.)	!temp
	  end if
	end if

	call nantest(nkn*nlvdi,uprv,'uprv')
	call nantest(nkn*nlvdi,vprv,'vprv')
	call nantest(nkn*nlvdi,tempv,'tempv')

	if( b1 ) then
	write(92,*) it,ndim
	do i=1,ndim
	  k = k1(i)
	  nlev = ilhkv(k)
	  write(92,*) k,nlev,znv(k)
	  write(92,*) (uprv(l,k),l=1,nlev)
	  write(92,*) (tempv(l,k),l=1,nlev)
	end do
	end if

	if( b2 ) then
	write(93,*) it,ndim
	do i=1,ndim
	  k = k2(i)
	  nlev = ilhkv(k)
	  write(93,*) k,nlev,znv(k)
	  write(93,*) (uprv(l,k),l=1,nlev)
	  write(93,*) (tempv(l,k),l=1,nlev)
	end do
	end if

	if( b3 ) then
	k = k3(1)
	nlev = ilhkv(k)
	call animdata(94,it,nlev,uprv(1,k))
	call animdata(96,it,nlev,tempv(1,k))
	!write(94,*) it,nlev
	!write(94,*) (uprv(l,k),l=1,nlev)
	!write(96,*) it,nlev,1.
	!write(96,*) (tempv(l,k),l=1,nlev)
	end if

	if( b4 ) then
	k = k4(1)
	nlev = ilhkv(k)
	write(94,*) it,nlev,1.
	write(94,*) (uprv(l,k),l=1,nlev)
	write(96,*) it,nlev,1.
	write(96,*) (tempv(l,k),l=1,nlev)
	end if

	end

c*****************************************************************

	subroutine concmass(it)

c set up zeta and conc for mass conservation test

	use mod_conz
	use mod_hydro
	use levels
	use basin

	implicit none

	include 'param.h'

	integer it



	integer k,ie,ii,n,ip,l
	integer ntot,nlev
	integer ibase,iw
	real y0,dy,dz,y,z
	real conz
	real res

	integer icall
	save icall
	integer kn(5)
	save kn
	data icall / 0 /
	data kn / 3186,2729,2507,2371,2180 /

	if( icall .eq. 0 ) then

	y0 = 28000.
	dy = 25000.
	dz = 0.5
	ntot = 0

	do k=1,nkn
	  y = ygv(k)
	  z = dz*(y-y0)/dy
	  znv(k) = z
	  nlev = ilhkv(k)
	  if( z .gt. 0. ) then
	    conz = 100.
	    ntot = ntot + 1
	  else
	    conz = 0.
	  end if
	  do l=1,nlev
	    cnv(l,k) = conz
	  end do
	end do

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    zenv(ii,ie) = znv(k)
	  end do
	end do

	call setweg(0,iw)

	icall = 1

	write(6,*) 'conz and level initialized... ',ntot,nkn

	end if

	call massconc(1,cnv,nlvdi,res)
	write(6,*) 'total dissolved mass: ',res
	write(6,*) 'values: ',(cnv(1,kn(ii)),ii=1,5)

	end

c*****************************************************************

	subroutine concmass1(it)

c mass conservation test

	use mod_conz
	use mod_geom_dynamic
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer it

	include 'param.h'


	real res
	integer i,ie,iweg

	iweg = 0
	do ie=1,nel
	  if(iwegv(ie) .ne. 0 ) iweg = iweg + 1
	end do

	call massconc(1,cnv,nlvdi,res)		!3D

	write(6,*) 'total dissolved mass: ',it,res,iweg
c	call debug_node(4179)
c	call debug_dry

	end

c*****************************************************************

	subroutine hakata(it)

	implicit none

        integer it

        integer ifemopa

        integer iunit
        save iunit
        data iunit / 0 /

        if( iunit .le. 0 ) then
          iunit = ifemopa('Opening Hakata file','.hak','f','new')
        end if

        call hakanode(iunit,it,1608)    !901
c        call hakanode(iunit,it,777)     !178
        call hakanode(iunit,it,585)     !1032 C-10
        call hakanode(iunit,it,383)     !231

	end

c*****************************************************************

	subroutine hakanode(iunit,it,k)

	use mod_ts
	use mod_diff_visc_fric
	use mod_hydro_print
	use levels

	implicit none

        integer iunit
        integer it
        integer k

        include 'param.h'


        integer nlev

        nlev = ilhkv(k)

        write(iunit,1000) k,'temp'
	call animdata(iunit,it,nlev,tempv(1,k))
        write(iunit,1000) k,'salt'
	call animdata(iunit,it,nlev,saltv(1,k))
        write(iunit,1000) k,'u'
	call animdata(iunit,it,nlev,uprv(1,k))
        write(iunit,1000) k,'v'
	call animdata(iunit,it,nlev,vprv(1,k))
        write(iunit,1000) k,'vis'
	call animdata(iunit,it,nlev+1,visv(0,k))
        write(iunit,1000) k,'dif'
	call animdata(iunit,it,nlev+1,difv(0,k))

        return
 1000   format(i8,2x,a)
	end

c*****************************************************************
c*****************************************************************

	subroutine animhead(iunit,lmax,ivert,xmin,xmax,ymin,ymax)

	implicit none

	integer iunit
	integer lmax,ivert
	real xmin,xmax,ymin,ymax

	integer magic,nvers
	parameter( magic = 46728645 , nvers = 1 )

	write(iunit,*) magic,nvers
	write(iunit,*) lmax,ivert
	write(iunit,*) xmin,xmax,ymin,ymax

	end

c*****************************************************************

	subroutine animdata(iunit,it,n,val)

	implicit none

	integer iunit
	integer it,n
	real val(n)

	integer i

	write(iunit,*) it,n
	write(iunit,*) (val(i),i=1,n)

	end

c*****************************************************************

	subroutine close_inlets

	implicit none

	include 'femtime.h'
c local
        integer ie,ii,k
c	integer ibc,nbc
	integer ibc,nbc
	integer ibctot
	real dz,dzcmh
	real zdate
	real dt
	integer nbnds,itybnd
	integer ifemopa

	integer iunit
	save iunit
	data iunit / 0 /

	if( iunit .eq. 0 ) then
	  iunit = ifemopa('','.ccc','form','new')
	end if

	nbc = nbnds()
	dz = 0.
	dzcmh = 0.
	call get_timestep(dt)

	ibctot=0
	do ibc=1,nbc
	  if( itybnd(ibc) .lt. 0 ) ibctot = ibctot + 1
	end do

	if( ibctot .eq. nbc ) then	!lagoon closed

	    dzcmh = 2.
	    dz = dt * dzcmh / ( 100. * 3600. )

	    call raise_zeta(dz)

	end if

	if( it .eq. 86400 ) call set_zdate(0,0.40)
	call get_zdate(0,zdate)
	write(iunit,*) it,dzcmh,zdate,(itybnd(ibc),ibc=1,nbc)
	
	end

c*****************************************************************

	subroutine close_inlets1

	use mod_hydro
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'
	include 'femtime.h'
c local
        integer ie,ii,k
c	integer ibc,nbc
	integer ibc
	integer ibctot
	real dz,dzmmh
	real zdate
	real dt
	integer nbnds,itybnd
	integer ifemopa
	real getpar

	logical bfill
	integer ndim
	parameter(ndim=200)
	real vals(ndim)
	integer icl,nbcc,iclose,n,ih,iclass,itcl
	integer itotal,ittot,iclose_old
	integer nbc
	real psv,dsv,zrise,zextra,zsalv,zbound,zfranco
	real zmax,wmax,rmax
	real windv,rainv,zlevel
	real t
	character*4 class
	real zvbnds

	save itotal,ittot,iclose_old

	integer iunit,i66
	save iunit,i66
	data iunit / 0 /

        integer icm
        icm(t) = nint(100.*t)

	bfill = .false.
	bfill = .true.		!fill lagoon with rain if closed

	nbc = nbnds()

	if( iunit .eq. 0 ) then
	  iunit = ifemopa('','.ccc','form','new')
	  i66 = ifemopa('','.cc1','form','new')
	  itotal = 0
	  ittot = 0
	  iclose_old = 0
	end if

c dsl  pts  btz  cfg  pbo  vgr  tso  fus  pov  ser  tre  gbo  vdg  lsl
c 2750 2449 19   919  533  1385 936  2174 2032 3140 3342 4086 4343 3971

c 1336 1717

	k = 1336		!2750 extern
	dsv = znv(k)
	k = 1717		!2449 extern
	psv = znv(k)

	zlevel = psv
	call get_timestep(dt)

	zrise   = 0.01*getpar('zrise')
	zfranco = 0.01*getpar('zfranc')
	zsalv   = 0.01*getpar('zsalv')
        zextra  = zsalv - 1.

c iclose=1      => lagoon is completely closed

	call is_closed(0,nbcc,icl)      !nbcc is total number of OBC
	iclose = 0
	if( icl .eq. nbcc ) iclose = 1

          call get_prev(ndim,it,zfranco,n,vals)
          call get_wind_rain(it,windv,rainv,wmax,rmax,dzmmh)

          call class12new(it,n,vals,zsalv,zextra,zrise,wmax,rmax,psv
     +                  ,zmax,zdate,ih,iclass,class)

          call set_zdate(0,zdate)

          if( iclose .eq. 1 ) then      !already closed
            dz = dzmmh*dt/(1000.*3600.)   !convert to [m/timestep]
            if( bfill ) call raise_zeta(dz)
          end if

	  zlevel = 0
	  if( ih .gt. 0 ) zlevel = zrise+vals(ih)       !level at time ih
	  zbound = zvbnds(1)                            !level outside
          write(i66,2300) it,iclass,class,icl,ih
     +          ,icm(zmax),icm(zlevel),icm(psv),icm(dsv)
     +		,icm(zbound),icm(zdate),(itybnd(ibc),ibc=1,nbc)

c ittot         number of time steps inlets are completely closed
c itotal        number of total closures

          if( iclose .eq. 1 ) then
            ittot = ittot + 1
	    if( iclose_old .eq. 0 ) itotal = itotal + 1
          end if

	iclose_old = iclose

	write(iunit,*) it,dzmmh,zdate
     +			,(itybnd(ibc),ibc=1,nbc),itotal,ittot
	
 2200       format(20i4)
 2100     format(i10,i2,a5,i2,i3,4f6.2)
 2300     format(i10,i2,a5,i2,i3,3x,10i4)
 2301     format(i10,i2,a5,i2,i3,3x,10i4,2x,3i4)
 2000   format(i10,i4,7f9.2)
 2500   format(a,a4,i10,2f9.2,i10)

	end

c*****************************************************************

	subroutine test_hakata(it)

	implicit none

        integer it

	logical berror
	integer i,k
        integer ifemopa

	integer ndim
	parameter(ndim=4)
	integer nodes(ndim)
        integer iunit
        save iunit
        save nodes
        data iunit / 0 /
        data nodes / 901, 178, 1032, 231 /

        if( iunit .le. 0 ) then
          iunit = ifemopa('Opening out3d','.o3d','f','new')
	  write(6,*) 'Nodes from test_hakata: '
          call n2int(ndim,nodes,berror)
	  do i=1,ndim
	    write(6,*) i,nodes(i)
	  end do
	  if( berror ) stop 'error stop test_hakata'
        end if

	do i=1,ndim
	  k = nodes(i)
	  call node3d(iunit,it,k)
	end do

c        call hakanode(iunit,it,1608)    !901
c        call hakanode(iunit,it,777)     !178
c        call hakanode(iunit,it,585)     !1032 C-10
c        call hakanode(iunit,it,383)     !231

	end

c*****************************************************************

	subroutine node3d(iunit,it,k)

	use mod_ts
	use mod_hydro_print
	use mod_hydro_vel
	use mod_hydro
	use levels

	implicit none

        integer iunit
        integer it
        integer k

        include 'param.h'


        integer nlev
	integer l

        nlev = ilhkv(k)

	write(iunit,1000) 'header (z): ',it,k
	write(iunit,*) 1,znv(k)

	write(iunit,1000) 'header (u): ',it,k
	write(iunit,*) nlev,(uprv(l,k),l=1,nlev)
	write(iunit,1000) 'header (v): ',it,k
	write(iunit,*) nlev,(vprv(l,k),l=1,nlev)
	write(iunit,1000) 'header (w): ',it,k
	write(iunit,*) nlev+1,(wlnv(l,k),l=0,nlev)

	write(iunit,1000) 'header (temp): ',it,k
	write(iunit,*) nlev,(tempv(l,k),l=1,nlev)
	write(iunit,1000) 'header (salt): ',it,k
	write(iunit,*) nlev,(saltv(l,k),l=1,nlev)

        return
 1000   format(a,2i10)
	end

c*****************************************************************

	subroutine lago(it)

	use mod_hydro_print
	use mod_hydro_vel
	use levels

	implicit none

	integer it

	include 'param.h'


	logical berror
	integer i,k
	real u,v
 
	integer ndim
	parameter (ndim=3)
	real vel(ndim)
	integer icall
	integer nodes(ndim)
	save nodes,icall
	data nodes /1020,2021,2527/
	data icall /0/

	if( icall .eq. 0 ) then
	  icall = 1
	  call n2int(ndim,nodes,berror)
	  if( berror ) stop 'error stop'
	end if

	do i=1,ndim
	  k = nodes(i)
	  u=uprv(1,k)
	  v=vprv(1,k)
	  !write(8,*) i,k,u,v
	  vel(i) = sqrt(u*u+v*v)
	end do
	  
	write(8,*) it,(vel(i),i=1,ndim)
	write(9,*) (vel(i),i=1,ndim)

	end

c*****************************************************************

	subroutine aldo(it)

	use mod_conz
	use levels
	use basin

	implicit none

	integer it

	include 'param.h'


	integer l,k,lmax
        real x0,y0,x1,y1,x2,y2,x3,y3
        real xp,yp
        logical is_left

        x0 = 2427.180420
        y0 = 2443.140137
        x1 = 4201.409180
        y1 = 235.025009
        x2 = 4201.409180
        y2 = 235.025009
        x3 = 5664.500000
        y3 = 1452.153198

        if( it .eq. 86400-300 ) then
         do k=1,nkn
          xp = xgv(k)
          yp = ygv(k)
          if( .not. is_left(xp,yp,x0,y0,x1,y1) .and.
     +                  is_left(xp,yp,x2,y2,x3,y3) ) then
            lmax = ilhkv(k)
            do l=1,lmax
              cnv(l,k) = 100.
            end do
          end if
         end do
        end if

        end

c*****************************************************************

        function is_left(xp,yp,x0,y0,x1,y1)

        implicit none

        logical is_left
        real xp,yp,x0,y0,x1,y1
        real dx,dy,dxn,dyn,dxp,dyp
        real scal

        dx = x1-x0
        dy = y1-y0
        dxn = -dy
        dyn = dx

        dxp = xp - x0
        dyp = yp - y0

        scal = dxp*dxn + dyp*dyn

        is_left = scal .gt. 0.

        end

c*****************************************************************
c*****************************************************************
c*****************************************************************
c*****************************************************************
c*****************************************************************

        subroutine jamal

c computes residence time online - one value for whole lagoon

	use mod_conz, only : cnv
	use levels
	use basin

        implicit none

        include 'param.h'

	include 'femtime.h'


        integer ie,ii,k,lmax,l,ia
        logical bnoret,breset,bstir
        real vol,conz,perc,percmin
        double precision mass,volume
        double precision massaux,volaux
        real volnode
	integer ifemop
	integer iaout,itmin,idtreset
	integer secsmonth
	integer iconz
	real dt,c0
	real restime,restime1,restimec
	real remnant,rlast
	real resmed,resstd
	real v1v(nkn)

	real getpar

	integer iu,it0,ndata
	save iu,it0,ndata
	double precision remint,remlog,remtim
	save remint,remlog,remtim
	double precision rsum,rsumsq
	save rsum,rsumsq
        double precision mass0
        save mass0

        integer icall
        save icall
        data icall / 0 /

c------------------------------------------------------------
c parameters
c------------------------------------------------------------
c
c bnoret	true if no return flow is used (conzentrations outside
c		are explicitly set to 0)
c bstir		simulates completely stirred tank
c		(replaces at every time step conz with average conz)
c percmin	percentage to reach -> after this stop computation
c		use 0 if no premature end is desired
c iaout		area code of elements out of lagoon (used for init and retflow)
c		use -1 to if no outside areas exist
c c0		initial concentration of tracer
c itmin		time from when to compute residence time
c idtreset	time step to reset concentration to c0
c		use 0 if no reset is wanted
c
c secsmonth	how many seconds in a month
c--------------------------
c default settings
c--------------------------
        bnoret = .false.
	bstir = .false.
	percmin = 0.
	iaout = -1
	c0 = 1.
	itmin = 0
	idtreset = 0
c--------------------------
c set how many seconds are in a month
c--------------------------
	secsmonth = nint(30.5 * 86400)		! 366 days per year
	!secsmonth = secsmonth - 7200		! 365 days per year
c--------------------------
c nador
c--------------------------
c	bnoret = .true.
c	percmin = 1.
c	iaout = 0
c--------------------------
c alimini
c--------------------------
c	idtreset = 1 * secsmonth
c--------------------------
c mar menor
c--------------------------
c 	just use default settings
c--------------------------
c taranto
c--------------------------
c	itmin = -1
c	itmin = 0
c	idtreset = 3 * secsmonth
c--------------------------
c curonian lagoon
c--------------------------
	secsmonth = secsmonth - 7200		! 365 days per year
	c0 = 0.
	idtreset = 3 * secsmonth
	idtreset = 86400

c------------------------------------------------------------
c do not change anything after this point
c------------------------------------------------------------

	if( itmin .eq. -1 ) itmin = itanf

c------------------------------------------------------------
c can we run the routine?
c------------------------------------------------------------

	iconz = nint(getpar('iconz'))
	if( iconz <= 0 ) then
	  write(6,*) 'iconz = ',iconz
	  stop 'error stop jamal: iconz must be 1'
	end if

c------------------------------------------------------------
c is it time to run the routine?
c------------------------------------------------------------

        if( it .le. itmin ) return

c------------------------------------------------------------
c initialization -> decide on reset
c------------------------------------------------------------

	breset = .false.		!normally do not reset

        if( icall .eq. 0 ) then
          write(6,*) 'initialization of routine jamal'
	  iu = ifemop('.jam','formatted','new')
	  breset = .true.		!always reset at first call
	  it0 = it
        end if

	if( idtreset .gt. 0 ) then
	  if( it-it0 .ge. idtreset ) breset = .true.
	end if

c------------------------------------------------------------
c flag nodes that are inside lagoon (v1v(k)=1)
c------------------------------------------------------------

        do k=1,nkn
          v1v(k) = 0.
        end do

        do ie=1,nel
          ia = iarv(ie)
          if( ia .ne. iaout ) then
              do ii=1,3
                k = nen3v(ii,ie)
                v1v(k) = 1.
              end do
          end if
        end do

c------------------------------------------------------------
c reset concentrations
c------------------------------------------------------------

        if( breset ) then		!reset concentrations to c0
          write(6,*) 'resetting concentrations in jamal at time ',it
          do k=1,nkn
            lmax = ilhkv(k)
	    conz = 0.
            if( v1v(k) .ne. 0. ) conz = c0
            do l=1,lmax
              cnv(l,k) = conz
            end do
          end do
	end if

c------------------------------------------------------------
c total mass (only for nodes inside lagoon)
c------------------------------------------------------------

        mass = 0.
        volume = 0.
        do k=1,nkn
          if( v1v(k) .ne. 0. ) then
            lmax = ilhkv(k)
            do l=1,lmax
              vol = volnode(l,k,+1)
              conz = cnv(l,k)
              mass = mass + vol*conz
              volume = volume + vol
            end do
          end if
        end do

c------------------------------------------------------------
c stirred tank?
c------------------------------------------------------------

	if( bstir ) then
	  conz = mass / volume
          do k=1,nkn
            if( v1v(k) .ne. 0. ) then
              lmax = ilhkv(k)
              do l=1,lmax
                cnv(l,k) = conz
	      end do
	    end if
	  end do
	end if

        massaux = 0.
        volaux = 0.
        do k=1,nkn
          if( v1v(k) .ne. 0. ) then
            lmax = ilhkv(k)
            do l=1,lmax
              vol = volnode(l,k,+1)
              conz = cnv(l,k)
              massaux = massaux + vol*conz
              volaux = volaux + vol
            end do
          end if
        end do

	!write(67,*) it,mass,massaux,volume,volaux
	!write(67,*) it,mass,massaux

c------------------------------------------------------------
c reset variables to compute residence time
c------------------------------------------------------------

        if( breset ) then
	  mass0 = mass 
	  remint = 0.
	  remlog = 0.
	  remtim = 0.
	  it0 = it
	  ndata = 0
	  rsum = 0.
	  rsumsq = 0.
	end if

c------------------------------------------------------------
c write to file
c------------------------------------------------------------

	call get_timestep(dt)

	remnant = 0.
        if( mass0 .gt. 0. ) remnant = mass/mass0
        perc = 100.*remnant

	remint = remint + remnant*dt	!integrated remnant function
	restime = remint/86400.		!residence time in days

	rlast = remnant
	if( rlast .ge. 1. ) rlast = 0.
	restimec = restime/(1.-rlast)	!corrected residence time

	if( remnant .gt. 0. ) remlog = remlog - log(remnant)
	remtim = remtim + (it-it0)
	restime1 = 0.
	if( remlog .gt. 0. ) restime1 = ( remtim / remlog ) / 86400.

	ndata = ndata + 1
	rsum = rsum + restime1
	rsumsq = rsumsq + restime1*restime1
	resmed = rsum / ndata
	resstd = sqrt( rsumsq/ndata - resmed*resmed )

c it		time in seconds
c perc		percentage of mass still in domain
c restime	residence time computed by integrating
c restimec	residence time computed by integrating with correction
c restime1	residence time computed by fitting regression curve
c resmed	average of residence times computed
c resstd	standard deviation of residence time

        write(iu,1000) it,perc,restime,restimec,restime1,resmed,resstd

c------------------------------------------------------------
c finish computation if mass is below threshold
c------------------------------------------------------------

        if( mass0 .ne. 0. .and. perc .lt. percmin ) then
                stop 'finished computing'
        end if

c------------------------------------------------------------
c no return flow -> set outside areas to 0
c------------------------------------------------------------

        if( bnoret ) then
          do k=1,nkn
            if( v1v(k) .eq. 0. ) then
              lmax = ilhkv(k)
              do l=1,lmax
                cnv(l,k) = 0.
              end do
            end if
          end do
        end if

c------------------------------------------------------------
c remember initialization
c------------------------------------------------------------

        icall = icall + 1

c------------------------------------------------------------
c end of routine
c------------------------------------------------------------

	return
 1000	format(i10,6f10.2)
        end

c*****************************************************************

        subroutine jamal_fra

c reset conz for fdp

	use mod_conz
	use levels
	use basin

        implicit none

        include 'param.h'

	include 'femtime.h'


        integer ie,ii,k,lmax,l,ia
        logical bnoret,breset,bstir
        real vol,conz,perc,percmin
        double precision mass,volume
        double precision massaux,volaux
        real volnode
	integer ifemop
	integer iaout,itmin,idtreset
	real dt,c0
	real restime,restime1,restimec
	real remnant,rlast
	real resmed,resstd

	integer iu,it0,ndata
	save iu,it0,ndata
	double precision remint,remlog,remtim
	save remint,remlog,remtim
	double precision rsum,rsumsq
	save rsum,rsumsq
        double precision mass0
        save mass0

        integer icall
        save icall
        data icall / 0 /

	itmin = 0
	idtreset = 0
	idtreset = nint( 1 * 30.5 * 86400 )		!one month is 30.5 days

c------------------------------------------------------------
c is it time to run the routine?
c------------------------------------------------------------

        if( it .le. itmin ) return

c------------------------------------------------------------
c initialization -> decide on reset
c------------------------------------------------------------

	breset = .false.		!normally do not reset

        if( icall .eq. 0 ) then
          write(6,*) 'initialization of routine jamal'
	  breset = .true.		!always reset at first call
	  it0 = it
	  it0 = itanf
        end if

	if( idtreset .gt. 0 ) then
	  if( it-it0 .ge. idtreset ) breset = .true.
	end if

c------------------------------------------------------------
c reset concentrations
c------------------------------------------------------------

        if( breset ) then		!reset concentrations to c0
          write(6,*) 'resetting concentrations in jamal at time ',it
          write(87,*) 'resetting concentrations in jamal at time ',it
          do k=1,nkn
            lmax = ilhkv(k)
            do l=1,lmax
              cnv(l,k) = 0
            end do
          end do
	  it0 = it
	end if

c------------------------------------------------------------
c end of routine
c------------------------------------------------------------

	end

c*****************************************************************

        subroutine sedimt

	use mod_conz, only : cnv
	use levels
	use basin

        implicit none

	real, save, allocatable :: conzs(:)
	real, save, allocatable :: conza(:)
	real, save, allocatable :: conzh(:)

        integer ie,ii,k,lmax,l,ia
	integer iunit,ierr
        real vol,conz,perc,wsink,dt,sed,h,r,cnew,rhos
	real v1v(nkn)
        double precision mass,masss
        double precision dtime

        real volnode,depnode
	double precision dgetpar
        logical next_output_d,is_in_output_d

	double precision, save :: da_out(4)
        integer, save :: icall = 0

	if( icall < 0 ) return

c------------------------------------------------------------
c parameters
c------------------------------------------------------------

	wsink = 0.
	wsink = 1.e-4
	wsink = 1.e-5
	wsink = 5.e-5
	rhos = 2500.

c------------------------------------------------------------
c initialization
c------------------------------------------------------------

        if( icall .eq. 0 ) then

          icall = 1

          write(6,*) 'initialization of routine sedimt: ',wsink

	  allocate(conzs(nkn))
	  allocate(conza(nkn))
	  allocate(conzh(nkn))
	  conzs = 0.
	  conza = 0.
	  conzh = 0.
	  cnv = 0.

	  call init_output_d('itmcon','idtcon',da_out)
	  call scalar_output_init(da_out,1,3,'set',ierr)
	  if( ierr > 0 ) goto 99
	  if( ierr < 0 ) icall = -1
	  if( icall < 0 ) return

        end if

c------------------------------------------------------------
c is it time ?
c------------------------------------------------------------

        if( .not. is_in_output_d(da_out) ) return

c------------------------------------------------------------
c sinking
c------------------------------------------------------------

	call get_act_dtime(dtime)
	call get_timestep(dt)
	call getinfo(iunit)

	if( wsink .gt. 0. ) then
	  l = 1
          do k=1,nkn
              h = depnode(l,k,+1)
              vol = volnode(l,k,+1)
	      r = 0.
	      if( h .gt. 0. ) r = wsink/h
              conz = cnv(l,k)
              conz = max(0.,conz)
	      cnew = conz * exp(-r*dt)
              cnv(l,k) = cnew
              sed = vol * (conz-cnew) 
	      conzs(k) = conzs(k) + sed
              sed = h * (conz-cnew) 
	      conza(k) = conza(k) + sed
              sed = sed / rhos
	      conzh(k) = conzh(k) + sed
          end do
	end if

c------------------------------------------------------------
c compute total mass
c------------------------------------------------------------

        mass = 0.
        masss = 0.
        do k=1,nkn
            lmax = ilhkv(k)
            do l=1,lmax
              vol = volnode(l,k,+1)
              conz = cnv(l,k)
              mass = mass + vol*conz
            end do
	    masss = masss + conzs(k)
        end do

c------------------------------------------------------------
c write total mass
c------------------------------------------------------------

        write(6,*) 'sedimt: ',dtime,mass,masss,mass+masss
        write(iunit,*) 'sedimt: ',dtime,mass,masss,mass+masss

c------------------------------------------------------------
c write file
c------------------------------------------------------------

        if( .not. next_output_d(da_out) ) return

	call scalar_output_write(dtime,da_out,22,1,conzs)	![kg]
	call scalar_output_write(dtime,da_out,23,1,conzs)	![kg/m**2]
	call scalar_output_write(dtime,da_out,24,1,conzs)	![m]

c------------------------------------------------------------
c end of routine
c------------------------------------------------------------

	return
   99	continue
	stop 'error stop sedimt: error opening file'
        end

c*****************************************************************

	subroutine oscillation

	implicit none

	include 'param.h'

	include 'femtime.h'

	real kenerg,penerg,tenerg

	integer icall
	data icall /0/
	save icall

	if( icall .eq. 0 ) then
	  icall = 1
	end if

	call energ(0,kenerg,penerg)

	tenerg = kenerg + penerg
	write(9,*) it,kenerg,penerg,tenerg

	end

c*****************************************************************

	subroutine kreis

	implicit none

	include 'param.h'

	include 'femtime.h'

	real kenerg,penerg,tenerg

	integer icall
	data icall /0/
	save icall

	if( icall .eq. 0 ) then
	  icall = 1
	  call init_kreis()
	end if

	call energ(0,kenerg,penerg)

	tenerg = kenerg + penerg
	write(9,*) it,kenerg,penerg,tenerg

	end

c*****************************************************************

	subroutine init_kreis

	use mod_depth
	use mod_area
	use mod_hydro_vel
	use mod_hydro
	use levels, only : nlvdi,nlv
	use basin

	implicit none

	integer k,ie,ii,l
	real pi,dcori,z0,r0,f,omega,grav
	real aux,r02,r2
	real x,y,z,u,v
	!real v1v(nkn)

	pi = 4.*atan(1.)
	dcori = 45.
	z0 = 0.
	r0 = 0.
	f = 2.0 * 0.729E-4 * sin(dcori*pi/180.)
	omega = 0.5E-5
	grav = 9.81

	aux = (omega*f)/(2.*grav)
	r02 = r0*r0

	do k=1,nkn
	  x = xgv(k)
	  y = ygv(k)
	  r2 = x*x+y*y
	  z = z0 + aux * ( r2 - r02 )

	  znv(k) = z
	end do

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    zenv(ii,ie) = znv(k)
	  end do
	end do

	call make_new_depth

	do ie=1,nel
	  call baric(ie,x,y)
	  u = -omega * y
	  v = +omega * x
	  do l=1,nlv
	    ulnv(l,ie) = u
	    vlnv(l,ie) = v
	  end do
	end do

	call vtot
	call uvint
	call make_prvel

	end

c*****************************************************************

c*****************************************************************
        subroutine bclevvar

	use mod_depth
	use mod_layer_thickness
	use mod_hydro_vel
	use mod_hydro
	use basin, only : nkn,nel,ngr,mbw

        implicit none 

        include 'param.h' 


	include 'femtime.h'

        integer edim
        parameter(edim = 20) !20!number of element near the open boundary
        integer elese(edim)!elements near the boundary
c  !      data elese / 290,321,326,335,308,2112,324,336,334,325,
c !    +     2118,2113,328,295,2119,2116,333,2121,2120,2114
c!     +      /                                                !grid chao0_primapubb.grd
        data elese / 3002,3001,2902,2901,2802,2801,2702,2701,
     +  2602,2601,3004,3003,2904,2903,2804,2803,2704,2703,2604,2603
     +      /    !grid chao0_ordinato.grd
c     !   data elese /382,381,380,379,378,377,376,375,374,373,372,371,
c     !+  369,368,367,366,365,364,363,362,361,360,359,358,357,356,
c     !+  355,354,353
c     !+         /!grid bati_gradino.grd
        integer eleint(edim)
        integer ie,l,i
        logical berror          !true on return if error
        

        real u,v
        real ut,vt
        real h

        integer nlev
        parameter (nlev = 5)
      
        real upresc(nlev),vpresc(nlev)
        real umax
        real alpha

        umax = 0.1 !unity m/s
        upresc(1) = umax !vertical profile u
        upresc(2) = umax/2
        upresc(3) = -umax/2!-umax/2
        upresc(4) = -umax/2!-umax/2
        upresc(5) = -umax/2!-umax/2!deb220307
        
        vpresc(1) = 0. !vertical profile v
        vpresc(2) = 0.
        vpresc(3) = 0.
        vpresc(4) = 0.
        vpresc(5) = 0.

        do i=1,edim 
            eleint(i)=elese(i)
        enddo
                call e2int(edim,eleint,berror)
            
        !print*,eleint 
        if (it .le. 43200)then
                alpha=it/43200.
                !write(6,*)'funziona',it,alpha
        else
                alpha=1.
        endif

        do i=1,edim
                ie=eleint(i)
                do l=1,nlev
                       h = hdeov(l,ie)
                       ulnv(l,ie) = upresc(l)*alpha
                       vlnv(l,ie) = vpresc(l)*alpha
                !   write(6,*)'funzionaDEB',ie,l,ulnv(l,ie),vlnv(l,ie)
                       utlnv(l,ie) = ulnv(l,ie)*h
                       vtlnv(l,ie) = vlnv(l,ie)*h
         !  write(6,*) l,' h ',h,' ulnv ',ulnv(l,ie),' vlnv ',vlnv(l,ie)
                   !write(6,*) l,'utlnv',utlnv(l,ie),'vtlnv',vtlnv(l,ie)
                enddo
        enddo

        end
c*****************************************************************
        subroutine bclevvar_ini !dbf

	use mod_internal
	use mod_depth
	use mod_layer_thickness
	use mod_hydro_vel
	use mod_hydro
	use basin, only : nkn,nel,ngr,mbw

        implicit none 

        include 'param.h' 


        integer edim
        parameter(edim = 20) !number of element near the open boundary
        integer elese(edim)!elements near the boundary
c!        data elese / 290,321,326,335,308,2112,324,336,334,325,
c!     +     2118,2113,328,295,2119,2116,333,2121,2120,2114
c!     +      /                                                !grid chao0_primapubb.grd
          data elese / 3002,3001,2902,2901,2802,2801,2702,2701,
     +  2602,2601,3004,3003,2904,2903,2804,2803,2704,2703,2604,2603 
     +  /    !grid chao0_ordinato.grd
c    !   data elese /382,381,380,379,378,377,376,375,374,373,372,371,
c    ! +  369,368,367,366,365,364,363,362,361,360,359,358,357,356,
c    ! +  355,354,353
c    ! +     /!grid bati_gradino.grd

        integer eleint(edim)
        integer ie,l,i,ii
        logical berror          !true on return if error
       




        real u,v
        real ut,vt
        real h

        integer nlev
        parameter (nlev = 5)
      
        real upresc(nlev),vpresc(nlev)
        real umax

        do ie=1,nel
                iuvfix(ie)= 0
        enddo

        umax = 0.!0.1 !unity m/s
        upresc(1) = umax !vertical profile u
        upresc(2) = umax/2
        upresc(3) = -umax/2
        upresc(4) = -umax/2
        upresc(5) = -umax/2
        
        vpresc(1) = 0. !vertical profile v
        vpresc(2) = 0.
        vpresc(3) = 0.
        vpresc(4) = 0.
        vpresc(5) = 0.

        do i=1,edim 
            eleint(i)=elese(i)
        enddo
                call e2int(edim,eleint,berror)
            
        !print*,eleint 
        do i=1,edim
                ie=eleint(i)
                iuvfix(ie)=1
                do l=1,nlev
                       h = hdeov(l,ie)
                       ulov(l,ie) = upresc(l)!dbf 12feb2007
                       vlov(l,ie) = vpresc(l)!dbf 12feb2007
                       ulnv(l,ie) = upresc(l)
                       vlnv(l,ie) = vpresc(l)
                       utlov(l,ie) = ulov(l,ie)*h
                       utlnv(l,ie) = ulnv(l,ie)*h
                       vtlov(l,ie) = vlov(l,ie)*h
                       vtlnv(l,ie) = vlnv(l,ie)*h
           !write(6,*)'inizializzo velocita bclevvar_ini'
           !write(6,*) l,' h ',h,' ulnv ',ulov(l,ie),' vlnv ',vlov(l,ie)
           !write(6,*) l,' h ',h,' ulnv ',ulnv(l,ie),' vlnv ',vlnv(l,ie)
                   !write(6,*) l,'utlnv',utlnv(l,ie),'vtlnv',vtlnv(l,ie)
                enddo
        enddo

        end


c****************************************************************************


c*****************************************************************

	subroutine debora(it)

	implicit none

	integer it

	integer icall
	save icall
	data icall /0/

	if( icall .gt. 0 ) return
	if( it .lt. 258000 ) return
	!if( it .lt. 3000 ) return

	icall = 1

	write(6,*) 'changed in debora...'
	call putpar('ilin',0.)
	call putpar('itsplt',2.)

	end

c*****************************************************************

        subroutine tsinitdebora(it)

	use mod_ts
	use levels, only : nlvdi,nlv
	use basin

	implicit none

        include 'param.h'

	integer ie,k,i,l,it
        integer itype

	integer icall
	save icall
	data icall /0/

        if(icall .gt. 0) return

        do l=1,nlv
             do ie=1,nel
                itype = iarv(ie)        
               do i=1,3
                        k=nen3v(i,ie)
                        if( itype .lt. 90) then
                          saltv(l,k)=30.
                          tempv(l,k)=14.
                        else
                          saltv(l,k)=35.
                          tempv(l,k)=14.
                        endif
                end do
             end do
        end do

        icall = 1

        end

c*****************************************************************

	subroutine zinit

	use mod_hydro
	use basin

	implicit none

	include 'param.h'

	integer k,ie,ii
	real xmin,xmax,ymin,ymax
	real zav,dz,pi
	real x,y,r,z

	integer icall
	save icall
	data icall / 0 /

	if( icall .gt. 0 ) return

	icall = 1

	zav = 0.
	dz = .20

	pi = 4.*atan(1.)

	xmax = xgv(1)
	xmin = xgv(1)
	ymax = ygv(1)
	ymin = ygv(1)

	do k=1,nkn
	  xmax = max(xmax,xgv(k))
	  xmin = min(xmin,xgv(k))
	  ymax = max(ymax,ygv(k))
	  ymin = min(ymin,ygv(k))
	end do

	do k=1,nkn
	  x = xgv(k)
	  y = ygv(k)
	  r = (y - ymin)/(ymax-ymin)	! r is in [0:+1]
	  z = zav + (2*r-1.)*dz		! linear
	  z = zav - dz * cos(r*pi)	! sinusoidal
	  znv(k) = z
	end do

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    zenv(ii,ie) = znv(k)
	  end do
	end do

	end

c****************************************************************

	subroutine cprint(it)

	use mod_conz
	use mod_ts
	use levels

	implicit none

	include 'param.h'

	integer it

	integer i,k,lmax,l
	logical berror


	integer icall
	save icall
	data icall / 0 /

	integer ndim
	parameter (ndim=5)
	integer nodes(ndim)
	save nodes
	!data nodes / 984, 4860, 4636, 4585 /
	data nodes / 984, 4860, 4636, 4585 , 3452 /

	if( icall .eq. 0 ) then
	  icall = 1
	  call n2int(ndim,nodes,berror)
	  if( berror ) stop 'error stop cprint'
	end if

	write(87,*)
	write(87,*) 'time = ',it
	write(87,*)

	do i=1,ndim
	  k = nodes(i)
	  lmax = ilhkv(k)
	  write(87,*) i,k,lmax,(cnv(l,k),l=1,lmax)
	  write(81,*) i,k,lmax,(saltv(l,k),l=1,lmax)
	  write(82,*) i,k,lmax,(tempv(l,k),l=1,lmax)
	end do

	write(86,*) it,(cnv(1,nodes(i)),i=1,ndim)

	end

c****************************************************************
c****************************************************************
c****************************************************************
c****************************************************************
c****************************************************************

	subroutine diffus2d

c test for 2D horizontal diffusion algorithm

c dC/dt = a*laplace(C)    with    c(x,0+)=delta(x)
c
c C(x,t) =  (4*pi*a*t)**(-n/2) * exp( -|x|**2/(4*a*t) )		(n=n)
c
c C(x,t) =  1/sqrt(4*pi*a*t) * exp( -x**2/(4*a*t) )		(n=1)
c
c C(x,t) =  1/(4*pi*a*t) * exp( -|x|**2/(4*a*t) )		(n=2)

	use mod_conz
	use levels, only : nlvdi,nlv

	implicit none

	include 'param.h'

	include 'femtime.h'

	integer ndim
	parameter (ndim=101)

	integer k,k0,i,kin,n
	integer kmin,kmax,kstep
	integer itt,idtfrq,it0
	real c0,ctot,ct0,c
	real dx,rk
	real depth
	real caux(ndim)
	real xv(2*ndim)
	real yv(2*ndim)
	character*80 file

	integer ipint
	real getpar

	save ct0

	integer icall
	save icall
	data icall / 0 /

c----------------------------------------------------------
c initialization
c----------------------------------------------------------

	idtfrq = 500
	it0 = 1000		!initialize with this time

	c0 = 100.
	k0 = 51051
	dx = 100.
	depth = 10.

	if( icall .eq. 0 ) then
	  rk = getpar('chpar')
	  kin = ipint(k0)
	  do k=kin,kin
	    cnv(1,k) = c0
	  end do
	  call tsmass(cnv,+1,nlvdi,ctot)
	  write(6,*) 'diffus2d initialized (1): ',c0,ctot
	  call cv_init(it0,k0,rk,c0,ct0,cnv)
	  call tsmass(cnv,+1,nlvdi,ctot)
	  write(6,*) 'diffus2d initialized: ',it0,it,itanf,idt
	  write(6,*) 'diffus2d initialized: ',c0,ctot,ct0
	  ct0 = ctot/depth
	end if

	icall = icall + 1
	itt = it + it0 - idt

c----------------------------------------------------------
c compute average over section
c----------------------------------------------------------

	if( mod(itt,idtfrq) .ne. 0 ) return

	rk = getpar('chpar')
	call tsmass(cnv,+1,nlvdi,ctot)
	kin = ipint(k0)
	c = cnv(1,kin)
	write(6,*) 'diffus2d: ',c0,c,ctot,ct0

	kmin = 51001
	kmax = 51101
	kstep = 1
	call extract_conz(cnv,ndim,xv,yv,k0,kmin,kmax,kstep)
        call make_filename(itt,file,'hor','.txt')
	call write_file(file,ndim,xv,yv)

	kmin =   1051
	kmax = 101051
	kstep = 1000
	call extract_conz(cnv,ndim,xv,yv,k0,kmin,kmax,kstep)
        call make_filename(itt,file,'ver','.txt')
	call write_file(file,ndim,xv,yv)

	kmin =   1001
	kmax = 101101
	kstep = 1001
	call extract_conz(cnv,ndim,xv,yv,k0,kmin,kmax,kstep)
        call make_filename(itt,file,'dia','.txt')
	call write_file(file,ndim,xv,yv)

	call make_anal(itt,ndim,xv,yv,ct0,rk,dx)
        call make_filename(itt,file,'ana','.txt')
	call write_file(file,ndim,xv,yv)

	call extract_irreg(cnv,2*ndim,n,xv,yv,k0,dx,60.)
        call make_filename(itt,file,'p60','.txt')
	call write_file(file,n,xv,yv)

c----------------------------------------------------------
c end of routine
c----------------------------------------------------------

	end

c****************************************************************

	subroutine write_file(file,ndim,xv,yv)

	implicit none

	character*(*) file
	integer ndim
	real xv(ndim)
	real yv(ndim)

	integer i,iu

	iu = 88

	open(iu,file=file,status='unknown',form='formatted')
	do i=1,ndim
	  write(iu,*) xv(i),yv(i)
	end do
	close(iu)

	end

c****************************************************************

	subroutine make_anal(it,ndim,xv,yv,c0,rk,dx)

c C(x,t) =  1/(4*pi*a*t) * exp( -|x|**2/(4*a*t) )		(n=2)

	implicit none

	include 'param.h'

	integer it
	integer ndim
	real xv(ndim),yv(ndim)
	real c0,rk,dx

	integer i,nc
	real pi,aux,x,y,r
	real ctot,dxx,daux,fact
	double precision dtot

	pi = 4. * atan(1.)
	aux = 4. * rk * it
	nc = 1 + ndim / 2

	do i=1,ndim
	  x = (i-nc)*dx
	  y = c0 * exp(-x*x/aux) / (pi*aux)
	  xv(i) = x
	  yv(i) = y
	end do

	fact = 0.9
	fact = 0.95
	r = nc * dx
	dtot = 0.
	do i=1,1000
	  x = r
	  y = c0 * exp(-x*x/aux) / (pi*aux)
	  dxx = (1.-fact) * r
	  daux = 2.*pi*r*y*dxx
	  dtot = dtot + daux
	  if( dtot .gt. 100 .and. daux .lt. 1.e-10 ) goto 1
	  !write(6,*) i,x,y,daux,real(dtot)
	  r = fact * r
	end do
	write(6,*) i,daux
	stop 'error stop make_anal: error too high'

    1	continue
	ctot = dtot
	write(6,*) 'analytic solution: ',i,y,daux,ctot

	end

c****************************************************************

	subroutine cv_init(it,k0,rk,c0,ct0,cv)

c initializes cnv with analytic solution (2D)

	use levels, only : nlvdi,nlv
	use basin

	implicit none

	include 'param.h'

	integer it,k0
	real rk,c0,ct0
        real cv(nlvdi,nkn)

	integer i,k,kin
	real x,y,dx,dy,r
	real x0,y0
	real pi,aux,ctot

	integer ipint

	pi = 4. * atan(1.)
	aux = 4. * rk * it
	ctot = c0 * pi * aux
	ct0 = ctot

	kin = ipint(k0)
	x0 = xgv(kin)
	y0 = ygv(kin)

	do k=1,nkn
	  x = xgv(k)
	  y = ygv(k)
	  dx = x-x0
	  dy = y-y0
	  r = sqrt( dx*dx + dy*dy )
	  y = ctot * exp(-r*r/aux) / (pi*aux)
	  cv(1,k) = y
	end do

	end

c****************************************************************

	subroutine extract_irreg(cv,ndim,n,xv,yv,k0,dx,phi)

c extracts conz from array

	use levels, only : nlvdi,nlv
	use basin

	implicit none

	include 'param.h'

        real cv(nlvdi,nkn)
	integer ndim
	integer n
	real xv(ndim)
	real yv(ndim)
	integer k0
	real dx,phi

	integer i,kin0,ie,ip
	real x0,y0
	real pi,rad,phi0
	real r,x,y
	real zp

	integer ipint

	write(6,*) 'starting to extract... phi = ',phi

	kin0 = ipint(k0)
	x0 = xgv(kin0)
	y0 = ygv(kin0)
	write(6,*) 'center... ',k0,kin0,x0,y0

c------------------------------------------------------------
c prepare parameters
c------------------------------------------------------------

	xv(1) = 0.
	yv(1) = cv(1,kin0)

	pi = 4.*atan(1.)
	rad = pi/180.
	phi0 = phi * rad

c------------------------------------------------------------
c go backward
c------------------------------------------------------------

	ie = 0
	do i=2,ndim
	  r = -(i-1)*dx
	  x = x0 + r * cos(phi0)
	  y = y0 + r * sin(phi0)
	  call intp_on_node(ie,x,y,cv,zp)
	write(6,*) ie,r,x,y,zp
	  if( ie .le. 0 ) goto 1
	  xv(i) = r
	  yv(i) = zp
	end do
	stop 'error stop extract_irreg: too many nodes (1)'
    1	continue
	n = i-1
	write(6,*) 'backward... ',n
	
c------------------------------------------------------------
c invert arrays
c------------------------------------------------------------

	do i=1,n/2
	  ip = n+1-i
	  call swap_ggu(xv(i),xv(ip))
	  call swap_ggu(yv(i),yv(ip))
	end do

c------------------------------------------------------------
c go foreward
c------------------------------------------------------------

	ie = 0
	do i=n+1,ndim
	  r = (i-n)*dx
	  x = x0 + r * cos(phi0)
	  y = y0 + r * sin(phi0)
	  call intp_on_node(ie,x,y,cv,zp)
	  if( ie .le. 0 ) goto 2
	  xv(i) = r
	  yv(i) = zp
	end do
	stop 'error stop extract_irreg: too many nodes (2)'
    2	continue
	n = i-1
	write(6,*) 'foreward... ',n
	
	write(6,*) 'finished to extract... n = ',n

c------------------------------------------------------------
c end of routine
c------------------------------------------------------------

	end

c****************************************************************

	subroutine swap_ggu(a,b)

	implicit none

	real a,b

	real aux

	aux = a
	a = b
	b = aux

	end

c****************************************************************

	subroutine intp_on_node(ie,x,y,cv,zp)

	use levels, only : nlvdi,nlv
	use basin

	implicit none

	include 'param.h'

	integer ie
	real x,y
	real cv(nlvdi,nkn)
	real zp

	integer ii,k
	real z(3)

	call find_elem_from_old(ie,x,y,ie)

	if( ie .le. 0 ) return

	do ii=1,3
	  k = nen3v(ii,ie)
	  z(ii) = cv(1,k)
	end do

	call femintp(ie,z,x,y,zp)

	end

c****************************************************************

	subroutine extract_conz(cv,ndim,xv,yv,k0,kmin,kmax,kstep)

c extracts conz from array

	use levels, only : nlvdi,nlv
	use basin

	implicit none

	include 'param.h'

        real cv(nlvdi,nkn)
	integer ndim
	real xv(ndim)
	real yv(ndim)
	integer k0,kmin,kmax,kstep

	integer i,k,kin
	real x,y,dx,dy,dxy
	real x0,y0
	real fact

	integer ipint

	kin = ipint(k0)
	x0 = xgv(kin)
	y0 = ygv(kin)

	i = 0
	fact = -1
	do k=kmin,kmax,kstep
	  kin = ipint(k)
	  if( k .eq. k0 ) fact = +1
	  x = xgv(kin)
	  y = ygv(kin)
	  dx = x-x0
	  dy = y-y0
	  dxy = sqrt( dx*dx + dy*dy )
	  i = i + 1
	  if( i .gt. ndim ) stop 'error stop extract_conz: ndim'
	  xv(i) = fact*dxy
	  yv(i) = cv(1,kin)
	end do

	end

c****************************************************************
c****************************************************************
c****************************************************************
c****************************************************************
c****************************************************************

	subroutine diffus

c test for horizontal diffusion algorithm

	use mod_conz

	implicit none

	include 'param.h'

	include 'femtime.h'

	integer k,k0,k1,i
	real c0,cmed
	real caux(25)

	integer icall
	save icall
	data icall / 0 /

c----------------------------------------------------------
c initialization
c----------------------------------------------------------

	c0 = 100.

	if( icall .eq. 0 ) then
	  do k=109,117
	    cnv(1,k) = c0
	  end do
	end if

	icall = icall + 1

c----------------------------------------------------------
c compute average over section
c----------------------------------------------------------

	do i=1,25
	  k1 = i*9
	  k0 = k1 - 8
	  cmed = 0.
	  do k=k0,k1
	    cmed = cmed + cnv(1,k)
	  end do
	  cmed = cmed / 9
	  caux(i) = cmed
	end do

	write(88,*) it,caux

c----------------------------------------------------------
c end of routine
c----------------------------------------------------------

	end

c****************************************************************

	subroutine viscos

c viscosity algorithm

	use mod_hydro_print
	use levels, only : nlvdi,nlv

	implicit none

	include 'param.h'

	include 'femtime.h'


	integer k,i,l

	write(88,*) it
	do k=109,117
	  write(88,'(12i5)') k,nlv,(nint(1000*vprv(l,k)),l=1,nlv)
	end do

	write(89,*) it
	do l=1,nlv
	  write(89,'(12i5)') l,(nint(1000*vprv(l,k)),k=109,117)
	end do

	end

c****************************************************************

	subroutine vdiffus(mode)

c diffusion algorithm
c
c solution of purely diffusional part :
c
c dC/dt = a*laplace(C)    with    c(x,0+)=delta(x)
c
c C(x,t) =  (4*pi*a*t)**(-n/2) * exp( -|x|**2/(4*a*t) )
c
c for n-dimensions and
c
c C(x,t) =  1/sqrt(4*pi*a*t) * exp( -x**2/(4*a*t) )
c
c for 1 dimension
c
c the solution is normalized, i.e.  int(C(x,t)dx) = 1 over the whole area

	use mod_conz
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer mode

	include 'param.h'

	include 'femtime.h'

	integer k,i,l,lc,lmax
	integer it0
	real c0,cmed,pi,a,t,aux,x2,cc,dc
	real caux(nlvdi)
	save it0

	real getpar

	integer icall
	save icall
	data icall / 0 /

c C(x,t) =  1/sqrt(4*pi*a*t) * exp( -x**2/(4*a*t) )

	lc = (nlv+1)/2
	c0 = 100.
	pi = 4. * atan(1.)
	a = getpar('diftur')
	dc = -0.5

	if( icall .eq. 0 ) then
	  if( mode .eq. 1 ) then
	    it0 = it
	    do k=1,nkn
	      cnv(lc,k) = c0
	    end do
	  else
	    do k=1,nkn
	      lmax = ilhkv(k)
	      do l=1,lmax
	        cnv(l,k) = c0 + l * dc
	      end do
	    end do
	  end if
	  icall = 1
	  write(68,*) 'vdiffus: ',mode,lc,c0,pi,a,it0  
	end if

	t = it - it0
	aux = 4. * a * t
	do l=1,nlv
	  x2 = (l-lc)**2
	  cc = 0.
	  if( aux .gt. 0. ) cc = c0 * exp( -x2 / aux ) / sqrt( pi * aux )
	  caux(l) = cc
	end do
	if( aux .eq. 0 ) caux(lc) = c0

	write(68,*) it,nlv
	do l=1,nlv
	  !k = i*9 - 4
	  write(68,'(i4,6f12.4)') l,(cnv(l,i*9-4),i=1,25,6),caux(l)
	end do

	end

c****************************************************************

        subroutine uv_bottom

	use levels
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        include 'param.h'


        integer k,l,m
        real u,v
        double precision uz,cdir

        real rdebug(nkn)

        integer iudeb
        save iudeb
        data iudeb /0/

        do k = 1,nkn

          l = ilhkv(k)                          !bottom layer
          m = l

          call getuv(m,k,u,v)                   !gets velocities u/v
          call getmd(u,v,UZ,CDIR)               !gets UZ, CDIR

          rdebug(k) = uz
        end do

        write(6,*) 'debug value written... ',iudeb
        call conwrite(iudeb,'.ggu',1,888,1,rdebug)

        end

c********************************************************************

	subroutine joel

c writes node list for Malta Coastal Model

	use mod_depth
	use levels
	use basin

	implicit none

        include 'param.h'



        integer k,l,ie,ii,i
	integer nbc
        real x,y,h

	real x0,y0,phi0
	real pi,rad
	real xfact,yfact
	real v1v(nkn)

	integer ibc,nnodes,kext
	integer nkbnds,kbnds,ipext,nbnds
	logical bproj

	integer icall
	save icall
	data icall / 0 /

c----------------------------------------------------------
c do only at first call
c----------------------------------------------------------

	if( icall .ne. 0 ) return
	icall = 1

c----------------------------------------------------------
c initialize geographic projection
c----------------------------------------------------------

	bproj = .false.

	x0 = 0.
	y0 = 0.
	xfact = 1.
	yfact = 1.

	if( bproj ) then
	  x0 = 14.047
	  y0 = 35.672
	  phi0 = 36.

	  pi = 4. * atan(1.)
	  rad = pi/180.

	  yfact = 60*1852
	  xfact = yfact * cos(phi0*rad)
	end if

c----------------------------------------------------------
c get deepest depth on node
c----------------------------------------------------------

        do k = 1,nkn
	  v1v(k) = -999.
	end do

        do ie = 1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    v1v(k) = max(v1v(k),hev(ie))
	  end do
	end do

c----------------------------------------------------------
c write out values for total domain
c----------------------------------------------------------

	open(1,file='joel_total.dat')

	write(1,*) nkn
        do k = 1,nkn
	  kext = ipext(k)
	  l = ilhkv(k)
	  h = v1v(k)
	  x = xgv(k) / xfact + x0
	  y = ygv(k) / yfact + y0
	  write(1,1000) k,kext,l,x,y,h
	end do

	close(1)

c----------------------------------------------------------
c write out values for open boundary
c----------------------------------------------------------

	nbc = nbnds()

	if( nbc .gt. 0 ) then
	open(1,file='joel_bound.dat')

	ibc = 1
	nnodes = nkbnds(ibc)

	write(1,*) nnodes
        do i = 1,nnodes
	  k = kbnds(ibc,i)			!get boundary nodes
	  kext = ipext(k)
	  l = ilhkv(k)
	  h = v1v(k)
	  x = xgv(k) / xfact + x0
	  y = ygv(k) / yfact + y0
	  write(1,1000) k,kext,l,x,y,h
	end do

	close(1)
	end if

c----------------------------------------------------------
c write out values z-levels
c----------------------------------------------------------

	open(1,file='joel_zlevels.dat')

	write(1,*) nlv
        do l = 1,nlv
	  h = hlv(l)
	  write(1,*) l,h
	end do

	close(1)

c----------------------------------------------------------
c write boundary nodes
c----------------------------------------------------------

        call print_bound_nodes(nkn,v1v)

c----------------------------------------------------------
c end of routine
c----------------------------------------------------------

 1000	format(2i10,i5,3f16.8)
	end

c********************************************************************

        subroutine andreac

c introduce le condizioni iniziali per la concentrazione

	use mod_conz
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        include 'param.h'


	integer ks

        integer icall
        save icall
        data icall / 0 /

	ks = 4312
        if( icall .eq. 0 .and. nkn .ge. ks .and. ks > 0 ) then
          cnv(1,ks-1)=100.
          cnv(1,ks)=100.
       	  icall=1
        end if

	end

c********************************************************************

        subroutine tvd_test(it)

c introduce le condizioni iniziali per la concentrazione

	use mod_conz
	use basin

        implicit none

	integer it

        include 'param.h'

	integer ndim
	parameter (ndim=2000)

	integer nelems
	integer iel(ndim)
	real conz(ndim)
	real yco(ndim)
	save nelems,iel,conz,yco


	integer i,ii,ie,k,ien
	integer igrad
	real grad,gradaux,gradmax
	real dy,x0,y,ymin,ymax
	real cc(3)
	real ccc,cmin,cmax

        integer icall
        save icall
        data icall / 0 /

	dy = 100.
	x0 = 1100.
	call mima(ygv,nkn,ymin,ymax)

        if( icall .eq. 0 ) then
       	  icall = 1

	  y = ymin
	  i = 0
	  ie = 0
	  do while( y+dy .le. ymax )
	    i = i + 1
	    if( i .gt. ndim ) stop 'error stop tvd_test: ndim'
	    y = ymin + (i-1) * dy + dy/2.
	    call find_elem_from_old(ie,x0,y,ien)
	    !write(6,*) ymin,ymax,x0,y,ie,ien
	    if( ien .le. 0 ) stop 'error stop tvd_test (1): iel'
	    iel(i) = ien
	    yco(i) = y
	    ie = ien
	  end do

	  nelems = i
	  write(6,*) 'set up of tvd test: ',ymin,ymax,nelems

	  write(68,*) 46728645,1
	  write(68,*) nelems,0
	  write(68,*) 0,nelems+1,-10,110.

        end if

	do i=1,nelems

	  y = yco(i)
	  ie = iel(i)

	  do ii=1,3
	    k = nen3v(ii,ie)
	    cc(ii) = cnv(1,k)
	  end do

	  call femintp(ie,cc,x0,y,ccc)

	  conz(i) = ccc
	end do

	cmin = +1.e+30
	cmax = -1.e+30
	do k=1,nkn
	  cmin = min(cmin,cnv(1,k))
	  cmax = max(cmax,cnv(1,k))
	end do
	write(71,*) it,cmin,cmax	!only for check
	cmin = min(cmin,0.)
	cmax = max(cmax,100.)
	write(70,*) it,cmin,cmax-100.

	write(68,*) it,nelems
	write(68,*) (conz(i),i=1,nelems)

	grad = 0.
	igrad = 0
	do i=4,nelems-1			!really start from 2 -> bug
	  gradaux = abs(conz(i+1)-conz(i-1))
	  if( gradaux .gt. grad ) then
	    grad = gradaux
	    igrad = i
	  end if
	end do

	gradmax = 100./(2.*dy)
	grad = grad / (2.*dy)
	grad = 100. * grad / gradmax	!is fraction of max grad in percentage
	write(69,*) it,grad,igrad

        end

c********************************************************************

        subroutine make_filename(it,file,pre,ext)

        implicit none

        integer it
        character*(*) file,pre,ext

        integer i,n
        integer ichanm

        write(file,'(a,i10,a)') pre,it,ext
        n = ichanm(file)
        do i=n,1,-1
          if( file(i:i) .eq. ' ' ) file(i:i) = '0'
        end do

        end

c**********************************************************************

	subroutine black_sea_nudge

	use mod_ts
	use levels, only : nlvdi,nlv
	use basin

	implicit none

	include 'param.h'

	integer k,l
	real xc1,yc1,xc2,yc2
	real xm1,ym1,xm2,ym2
	real tc,tm
	real scalc,scalm
	real ncx,ncy,nmx,nmy
	real x,y
	real tau

	integer icall
	save icall
	data icall / 0 /

	if( icall .ne. 0 ) return
	icall = 1

	xc1 = 29.
	yc1 = 43.5
	xc2 = 30.5
	yc2 = 45.5

	xm1 = 29.5
	ym1 = 43.5
	xm2 = 31.
	ym2 = 45.5

	tc = 259200.
	tm = 86400.

	ncx = -(yc2-yc1)
	ncy = +(xc2-xc1)
	nmx = -(ym2-ym1)
	nmy = +(xm2-xm1)

	do k=1,nkn
	  x = xgv(k)
	  y = ygv(k)

	  scalc = ncx * (x-xc1) + ncy * (y-yc1)
	  scalm = nmx * (x-xm1) + nmy * (y-ym1)

	  if( scalc .gt. 0. ) then		!coast
	    tau = tc
	  else if( scalm .lt. 0. ) then		!open sea
	    tau = tm
	  else					!intermediate
	    tau = 0.5 * (tc+tm)
	  end if

	  if( tau .gt. 0. ) tau = 1. / tau

	  do l=1,nlvdi
	    ttauv(l,k) = tau
	    stauv(l,k) = tau
	  end do
	end do

	write(6,*) 'relaxation time for nudging set up'

	end

!**********************************************************************
! nudging in funzione della dimensione degli elementi in modo da ottenere
! un relaxation factor maggiore a costa rispettio al mare aperto
! Discrimina a seconda dei type degli elementi

        subroutine marmara_sea_nudge

        use mod_ts
        use levels, only : nlvdi,nlv
        use basin
        use evgeom

        implicit none

!       include 'param.h'
!       include 'basin.h'
!        include 'ev.h'
!        include 'ts.h'

        integer k,l,ii,ie
        real rgil,rmed,rbla,rmar     !relaxaxion factors
        real rop,rco            !relaxaxion factor for open sea and coast

        real, save, allocatable :: tau(:)
        real, save, allocatable :: rtd(:)
        real, save, allocatable :: ark(:)
        real area,areal,areak
        real ttt
	real xe,ye
        integer iet

        integer ia_out(4)
        save ia_out
        logical has_output,next_output
        logical surface

        integer icall
        save icall
        data icall / 0 /

        !-------------------------------------------
        ! Set surface or 3D nudging
        !-------------------------------------------
        surface = .true.

        !-------------------------------------------
        ! Set relaxation times
        !-------------------------------------------
	! For surface water temperature the value of the relaxation
	! coefficient is spatially varying over the model domain from
	! 8 days in the open sea and increasing, thus diminishing the
	! restoration contribution, toward the coast, up to 40 days 
	! in the TSS. For surface salinity the relaxation time is double
	! the one for temperature. We impose strong nudging of T/S at
	! Gibraltar over the whole water column as boundary conditions.

        rmed = 8.*86400.       !type = 5
        rbla = 8.*86400.       !type = 1
        rmar = 40.*86400.      !type = 2,3,4
        rgil = 3.*86400.       !type = 11

        !-------------------------------------------
        ! Initialization
        !-------------------------------------------
        if( icall .ne. 0 ) return
        icall = 1

        allocate(tau(nkn))
        allocate(ark(nkn))
        allocate(rtd(nkn))

        do k = 1,nkn
          tau(k) = 0.
          ark(k) = 0.
        end do

        !-------------------------------------------
        ! Set relaxation based on type and element size
        !-------------------------------------------
        do ie = 1,nel
          iet = iarv(ie)
          area = 12. * ev(10,ie)
	  call baric(ie,xe,ye)
          if ( iet .eq. 1 ) then        !black sea
            areal = 1E+07               !average element are in the open sea
            rop = rbla
            rco = rmar
            ttt = rop + (areal - area)*(rco-rop)/areal
            if (area .gt. areal) ttt  = rop
          elseif (iet .ge. 2 .and. iet .le. 4 ) then    !marmara sea
            ttt  = rmar
          elseif (iet .eq. 5 ) then     !mediterranean sea
            areal = 1E+07               !average element are in the open sea
            rop = rmed
            rco = rmar
            ttt = rop + (areal - area)*(rco-rop)/areal
            if (area .gt. areal) ttt  = rop
          elseif (iet .eq. 6) then      !out bosphorous
            areal = 1E+07               !average element are in the open sea
            rop = rbla
            rco = rmar
            ttt = rop + (areal - area)*(rco-rop)/areal
            ttt = max(ttt,rop)
          elseif (iet .eq. 7) then      !out dardanelles
            areal = 1E+07               !average element are in the open sea
            rop = rbla
            rco = rmar
            ttt = rop + (areal - area)*(rco-rop)/areal
            ttt = max(ttt,rop)
          elseif ( iet .eq. 10) then    !Azov sea
            areal = 1E+07               !average element are in the open sea
            rop = rbla
            rco = rmar
            ttt = rop + (areal - area)*(rco-rop)/areal
            if (area .gt. areal) ttt  = rop
          elseif ( iet .eq. 11) then    !Gibraltar Strait
            if ( xe .lt. -6.0 ) then
              areal = 2.2E+07               !average element are in the open sea
              rop = rgil
              rco = rmar
              ttt = rop + (areal - area)*(rco-rop)/areal
              if (area .gt. areal) ttt  = rop
	      !ttt = rgil
            else
              areal = 2.2E+07               !average element are in the open sea
              rop = rmed
              rco = rmar
              ttt = rop + (areal - area)*(rco-rop)/areal
              if (area .gt. areal) ttt  = rop
            end if
          else
            stop 'error: no element with this type'
          end if
          !write(222,*)ie,iet,area,areal,ttt

          do ii = 1,3
            k = nen3v(ii,ie)
            areak = area / 3.
            ark(k) = ark(k) + areak
            tau(k) = tau(k) + ttt*areak
          end do
        end do

        !-------------------------------------------
        ! Set variables ttauv and stauv
        !-------------------------------------------
        do k = 1,nkn
          tau(k) = tau(k) / ark(k)
          rtd(k) = tau(k)/86400.
          if ( surface ) then
            ttauv(1,k) = 1. / tau(k)
            stauv(1,k) = 1. / (tau(k)*2.)
            if ( xgv(k) .lt. -6.0 ) then	!gibraltar
              do l=2,nlvdi
                ttauv(l,k) = ttauv(1,k)
                stauv(l,k) = ttauv(1,k)
              end do
	    else
              do l=2,nlvdi
                ttauv(l,k) = 0.
                stauv(l,k) = 0.
              end do
            end if
          else
            do l=1,nlvdi
              ttauv(l,k) = 1. / tau(k)
              stauv(l,k) = 1. / (tau(k)*2.)
            end do
          end if
        end do

        !-------------------------------------------
        ! Write output for plotting
        !-------------------------------------------
        call init_output('itmcon','idtcon',ia_out)
        if( has_output(ia_out) ) then
           call open_scalar_file(ia_out,1,1,'rlx')
        end if
        call write_scalar_file(ia_out,41,1,rtd)

        write(6,*) 'relaxation time for nudging set up'

        end

c**********************************************************************

        subroutine ggu_ginevra

	use mod_ts
	use mod_diff_visc_fric
	use mod_hydro_print
	use mod_hydro
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'

	include 'femtime.h'



	integer ks,lmax,l

	integer icall
	save icall
	data icall / 0 /

	if( icall .eq. -1 ) return

	ks = 678

	if( icall .eq. 0 ) then
	  icall = -1
	  if( ks .le. 0 ) return
	  lmax = ilhkv(ks)
	  write(177,*) lmax
	  do l=1,lmax
	    write(177,*) l,hlv(l)
	  end do
	  icall = 1
	end if

	lmax = ilhkv(ks)
	write(185,*) it,lmax,(uprv(l,ks),l=1,lmax)
	write(186,*) it,lmax,(vprv(l,ks),l=1,lmax)
	write(187,*) it,lmax,(tempv(l,ks),l=1,lmax)
	write(188,*) it,lmax,(visv(l,ks),l=1,lmax)
	write(189,*) it,lmax,(difv(l,ks),l=1,lmax)
	write(190,*) it,znv(ks)

	end

c**********************************************************************
c**********************************************************************
c**********************************************************************

        subroutine skadar_debug

c debugging skadar application

	implicit none

	include 'femtime.h'

	integer ie,it0

	it0 = 16917074
	ie = 8964

	it0 = 15551273
	ie = 4426

	if( it .lt. it0 ) return
	return

	call check_set_unit(345)
        call check_elem(ie)
        call check_nodes_in_elem(ie)

	end

c**********************************************************************

        subroutine skadar

	implicit none

	call limit_temp
	call write_meteo

	end

c**********************************************************************

        subroutine write_meteo

	use mod_ts
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'

	include 'femtime.h'


	integer l,k,lmax,ks
	real qs,ta,rh,wb,uw,cc,p

	integer icall
	save icall
	data icall /0/

	integer ipint

	ks = 7741
	if( icall .eq. 0 ) ks = ipint(ks)

	call meteo_get_heat_values(k,qs,ta,rh,wb,uw,cc,p)

	if( mod(it,21600) .eq. 0 ) write(565,*) it,qs,ta,rh,uw,cc
	write(566,*) it,qs,ta,rh,uw,cc
	write(567,*) it,tempv(1,ks)

	icall = 1

	end

c**********************************************************************

        subroutine limit_temp

	use mod_ts
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'

	include 'femtime.h'


	integer l,k,lmax
	real tmin

	tmin = 0.

	do k=1,nkn
	  lmax = ilhkv(k)
	  do l=1,lmax
	    if( tempv(l,k) .lt. tmin ) tempv(l,k) = tmin
	  end do
	end do

	end

c**********************************************************************

	subroutine wet_dry

c time of inundation for theseus

	use mod_geom_dynamic
	use mod_depth
	use mod_layer_thickness
	use mod_ts
	use evgeom
	use basin

	implicit none

	include 'femtime.h'

	logical binit,blast
	integer ie,n,k,ii
	integer idtwrite
	real area,dt,h
	real sedim,smed,s
	real sedim_reed,sedim_salt
	real sfact
	character*80 file

	integer, save, allocatable :: idry(:)
	double precision, save, allocatable :: salt_aver_k(:)
	real, save, allocatable :: salt_aver_e(:)

	real, save :: sedim_save = 0.
	integer, save :: isum = 0
	integer, save :: iu122 = 0
	integer, save :: iu123 = 0
	integer, save :: iu124 = 0
	integer, save :: iu127 = 0
	integer, save :: iu128 = 0

	integer, save :: icall = 0

	real getpar
	integer ifemop

	binit = .false.
	blast = it .eq. itend
	idtwrite = 86400
	idtwrite = 86400*30.5
	smed = 5.

c-----------------------------------------
c first call
c-----------------------------------------

	if( icall .eq. 0 ) then
	  allocate(idry(nel))
	  allocate(salt_aver_k(nkn))
	  allocate(salt_aver_e(nel))

	  idry = 0.
	  salt_aver_k = 0.
	  salt_aver_e = smed

	  area = 0.
	  do ie=1,nel
	    area = area + ev(10,ie)
	  end do
	  area = area * 12.
	  if( iu122 == 0 ) iu122 = ifemop('.122','form','new')
	  write(iu122,*) icall,nel,area
	  write(iu122,*) (12.*ev(10,ie),ie=1,nel)
	  close(iu122)
	  sedim_save = getpar('sedim')	! sedimentation in [mm/y]
	  if( sedim_save .lt. 0. ) then
	    sedim_save = -sedim_save
	    binit = .true.		! initialize hev from file
	  end if
	  isum = 0
	end if

c-----------------------------------------
c accumulate
c-----------------------------------------

	icall = icall + 1

	do ie=1,nel
	  if( iwegv(ie) .ne. 0 ) then
	    idry(ie) = idry(ie) + 1	! changed !!! -> was 1, should be idt
	  end if
	end do

	isum = isum + 1
	do k=1,nkn
	  salt_aver_k(k) = salt_aver_k(k) + saltv(1,k)
	end do

c-----------------------------------------
c write to file
c-----------------------------------------

	if( mod(it,idtwrite) .eq. 0 ) then
	  if( iu124 == 0 ) iu124 = ifemop('.124','form','new')
	  write(iu124,*) icall,nel,it,idtwrite
	  write(iu124,*) (idry(ie),ie=1,nel)

	  do k=1,nkn
	    salt_aver_k(k) = salt_aver_k(k) / isum
	  end do
	  do ie=1,nel
	    s = 0.
	    do ii=1,3
	      k = nen3v(ii,ie)
	      s = s + salt_aver_k(k)
	    end do
	    salt_aver_e(ie) = s / 3.
	  end do
	  if( iu127 == 0 ) iu127 = ifemop('.127','form','new')
	  if( iu128 == 0 ) iu128 = ifemop('.128','form','new')
	  write(iu127,*) icall,nel,it,idtwrite
	  write(iu127,*) (salt_aver_e(ie),ie=1,nel)
	  write(iu128,*) icall,nkn,it,idtwrite
	  write(iu128,*) (salt_aver_k(k),k=1,nkn)
	  do k=1,nkn
	    salt_aver_k(k) = 0.
	  end do
	  isum = 0
	end if

	if( blast ) then
	  if( iu123 == 0 ) iu123 = ifemop('.123','form','new')
	  write(iu123,*) icall,nel,it,0
	  write(iu123,*) (idry(ie),ie=1,nel)
	end if

c-----------------------------------------
c modify depth
c-----------------------------------------

	sedim = sedim_save

	if( sedim .gt. 0. ) then

	call get_timestep(dt)
	sfact = dt * 0.001/(365*86400)			! [m/timestep]
	!sedim = sedim*0.001/(365*86400)		! [m/s]
	!sedim = sedim * dt
	sedim_reed = sedim * sfact
	sedim_salt = sedim * sfact * 0.66

	do ie=1,nel
	  if( iwegv(ie) .eq. 0 ) then	!only for flooded elements
	    h = hdenv(1,ie)
	    s = salt_aver_e(ie)
	    if( h .lt. 0.50 ) then	!only for shallow areas
	      if( s .lt. 5. ) then
		sedim = sedim_reed
	      else if( s .gt. 10. ) then
		sedim = sedim_salt
	      else
		sedim = sedim_reed+(sedim_salt-sedim_reed)*(s-5.)/(10.-5.)
	      end if
	      hev(ie) = hev(ie) - sedim
	    end if
	  end if
	end do

	call adjourn_depth_from_hev
	!call set_last_layer

	end if

c-----------------------------------------
c read or write hev
c-----------------------------------------

	if( binit ) then	! initialize hev
	  write(6,*) 'hev initialized'
	  call getfnm('hsedim',file)
	  if( file == ' ' ) stop 'error stop wet_dry: no init file'
	  call read_in_hev(file)
	  call adjourn_depth_from_hev
	  call set_last_layer
          call setweg(0,n)
          call setznv           ! -> change znv since zenv has changed
	end if

	if( blast ) then	! last time step -> write out
	  call def_make('.120',file)
	  call write_out_hev(file)
	  write(6,*) 'sedim used: ',sedim_save,binit
	end if

c-----------------------------------------
c end of routine
c-----------------------------------------

	end

c**********************************************************************

        subroutine init_ts

	use mod_ts
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'

	include 'femtime.h'


c	integer ndim
c	parameter (ndim=10)	!must be at least the number of layers used
c	real tobs(ndim)
c	real sobs(ndim)
c	data tobs /20.,20.,15.,13.,12.,10.,8.,6.,6.,6./
c	data sobs /35.,35.,35.,36.,36.,36.,37.,37.,37.,37./

	integer ndim
        parameter (ndim=22)     !must be at least the number of layers used
        real tobs(ndim)
        real sobs(ndim)
        data tobs /7.4,7.35,7.3,7.4,7.35,7.2,7.3,7.25,7.25,7.25,
     *  7.15,7.15,7.15,7.,7.55,7.5,7.55,7.45,7.5,7.55,7.55,7.85/
        data sobs /35.278,35.452,35.626,35.8,35.974,36.148,36.2995,
     *  36.451,36.6025,36.754,36.9055,37.057,37.1324,37.2078,37.2832,
     *  37.3586,37.434,37.5094,37.5848,37.6602,37.7356,37.811/

	integer l,k,lmax

	integer icall
	save icall
	data icall /0/

	if( icall .ne. 0 ) return
	icall = 1

        write(6,*) '========================================='
        write(6,*) '========================================='
        write(6,*) '========================================='
        write(6,*) 'initializing T/S data with profile'
        write(6,*) '========================================='
        write(6,*) '========================================='
        write(6,*) '========================================='

	do k=1,nkn
	  lmax = ilhkv(k)
	  if( lmax .gt. ndim ) goto 99
	  do l=1,lmax
	    tempv(l,k) = tobs(l)
	    saltv(l,k) = sobs(l)
	  end do
	end do

	return
   99	continue
	write(6,*) 'lmax,ndim: ',lmax,ndim
	stop 'error stop init_ts: ndim must be at least lmax'
	end

c**********************************************************************

	subroutine fluid_mud

	use mod_meteo
	use mod_fluidmud
	use levels, only : nlvdi,nlv
	use basin

	implicit none

	include 'param.h'

	integer k,l,lmax

	integer icall
	save icall
	data icall / 0 /

	if( icall .eq. 10 ) then
	  write(6,*) '********** fluid mud concentration initialized'
	  do k=1,nkn
	    lmax = nlvdi
	    do l=1,lmax
	      mudc(l,k) = l*10.
	      if( 2*(l/2) .eq. l ) mudc(l,k) = 0.
	    end do
	  end do
	end if

	icall = icall + 1

        do k = 1, nkn
          if (xgv(k).gt.1..and.xgv(k).lt.2.) then
            if (ygv(k).gt.1..and.ygv(k).lt.2.) then
              tauxnv(k) = 0.1 
            endif
          endif
        end do

	end

c**********************************************************************

	subroutine test_par

	implicit none

	integer iunit

	integer icall
	save icall
	data icall /0/

	if( icall .ne. 0 ) return
	icall = 1

	iunit = 99

	call pripar(iunit)
	call prifnm(iunit)

	end

c**********************************************************************

	subroutine conz_decay_curonian

	use mod_conz
	use levels
	use basin

	implicit none

	include 'param.h'



	integer k,l,lmax
	real tau1,tau2,aux1,aux2,aux
	real dt
	real x1,y1,x2,y2,x,y
	real qx,qy,d,d0
	logical is_north

c sim 5:	+10 -10 0.5
c sim 6:	+8 -8 0.4

	tau1 = +8.0		!growth rate - efolding time in days
	tau2 = -8.0		!decay rate - efolding time in days
	d0 = 0.4		!distance at which to apply tau1/tau2

	x1 = 20.99
	y1 = 55.3
	x2 = 21.2
	y2 = 55.18

	call get_timestep(dt)
	aux1 = exp(dt/(tau1*86400))
	aux2 = exp(dt/(tau2*86400))

        do k=1,nkn
          lmax = ilhkv(k)
	  x = xgv(k)
	  y = ygv(k)
          do l=1,lmax
	    call dist_point_to_line(x,y,x1,y1,x2,y2,qx,qy,d)
	    !d = 0.	! no decay
	    !d = d0	! constant decay - use above values
	    if( is_north(x,y,x1,y1,x2,y2) ) then
	      aux = 1.+(aux2-1.)*d/d0
              cnv(l,k) = aux * cnv(l,k)
	    else
	      aux = 1.+(aux1-1.)*d/d0
              cnv(l,k) = aux * cnv(l,k)
	    end if
          end do
        end do

	end

c**********************************************************************

	function is_north(x,y,x1,y1,x2,y2)

	implicit none

	logical is_north
	real x,y,x1,y1,x2,y2

	logical left

	is_north = left(x1,y1,x2,y2,x,y)

	end

c**********************************************************************

        subroutine set_yaron

c momentum input for yaron

        use mod_internal
        use mod_layer_thickness
        use evgeom
        use levels
        use basin

        implicit none

        integer kin,lin,ie,ii,k,lmax,nelem
        real rnx,rny,rfact,q,area,h,fact

        kin = 3935
        kin = 2088
        kin = 0
        lin = 8
        nelem = 6
        nelem = 4
        rnx = -1
        rny = 0.
        rfact = 1.1
        q = 10.

        if( kin .le. 0 ) return

        do ie=1,nel
          lmax = ilhv(ie)
          area = 12. * ev(10,ie)
          do ii=1,3
            k = nen3v(ii,ie)
            if( k .eq. kin .and. lmax .le. lin ) then
              h = hdeov(lin,ie)
              fact = rfact * q*q / (h*area*sqrt(area)*nelem)
              write(17,*) 'yaron: ',ie,k,fact,fxv(lin,ie)
              fxv(lin,ie) = fxv(lin,ie) - fact * rnx
              fyv(lin,ie) = fyv(lin,ie) - fact * rny
            end if
          end do
        end do

        end

c*******************************************************************

        subroutine lock_exchange

c lock exchange experiment to test non-hydrostatic dynamics

	use mod_hydro
	use mod_ts
	use mod_layer_thickness
        use levels
        use basin
	use mod_diff_visc_fric

        implicit none

	integer k,l,lmax
        integer icall
        save icall
        data icall / 0 /
        real salval /1.3592/

        if( icall .gt. 0 ) return
        icall = icall+1

	do k=1,nkn
	  lmax = ilhkv(k)
	  do l=1,lmax
	    if(ipv(k).lt.3259)then
	      saltv(l,k) = 0.
	      znv(k) = 0.
	    else
	      saltv(l,k) = salval
	      znv(k) = 0.	   
	    endif
	  enddo
	enddo

	end

c*******************************************************************

	subroutine mpi_test_basin(mode)

c mode == 1	initialize for wind
c mode == 2	initialize for zeta (oscillation)
c mode == 3	initialize only output
c mode == 4	initialize only output for mmba

	use mod_meteo
	use basin
	use levels
	use mod_hydro
	use mod_hydro_print
	use shympi
        use mod_conz

	implicit none

	integer mode

	logical :: bwind
	logical :: bzeta
	logical :: bconz
	logical :: bconzinit
	logical :: breinit
	integer k,ie,ii,i,kk
        integer lsup,lbot
	real z,dz,y,dy,pi,dc,x,dx
	double precision dtime
	real, save :: cd = 2.5e-3
	real, save :: wind = 10.
	real, save :: wfact = 1.025/1000.
	real :: stress
	integer, save :: icall = 0

	integer, parameter :: ndim = 5
	integer, save :: nodes(ndim) = (/5,59,113,167,221/)
	integer, parameter :: nmdim = 9
	integer, save :: nmodes(nmdim) = (/132539,131111,132453
     +		,202,201,204,203,129327,130459/)

	integer ipint

!--------------------------------------------------------------
! preliminary
!--------------------------------------------------------------

	if( mode < 1 .or. mode > 4 ) then
	 stop 'error stop mpi_test_basin: mode'
	end if
	bwind = mode == 1
	bzeta = mode == 2

	if( mode <= 3 ) then
	if( nkn_global /= 225 .or. nel_global /= 384 ) then
	  write(6,*) 'nkn,nel: ',nkn,nel
	  stop 'error stop mpi_test_basin: basin'
	end if
	end if

	call get_act_dtime(dtime)

        bconz = mod_conz_is_initialized()
        breinit = ( dtime == 43200. )
	breinit = .false.
	bconzinit = .false.

	dz = 0.1
	dy = 3000.
	dx = 1000.
	dc = 10.
	pi = 4.*atan(1.)

!--------------------------------------------------------------
! initialization
!--------------------------------------------------------------

	if( icall .eq. 0 ) then

	  write(6,*) 'initializing basin for MPI: ',bwind,bzeta,bconz

	  if( bwind ) then
	    stress = wfact * cd * wind * wind
	    wxv = 0.
	    wyv = wind
	    metws = wind
	    ppv = 0.
	    tauxnv = 0.
	    tauynv = stress
	  else if( bzeta ) then
	    do k=1,nkn
	      y = ygv(k) - dy
	      z = dz * (y/dy)
	      z = dz * sin( (pi/2.) * (y/dy) )
	      znv(k) = z
	    end do
	    call setzev
	    call make_new_depth
          end if

          if( bconzinit ) then
            write(6,*) 'initializing concentration from subcus...'
	    do k=1,nkn
              x = xgv(k)
	      y = ygv(k) - dy
	      z = dc * (1.+y/dy) + dc * (1.+x/dx)
              z = 10.
              if( y/dy > 0 ) z = 20.
	      cnv(:,k) = z
	    end do
	  end if

          write(6,*) 'initializing node output from subcus...'
	  if( mode <= 3 ) then
	    call register_nodes(ndim,nodes)
	  else
	    call register_nodes(nmdim,nmodes)
	  end if

	end if

	icall = icall + 1

!--------------------------------------------------------------
! re-initialization
!--------------------------------------------------------------

        if( bconz .and. breinit ) then
          write(6,*) 're-initializing concentration from subcus...'
	  do k=1,nkn
            x = xgv(k)
	    y = ygv(k) - dy
	    z = dc * (1.+y/dy) + dc * (1.+x/dx)
            z = 10.
            if( y/dy > 0 ) z = 20.
	    cnv(:,k) = z
	  end do
	end if

!--------------------------------------------------------------
! output results
!--------------------------------------------------------------

        !write(6,*) 'nodes: ',ndim,nodes
        !flush(6)

        lsup = min(2,nlv)
        lbot = max(nlv-1,1)

	if( mode <= 3 ) then
	do i=1,ndim
	  k = nodes(i)
	  if( k .le. 0 ) cycle
	  write(500+i,*) dtime,znv(k)
	  write(600+i,*) dtime,vprv(lsup,k)
	  write(700+i,*) dtime,vprv(lbot,k)
          if( bconz ) then
	  write(800+i,*) dtime,cnv(lsup,k)
	  write(900+i,*) dtime,cnv(lbot,k)
          end if
	end do
	else
	do i=1,nmdim
	  k = nmodes(i)
	  if( k .le. 0 ) cycle
	  write(500+i,*) dtime,znv(k)
	end do
	end if
	
!--------------------------------------------------------------
! end of routine
!--------------------------------------------------------------

	end

c*******************************************************************

	subroutine register_nodes(n,nodes)

	use shympi

	implicit none

	integer n
	integer nodes(n)

	integer i,k,kk
	integer ipint

	  do i=1,n
	    k = nodes(i)
	    kk = ipint(k)
	    if( kk .le. 0 ) then
	      write(6,*) '**** ignoring not existing node ',k,kk,my_id
	      !stop 'error stop mpi_test_basin: no such node'
	    else if( .not. shympi_is_inner_node(kk) ) then
	      write(6,*) '**** ignoring ghost node ',k,kk,my_id
	      kk = 0
	    else
	      write(6,*) '**** register node ',k,kk,my_id
	    end if
	    nodes(i) = kk
	  end do

	end

c*******************************************************************

	subroutine cyano_diana

! write out values for cyanobacteria scum simulations

	use basin
	use mod_waves
	use mod_meteo
	use mod_hydro_print

	implicit none

	integer, save :: nvar = 5			!number of vars written
	character*80, save :: file = 'cdiana.txt'	!file name of dates
	double precision, save :: dtout = 1800.		!output frequence
	double precision, save :: dtwin = 3.*3600.	!output window

	character*20 :: aline
	integer ierr,k
	integer, save :: id = 0
	logical, save :: bfinish = .false.
	double precision :: dtime,atime
	double precision, save :: astart,anext
	double precision, save :: aend = 0		!impossible time
	real s,d
	real uvmed(nkn),uvdir(nkn),windir(nkn)

	if( bfinish ) return	!nothing more to read

        call get_absolute_act_time(atime)

	if( aend < atime ) then	!read new output time
	  ierr = 0
	  call get_new_cyano_time(file,dtout,dtwin,astart,aend,ierr)
	  if( ierr /= 0 ) bfinish = .true.
	  if( bfinish ) return	!nothing more to read
	  anext = astart
	end if

	if( atime < anext ) return	!no output yet
	anext = anext + dtout

	if( id == 0 ) then
          call shyfem_init_scalar_file('cyano',nvar,.true.,id)
	end if

	do k=1,nkn
	  call convert_wind_sd(uprv(1,k),vprv(1,k),s,d)
	  uvmed(k) = s
	  d = d + 180.
	  if( d > 360. ) d = d - 360.
	  uvdir(k) = d
	  call convert_wind_sd(wxv,wyv,s,d)
	  windir(k) = d
	end do
	!uvmed = sqrt( uprv(1,:)**2 + vprv(1,:)**2 )

	call get_act_dtime(dtime)
	call get_act_timeline(aline)
	write(6,*) 'new cyano output written: ',aline
        call shy_write_scalar_record(id,dtime,231,1,waveh)
        call shy_write_scalar_record(id,dtime,28,1,metws)
        call shy_write_scalar_record(id,dtime,29,1,windir)
        call shy_write_scalar_record(id,dtime,6,1,uvmed)
        call shy_write_scalar_record(id,dtime,7,1,uvdir)

	end

c*******************************************************************

	subroutine get_new_cyano_time(file,dtout,dtwin,astart,aend,ierr)

! reads date and determines start and end date for output

	use iso8601

	implicit none

	character*(*) file
	double precision :: dtout,dtwin,astart,aend
	integer ierr

	integer date,time
	integer ios
	integer, save :: iud = 0
	double precision atime,aout
	character*80 line
	character*20 aline

	if( iud == 0 ) then
	  open(iud,file=file,status='old',form='formatted',iostat=ios)
	  if( ios /= 0 ) then
	    write(6,*) 'cannot open file: ',trim(file)
	    stop 'error stop get_new_cyano_time: error opening file'
	  end if
	end if

	ierr = -1
	read(iud,'(a)',iostat=ios) line
	if( ios > 0 ) then
	  stop 'error stop get_new_cyano_time: read error'
	end if
	if( ios < 0 ) return

	call string2date(line,date,time,ierr)
	call dts_to_abs_time(date,time,aout)

	atime = dnint(aout/dtout)
	atime = atime * dtout
	astart = atime - dtwin
	aend = atime + dtwin
	  
	write(6,*) '------------------------------------------'
	call dts_format_abs_time(aout,aline)
	write(6,*) 'new output date for cyano read: ',aline
	call dts_format_abs_time(atime,aline)
	write(6,*) 'central output date for cyano: ',aline
	call dts_format_abs_time(astart,aline)
	write(6,*) 'start output date for cyano: ',aline
	call dts_format_abs_time(aend,aline)
	write(6,*) 'end output date for cyano: ',aline
	write(6,*) '------------------------------------------'

	end

c*******************************************************************

