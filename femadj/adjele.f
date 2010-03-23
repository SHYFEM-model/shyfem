c
c $Id: adjele.f,v 1.21 2010-03-11 15:31:28 georg Exp $
c
c description :
c
c main for mesh adjustment
c
c contents :
c
c adjele			main
c
c revision log :
c
c 19.05.2003    ggu     some more checks
c 02.12.2004    ggu     some more debug checks (version 1.41)
c 02.12.2004    ggu     copyright statement added
c 20.03.2007    ggu     fake version
c 16.04.2008    ggu     new Makefile structure
c 09.12.2008    ggu     small changes to Makefile
c 26.01.2009    ggu     makedepend, compiler warnings
c 06.04.2009    ggu     new param.h, include basin.h, better error handling
c 24.04.2009    ggu     new call to rdgrd()
c
c***************************************************************

	program adjele

c adjusts elements after automatic mesh generator

	implicit none

        character*20 version
        parameter (version='1.56')

c------------------------------------------------------- parameters
	include 'param.h'
c------------------------------------------------------- grade index
	include 'grade.h'
c------------------------------------------------------- basin
	include 'basin.h'
c------------------------------------------------------- local
	integer kspecial
	integer nlidim,nlndim
	integer nmax,k
	integer ner
	integer nco,nknh,nli
	logical bstop, bplot
	character*80 file
	real dx(nkndim),dy(nkndim)
	integer ic(nkndim)
	integer ianv(nkndim)
	integer ipaux(nkndim)
	integer iaux(neldim)
	real raux(neldim)
c------------------------------------------------------- end declaration

	kspecial = 0		!set to 0 for no debug

	file = 'old.grd'
	file = ' '
	ner = 6
	bstop = .false.
	bplot = .false.		!for plotting of basin and grades
	bplot = .true.		!for plotting of basin and grades

c-------------------------------------------------------------------

        call copyright(version)

c-------------------------------------------------------------------

c read grid file with nodes and elements

        nlidim = 0
        nlndim = 0
        call rdgrd(
     +                   file
     +                  ,bstop
     +                  ,nco,nkn,nel,nli
     +                  ,nkndim,neldim,nlidim,nlndim
     +                  ,ipv,ipev,iaux
     +                  ,ianv,iarv,iaux
     +                  ,hkv,hev,raux
     +                  ,xgv,ygv
     +                  ,nen3v
     +                  ,iaux,iaux
     +                  )

        if( bstop ) goto 97

	write(6,'(a,5i7)') 'grid parameters: ',nkn,nel,nli

c external to internal numbering

        call ex2in(nkn,3*nel,nlidim,ipv,ipaux,nen3v,iaux,bstop)
        if( bstop ) goto 99

c save depth information in elements to nodes

	nknh = 0
	do k=1,nkn
	  if( hkv(k) .ne. -999. ) nknh = nknh + 1
	end do

	if( nknh .eq. 0 ) then
	  write(6,*) 'copying element depth to nodes...'
	  call hev2hkv(nkn,nel,nen3v,hev,hkv,ic)
	end if

c determine grade

	call maxgrd(nkn,nel,nen3v,ngrade,nmax)

	write(6,*) 'maximum grade: ',nmax
	if( nmax .gt. ngrdim ) goto 98

	call statgrd('first call',nkn,nmax,ngrade,nbound)

c make boundary nodes (flag nbound)

	call mkbound(nkn,nel,ngrdim,nen3v,ngrade,nbound,ngri)
        call mkstatic(nkn,ianv,nbound)
	call statgrd('boundary nodes',nkn,nmax,ngrade,nbound)

	call nodeinfo(kspecial,nkn,nel,ngrdim,ngri)

c plot grade

	call qopen

	if( bplot) call plobas(nkn,nel,nen3v,ngrade,xgv,ygv)

c eliminate 4- grades

        write(6,*) 'start eliminating nodes ...'

	call chkgrd
        write(6,*) 'chkgrd ok ...'
c	call pltsgrd(4,nkn,nel,nen3v,ngrade,xgv,ygv)
	call eliml(nkn,nel,ngrdim,ngrade,nbound,ngri,nen3v)
	if( bplot) call plobas(nkn,nel,nen3v,ngrade,xgv,ygv)
	call statgrd('4- grades',nkn,nmax,ngrade,nbound)
	call nodeinfo(kspecial,nkn,nel,ngrdim,ngri)

c eliminate 8+ grades

	call chkgrd
	call elimh(8,nkn,nel,ngrdim,ngrade,nbound,ngri,nen3v)
	if( bplot) call plobas(nkn,nel,nen3v,ngrade,xgv,ygv)
	call statgrd('8+ grades',nkn,nmax,ngrade,nbound)

	write(6,*) 'checking after 8+'
	call chkgrd
	call nodeinfo(kspecial,nkn,nel,ngrdim,ngri)

c eliminate 7+ grades

	call chkgrd
	call elimh(7,nkn,nel,ngrdim,ngrade,nbound,ngri,nen3v)
	if( bplot) call plobas(nkn,nel,nen3v,ngrade,xgv,ygv)
	call statgrd('7+ grades',nkn,nmax,ngrade,nbound)

	write(6,*) 'checking after 7+'
	call chkgrd
	call nodeinfo(kspecial,nkn,nel,ngrdim,ngri)

c smoothing

	call wrgrd('new_nosmooth.grd',nkn,nel,xgv,ygv,nen3v)

        call smooth(50,0.1,nkn,nel
     +                          ,nen3v,nbound
     +                          ,xgv,ygv,dx,dy,ic)
	if( bplot) call plobas(nkn,nel,nen3v,ngrade,xgv,ygv)
	call nodeinfo(kspecial,nkn,nel,ngrdim,ngri)

c again ...

	call chkgrd
        call eliml(nkn,nel,ngrdim,ngrade,nbound,ngri,nen3v)
	call elimh(8,nkn,nel,ngrdim,ngrade,nbound,ngri,nen3v)
	call elimh(7,nkn,nel,ngrdim,ngrade,nbound,ngri,nen3v)
	if( bplot) call plobas(nkn,nel,nen3v,ngrade,xgv,ygv)
	call statgrd('one more time',nkn,nmax,ngrade,nbound)
	call nodeinfo(kspecial,nkn,nel,ngrdim,ngri)

c eliminate 5 grades

	call wrgrd('new_help.grd',nkn,nel,xgv,ygv,nen3v)

	call chkgrd
	call elim5(nkn,nel,ngrdim,ngrade,nbound,ngri,nen3v)
	if( bplot) call plobas(nkn,nel,nen3v,ngrade,xgv,ygv)
	call statgrd('5 grades',nkn,nmax,ngrade,nbound)

	write(6,*) 'checking after 5'
	call chkgrd
	call nodeinfo(kspecial,nkn,nel,ngrdim,ngri)

c eliminate 5-5 grades

	call elim57(nkn,nel,ngrdim,ngrade,nbound,ngri,nen3v)
	if( bplot) call plobas(nkn,nel,nen3v,ngrade,xgv,ygv)
	call statgrd('5-5 grades',nkn,nmax,ngrade,nbound)

	write(6,*) 'checking after 5-5'
	call chkgrd
	call nodeinfo(kspecial,nkn,nel,ngrdim,ngri)

c one more time

	write(6,*) 'one more round...'
	call elimh(8,nkn,nel,ngrdim,ngrade,nbound,ngri,nen3v)
        call elimh(7,nkn,nel,ngrdim,ngrade,nbound,ngri,nen3v)
        call elim5(nkn,nel,ngrdim,ngrade,nbound,ngri,nen3v)
        call elim57(nkn,nel,ngrdim,ngrade,nbound,ngri,nen3v)
        call statgrd('all again',nkn,nmax,ngrade,nbound)
	call chkgrd
	call nodeinfo(kspecial,nkn,nel,ngrdim,ngri)

c smoothing

        call smooth(50,0.1,nkn,nel
     +                          ,nen3v,nbound
     +                          ,xgv,ygv,dx,dy,ic)
	if( bplot) call plobas(nkn,nel,nen3v,ngrade,xgv,ygv)
	call nodeinfo(kspecial,nkn,nel,ngrdim,ngri)

c write to grd file

	call chkgrd
	call nodeinfo(kspecial,nkn,nel,ngrdim,ngri)
	call wrgrd('new.grd',nkn,nel,xgv,ygv,nen3v)

	call qclose

	stop
   97	continue
	write(6,*) 'error reading grd file'
	stop 'error stop rdgrd'
   98	continue
	write(6,*) 'Grade is too high. Please increase ngrdim.'
	write(6,*) '(you can set it in Rules.make)'
	stop 'error stop: ngrdim'
   99	continue
	write(6,*) 'error external to internal numbering'
	stop 'error stop ex2in'
	end

c***********************************************************

	subroutine hev2hkv(nkn,nel,nen3v,hev,hkv,ic)

c saves information about depth to nodes

	implicit none

	integer nkn,nel
	integer nen3v(3,1)
	real hev(1),hkv(1)
	integer ic(1)

	integer k,ie,ii

	do k=1,nkn
	  hkv(k) = 0.
	  ic(k) = 0
	end do

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    hkv(k) = hkv(k) + hev(ie)
	    ic(k) = ic(k) + 1
	  end do
	end do

	do k=1,nkn
	  if( ic(k) .gt. 0 ) then
	    hkv(k) = hkv(k) / ic(k)
	  end if
	end do

	end

c***********************************************************

        subroutine copyright(version)

c writes copyright 

        implicit none

        character*(*) version


        write(6,*)
        write(6,*) ' ----------------------------------------'
        write(6,*)
        write(6,*) ' ADJELE - Regularize finite element grids'
        write(6,*) ' Copyright (c)  Georg Umgiesser 1985-2010'
        write(6,*)
        write(6,*) ' version ',version
        write(6,*)
        write(6,*) ' ----------------------------------------'
        write(6,*)

        end

c***********************************************************

