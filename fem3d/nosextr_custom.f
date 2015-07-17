c
c $Id: nosextr_custom.f,v 1.2 2009-04-07 10:43:57 georg Exp $
c
c revision log :
c
c 18.11.1998    ggu     check dimensions with dimnos
c 24.02.1999    ggu     use n2int for node number translation
c 03.12.2001    ggu     cleaned up, hakata bay
c
c****************************************************************

	program nosextr_custom

c extracts single nodes from nos file -> creates time series

	use basin !COMMON_GGU_SUBST

	implicit none

	include 'param.h'

c--------------------------------------------------
COMMON_GGU_DELETED	include 'basin.h'


c--------------------------------------------------

	character*80 title
	real cv(nkndim)
	real cv3(nlvdim,nkndim)

	integer ilhkv(nkndim)
	real hlv(nlvdim)
	real hev(neldim)

	logical berror
        integer i,n,k,ke,l
        integer nread,nunit
        integer nvers
        integer nlv,nvar,ivar,ierr
        integer nin,it

	integer iapini,ideffi,ialfa

c--------------------------------------------------
	integer nudim
	parameter(nudim=300)
	integer iunit(nudim)
	character*50 file

c---------------------------------------------------------------
c nodes for extraction
c---------------------------------------------------------------

	integer ndim
c	parameter( ndim = 11 )          !taranto
c	parameter( ndim = 21 )
c	parameter( ndim = 4 )           !curonian
c	parameter( ndim = 5 )
c	parameter( ndim = 27 )          !francesca
c	parameter( ndim = 7 )           !edison
	parameter( ndim = 3 )           !parallel

	integer nodese(ndim)	!external numbers
	integer nodes(ndim)	!internal numbers

c	data nodese /1601,1762,1216,1523,726
c     +                  ,1075,534,1043,1032,747,18,1256,959
c     +                  ,258,1204,880,1280,45,8,135,1458/
c	data nodese /938,1483,2104,3306,3226,2162,859/
c	data nodese /328,1970,2784,2683,2435/
c        data nodese /1017,726,3133,4080,464,
c     +    308,354,847,2934,1813,3494/ !extr nodes staz. CTD Taranto
c        data nodese /526,1026,1526,2026,2526,
c     +               3026,3526,4026,4526/
c        data nodese /1316,1843,1623,1871/ !curonian
c        data nodese /1512,484,928,326,1241,
c     +              243,450,384/ !curonian stat
c
c         data nodese /2703,2765,2766,2754,2536,2744,3374,3373,
c     +               3372,3270,3293,3597,3417,3416,3415,3340,
c     +               3296,3601,2767,2770,2773,3369,3371,3377,
c     +               3380,3411,3414/  !nodi francesca
c
c	data nodese /6263,6273,3545,6262,6261,6274,466/	!nodi edison
c	data nodese /6260,6270,6269,6257,6256,3540,6286/	!nodi edison

	data nodese /5561,2527,3604/				!nodi parallel

c---------------------------------------------------------------
c open simulation and basin
c---------------------------------------------------------------

	if(iapini(3,nkndim,neldim,0).eq.0) then
		stop 'error stop : iapini'
	end if

c---------------------------------------------------------------
c initializing units and nodes to be extracted
c---------------------------------------------------------------

	nunit = 60
	do i=1,nudim
	  iunit(i) = 0
	end do

	nread=0

	do i=1,ndim
	  nodes(i) = nodese(i)
	end do

	call n2int(ndim,nodes,berror)

	if( berror ) stop 'error stop nosextr'

	write(6,*) 'Extracting ',ndim,' nodes :'
	write(6,*) (nodese(i),i=1,ndim)

c---------------------------------------------------------------
c open files and read headers
c---------------------------------------------------------------

	nin=ideffi('datdir','runnam','.nos','unform','old')
	if(nin.le.0) goto 100

        nvers=3
	call rfnos(nin,nvers,nkn,nel,nlv,nvar,title,ierr)
        if(ierr.ne.0) goto 100

        write(6,*) 'nvers    : ',nvers
        write(6,*) 'nkn,nel  : ',nkn,nel
        write(6,*) 'nlv,nvar : ',nlv,nvar
        write(6,*) 'title    : ',title

        call dimnos(nin,nkndim,neldim,nlvdim)

	call rsnos(nin,ilhkv,hlv,hev,ierr)
        if(ierr.ne.0) goto 100

	write(6,*) 'Available levels: ',nlv
	write(6,*) (hlv(l),l=1,nlv)

c---------------------------------------------------------------
c time loop
c---------------------------------------------------------------

  300   continue

	call rdnos(nin,it,ivar,nlvdim,ilhkv,cv3,ierr)

        if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
        if(ierr.ne.0) goto 100

	nread=nread+1
	write(6,*) 'time : ',it,ivar

	if( ivar .gt. nudim ) then              !ivar too high
          write(6,*) ivar,nudim
          stop 'error stop: nudim'
        else if( iunit(ivar) .eq. 0 ) then      !not yet initialized
	  nunit = nunit + 1
	  file = ' '
	  n = ialfa(float(ivar),file,-1,-1)
	  file(n+1:) = '.dat'
	  write(6,*) file
	  open(nunit,file=file,status='unknown',form='formatted')
	  iunit(ivar) = nunit
c	  write(nunit,'(2i10)') 0,ndim
c	  write(nunit,'(6i12)') (nodese(i),i=1,ndim)
	else                                    !already initialized
	  nunit = iunit(ivar)
	end if

	!write(nunit,'(i10)') it
	!write(nunit,'(6e12.4)') (cv3(1,nodes(i)),i=1,ndim)
        write(nunit,'(i10,30e12.4)') it,(cv3(1,nodes(i)),i=1,ndim)

	do i=1,ndim
	  k = nodes(i)
	  ke = nodese(i)
	  write(79,*) it,ke,ilhkv(k),ivar,k,i
	  write(79,*) (cv3(l,k),l=1,ilhkv(k))
	end do

	goto 300

c---------------------------------------------------------------
c end of time loop
c---------------------------------------------------------------

  100	continue

	write(6,*)
	write(6,*) nread,' records read'
	write(6,*)
	write(6,*) 'data written to file 79'
	write(6,*)

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

	end

c***************************************************************

