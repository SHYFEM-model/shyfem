c
c $Id: nosextr_nodes.f,v 1.3 2009-04-07 10:43:57 georg Exp $
c
c extract nodes from NOS file
c
c revision log :
c
c 18.11.1998    ggu     check dimensions with dimnos
c 24.02.1999    ggu     use n2int for node number translation
c 03.12.2001    ggu     cleaned up, hakata bay
c 03.06.2011    ggu     routine adjourned
c 25.01.2013    ggu     routines cleaned
c 17.05.2013    ggu     prepared for write_profile
c 05.09.2013    ggu     new call to get_layer_thickness()
c
c***************************************************************

	program nosextr_nodes

c extracts single nodes from NOS file -> creates time series

	use mod_depth
	use evgeom
	use levels
	use basin

	implicit none

	include 'param.h'

c--------------------------------------------------

	character*80 title



	real, allocatable :: hl(:)
	real, allocatable :: cv(:)
	real, allocatable :: cv3(:,:)

	character*80 format
	integer nread,ierr
	integer nvers,nin
	integer nknnos,nelnos

	integer it
	integer k,ke,i
	integer l,lmax
	real z,h

	integer nvar,ivar,iu

	integer iapini,ideffi

c---------------------------------------------------------------
c nodes for extraction
c---------------------------------------------------------------

	integer ndim
	integer nnodes
	integer, allocatable :: nodes(:)	!node numbers
	integer, allocatable :: nodese(:)	!external node numbers

c---------------------------------------------------------------
c initialize params
c---------------------------------------------------------------

        nread=0

c---------------------------------------------------------------
c open simulation and basin
c---------------------------------------------------------------

	if(iapini(3,0,0,0).eq.0) then
		stop 'error stop : iapini'
	end if

	call ev_init(nel)
	call mod_depth_init(nkn,nel)

	call set_ev

c---------------------------------------------------------------
c open NOS file and read header
c---------------------------------------------------------------

	nin=ideffi('datdir','runnam','.nos','unform','old')
	if(nin.le.0) goto 100

        nvers=3
	call rfnos(nin,nvers,nknnos,nelnos,nlv,nvar,title,ierr)
        if(ierr.ne.0) goto 100

	write(6,*)
        write(6,*) 'title        : ',title
	write(6,*)
        write(6,*) 'nvers        : ',nvers
        write(6,*) 'nkn,nel      : ',nknnos,nelnos
        write(6,*) 'nlv          : ',nlv
        write(6,*) 'nvar         : ',nvar
	write(6,*)

	if( nkn .ne. nknnos .or. nel .ne. nelnos ) goto 94
	nlvdi = nlv

	allocate(hl(nlvdi))
	allocate(cv(nkn))
	allocate(cv3(nlvdi,nkn))

        call dimnos(nin,nkn,nel,nlvdi)
	call levels_init(nkn,nel,nlv)

	call rsnos(nin,ilhkv,hlv,hev,ierr)
        if(ierr.ne.0) goto 100

        call init_sigma_info(nlv,hlv)
        call makehev(hev)
	call makehkv_minmax(hkv,+1)

	write(6,*) 'Available levels: ',nlv
	write(6,*) (hlv(l),l=1,nlv)

c---------------------------------------------------------------
c get nodes to extract from STDIN
c---------------------------------------------------------------

	ndim = nkn
	allocate(nodes(ndim))
	allocate(nodese(ndim))
        call get_nodes_from_stdin(ndim,nnodes,nodes,nodese)

	if( nnodes .le. 0 ) goto 95

	write(6,*) nnodes
	write(6,*) (nodes(i),i=1,nnodes)
	write(6,*) (nodese(i),i=1,nnodes)

c---------------------------------------------------------------
c loop on input records
c---------------------------------------------------------------

  300   continue

	call rdnos(nin,it,ivar,nlvdi,ilhkv,cv3,ierr)

        if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
        if(ierr.ne.0) goto 100

	nread=nread+1
	write(6,*) 'time : ',it,nread,ivar

	call get_unit(ivar,iu)		!gets unit number for variable

c	---------------------------------------------------------
c	write to file
c	---------------------------------------------------------

	do i=1,nnodes
	  k = nodes(i)
	  ke = nodese(i)
	  lmax = ilhkv(k)
	  write(4,*) it,i,ke,k,lmax,ivar
	  write(4,*) (cv3(l,k),l=1,lmax)
	  write(3,*) it,i,ke,k,lmax,ivar
          write(format,'(a,i4,a)') '(',lmax,'(f12.4))'
          write(3,format) (cv3(l,k),l=1,lmax)
	  !write(3,'((6f10.2))') (cv3(l,k),l=1,lmax)
	  !write(3,'(*(f10.2))') (cv3(l,k),l=1,lmax)

	  z = 0.
	  h = hkv(k)
	  call write_profile_c(it,i,k,ke,lmax,ivar,h,z
     +				,cv3(1,k),hlv,hl)
	end do

        write(format,'(a,i4,a)') '(i12,',nnodes,'(f12.4))'
        write(iu,format) it,(cv3(1,nodes(i)),i=1,nnodes)
        !write(iu,'(i10,30e12.4)') it,(cv3(1,nodes(i)),i=1,nnodes)
        !write(iu,'(i10,*(e12.4))') it,(cv3(1,nodes(i)),i=1,nnodes)

	goto 300

  100	continue

c---------------------------------------------------------------
c end of loop
c---------------------------------------------------------------

	write(6,*)
	write(6,*) nread,' records read'
	write(6,*)
	write(6,*) 'data written to file 3 and 4'
	write(6,*)

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

	stop
   95   continue
        write(6,*) 'no nodes read'
        stop 'error stop nosextr_nodes: no nodes'
   94   continue
        write(6,*) 'incompatible simulation and basin'
        write(6,*) 'nkn: ',nkn,nknnos
        write(6,*) 'nel: ',nel,nelnos
        stop 'error stop nosextr_nodes: nkn,nel'
	end

c***************************************************************

	subroutine get_unit(ivar,iu)

c handles unit numbers for files of single variables

	implicit none

	integer ivar		!actual variable number
	integer iu		!unit number to be used (return)

	integer nudim
	parameter(nudim=300)
	integer iunit(nudim)
	save iunit

	integer n,i
	character*50 file
	integer ialfa

	integer nunit
	save nunit
	data nunit / 0 /

	if( nunit .eq. 0 ) then
	  nunit = 60
	  do i=1,nudim
	    iunit(i) = 0
	  end do
	end if

	if( ivar .gt. nudim ) then              !ivar too high
          write(6,*) 'nudim is too low ',ivar,nudim
          write(6,*) 'please increase nudim'
          stop 'error stop get_unit: nudim'
	end if

        if( iunit(ivar) .eq. 0 ) then      !not yet initialized
	  nunit = nunit + 1
	  file = ' '
	  n = ialfa(float(ivar),file,-1,-1)
	  file(n+1:) = '.dat'
	  write(6,*) 'opening file ',file
	  open(nunit,file=file,status='unknown',form='formatted')
	  iunit(ivar) = nunit
	end if

	iu = iunit(ivar)

	end

c***************************************************************

