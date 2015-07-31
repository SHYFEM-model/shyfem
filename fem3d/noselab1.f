c
c $Id: noselab.f,v 1.8 2008-11-20 10:51:34 georg Exp $
c
c revision log :
c
c 18.11.1998    ggu     check dimensions with dimnos
c 06.04.1999    ggu     some cosmetic changes
c 03.12.2001    ggu     some extra output -> place of min/max
c 09.12.2003    ggu     check for NaN introduced
c 07.03.2007    ggu     easier call
c 08.11.2008    ggu     do not compute min/max in non-existing layers
c 07.12.2010    ggu     write statistics on depth distribution (depth_stats)
c 06.05.2015    ggu     noselab started
c 05.06.2015    ggu     many more features added
c
c**************************************************************

	subroutine noselab

	use clo
	use elabutil

        use basin
        use mod_depth
        use evgeom
        use levels

c elaborates nos file

	implicit none

	include 'param.h'

	integer, parameter :: ndim = 1000
	integer iusplit(ndim)

	real, allocatable :: cv2(:)
	real, allocatable :: cv3(:,:)
	real, allocatable :: cv3all(:,:,:)
	real, allocatable :: vol3(:,:)

	integer, allocatable :: ivars(:)
	integer, allocatable :: nodes(:)

	integer, allocatable :: naccu(:)
	double precision, allocatable :: accum(:,:,:)

	real, allocatable :: hl(:)

	integer nread,nelab,nrec,nin
	integer nvers
	integer nknnos,nelnos,nvar
	integer ierr
	integer it,ivar,itvar,itnew,itold,iaux
	integer i,j,l,k,lmax,nnodes,node
	integer ip,nb,naccum
	character*80 title,name
	character*20 dline
	character*80 basnam,simnam
	real rnull
	real cmin,cmax,cmed,vtot

	integer iapini
	integer ifem_open_file

c--------------------------------------------------------------

	nread=0
	nelab=0
	nrec=0
	rnull=0.
	rnull=-1.
	bopen = .false.

c--------------------------------------------------------------
c set command line parameters
c--------------------------------------------------------------

	call elabutil_init('NOS')

	!--------------------------------------------------------------
	! open input files
	!--------------------------------------------------------------

	call ap_init(bask,modeb,0,0)

	call open_nos_type('.nos','old',nin)

	call peek_nos_header(nin,nknnos,nelnos,nlv,nvar)

	if( bneedbasin ) then
	  if( nkn /= nknnos .or. nel /= nelnos ) goto 92
	else
	  nkn = nknnos
	  nel = nelnos
	end if

        call mod_depth_init(nkn,nel)
        call levels_init(nkn,nel,nlv)

	allocate(cv2(nkn))
	allocate(cv3(nlv,nkn))
	allocate(vol3(nlv,nkn))
	allocate(cv3all(nlv,nkn,nvar))
        allocate(hl(nlv))
	allocate(ivars(nvar))

	nlvdi = nlv
	call read_nos_header(nin,nkn,nel,nlvdi,ilhkv,hlv,hev)
	call nos_get_params(nin,nkn,nel,nlv,nvar)

	call init_sigma_info(nlv,hlv)

	if( bneedbasin ) then
	  call outfile_make_hkv(nkn,nel,nen3v,hev,hkv)
	  call init_volume(nlvdi,nkn,nel,nlv,nen3v,ilhkv
     +                          ,hlv,hev,hl,vol3)
	end if

	if( bverb ) call depth_stats(nkn,nlvdi,ilhkv)

	if( bnodes .or. bnode ) then
	  if( bnodes ) then
	    nnodes = 0
	    call get_node_list(nodefile,nnodes,nodes)
	    allocate(nodes(nnodes))
	    call get_node_list(nodefile,nnodes,nodes)
	  else
	    nnodes = 1
	    allocate(nodes(nnodes))
	    nodes(1) = nodesp
	  end if
	  write(6,*) 'nodes: ',nnodes,(nodes(i),i=1,nnodes)
	  call convert_internal_nodes(nnodes,nodes)
	  bnodes = .true.
	end if

	!--------------------------------------------------------------
	! time management
	!--------------------------------------------------------------

	call nos_get_date(nin,date,time)
	call elabutil_date_and_time

	!--------------------------------------------------------------
	! averaging
	!--------------------------------------------------------------

	call elabutil_set_averaging(nvar)

	if( btrans ) then
	  allocate(naccu(istep))
	  allocate(accum(nlvdi,nkn,istep))
	  naccum = 0
	  naccu = 0
	  accum = 0.
	end if

	!write(6,*) 'mode: ',mode,ifreq,istep

	!--------------------------------------------------------------
	! open output file
	!--------------------------------------------------------------

	iusplit = 0

	boutput = boutput .or. btrans
	bopen = boutput .and. .not. bsplit

	if( bopen ) then
          call open_nos_file('out','new',nb)
          call nos_init(nb,0)
          call nos_clone_params(nin,nb)
	  if( b2d ) then
	    call nos_set_params(nb,0,0,1,0)
	  end if
          call write_nos_header(nb,ilhkv,hlv,hev)
	end if

c--------------------------------------------------------------
c loop on data
c--------------------------------------------------------------

	it = 0
	if( .not. bquiet ) write(6,*)

	cv3 = 0.

	do

	 itold = it

	 do i=1,nvar
	  call nos_read_record(nin,it,ivar,nlvdi,ilhkv,cv3,ierr)
          if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
          if(ierr.ne.0) exit
	  if( i == 1 ) itvar = it
	  if( itvar /= it ) goto 85
	  ivars(i) = ivar
	  cv3all(:,:,i) = cv3(:,:)
	  nread=nread+1
	 end do

         if(ierr.ne.0) exit
	 nrec = nrec + 1

	 if( nrec == 1 ) itold = it
	 call nos_peek_record(nin,itnew,iaux,ierr)
	 !write(6,*) 'peek: ',it,itnew,ierr
	 if( ierr .ne. 0 ) itnew = it

	 if( .not. elabutil_check_time(it,itnew,itold) ) cycle

	 do i=1,nvar

	  ivar = ivars(i)
	  cv3(:,:) = cv3all(:,:,i)

	  nelab=nelab+1

	  if( .not. bquiet ) then
	    dline = ' '
	    if( bdate ) call dtsgf(it,dline)
	    write(6,*) 'time : ',it,'  ',dline,'   ivar : ',ivar
	  end if

	  if( bwrite ) then
	    do l=1,nlv
	      do k=1,nkn
	        cv2(k)=cv3(l,k)
	        if( l .gt. ilhkv(k) ) cv2(k) = rnull
	      end do
	      call mimar(cv2,nkn,cmin,cmax,rnull)
              call aver(cv2,nkn,cmed,rnull)
              call check1Dr(nkn,cv2,0.,-1.,"NaN check","cv2")
	      write(6,*) 'l,min,max,aver : ',l,cmin,cmax,cmed
	    end do
	  end if

	  if( btrans ) then
	    call nos_time_aver(mode,i,ifreq,istep,nkn,nlvdi
     +					,naccu,accum,cv3,boutput)
	  end if

	  if( baverbas ) then
	    call make_aver(nlvdi,nkn,ilhkv,cv3,vol3
     +                          ,cmin,cmax,cmed,vtot)
	    call write_aver(it,ivar,cmin,cmax,cmed,vtot)
	  end if

	  if( b2d ) then
	    call make_vert_aver(nlvdi,nkn,ilhkv,cv3,vol3,cv2)
	  end if

	  if( bsplit ) then
            call get_split_iu(ndim,iusplit,ivar,nin,ilhkv,hlv,hev,nb)
	  end if

	  if( boutput ) then
	    if( bverb ) write(6,*) 'writing to output: ',ivar
	    if( bsumvar ) ivar = 30
	    if( b2d ) then
              call nos_write_record(nb,it,ivar,1,ilhkv,cv2,ierr)
	    else
              call nos_write_record(nb,it,ivar,nlvdi,ilhkv,cv3,ierr)
	    end if
            if( ierr .ne. 0 ) goto 99
	  end if

	  if( bnodes ) then
	    do j=1,nnodes
	      node = nodes(j)
	      call write_node(j,node,cv3,it,ivar)
	    end do
	  end if

	 end do		!loop on ivar
	end do		!time do loop

c--------------------------------------------------------------
c end of loop on data
c--------------------------------------------------------------

c--------------------------------------------------------------
c final write of variables
c--------------------------------------------------------------

	if( btrans ) then
	  !write(6,*) 'istep,naccu: ',istep,naccu
	  do ip=1,istep
	    naccum = naccu(ip)
	    !write(6,*) 'naccum: ',naccum
	    if( naccum > 0 ) then
	      write(6,*) 'final aver: ',ip,naccum
	      call nos_time_aver(-mode,ip,ifreq,istep,nkn,nlvdi
     +					,naccu,accum,cv3,boutput)
	      if( bsumvar ) ivar = 30
              call nos_write_record(nb,it,ivar,nlvdi,ilhkv,cv3,ierr)
              if( ierr .ne. 0 ) goto 99
	    end if
	  end do
	end if

c--------------------------------------------------------------
c write final message
c--------------------------------------------------------------

	write(6,*)
	write(6,*) nread,' records read'
	write(6,*) nrec ,' unique time records read'
	write(6,*) nelab,' records elaborated'
	write(6,*)

	if( bsplit ) then
	  write(6,*) 'output written to following files: '
	  do ivar=1,ndim
	    if( iusplit(ivar) .gt. 0 ) then
              write(name,'(i4)') ivar
	      write(6,*) trim(adjustl(name))//'.nos'
	    end if
	  end do
	else if( boutput ) then
	  write(6,*) 'output written to file out.nos'
	end if

	call ap_get_names(basnam,simnam)
	write(6,*) 'names used: '
	write(6,*) 'basin: ',trim(basnam)
	write(6,*) 'simul: ',trim(simnam)

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	stop
   85	continue
	write(6,*) 'it,itvar,i,ivar,nvar: ',it,itvar,i,ivar,nvar
	stop 'error stop noselab: time mismatch'
   92	continue
	write(6,*) 'incompatible basin: '
	write(6,*) 'nkn,nknnos: ',nkn,nknnos
	write(6,*) 'nel,nelnos: ',nel,nelnos
	stop 'error stop noselab: parameter mismatch'
   99	continue
	write(6,*) 'error writing to file unit: ',nb
	stop 'error stop noselab: write error'
	end

c***************************************************************
c***************************************************************
c***************************************************************

	subroutine nos_time_aver(mode,nread,ifreq,istep,nkn,nlvddi
     +					,naccu,accum,cv3,bout)

c mode:  1:aver  2:sum  3:min  4:max
c
c mode negative: only transform, do not accumulate

	implicit none

	integer mode
	integer nread,ifreq,istep
	integer nkn,nlvddi
	integer naccu(istep)
	double precision accum(nlvddi,nkn,istep)
	real cv3(nlvddi,nkn)
	logical bout

	integer ip,naccum
	integer k,l

	bout = .false.
	ip = mod(nread,istep)
	if( ip .eq. 0 ) ip = istep

	!write(6,*) 'ip: ',ip,istep,nread,mode

	if( mode == 1 .or. mode == 2 ) then
	  naccu(ip) = naccu(ip) + 1
	  accum(:,:,ip) = accum(:,:,ip) + cv3(:,:)
	else if( mode == 3 ) then
	  do k=1,nkn
	    do l=1,nlvddi
	      accum(l,k,ip) = min(accum(l,k,ip),cv3(l,k))
	    end do
	  end do
	else if( mode == 3 ) then
	  do k=1,nkn
	    do l=1,nlvddi
	      accum(l,k,ip) = max(accum(l,k,ip),cv3(l,k))
	    end do
	  end do
	end if

	if( naccu(ip) == ifreq .or. mode < 0 ) then	!here ip == 1
	  naccum = naccu(ip)
	  if( abs(mode) > 1 ) naccum = 1
	  write(6,*) 'averaging: ',ip,naccum
	  bout = .true.
	  if( naccum > 0 ) then
	    cv3(:,:) = accum(:,:,ip)/naccum
	  end if
	  naccu(ip) = 0
	  accum(:,:,ip) = 0.
	end if

	end

c***************************************************************

        subroutine get_split_iu(ndim,iu,ivar,nin,ilhkv,hlv,hev,nb)

        implicit none

        integer ndim
        integer iu(ndim)
        integer ivar
        integer nin
        integer ilhkv(1)
        real hlv(1)
        real hev(1)
	integer nb		!unit to use for writing (return)

        integer nkn,nel,nlv,nvar
        integer ierr
        character*80 name

        if( ivar > ndim ) then
          write(6,*) 'ndim,ivar: ',ndim,ivar
          stop 'error stop: ndim'
        end if

        if( iu(ivar) .le. 0 ) then      !open file
          write(name,'(i4)') ivar
          call open_nos_file(name,'new',nb)
          call nos_init(nb,0)
          call nos_clone_params(nin,nb)
          call nos_get_params(nb,nkn,nel,nlv,nvar)
          call nos_set_params(nb,nkn,nel,nlv,1)
          call write_nos_header(nb,ilhkv,hlv,hev)
          iu(ivar) = nb
        end if

        nb = iu(ivar)

        end

c***************************************************************

