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
c 10.09.2015    ggu     std and rms for averaging implemented
c 11.09.2015    ggu     write in gis format
c 23.09.2015    ggu     handle more than one file (look for itstart)
c 19.02.2016    ggu     bug fixes for bsumvar and mode==2 (sum)
c 22.02.2016    ggu     handle catmode
c 08.09.2016    ggu     some minor changes, map_influence, custom dates
c
c**************************************************************

	subroutine noselab

	use clo
	use elabutil
	use elabtime
	use custom_dates

        use basin
        use mod_depth
        use evgeom
        use levels

c elaborates nos file

	implicit none

	integer, parameter :: ndim = 1000
	integer iusplit(ndim)

	real, allocatable :: cv2(:)
	real, allocatable :: cv3(:,:)
	real, allocatable :: cv3all(:,:,:)
	real, allocatable :: vol3(:,:)

	integer, allocatable :: ivars(:)
	!integer, allocatable :: nodes(:)

	integer, allocatable :: naccu(:)
	double precision, allocatable :: accum(:,:,:)
	double precision, allocatable :: std(:,:,:,:)

	real, allocatable :: hl(:)

	logical bforce
	integer nwrite,nread,nelab,nrec,nin,nold
	integer nvers
	integer nknnos,nelnos,nvar
	integer ierr
	integer it,ivar,itvar,itnew,itold,iaux,itstart,iv,itfirst
	integer i,j,l,k,lmax,node
	integer ip,nb,naccum
	integer ifile
	integer date,time
	character*80 title,name,file
	character*20 dline
	character*80 basnam,simnam
	real rnull
	real cmin,cmax,cmed,vtot
	double precision dtime

	!logical, parameter :: bmap = .false.
	real, parameter :: pthresh = 30.
	real, parameter :: cthresh = 20.

	integer iapini
	integer ifem_open_file
	logical concat_cycle

c--------------------------------------------------------------

	nread=0
	nwrite=0
	nelab=0
	nrec=0
	rnull=0.
	rnull=-1.
	bopen = .false.
	bforce = .false.

c--------------------------------------------------------------
c set command line parameters
c--------------------------------------------------------------

	call elabutil_init('NOS')

	!--------------------------------------------------------------
	! open input files
	!--------------------------------------------------------------

	call ap_init(bask,modeb,0,0)

	!call open_nos_type('.nos','old',nin)
	ifile = 1
	call clo_get_file(ifile,file)
	call open_shy_file(file,'old',nin)
	write(6,*) '================================'
	write(6,*) 'reading file: ',trim(file)
	write(6,*) '================================'
	call clo_get_file(ifile+1,file)
	itstart = -1
	if( file /= ' ' ) call nos_get_it_start(file,itstart)

	call nos_is_nos_file(nin,nvers)
	if( nvers .le. 0 ) then
	  write(6,*) 'nvers: ',nvers
	  stop 'error stop noselab: not a valid nos file'
	end if

	call peek_nos_header(nin,nknnos,nelnos,nlv,nvar)

	if( bneedbasin ) then
	  if( nkn /= nknnos .or. nel /= nelnos ) goto 92
	else
	  nkn = nknnos
	  nel = nelnos
	end if

        call mod_depth_init(nkn,nel)
        call levels_init(nkn,nel,nlv)
	if( bneedbasin ) then
          call ev_init(nel)
          call set_ev
	end if

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
	  call outfile_make_hkv(nkn,nel,nen3v,hm3v,hev,hkv)
	  call ilhk2e(nkn,nel,nen3v,ilhkv,ilhv)
	  call adjust_layer_index(nel,nlv,hev,hlv,ilhv)
	  call init_volume(nlvdi,nkn,nel,nlv,nen3v,ilhkv
     +                          ,hlv,hev,hl,vol3)
	  !vol3=1.
	end if

	!if( bverb ) write(6,*) 'hlv: ',nlv,hlv
	if( bverb ) call depth_stats(nkn,nlvdi,ilhkv)
	if( binfo ) return

	call handle_nodes

	!--------------------------------------------------------------
	! time management
	!--------------------------------------------------------------

	call nos_get_date(nin,date,time)
	call elabtime_date_and_time(date,time)
	bdate = date > 0

	!--------------------------------------------------------------
	! averaging
	!--------------------------------------------------------------

	call elabutil_set_averaging(nvar)
	call custom_dates_init(itstart,datefile)
        !write(6,*) 'ggu: ',ifreq,istep

        if( btrans .and. nvar > 1 ) then
	  if( .not. bsumvar ) then
            stop 'error stop noselab: only one variable with averaging'
	  end if
        end if

	if( btrans ) then
	  allocate(naccu(istep))
	  allocate(accum(nlvdi,nkn,istep))
	  allocate(std(nlvdi,nkn,16,istep))	!also used for directions
	else
	  allocate(naccu(1))
	  allocate(accum(1,1,1))
	  allocate(std(1,1,1,1))	!also used for directions
	end if
	naccum = 0
	naccu = 0
	accum = 0.
	std = 0.

	!write(6,*) 'mode: ',mode,ifreq,istep

	!--------------------------------------------------------------
	! open output file
	!--------------------------------------------------------------

	iusplit = 0

	boutput = boutput .or. btrans .or. bsplit
	bopen = boutput .and. .not. bsplit
	bopen = boutput .or. bsumvar

	if( bopen ) then
          call open_nos_file('out','new',nb)
          call nos_init(nb,0)
          call nos_clone_params(nin,nb)
	  if( b2d ) then
	    call nos_set_params(nb,0,0,1,0)
	  end if
	  if( bsumvar ) then
	    call nos_set_params(nb,0,0,0,1)
	  end if
          call write_nos_header(nb,ilhkv,hlv,hev)
	end if

	if( outformat == 'gis' ) call gis_write_connect

c--------------------------------------------------------------
c loop on data
c--------------------------------------------------------------

	it = 0
	!if( .not. bquiet ) write(6,*)

	cv3 = 0.

	do

	 itold = it

	 do iv=1,nvar
	  call nos_read_record(nin,it,ivar,nlvdi,ilhkv,cv3,ierr)
          if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
          if(ierr.ne.0) exit
	  if( iv == 1 ) itvar = it
	  if( itvar /= it ) goto 85
	  ivars(iv) = ivar
	  cv3all(:,:,iv) = cv3(:,:) * fact
	  nread=nread+1
	 end do

         if(ierr.ne.0) then	!EOF - see if we have to read another file
	   it = itold
	   if( itstart == -1 ) exit
	   nold = nin
	   ifile = ifile + 1
	   call clo_get_file(ifile,file)
	   call open_shy_file(file,'old',nin)
	   write(6,*) '================================'
	   write(6,*) 'reading file: ',trim(file)
	   write(6,*) '================================'
	   call read_nos_header(nin,nkn,nel,nlvdi,ilhkv,hlv,hev)
	   call nos_check_compatibility(nin,nold)
	   call nos_peek_record(nin,itnew,iaux,ierr)
	   call nos_get_date(nin,date,time)
	   call elabtime_date_and_time(date,time)
	   bdate = date > 0
	   it = itold		!reset time of last successfully read record
	   call nos_close(nold)
	   close(nold)
	   call clo_get_file(ifile+1,file)
	   itstart = -1
	   if( file /= ' ' ) call nos_get_it_start(file,itstart)
	   cycle
	 end if

	 nrec = nrec + 1

	 if( nrec == 1 ) itold = it
	 call nos_peek_record(nin,itnew,iaux,ierr)
	 !write(6,*) 'peek: ',it,itnew,ierr
	 if( ierr .ne. 0 ) itnew = it

	 if( concat_cycle(it,itold,itstart,nrec) ) cycle
	 !if( elabutil_over_time_a(atime,atnew,atold) ) exit
	 if( .not. elabtime_check_time(it,itnew,itold) ) cycle

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
	      write(6,1000) 'l,min,max,aver : ',l,cmin,cmax,cmed
 1000	      format(a,i5,3g16.6)
	    end do
	  end if

	  if( btrans ) then
	    call custom_dates_over(it,bforce)
	    call nos_time_aver(bforce,avermode,nread,ifreq,istep,nkn,nlvdi
     +				,naccu,accum,std,threshold,cv3,boutput)
	  end if

	  if( baverbas ) then
	    call make_basin_aver(nlvdi,nkn,ilhkv,cv3,vol3
     +                          ,cmin,cmax,cmed,vtot)
	    call write_aver(it,ivar,cmin,cmax,cmed,vtot)
	  end if

	  if( b2d ) then
	    call make_vert_aver(nlvdi,nkn,ilhkv,cv3,vol3,cv2)
	  end if

	  if( bsplit ) then
	    !write(6,*) 'splitting ',ivar,iusplit(ivar),boutput
            call get_split_iu(ndim,iusplit,ivar,nin,ilhkv,hlv,hev,nb)
	  end if

	  if( bmap ) then	!creates influence map
	    boutput = .false.
	    if( i == nvar ) then
	      boutput = .true.
	      ivar = 75
	      call comp_map0(nlvdi,nkn,nvar,pthresh,cthresh,cv3all,cv3)
	    end if
	  end if

	  if( bsumvar ) then
	    boutput = .false.
	    if( i == nvar ) then
	      boutput = .true.
	      ivar = 10
	      cv3 = sum(cv3all,3)
	    end if
	  end if

	  if( boutput ) then
	    nwrite = nwrite + 1
	    if( bverb ) write(6,*) 'writing to output: ',ivar
	    if( bthreshold ) ivar = 199
	    if( b2d ) then
              call noselab_write_record(nb,it,ivar,1,ilhkv,cv2,ierr)
	    else
              call noselab_write_record(nb,it,ivar,nlvdi,ilhkv,cv3,ierr)
	    end if
            if( ierr .ne. 0 ) goto 99
	  end if

	  if( bnodes ) then
	    dtime = it
	    call write_nodes(dtime,ivar,cv3)
	    if( .not. b2d ) then	!still to average
	      call make_vert_aver(nlvdi,nkn,ilhkv,cv3,vol3,cv2)
	    end if
	    call write_2d_all_nodes(nnodes,nodes,cv2,it,ivar)
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
	      nwrite = nwrite + 1
	      !write(6,*) 'final aver: ',ip,naccum
	      call nos_time_aver(bforce,-avermode,ip,ifreq,istep,nkn,nlvdi
     +				,naccu,accum,std,threshold,cv3,boutput)
	      if( bsumvar ) ivar = 10
	      if( bthreshold ) ivar = 199
              call nos_write_record(nb,it,ivar,nlvdi,ilhkv,cv3,ierr)
              if( ierr .ne. 0 ) goto 99
	    end if
	  end do
	end if

c--------------------------------------------------------------
c write final message
c--------------------------------------------------------------

	write(6,*)
	write(6,*) nread, ' records read'
	write(6,*) nrec , ' unique time records read'
	write(6,*) nelab, ' records elaborated'
	write(6,*) ifile, ' files read'
	write(6,*) nwrite,' records written'
	write(6,*)

	if( bsplit ) then
	  write(6,*) 'output written to following files: '
	  do ivar=1,ndim
	    nb = iusplit(ivar)
	    if( nb .gt. 0 ) then
              write(name,'(i4)') ivar
	      write(6,*) trim(adjustl(name))//'.nos'
	      close(nb)
	    end if
	  end do
	else if( boutput ) then
	  write(6,*) 'output written to file out.nos'
	  close(nb)
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
	write(6,*) 'it,itvar,iv,ivar,nvar: ',it,itvar,iv,ivar,nvar
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

	subroutine nos_time_aver(bforce
     +				,mode,nread,ifreq,istep,nkn,nlvddi
     +				,naccu,accum,std,threshold,cv3,bout)

c mode:  1:aver  2:sum  3:min  4:max  5:std  6:rms  7:thres  8:averdir
c
c mode negative: only transform, do not accumulate

	implicit none

	logical bforce
	integer mode
	integer nread,ifreq,istep
	integer nkn,nlvddi
	integer naccu(istep)
	double precision accum(nlvddi,nkn,istep)
	double precision std(nlvddi,nkn,16,istep)
	double precision threshold
	real cv3(nlvddi,nkn)
	logical bout

	integer ip,naccum,mmode
	integer k,l,id,idmax
	double precision dmax

	if( mode .eq. 0 ) return

	bout = .false.
	ip = mod(nread,istep)
	if( ip .eq. 0 ) ip = istep

	!write(6,*) 'ip1: ',mode,ifreq,istep,nread,ip,naccu(ip)

	if( mode == 1 .or. mode == 2 ) then
	  accum(:,:,ip) = accum(:,:,ip) + cv3(:,:)
	else if( mode == 3 ) then
	  accum(:,:,ip) = min(accum(:,:,ip),cv3(:,:))
	  !do k=1,nkn
	  !  do l=1,nlvddi
	  !    accum(l,k,ip) = min(accum(l,k,ip),cv3(l,k))
	  !  end do
	  !end do
	else if( mode == 4 ) then
	  accum(:,:,ip) = max(accum(:,:,ip),cv3(:,:))
	  !do k=1,nkn
	  !  do l=1,nlvddi
	  !    accum(l,k,ip) = max(accum(l,k,ip),cv3(l,k))
	  !  end do
	  !end do
	else if( mode == 5 ) then
	  accum(:,:,ip) = accum(:,:,ip) + cv3(:,:)
	  std(:,:,1,ip) = std(:,:,1,ip) + cv3(:,:)**2
	else if( mode == 6 ) then
	  accum(:,:,ip) = accum(:,:,ip) + cv3(:,:)**2
	else if( mode == 7 ) then
	  where( cv3(:,:) >= threshold )
	    accum(:,:,ip) = accum(:,:,ip) + 1.
	  end where
	else if( mode == 8 ) then
	  do k=1,nkn
	    do l=1,nlvddi
	      id = nint( cv3(l,k)/22.5 )
	      if( id == 0 ) id = 16
	      if( id < 0 .or. id > 16 ) stop 'error stop: direction'
	      std(l,k,id,ip) = std(l,k,id,ip) + 1.
	    end do
	  end do
	end if

	if( mode > 0 ) naccu(ip) = naccu(ip) + 1
	!write(6,*) 'ip2: ',mode,ifreq,istep,nread,ip,naccu(ip)

	if( naccu(ip) == ifreq .or. mode < 0 .or. bforce ) then
	  naccum = max(1,naccu(ip))
	  mmode = abs(mode)
	  if( mmode == 2 ) naccum = 1			!sum
	  if( mmode == 3 ) naccum = 1			!min
	  if( mmode == 4 ) naccum = 1			!max
	  if( mmode == 7 ) naccum = 1			!threshold
	  if( naccum > 0 ) cv3(:,:) = accum(:,:,ip)/naccum
	  if( mmode == 5 ) then
	    cv3(:,:) = sqrt( std(:,:,1,ip)/naccum - cv3(:,:)**2 )
	  else if( mmode == 6 ) then
	    cv3(:,:) = sqrt( cv3(:,:) )
	  else if( mmode == 8 ) then
	    do k=1,nkn
	      do l=1,nlvddi
		dmax = 0.
		idmax = 0
	        do id=1,16
		  if( std(l,k,id,ip) > dmax ) then
		    idmax = id
		    dmax = std(l,k,id,ip)
		  end if
		end do
		if( idmax == 16 ) idmax = 0
		cv3(l,k) = idmax * 22.5
	      end do
	    end do
	  end if
	  write(6,*) 'averaging: ',ip,naccum,naccu(ip)
	  bout = .true.
	  naccu(ip) = 0
	  accum(:,:,ip) = 0.
	  std(:,:,:,ip) = 0.
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

	subroutine noselab_write_record(nb,it,ivar,nlvdi,ilhkv,cv,ierr)

	use basin
        use elabutil

	implicit none

	integer nb,it,ivar,nlvdi
	integer ilhkv(nlvdi)
	real cv(nlvdi,*)
	integer ierr

	double precision dtime

	ierr = 0

	if( outformat == 'nos' .or. outformat == 'native') then
          call nos_write_record(nb,it,ivar,nlvdi,ilhkv,cv,ierr)
	else if( outformat == 'gis' ) then
	  dtime = it
          call gis_write_record(dtime,ivar,nkn,nlvdi
     +				,ilhkv,cv,xgv,ygv)
	else
	  write(6,*) 'output format unknown: ',outformat
	  stop 'error stop noselab_write_record: output format'
	end if

        end

c***************************************************************

        subroutine comp_map0(nlvdi,nkn,nvar,pt,ct,cvv,valri)

c compute dominant discharge and put index in valri

        implicit none

        integer nlvdi,nkn,nvar
        real pt,ct
        real cvv(nlvdi,nkn,nvar)
        real valri(nlvdi,nkn)

        integer iv,k,ismax,l
        real conz, pconz
        real sum,rmax
        real cthresh,pthresh

        pthresh = pt     !threshold on percentage
        cthresh = ct     !threshold on concentration - 0 for everywhere

	!pthresh = 30.
	!cthresh = 20.

        do l=1,nlvdi
          do k=1,nkn
                sum = 0.
                rmax = 0.
                ismax = 0
                do iv=1,nvar
                   conz = cvv(l,k,iv)
                   sum = sum + conz
                   if( conz .gt. rmax ) then
                        rmax = conz
                        ismax = iv
                   end if
                end do

                conz = 0.
                if( ismax .gt. 0 ) conz = cvv(l,k,ismax)
                pconz = 0.
                if( sum .gt. 0. ) pconz = (conz/sum)*100

                valri(l,k) = 0.
                if( conz .gt. cthresh ) then
                  if( pconz .gt. pthresh ) then
                    valri(l,k) = ismax
                  end if
                end if
          end do
        end do

        end

c***************************************************************


