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

	program noselab

	use basin !COMMON_GGU_SUBST
	use clo

c elaborates nos file

	implicit none

	include 'param.h'
COMMON_GGU_DELETED	include 'basin.h'

	integer, parameter :: ndim = 1000
	integer iusplit(ndim)

	real, allocatable :: cv2(:)
	real, allocatable :: cv3(:,:)
	real, allocatable :: cv(:)
	real, allocatable :: vol3(:,:)

	integer, allocatable :: naccu(:)
	double precision, allocatable :: accum(:,:,:)

	integer, allocatable :: ilhkv(:)
	real, allocatable :: hlv(:)
	real, allocatable :: hl(:)
	real, allocatable :: hev(:)
	real, allocatable :: hkv(:)

	logical bwrite,bquiet,bask,bmem,bverb
	logical btrans,bout,bopen,btimew
	logical baver,bsum,bmin,bmax,bsumvar,b2d,baverbas
	logical bdate,bsplit,bneedbasin,boutput,bnode
	logical btmin,btmax
	integer mode,modeb
	integer date,time
	integer nread,nin
	integer nvers
	integer nknnos,nelnos,nlv,nvar,nlvdi
	integer nodesp,invar,lmaxsp
	integer ierr
	integer it,ivar
	integer l,k,lmax
	integer ifreq,istep,ip,nb,naccum
	integer datetime(2)
	character*80 title,name
	character*20 dline
	character*80 infile
	character*80 basnam,simnam
	character*80 stmin,stmax
	real rnull
	real cmin,cmax,cmed,vtot
	double precision atmin,atmax,atime,dtime

	integer iapini
	integer ifem_open_file

c--------------------------------------------------------------

	nread=0
	rnull=0.
	rnull=-1.
	bopen = .false.

c--------------------------------------------------------------
c set command line parameters
c--------------------------------------------------------------

        call clo_init('noselab','nos-file','2.1')

        call clo_add_info('returns info on or elaborates a nos file')

        call clo_add_sep('what to do (only one of these may be given)')

        call clo_add_option('out',.false.,'writes new nos file')
        call clo_add_option('averbas',.false.,'average over basin')
        call clo_add_option('aver',.false.,'average over records')
        call clo_add_option('sum',.false.,'sum over records')
        call clo_add_option('min',.false.,'minimum of records')
        call clo_add_option('max',.false.,'maximum of records')
        call clo_add_option('sumvar',.false.,'sum over variables')
        call clo_add_option('split',.false.,'split file for variables')
	call clo_add_option('2d',.false.,'average vertically to 2d field')

        call clo_add_sep('options in/output')

        !call clo_add_option('basin name',' ','name of basin to be used')
        call clo_add_option('mem',.false.,'if no file given use memory')
        call clo_add_option('ask',.false.,'ask for simulation')
        call clo_add_option('verb',.false.,'be more verbose')
        call clo_add_option('write',.false.,'write min/max of values')
        call clo_add_option('quiet',.false.,'do not be verbose')

        call clo_add_sep('additional options')

        call clo_add_option('node n',0,'extract vars of node number n')
	call clo_add_option('freq n',0.,'frequency for aver/sum/min/max')
        call clo_add_option('tmin time',' '
     +                          ,'only process starting from time')
        call clo_add_option('tmax time',' '
     +                          ,'only process up to time')

c--------------------------------------------------------------
c get command line parameters
c--------------------------------------------------------------

        call clo_parse_options

        call clo_get_option('out',bout)
        call clo_get_option('averbas',baverbas)
        call clo_get_option('aver',baver)
        call clo_get_option('sum',bsum)
        call clo_get_option('min',bmin)
        call clo_get_option('max',bmax)
        call clo_get_option('sumvar',bsumvar)
        call clo_get_option('split',bsplit)
        call clo_get_option('2d',b2d)

        call clo_get_option('node',nodesp)

        call clo_get_option('mem',bmem)
        call clo_get_option('ask',bask)
        call clo_get_option('verb',bverb)
        call clo_get_option('write',bwrite)
        call clo_get_option('quiet',bquiet)

        call clo_get_option('freq',ifreq)
        call clo_get_option('tmin',stmin)
        call clo_get_option('tmax',stmax)

	if( .not. bask .and. .not. bmem ) call clo_check_files(1)
	call clo_get_file(1,infile)
	call ap_set_names(' ',infile)

	if( .not. bquiet ) then
	  call shyfem_copyright('noselab - Elaborate NOS files')
	end if

	!--------------------------------------------------------------
	! open input files
	!--------------------------------------------------------------

	bnode = nodesp > 0

	boutput = bout .or. bsplit .or. b2d
	!btrans is added later

	bneedbasin = b2d .or. baverbas .or. bnode

	modeb = 2
	if( bneedbasin ) modeb = 3

	call ap_init(bask,modeb,0,0)

	call open_nos_type('.nos','old',nin)

	call peek_nos_header(nin,nknnos,nelnos,nlv,nvar)

	nlvdi = nlv

	allocate(cv2(nknnos))
	allocate(cv3(nlv,nknnos))
	allocate(vol3(nlv,nknnos))

	allocate(ilhkv(nknnos))
	allocate(hlv(nlv))
	allocate(hl(nlv))
	allocate(hev(nelnos))
	allocate(hkv(nknnos))

	call read_nos_header(nin,nknnos,nelnos,nlvdi,ilhkv,hlv,hev)
	call nos_get_params(nin,nknnos,nelnos,nlv,nvar)

	call init_sigma_info(nlv,hlv)

	call nos_get_date(nin,date,time)
	bdate = date .gt. 0
	if( bdate ) call dtsini(date,time)
	datetime(1) = date
	datetime(2) = time

	if( bneedbasin ) then
	  if( nkn /= nknnos .or. nel /= nelnos ) goto 92
	  call nos_make_hkv(nkn,nel,nen3v,hev,hkv)
	  call init_volume(nlvdi,nkn,nel,nlv,nen3v,ilhkv
     +                          ,hlv,hev,hl,vol3)
	else
	  nkn = nknnos
	  nel = nelnos
	end if

	call depth_stats(nkn,nlvdi,ilhkv)

	if( bnode ) then
	  call convert_internal_node(nodesp)
	  invar = 0
	  lmaxsp = ilhkv(nodesp)
	  allocate(cv(nvar))
	end if

	!--------------------------------------------------------------
	! time management
	!--------------------------------------------------------------

        atmin = 0.
        atmax = 0.
        btmin = stmin .ne. ' '
        btmax = stmax .ne. ' '
        if( btmin ) call fem_file_string2time(stmin,atmin)
        if( btmax ) call fem_file_string2time(stmax,atmax)

	if( bverb ) then
	  write(6,*) 'time limits: '
          write(6,*) stmin(1:len_trim(stmin)),btmin,atmin
          write(6,*) stmax(1:len_trim(stmax)),btmax,atmax
	end if
	if( .not. bquiet ) write(6,*) 

	!--------------------------------------------------------------
	! averaging
	!--------------------------------------------------------------

	mode = 0
	btrans = .false.
	if( baver ) mode = 1
	if( bsum )  mode = 2
	if( bmin )  mode = 3
	if( bmax )  mode = 4
	if( bsumvar )  mode = 2

	if( mode > 0 ) then	!prepare for averaging
	  btrans = .true.
	  if( bsumvar ) then		!sum over variables
	    istep = 1
	    ifreq = nvar
	  else if( nvar > 1 ) then
	    goto 91
	  else if( ifreq .ge. 0 ) then	!normal averaging
	    istep = 1
	  else				!accumulate every -ifreq record
	    istep = -ifreq
	    ifreq = 0
	  end if
	  allocate(naccu(istep))
	  allocate(accum(nlvdi,nknnos,istep))
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

	cv3 = 0.

	do

	  call nos_read_record(nin,it,ivar,nlvdi,ilhkv,cv3,ierr)

          if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
          if(ierr.ne.0) exit

	  nread=nread+1

	  dtime = it
	  call fem_file_convert_time(datetime,dtime,atime)
          btimew = .true.
          if( btmin ) btimew = btimew .and. atime >= atmin
          if( btmax ) btimew = btimew .and. atime <= atmax

	  if( .not. btimew ) cycle	!outside of time window

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
	    call nos_time_aver(mode,nread,ifreq,istep,nknnos,nlvdi
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

	  if( bnode ) then
	    invar = invar + 1
	    cv(invar) = cv3(1,nodesp)
	    if( invar .eq. nvar ) then
	      call write_node_2d(it,nvar,cv)
	      invar = 0
	    end if
	    call write_node(nodesp,nlvdi,cv3,it,ivar,lmaxsp)
	  end if

	end do	!do loop

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
	      call nos_time_aver(-mode,ip,ifreq,istep,nknnos,nlvdi
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
   92	continue
	write(6,*) 'incompatible basin: '
	write(6,*) 'nkn,nknnos: ',nkn,nknnos
	write(6,*) 'nel,nelnos: ',nel,nelnos
	stop 'error stop noselab: parameter mismatch'
   91	continue
	write(6,*) 'file contains different variables: ',nvar
	stop 'error stop noselab: averaging only with one variable'
   99	continue
	write(6,*) 'error writing to file unit: ',nb
	stop 'error stop noselab: write error'
	end

c***************************************************************
c***************************************************************
c***************************************************************

        subroutine aver(xx,n,xaver,rnull)

c computes min/max of vector
c
c xx            vector
c n             dimension of vector
c xmin,xmax     min/max value in vector
c rnull		invalid value

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

c***************************************************************

        subroutine mimar(xx,n,xmin,xmax,rnull)

c computes min/max of vector (2d)
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

        subroutine mimar_s(xx,nlvddi,n,xmin,xmax,rnull)

c computes min/max of vector (3d)
c
c xx            vector
c n             dimension of vector
c xmin,xmax     min/max value in vector
c rnull		invalid value

        implicit none

        integer n
        integer nlvddi
        real xx(nlvddi,n)
        real xmin,xmax,rnull

        integer k,l
        real x

        do k=1,n
          do l=1,nlvddi
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

	subroutine depth_stats(nkn,nlvddi,ilhkv)

c	computes statistics on levels

	implicit none

	include 'param.h'

	integer nkn
	integer nlvddi
	integer ilhkv(1)

	integer count(nlvddi)
	integer ccount(nlvddi)

	integer nlv,lmax,l,k,nc,ll

	nlv = 0
	do l=1,nlvddi
	  count(l) = 0
	  ccount(l) = 0
	end do

	do k=1,nkn
	  lmax = ilhkv(k)
	  if( lmax .gt. nlvddi ) stop 'error stop depth_stats: lmax'
	  count(lmax) = count(lmax) + 1
	  nlv = max(nlv,lmax)
	end do

	do l=nlv,1,-1
	  nc = count(l)
	  do ll=1,l
	    ccount(ll) = ccount(ll) + nc
	  end do
	end do

	nc = 0
	write(6,*) 'statistics for layers: ',nlv
	do l=1,nlv
	  if( count(l) > 0 ) then
	    write(6,*) l,count(l),ccount(l)
	    nc = nc + count(l)
	  end if
	end do
	write(6,*) 'total count: ',nc

	end

c***************************************************************

	subroutine nos_time_aver(mode,nread,ifreq,istep,nknnos,nlvddi
     +					,naccu,accum,cv3,bout)

c mode:  1:aver  2:sum  3:min  4:max
c
c mode negative: only transform, do not accumulate

	implicit none

	integer mode
	integer nread,ifreq,istep
	integer nknnos,nlvddi
	integer naccu(istep)
	double precision accum(nlvddi,nknnos,istep)
	real cv3(nlvddi,nknnos)
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
	  do k=1,nknnos
	    do l=1,nlvddi
	      accum(l,k,ip) = min(accum(l,k,ip),cv3(l,k))
	    end do
	  end do
	else if( mode == 3 ) then
	  do k=1,nknnos
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

	subroutine nos_make_hkv(nkn,nel,nen3v,hev,hkv)

c averages vertically

	implicit none

	integer nkn,nel
	integer nen3v(3,nel)
	real hev(nel)
	real hkv(nkn)

	integer k,ie,ii
	real h

	do ie=1,nel
	  h = hev(ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    hkv(k) = max(hkv(k),h)
	  end do
	end do

	end

c***************************************************************

	subroutine write_aver(it,ivar,cmin,cmax,cmed,vtot)

c writes basin average to file

	implicit none

	integer it,ivar
	real cmin,cmax,cmed,vtot

	real totmass

	totmass = cmed * vtot

        write(6,1234) it,ivar,cmin,cmed,cmax,totmass
        write(100+ivar,1235) it,cmin,cmed,cmax,totmass
        write(100,'(i10,e14.6)') it,vtot

 1234   format(2i10,3f12.4,e14.6)
 1235   format(i10,3f12.4,e14.6)
	end

c***************************************************************

	subroutine write_node(node,nlvddi,cv3,it,ivar,lmax)

	implicit none

	integer node
	integer nlvddi
	real cv3(nlvddi,*)
	integer it
	integer ivar
	integer lmax

	integer l

	write(88,*) it,node,ivar,lmax
	write(88,*) (cv3(l,node),l=1,lmax)

	end

c***************************************************************

	subroutine write_node_2d(it,nvar,cv)

	implicit none

	integer it
	integer nvar
	real cv(nvar)

	integer i

	write(89,*) it,(cv(i),i=1,nvar)

	end

c***************************************************************

	subroutine convert_internal_node(node)

	use basin !COMMON_GGU_SUBST

	implicit none

	include 'param.h'
COMMON_GGU_DELETED	include 'basin.h'

	integer node

	integer ni
	integer ipint

	if( node <= 0 ) return

	ni = ipint(node)

	if( ni <= 0 ) then
	  write(6,*) 'cannot find node: ',node
	  stop 'error stop convert_internal_node: no such node'
	end if

	node = ni

	end

c***************************************************************
