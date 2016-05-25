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
c
c**************************************************************

	program nos2shy

	use clo
	use elabutil
	use elabtime
	use shyfile

        use basin
        use mod_depth
        use evgeom
        use levels

c elaborates nos file

	implicit none

	include 'param.h'

	integer, parameter :: ndim = 1000
	integer iusplit(ndim)

	real, allocatable :: cv3(:,:)
	real, allocatable :: cv3all(:,:,:)
	real, allocatable :: vol3(:,:)

	integer, allocatable :: ivars(:)
	real, allocatable :: hl(:)

	integer nwrite,nread,nelab,nrec,nin,nold
	integer nvers
	integer nknnos,nelnos,nvar
	integer ierr
	integer it,ivar,itvar,itnew,iaux,itstart
	integer i,j,l,k,lmax,node
	integer ip,nb,naccum
	integer ifile
	integer id,ftype
	integer date,time
	character*80 title,name,file
	character*80 sfile
	character*20 dline
	character*80 basnam,simnam
	real rnull
	real cmin,cmax,cmed,vtot
	double precision dtime

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

c--------------------------------------------------------------
c set command line parameters
c--------------------------------------------------------------

	call elabutil_init('NOS')

	!--------------------------------------------------------------
	! open input files
	!--------------------------------------------------------------

	boutput = .true.
	bneedbasin = .true.
	!call ap_init(bask,modeb,0,0)
	call ap_init_basin

	ifile = 1
	call clo_get_file(ifile,file)
	call open_shy_file(file,'old',nin)
	write(6,*) '================================'
	write(6,*) 'reading file: ',trim(file)
	write(6,*) '================================'

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
          call ilhk2e(nkn,nel,nen3v,ilhkv,ilhv)
          call adjust_layer_index(nel,nlv,hev,hlv,ilhv)
	  call init_volume(nlvdi,nkn,nel,nlv,nen3v,ilhkv
     +                          ,hlv,hev,hl,vol3)
	end if

	if( bverb ) call depth_stats(nkn,nlvdi,ilhkv)

	!--------------------------------------------------------------
	! time management
	!--------------------------------------------------------------

	call nos_get_date(nin,date,time)
	call elabtime_date_and_time(date,time)

	!--------------------------------------------------------------
	! open output file
	!--------------------------------------------------------------

	bopen = boutput

	if( bopen ) then
	  sfile = 'out.scal.shy'
	  ftype = 2
	  call shy_open_output_file(sfile,1,nlv,nvar,ftype,id)
	  call nos_transfer_simul_params(nin,id)
          call shy_make_header(id)
	end if

c--------------------------------------------------------------
c loop on data
c--------------------------------------------------------------

	it = 0
	cv3 = 0.

	do

	 do i=1,nvar
	  call nos_read_record(nin,it,ivar,nlvdi,ilhkv,cv3,ierr)
          if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
          if(ierr.ne.0) exit
	  if( i == 1 ) itvar = it
	  if( itvar /= it ) goto 85
	  ivars(i) = ivar
	  nread=nread+1
	  cv3all(:,:,i) = cv3(:,:)
	 end do

         if(ierr.ne.0) exit

	 dtime = it
	 nrec = nrec + 1
	 nelab = nelab + 1

	 do i=1,nvar

	  ivar = ivars(i)
	  cv3(:,:) = cv3all(:,:,i)

	  if( .not. bquiet ) then
	    dline = ' '
	    if( bdate ) call dtsgf(it,dline)
	    write(6,*) 'time : ',it,'  ',dline,'   ivar : ',ivar
	  end if

	  if( boutput ) then
	    nwrite = nwrite + 1
	    if( bverb ) write(6,*) 'writing to output: ',ivar
	    call shy_write_scalar_record(id,dtime,ivar,nlv,cv3)
	  end if

	 end do		!loop on ivar
	end do		!time do loop

c--------------------------------------------------------------
c end of loop on data
c--------------------------------------------------------------

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

	if( boutput ) then
	  write(6,*) 'output written to file ',trim(sfile)
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

        subroutine nos_transfer_simul_params(nin,id)

	use shyfile

	implicit none

	integer nin,id

	integer date,time
	character*80 title
	character*80 femver

	call nos_get_date(nin,date,time)
        call nos_get_title(nin,title)
        call nos_get_femver(nin,femver)

        call shy_set_date(id,date,time)
        call shy_set_title(id,title)
        call shy_set_femver(id,femver)

	end

c***************************************************************

