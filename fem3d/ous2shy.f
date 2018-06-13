c
c $Id: ouselab.f,v 1.8 2008-11-20 10:51:34 georg Exp $
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
c 11.09.2015    ggu     split feature added
c 17.03.2016    ggu     outformat git added
c
c**************************************************************

	program ous2shy

	use clo
	!use elabutil
	use elabtime

	use basin
        use mod_depth
        use mod_hydro
        use mod_hydro_baro
        use evgeom
        use levels
        use shyutil
        use shympi
        use iso8601

c elaborates ous file

	implicit none

	real, allocatable :: zv(:)
	real, allocatable :: ze(:,:)
	real, allocatable :: uv(:,:)
	real, allocatable :: vv(:,:)

	real, allocatable :: vars(:,:,:)
	real, allocatable :: vars2d(:,:,:)
	integer, allocatable :: idims(:,:)

	real, allocatable :: hl(:)

	logical bopen,bneedbasin,boutput,bquiet,bverb,bsilent
	logical bdate,bdstring
	integer nread,nelab,nrec,nin,nwrite
	integer nndim,nvar,iv
	integer nvers
	integer nknous,nelous
	integer invar
	integer ierr
	integer it,ivar,itvar,itnew,itold,iaux
	integer i,l,k,lmax
	integer iano,ks
	integer ip,nb,naccum
	integer ftype,id
	integer date,time
	character*80 title,name
	character*80 sfile
	character*20 dline
	character*20 vers
	character*80 version
	character*80 dstring
	character*80 basnam,simnam
	character*80 basfile,simfile
	real rnull
	real zmin,zmax
	real umin,umax,vmin,vmax
	real volume,area
	double precision dtime,dtfirst,dtlast

	integer iapini
	integer ifem_open_file

c--------------------------------------------------------------

	!ifile = 1
        nwrite = 0
	nread=0
	nelab=0
	nrec=0
	rnull=0.
	rnull=-999.
	bopen = .false.
	dtfirst = 0.
	dtlast = 0.

	ks = -1			!write special node
	iano = -1		!no computation for this area code

c--------------------------------------------------------------
c set command line parameters
c--------------------------------------------------------------

        call get_shyfem_version(vers)
        version = '2.0' // ' (SHYFEM version ' // trim(vers) // ')'

        call clo_init('ous2shy','ous-file bas-file',version)
        call clo_add_info('transforms ous file into shy file')

        call clo_add_sep('general options')
        call clo_add_option('verbose',.false.,'be more verbose')
        call clo_add_option('quiet',.false.,'do not write time records')
        call clo_add_option('silent',.false.,'do not write anything')

        call clo_add_sep('date options')
        call clo_add_option('date date',' '
     +		,'give reference date if not contained in file')

        call clo_parse_options
 
        call clo_get_option('verbose',bverb)
        call clo_get_option('quiet',bquiet)
        call clo_get_option('silent',bsilent)
        call clo_get_option('date',dstring)

	if( .not. bsilent ) then
	 call shyfem_copyright('ous2shy - Transforms OUS to SHY files')
	end if

	!--------------------------------------------------------------
	! open input files
	!--------------------------------------------------------------

	boutput = .true.
	bneedbasin = .true.
	if( bsilent ) bquiet = .true.

	call compat_get_files(simfile,basfile)
	!call ap_init(bask,modeb,0,0)
	!call ap_init_basin

	call basin_read(basfile)

	call open_shy_file(simfile,'old',nin)
	if( .not. bquiet ) then
	write(6,*) '================================'
	write(6,*) 'reading file: ',trim(simfile)
	write(6,*) '================================'
	end if

        call ous_is_ous_file(nin,nvers)
        if( nvers .le. 0 ) then
          write(6,*) 'nvers: ',nvers
          stop 'error stop ouselab: not a valid ous file'
        end if

	call peek_ous_header(nin,nknous,nelous,nlv)

	if( bneedbasin ) then
	  if( nkn /= nknous .or. nel /= nelous ) goto 92
	else
	  nkn = nknous
	  nel = nelous
	end if

	call mod_hydro_baro_init(nel)
	call mod_depth_init(nkn,nel)
	call levels_init(nkn,nel,nlv)
	call mod_hydro_init(nkn,nel,nlv)

	allocate(hl(nlv))
	allocate(zv(nkn))
	allocate(ze(3,nel))
	allocate(uv(nlvdi,nel))
	allocate(vv(nlvdi,nel))

	nndim = max(3*nel,nkn)
	nvar = 4
	allocate(vars(nlvdi,nndim,nvar))

	nlvdi = nlv
        call read_ous_header(nin,nkn,nel,nlvdi,ilhv,hlv,hev)
        call ous_get_params(nin,nkn,nel,nlv)

	if( .not. bquiet ) then
	  call write_ous_info(nin)
          write(6,*) 'levels: '
          write(6,'(5g14.6)') (hlv(l),l=1,nlv)
	end if

	call init_sigma_info(nlv,hlv)

	if( bneedbasin ) then
	  call outfile_make_hkv(nkn,nel,nen3v,hm3v,hev,hkv)
	  call ilhe2k(nkn,nel,nen3v,ilhv,ilhkv)
	end if

        call shympi_init(.false.)               !call after basin has been read
	call shympi_set_hlv(nlv,hlv)

	if( bverb ) call depth_stats(nkn,nlvdi,ilhkv)

	!--------------------------------------------------------------
	! time management
	!--------------------------------------------------------------

	call ous_get_date(nin,date,time)
	if( date < 10000 ) date = date * 10000 + 101
	bdate = ( date > 0 ) 
	bdstring = ( dstring /= ' ' )
	if( .not. bdate .and. .not. bdstring ) then
	  write(6,*) 'no date in file... please specify with -date'
	  stop 'error stop ous2shy: no date'
	else if( bdate .and. bdstring ) then
	  write(6,*) 'date in file... cannot specify -date'
	  stop 'error stop ous2shy: no -date possible'
	else if( bdstring ) then
	  call string2date(dstring,date,time,ierr)
	  if( ierr /= 0 ) then
	    write(6,*) 'cannot parse date string: ',trim(dstring)
	    stop 'error stop ous2shy: error in date'
	  end if
	end if
	call ous_set_date(nin,date,time)
	call elabtime_date_and_time(date,time)
	call date2string(date,time,dline)
	if( .not. bquiet ) write(6,*) 'reference date used: ',dline
	bdate = .true.

	!--------------------------------------------------------------
	! open output file
	!--------------------------------------------------------------

	bopen = boutput

	if( bopen ) then
          sfile = 'out.hydro.shy'
          ftype = 1
          call shy_open_output_file(sfile,3,nlv,nvar,ftype,id)
          call ous_transfer_simul_params(nin,id)
          call shy_make_header(id)
	end if

c--------------------------------------------------------------
c loop on data
c--------------------------------------------------------------

	area = 0.
	volume = 0.
	it = 0
	if( .not. bquiet ) write(6,*)

	do

	  itold = it

	  call new_read_record2(nin,it,nlvdi,nndim,nvar,vars,ierr)

          if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
          if(ierr.ne.0) exit

	  dtime = it
	  if( nrec == 0 ) dtfirst = dtime
	  nread=nread+1
	  nrec = nrec + 1

	  nelab=nelab+1

	  if( bverb ) then
	    call dtsgf(it,dline)
	    write(6,*) 'time : ',it,'  ',dline
	  end if

	  if( boutput ) then
            nwrite = nwrite + 1
	    if( bverb ) write(6,*) 'writing to output '
	    call transfer_uvz2(nlvdi,nndim,nvar,vars,zv,ze,uv,vv)
	    call shy_write_hydro_records(id,dtime,nlvdi,zv,ze,uv,vv)
	  end if

	  dtlast = dtime

	end do		!time do loop

c--------------------------------------------------------------
c end of loop on data
c--------------------------------------------------------------

c--------------------------------------------------------------
c write final message
c--------------------------------------------------------------

	if( .not. bsilent ) then
	write(6,*)
	write(6,*) nread,' total records read'
	!write(6,*) nrec ,' unique time records read'
	write(6,*) nelab,' records elaborated'
	write(6,*) nwrite,' records written'
	write(6,*)
	it = dtfirst
	call dtsgf(it,dline)
	write(6,*) 'first time record: ',dline
	it = dtlast
	call dtsgf(it,dline)
	write(6,*) 'last time record:  ',dline
	write(6,*)
	if( boutput ) then
	  write(6,*) 'output written to file ',trim(sfile)
	end if
	end if

	!call ap_get_names(basnam,simnam)
	!write(6,*) 'names used: '
	!write(6,*) 'basin: ',trim(basnam)
	!write(6,*) 'simul: ',trim(simnam)

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	stop
   92	continue
	write(6,*) 'incompatible basin: '
	write(6,*) 'nkn,nknous: ',nkn,nknous
	write(6,*) 'nel,nelous: ',nel,nelous
	stop 'error stop ouselab: parameter mismatch'
   99	continue
	write(6,*) 'error writing to file unit: ',nb
	stop 'error stop ouselab: write error'
	end

c***************************************************************
c***************************************************************
c***************************************************************


c***************************************************************
c***************************************************************
c***************************************************************

	subroutine new_read_record2(nb,it,nlvddi,nndim,nvar,vars,ierr)

	use basin
	use levels

	implicit none

	integer nb
	integer it
	integer nlvddi,nndim,nvar
	real vars(nlvddi,nndim,nvar)
	integer ierr

	real znv(nkn)
	real zenv(3*nel)
	real utlnv(nlvddi,nel)
	real vtlnv(nlvddi,nel)

        call ous_read_record(nb,it,nlvddi,ilhv,znv,zenv
     +                          ,utlnv,vtlnv,ierr)

	vars(1,1:nkn,1)   = znv(1:nkn)
	vars(1,1:3*nel,2) = zenv(1:3*nel)
	vars(:,1:nel,3)   = utlnv(:,1:nel)
	vars(:,1:nel,4)   = vtlnv(:,1:nel)

	end

c***************************************************************

	subroutine transfer_uvz2(nlvddi,nndim,nvar,vars,znv,zenv,uv,vv)

	use basin

	implicit none

	integer nlvddi,nndim,nvar
	real vars(nlvddi,nndim,nvar)
	real znv(nkn)
	real zenv(3*nel)
	real uv(nlvddi,nel)
	real vv(nlvddi,nel)

	znv(1:nkn)     = vars(1,1:nkn,1)
	zenv(1:3*nel)  = vars(1,1:3*nel,2)
	uv(:,1:nel)    = vars(:,1:nel,3)
	vv(:,1:nel)    = vars(:,1:nel,4)

	end

c***************************************************************

        subroutine ous_transfer_simul_params(nin,id)

        use shyfile

        implicit none

        integer nin,id

        integer date,time
        character*80 title
        character*80 femver

        call ous_get_date(nin,date,time)
        call ous_get_title(nin,title)
        call ous_get_femver(nin,femver)

        call shy_set_date(id,date,time)
        call shy_set_title(id,title)
        call shy_set_femver(id,femver)

        end

c***************************************************************
c***************************************************************
c***************************************************************

	subroutine compat_get_files(simfile,basfile)

! gets exactly one simfile and one basfile from command line

	use clo
	use basin

	implicit none

	character*(*) simfile
	character*(*) basfile

	integer ifile,nfile,ios
	integer nin,nvers
	character*80 file

	simfile = ' '
	basfile = ' '

	nfile = clo_number_of_files()

	do ifile=1,nfile

	  call clo_get_file(ifile,file)

	  call open_shy_file(file,'old',nin)
	  call ous_is_ous_file(nin,nvers)
	  close(nin,iostat=ios)
	  if( nvers > 0 ) then
	    if( simfile /= ' ' ) goto 99
	    simfile = file
	    cycle
	  end if

	  if( basin_is_basin(file) ) then
	    if( basfile /= ' ' ) goto 98
	    basfile = file
	    cycle
	  end if

	  goto 97

	end do

	if( simfile == ' ' ) then
	  write(6,*) 'no OUS file given on command line'
	  call clo_usage
	  stop 'error stop compat_get_files: no file'
	end if

	if( basfile == ' ' ) then
	  write(6,*) 'no BAS file given on command line'
	  call clo_usage
	  stop 'error stop compat_get_files: no file'
	end if

	return
   97	continue
	write(6,*) 'not recognized file type: ',trim(file)
	stop 'error stop compat_get_files: unknown file'
   98	continue
	write(6,*) 'more than one BAS file on command line'
	write(6,*) 'can handle only one file'
	stop 'error stop compat_get_files: BAS files'
   99	continue
	write(6,*) 'more than one OUS file on command line'
	write(6,*) 'can handle only one file'
	stop 'error stop compat_get_files: OUS files'
	end

c***************************************************************

