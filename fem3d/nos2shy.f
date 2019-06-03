
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

c revision log :
c
c 18.11.1998	ggu	check dimensions with dimnos
c 06.04.1999	ggu	some cosmetic changes
c 03.12.2001	ggu	some extra output -> place of min/max
c 09.12.2003	ggu	check for NaN introduced
c 07.03.2007	ggu	easier call
c 08.11.2008	ggu	do not compute min/max in non-existing layers
c 07.12.2010	ggu	write statistics on depth distribution (depth_stats)
c 06.05.2015	ggu	noselab started
c 05.06.2015	ggu	many more features added
c 10.09.2015	ggu	std and rms for averaging implemented
c 11.09.2015	ggu	write in gis format
c 23.09.2015	ggu	handle more than one file (look for itstart)
c 19.02.2016	ggu	bug fixes for bsumvar and mode==2 (sum)
c 22.02.2016	ggu	handle catmode
c 28.04.2016	ggu	changed VERS_7_5_9
c 25.05.2016	ggu	changed VERS_7_5_10
c 07.06.2016	ggu	changed VERS_7_5_12
c 14.11.2017	ggu	changed VERS_7_5_36
c 22.02.2018	ggu	changed VERS_7_5_42
c 06.07.2018	ggu	changed VERS_7_5_48
c 16.02.2019	ggu	changed VERS_7_5_60
c
c**************************************************************

	program nos2shy

	use clo
	!use elabutil
	use elabtime
	use shyfile
	use shympi

        use basin
        use mod_depth
        use evgeom
        use levels
	use iso8601

c elaborates nos file

	implicit none

	integer, parameter :: ndim = 1000
	integer iusplit(ndim)

	real, allocatable :: cv3(:,:)
	real, allocatable :: cv3all(:,:,:)
	real, allocatable :: vol3(:,:)

	integer, allocatable :: ivars(:)
	real, allocatable :: hl(:)

	logical bopen,bneedbasin,boutput,bquiet,bverb,bsilent
	logical bdate,bdstring
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
	character*20 vers
	character*80 version
	character*80 dstring
	character*80 basnam,simnam
	character*80 basfile,simfile
	real rnull
	real cmin,cmax,cmed,vtot
	double precision dtime,dtfirst,dtlast

	integer iapini
	integer ifem_open_file
	logical concat_cycle

c--------------------------------------------------------------

	ifile = 1
	nread=0
	nwrite=0
	nelab=0
	nrec=0
	rnull=0.
	rnull=-1.
	bopen = .false.
	dtfirst = 0.
	dtlast = 0.

c--------------------------------------------------------------
c set command line parameters
c--------------------------------------------------------------

        call get_shyfem_version(vers)
        version = '2.0' // ' (SHYFEM version ' // trim(vers) // ')'

        call clo_init('nos2shy','nos-file bas-file',version)
        call clo_add_info('transforms nos file into shy file')

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
	 call shyfem_copyright('nos2shy - transforms NOS to SHY files')
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

	if( .not. bquiet ) then
	  call write_nos_info(nin)
          write(6,*) 'levels: '
          write(6,'(5g14.6)') (hlv(l),l=1,nlv)
	end if

	call init_sigma_info(nlv,hlv)

	if( bneedbasin ) then
	  call outfile_make_hkv(nkn,nel,nen3v,hm3v,hev,hkv)
          call ilhk2e(nkn,nel,nen3v,ilhkv,ilhv)
          call adjust_layer_index(nel,nlv,hev,hlv,ilhv)
	  call init_volume(nlvdi,nkn,nel,nlv,nen3v,ilhkv
     +                          ,hlv,hev,hl,vol3)
	end if

	if( bverb ) call depth_stats(nkn,nlvdi,ilhkv)

        call shympi_init(.false.)               !call after basin has been read
        call shympi_set_hlv(nlv,hlv)

	!--------------------------------------------------------------
	! time management
	!--------------------------------------------------------------

	call nos_get_date(nin,date,time)
	if( date < 10000 ) date = date * 10000 + 101
	bdate = ( date > 0 ) 
	bdstring = ( dstring /= ' ' )
	if( .not. bdate .and. .not. bdstring ) then
	  write(6,*) 'no date in file... please specify with -date'
	  stop 'error stop nos2shy: no date'
	else if( bdate .and. bdstring ) then
	  write(6,*) 'date in file... cannot specify -date'
	  stop 'error stop nos2shy: no -date possible'
	else if( bdstring ) then
	  call string2date(dstring,date,time,ierr)
	  if( ierr /= 0 ) then
	    write(6,*) 'cannot parse date string: ',trim(dstring)
	    stop 'error stop nos2shy: error in date'
	  end if
	end if
	call nos_set_date(nin,date,time)
	call elabtime_date_and_time(date,time)
	call date2string(date,time,dline)
	if( .not. bquiet ) write(6,*) 'reference date used: ',dline
	bdate = .true.

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
	 if( nrec == 0 ) dtfirst = dtime
	 nrec = nrec + 1
	 nelab = nelab + 1

	 do i=1,nvar

	  ivar = ivars(i)
	  cv3(:,:) = cv3all(:,:,i)

	  if( bverb ) then
	    call dtsgf(it,dline)
	    write(6,*) 'time : ',it,'  ',dline,'   ivar : ',ivar
	  end if

	  if( boutput ) then
	    nwrite = nwrite + 1
	    if( bverb ) write(6,*) 'writing to output: ',ivar
	    call shy_write_scalar_record(id,dtime,ivar,nlv,cv3)
	  end if

	 end do		!loop on ivar
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
	write(6,*) nread, ' records read'
	write(6,*) nrec , ' unique time records read'
	write(6,*) nelab, ' records elaborated'
	write(6,*) ifile, ' files read'
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
	  call nos_is_nos_file(nin,nvers)
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
	  write(6,*) 'no NOS file given on command line'
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
	write(6,*) 'more than one NOS file on command line'
	write(6,*) 'can handle only one file'
	stop 'error stop compat_get_files: NOS files'
	end

c***************************************************************
