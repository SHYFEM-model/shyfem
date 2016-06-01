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
c 11.09.2015    ggu     split feature added
c 17.03.2016    ggu     outformat git added
c
c**************************************************************

	program ous2shy

	use clo
	use elabutil
	use elabtime

	use basin
        use mod_depth
        use mod_hydro
        use mod_hydro_baro
        use evgeom
        use levels
        use shyutil

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
	character*80 basnam,simnam
	real rnull
	real zmin,zmax
	real umin,umax,vmin,vmax
	real volume,area
	double precision dtime

	integer iapini
	integer ifem_open_file

c--------------------------------------------------------------

        nwrite = 0
	nread=0
	nelab=0
	nrec=0
	rnull=0.
	rnull=-999.
	bopen = .false.

	ks = -1			!write special node
	iano = -1		!no computation for this area code

c--------------------------------------------------------------
c set command line parameters
c--------------------------------------------------------------

	call elabutil_init('OUS')

	!--------------------------------------------------------------
	! open input files
	!--------------------------------------------------------------

	boutput = .true.
	bneedbasin = .true.

	call ap_init_basin

        call open_ous_type('.ous','old',nin)

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

	call init_sigma_info(nlv,hlv)

	if( bneedbasin ) then
	  call outfile_make_hkv(nkn,nel,nen3v,hm3v,hev,hkv)
	  call ilhe2k(nkn,nel,nen3v,ilhv,ilhkv)
	end if

	if( bverb ) call depth_stats(nkn,nlvdi,ilhkv)

	!--------------------------------------------------------------
	! time management
	!--------------------------------------------------------------

	call ous_get_date(nin,date,time)
	call elabtime_date_and_time(date,time)

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

	if( outformat == 'gis' ) call gis_write_connect

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
	  nread=nread+1
	  nrec = nrec + 1

	  nelab=nelab+1

	  if( .not. bquiet ) then
	    dline = ' '
	    if( bdate ) call dtsgf(it,dline)
	    write(6,*) 'time : ',it,'  ',dline
	  end if

	  if( boutput ) then
            nwrite = nwrite + 1
	    if( bverb ) write(6,*) 'writing to output '
	    call transfer_uvz2(nlvdi,nndim,nvar,vars,zv,ze,uv,vv)
	    call shy_write_hydro_records(id,dtime,nlvdi,zv,ze,uv,vv)
	  end if

	end do		!time do loop

c--------------------------------------------------------------
c end of loop on data
c--------------------------------------------------------------

c--------------------------------------------------------------
c write final message
c--------------------------------------------------------------

	write(6,*)
	write(6,*) nread,' total records read'
	!write(6,*) nrec ,' unique time records read'
	write(6,*) nelab,' records elaborated'
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
