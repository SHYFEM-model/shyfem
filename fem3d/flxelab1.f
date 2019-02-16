
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
c 18.11.1998    ggu     check dimensions with dimnos
c 06.04.1999    ggu     some cosmetic changes
c 03.12.2001    ggu     some extra output -> place of min/max
c 09.12.2003    ggu     check for NaN introduced
c 07.03.2007    ggu     easier call
c 08.11.2008    ggu     do not compute min/max in non-existing layers
c 07.12.2010    ggu     write statistics on depth distribution (depth_stats)
c 06.05.2015    ggu     noselab started
c 05.06.2015    ggu     many more features added
c 05.10.2015    ggu     started flxelab
c 26.10.2017    ggu     various user related improvements, new flx subroutines
c
c**************************************************************

	subroutine flxelab

	use clo
	use elabutil
	use elabtime
	use shyfem_strings

        use basin
        use mod_depth
        use evgeom
        use levels

c elaborates flx file

	implicit none

	integer, allocatable :: kflux(:)
	integer, allocatable :: nlayers(:)
	integer, allocatable :: nsnodes(:,:)
	integer, allocatable :: ivars(:)
	real, allocatable :: fluxes(:,:,:)
	real, allocatable :: fluxes0(:)
	character*80, allocatable :: strings(:)

	integer, allocatable :: naccu(:)
	double precision, allocatable :: accum(:,:,:)

	logical b3d
	integer nread,nelab,nrec,nin,nout
	integer nvers
	integer nknnos,nelnos,nvar
	integer ierr
	integer ivar,iv,ivarnew,ivarfirst
	integer ii,i,j,l,k,lmax,node,nn,n1,n2
	integer ip,nb,naccum
	integer kfluxm
	integer idtflx,nlmax,nsect
	integer nscdi,nfxdi
	integer date,time
	character*80 title,name,femver
	character*20 dline
	character*80 basnam,simnam
        character*10 :: format,range
        character*10 :: short
        character*40 :: full
	real rnull
	real cmin,cmax,cmed,vtot
	real zmin,zmax,zmed
	real href,hzmin
	real vmin,vmax,vmed
	double precision atime,atime0,atfirst,atold,atnew,atlast,atwrite
	character*1 :: what(5) = (/'u','v','z','m','a'/)
	character*80 file

	integer iapini,ifileo
	integer ifem_open_file

c--------------------------------------------------------------

	nread=0
	nelab=0
	nrec=0
	nout=0
	rnull=0.
	rnull=-1.
	bopen = .false.

c--------------------------------------------------------------
c set command line parameters
c--------------------------------------------------------------

	call elabutil_init('FLX','flxelab')

	!--------------------------------------------------------------
	! open input files
	!--------------------------------------------------------------

        call clo_reset_files
        call clo_get_next_file(file)
        nin = ifileo(0,file,'unform','old')

        call flx_is_flx_file(nin,nvers)
        if( nvers .le. 0 ) then
          write(6,*) 'nvers: ',nvers
          stop 'error stop flxelab: not a valid flx file'
        end if

	call flx_read_header(nin,nvers,nsect,kfluxm,idtflx
     +                                  ,nlmax,nvar,ierr)
	if( ierr /= 0 ) goto 93

	nscdi = nsect
	nfxdi = kfluxm
	nlvdi = nlmax
	nlv = nlmax

	allocate(kflux(kfluxm))
	allocate(nlayers(nscdi))
	allocate(ivars(nvar))
	allocate(nsnodes(2,nscdi))
	allocate(strings(nscdi))
	allocate(fluxes(0:nlvdi,3,nscdi))
	allocate(fluxes0(nscdi))

        call flx_read_header2(nin,nvers,nsect,kfluxm
     +                          ,kflux,nlayers
     +                          ,atime0,title,femver,strings,ierr)
	if( ierr /= 0 ) goto 92

	call set_nodes(nsect,kfluxm,kflux,nsnodes)

	if( .not. bquiet ) then
          write(6,*) 'nvers      : ',nvers
          write(6,*) 'nsect      : ',nsect
          write(6,*) 'kfluxm     : ',kfluxm
          write(6,*) 'idtflx     : ',idtflx
          write(6,*) 'nlmax      : ',nlmax
          write(6,*) 'nvar       : ',nvar 
          write(6,*) 'title      : ',trim(title)
          write(6,*) 'femver     : ',trim(femver)
          write(6,*) 'Sections contained in file:'
          write(6,*) ' i   nodes  layers  description'
          do i=1,nsect
	    nn = 1 + nsnodes(2,i) - nsnodes(1,i)
            write(6,1000) i,nn,nlayers(i),'  ',trim(strings(i))
 1000       format(i3,i8,i8,a,a)
          end do

	  if( bverb ) then
           write(6,*) 'Extra information on sections: '
           do i=1,nsect
	    n1 = nsnodes(1,i)
	    n2 = nsnodes(2,i)
	    nn = 1 + n2 - n1
            write(6,*) i,nn,(kflux(j),j=n1,n2)
	   end do
	  end if
	end if

	if( binfo ) return

	b3d = nlmax > 1

	!--------------------------------------------------------------
	! time management
	!--------------------------------------------------------------

	date = 0
	time = 0
	call elabtime_date_and_time(date,time)
	call elabtime_set_minmax(stmin,stmax)

	call flx_peek_record(nin,nvers,atime,ivar,ierr)
	if( ierr /= 0 ) goto 91
        atfirst = atime
        atold = atime

        !--------------------------------------------------------------
        ! see what is in the file
        !--------------------------------------------------------------

        if( .not. bquiet .and. nvar > 0 ) then
          write(6,*) 'Variables contained in file:'
          write(6,*) ' i ivar  short     full'
	end if

	ivarfirst = 0
        do iv=1,nvar
          call flx_read_record(nin,nvers,atime
     +                  ,nlvdi,nsect,ivar
     +                  ,nlayers,fluxes,ierr)
          if( ierr /= 0 ) goto 91
	  ivars(iv) = ivar
	  if( iv == 1 ) ivarfirst = ivar
          call strings_get_short_name(ivar,short)
          call strings_get_full_name(ivar,full)
          if( .not. bquiet ) then
	    write(6,'(i3,i5,a,a,a)') iv,ivar,'  ',short,full
	  end if
        end do

        do iv=1,nvar
          backspace(nin)
        end do

	if( binfo ) return

	!--------------------------------------------------------------
	! averaging
	!--------------------------------------------------------------

	!call elabutil_set_averaging(nvar)

	btrans = .false.
	if( btrans ) then
	  allocate(naccu(istep))
	  allocate(accum(nlvdi,nkn,istep))
	else
	  allocate(naccu(1))
	  allocate(accum(1,1,1))
	end if
	naccum = 0
	naccu = 0
	accum = 0.

	!write(6,*) 'mode: ',mode,ifreq,istep

	!--------------------------------------------------------------
	! initialize variables
	!--------------------------------------------------------------

	fluxes = 0
	fluxes0 = 0

	!--------------------------------------------------------------
	! open output file
	!--------------------------------------------------------------

	boutput = boutput .or. btrans
	bopen = boutput
	if( bsplit ) then
	  boutput = .false.
	  bopen = .false.
	end if

	if( bopen ) then
	  nb = ifileo(0,'out.flx','unform','new')
	  call flx_write_header(nb,0,nsect,kfluxm
     +			,idtflx,nlmax,nvar,ierr)
	  if( ierr /= 0 ) goto 99
	  call flx_write_header2(nb,0,nsect,kfluxm
     +			,kflux,nlayers
     +			,atime0,title,femver,strings,ierr)
	  if( ierr /= 0 ) goto 99
	end if

c--------------------------------------------------------------
c loop on data
c--------------------------------------------------------------

	atime = atold
	atwrite = -1
	if( .not. bquiet ) write(6,*)

	do

	 atold = atime

	 call flx_read_record(nin,nvers,atime,nlvdi,nsect,ivar
     +			,nlayers,fluxes,ierr)
         if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
         if(ierr.ne.0) exit

	 nread = nread + 1
	 if( ivar == ivarfirst ) then
	   nrec = nrec + 1
	   iv = 0
	 end if
	 iv = iv + 1

	 atlast = atime
	 call flx_peek_record(nin,nvers,atnew,ivarnew,ierr)
	 if( ierr .ne. 0 ) atnew = atime

         if( elabtime_over_time(atime,atnew,atold) ) exit
         if( .not. elabtime_in_time(atime,atnew,atold) ) cycle
	 !if( .not. elabtime_check_time(atime,atnew,atold) ) cycle

	  if( ivar == 0 ) nelab=nelab+1

	  if( bverb .and. atwrite /= atime ) then
            call dts_format_abs_time(atime,dline)
            write(6,*) 'time : ',atime,'  ',dline
            atwrite = atime
	  end if

	  if( bwrite ) then
	    fluxes0(:) = fluxes(0,1,:)
            call strings_get_short_name(ivar,short)
	    call minmaxmed(fluxes0,nsect,vmin,vmax,vmed)
            write(6,*) trim(short)//' min,max,aver : ',vmin,vmax,vmed
	  end if

	  if( bsplit ) then
            call split_fluxes_0d(atime,nlvdi,nsect,nvar,iv,ivar,fluxes)
            call split_fluxes_2d(atime,nlvdi,nsect,nvar,iv,ivar
     +				,fluxes,bsplitall)
	    if( b3d ) then
	      call split_fluxes_3d(atime,nlvdi,nsect,nvar,iv,ivar
     +				,nlayers,fluxes,bsplitall)
	    end if
	  end if

	  if( boutput ) then
	    if( ivar == 0 ) nout = nout + 1
	    if( bverb ) write(6,*) 'writing output var: ',ivar
            call flx_write_record(nb,0,atime
     +			,nlvdi,nsect,ivar,nlayers,fluxes,ierr)
	    if( ierr /= 0 ) goto 99
	  end if

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
!	      call nos_time_aver(-mode,ip,ifreq,istep,nkn,nlvdi
!     +					,naccu,accum,cv3,boutput)
!	      if( bsumvar ) ivar = 30
!              call nos_write_record(nb,it,ivar,nlvdi,ilhkv,cv3,ierr)
              if( ierr .ne. 0 ) goto 99
	    end if
	  end do
	end if

c--------------------------------------------------------------
c write final message
c--------------------------------------------------------------

	if( .not. bsilent ) then
          write(6,*)
	  call dts_format_abs_time(atfirst,dline)
          write(6,*) 'first time record: ',dline
	  call dts_format_abs_time(atlast,dline)
          write(6,*) 'last time record:  ',dline
	  write(6,*)
	  write(6,*) nread,' data records read'
	  write(6,*) nrec ,' time records read'
	  write(6,*) nelab,' time records elaborated'
	  write(6,*) nout ,' time records written to file'
	  write(6,*)
	end if

	if( .not. bquiet ) then
	 if( bsplit ) then
          write(6,*) 'output written to following files: '
          write(6,*) '  disch.what.dim.sect'
	  if( bsplitall ) then
	   write(6,*) 'disch is: '
	   write(6,*) '  disch or disch_total    total discharge'
	   write(6,*) '  disch_positive          positive discharge'
	   write(6,*) '  disch_negative          negative discharge'
	   write(6,*) '  disch_absolute          absolute discharge'
	  end if
          write(6,*) 'what is one of the following:'
	  call write_vars(nvar,ivars)
          write(6,*) 'dim is 0d, 2d or 3d'
          write(6,*) '  0d for depth averaged variables (one column)'
          write(6,*) '  2d for depth averaged variables'
          write(6,*) '  3d for output at each layer'
	  call compute_range(nsect,range)
	  write(6,'(a,a)') ' sect is consecutive section numbering: '
     +				,trim(range)
	  write(6,*) 'discharge data is normally total discharge'
	  write(6,*) 'for 2d files the 4 columns are: '//
     +			'total,positive,negative,absolute'
	  write(6,*) 'for 3d files the columns are the total '//
     +			'discharge for each layer'
	 else if( boutput ) then
	  write(6,*) 'output written to file out.flx'
	 end if
	end if

	!call ap_get_names(basnam,simnam)
	!write(6,*) 'names used: '
	!write(6,*) 'simul: ',trim(simnam)

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	stop
   91	continue
	write(6,*) 'error reading first data record'
	write(6,*) 'maybe the file is empty'
	stop 'error stop flxelab: empty record'
   92	continue
	write(6,*) 'error reading second header'
	stop 'error stop flxelab: error in header'
   93	continue
	write(6,*) 'error reading first header'
	stop 'error stop flxelab: error in header'
   99	continue
	write(6,*) 'error writing to file unit: ',nb
	stop 'error stop flxelab: write error'
	end

c***************************************************************
c***************************************************************
c***************************************************************

        subroutine split_fluxes_2d(atime,nlvddi,nsect
     +				,nvar,iv,ivar,fluxes,bext)

c writes 2d fluxes to file (only for ivar=0)

        implicit none

	double precision atime
        integer nlvddi                  !vertical dimension
        integer nsect                   !total number of sections
	integer nvar
	integer iv
        integer ivar                    !type of variable (0: water fluxes)
        real fluxes(0:nlvddi,3,nsect)   !fluxes
	logical bext			!extended fluxes

        integer i,j,iu
        real ptot(4,nsect)
	character*80 name,filename
	character*20 dline
        integer, save, allocatable :: iusplit(:,:)
        integer, save :: iubox = 0
        integer, save :: icall = 0

        if( icall == 0 ) then
	  allocate(iusplit(nsect,nvar))
	  iusplit = 0
	  iubox = 0
	  if( bext ) then
            name = 'disch.box.2d.txt'
	    call get_new_unit(iubox)
            open(iubox,file=name,form='formatted',status='unknown')
	  end if
	end if

	if( iusplit(1,iv) == 0 ) then
	  call ivar2filename(ivar,filename)
	  do j=1,nsect
            call make_iunit_name('disch.',filename,'2d',j,iu)
            iusplit(j,iv) = iu
          end do
        end if

	icall = icall + 1

        call dts_format_abs_time(atime,dline)

        do j=1,nsect
          ptot(1,j) = fluxes(0,1,j)                     !total
          ptot(2,j) = fluxes(0,2,j)                     !positive
          ptot(3,j) = fluxes(0,3,j)                     !negative
          ptot(4,j) = fluxes(0,2,j) + fluxes(0,3,j)     !absolute
	  iu = iusplit(j,iv)
	  call iusplit_info(iu,j,iv,nsect,nvar)
	  write(iu,'(a20,4f14.3)') dline,(ptot(i,j),i=1,4)
        end do

c next is box format for Ali

	if( bext .and. ivar == 0 ) then
          write(iubox,*) atime
          write(iubox,*) 0
          write(iubox,*) nsect
          do j=1,nsect
            write(iubox,*) 0,0,(ptot(i,j),i=1,3)
          end do
	end if

        end

c****************************************************************

        subroutine split_fluxes_3d(atime,nlvddi,nsect,nvar,iv,ivar
     +				,nlayers,fluxes,bext)

c writes 3d fluxes to file

	use shyfem_strings

        implicit none

	double precision atime
        integer nlvddi                  !vertical dimension
        integer nsect                   !total number of sections
	integer nvar
	integer iv
        integer ivar                    !type of variable (0: water fluxes)
        integer nlayers(nsect)          !max layers for section
        real fluxes(0:nlvddi,3,nsect)   !fluxes
	logical bext			!extended write

        integer lmax,l,j,i
	integer iu,iunit
        integer, save, allocatable :: iusplit(:,:,:)
        integer, save :: icall = 0
	integer, save :: imin,imax
	real port(0:nlvddi,5,nsect)
	character*80 format,filename
	character*10 short
	character*20 dline
	character*14, save :: what(5) = (/ 
     +				 'disch_total   '
     +				,'disch_positive'
     +				,'disch_negative'
     +				,'disch_absolute'
     +				,'disch         '
     +				/)

        if( icall == 0 ) then
	  imin = 5
	  imax = 5
	  if( bext ) then		!write all fluxes, not only total
	    imin = 1
	    imax = 4
	  end if
	  allocate(iusplit(nsect,nvar,5))
	  iusplit = 0
	end if

	if( iusplit(1,iv,imin) == 0 ) then
	  call ivar2filename(ivar,filename)
	  do j=1,nsect
	    do i=imin,imax
              call make_iunit_name(what(i),'.'//filename,'3d',j,iu)
	      iusplit(j,iv,i) = iu
	    end do
	  end do
	end if

	icall = icall + 1

        call dts_format_abs_time(atime,dline)

	if( bext ) then
	  port(:,1:3,:) = fluxes(:,1:3,:)
	  port(:,4,:) = fluxes(:,2,:)+fluxes(:,3,:)	!absolute
	else
	  port(:,5,:) = fluxes(:,1,:)			!equal to 1 (total)
	end if

        do j=1,nsect
          lmax = nlayers(j)
	  write(format,'(a,i3,a)') '(a20,',lmax,'f12.3)'
	  do i=imin,imax
	    iu = iusplit(j,iv,i)
	    write(iu,format) dline,(port(l,i,j),l=1,lmax)
	  end do
        end do

	return
   99	continue
        end

c****************************************************************

        subroutine split_fluxes_0d(atime,nlvddi,nsect,nvar,iv,ivar
     +					,fluxes)

c splits barotropic values of discharges - old version

        implicit none

	double precision atime
	integer nlvddi
	integer nsect
	integer nvar
	integer iv
	integer ivar
        real fluxes(0:nlvddi,3,nsect)   !fluxes

        integer j,iu
	integer, save :: icall = 0
        integer, save, allocatable :: iusplit(:)
	character*20 dline

        if( ivar .ne. 0 ) return
        if( iv .ne. 1 ) return

	if( icall == 0 ) then
	  allocate(iusplit(nsect))
	  iusplit = 0
	  do j=1,nsect
	    call make_iunit_name('disch','.mass','0d',j,iu)
	    iusplit(j) = iu
	  end do
	end if

	icall = icall + 1

        call dts_format_abs_time(atime,dline)

	do j=1,nsect
	  iu = iusplit(j)
	  !write(iu,*) it,(fluxes(0,ii,j),ii=1,3)
	  write(iu,*) dline,fluxes(0,1,j)
	end do

        end

c*******************************************************************

	subroutine set_nodes(nsect,kfluxm,kflux,nsnodes)

	implicit none

	integer nsect,kfluxm
	integer kflux(kfluxm)
	integer nsnodes(2,nsect)

	logical bsect
	integer i,k,is,nn

	nsnodes = 0

	is = 0
	bsect = .false.

	do i=1,kfluxm
	  k = kflux(i)
	  if( bsect ) then
	    if( k <= 0 ) then
	      bsect = .false.
	      nsnodes(2,is) = i-1
	    end if
	  else
	    if( k > 0 ) then
	      is = is + 1
	      if( is > nsect ) goto 99
	      bsect = .true.
	      nsnodes(1,is) = i
	    end if
	  end if
	end do

	if( bsect ) nsnodes(2,is) = i-1

	if( is /= nsect ) goto 99

	!do is=1,nsect
	!  nn = 1 + nsnodes(2,is) - nsnodes(1,is)
	!  write(6,*) nn,nsnodes(:,is)
	!end do

	return
   99	continue
	write(6,*) 'kflux does not contain expected sections: '
	write(6,*) 'is,nsect: ',is,nsect
	stop 'error stop set_nodes: nsect'
	end

c*******************************************************************

	subroutine iusplit_info(iu,j,iv,nsect,nvar)

	implicit none

	integer iu
	integer j,iv
	integer nsect,nvar

	if( iu < 10 .or. iu > 1000 ) then
	  write(6,*) 'error in unit number: ',iu
	end if

	if( j < 1 .or. j > nsect ) then
	  write(6,*) 'error in section number: ',j,nsect
	end if

	if( iv < 1 .or. iv > nvar ) then
	  write(6,*) 'error in variable number: ',iv,nvar
	end if

	end

c*******************************************************************

