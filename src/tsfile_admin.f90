! $Id: newpr.f,v 1.24 2010-02-22 15:38:36 georg Exp $
!
! reading and interpolation of external files
!
! contents :
!
! revision log :
!
! 29.10.2012    ggu     created from scratch
! 17.06.2013    ggu     do not pass function into subroutine
! 02.07.2014    ggu     new framework finished
! 10.07.2014    ggu     only new file format allowed
!
!*******************************************************************	
!*******************************************************************	
!*******************************************************************	
!-------------------------------------------------------------------
        module tsfile_admin
!-------------------------------------------------------------------
        contains
!-------------------------------------------------------------------

	subroutine ts_file_open(name,it,nkn,nlv,iunit)
	integer iunit(3)
	character*(*) name
	  call ts_file_open_1(name,it,nkn,nlv,iunit)
	end

	subroutine ts_file_descrp(iunit,name)
	use intp_fem_file
	integer iunit(3)
	character*(*) name
	  call iff_set_description(iunit(1),0,name)
	end

	subroutine ts_next_record(it,iunit,nlvddi,nkn,nlv,value)
	include 'param.h'
	integer iunit(3)
	double precision value(nlvddi,nkn)
	  call ts_next_record_1(it,iunit,nlvddi,nkn,nlv,value)
	end

	subroutine ts_file_close(info)
	integer info(3)
	  call ts_file_close_1(info)
	end

!*******************************************************************	
!*******************************************************************	
!*******************************************************************	
!*******************************************************************	

	subroutine ts_file_open_1(file,it,np,nlv,iunit)

! opens T/S file

	use intp_fem_file

	implicit none

	include 'param.h'

	character*(*) file		!name of file
	integer it
	integer np			!number of points expected
	integer nlv
	integer iunit(3)		!unit number (return)

	integer nvar,nexp,lexp,nintp
	integer id
	double precision dtime
	integer nodes(1)
	double precision vconst(1)

	dtime = it
	nvar = 1
	nexp = np
	lexp = nlv
	nintp = 2
	nodes = 0
	vconst = 0.

!$OMP CRITICAL
	call iff_init(dtime,file,nvar,nexp,lexp,nintp,nodes,vconst,id)
!$OMP END CRITICAL

	iunit(1) = id

	return
   99	continue
	write(6,*) 'Cannot open file: ',file
	stop 'error stop ts_file_open: error file open'
	end

!*******************************************************************

        subroutine iff_ts_mpi(dtime,file,nvar,nexp,nmpi,lexp,nintp,vconst,id)

        use shympi
        use femfile
        use intp_fem_file
        use tsfile

        implicit none

        double precision dtime  !initial time
        character*(*) file      !file name
        integer nvar            !expected number of variables (might change)
        integer nexp            !expected number of points
        integer nmpi            !expected number of points mpi
        integer lexp            !expected max vertical data (0 for 2D)
        integer nintp           !requested interpolation (2 linear, 4 cubic)
        double precision vconst(nvar)       !constant values in case no file is given
        integer id              !identification of file info (return)

! nvar is return value (if file can be read)

        integer iformat,iunit
        integer ierr,np,i
        integer nvar_orig,nvar_read
        integer datetime(2)
        integer ntype,itype(2)
        integer id0,ibc
        logical breg
        logical bok
        logical bts,bfem,bnofile,bfile,berror,bnosuchfile,boperr
        type(info), pointer :: p

        !---------------------------------------------------------
        ! get new id for file
        !---------------------------------------------------------

        idlast = idlast + 1
        if( idlast > 300 ) then
          stop 'error stop iff_ts_mpi: too many files opened'
        end if
        id = idlast

        call iff_init_entry(id)

        !---------------------------------------------------------
        ! check input parameters
        !---------------------------------------------------------

        if( nvar < 1 ) goto 97
        if( nexp < 1 ) goto 97
        if( lexp < 0 ) goto 97
        if( nintp < 1 ) goto 97

        !---------------------------------------------------------
        ! get info on file
        !---------------------------------------------------------

        nvar_orig = nvar

        call iff_get_file_info(file,np,nvar_read,ntype,iformat)

        bnofile = iformat == iform_none                 !no file given
        bfile = .not. bnofile                           !file has been given
        bnosuchfile = iformat == iform_no_such_file     !file not existing
        !bnofile = iformat < 0
        !berror = iformat < 0
        boperr = iformat == iform_error_opening         !error opening
        berror = iformat == iform_error                 !error file

        id0 = 0
        if( boperr ) then       !see if we can take data from other file
          id0 = iff_find_id_to_file(file)
          id0 = 0               !not working
          if( id0 > 0 ) then
            boperr = .false.
            iformat = pinfo(id0)%iformat
            ntype = pinfo(id0)%ntype
          end if
        end if

        bts = iformat == iform_ts                       !file is TS
        bfem = iformat >= 0 .and. iformat <= 2

        call fem_file_make_type(ntype,2,itype)
        breg = itype(2) > 0

        if( file /= ' ' .and. bnofile ) goto 99
        if( bnosuchfile ) goto 99
        if( bfile .and. berror ) goto 93
        if( bfile .and. boperr ) goto 91
        if( bfile .and. np < 1 ) goto 96
        if( .not. breg .and. np > 1 .and. np /= nexp ) goto 96

        if( bfile ) then
          if( nvar_orig /= nvar_read ) goto 92
        end if
        nvar = nvar_orig

        !---------------------------------------------------------
        ! store information
        !---------------------------------------------------------

        p => pinfo(id)

        pinfo(id)%iunit = -1
        pinfo(id)%nvar = nvar
        pinfo(id)%nintp = nintp
        pinfo(id)%file = file
        pinfo(id)%nexp = nexp
        pinfo(id)%nmpi = nmpi
        pinfo(id)%lexp = lexp

        pinfo(id)%id0 = id0
        pinfo(id)%iformat = iformat
        pinfo(id)%ntype = ntype
        pinfo(id)%ireg = itype(2)

        if( p%nexp /= nexp ) then       !dummy test - delete later
          stop 'error stop iff_ts_mpi: internal error (99)'
        end if

        !---------------------------------------------------------
        ! get data description and allocate data structure
        !---------------------------------------------------------

        allocate(pinfo(id)%strings_file(nvar))

        if( id0 > 0 ) then
          pinfo(id)%strings_file = pinfo(id0)%strings_file
          pinfo(id)%bonepoint = pinfo(id0)%bonepoint
        else if( bfem ) then
          call fem_file_get_data_description(file,pinfo(id)%strings_file,ierr)
          if( ierr /= 0 ) goto 98
        else if( bts ) then
          pinfo(id)%bonepoint = .true.
          pinfo(id)%strings_file = ' '
        else if( bnofile ) then         !constant
          pinfo(id)%nintp = 0
          pinfo(id)%ilast = 1
          pinfo(id)%bonepoint = .true.
          pinfo(id)%strings_file = ' '
          call iff_allocate_fem_data_structure(id)
          pinfo(id)%time(1) = 0.
          do i=1,nvar
            pinfo(id)%data(1,1,i,1) = vconst(i)
          end do
        else
          stop 'error stop iff_ts_mpi: internal error (3)'
        end if

        !---------------------------------------------------------
        ! finally open files
        !---------------------------------------------------------

        if( bnofile ) then
          return
        else if( bfem ) then
          !call fem_file_read_open(file,nexp,iformat,iunit)
          call fem_file_read_open(file,np,iformat,iunit)
        else if( bts ) then
          call ts_open_file(file,nvar,datetime,iunit)
          pinfo(id)%datetime = datetime
        else
          stop 'error stop iff_ts_mpi: internal error (3)'
        end if

        if( iunit < 0 ) goto 99
        if( iunit == 0 ) goto 90
        pinfo(id)%iunit = iunit

        if(my_id.eq.0) write(6,*) 'file opened: ',id,trim(file)

        !---------------------------------------------------------
        ! populate data base
        !---------------------------------------------------------

        call iff_fill_records_mpi(id)

        !---------------------------------------------------------
        ! end of routine
        !---------------------------------------------------------

        return
   90   continue
        write(6,*) 'error in opening file: ',trim(file)
        write(6,*) 'iformat = ',iformat
        stop 'error stop iff_ts_mpi, 90'
   91   continue
        write(6,*) 'error opening file: ',trim(file)
        id0 = iff_find_id_to_file(file)
        if( id0 > 0 ) then
          ibc = pinfo(id0)%ibc
          write(6,*) 'the file is already open on boundary ',ibc
        else
          write(6,*) 'iformat = ',iformat
        end if
        stop 'error stop iff_ts_mpi, 91'
   92   continue
        write(6,*) 'error in file: ',trim(file)
        write(6,*) 'file does not contain correct number of variables'
        write(6,*) 'nvar file = ',nvar_read
        write(6,*) 'nvar expected = ',nvar_orig
        stop 'error stop iff_ts_mpi, 92'
   93   continue
        write(6,*) 'error in file: ',trim(file)
        write(6,*) 'iformat = ',iformat
        stop 'error stop iff_ts_mpi, 93'
   96   continue
        write(6,*) 'file does not contain expected data size'
        write(6,*) 'nexp,np: ',nexp,np
        call iff_print_file_info(id)
        stop 'error stop iff_ts_mpi, 96'
   97   continue
        write(6,*) 'error in input parameters of routine: '
        write(6,*) 'file: ',trim(file)
        write(6,*) 'nvar: ',nvar
        write(6,*) 'nexp,lexp: ',nexp,lexp
        write(6,*) 'nintp: ',nintp
        write(6,*) 'nkn_fem: ',nkn_fem
        call iff_print_file_info(id)
        stop 'error stop iff_ts_mpi, 97'
   98   continue
        write(6,*) 'error reading data description of file: ',trim(file)
        call iff_print_file_info(id)
        stop 'error stop iff_ts_mpi, 98'
   99   continue
        write(6,*) 'no such file: ',trim(file)
        write(6,*) 'iformat = ',iformat
        stop 'error stop iff_ts_mpi, 99'

        end subroutine iff_ts_mpi

!****************************************************************

        function iff_read_ts_record_mpi(id,dtime)

        use intp_fem_file
        use shympi
        use timing

        logical iff_read_ts_record_mpi
        integer id
        double precision dtime, time1

        if( iff_read_header(id,dtime) ) then

          if(ln_timing) time1 = shympi_wtime()

          call iff_read_ts_mpi(id,dtime)

          if(ln_timing) io_time = io_time + shympi_wtime() - time1

          iff_read_ts_record_mpi = .true.
        else
          call iff_close_file(id)
          iff_read_ts_record_mpi = .false.
        end if

        end function iff_read_ts_record_mpi

!****************************************************************

        subroutine iff_fill_records_mpi(id)

        use intp_fem_file
        integer id

        double precision dtime,dtime2
        logical bok

        if( .not. iff_read_ts_record_mpi(id,dtime) ) goto 99

        bok = iff_peek_next_record(id,dtime2)

        pinfo(id)%nintp = 0
        pinfo(id)%ilast = 1
        call iff_allocate_fem_data_structure(id)
        call iff_space_interpolate(id,1,dtime)
        call iff_close_file(id)

        return
   99   continue
        call iff_print_file_info(id)
        write(6,*) 'error reading first record of file'
        stop 'error stop iff_fill_records_mpi: read error'

        end  subroutine iff_fill_records_mpi

!*******************************************************************

        subroutine ts_file_open_mpi(file,it,np,nmpi,nlv,iunit,nvar)

! opens T/S file

        implicit none

        include 'param.h'

        character*(*) file              !name of file
        integer it
        integer np			!number of points expected
	integer nmpi                    !number of points subdomain
        integer nlv
        integer iunit(3)                !unit number (return)

        integer nvar,nexp,lexp,nintp
        integer id
        double precision dtime
        !integer nodes(1)
        double precision vconst(1)

        dtime = it
        !nvar = 1
        nexp = np
        lexp = nlv
        nintp = 2
        !nodes = 0
        vconst = 0.

!$OMP CRITICAL
        call iff_ts_mpi(dtime,file,nvar,nexp,nmpi,lexp,nintp,vconst,id)
!$OMP END CRITICAL

        iunit(1) = id

        return

        end subroutine

!****************************************************************

        subroutine iff_read_ts_mpi(id,dtime)

        use intp_fem_file

        integer id
        double precision dtime

        integer iunit,nvers,np,lmax,nmpi
        integer nlvddi,nvar
        integer ierr,i,iformat
        logical bnofile,bts
        character*60 string

        iformat = pinfo(id)%iformat
        bts = iformat == iform_ts
        bnofile = iformat < 0

        if( bnofile ) return

        iunit = pinfo(id)%iunit
        nvers = pinfo(id)%nvers
        np = pinfo(id)%np
        lmax = pinfo(id)%lmax
        nmpi = pinfo(id)%nmpi
        nvar = pinfo(id)%nvar

        nlvddi = lmax

        if( bts ) then
          ! ts data has already been read
        else
          do i=1,nvar
            call fem_file_read_ts_mpi(iformat,iunit,nvers,np,lmax,nmpi,string,pinfo(id)%ilhkv_file      &
     &                                   ,pinfo(id)%hd_file,nlvddi,pinfo(id)%data_file(1,1,i),ierr)
            if( ierr /= 0 ) goto 99
            if( string /= pinfo(id)%strings_file(i) ) goto 98
          end do
        end if

        return
   98   continue
        write(6,*) 'string description has changed for var ',i
        write(6,*) 'time: ',dtime
        write(6,*) 'old: ',pinfo(id)%strings_file(i)
        write(6,*) 'new: ',string
        call iff_print_file_info(id)
        stop 'error stop iff_read_ts_mpi'
   99   continue
        write(6,*) 'error reading data: ',ierr
        write(6,*) 'time: ',dtime
        call iff_print_file_info(id)
        stop 'error stop iff_read_ts_mpi'

        end subroutine iff_read_ts_mpi

!******************************************************************

      subroutine fem_file_read_ts_mpi(iformat,iunit,nvers,np,lmax,nmpi,string,ilhkv,hd_par,nlvddi,data_par,ierr)

! reads data of the file

        use shympi
        use bnd_geom

        implicit none

        integer iformat         !formatted or unformatted
        integer iunit           !file unit
        integer nvers           !version of file format
        integer np              !size of data (horizontal, nodes or elements)
        integer lmax            !vertical values
        integer nmpi            !size of mpi data (horizontal, nodes or elements)
        character*(*) string    !string explanation
        integer ilhkv(nmpi)     !number of layers in point k (node)
        double precision hd_par(nmpi)               !total depth
        integer nlvddi          !vertical dimension of data
        double precision data_par(nlvddi,nmpi)
        integer ierr            !return error code
        integer h
        logical b2d
        integer k,lm,l,npp
        double precision hdepth
        character*80 text
        integer p,ipext

        ierr = 0
        npp = np

	if( iformat == 1 ) then
          read(iunit,'(a)',err=13) text
          if( nvers >= 3 ) read(iunit,*,err=11) npp,lmax
          if( lmax > nlvddi ) goto 97
          if( np /= npp ) goto 96
          b2d = lmax .le. 1
          if( b2d ) then
            do k=1,np
	      do h=1,nmpi
                if (domain%nodes%globalID(h).eq.k) then
                  read(iunit,*,err=15) data_par(1,h)
                  exit
                end if
                if (h.eq.nmpi) read(iunit,*,err=15)
	      end do
	    end do
          else
            do k=1,np
	      do h=1,nmpi
	        if(k .eq. domain%nodes%globalID(h)) then
        	  read(iunit,*,err=15) lm,hd_par(h),(data_par(l,h),l=1,min(lm,lmax))
		  ilhkv(h) = lm
		  exit
        	  if( lm .gt. lmax ) goto 99
	        else
		  if(h.eq.nmpi) read(iunit,*,err=15) 
	        end if
	      end do
            end do
          end if
        else
          read(iunit,err=13) text
          if( nvers >= 3 ) read(iunit,err=11) npp,lmax
          if( lmax > nlvddi ) goto 97
          if( np /= npp ) goto 96
          b2d = lmax .le. 1
          if( b2d ) then
            do k=1,np
	      do h=1,nmpi
                if (domain%nodes%globalID(h).eq.k) then
                  read(iunit,err=15) data_par(1,h)
                  exit
                end if
                if (h.eq.nmpi) read(iunit,err=15)
              end do
            end do 
          else
            do k=1,np
	      do h=1,nmpi
	        if(k .eq. domain%nodes%globalID(h)) then
        	  read(iunit,err=15) lm,hd_par(h),(data_par(l,h),l=1,min(lm,lmax))
		  ilhkv(h) = lm
		  exit
        	  if( lm .gt. lmax ) goto 99
	        else
		  if(h.eq.nmpi) read(iunit,*,err=15) 
	        end if
	      end do
            end do
          end if
        end if
        if( b2d ) then
          ilhkv = 1
          hd_par = 10000.
        end if
        string = trim(text)

        return
   11   continue
        write(6,*) 'error reading data size'
        ierr = 11
        return
   13   continue
        write(6,*) 'error reading string description'
        ierr = 13
        return
   15   continue
        write(6,*) 'error reading data record'
        ierr = 15
        return
   96   continue
        write(6,*) 'error reading data record: size mismatch'
        write(6,*) 'np,npp: ',np,npp
        ierr = 96
        return
   97   continue
        write(6,*) 'error reading data record: lmax > nlvddi'
        write(6,*) 'lmax,nlvddi: ',lmax,nlvddi
        ierr = 97
        return
   99   continue
        write(6,*) 'error reading data record: too much vertical data'
        write(6,*) 'k,lm,lmax: ',k,lm,lmax
        ierr = 99
        return
      end subroutine

!*******************************************************************	

	subroutine ts_next_record_mpi(it,iunit,nlvddi,nkn,nlv,value,ivar)

        use  conz_util
        use intp_fem_file

	implicit none

	include 'param.h'

	integer it
	integer iunit(3)
	integer nlvddi
	integer nkn
	integer nlv
	double precision value(nlvddi,nkn)

	integer id,ldim,ndim,ivar
        double precision vmin,vmax
	double precision dtime
	character*80 string

!--------------------------------------------------------------
! initialize
!--------------------------------------------------------------

	id = iunit(1)

!--------------------------------------------------------------
! read new data
!--------------------------------------------------------------

	dtime = it
	ndim = nkn
	ldim = nlvddi

	!write(6,*)'reading T/S values', it,dtime

	!call iff_read_and_interpolate(id,dtime)
	call iff_time_interpolate(id,dtime,ivar,ndim,ldim,value)

!--------------------------------------------------------------
! some statistics
!--------------------------------------------------------------

        call conmima(nlvddi,value,vmin,vmax)

        !write(6,*) 'min/max: ',vmin,vmax

!--------------------------------------------------------------
! end of routine
!--------------------------------------------------------------

	end



!*******************************************************************	

	subroutine ts_next_record_1(it,iunit,nlvddi,nkn,nlv,value)

	use intp_fem_file
        use conz_util

	implicit none

	include 'param.h'

	integer it
	integer iunit(3)
	integer nlvddi
	integer nkn
	integer nlv
	double precision value(nlvddi,nkn)

	integer id,ldim,ndim,ivar
        double precision vmin,vmax
	double precision dtime
	character*80 string

!--------------------------------------------------------------
! initialize
!--------------------------------------------------------------

	id = iunit(1)

!--------------------------------------------------------------
! read new data
!--------------------------------------------------------------

	dtime = it
	ivar = 1
	ndim = nkn
	ldim = nlvddi

	!write(6,*)'reading T/S values', it,dtime

	call iff_read_and_interpolate(id,dtime)
	call iff_time_interpolate(id,dtime,ivar,ndim,ldim,value)

!--------------------------------------------------------------
! some statistics
!--------------------------------------------------------------

        call conmima(nlvddi,value,vmin,vmax)

        !write(6,*) 'min/max: ',vmin,vmax

!--------------------------------------------------------------
! end of routine
!--------------------------------------------------------------

	end

!*******************************************************************	

	subroutine ts_file_close_1(iunit)

! closes T/S file

	use intp_fem_file

	implicit none

	integer iunit(3)

	integer id

	id = iunit(1)
	call iff_forget_file(id)

	end

!*******************************************************************	
!*******************************************************************	
!*******************************************************************	
!****** old routines ***********************************************	
!*******************************************************************	
!*******************************************************************	
!*******************************************************************	


!*******************************************************************	
!*******************************************************************	
!*******************************************************************	

!-------------------------------------------------------------------
        end module tsfile_admin
!-------------------------------------------------------------------
