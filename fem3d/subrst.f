c
c $Id: subrst.f,v 1.11 2010-03-11 15:36:39 georg Exp $
c
c restart routines
c
c contents :
c
c subroutine inirst		reads and initializes values from restart
c subroutine admrst		administers writing of restart file

c subroutine wrrst(it,iunit)	writes one record of restart data
c subroutine rdrst(itrst,iunit)	reads one record of restart data
c
c subroutine skip_rst(iunit,atime,it,nvers,nrec,nkn,nel,nlv,iflag,ierr)
c				returns info on record in restart file
c
c revision log :
c
c 02.09.2004    ggu     started with routine rdrst()
c 18.10.2006    ggu     included hm3v in restart file -> nvers = 4
c 13.06.2008    ggu     new version 5 -> S/T/rho
c 09.01.2009    ggu     bugfix in inirst - file opened with status=new
c 23.03.2009    ggu     bugfix in admrst - itmrst=itanf if itmrst<itanf
c 07.05.2009    ggu     new parameter ityrst
c 29.05.2009    ggu     use closest record for restart (if ityrst=2)
c 13.11.2009    ggu     keep track of restart: /rstrst/ and has_restart()
c 27.11.2009    ggu     deal with empty file, rdrst() restructured
c 19.01.2010    ggu     initialize also conz, has_restart() is function
c 11.03.2010    ggu     write also vertical velocity
c 10.02.2012    ggu     write only last record, restart from last record
c 16.02.2012    aac     write also ecological varibles
c 28.08.2012    aac     bug fix for restart time = -1 (rdrst_record)
c 11.12.2014    ccf     bug fix for atime
c
c
c notes :
c
c  |--------|---------|------------|--------------|--------------|
c  | ityrst | no file | file empty ! itrec/=itrst ! itrec==itrst !
c  |--------|---------|------------|--------------|--------------|
c  |   0    |  error  |   error    |     error    |    warm      |
c  |--------|---------|------------|--------------|--------------|
c  |   1    |  cold   |   error    |     error    |    warm      |
c  |--------|---------|------------|--------------|--------------|
c  |   2    |  cold   |   cold     |     warm     |    warm      |
c  |--------|---------|------------|--------------|--------------|
c
c*********************************************************************

        subroutine inirst

c reads and initializes values from restart

        implicit none

	integer iokrst,nvers,ibarcl,iconz,iwvert,ieco
	common /rstrst/ iokrst,nvers,ibarcl,iconz,iwvert,ieco
	save /rstrst/

        integer itrst,iunit,ierr,ityrst,it
	integer itanf,itend
	integer date,time
	double precision dit
	double precision atime,atime0,atrst
        character*80 name
        character*20 line

        real getpar
	double precision dgetpar
        integer ifileo

c-----------------------------------------------------------------
c get parameters
c-----------------------------------------------------------------

	call init_restart	!sets params to 0

        !itrst = nint(dgetpar('itrst'))
	call convert_date('itrst',itrst)
        ityrst = nint(dgetpar('ityrst'))
        call getfnm('restrt',name)
        if(name.eq.' ') return

        date = nint(dgetpar('date'))
        time = nint(dgetpar('time'))

	call dts_to_abs_time(date,time,atime0)
	atrst = atime0 + itrst
	if( itrst == -1 ) atrst = -1.

c-----------------------------------------------------------------
c name of restart file given -> open and read
c-----------------------------------------------------------------

        write(6,*) '---------------------------------------------'
        write(6,*) '... performing restart from file'
        write(6,*) name
        write(6,*) '---------------------------------------------'

        iunit = ifileo(1,name,'unformatted','old')
        if( iunit .le. 0 ) then
          if( ityrst .le. 0 ) goto 98
          write(6,*) '*** Cannot find restart file ...'
          write(6,*) '*** Continuing with cold start...'
          return
        end if

	atime = atrst
        call rdrst(atime,iunit,ierr)
        if( ierr .gt. 0 ) then
          if( ityrst .le. 1 ) goto 97
          if( ierr .eq. 95 ) then	!file is empty
            write(6,*) '*** No data in restart file ...'
            write(6,*) '*** Continuing with cold start...'
            return
          end if
          write(6,*) '*** Another time record is used for restart'
	  call dts_format_abs_time(atrst,line)
          write(6,*) '*** Looking for time = ',itrst,line
	  call dts_format_abs_time(atime,line)
	  it = nint(atime-atime0)
          write(6,*) '*** Finding time = ',it,line
          write(6,*) '*** Continuing with hot start...'
        end if

	close(iunit)

	if( itrst .eq. -1 ) then	!reset initial time
	  it = nint(atime-atime0)
	  call dts_format_abs_time(atime,line)
	  write(6,*) 'setting new initial time: ',it,line
	  dit = it
	  call dputpar('itanf',dit)
	  write(6,*) 'new itanf: ',it,line
	end if

	itanf = nint(dgetpar('itanf'))
	itend = nint(dgetpar('itend'))

        write(6,*) '---------------------------------------------'
        write(6,*) 'A restart has been performed'
	call dtsgf(itrst,line)
        write(6,*) ' requested restart time = ',itrst,line
	it = nint(atime-atime0)
	call dtsgf(it,line)
        write(6,*) ' used restart time = ',it,line
        write(6,*) ' nvers,ibarcl,iconz,ieco ',nvers,ibarcl,iconz,ieco
        write(6,*) ' iwvert ',iwvert
	call dtsgf(itanf,line)
        write(6,*) ' itanf = ',itanf,line
	call dtsgf(itend,line)
        write(6,*) ' itend = ',itend,line
        write(6,*) '---------------------------------------------'

	iokrst = 1

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

        return
   97   continue
        write(6,*) 'no record found for time = ',itrst
        stop 'error stop inirst: Cannot find time record'
   98   continue
        write(6,*) 'no such file : ',name
        stop 'error stop inirst: Cannot read restart file'
        end

c*******************************************************************

	subroutine init_restart

	implicit none

	integer iokrst,nvers,ibarcl,iconz,iwvert,ieco
	common /rstrst/ iokrst,nvers,ibarcl,iconz,iwvert,ieco
	save /rstrst/

	iokrst = 0
	nvers = 0
	ibarcl = 0
	iconz = 0
	iwvert = 0
	ieco = 0

	end

c*******************************************************************

	function has_restart(icode)

c gives indication if restart has been performed
c
c icode indicates what information is requested
c
c icode = 0	general restart
c icode = 1	basic restart (hydro values)
c icode = 2	depth values
c icode = 3	t/s/rho values
c icode = 4	conz values
c icode = 5	vertical velocity
c icode = 6	ecological variables

	implicit none

	logical has_restart
	integer icode

	integer iokrst,nvers,ibarcl,iconz,iwvert,ieco
	common /rstrst/ iokrst,nvers,ibarcl,iconz,iwvert,ieco
	save /rstrst/

	if( iokrst .le. 0 ) then
	  has_restart = .false.
	else if( icode .eq. 0 .or. icode .eq. 1 ) then
	  has_restart = .true.
	else if( icode .eq. 2 ) then
	  has_restart = nvers .ge. 4
	else if( icode .eq. 3 ) then
	  has_restart = ibarcl .gt. 0
	else if( icode .eq. 4 ) then
	  has_restart = iconz .gt. 0
	else if( icode .eq. 5 ) then
	  has_restart = iwvert .gt. 0
	else if( icode .eq. 6 ) then
	  has_restart = ieco .gt. 0
	else
	  has_restart = .false.
	end if

	end

c*******************************************************************

        subroutine admrst

c administers writing of restart file

        implicit none

	include 'femtime.h'

	integer ierr
        integer idtrst,itmrst,iunit

        real getpar
        double precision dgetpar
        integer ifemop
	!integer fsync
	logical has_output,next_output

	logical bonce
	save bonce
	integer ia_out(4)
	save ia_out
        integer icall
        save icall
        data icall / 0 /

        if( icall .le. -1 ) return

c-----------------------------------------------------
c initializing
c-----------------------------------------------------

        if( icall .eq. 0 ) then

          call convert_date('itmrst',itmrst)
          call convert_time('idtrst',idtrst)

	  if( idtrst .lt. 0 ) then	!only last record saved
	    if( idtrst .eq. -1 ) then	!only at the end of the simulation
	      itmrst = itanf
	      idtrst = -(itend-itanf)
	    end if
	    bonce = .true.
	    idtrst = -idtrst
	  else
	    bonce = .false.
	  end if

          icall = -1
          call set_output_frequency(itmrst,idtrst,ia_out)
          if( .not. has_output(ia_out) ) return
          icall = 1

	  if( .not. bonce ) then
            iunit = ifemop('.rst','unformatted','new')
            if( iunit .le. 0 ) goto 98
	    ia_out(4) = iunit
	  end if

        end if

c-----------------------------------------------------
c normal call and writing
c-----------------------------------------------------

        if( .not. next_output(ia_out) ) return

	if( bonce ) then
          iunit = ifemop('.rst','unformatted','new')
          if( iunit .le. 0 ) goto 98
          call wrrst(it,iunit)
	  close(iunit)
	else
	  iunit = ia_out(4)
          call wrrst(it,iunit)
	  flush(iunit)
          !ierr = fsync(fnum(iunit))
	end if

c-----------------------------------------------------
c end of routine
c-----------------------------------------------------

        return
   98   continue
        stop 'error stop admrst: Cannot open rst file'
        end

c*******************************************************************

        subroutine wrrst(it,iunit)

c writes one record of restart data

	use mod_conz
	use mod_geom_dynamic
	use mod_ts
	use mod_hydro_vel
	use mod_hydro
	use levels, only : nlvdi,nlv
	use basin

        implicit none

        integer it
        integer iunit

        include 'param.h'




        integer ii,l,ie,k,i
	integer ibarcl,iconz,ieco
        integer nvers
	integer date,time

	real getpar
        double precision dgetpar

        nvers = 9

	ibarcl = nint(getpar('ibarcl'))
	iconz = nint(getpar('iconz'))
	ieco = nint(getpar('ibfm'))
        date = nint(dgetpar('date'))
        time = nint(dgetpar('time'))

        write(iunit) it,nvers,1
        write(iunit) date,time
        write(iunit) nkn,nel,nlv

        write(iunit) (iwegv(ie),ie=1,nel)
        write(iunit) (znv(k),k=1,nkn)
        write(iunit) ((zenv(ii,ie),ii=1,3),ie=1,nel)
        write(iunit) ((utlnv(l,ie),l=1,nlv),ie=1,nel)
        write(iunit) ((vtlnv(l,ie),l=1,nlv),ie=1,nel)

        write(iunit) ((hm3v(ii,ie),ii=1,3),ie=1,nel)

        write(iunit) ibarcl
	if( ibarcl .gt. 0 ) then
          write(iunit) ((saltv(l,k),l=1,nlv),k=1,nkn)
          write(iunit) ((tempv(l,k),l=1,nlv),k=1,nkn)
          write(iunit) ((rhov(l,k),l=1,nlv),k=1,nkn)
	end if

        write(iunit) iconz
	if( iconz .eq. 1 ) then
          write(iunit) ((cnv(l,k),l=1,nlv),k=1,nkn)
	else if( iconz .gt. 1 ) then
          write(iunit) (((conzv(l,k,i),l=1,nlv),k=1,nkn),i=1,iconz)
	end if
	
        write(iunit) nlv-1
	if( nlv .gt. 1 ) then
          write(iunit) ((wlnv(l,k),l=0,nlv),k=1,nkn)
	end if

	write(iunit) ieco
	if( ieco .gt. 0 ) then
	  call write_restart_eco(iunit)
        end if

        end

c*******************************************************************

        subroutine rdrst(atrst,iunit,ierr)

c reads restart file until it finds itrst

        implicit none

        double precision atrst
        integer iunit
        integer ierr            !error code - different from 0 if error

        integer ii,l,ie,k
        integer it
        integer irec
        logical bloop,blast,bnext
        double precision atime,alast
	character*20 line

        irec = 0
	ierr = 0
        bloop = .true.
	blast = atrst .eq. -1		! take last record

        do while( bloop )
          call rdrst_record(atime,it,iunit,ierr)
          if( ierr .eq. 0 ) irec = irec + 1
	  bnext = atime .lt. atrst .or. blast	!look for more records
          bloop = ierr .eq. 0 .and. bnext
	  alast = atime
        end do

	if( irec .gt. 0 ) then
	  if( blast ) then
	    ierr = 0
	    atrst = alast
	    return
          else if( atime .eq. atrst ) then
	    return
	  end if
	end if

        if( ierr .ne. 0 ) then          !EOF
          if( irec .eq. 0 ) then        !no data found
                goto 95
          else                          !last record in file has smaller time
                goto 97
          end if
        else                            !read past desired time
                goto 97
        end if

        return
   95   continue
        write(6,*) 'reading restart file... '
        write(6,*) 'no records in file (file is empty)'
        ierr = 95
        return
   97   continue
        write(6,*) 'reading restart file... '
	call dts_format_abs_time(atime,line)
        write(6,*) 'last record read at time = ',atime,line
	call dts_format_abs_time(atrst,line)
        write(6,*) 'no record found for time = ',atrst,line
        ierr = 97
        return
        end

c*******************************************************************

	subroutine skip_rst(iunit,atime,it,nvers,nrec
     +				,nkn,nel,nlv,iflag,ierr)

c returns info on record in restart file and skips data records

	implicit none

	integer iunit,it,nvers,nrec,nkn,nel,nlv,iflag,ierr
	double precision atime
	integer ibarcl,iconz,iwvert,ieco
	integer date,time

	iflag = 0
	date = 0
	time = 0

	read(iunit,end=2,err=3) it,nvers,nrec
        if( nvers > 8 ) read(iunit) date,time
        read(iunit) nkn,nel,nlv

	atime = 0.
	if( date > 0 ) call dts_to_abs_time(date,time,atime)
	atime = atime + it

	read(iunit)
	read(iunit)
	read(iunit)
	read(iunit)
	read(iunit)
	if( nvers .ge. 4 ) then
	  iflag = iflag + 1
	  read(iunit)
	end if

	if( nvers .ge. 5 ) then
	  read(iunit) ibarcl
	  if( ibarcl .gt. 0 ) then
	    iflag = iflag + 10
	    read(iunit)
	    read(iunit)
	    read(iunit)
	  end if
	end if

	if( nvers .ge. 6 ) then
	  read(iunit) iconz
	  if( iconz .gt. 0 ) then
	    iflag = iflag + 100
	    read(iunit)
	  end if
	end if

	if( nvers .ge. 7 ) then
	  read(iunit) iwvert
	  if( iwvert .gt. 0 ) then
	    iflag = iflag + 1000
	    read(iunit)
	  end if
	end if

	if( nvers .ge. 8 ) then
          read(iunit) ieco
          if( ieco .gt. 0 ) then
	    iflag = iflag + 10000
	    call skip_restart_eco(iunit)
          end if
        end if

	ierr = 0
	return

    2	continue
	ierr = -1
	return
    3	continue
	write(6,*) 'skip_rst: error in reading restart file'
	ierr = 1
	return
	end

c*******************************************************************

        subroutine rdrst_record(atime,it,iunit,ierr)

c reads one record of restart data

	use mod_conz
	use mod_geom_dynamic
	use mod_ts
	use mod_hydro_vel
	use mod_hydro
	use levels, only : nlvdi,nlv
	use basin

        implicit none

        double precision atime
	integer it
        integer iunit
        integer ierr            !error code - different from 0 if error

        include 'param.h'




	integer iokrst,nvers,ibarcl,iconz,iwvert,ieco
	common /rstrst/ iokrst,nvers,ibarcl,iconz,iwvert,ieco
	save /rstrst/

        integer ii,l,ie,k,i
        integer nversaux,nrec
        integer nknaux,nelaux,nlvaux
	integer date,time

        read(iunit,end=97) it,nvers,nrec
        if( nvers .lt. 3 ) goto 98

        ierr = 0
	ibarcl = 0
	iconz = 0
	iwvert = 0
	ieco = 0
	date = 0
	time = 0

        if( nvers > 8 ) read(iunit) date,time

	atime = 0.
	if( date > 0 ) call dts_to_abs_time(date,time,atime)
	atime = atime + it

          read(iunit) nknaux,nelaux,nlvaux
          if( nknaux .ne. nkn ) goto 99
          if( nelaux .ne. nel ) goto 99
          if( nlvaux .ne. nlv ) goto 99

          read(iunit) (iwegv(ie),ie=1,nel)
          read(iunit) (znv(k),k=1,nkn)
          read(iunit) ((zenv(ii,ie),ii=1,3),ie=1,nel)
          read(iunit) ((utlnv(l,ie),l=1,nlv),ie=1,nel)
          read(iunit) ((vtlnv(l,ie),l=1,nlv),ie=1,nel)

          if( nvers .ge. 4 ) then
            read(iunit) ((hm3v(ii,ie),ii=1,3),ie=1,nel)
          end if

          if( nvers .ge. 5 ) then
            read(iunit) ibarcl
            if( ibarcl .gt. 0 ) then
              read(iunit) ((saltv(l,k),l=1,nlv),k=1,nkn)
              read(iunit) ((tempv(l,k),l=1,nlv),k=1,nkn)
              read(iunit) ((rhov(l,k),l=1,nlv),k=1,nkn)
            end if
          end if

          if( nvers .ge. 6 ) then
            read(iunit) iconz
	    if( iconz .eq. 1 ) then
              read(iunit) ((cnv(l,k),l=1,nlv),k=1,nkn)
	    else if( iconz .gt. 1 ) then
              read(iunit) (((conzv(l,k,i),l=1,nlv),k=1,nkn),i=1,iconz)
	    end if
	  end if

	  if( nvers .ge. 7 ) then
	    read(iunit) iwvert
	    if( iwvert .gt. 0 ) then
              read(iunit) ((wlnv(l,k),l=0,nlv),k=1,nkn)
	    end if
	  end if

	  if( nvers .ge. 8 ) then
	    read(iunit) ieco
            if( ieco .gt. 1 ) then
	      call read_restart_eco(iunit)
	    end if
          end if

        return
   97   continue
        ierr = -1
        return
   98   continue
        write(6,*) 'reading restart file...'
        write(6,*) 'nvers: ',nvers
        stop 'error stop rdrst: cannot read this version'
   99   continue
        write(6,*) 'reading restart file...'
        write(6,*) nkn,nel,nlv
        write(6,*) nknaux,nelaux,nlvaux
        stop 'error stop rdrst: incompatible parameters'
        end

c*******************************************************************

