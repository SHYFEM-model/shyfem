
c info on restart file

c******************************************************************

	program rstinf

	implicit none

	integer iunit,it,nvers,nrec,nkn,nel,nlv,iflag,ierr
	integer nread
	double precision atime
	double precision atime_anf
	double precision atime_end
	character*20 aline
	character*60 file
	character*52 title1
	character*72 title2

!-------------------------------------------------------------------
! create title strings
!-------------------------------------------------------------------

	title1 = 'version nrec       nkn       nel       nlv     iflag'
c                 12345678901234567890123456789012345678901234567890123
	title2 = '   irec                         atime     date'

!-------------------------------------------------------------------
! initialize and open file
!-------------------------------------------------------------------

	nread = 0
	iunit = 1
	file = ' '

        call rst_init(file)

	open(iunit,file=file,status='old',form='unformatted')

!-------------------------------------------------------------------
! loop on records
!-------------------------------------------------------------------

	do

	call rst_skip_record(iunit,atime,nvers,nrec
     +					,nkn,nel,nlv,iflag,ierr)
	if( ierr .ne. 0 ) exit

	if( nread == 0 ) then
	  write(6,1000) title1
	  write(6,1010) nvers,nrec,nkn,nel,nlv,iflag
	  write(6,*)
	  write(6,1001) title2
	  atime_anf = atime
	end if

	nread = nread + 1
	call dts_format_abs_time(atime,aline)
	write(6,1011) nread,atime,aline
	atime_end = atime

	end do

	if( ierr > 0 ) stop 'error stop rstinf: error reading record'

!-------------------------------------------------------------------
! final message
!-------------------------------------------------------------------

	write(6,1001) title2
	write(6,*)
	write(6,1000) title1
	write(6,1010) nvers,nrec,nkn,nel,nlv,iflag
	write(6,*)
	write(6,*) 'Number of records read: ',nread
	call dts_format_abs_time(atime_anf,aline)
	write(6,*) 'Initial time in file:   ',atime_anf,aline
	call dts_format_abs_time(atime_end,aline)
	write(6,*) 'Final time in file:     ',atime_end,aline
	write(6,*)
	write(6,*) 'Meaning of iflag:'
	write(6,*) '         1          hydro'
	write(6,*) '        10          depth'
	write(6,*) '       100          ibarcl (T/S/rho)'
	write(6,*) '      1000          iconz (cnv/conzv)'
	write(6,*) '     10000          iwvert (wlnv)'
	write(6,*) '    100000          ieco (ecological variables)'

!-------------------------------------------------------------------
! end of routine
!-------------------------------------------------------------------

	stop
 1000	format(a52)
 1001	format(a72)
 1010	format(i7,i5,4i10)
 1011	format(i7,f30.2,5x,a20)
	end

c******************************************************************

        subroutine rst_init(rstfile)

        use clo

        implicit none

        character*(*) rstfile

        call shyfem_copyright('rstinf - info on restart file')

        call clo_init('rstinf','rstfile','1.2')

        call clo_add_info('returns info on records of restart file')

        call clo_parse_options

        call clo_check_files(1)
        call clo_get_file(1,rstfile)

        end

c******************************************************************

