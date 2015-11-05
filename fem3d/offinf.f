
c**********************************************************

	program offinf

c shows content of offline data file

	implicit none

	integer it,nkn,nel,nrec,iu,i,type
	integer ios
	character*60 name
	double precision buffer(5)

        call off_init(name)

	nrec = 0
	iu = 1

	open(iu,file=name,status='old',form='unformatted',iostat=ios)
	if( ios /= 0 ) stop 'error stop offinf: opening file'

	do
          read(iu,iostat=ios) it,nkn,nel,type
	  if( ios /= 0 ) exit
	  nrec = nrec + 1
	  write(6,*) nrec,it,nkn,nel,type
	  if( type /= 3 ) stop 'error stop offinf: type /= 3'

	  do i=1,9
	    read(iu) buffer
	  end do

	  !write(6,'(5g14.6)') buffer
	end do

	if( ios > 0 ) stop 'error stop offinf: reading file'

	close(iu)

	end

c**********************************************************

        subroutine off_init(offfile)

        use clo

        implicit none

        character*(*) offfile

        call shyfem_copyright('offinf - info on offline file')

        call clo_init('offinf','offfile','1.0')

	call clo_add_info('returns info on records of offline file')

        call clo_parse_options

        call clo_check_files(1)
        call clo_get_file(1,offfile)

        end

c**********************************************************

