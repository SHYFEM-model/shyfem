
	program test_unformatted

	implicit none

	integer iunit,nl
	integer fields(3)

	open(1,file='unform.dat',form='unformatted',status='unknown')
	write(1) 1
	close(1)

	open(1,file='stream.dat',form='unformatted',status='unknown'
     +			,access='stream')
	write(1) 1
	close(1)

	iunit = 1
	open(iunit,file='unform.dat',form='unformatted',status='old'
     +			,access='stream')
	call skip_record(iunit,nl)
	close(iunit)

	write(6,*) nl,' bytes contained in record'

	end

!**************************************************************

	subroutine skip_record(iunit,nl)

	implicit none

	integer iunit
	integer nl

	integer nle

	call read_length(iunit,nl)
	call skip_bytes(iunit,nl)
	call read_length(iunit,nle)

	if( nl /= nle ) then
	  write(6,*) 'cannot read record: ',nl,nle
	  stop 'error stop read_record: nl/=nle'
	end if

	end

!**************************************************************

	subroutine read_length(iunit,nl)

	implicit none

	integer iunit
	integer nl

	read(iunit) nl

	end

!**************************************************************

	subroutine skip_bytes(iunit,nb)

	implicit none

	integer iunit
	integer nb

	integer (kind=1), allocatable :: aux(:)

	allocate(aux(nb))

	read(iunit) aux(1:nb)

	end

!**************************************************************

	subroutine read_words(iunit,nb,words)

	implicit none

	integer iunit
	integer nb
	integer words(*)

	integer nw

	nw = nb / 4
	if( 4*nw /= nb ) then
	  write(6,*) 'nb = ',nb,'  nw = ',nw
	  stop 'error stop read_words: number of bytes no multiple of 4'
	end if

	read(iunit) words(1:nw)

	end

!**************************************************************

