
	module shyelab_unit

	integer, save :: iunit = 100

	end module shyelab_unit

!***************************************************************

	subroutine get_new_unit(iu)

	use shyelab_unit

	implicit none

	integer iu
	logical bopen

	do
	  iunit = iunit + 1
	  inquire(unit=iunit,opened=bopen)
	  if( .not. bopen ) exit
	end do

	iu = iunit

	end

!***************************************************************
!***************************************************************
!***************************************************************

	subroutine make_iunit_name(short,modi,dim,j,iu)

	implicit none

	character*(*) short,modi,dim
	integer j
	integer iu

	logical bopen
	character*5 numb
	character*80 dimen
	character*80 name

        write(numb,'(i5)') j
        numb = adjustl(numb)

	dimen = '.' // trim(dim) // '.'
	name = trim(short)//trim(modi)//trim(dimen)//numb

	call get_new_unit(iu)

	inquire(file=name,opened=bopen)
	if( bopen ) goto 99
	inquire(unit=iu,opened=bopen)
	if( bopen ) goto 98

        !write(6,*) 'opening file : ',iu,'  ',trim(name)
        open(iu,file=name,form='formatted',status='unknown')

	return
   99	continue
	write(6,*) 'file already open: ',trim(name)
	stop 'error stop make_iunit_name: internal error (1)'
   98	continue
	write(6,*) 'unit already open: ',iu
	stop 'error stop make_iunit_name: internal error (2)'
	end

!***************************************************************

