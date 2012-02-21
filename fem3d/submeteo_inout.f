
c***************************************************************
c***************************************************************
c***************************************************************

	subroutine write_heat_unformatted(iunit,it,n,s,t,h,c)

	implicit none

	integer iunit,it,n
	real s(n),t(n),h(n),c(n)

	integer i,nvar,id

	nvar = 4
	id = 1002

	write(iunit) it,id,n,nvar
	write(iunit) (s(i),t(i),h(i),c(i),i=1,n)

	end

c***************************************************************

	subroutine read_heat_unformatted(iunit,it,n,s,t,h,c,ierr)

	implicit none

	integer iunit,it,n,ierr
	real s(1),t(1),h(1),c(1)

	integer i,n_read,nvar,nvar_read,id,id_read

	nvar = 4
	id = 1002
	ierr = 0

	read(iunit,iostat=ierr) it,id_read,n_read,nvar_read

	if( ierr .lt. 0 ) return	!EOF
	if( ierr .gt. 0 ) goto 99
	if( nvar .ne. nvar_read ) goto 98
	if( id .ne. id_read ) goto 98
	if( n .ne. 0 .and. n .ne. n_read ) goto 98

	read(iunit) (s(i),t(i),h(i),c(i),i=1,n)

	n = n_read

	return
   98	continue
	write(6,*) 'it: ',it
	write(6,*) 'n,nvar,id (read): ',n_read,nvar_read,id_read
	write(6,*) 'n,nvar,id (expected): ',n,nvar,id
	stop 'errro stop read_heat_unformatted: parameters in header record'
   99	continue
	stop 'errro stop read_heat_unformatted: cannot read header record'
	end

c***************************************************************
c***************************************************************
c***************************************************************

	subroutine write_rain_unformatted(iunit,it,n,r)

	implicit none

	integer iunit,it,n
	real r(n)

	integer i,nvar,id

	nvar = 1
	id = 1003

	write(iunit) it,id,n,nvar
	write(iunit) (r(i),i=1,n)

	end

c***************************************************************

	subroutine read_rain_unformatted(iunit,it,n,r,ierr)

	implicit none

	integer iunit,it,n,ierr
	real r(1)

	integer i,n_read,nvar,nvar_read,id,id_read

	nvar = 1
	id = 1003
	ierr = 0

	read(iunit,iostat=ierr) it,id_read,n_read,nvar_read

	if( ierr .lt. 0 ) return	!EOF
	if( ierr .gt. 0 ) goto 99
	if( nvar .ne. nvar_read ) goto 98
	if( id .ne. id_read ) goto 98
	if( n .ne. 0 .and. n .ne. n_read ) goto 98

	read(iunit) (r(i),i=1,n)

	n = n_read

	return
   98	continue
	write(6,*) 'it: ',it
	write(6,*) 'n,nvar,id (read): ',n_read,nvar_read,id_read
	write(6,*) 'n,nvar,id (expected): ',n,nvar,id
	stop 'errro stop read_rain_unformatted: parameters in header record'
   99	continue
	stop 'errro stop read_rain_unformatted: cannot read header record'
	end

c***************************************************************
c***************************************************************
c***************************************************************

	subroutine write_wind_unformatted(iunit,it,n,wx,wy,p)

	implicit none

	integer iunit,it,n
	real wx(n),wy(n),p(n)

	integer i,nvar,id

	nvar = 3
	id = 1001

	write(iunit) it,id,n,nvar
	write(iunit) (wx(i),wy(i),p(i),i=1,n)

	end

c***************************************************************

	subroutine read_wind_unformatted(iunit,it,n,wx,wy,p,ierr)

	implicit none

	integer iunit,it,n,ierr
	real wx(1),wy(1),p(1)

        real pstd
        parameter ( pstd = 1013.25 )

	integer i,n_read,nvar,nvar_read,mode,id,id_read

	nvar = 3
	id = 1001
	ierr = 0
	mode = 0

	read(iunit,iostat=ierr) it,id_read,n_read,nvar_read
	if( ierr .lt. 0 ) return	!EOF

	if( ierr .gt. 0 ) then
	  backspace(iunit)
	  read(iunit,iostat=ierr) it,n_read		!old read
	  if( ierr .gt. 0 ) goto 99
	  mode = 1
	  id_read = id
	  nvar_read = nvar
	  if( n_read .lt. 0 ) then
	    mode = 2
	    n_read = -n_read
	  end if
	end if

	if( nvar .ne. nvar_read ) goto 98
	if( id .ne. id_read ) goto 98
	if( n .ne. 0 .and. n .ne. n_read ) goto 98

	if( mode .eq. 0 ) then		!new read
	  read(iunit) (wx(i),wy(i),p(i),i=1,n)
	else if( mode .eq. 1 ) then	!old read - no pressure
	  read(iunit) (wx(i),wy(i),i=1,n)
	  do i=1,n
	    p(i) = 100.*pstd
	  end do
	else if( mode .eq. 2 ) then	!old read - with pressure
	  read(iunit) (wx(i),wy(i),i=1,n),(p(i),i=1,n)
	end if

	n = n_read

	return
   98	continue
	write(6,*) 'it: ',it
	write(6,*) 'n,nvar,id (read): ',n_read,nvar_read,id_read
	write(6,*) 'n,nvar,id (expected): ',n,nvar,id
	stop 'errro stop read_wind_unformatted: parameters in header record'
   99	continue
	stop 'errro stop read_wind_unformatted: cannot read header record'
	end

c***************************************************************
c***************************************************************
c***************************************************************

	subroutine test_meteo_unformatted(file,it,id,n,nvar)

	implicit none

	character*(*) file
	integer it,id,n,nvar
	integer ierr

	integer iunit

	open(iunit,file=file,form='unformatted',status='old')

	it = 0
	id = 0
	n = 0
	nvar = 0

	read(iunit,iostat=ierr) it,id,n,nvar
	if( ierr .lt. 0 ) goto 1	!EOF

	if( ierr .gt. 0 ) then
	  backspace(iunit)
	  read(iunit,iostat=ierr) it,n		!old wind read
	  if( ierr .ne. 0 ) goto 1
	  id = 1001
	  nvar = 2
	  if( n .lt. 0 ) then
	    n = -n
	    nvar = 3
	  end if
	end if

    1	continue
	close(iunit)

	end

c***************************************************************
c***************************************************************
c***************************************************************

	subroutine test_wind

	implicit none
	
	integer ndim
	parameter (ndim=10000)

	character*80 file
	integer iunit,irec,ierr,it,id
	integer n,nvar
	real wx(ndim)
	real wy(ndim)
	real p(ndim)

	irec = 0
	iunit = 1
	file = 'meteodat.win'
	file = 'wind.win'

	call test_meteo_unformatted(file,it,id,n,nvar)
	write(6,*) 'file info: ',it,id,n,nvar
	if( id .eq. 0 ) stop

	open(iunit,file=file,form='unformatted',status='old')

	do while(.true.)
	  n = 0
	  call read_wind_unformatted(iunit,it,n,wx,wy,p,ierr)
	  if( ierr .lt. 0 ) goto 1
	  irec = irec + 1
	  write(6,*) irec,it,n
	end do

    1	continue

	end

c***************************************************************
c	program test_wind_main
c	call test_wind
c	end
c***************************************************************

