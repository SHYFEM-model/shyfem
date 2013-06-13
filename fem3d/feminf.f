
	program feminf

c shows content of fem boundary/initial files

	implicit none

        include 'param.h'

	logical bformat
	integer np,nunit
	integer it,nvers,lmax,nvar,ntype
	integer nrec,i
	integer ierr
	real hlv(nlvdim)
	character*10 type
	character*60 string,file

	integer iapini

        if(iapini(2,nkndim,neldim,0).eq.0) then
                stop 'error stop : iapini'
        end if

	type = ' '
        call def_make(type,file)

	np = 0
        call fem_file_read_open(file,np,nunit,bformat)
	if( nunit .eq. 0 ) goto 99

        call fem_file_get_params(bformat,nunit,it
     +                          ,nvers,np,lmax,nvar,ntype,ierr)
	if( ierr .ne. 0 ) goto 98

        write(6,*)
        write(6,*) 'bformat = ',bformat
        write(6,*) ' nvers  = ',nvers
        write(6,*) '    np  = ',np,   ' ntype = ',ntype
        write(6,*) '  lmax  = ',lmax, '  nvar = ',nvar
        write(6,*)

	nrec = 0
    1	continue
          call fem_file_read_header(bformat,nunit,it
     +                  ,nvers,np,lmax,nvar,ntype,nlvdim,hlv,ierr)
	  if( ierr .ne. 0 ) goto 2

          do i=1,nvar
            call fem_file_skip_data(bformat,nunit
     +                          ,nvers,np,lmax
     +                          ,string)
	    if( nrec .eq. 0 ) write(6,*) i,string
	  end do
	  if( nrec .eq. 0 ) write(6,*)

	  nrec = nrec + 1
	  write(6,*) nrec,it

	  goto 1
    2	continue
	if( ierr .gt. 0 ) goto 97

	stop
   97	continue
	stop 'error stop feminf: reading header'
   98	continue
	stop 'error stop feminf: getting params'
   99	continue
	stop 'error stop feminf: opening file'
	end

