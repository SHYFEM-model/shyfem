
	program feminf

c shows content of fem boundary/initial files

	implicit none

        include 'param.h'

	logical bformat,bread
	integer np,nunit
	integer it,nvers,lmax,nvar,ntype
	integer nrec,i
	integer ierr
	integer year0
	real hlv(nlvdim)
	character*10 type
	character*60 string,file
	character*30 line

	integer ilhkv(nkndim)
	real hd(nkndim)
	real data(nlvdim,nkndim)
	integer ivalue,mima

	integer iapini

	ivalue = 1		!value to be written
	bread = ivalue.gt. 0

	year0=0
	year0=2007

	if( year0 .gt. 0 ) then
	  call dtsyear(year0)
	end if

        if(iapini(2,nkndim,neldim,0).eq.0) then
                stop 'error stop : iapini'
        end if

	type = ' '
        call def_make(type,file)

	np = 0
        call fem_file_read_open(file,np,nunit,bformat)
	if( nunit .eq. 0 ) goto 99

	write(6,*) '+++++++++++++== ',bformat,np
	if( .false. ) then
        call fem_file_get_params(bformat,nunit,it
     +                          ,nvers,np,lmax,nvar,ntype,ierr)
	if( ierr .ne. 0 ) goto 98

        write(6,*)
        write(6,*) 'bformat = ',bformat
        write(6,*) ' nvers  = ',nvers
        write(6,*) '    np  = ',np,   ' ntype = ',ntype
        write(6,*) '  lmax  = ',lmax, '  nvar = ',nvar
	end if

        write(6,*)
        write(6,*) 'bformat = ',bformat
	if( year0 .gt. 0 ) write(6,*) 'reference year: ',year0
        write(6,*)

	nrec = 0
    1	continue
          call fem_file_read_header(bformat,nunit,it
     +                  ,nvers,np,lmax,nvar,ntype,nlvdim,hlv,ierr)
	  if( ierr .ne. 0 ) goto 2

	  if( bread ) then
	    if( np .gt. nkndim ) goto 89
	    if( lmax .gt. nlvdim ) goto 89
	  end if

          do i=1,nvar
	   if( bread ) then
            call fem_file_read_data(bformat,nunit
     +                          ,nvers,np,lmax
     +                          ,ilhkv,hd
     +                          ,string,nlvdim,data)
	    if( ivalue .eq. i ) then
	      call minmax(it,nlvdim,np,ilhkv,data)
	      call write_value(it,nlvdim,np,data)
	      call write_node(it,nlvdim,np,data)
	    end if
	   else
            call fem_file_skip_data(bformat,nunit
     +                          ,nvers,np,lmax
     +                          ,string)
	   end if
	   if( nrec .eq. 0 .and. i .eq. 1 ) then
	     write(6,*) nvers,np,lmax,nvar
	     write(6,*)
	   end if

	   if( nrec .eq. 0 ) write(6,*) i,string
	  end do
	  if( nrec .eq. 0 ) write(6,*)

	  nrec = nrec + 1
	  call make_time(it,year0,line)
	  write(6,'(i6,i12,a,a)') nrec,it,'  ',line

	  goto 1
    2	continue
	if( ierr .gt. 0 ) goto 97

	stop
   89	continue
	write(6,*) 'nlvdim,lmax: ',nlvdim,lmax
	write(6,*) 'nkndim,np: ',nkndim,np
	stop 'error stop feminf: reading header'
   97	continue
	stop 'error stop feminf: reading header'
   98	continue
	stop 'error stop feminf: getting params'
   99	continue
	stop 'error stop feminf: opening file'
	end

c*****************************************************************

	subroutine make_time(it,year0,line)

	implicit none

	integer it
	integer year0
	character*(*) line

	integer year,month,day,hour,min,sec

	line = ' '
	if( year0 .le. 0 ) return

	call dts2dt(it,year,month,day,hour,min,sec)
	call dtsform(year,month,day,hour,min,sec,line)

	end

c*****************************************************************

	subroutine write_node(it,nlvdim,np,data)

	implicit none

	integer it
	integer nlvdim,np
	real data(nlvdim,1)

	integer nnodes
	parameter(nnodes=4)
	integer nodes(nnodes)
	save nodes
	data nodes /9442,10770,13210,14219/

	integer n,i

	n = nnodes
	write(90,'(i10,10i6)') it,(ifix(data(1,nodes(i))),i=1,n)

	end

c*****************************************************************

	subroutine write_value(it,nlvdim,np,data)

	implicit none

	integer it
	integer nlvdim,np
	real data(nlvdim,1)

	integer n,nskip,i

	n = 10
	nskip = np/n

	!write(89,*) np,n,nskip,n*nskip
	write(89,'(i10,10i6)') it,(ifix(data(1,i*nskip)),i=1,n)

	end

c*****************************************************************

	subroutine minmax(it,nlvdim,np,ilhkv,data)

	implicit none

	integer it
	integer nlvdim,np
	integer ilhkv(1)
	real data(nlvdim,1)

	integer k,l,lmax
	real vmin,vmax,v

	vmin = data(1,1)
	vmax = data(1,1)

	do k=1,np
	  lmax = ilhkv(k)
	  do l=1,lmax
	    v = data(l,k)
	    vmax = max(vmax,v)
	    vmin = min(vmin,v)
	  end do
	end do

	write(86,*) 'min/max: ',it,vmin,vmax

	end

c*****************************************************************
