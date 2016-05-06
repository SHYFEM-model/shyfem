
c utilities for NOS files
c
c revision log :
c
c 29.04.2010    ggu     new file from nosaver
c 07.05.2010    ggu     new routines qopen_nos_file(), checks ev initialization
c 15.12.2010    ggu     volume computation also for sigma layers
c 10.11.2011    ggu     new routines for hybrid levels, init_volume() changed
c 02.12.2011    ggu     bug fix for call to get_sigma_info() (missing argument)
c 10.02.2012    ggu     new routines to get initial/final time of records
c 25.01.2013    ggu     new routines nos_get_vars()
c 05.09.2013    ggu     new call to get_layer_thickness()
c 20.01.2014    ggu     new helper routines
c 23.09.2015    ggu     close files in nos_get_it_start() nos_get_it_end()
c
c***************************************************************

	subroutine make_vert_aver(nlvddi,nkn,ilhkv,cv3,vol3,cv2)

	implicit none

	integer nlvddi
	integer nkn
	integer ilhkv(nkn)
	real cv3(nlvddi,nkn)
	real vol3(nlvddi,nkn)
	real cv2(nkn)

	integer k,l,lmax
	double precision c,v
	double precision cctot,vvtot
	integer :: ks = 3096
	logical bdebug

	do k=1,nkn
	  cctot = 0.
	  vvtot = 0.
	  lmax = ilhkv(k)
	  bdebug = k == ks
	  do l=1,lmax
	    c = cv3(l,k)
	    v = vol3(l,k)
	    cctot = cctot + c*v
	    vvtot = vvtot + v
	    if( bdebug ) write(42,*) l,v,c
	  end do
	  if( bdebug ) write(42,*) lmax,vvtot,cctot,real(cctot/vvtot)
	  cctot = cctot / vvtot
	  cv2(k) = cctot
	end do

	end

c***************************************************************

	subroutine make_basin_aver(nlvddi,nkn,ilhkv,cv3,vol3
     +				,cmin,cmax,cmed,vtot)

	implicit none

	integer nlvddi
	integer nkn
	integer ilhkv(nkn)
	real cv3(nlvddi,nkn)
	real vol3(nlvddi,nkn)
	real cmin,cmax,cmed,vtot

	integer k,l,lmax
	double precision c,v
	double precision cctot,vvtot

	cmin = cv3(1,1)
	cmax = cv3(1,1)
	cctot = 0.
	vvtot = 0.

	do k=1,nkn
	  lmax = ilhkv(k)
	  do l=1,lmax
	    c = cv3(l,k)
	    v = vol3(l,k)
	    cmin = min(cmin,c)
	    cmax = max(cmax,c)
	    cctot = cctot + c*v
	    vvtot = vvtot + v
	  end do
	end do

	cmed = cctot / vvtot
	vtot = vvtot

	end

c***************************************************************

	subroutine make_acumulate(nlvddi,nkn,ilhkv,cv3,cvacu)

	implicit none

	integer nlvddi
	integer nkn
	integer ilhkv(nkn)
	real cv3(nlvddi,nkn)
	real cvacu(nlvddi,nkn)

	integer k,l,lmax

        do k=1,nkn
          lmax = ilhkv(k)
          do l=1,lmax
            cvacu(l,k) = cvacu(l,k) + cv3(l,k)
          end do
        end do

	end

c***************************************************************
c***************************************************************
c***************************************************************

	subroutine write_nos_header(iu,ilhkv,hlv,hev)

c other variables are stored internally
c
c must have been initialized with nos_init
c all other variables must have already been stored internally (title,date..)

	implicit none

	integer iu
	integer ilhkv(1)
	real hlv(1)
	real hev(1)

	integer nkn,nel,nlv,nvar
	integer ierr

	call nos_get_params(iu,nkn,nel,nlv,nvar)
	call nos_write_header(iu,nkn,nel,nlv,nvar,ierr)
	if( ierr .ne. 0 ) goto 99
	call nos_write_header2(iu,ilhkv,hlv,hev,ierr)
	if( ierr .ne. 0 ) goto 99

	return
   99	continue
	write(6,*) 'error in writing header of NOS file'
	stop 'error stop write_nos_header: writing header'
	end

c***************************************************************

	subroutine peek_nos_header(iu,nkn,nel,nlv,nvar)

c get size of data

	implicit none

	integer iu
	integer nkn,nel,nlv,nvar

	integer nvers
	integer ierr

	nvers = 5
	call nos_init(iu,nvers)

	call nos_read_header(iu,nkn,nel,nlv,nvar,ierr)
	if( ierr .ne. 0 ) goto 99

	call nos_close(iu)
	rewind(iu)

	return
   99	continue
	write(6,*) 'error in reading header of NOS file'
	stop 'error stop peek_nos_header: reading header'
	end

c***************************************************************

	subroutine read_nos_header(iu,nknddi,nelddi,nlvddi,ilhkv,hlv,hev)

c other variables are stored internally

	implicit none

	integer iu
	integer nknddi,nelddi,nlvddi
	integer ilhkv(nknddi)
	real hlv(nlvddi)
	real hev(nelddi)

	integer nvers
	integer nkn,nel,nlv,nvar
	integer ierr
	integer l
	integer date,time
	character*50 title,femver

	nvers = 5

	call nos_init(iu,nvers)

	call nos_read_header(iu,nkn,nel,nlv,nvar,ierr)
	if( ierr .ne. 0 ) goto 99

	call dimnos(iu,nknddi,nelddi,nlvddi)

	call getnos(iu,nvers,nkn,nel,nlv,nvar)
	call nos_get_date(iu,date,time)
	call nos_get_title(iu,title)
	call nos_get_femver(iu,femver)

        write(6,*) 'nvers     : ',nvers
        write(6,*) 'nkn,nel   : ',nkn,nel
        write(6,*) 'nlv,nvar  : ',nlv,nvar
        write(6,*) 'title     : ',title
        write(6,*) 'femver    : ',femver
        write(6,*) 'date,time : ',date,time

	call nos_read_header2(iu,ilhkv,hlv,hev,ierr)
	if( ierr .ne. 0 ) goto 99

        write(6,*) 'Available levels: ',nlv
        write(6,*) (hlv(l),l=1,nlv)

	return
   99	continue
	write(6,*) 'error in reading header of NOS file'
	stop 'error stop read_nos_header: reading header'
	end

c***************************************************************
c***************************************************************
c***************************************************************

	subroutine open_nos_type(type,status,nunit)

c open NOS file with default simulation name and given extension

	implicit none

	character*(*) type,status
	integer nunit

	integer nb
	character*80 file

        integer ifileo

	call def_make(type,file)
	nb = ifileo(0,file,'unform',status)

	if( nb .le. 0 ) then
	  write(6,*) 'file: ',file
	  stop 'error stop open_nos_type: opening file'
	end if

	nunit = nb

	end

c***************************************************************

	subroutine open_nos_file(name,status,nunit)

	implicit none

	character*(*) name,status
	integer nunit

	integer nb
	character*80 file

        integer ifileo

	call mkname(' ',name,'.nos',file)
	nb = ifileo(0,file,'unform',status)

	if( nb .le. 0 ) then
	  write(6,*) 'file: ',file
	  stop 'error stop open_nos_file: opening file'
	end if

	nunit = nb

	end

c***************************************************************

        subroutine qopen_nos_file(text,status,nunit)

c asks for name and opens nos file

        implicit none

        character*(*) text,status
        integer nunit

        character*80 name

        write(6,*) text
        read(5,'(a)') name
        write(6,*) name

        call open_nos_file(name,status,nunit)

        end

c***************************************************************
c***************************************************************
c***************************************************************

	subroutine init_volume(nlvddi,nkn,nel,nlv,nen3v,ilhkv
     +				,hlv,hev,hl,vol3)

c initializes volumes just in case no volume file is found
c
c we just set everything to 1.
c we could do better using information on node area and depth structure

	implicit none

	integer nlvddi
	integer nkn,nel,nlv
	integer nen3v(3,nel)
	integer ilhkv(nkn)
	real hlv(nlv)
	real hev(nel)
	real hl(nlv)		!aux vector for layer thickness
	real vol3(nlvddi,nkn)

	logical bvolwrite,bdebug
	integer ie,ii,k,l,lmax,nsigma,nlvaux,ks
	real z,h,hsigma
	double precision ak,vk,ve
	double precision, allocatable :: volk(:,:)

	real weight_elem

        bvolwrite = .true.
        bvolwrite = .false.
	ks = 2985
	ks = 3096
	ks = 0

        call get_sigma_info(nlvaux,nsigma,hsigma)
	z = 0.			!do not use water level

	vol3 = 0.
        allocate(volk(nlvddi,nkn))
	volk = 0.

	do ie=1,nel
	  ak = 4. * weight_elem(ie)	!area of vertex
	  h = hev(ie)
	  call get_layer_thickness(nlv,nsigma,hsigma,z,h,hlv,hl)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    bdebug = k == ks .and. nlv > 1
	    lmax = ilhkv(k)
	    do l=1,lmax
	      vk = ak * hl(l)
	      volk(l,k) = volk(l,k) + vk
	      if( bdebug ) write(82,*) ie,ii,lmax,l,vk,ak,hl(l),h
	    end do
	  end do
          if( bvolwrite .and. mod(ie,nel/nel) == -1 ) then
	   if( nlv > 1 ) then
            write(62,*) ie,h
            write(62,*) hl
           end if
          end if
	end do

	vol3 = volk

	deallocate(volk)

	if( nlv <= 1 ) bvolwrite = .false.
	if( bvolwrite ) then
	write(72,*) 'volume from init_volume'
        do k=1,nkn,nkn/10
          write(72,*) k
          write(72,*) vol3(:,k)
        end do
	end if

	end

c***************************************************************

	subroutine get_volume(nvol,it,nlvddi,ilhkv,vol3)

c reads volumes

	implicit none

	integer nvol
	integer it
	integer nlvddi
	integer ilhkv(1)
	real vol3(nlvddi,1)

	integer ivar,ierr

	integer icall,itold
	save icall,itold
	data icall,itold /0,0/

	if( icall .gt. 0 .and. it .eq. itold ) return	!already read
	if( icall .le. 0 ) itold = it - 1		!force read
	if( it .lt. itold ) goto 95			!error - it too small

	do while( itold .lt. it )
	  call rdnos(nvol,itold,ivar,nlvddi,ilhkv,vol3,ierr)
          if(ierr.gt.0) goto 94			!read error
          if(ierr.ne.0) goto 93			!EOF
          if(ivar.ne.66) goto 92		!ivar should be 66
	end do

	icall = 1
	if( it .lt. itold ) goto 95		!HACK - just temporary
	if( itold .ne. it ) goto 96

	return
   92	continue
	write(6,*) ivar,66
	stop 'error stop get_volume: wrong variable'
   93	continue
	stop 'error stop get_volume: EOF found reading nos file'
   94	continue
	stop 'error stop get_volume: error reading nos file'
   95	continue
	write(6,*) 'it in vol file is higher than requested: ',itold,it
	return		!FIXME -> should signal error
	stop 'error stop get_volume: it too small'
   96	continue
	write(6,*) it,itold
	stop 'error stop get_volume: no vol record for it'
	end

c***************************************************************
c***************************************************************
c***************************************************************

	subroutine check_equal_i(text,n,a1,a2)

c tests array to be equal

	implicit none

	character*(*) text
	integer n
	integer a1(n)
	integer a2(n)

	integer i,imin,imax

	do i=1,n
	  if( a1(i) .ne. a2(i) ) goto 99
	end do

	return
   99	continue
	imin = max(1,i-5)
	imax = min(n,i-5)
	write(6,*) 'first array: ',imin,i,imax
	write(6,*) (a1(i),i=imin,imax)
	write(6,*) 'second array: ',imin,i,imax
	write(6,*) (a2(i),i=imin,imax)
	write(6,*) 'arrays are not equal: ',text
	stop 'error stop check_iqual_i: arrays differ'
	end

c***************************************************************

	subroutine check_equal_r(text,n,a1,a2)

c tests array to be equal

	implicit none

	character*(*) text
	integer n
	real a1(n)
	real a2(n)

	integer i,imin,imax

	do i=1,n
	  if( a1(i) .ne. a2(i) ) goto 99
	end do

	return
   99	continue
	imin = max(1,i-5)
	imax = min(n,i-5)
	write(6,*) 'first array: ',imin,i,imax
	write(6,*) (a1(i),i=imin,imax)
	write(6,*) 'second array: ',imin,i,imax
	write(6,*) (a2(i),i=imin,imax)
	write(6,*) 'arrays are not equal: ',text
	stop 'error stop check_iqual_r: arrays differ'
	end

c***************************************************************
c***************************************************************
c***************************************************************

	subroutine nos_get_it_start(file,itstart)

c gets it of first record

	implicit none

	character*(*) file
	integer itstart

	integer nunit,nvers,nvar
	integer it,ivar,ierr
	character*80 title

	nvers = 5
	itstart = -1

	call open_nos_file(file,'old',nunit)
	if( nunit .le. 0 ) return
	call nos_init(nunit,nvers)
	call nos_skip_header(nunit,nvar,ierr)
	if( ierr .ne. 0 ) goto 1
	call nos_skip_record(nunit,it,ivar,ierr)
	if( ierr .ne. 0 ) goto 1
	itstart = it

    1	continue
	call nos_close(nunit)
	close(nunit)

	end

c***************************************************************

	subroutine nos_get_it_end(file,itend)

c gets it of last record

	implicit none

	character*(*) file
	integer itend

	integer nunit,nvers,nvar
	integer it,itlast,ivar,ierr
	character*80 title

	nvers = 5
	itend = -1
	itlast = -1

	call open_nos_file(file,'old',nunit)
	if( nunit .le. 0 ) return
	call nos_init(nunit,nvers)
	call nos_skip_header(nunit,nvar,ierr)
	if( ierr .ne. 0 ) goto 1

	do
	  call nos_skip_record(nunit,it,ivar,ierr)
	  if( ierr .gt. 0 ) goto 1
	  if( ierr .lt. 0 ) exit
	  itlast = it
	end do
	itend = itlast

    1	continue
	call nos_close(nunit)
	close(nunit)

	end

c***************************************************************

	subroutine nos_get_vars(nin,nvar,ivars)

	implicit none

	integer nin
	integer nvar
	integer ivars(nvar)

	integer i,ivar,it,ierr

	do i=1,nvar
	  call nos_skip_record(nin,it,ivar,ierr)
	  if( ierr .ne. 0 ) goto 99
	  ivars(i) = ivar
	end do

	do i=1,nvar
	  call nos_back_record(nin)
	end do

	return
   99	continue
	write(6,*) 'not enough variables: ',nvar,i,ierr
	stop 'error stop nos_get_vars: nvar'
	end

c***************************************************************

	function check_nos_file(file)

	implicit none

	logical check_nos_file
	character*(*) file

	integer nb,nvers
        integer ifileo

	check_nos_file = .false.

	nb = ifileo(0,file,'unform','old')
	if( nb .le. 0 ) return
	call nos_is_nos_file(nb,nvers)
	close(nb)

	check_nos_file = nvers > 0

	end

c***************************************************************

