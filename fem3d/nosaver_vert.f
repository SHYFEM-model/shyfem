c
c $Id: nosaverf.f,v 1.1 2008-04-11 16:05:37 georg Exp $
c
c revision log :
c
c 18.11.1998    ggu     check dimensions with dimnos
c 07.03.2007    ggu     easier calls
c 10.04.2008    ggu     copied from nosaver -> frequency introduced
c 29.04.2010    ggu     new from nosaver_basin (using volumes)
c 07.05.2010    ggu     call whnos with nlv=1 (bug)
c 10.11.2011    ggu     call to init_volume() changed for hybrid levels
c
c**********************************************************

	program nosaver_vert

c averages vertically records

        implicit none

	include 'param.h'
	include 'basin.h'

	character*80 title
	real cv3(nlvdim,nkndim)
	real cv2(nkndim)
	real vol3(nlvdim,nkndim)

	integer ilhkv(nkndim)
	real hlv(nlvdim)
	integer ilhkv2(nkndim)
	real hlv2(nlvdim)
	real hev2(neldim)
	real hl(nlvdim)

	include 'depth.h'

	logical bvol
        integer nread
        integer l,k,nin,nvol,nb2
        integer nlv,nvar
	integer nkn1,nkn2,nel1,nel2,nlv2
        integer it,ivar
        integer ierr
        integer nvers
	real cmin,cmax,cmed,vtot

        integer iapini,ideffi,ifileo,ifem_open_file,ifem_choose_file

c-----------------------------------------------------------------
c initialize params
c-----------------------------------------------------------------

	nread = 0

c-----------------------------------------------------------------
c open input files
c-----------------------------------------------------------------

	if(iapini(3,nkndim,neldim,0).eq.0) then
		stop 'error stop : iapini'
	end if

        nvers=3

c	----------------------------------------------------------
c	file containing volumes
c	----------------------------------------------------------

        nvol = ifem_choose_file('.fvl','old')
	bvol = nvol .gt. 0

	if( bvol ) then
	  write(6,*) 'volume file opened... using it'
          call rhnos(nvol,nvers,nkndim,neldim,nlvdim,nkn2,nel2,nlv2,nvar
     +                          ,ilhkv2,hlv2,hev2,title)
	else
	  write(6,*) 'cannot open volume file... doing without'
	end if

c	----------------------------------------------------------
c	file containing variables to average
c	----------------------------------------------------------

        nin = ifem_open_file('.nos','old')
        call rhnos(nin,nvers,nkndim,neldim,nlvdim,nkn1,nel1,nlv,nvar
     +                          ,ilhkv,hlv,hev,title)

        call init_sigma_info(nlv,hlv)

        call init_volume(nlvdim,nkn,nel,nlv,nen3v,ilhkv,hlv,hev,hl,vol3)

c	----------------------------------------------------------
c	output file
c	----------------------------------------------------------

        call open_nos_file('nos2d','new',nb2)
        call whnos(nb2,nvers,nkn,nel,1,nvar,ilhkv,hlv,hev,title)

c-----------------------------------------------------------------
c check compatibility
c-----------------------------------------------------------------

	if( nkn .ne. nkn1 ) goto 96
	if( nel .ne. nel1 ) goto 96

	if( bvol ) then
	  if( nkn .ne. nkn2 ) goto 95
	  if( nel .ne. nel2 ) goto 95
	  call check_equal_i('ilhkv',nkn,ilhkv,ilhkv2)
	  call check_equal_r('hlv',nlv,hlv,hlv2)
	  call check_equal_r('hev',nel,hev,hev2)
	end if

c-----------------------------------------------------------------
c loop on input records
c-----------------------------------------------------------------

	do while(.true.)

c	  --------------------------------------------------------
c	  read next record
c	  --------------------------------------------------------

	  call rdnos(nin,it,ivar,nlvdim,ilhkv,cv3,ierr)
          if(ierr.gt.0) goto 94
          if(ierr.ne.0) goto 100

	  if( bvol ) call get_volume(nvol,it,nlvdim,ilhkv,vol3)

	  nread=nread+1
	  write(6,*) 'time : ',it,ivar

c	  --------------------------------------------------------
c	  compute vertical average
c	  --------------------------------------------------------

	  call make_vert_aver(nlvdim,nkn,ilhkv,cv3,vol3,cv2)

c	  --------------------------------------------------------
c	  write output
c	  --------------------------------------------------------

          call wrnos(nb2,it,ivar,1,ilhkv,cv2,ierr)

	end do

c-----------------------------------------------------------------
c end of loop on input records
c-----------------------------------------------------------------

  100	continue

        if( .not. bvol ) then
          write(6,*)
          write(6,*) 'no volume file found: average done without'
        end if

	write(6,*)
	write(6,*) 'new 2D file written to nos2d.nos'
	write(6,*)
	write(6,*) nread,' records read in total'
	write(6,*)

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	stop
   94	continue
	stop 'error stop nosaver: error reading nos file'
   95	continue
	write(6,*) 'error parameters in fvl file : '
	write(6,*) 'nkn: ',nkn,nkn2
	write(6,*) 'nel: ',nel,nel2
	stop 'error stop nosaver: nkn,nel'
   96	continue
	write(6,*) 'error parameters in nos file: '
	write(6,*) 'nkn: ',nkn,nkn1
	write(6,*) 'nel: ',nel,nel1
	stop 'error stop nosaver: nkn,nel'
	end

c***************************************************************

