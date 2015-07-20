c
c $Id: nosaverf.f,v 1.1 2008-04-11 16:05:37 georg Exp $
c
c revision log :
c
c 18.11.1998    ggu     check dimensions with dimnos
c 07.03.2007    ggu     easier calls
c 10.04.2008    ggu     copied from nosaver -> frequency introduced
c 29.04.2010    ggu     rewritten (using volumes)
c 10.11.2011    ggu     call to init_volume() changed for hybrid levels
c
c**********************************************************

	program nosaver_basin

c averages records over basin
c
c creates 4 values: aver, min, max, sum

	use mod_depth
	use basin

        implicit none

	include 'param.h'

	character*80 title
	real cv3(nlvdim,nkndim)
	real vol3(nlvdim,nkndim)

	integer ilhkv(nkndim)
	real hlv(nlvdim)
	integer ilhkv2(nkndim)
	real hlv2(nlvdim)
	real hev2(neldim)
	real hl(nlvdim)


	logical bvol
        integer nread
        integer l,k,nin,nvol
        integer nlv,nvar
	integer nkn1,nkn2,nel1,nel2,nlv2
        integer it,ivar
        integer ierr
        integer nvers
	real cmin,cmax,cmed,vtot,totmass

        integer iapini,ideffi,ifileo
        integer ifem_open_file,ifem_test_file,ifem_choose_file

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
	  !write(6,*) 'time : ',it,ivar

c	  --------------------------------------------------------
c	  compute average over basin
c	  --------------------------------------------------------

	  call make_aver(nlvdim,nkn,ilhkv,cv3,vol3,cmin,cmax,cmed,vtot)

c	  --------------------------------------------------------
c	  write output
c	  --------------------------------------------------------

	  totmass = cmed * vtot
	  write(6,1234) it,ivar,cmin,cmed,cmax,totmass
	  write(100+ivar,1235) it,cmin,cmed,cmax,totmass
	  write(100,'(i10,e14.6)') it,vtot

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
 1234	format(2i10,3f12.4,e14.6)
 1235	format(i10,3f12.4,e14.6)
	end

c***************************************************************

