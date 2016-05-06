c
c $Id: nosresidence.f,v 1.9 2010-03-22 15:29:31 georg Exp $
c
c analyzes NOS file for residence time
c
c revision log :
c
c 05.03.2004    ggu     copied from nosextr
c 18.10.2006    ggu     heavily commented
c 13.02.2009    ggu     compute min/max of residence time
c 29.10.2009    ggu     compute and write 2D fields, use log to compute restime
c 12.01.2010    ggu     lots of changes
c 29.04.2010    ggu     completely changed
c 10.11.2011    ggu     call to init_volume() changed for hybrid levels
c
c****************************************************************

	program nosresidence

c analyzes CON file for residence time

	use mod_depth
	use evgeom
	use levels
	use basin

	implicit none

	include 'param.h'

c--------------------------------------------------

	character*80 title

	real cv3(nlvdim,nkndim)		!conz read
	real cv3d(nlvdim,nkndim)	!conz of last step
	real cv2(nkndim)		!conz 2d
	real cvres3(nlvdim,nkndim)	!computed RT 3D
	real cvres2(nkndim)		!computed RT 2D
	real cvrestot(nlvdim,nkndim)	!average RT of all replica
	real vol3(nlvdim,nkndim)	!volume
        double precision cvacu(nlvdim,nkndim)	!conz integrated

	integer ilhkv2(nkndim)
	real hev2(neldim)
	real hlv2(nlvdim)
	real hl(nlvdim)

	logical bminmax,balways,breset
	logical blog,badj,bvol
	integer nvol,nb3,nb2
	integer nkn1,nkn2,nel1,nel2,nlv1,nlv2
        integer l,lmax,k
	integer itstart,itend
	integer nin,nvers,nvar
	integer it,ivar,ierr,it0
	integer nread,nused,nrepl
	integer minused
        real conz,c0,conze,conzold
	real dt,ddt,secs_in_day,eps
	real xmin,xmax
	real rl,clog,cl
	real ctop,ccut
	real res
	real cmin,cmax,cmed,vtot
        double precision tacu

	integer iapini,ideffi
	integer ifem_test_file,ifem_open_file,ifem_choose_file

c---------------------------------------------------------------
c parameters to be changed
c---------------------------------------------------------------

c dt		time step used (output)
c c0		initial concentration
c it0		start of simulation (where conz = c0)
c bminmax	true  -> compute only between itstart and itend
c		false -> compute always
c breset	true  -> concentration gets reset during simulation
c itstart	start of computation if bminmax is true
c itend		end of computation if bminmax is true
c
c blog		use logarithmic regression to compute residence time
c badj		adjust residence time for tail of distribution
c
c ctop		maximum to be used for frequency curve
c ccut		cut residence time at this level (cor res time computation)

c        dt = 600
c        dt = 3600
c	c0 = 100.

        dt = 3600
        dt = 86400
	c0 = 1.
	it0 = 0
	ctop = 150.
	ccut = 150.
	ctop = 7000.
	ccut = 7000.
	ctop = 400.
	ccut = 500.
	minused = 5

	bminmax = .true.
	bminmax = .false.
	itstart = 86400
	itend = 30000000

	blog = .false.
	badj = .true.
	badj = .false.

c---------------------------------------------------------------
c do not change anything beyond this point
c---------------------------------------------------------------

	secs_in_day = 86400.
	eps = 0.3
	conzold = 1.e+30

c---------------------------------------------------------------
c open simulation and basin
c---------------------------------------------------------------

	if(iapini(3,nkndim,neldim,0).eq.0) then
		stop 'error stop : iapini'
	end if

	if( nkn > nkndim .or. nel > neldim ) then
	  write(6,*) nkn,nkndim,nel,neldim
	  stop 'error stop: dimensions'
	end if

	call set_ev

        nvers=3

c       ----------------------------------------------------------
c       file containing volumes
c       ----------------------------------------------------------

        nvol = ifem_choose_file('.fvl','old')
        bvol = nvol .gt. 0

        if( bvol ) then
	  write(6,*) 'volume file opened... using it'
          call rhnos(nvol,nvers,nkndim,neldim,nlvdim,nkn2,nel2,nlv2,nvar
     +                          ,ilhkv2,hlv2,hev2,title)
        else
          write(6,*) 'cannot open volume file... doing without'
        end if

c       ----------------------------------------------------------
c	file containing concentrations
c       ----------------------------------------------------------

	nin = ifem_open_file('.con','old')
	if(nin.le.0) goto 100

	write(6,*) 'reading file with concentrations: '
        call rhnos(nin,nvers,nkndim,neldim,nlvdim,nkn1,nel1,nlv,nvar
     +                          ,ilhkv,hlv,hev,title)

	call init_sigma_info(nlv,hlv)

	write(6,*) 'initializing volumes...'
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

	write(6,*) 'compatibility check passed...'

c---------------------------------------------------------------
c initialize variables and arrays
c---------------------------------------------------------------

        it = 0
        ivar = 99

        call open_nos_file('nosres','new',nb3)
        call whnos(nb3,nvers,nkn,nel,nlv,1,ilhkv,hlv,hev,title)

        call open_nos_file('nosres2d','new',nb2)
        call whnos(nb2,nvers,nkn,nel,1,1,ilhkv,hlv,hev,title)

	balways = .not. bminmax

	nread = 0
	nused = 0
	nrepl = 0

	call acu_reset_0(tacu,cvacu,ilhkv)

        do k=1,nkn
          lmax = ilhkv(k)
          do l=1,lmax
            cv3d(l,k) = 0.
            cvrestot(l,k) = 0.
          end do
        end do

        open(76,file='resi.txt',status='unknown',form='formatted')

c---------------------------------------------------------------
c time loop -> accumulate concentrations
c---------------------------------------------------------------

  300   continue

	call rdnos(nin,it,ivar,nlvdim,ilhkv,cv3,ierr)

        if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
        if(ierr.ne.0) goto 100

        if( bvol ) call get_volume(nvol,it,nlvdim,ilhkv,vol3)

	nread=nread+1
	!write(6,*) 'time : ',it,ivar

	if( balways .or. it .ge. itstart .and. it .le. itend ) then

	  nused = nused + 1

	  call make_basin_aver(nlvdim,nkn,ilhkv,cv3,vol3
     +				,cmin,cmax,cmed,vtot)
	  conz = cmed
	  write(68,*) it,cmin,cmed,cmax

c	  -----------------------------------------
c	  compute residence times
c	  -----------------------------------------

	  if( conz - conzold .gt. eps ) then	!reset
	    write(6,*) '-------------------------------------------'
	    write(6,*) 'computing res time: ',it,conzold,conz,nused
	    write(6,*) '-------------------------------------------'
	    if( nused .gt. minused ) then
	      call acu_comp_0(blog,badj,it,dt,c0,ccut,ilhkv,tacu,cvacu
     +				,nb3,nb2
     +				,cv3d,vol3,cvres3,cvres2)
	      call acu_freq_0(it,ctop,ilhkv,cvres3,vol3)
	      nrepl = nrepl + 1
	      call make_acumulate(nlvdim,nkn,ilhkv,cvres3,cvrestot)
	    end if
	    call acu_reset_0(tacu,cvacu,ilhkv)
	    it0 = it
	    nused = 0
	  end if
	  !write(68,*) it,conz,conzold
	  conzold = conz

c	  -----------------------------------------
c	  accumulate for residence time
c	  -----------------------------------------

	  tacu = tacu + (it-it0)

          do k=1,nkn
            lmax = ilhkv(k)
            do l=1,lmax
              conz = cv3(l,k)
              if ( conz .lt. 0. ) conz = 0.
	      if( blog ) then
                if ( conz .gt. c0 ) conz = c0
	        rl = - log(conz/c0)
                cvacu(l,k) = cvacu(l,k) + rl
	      else
                cvacu(l,k) = cvacu(l,k) + conz
	      end if
	      cv3d(l,k) = conz
            end do
          end do

	end if

	goto 300

c---------------------------------------------------------------
c end of time loop
c---------------------------------------------------------------

  100	continue

c---------------------------------------------------------------
c compute last residence time
c---------------------------------------------------------------

	    conz = 0.

	    write(6,*) '-------------------------------------------'
	    write(6,*) 'computing res time: ',it,conzold,conz,nused
	    write(6,*) '-------------------------------------------'

	    if( nused .gt. minused ) then
	      nrepl = nrepl + 1
	      call acu_comp_0(blog,badj,it,dt,c0,ccut,ilhkv,tacu,cvacu
     +				,nb3,nb2
     +				,cv3d,vol3,cvres3,cvres2)
	      call acu_freq_0(it,ctop,ilhkv,cvres3,vol3)
	      call make_acumulate(nlvdim,nkn,ilhkv,cvres3,cvrestot)
	    end if

	    write(6,*) '-------------------------------------------'
	    write(6,*) 'computing final res time: ',nrepl
	    write(6,*) '-------------------------------------------'

	    if( nrepl .gt. 0 ) then
	      it = 0
	      call acu_final(it,ilhkv,nrepl,nb3,nb2,cvrestot,vol3,cv2)
	      call acu_freq_0(it,ctop,ilhkv,cvrestot,vol3)
	    end if

c---------------------------------------------------------------
c finish up
c---------------------------------------------------------------

	write(6,*)
	write(6,*) 'parameters   dt = ',dt,'     c0 = ',c0
	write(6,*) nread,' records read'
	write(6,*) nused,' records used'
        write(6,*) 'levels    : ',lmax
	if( bminmax ) write(6,*) 'min/max time used : ',itstart,itend
	write(6,*)
	write(6,*) 'output written to nosres.nos and nosres2d.nos'
	write(6,*)

        if( .not. bvol ) then
          write(6,*) 'no volume file found: average done without'
          write(6,*)
        end if

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

	stop
   95   continue
        write(6,*) 'error parameters in fvl file : '
        write(6,*) 'nkn: ',nkn,nkn2
        write(6,*) 'nel: ',nel,nel2
        stop 'error stop nosaver: nkn,nel'
   96   continue
        write(6,*) 'error parameters in nos file: '
        write(6,*) 'nkn: ',nkn,nkn1
        write(6,*) 'nel: ',nel,nel1
        stop 'error stop nosaver: nkn,nel'
	end

c***************************************************************
c***************************************************************
c***************************************************************
c***************************************************************
c***************************************************************

	subroutine acu_reset_0(tacu,cvacu,ilhkv)

	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'

	double precision tacu
	double precision cvacu(nlvdim,nkndim)
	integer ilhkv(nkndim)


	integer k,lmax,l

	tacu = 0.
        do k=1,nkn
          lmax = ilhkv(k)
          do l=1,lmax
            cvacu(l,k) = 0.
          end do
        end do

	end

c***************************************************************

	subroutine acu_comp_0(blog,badj,it,dt,c0,ccut,ilhkv,tacu,cvacu
     +				,nb3,nb2
     +				,cv3d,vol3,cv3,cv2)

c compute residence time

	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'

	logical blog,badj
	integer it
	real dt
	real c0
	real ccut  !cut residence time at this level (cor res time computation)
	integer ilhkv(nkndim)
	double precision tacu
	double precision cvacu(nlvdim,nkndim)		!accumulated conz
	integer nb3,nb2
	real cv3d(nlvdim,nkndim)			!last concentration
	real vol3(nlvdim,nkndim)			!volumes
	real cv3(nlvdim,nkndim)				!computed RT 3D
	real cv2(nkndim)				!computed RT 2D


	integer k,lmax,l,ivar,ierr
	real conz,conze,res,rese
	real cmin,cmax,cmed,vtot
	real secs_in_day,ddt

c---------------------------------------------------------------
c set parameters
c---------------------------------------------------------------

	secs_in_day = 86400.

	tacu = tacu / secs_in_day
	ddt = dt / secs_in_day

c---------------------------------------------------------------
c compute residence times -> put in cv3
c---------------------------------------------------------------

        do k=1,nkn
          lmax = ilhkv(k)
          do l=1,lmax
            conz = cvacu(l,k)
	    if( blog ) then
	      conz = tacu / conz
	    else
              conz = ddt * conz / c0		!convert to res time
	      if( badj ) then
	        conze = cv3d(l,k) / c0
	        if( conze .ge. 1 ) conze = 0.
	        if( conze .le. 0 ) conze = 0.
                conz = conz / ( 1. - conze )	!adjusted res time
	      end if
	    end if
	    if ( ccut .gt. 0. .and. conz .gt. ccut ) conz = ccut
            cv3(l,k) = conz
          end do
        end do

c---------------------------------------------------------------
c compute average residence times
c---------------------------------------------------------------

	write(6,*) 'basin wide residence times:'

	call make_basin_aver(nlvdim,nkn,ilhkv,cv3,vol3
     +				,cmin,cmax,cmed,vtot)
	call make_vert_aver(nlvdim,nkn,ilhkv,cv3,vol3,cv2)

	write(6,*) ' (aver/min/max): ',cmed,cmin,cmax

c---------------------------------------------------------------
c write to file and terminal
c---------------------------------------------------------------

	ivar = 99

        call wrnos(nb3,it,ivar,nlvdim,ilhkv,cv3,ierr)
        call wrnos(nb2,it,ivar,1,ilhkv,cv2,ierr)

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

	end

c**********************************************************************

	subroutine acu_final(it,ilhkv,nrepl,nb3,nb2,cv3,vol3,cv2)

	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'

	integer it
	integer ilhkv(nkndim)
	integer nrepl
	integer nb3,nb2
	real cv3(nlvdim,nkndim)
	real vol3(nlvdim,nkndim)			!volumes
	real cv2(nkndim)


	integer k,lmax,l,ivar,ierr
	real cmin,cmax,cmed,vtot

c---------------------------------------------------------------
c compute residence times
c---------------------------------------------------------------

        do k=1,nkn
          lmax = ilhkv(k)
          do l=1,lmax
	    cv3(l,k) = cv3(l,k) / nrepl
	  end do
	end do

c---------------------------------------------------------------
c compute average residence times
c---------------------------------------------------------------

	write(6,*) 'basin wide residence times (average):'

	call make_basin_aver(nlvdim,nkn,ilhkv,cv3,vol3
     +				,cmin,cmax,cmed,vtot)
	call make_vert_aver(nlvdim,nkn,ilhkv,cv3,vol3,cv2)

	write(6,*) ' (aver/min/max): ',cmed,cmin,cmax

c---------------------------------------------------------------
c write to file and terminal
c---------------------------------------------------------------

	ivar = 99

        call wrnos(nb3,it,ivar,nlvdim,ilhkv,cv3,ierr)
        call wrnos(nb2,it,ivar,1,ilhkv,cv2,ierr)

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

	end

c**********************************************************************
c**********************************************************************
c**********************************************************************

	subroutine acu_freq_0(it,ctop,ilhkv,cv3,vol3)

c write histogram

	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'

	integer it
	real ctop			!cut at this value of residence time
	integer ilhkv(nkndim)
	real cv3(nlvdim,nkndim)
	real vol3(nlvdim,nkndim)

	integer ndim
	parameter (ndim=100)


	logical bdebug
	integer k,lmax,l,i,ic
	integer icount(0:ndim)
	double precision dcount(0:ndim)
	real conz,c,amax,dc,val,tot,cmax
	real v,dw
	character*50 file

	bdebug = .true.
	bdebug = .false.

c---------------------------------------------------------------
c compute maximum
c---------------------------------------------------------------

	cmax = 0.
        do k=1,nkn
          lmax = ilhkv(k)
          do l=1,lmax
	    cmax = max(cmax,cv3(l,k))
	  end do
	end do

	amax = cmax
	if( ctop .gt. 0. .and. amax .gt. ctop ) amax = ctop

c---------------------------------------------------------------
c compute frequency curve
c---------------------------------------------------------------

	ic = 0
	dc = 0.
	do i=0,ndim
	  icount(i) = 0
	  dcount(i) = 0.
	end do
	conz = 0.

        do k=1,nkn
          lmax = ilhkv(k)
          do l=1,lmax
	    conz = cv3(l,k)
	    v = vol3(l,k)

	    i = nint(ndim*conz/amax)
	    if (i .lt. 0) i = 0
	    if (i .gt. ndim) i = ndim

	    icount(i) = icount(i) + 1
	    ic = ic + 1
	    dcount(i) = dcount(i) + v
	    dc = dc + v
	  end do
	end do

c---------------------------------------------------------------
c write frequency curve to file
c---------------------------------------------------------------

	dw = 1.
	tot = 0.
	call make_name(it,file,'freq_by_bin_','.his')
	open(11,file=file,status='unknown',form='formatted')
	do i=0,ndim
	  c = i*amax/ndim
	  val = 100.*icount(i)/(float(ic)*dw)
	  tot = tot + val*dw
	  write(11,*) i,val,icount(i)
	end do
	close(11)
	if( bdebug ) write(6,*) 'writing his: ',tot,file

	dw = amax/100.
	tot = 0.
	call make_name(it,file,'freq_by_res_','.his')
	open(11,file=file,status='unknown',form='formatted')
	do i=0,ndim
	  c = i*amax/ndim
	  val = 100.*icount(i)/(float(ic)*dw)
	  tot = tot + val*dw
	  write(11,*) c,val,icount(i)
	end do
	close(11)
	if( bdebug ) write(6,*) 'writing his: ',tot,file

	dw = amax/100.
	tot = 0.
	call make_name(it,file,'freq_by_vol_','.his')
	open(11,file=file,status='unknown',form='formatted')
	do i=0,ndim
	  c = i*amax/ndim
	  val = 100.*dcount(i)/(dc*dw)
	  tot = tot + val*dw
	  write(11,*) c,val,icount(i)
	end do
	close(11)
	if( bdebug ) write(6,*) 'writing his: ',tot,file

c---------------------------------------------------------------
c write out all data to file (for debug and median)
c---------------------------------------------------------------

	write(76,*) nkn,it
        do k=1,nkn
	  conz = cv3(1,k)
	  write(76,*) conz
	end do

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

	end

c**********************************************************************

	subroutine make_name(it,file,pre,ext)

	implicit none

	integer it
	character*(*) file,pre,ext

	integer i,n
	integer ichanm

	write(file,'(a,i10,a)') pre,it,ext
	n = ichanm(file)
	do i=n,1,-1
	  if( file(i:i) .eq. ' ' ) file(i:i) = '0'
	end do

	end

c**********************************************************************
c**********************************************************************
c**********************************************************************
c**********************************************************************
c**********************************************************************

