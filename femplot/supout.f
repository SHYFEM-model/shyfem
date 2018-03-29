c
c $Id: supout.f,v 1.15 2010-02-26 15:29:19 georg Exp $
c
c routines for reading data files
c
c revision log :
c
c 31.10.2003  ggu     new routines velclose(), resetsim()
c 31.10.2003  ggu     new routines for handling wind
c 22.09.2004  ggu     new routines for handling ous file
c 22.09.2004  ggu     use OUS file for 3D, try to use alsways 3D (level=0)
c 05.10.2004  ggu     adjustments -> use always 3D data structure
c 14.03.2007  ggu     new routines for wave plotting
c 17.09.2008  ggu     new routine level_e2k to compute ilhkv from ilhv
c 09.10.2009  ggu     read also pressure from wind file
c 13.10.2009  ggu     set nlv once file is read
c 23.02.2010  ggu     change in reading wind file
c 30.03.2011  ggu     new routines to handle fvl files (not yet integrated)
c 31.03.2011  ggu     use fvl routines to exclude areas from plot
c 12.07.2011  ggu     in prepsim use what is available for dry areas
c 18.08.2011  ggu     bug fix in nosopen() -> extra comma eliminated 
c 31.08.2011  ggu     new routines for handling EOS files
c 14.11.2011  ggu     call to init_sigma_info() to setup layer info
c 19.12.2011  ggu     new routine level_k2e -> called for nos files
c 13.06.2013  ggu     new routines for handling FEM files
c 03.09.2013  ggu     level_k2e -> level_k2e_sh, level_e2k -> level_e2k_sh
c 05.03.2014  ggu     new read for ous and nos files (use date)
c 20.10.2014  ggu     deleted is2d() and out reading routines
c 10.02.2015  ggu     use different file units (more than one can be opened)
c 14.09.2015  ggu     introduced bwind and bvel (for velocities)
c
c**********************************************************
c**********************************************************
c**********************************************************

        subroutine velopen

c opens velocity file (2d and 3d)

        implicit none

        call ousopen

        end

c**********************************************************

        function velnext(it)

c gets next velocity field (2d and 3d)

        implicit none

        logical velnext
        integer it

        logical outnext, ousnext

        velnext = ousnext(it)

        end

c**********************************************************

        subroutine velclose

c closes velocity file

        implicit none

        call ousclose

        end

c******************************************************
c******************************************************
c******************************************************

        subroutine reset_dry_mask

c resets mask data structure

	use mod_hydro_plot

	implicit none

        bwater = .true.
        bkwater = .true.
        bplot = .true.

        end

c******************************************************

        subroutine adjust_no_plot_area

        use basin
	use mod_hydro_plot

        implicit none

        integer ie
        integer ianopl
	real getpar

        ianopl = nint(getpar('ianopl'))
	if( ianopl < 0 ) return

        !write(6,*) 'applying no plot mask for area code = ',ianopl

        do ie=1,nel
          if( iarv(ie) == ianopl ) then
	    bwater(ie) = .false.
	    bplot(ie) = .false.
	  end if
        end do

        end

c******************************************************

	subroutine prepare_dry_mask

c prepares simulation for use - computes wet and dry areas

	use mod_hydro_plot
	use mod_hydro
	use levels

	implicit none

	logical bshowdry
	integer level
	real href,hzmin,hdry

	logical fvl_is_available
	logical ous_is_available
	integer getlev
	real getpar

c---------------------------------------------------
c set up mask of water points
c---------------------------------------------------

c set bshowdry = .true.  if you want to show dried out areas
c set bshowdry = .false. if you want to plot all areas

	hdry = 0.05				!level for drying
	bshowdry = .true.			!show dry areas, else plot all
	!bshowdry = .false.			!show dry areas, else plot all

	if( .not. bshowdry ) hdry = -1.e+30	!no drying

        href = getpar('href')
        hzmin = getpar('hzmin')
	level = getlev()

	call reset_dry_mask

	if( ous_is_available() ) then			!...handle on elements

	  write(6,*) 'using zeta for dry areas'
	  if( bshowdry ) then
            call set_dry_mask(bwater,znv,zenv,href,hzmin) !false if znv/=zenv
	  end if
          call set_level_mask(bwater,ilhv,level)	!element has this level
	  call make_dry_node_mask(bwater,bkwater)	!copy elem to node mask

	else if( fvl_is_available() ) then		!...handle on nodes

	  write(6,*) 'using fvl file for dry areas: ',hdry
	  call set_dry_volume_mask(bkwater,hdry)	!guess if dry using vol
	  call make_dry_elem_mask(bwater,bkwater)	!copy node to elem mask
          call set_level_mask(bwater,ilhv,level)	!element has this level
	  call make_dry_node_mask(bwater,bkwater)	!copy elem to node mask

	else

	  write(6,*) 'no information on dry areas: ',hdry
          call set_level_mask(bwater,ilhv,level)	!element has this level
	  call make_dry_node_mask(bwater,bkwater)	!copy elem to node mask

	end if

        call adjust_no_plot_area
	call make_dry_node_mask(bwater,bkwater)		!copy elem to node mask
        call info_dry_mask(bwater,bkwater)

c---------------------------------------------------
c end of routine
c---------------------------------------------------

	end

c******************************************************

	subroutine set_dry_volume_mask(bkw,hdry)

	use mod_hydro_plot
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	logical bkw(1)
	real hdry

	integer k,idry
	real vol,area,h
	real hhmin,hhmax

	hhmin = 1.e+30
	hhmax = -1.e+30
	idry = 0

	do k=1,nkn
	  vol = fvlv(1,k)
	  area = arfvlv(k)
	  h = vol/area
	  hhmin = min(hhmin,h)
	  hhmax = max(hhmax,h)
	  bkw(k) = h .ge. hdry
	  if( .not. bkw(k) ) idry = idry + 1
	end do

	write(6,*) 'min/max depth: ',hhmin,hhmax,hdry,idry

	end

c******************************************************
c******************************************************
c******************************************************

        subroutine waveini

        implicit none

	include 'supout.h'

	integer icall
	save icall
	data icall /0/

	if( bdebug_out ) then
	  write(6,*) 'debug_out: waveini'
	  write(6,*) icall,nunit_fvl
	end if

	if( icall .ne. 0 ) return

	icall = 1

	nunit_wave = 0
        iwave = 0	! 0 = wave height   1 = wave period

        end

c******************************************************

        subroutine waveclose

        implicit none

	include 'supout.h'

	if( nunit_wave .gt. 0 ) close(nunit_wave)
	nunit_wave = 0

        end

c******************************************************

	subroutine waveopen

	use mod_depth
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'simul.h'
	include 'supout.h'

	integer nknaux,nelaux,nlvaux
	integer nknddi,nelddi,nlvddi
	integer nvers,nvar
	integer date,time

	integer ifemop

	call waveini

c open file

        call open_nos_type('.wav','old',nunit)

	call peek_nos_header(nunit,nknaux,nelaux,nlvaux,nvar)
	if( nkn .ne. nknaux ) goto 99
	if( nel .ne. nelaux ) goto 99
        nlv = nlvaux
	call allocate_simulation(0)
	call get_dimension_post(nknddi,nelddi,nlvddi)
        call read_nos_header(nunit,nknddi,nelddi,nlvddi,ilhkv,hlv,hev)

	if( nvar .ne. 3 ) goto 99

	call level_k2e_sh
        call init_sigma_info(nlv,hlv)

c initialize time

        call nos_get_date(nunit,date,time)
        if( date .ne. 0 ) then
          call dtsini(date,time)
	  call ptime_set_date_time(date,time)
	  call ptime_min_max
	  call iff_init_global_date(date,time)
        end if

	nunit_wave = nunit

	return
   99	continue
	write(6,*) nkn,nel,nlv
	write(6,*) nknaux,nelaux,nlvaux
	write(6,*) nvar
	stop 'error stop waveopen: parameter mismatch'
        end

c******************************************************

	function wavenext(it)

	use mod_hydro_plot
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	logical wavenext
	integer it

	include 'supout.h'

        integer ierr,nlvddi,ivar
	real v1v(nkn)
	real v2v(nkn)
	real v3v(nkn)


	call waveini
	nunit = nunit_wave

	wavenext = .false.
	nlvddi = 1
	call rdnos(nunit,it,ivar,nlvddi,ilhkv,v1v,ierr)
	if( ierr .gt. 0 ) goto 99
	if( ierr .lt. 0 ) return
	if( ivar .ne. 31 ) goto 97
	call rdnos(nunit,it,ivar,nlvddi,ilhkv,v2v,ierr)
	if( ierr .gt. 0 ) goto 99
	if( ierr .lt. 0 ) goto 98
	if( ivar .ne. 32 ) goto 97
	call rdnos(nunit,it,ivar,nlvddi,ilhkv,v3v,ierr)
	if( ierr .gt. 0 ) goto 99
	if( ierr .lt. 0 ) goto 98
	if( ivar .ne. 33 ) goto 97

	wavenext = .true.

	if( iwave .eq. 0 ) then
	  call polar2xy(nkn,v1v,v3v,uvnode,vvnode)
	else
	  call polar2xy(nkn,v2v,v3v,uvnode,vvnode)
	end if

	call ptime_set_itime(it)

	return
   99	continue
	stop 'error stop wavenext: error reading data'
   98	continue
	stop 'error stop wavenext: not enough records'
   97	continue
	stop 'error stop wavenext: records not in order'
        end

c******************************************************

	subroutine polar2xy(n,speed,dir,uv,vv)

	implicit none

	integer n
	real speed(1), dir(1)
	real uv(1), vv(1)

	integer i
	real rad,a

	rad = atan(1.) / 45.

	do i=1,n
	  a = dir(i)
          a = 90. - a + 180.
          do while( a .lt. 0. )
            a = a + 360.
          end do
          a = mod(a,360.)

	  uv(i) = speed(i) * cos( rad * a )
	  vv(i) = speed(i) * sin( rad * a )
	end do

	end

c******************************************************
c******************************************************
c******************************************************

	subroutine ousini

c initializes internal data structure for OUS file

	implicit none

	include 'supout.h'

	integer icall
	save icall
	data icall /0/

	if( bdebug_out ) then
	  write(6,*) 'debug_out: ousini'
	  write(6,*) icall,nunit_fvl
	end if

	if( icall .ne. 0 ) return

	icall = 1

	nunit_ous = 0

	end

c******************************************************

	function ous_is_available()

c checks if OUS file is opened

	implicit none

	logical ous_is_available

	include 'supout.h'

	call ousini
	ous_is_available = nunit_ous .gt. 0

	end

c******************************************************

	subroutine ousclose

c closes OUS file

	implicit none

	include 'supout.h'

	call ousini
	if( nunit_ous .gt. 0 ) close(nunit_ous)
	nunit_ous = 0

	end

c******************************************************

        subroutine ousinfo(nvers,nkn,nel,nlv)

c returns info on OUS parameters

        implicit none

	include 'supout.h'

        integer nvers,nkn,nel,nlv

	call ousini
        call getous(nunit_ous,nvers,nkn,nel,nlv)

        end

c******************************************************

	subroutine ousopen

c opens OUS file and reads header

	use mod_depth
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'supout.h'
	include 'simul.h'

	character*80 file

	integer nvers
	integer nknaux,nelaux,nlvaux,nvar
	integer nknous,nelous,nlvous
	integer nknddi,nelddi,nlvddi
	integer ierr
        integer i,l
	integer date,time
        real href,hzmin
	integer ifemop

c initialize routines

	call ousini

c open and read header

        call open_ous_type('.ous','old',nunit)

	call peek_ous_header(nunit,nknous,nelous,nlvous)
	if( nkn .ne. nknous ) goto 99
	if( nel .ne. nelous ) goto 99
        nlv = nlvous
	call allocate_simulation(0)
	call get_dimension_post(nknddi,nelddi,nlvddi)
        call read_ous_header(nunit,nknddi,nelddi,nlvddi,ilhv,hlv,hev)

	call level_e2k_sh			!computes ilhkv
        call init_sigma_info(nlv,hlv)

c initialize time

	call ous_get_date(nunit,date,time)
	if( date .ne. 0 ) then
	  call dtsini(date,time)
	  call ptime_set_date_time(date,time)
	  call ptime_min_max
	  call iff_init_global_date(date,time)
	end if

	nunit_ous = nunit

c end

	return
   96	continue
	stop 'error stop ousopen: cannot open OUS file'
   97	continue
	stop 'error stop ousopen: error reading second header'
   98	continue
	stop 'error stop ousopen: error reading first header'
   99	continue
	write(6,*) 'error in parameters :'
	write(6,*) 'nkn : ',nkn,nknous
	write(6,*) 'nel : ',nel,nelous
	write(6,*) 'nlv : ',nlvdi,nlvous
	stop 'error stop ousopen'
	end

c******************************************************

	function ousnext(it)

c reads next OUS record - is true if a record has been read, false if EOF

	use mod_hydro
	use levels

	implicit none

	logical ousnext		!true if record read, flase if EOF
	integer it		!time of record

	include 'supout.h'

	integer ierr

	call ousini
	nunit = nunit_ous

	call rdous(nunit,it,nlvdi,ilhv,znv,zenv,utlnv,vtlnv,ierr)

c set return value

	if( ierr .gt. 0 ) then
		!stop 'error stop ousnext: error reading data record'
		write(6,*) '*** ousnext: error reading data record'
		ousnext = .false.
	else if( ierr .lt. 0 ) then
		ousnext = .false.
	else
		ousnext = .true.
	end if

	call ptime_set_itime(it)

	end

c******************************************************
c******************************************************
c******************************************************
c 3d-version of concentration
c******************************************************
c******************************************************
c******************************************************

	subroutine nosini

c initializes internal data structure for NOS file

	implicit none

	include 'supout.h'

	integer icall
	save icall
	data icall /0/

	if( bdebug_out ) then
	  write(6,*) 'debug_out: nosini'
	  write(6,*) icall,nunit_fvl
	end if

	if( icall .ne. 0 ) return

	icall = 1

	nunit_nos = 0

	end

c******************************************************

	subroutine nosclose

c closes NOS file

	implicit none

	include 'supout.h'

	call nosini
	if( nunit_nos .gt. 0 ) close(nunit_nos)
	nunit_nos = 0

	end

c******************************************************

	subroutine nosopen(type)

c opens NOS file and reads header

	use mod_depth
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	character*(*) type

	include 'supout.h'
	include 'simul.h'

	character*80 file

	integer nvers
	integer nknaux,nelaux,nlvaux,nvar
	integer nknnos,nelnos,nlvnos
	integer nknddi,nelddi,nlvddi
	integer date,time
	integer ierr,l
	integer ifemop

c initialize routines

	call nosini

c open file

	if( type /= ' ' ) then
          call open_nos_type(type,'old',nunit)
	else
          call open_nos_type('.nos','old',nunit)
	end if

	call peek_nos_header(nunit,nknnos,nelnos,nlvnos,nvar)
	if( nkn .ne. nknnos ) goto 99
	if( nel .ne. nelnos ) goto 99
        nlv = nlvnos
	call allocate_simulation(0)
	call get_dimension_post(nknddi,nelddi,nlvddi)
        call read_nos_header(nunit,nknddi,nelddi,nlvddi,ilhkv,hlv,hev)

	call level_k2e_sh
        call init_sigma_info(nlv,hlv)

c initialize time

        call nos_get_date(nunit,date,time)
        if( date .ne. 0 ) then
          call dtsini(date,time)
	  call ptime_set_date_time(date,time)
	  call ptime_min_max
	  call iff_init_global_date(date,time)
        end if

	nunit_nos = nunit

c end

	!write(6,*) 'gguuuuuu: (nosopen)',nunit,nunit_fvl

	return
   99	continue
	write(6,*) 'error in parameters : basin - simulation'
	write(6,*) 'nkn : ',nkn,nknnos
	write(6,*) 'nel : ',nel,nelnos
	write(6,*) 'nlv : ',nlvdi,nlvnos
	write(6,*) 'parameters are different between basin and simulation'
	stop 'error stop nosopen'
	end

c******************************************************

	function nosnext(it,ivar,nlvddi,array)

c reads next NOS record - is true if a record has been read, false if EOF

	use levels

	implicit none

	logical nosnext		!true if record read, flase if EOF
	integer it		!time of record
	integer ivar		!type of variable
	integer nlvddi		!dimension of vertical coordinate
	real array(nlvddi,1)	!values for variable

	include 'supout.h'

	integer ierr

	if( nlvddi .ne. nlvdi ) stop 'error stop nosnext: nlvddi'

	call nosini
	nunit = nunit_nos
	!write(6,*) 'gguuuuuu: (nosnext)',nunit,nunit_fvl

	call rdnos(nunit,it,ivar,nlvddi,ilhkv,array,ierr)

c set return value

	if( ierr .gt. 0 ) then
		!stop 'error stop nosnext: error reading data record'
		write(6,*) '*** nosnext: error reading data record'
		nosnext = .false.
	else if( ierr .lt. 0 ) then
		nosnext = .false.
	else
		nosnext = .true.
	end if

	call ptime_set_itime(it)

	end

c******************************************************
c******************************************************
c******************************************************
c routines to read fvl file
c******************************************************
c******************************************************
c******************************************************

	subroutine fvlini

c initializes internal data structure for FVL file

	implicit none

	include 'supout.h'

	integer icall
	save icall
	data icall /0/

	if( bdebug_out ) then
	  write(6,*) 'debug_out: fvlini'
	  write(6,*) icall,nunit_fvl
	end if

	if( icall .ne. 0 ) return

	icall = 1

	nunit_fvl = 0

	end

c******************************************************

	function fvl_is_available()

c checks if FVL file is opened

	implicit none

	logical fvl_is_available

	include 'supout.h'

	call fvlini
	fvl_is_available = nunit_fvl .gt. 0

	end

c******************************************************

	subroutine fvlclose

c closes FVL file

	implicit none

	include 'supout.h'

	call fvlini
	if( nunit_fvl .gt. 0 ) close(nunit_fvl)
	nunit_fvl = 0

	end

c******************************************************

	subroutine fvlopen(type)

c opens FVL file and reads header

	use mod_depth
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	character*(*) type

	include 'supout.h'
	include 'simul.h'

        real hlv1(nlv), hev1(nel)
        integer ilhkv1(nkn)

	character*80 file

	integer nvers
	integer nknaux,nelaux,nlvaux,nvar
	integer ierr
	integer ifem_choose_file

c initialize routines

	call fvlini

c open file

	nunit = ifem_choose_file(type,'old')
	if( nunit .le. 0 ) then
		write(6,*) 'Cannot open fvl file ... doing without...'
		nunit_fvl = 0
		return
		!stop 'error stop conopen: cannot open NOS file'
        else
                write(6,*) 'fvl file opened :'
                inquire(nunit,name=file)
                write(6,*) file
                write(6,*) 'Reading file ...'
	end if

c read first header

	nvers = 3
	call rfnos(nunit,nvers,nknaux,nelaux,nlvaux,nvar,descrp,ierr)

	if( ierr .ne. 0 ) then
		stop 'error stop nosopen: error reading first header'
	end if

        write(6,*)
        write(6,*)   descrp
        write(6,*)
        write(6,*) ' nvers = ', nvers
        write(6,*) '   nkn = ',nknaux, '   nel = ',nelaux
        write(6,*) '   nlv = ',nlvaux, '  nvar = ',  nvar
        write(6,*)

	if( nkn .ne. nknaux ) goto 99
	if( nelaux .ne. 0 .and. nel .ne. nelaux ) goto 99
	if( nlv .ne. nlvaux ) goto 99
	if( nvar .ne. 1 ) goto 99

c read second header

	call rsnos(nunit,ilhkv1,hlv1,hev1,ierr)

	if( ierr .ne. 0 ) then
		stop 'error stop nosopen: error reading second header'
	end if

	call array_i_check(nkn,ilhkv1,ilhkv,'ilhkv')
	call array_check(nlv,hlv1,hlv,'hlv')
	call array_check(nel,hev1,hev,'hev')

	nunit_fvl = nunit

c end

	return
   99	continue
	write(6,*) 'error in parameters :'
	write(6,*) 'nkn : ',nkn,nknaux
	write(6,*) 'nel : ',nel,nelaux
	write(6,*) 'nlv : ',nlv,nlvaux
	write(6,*) 'nvar: ',nvar
	write(6,*) 'not using FVL file'
	nunit_fvl = 0
	!stop 'error stop fvlopen'
	end

c******************************************************

	subroutine fvlnext(it,nlvddi,array)

c reads next FVL record

	use levels

	implicit none

	integer it		!time of record
	integer nlvddi		!dimension of vertical coordinate
	real array(nlvddi,1)	!values for variable

	include 'supout.h'

	integer ierr
	integer ivar		!type of variable

	integer itfvl
	save itfvl
	data itfvl / -1 /

	if( nlvddi .ne. nlvdi ) stop 'error stop fvlnext: nlvddi'

	call fvlini
	nunit = nunit_fvl

	if( nunit .eq. 0 ) return	!file closed
	if( it .eq. itfvl ) return	!already read

	if( itfvl .eq. -1 ) itfvl = it - 1	!force first read
	nunit = abs(nunit)
	ierr = 0
	ivar = 66

	do while( ierr .eq. 0 .and. itfvl .lt. it )
	  call rdnos(nunit,itfvl,ivar,nlvddi,ilhkv,array,ierr)
	  write(6,*) 'fvl file read...',itfvl,ivar
	end do

c set return value

	if( ierr .ne. 0 ) then
	  if( ierr .gt. 0 ) then
		!stop 'error stop fvlnext: error reading data record'
		write(6,*) '*** fvlnext: error reading data record'
	  else if( ierr .lt. 0 ) then
		write(6,*) '*** fvlnext: EOF encountered'
	  end if
	  itfvl = -1
	  nunit_fvl = 0
	  return
	end if

c check results

	if( it .ne. itfvl ) nunit_fvl = -abs(nunit_fvl)
	if( ivar .ne. 66 ) goto 99

c end

	return
   99	continue
	write(6,*) it,itfvl,ivar
	stop 'error stop fvlnext: time or variable'
	end

c******************************************************
c******************************************************
c******************************************************
c element values
c******************************************************
c******************************************************
c******************************************************

	subroutine eosini

c initializes internal data structure for EOS file

	implicit none

	include 'supout.h'

	integer icall
	save icall
	data icall /0/

	if( bdebug_out ) then
	  write(6,*) 'debug_out: eosini'
	  write(6,*) icall,nunit_fvl
	end if

	if( icall .ne. 0 ) return

	icall = 1

	nunit_eos = 0

	end

c******************************************************

	subroutine eosclose

c closes EOS file

	implicit none

	include 'supout.h'

	call eosini
	if( nunit_eos .gt. 0 ) close(nunit_eos)
	nunit_eos = 0

	end

c******************************************************

	subroutine eosopen(type)

c opens EOS file and reads header

	use mod_depth
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	character*(*) type

	include 'supout.h'
	include 'simul.h'

	character*80 file

	integer nvers
	integer nknaux,nelaux,nlvaux,nvar
	integer nknddi,nelddi,nlvddi
	integer ierr
	integer ifemop

c initialize routines

	call eosini

c open file

	nunit = ifemop(type,'unform','old')
	if( nunit .le. 0 ) then
		stop 'error stop conopen: cannot open EOS file'
        else
                write(6,*) 'File opened :'
                inquire(nunit,name=file)
                write(6,*) file
                write(6,*) 'Reading file ...'
	end if

c read first header

	nvers = 3
	call rfeos(nunit,nvers,nknaux,nelaux,nlvaux,nvar,descrp,ierr)

	if( ierr .ne. 0 ) then
		stop 'error stop eosopen: error reading first header'
	end if

        write(6,*)
        write(6,*)   descrp
        write(6,*)
        write(6,*) ' nvers = ', nvers
        write(6,*) '   nkn = ',nknaux, '   nel = ',nelaux
        write(6,*) '   nlv = ',nlvaux, '  nvar = ',  nvar
        write(6,*)

	if( nkn .ne. nknaux ) goto 99
	if( nelaux .ne. 0 .and. nel .ne. nelaux ) goto 99

	nlv = nlvaux
	call allocate_simulation(0)
	call get_dimension_post(nknddi,nelddi,nlvddi)
	if( nlvddi .lt. nlv ) goto 99

c read second header

	call rseos(nunit,ilhv,hlv,hev,ierr)

	if( ierr .ne. 0 ) then
		stop 'error stop eosopen: error reading second header'
	end if

	call level_e2k_sh		!computes ilhkv

	nunit_eos = nunit

	return
   99	continue
	write(6,*) 'error in parameters :'
	write(6,*) 'nkn : ',nkn,nknaux
	write(6,*) 'nel : ',nel,nelaux
	write(6,*) 'nlv : ',nlvdi,nlvaux
	stop 'error stop eosopen'
	end

c******************************************************

	function eosnext(it,ivar,nlvddi,array)

c reads next EOS record - is true if a record has been read, false if EOF

	use levels

	implicit none

	logical eosnext		!true if record read, flase if EOF
	integer it		!time of record
	integer ivar		!type of variable
	integer nlvddi		!dimension of vertical coordinate
	real array(nlvddi,1)	!values for variable

	include 'supout.h'

	integer ierr

	if( nlvddi .ne. nlvdi ) stop 'error stop eosnext: nlvddi'

	call eosini
	nunit = nunit_eos

	call rdeos(nunit,it,ivar,nlvddi,ilhv,array,ierr)

c set return value

	if( ierr .gt. 0 ) then
		!stop 'error stop eosnext: error reading data record'
		write(6,*) '*** eosnext: error reading data record'
		eosnext = .false.
	else if( ierr .lt. 0 ) then
		eosnext = .false.
	else
		eosnext = .true.
	end if

	call ptime_set_itime(it)

	end

c******************************************************
c******************************************************
c******************************************************
c aux routines
c******************************************************
c******************************************************
c******************************************************

	subroutine level_e2k_sh

c computes max level at nodes from elements

	use levels
	use basin

	implicit none

	call level_e2k(nkn,nel,nen3v,ilhv,ilhkv)

	end

c******************************************************

	subroutine level_k2e_sh

c computes level at elems from nodes (not exact)

	use levels
	use basin

	implicit none

	call level_k2e(nkn,nel,nen3v,ilhkv,ilhv)

	end

c******************************************************

	subroutine array_check(n,a1,a2,text)

c checks if arrays are equal

	implicit none

	integer n
	real a1(n)
	real a2(n)
	character*(*) text

	integer i

	do i=1,n
	  if( a1(i) .ne. a2(i) ) then
	    write(6,*) text,' : arrays differ'
	    write(6,*) n,i,a1(i),a2(i)
	    stop 'error stop array_check: arrays differ'
	  end if
	end do

	end

c******************************************************

	subroutine array_i_check(n,a1,a2,text)

c checks if arrays are equal

	implicit none

	integer n
	integer a1(n)
	integer a2(n)
	character*(*) text

	integer i

	do i=1,n
	  if( a1(i) .ne. a2(i) ) then
	    write(6,*) text,' : arrays differ'
	    write(6,*) n,i,a1(i),a2(i)
	    stop 'error stop array_check: arrays differ'
	  end if
	end do

	end

c******************************************************
c******************************************************
c******************************************************
c fem files
c******************************************************
c******************************************************
c******************************************************

	subroutine femini

c initializes internal data structure for FEM files

	implicit none

	include 'supout.h'

	integer icall
	save icall
	data icall /0/

	if( bdebug_out ) then
	  write(6,*) 'debug_out: femini'
	  write(6,*) icall,nunit_fvl
	end if

	if( icall .ne. 0 ) return

	icall = 1

	nunit_fem = 0
	iformat = 0

	end

c******************************************************

	subroutine femclose

c closes FEM file

	implicit none

	include 'supout.h'

	call femini
	if( nunit_fem .gt. 0 ) close(nunit_fem)
	nunit_fem = 0

	end

c******************************************************

	subroutine femopen(type)

c opens FEM file and reads header

	use mod_depth
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	character*(*) type

	include 'supout.h'
	include 'simul.h'

	character*80 file

	logical bformat,bregdata
	integer nvers,np,it,lmax,ntype
	integer nknaux,nelaux,nlvaux,nvar
	integer nknddi,nelddi,nlvddi
	integer ierr,l
	integer datetime(2)
	real regpar(7)
	double precision dtime
	integer ifemop,fem_file_regular

c initialize routines

	call femini

c open file

	np = nkn
	call def_make(type,file)
	call fem_file_read_open(file,np,iformat,nunit)

	if( nunit .le. 0 ) then
                write(6,*) file
		stop 'error stop femopen: cannot open FEM file'
        else
                write(6,*) 'File opened :'
                inquire(nunit,name=file)
                write(6,*) file
                write(6,*) 'Reading FEM file ...'
	end if

c read first header

        call fem_file_read_params(iformat,nunit,dtime
     +                          ,nvers,np,lmax,nvar,ntype,datetime,ierr)

	if( ierr .ne. 0 ) then
		write(6,*) 'ierr = ',ierr
		stop 'error stop femopen: error reading header'
	end if

	bregdata = fem_file_regular(ntype) > 0

        write(6,*)
        write(6,*) ' nvers = ', nvers
        write(6,*) '   nkn = ',np,   ' ntype = ',ntype
        write(6,*) '   nlv = ',lmax, '  nvar = ',nvar
        write(6,*)
	if( bregdata ) then
          write(6,*) 'the file is a regular grid'
          write(6,*)
	end if

	if( .not. bregdata .and. nkn .ne. np ) goto 99
	!if( nkn .lt. np ) goto 99
	if( bregdata .and. lmax > 1 ) goto 98

        nlv = lmax
	call allocate_simulation(nel)
	!if( bregdata ) call reallocate_2d_arrays(np) !re-allocate with minimum np
	call get_dimension_post(nknddi,nelddi,nlvddi)
	if( nlvddi .lt. lmax ) goto 99

	it = nint(dtime)

c read second header

	call fem_file_read_2header(iformat,nunit,ntype,nlv
     +				,hlv,regpar,ierr)

	if( ierr .ne. 0 ) then
		write(6,*) 'ierr = ',ierr
		stop 'error stop femopen: error reading hlv'
	end if

	call level_k2e_sh
	call init_sigma_info(nlv,hlv)		!sets up hlv

	write(6,*) 'hlv: ',nlv
	write(6,*) (hlv(l),l=1,nlv)

c rewind for a clean state

	nunit_fem = nunit
	rewind(nunit)

	return
   98	continue
	write(6,*) 'error in parameters : regular - lmax'
	write(6,*) 'the file is regular and 3D'
	write(6,*) 'Cannot handle interpolation yet'
	stop 'error stop femopen'
   99	continue
	write(6,*) 'error in parameters : basin - simulation'
	write(6,*) 'nkn : ',nkn,np
	write(6,*) 'nlv : ',nlvdi,lmax
	write(6,*) 'parameters are different between basin and simulation'
	stop 'error stop femopen'
	end

c******************************************************

	function femnext(atime,ivar,nlvddi,nkn,array)

c reads next FEM record - is true if a record has been read, false if EOF

	use levels

	implicit none

	logical femnext			!true if record read, flase if EOF
	integer it			!time of record
	double precision atime		!absolute time
	integer ivar			!type of variable
	integer nlvddi			!dimension of vertical coordinate
	integer nkn			!number of points needed
	real array(nlvddi,nkn,1)	!values for variable

	include 'supout.h'

	logical bfound,bformat,bregdata,bwind,bvel
	integer ierr
	integer i,iv,ip
	integer nvers,np,lmax,nvar,ntype
	integer datetime(2)
	real regpar(7)
	double precision dtime
	real, allocatable :: p3read(:,:,:)
	real, allocatable :: v1v(:)
	real fact
	real vmin,vmax
	character*80 string

	integer fem_file_regular

	if( nlvddi .ne. nlvdi ) stop 'error stop femnext: nlvddi'

	call femini
	nunit = nunit_fem
	bformat = iformat .eq. 1

	!write(6,*) 'femnext format: ',bformat,iformat

	np = 0
        call fem_file_read_params(iformat,nunit,dtime
     +                          ,nvers,np,lmax,nvar,ntype,datetime,ierr)
	if( ierr .ne. 0 ) goto 7
	nlv = lmax
	regpar = 0
	call fem_file_read_2header(iformat,nunit,ntype,nlv
     +				,hlv,regpar,ierr)
	if( ierr .ne. 0 ) goto 7

	bregdata = fem_file_regular(ntype) > 0
	if( .not. bregdata .and. np .ne. nkn ) goto 99
	if( bregdata .and. lmax > 1 ) goto 98

	write(6,*) 'ggggggguuuuu bregdata: ',bregdata,np

	if( bregdata ) then
	  write(6,*) 'plotting regular grid...'
	  !write(6,*) 'not yet ready for regular grid...'
	  !stop
	  allocate(p3read(nlvddi,np,2),v1v(np))
	else
	  allocate(p3read(nlvddi,nkn,2),v1v(nkn))
	end if
	write(6,*) 'p3read allocated: ',nlvddi,nkn,np,2
	p3read = 0.

	!it = nint(dtime)
	call ptime_set_date_time(datetime(1),datetime(2))
	call ptime_set_dtime(dtime)
	call ptime_get_atime(atime)

	ip = 1
	bfound = .false.
	do i=1,nvar
	  if( bfound ) then
            call fem_file_skip_data(iformat,nunit
     +                          ,nvers,np,lmax
     +                          ,string,ierr)
	    if( ierr .ne. 0 ) goto 98
	  else
	    write(6,*) 'reading data: ',bregdata,nlvddi,nkn,np,ip
            call fem_file_read_data(iformat,nunit
     +                          ,nvers,np,lmax
     +				,string
     +                          ,ilhkv,v1v
     +                          ,nlvddi,p3read(1,1,ip),ierr)
	    if( ierr .ne. 0 ) goto 98
	    call string2ivar(string,iv)
	    bfound = iv .eq. ivar
	    bwind = iv .eq. 21
	    bvel = iv .eq. 2
	    if( bfound .and. ( bwind .or. bvel ) ) then
	      ip = ip + 1
	      if( ip .eq. 2 ) bfound = .false.
	    end if
	  end if
	end do

	if( ip == 3 ) then
	  v1v = sqrt( p3read(1,:,1)**2 + p3read(1,:,2)**2 )
	  vmin = minval(v1v)
	  vmax = maxval(v1v)
	  write(6,*) 'speed: ',vmin,vmax
	end if

	if( bregdata ) then
	  write(6,*) 'interpolating from regular grid... ',ip
	  ip = min(2,ip)
	  call fem_interpolate(nlvddi,nkn,np,ip,regpar,ilhkv
     +					,p3read,array)
	  regp = regpar		!save for later
	else
	  array = p3read
	end if

	deallocate(p3read,v1v)

	if( bfound ) then
	  call level_k2e_sh
	else
	  ivar = 0
	end if

c set return value

    7	continue

	if( ierr .gt. 0 ) then
		!stop 'error stop femnext: error reading data record'
		write(6,*) '*** femnext: error reading data record'
		femnext = .false.
	else if( ierr .lt. 0 ) then
		femnext = .false.
	else
		femnext = .true.
	end if

c end

	return
   98	continue
	write(6,*) 'error in parameters : regular - lmax'
	write(6,*) 'the file is regular and 3D'
	write(6,*) 'Cannot handle interpolation yet'
	stop 'error stop femopen'
   99	continue
	write(6,*) 'nkn,np: ',nkn,np
	stop 'error stop femnext: np different from nkn'
	end

c******************************************************

	subroutine femscale(nlvddi,nkn,fact,array)

	implicit none

	integer nlvddi,nkn
	real fact
	real array(nlvddi,nkn)

	integer k,l

	do k=1,nkn
	  do l=1,nlvddi
	    array(l,k) = fact * array(l,k)
	  end do
	end do

	end

c******************************************************

	subroutine fem_interpolate(nlvddi,nkn,np,ip,regpar,ilhkv
     +			,p3reg,array)

c interpolates from a regular grid (only for 2D)

	implicit none

	integer nlvddi,nkn,np,ip
	real regpar(7)
	integer ilhkv(nkn)
	real p3reg(nlvddi,np,ip)
	real array(nlvddi,nkn,ip)

	integer i,ivar,j
	integer nx,ny
	real x0,y0,dx,dy,flag
	real areg(np)
	real afem(nkn)

	nx = nint(regpar(1))
	ny = nint(regpar(2))
	x0 = regpar(3)
	y0 = regpar(4)
	dx = regpar(5)
	dy = regpar(6)
	flag = regpar(7)
	ilhkv = 1		!works only for 2D

	call setgeo(x0,y0,dx,dy,flag)
	
	do ivar=1,ip
	  do i=1,np
	    areg(i) = p3reg(1,i,ivar)
	  end do
	  !write(6,*) 'reg: ',(areg(j),j=1,np,np/10)

	  call am2av(areg,afem,nx,ny)

	  do i=1,nkn
	    array(1,i,ivar) = afem(i)
	  end do
	  !write(6,*) 'fem: ',(afem(j),j=1,nkn,nkn/10)
	end do

	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine allocate_simulation(npd)

	use mod_hydro_plot
	use mod_hydro_print
	use mod_hydro_vel
	use mod_hydro
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer npd

	integer np
	real flag

	nlvdi = nlv
	np = max(2*nel,npd)

	call levels_init(nkn,nel,nlvdi)
	call mod_hydro_init(nkn,nel,nlvdi)
	call mod_hydro_vel_init(nkn,nel,nlvdi)
	call mod_hydro_print_init(nkn,nlvdi)
	call mod_hydro_plot_init(nkn,nel,nlvdi,np)

	call mkareafvl			!area of finite volumes

        call get_flag(flag)
	p3 = flag

	write(6,*) 'allocate_simulation: ',nkn,nel,nlvdi,np

	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine get_dimension_post(nknddi,nelddi,nlvddi)

	implicit none

	integer nknddi,nelddi,nlvddi

	call get_dimension_post_2d(nknddi,nelddi)
	call get_dimension_post_3d(nlvddi)

	end

c*****************************************************************

	subroutine get_dimension_post_2d(nknddi,nelddi)

	use basin

	implicit none

	integer nknddi,nelddi

	call basin_get_dimension(nknddi,nelddi)

	end

c*****************************************************************

	subroutine get_dimension_post_3d(nlvddi)

	use levels

	implicit none

	integer nlvddi

	call levels_get_dimension(nlvddi)

	end

c*****************************************************************

