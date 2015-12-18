!
! utility routines for shyelab: elabutil
!
! revision log :
!
! 15.07.2015	ggu	written from scratch
! 22.09.2015	ggu	new routine open_shy_file()
! 10.10.2015	ggu	code added to handle FLX routines
!
!************************************************************

!====================================================
	module elabutil
!====================================================

	implicit none

	logical, save, private :: binitialized = .false.
	double precision, parameter :: flag = -999.

	logical, save :: bout
	logical, save :: baverbas
	logical, save :: baver
	logical, save :: baverdir
	logical, save :: bsum
	logical, save :: bmin
	logical, save :: bmax
	logical, save :: bstd
	logical, save :: brms
	logical, save :: bsumvar
	logical, save :: bsplit
	logical, save :: b2d

	logical, save :: bmem
	logical, save :: bask
	logical, save :: bverb
	logical, save :: bwrite
	logical, save :: bquiet
	logical, save :: bdate

	integer, save :: ifreq
	integer, save :: tmin
	integer, save :: tmax

	logical, save :: bnode
	logical, save :: bnodes
	logical, save :: boutput
	logical, save :: bneedbasin
	logical, save :: btrans

	logical, save :: bopen

	logical, save :: btmin
	logical, save :: btmax
	logical, save :: binclusive
	double precision, save :: atmin
	double precision, save :: atmax

	logical, save :: bthreshold
	double precision, save :: threshold

	integer, save :: nodesp
	integer, save :: nnodes = 0
	integer, save, allocatable :: nodes(:)

	real, save :: fact

	integer, save :: istep
	integer, save :: mode
	integer, save :: modeb

	integer, save :: date = 0
	integer, save :: time = 0
	integer, save :: datetime(2)

        character*80, save :: infile
        character*80, save :: stmin,stmax
        character*80, save :: nodefile
        character*10, save :: outformat

        INTERFACE elabutil_check_time
        MODULE PROCEDURE elabutil_check_time_i,elabutil_check_time_d
        END INTERFACE

!====================================================
	contains
!====================================================

	subroutine elabutil_init(type)

	use clo

	character*(*) type

	call elabutil_set_options(type)
	call clo_parse_options
	call elabutil_get_options(type)

	binitialized = .true.

	end subroutine elabutil_init

!************************************************************

	subroutine elabutil_set_options(type)

	use clo

	character*(*) type

	if( binitialized ) return

	if( type == 'SHY' ) then
          call clo_init('shyelab','shy-file','3.0')
	else if( type == 'NOS' ) then
          call clo_init('noselab','nos-file','3.0')
	else if( type == 'OUS' ) then
          call clo_init('ouselab','ous-file','3.0')
	else if( type == 'EXT' ) then
          call clo_init('extelab','ext-file','3.0')
	else if( type == 'FLX' ) then
          call clo_init('extelab','flx-file','3.0')
	else
	  write(6,*) 'type : ',trim(type)
	  stop 'error stop elabutil_set_options: unknown type'
	end if

        call clo_add_info('returns info on or elaborates a shy file')

        call clo_add_sep('what to do (only one of these may be given)')

        call clo_add_option('out',.false.,'writes new nos file')
        call clo_add_option('averbas',.false.,'average over basin')
        call clo_add_option('aver',.false.,'average over records')
        call clo_add_option('averdir',.false.,'average for directions')
        call clo_add_option('sum',.false.,'sum over records')
        call clo_add_option('min',.false.,'minimum of records')
        call clo_add_option('max',.false.,'maximum of records')
	call clo_add_option('std',.false.,'standard deviation of records')
        call clo_add_option('rms',.false.,'root mean square of records')
        call clo_add_option('sumvar',.false.,'sum over variables')
        call clo_add_option('split',.false.,'split file for variables')
	call clo_add_option('2d',.false.,'average vertically to 2d field')

	call clo_add_option('threshold t',flag,'records over threshold t')
	call clo_add_option('fact fact',1.,'multiply values by fact')

        call clo_add_sep('options in/output')

        !call clo_add_option('basin name',' ','name of basin to be used')
	call clo_add_option('mem',.false.,'if no file given use memory')
        call clo_add_option('ask',.false.,'ask for simulation')
        call clo_add_option('verb',.false.,'be more verbose')
        call clo_add_option('write',.false.,'write min/max of values')
        call clo_add_option('quiet',.false.,'do not be verbose')

        call clo_add_sep('additional options')

        call clo_add_option('node n',0,'extract vars of node number n')
        call clo_add_option('nodes file',' '
     +				,'extract vars at nodes given in file')
	call clo_add_option('freq n',0.,'frequency for aver/sum/min/max')
        call clo_add_option('tmin time',' '
     +                          ,'only process starting from time')
        call clo_add_option('tmax time',' '
     +                          ,'only process up to time')
        call clo_add_option('inclusive',.false.,'include time period')

        call clo_add_option('outformat form','native','output format')

	end subroutine elabutil_set_options

!************************************************************

	subroutine elabutil_get_options(type)

	use clo

	character*(*) type

	if( binitialized ) return

        call clo_get_option('out',bout)
        call clo_get_option('averbas',baverbas)
        call clo_get_option('aver',baver)
        call clo_get_option('averdir',baverdir)
        call clo_get_option('sum',bsum)
        call clo_get_option('min',bmin)
        call clo_get_option('max',bmax)
        call clo_get_option('std',bstd)
        call clo_get_option('rms',brms)
        call clo_get_option('sumvar',bsumvar)
        call clo_get_option('split',bsplit)
        call clo_get_option('2d',b2d)

        call clo_get_option('threshold',threshold)
        call clo_get_option('fact',fact)

        call clo_get_option('node',nodesp)
        call clo_get_option('nodes',nodefile)

        call clo_get_option('mem',bmem)
        call clo_get_option('ask',bask)
        call clo_get_option('verb',bverb)
        call clo_get_option('write',bwrite)
        call clo_get_option('quiet',bquiet)

        call clo_get_option('freq',ifreq)
        call clo_get_option('tmin',stmin)
        call clo_get_option('tmax',stmax)
        call clo_get_option('inclusive',binclusive)

        call clo_get_option('outformat',outformat)

        if( .not. bask .and. .not. bmem ) call clo_check_files(1)
        call clo_get_file(1,infile)
        call ap_set_names(' ',infile)

        if( .not. bquiet ) then
	  if( type == 'SHY' ) then
            call shyfem_copyright('shyelab - Elaborate SHY files')
	  else if( type == 'NOS' ) then
            call shyfem_copyright('noselab - Elaborate NOS files')
	  else if( type == 'OUS' ) then
            call shyfem_copyright('ouselab - Elaborate OUS files')
	  else if( type == 'EXT' ) then
            call shyfem_copyright('extelab - Elaborate EXT files')
	  else if( type == 'FLX' ) then
            call shyfem_copyright('flxelab - Elaborate FLX files')
	  else
	    write(6,*) 'type : ',trim(type)
	    stop 'error stop elabutil_get_options: unknown type'
	  end if
        end if

        bnode = nodesp > 0
        bnodes = nodefile .ne. ' '

        boutput = bout .or. bsplit .or. b2d
	boutput = boutput .or. outformat /= 'native'
        !btrans is added later

        bneedbasin = b2d .or. baverbas .or. bnode .or. bnodes
	bneedbasin = bneedbasin .or. outformat == 'gis'
	bneedbasin = bneedbasin .or. ( type == 'OUS' .and. bsplit )

        modeb = 2
        if( bneedbasin ) modeb = 3

	bthreshold = ( threshold /= flag )

	end subroutine elabutil_get_options

!************************************************************
!************************************************************
!************************************************************

	subroutine elabutil_date_and_time

        bdate = date .gt. 0
        if( bdate ) call dtsini(date,time)
        datetime(1) = date
        datetime(2) = time

        atmin = 0.
        atmax = 0.
        btmin = stmin .ne. ' '
        btmax = stmax .ne. ' '
        if( btmin ) call fem_file_string2time(stmin,atmin)
        if( btmax ) call fem_file_string2time(stmax,atmax)

        if( bverb ) then
          write(6,*) 'time limits: '
          write(6,*) stmin(1:len_trim(stmin)),btmin,atmin
          write(6,*) stmax(1:len_trim(stmax)),btmax,atmax
        end if

	end subroutine elabutil_date_and_time

!************************************************************

	function elabutil_check_time_i(it,itnew,itold)

! integer version

	logical elabutil_check_time_i
	integer it,itnew,itold

	double precision dtime,dtimenew,dtimeold

        dtime = it
	dtimenew = itnew
	dtimeold = itold

	elabutil_check_time_i = 
     +		elabutil_check_time_d(dtime,dtimenew,dtimeold)

	end function elabutil_check_time_i

!************************************************************

	function elabutil_check_time_d(dtime,dtimenew,dtimeold)

! double version (relativ)

	logical elabutil_check_time_d
	double precision dtime,dtimenew,dtimeold

	double precision atime,atimenew,atimeold

        call fem_file_convert_time(datetime,dtime,atime)
        call fem_file_convert_time(datetime,dtimenew,atimenew)
        call fem_file_convert_time(datetime,dtimeold,atimeold)

	elabutil_check_time_d =
     +		elabutil_check_time_a(atime,atimenew,atimeold)

	end function elabutil_check_time_d

!************************************************************

	function elabutil_check_time_a(atime,atimenew,atimeold)

! double version (absolute)

	logical elabutil_check_time_a
	double precision atime,atimenew,atimeold

	logical bdebug
	logical btimew

	bdebug = .true.
	bdebug = .false.
        btimew = .true.

        if( btmin ) btimew = btimew .and. atime >= atmin
        if( btmax ) btimew = btimew .and. atime <= atmax

	elabutil_check_time_a = btimew

	if( bdebug ) then
	  write(6,*) 'exclusive..........',btimew,binclusive
	  write(6,*) 'exclusive..........',atmin,atime,atmax
	end if

	if( .not. binclusive ) return

        if( btmin ) then
	  btimew = btimew .or. (atime < atmin .and. atmin < atimenew)
	end if
        if( btmax ) then
	  btimew = btimew .or. (atimeold < atmax .and. atmax < atime)
	end if

	elabutil_check_time_a = btimew

	end function elabutil_check_time_a

!************************************************************
!************************************************************
!************************************************************

	subroutine elabutil_set_averaging(nvar)

	integer nvar

        mode = 0
        btrans = .false.

        if( baver ) mode = 1
        if( bsum )  mode = 2
        if( bmin )  mode = 3
        if( bmax )  mode = 4
        if( bsumvar )  mode = 2
        if( bstd )  mode = 5
        if( brms )  mode = 6
        if( bthreshold )  mode = 7
        if( baverdir ) mode = 8

        if( mode > 0 ) then     !prepare for averaging
          btrans = .true.
          if( bsumvar ) then            !sum over variables
            istep = 1
            ifreq = nvar
          else if( nvar > 1 ) then
            write(6,*) 'file contains different variables: ',nvar
	    stop 'error stop noselab: averaging only with one variable'
          else if( ifreq .ge. 0 ) then  !normal averaging
            istep = 1
          else                          !accumulate every -ifreq record
            istep = -ifreq
            ifreq = 0
          end if
	end if

	end subroutine elabutil_set_averaging

c***************************************************************

	subroutine handle_nodes

	integer i

          if( bnodes ) then
            nnodes = 0
            call get_node_list(nodefile,nnodes,nodes)
            allocate(nodes(nnodes))
            call get_node_list(nodefile,nnodes,nodes)
          else if( bnode ) then
            nnodes = 1
            allocate(nodes(nnodes))
            nodes(1) = nodesp
          end if

          write(6,*) 'nodes: ',nnodes,(nodes(i),i=1,nnodes)
          call convert_internal_nodes(nnodes,nodes)

          if( bnode ) bnodes = .true.

	end subroutine handle_nodes

c***************************************************************

	subroutine write_nodes(dtime,ivar,nlvddi,cv3)

	double precision dtime
	integer ivar
	integer nlvddi
	real cv3(nlvddi,*)

	integer j,node,it

        do j=1,nnodes
          node = nodes(j)
          it = dtime
          call write_node(j,node,cv3,it,ivar)
        end do

	end subroutine write_nodes

!====================================================
	end module elabutil
!====================================================

c***************************************************************

        subroutine outfile_make_depth(nkn,nel,nen3v,hm3v,hev,hkv)

c averages vertically

        implicit none

        integer nkn,nel
        integer nen3v(3,nel)
        real hm3v(3,nel)
        real hev(nel)
        real hkv(nkn)

        integer k,ie,ii
        real h,hm

	hkv = -huge(1.)

        do ie=1,nel
	  hm = 0.
          do ii=1,3
            k = nen3v(ii,ie)
            h = hm3v(ii,ie)
            hkv(k) = max(hkv(k),h)
	    hm = hm + h
          end do
	  hev(ie) = hm / 3.
        end do

        end

c***************************************************************

        subroutine outfile_make_hkv(nkn,nel,nen3v,hev,hkv)

c averages vertically

        implicit none

        integer nkn,nel
        integer nen3v(3,nel)
        real hev(nel)
        real hkv(nkn)

        integer k,ie,ii
        real h

        do ie=1,nel
          h = hev(ie)
          do ii=1,3
            k = nen3v(ii,ie)
            hkv(k) = max(hkv(k),h)
          end do
        end do

        end

c***************************************************************

        subroutine depth_stats(nkn,nlvddi,ilhkv)

c       computes statistics on levels

        implicit none

        integer nkn
        integer nlvddi
        integer ilhkv(nkn)

        integer count(nlvddi)
        integer ccount(nlvddi)

        integer nlv,lmax,l,k,nc,ll

        nlv = 0
        do l=1,nlvddi
          count(l) = 0
          ccount(l) = 0
        end do

        do k=1,nkn
          lmax = ilhkv(k)
          if( lmax .gt. nlvddi ) stop 'error stop depth_stats: lmax'
          count(lmax) = count(lmax) + 1
          nlv = max(nlv,lmax)
        end do

        do l=nlv,1,-1
          nc = count(l)
          do ll=1,l
            ccount(ll) = ccount(ll) + nc
          end do
        end do

        nc = 0
        write(6,*) 'statistics for layers: ',nlv
        do l=1,nlv
          if( count(l) > 0 ) then
            write(6,*) l,count(l),ccount(l)
            nc = nc + count(l)
          end if
        end do
        write(6,*) 'total count: ',nc

        end

c***************************************************************

        subroutine convert_internal_node(node)
	integer node
	stop 'error stop convert_internal_node: not supported'
	end

        subroutine convert_internal_nodes(n,nodes)

        use basin

        implicit none

        integer n
        integer nodes(n)

        integer ne,ni,i
        integer ipint

        if( n <= 0 ) return

	do i=1,n
	  ne = nodes(i)
          if( ne <= 0 ) goto 99
          ni = ipint(ne)
          if( ni <= 0 ) goto 98
	  nodes(i) = ni
        end do

	return
   98	continue
        write(6,*) 'cannot find node: ',ne
        stop 'error stop convert_internal_nodes: no such node'
   99	continue
        write(6,*) 'cannot convert node: ',ne
        stop 'error stop convert_internal_nodes: no such node'
        end

c***************************************************************

	subroutine get_node_list(file,n,nodes)

c for n == 0 only checks how many nodes to read
c for n > 0 reads nodes into nodes() (error if n is too small)

	implicit none

	character*(*) file
	integer n
	integer nodes(n)

	integer nin,ios,node,ndim
	logical btest

	integer ifileo

	nin = ifileo(0,file,'form','old')
	if( nin .le. 0 ) goto 99

	ndim = n
	btest = ndim == 0

	n = 0
	do
	  read(nin,*,iostat=ios) node
	  if( ios > 0 ) goto 98
	  if( ios < 0 ) exit
	  if( node .le. 0 ) exit
	  n = n + 1
	  if( .not. btest ) then
	    if( n > ndim ) goto 96
	    nodes(n) = node
	  end if
	end do

	if( n == 0 ) goto 97

	close(nin)

	return
   96	continue
	write(6,*) 'n,ndim :',n,ndim
	write(6,*) 'file: ',trim(file)
	stop 'error stop get_node_list: dimension error'
   97	continue
	write(6,*) 'no data in file ',trim(file)
	stop 'error stop get_node_list: read error'
   98	continue
	write(6,*) 'read error in record ',n
	write(6,*) 'in file ',trim(file)
	stop 'error stop get_node_list: read error'
   99	continue
	write(6,*) 'file: ',trim(file)
	stop 'error stop get_node_list: cannot open file'
	end

c***************************************************************
c***************************************************************
c***************************************************************

        subroutine aver(xx,n,xaver,rnull)

c computes min/max of vector
c
c xx            vector
c n             dimension of vector
c xmin,xmax     min/max value in vector
c rnull         invalid value

        implicit none

        integer n
        real xx(n)
        real xaver,rnull

        integer i,nacu
        double precision acu

        nacu = 0
        acu = 0.
        xaver = rnull

        do i=1,n
          if(xx(i).ne.rnull) then
            acu = acu + xx(i)
            nacu = nacu + 1
          end if
        end do

        if( nacu .gt. 0 ) xaver = acu / nacu

        end

c***************************************************************

        subroutine mimar(xx,n,xmin,xmax,rnull)

c computes min/max of vector (2d)
c
c xx            vector
c n             dimension of vector
c xmin,xmax     min/max value in vector
c rnull         invalid value

        implicit none

        integer n,i,nmin
        real xx(n)
        real xmin,xmax,x,rnull

        do i=1,n
          if(xx(i).ne.rnull) goto 1
        end do
    1   continue

        if(i.le.n) then
          xmax=xx(i)
          xmin=xx(i)
        else
          xmax=rnull
          xmin=rnull
        end if

        nmin=i+1

        do i=nmin,n
          x=xx(i)
          if(x.ne.rnull) then
            if(x.gt.xmax) xmax=x
            if(x.lt.xmin) xmin=x
          end if
        end do

        end

c***************************************************************

        subroutine mimar_s(xx,nlvddi,n,xmin,xmax,rnull)

c computes min/max of vector (3d)
c
c xx            vector
c n             dimension of vector
c xmin,xmax     min/max value in vector
c rnull         invalid value

        implicit none

        integer n
        integer nlvddi
        real xx(nlvddi,n)
        real xmin,xmax,rnull

        integer k,l
        real x

        do k=1,n
          do l=1,nlvddi
            x=xx(l,k)
            if(x.ne.rnull) then
              if( x .lt. xmin .or. x .gt. xmax ) then
                write(6,*) l,k,x
              end if
            end if
          end do
        end do

        end

c***************************************************************

        subroutine write_aver(it,ivar,cmin,cmax,cmed,vtot)

c writes basin average to file

        implicit none

        integer it,ivar
        real cmin,cmax,cmed,vtot

        real totmass

        totmass = cmed * vtot

        write(6,1234) it,ivar,cmin,cmed,cmax,totmass
        write(100+ivar,1235) it,cmin,cmed,cmax,totmass
        write(100,'(i10,e14.6)') it,vtot

 1234   format(2i10,3f12.4,e14.6)
 1235   format(i10,3f12.4,e14.6)
        end

c***************************************************************

        subroutine write_node(i,node,cv3,it,ivar)

	use levels
	use mod_depth

        implicit none

        integer i
        integer node
        real cv3(nlvdi,*)
        integer it
        integer ivar

        integer ki,ke
        integer l,lmax
	real z,h
	character*40 format
	real hl(nlvdi)

	integer ipext

	ki = node
	ke = ipext(ki)
	lmax = ilhkv(ki)

        write(4,*) it,i,ke,ki,lmax,ivar
        write(4,*) (cv3(l,ki),l=1,lmax)

        write(3,*) it,i,ke,ki,lmax,ivar
        write(format,'(a,i4,a)') '(',lmax,'(f12.4))'
        write(3,format) (cv3(l,ki),l=1,lmax)

        z = 0.
        h = hkv(ki)
        call write_profile_c(it,i,ki,ke,lmax,ivar,h,z
     +				,cv3(1,ki),hlv,hl)

        end

c***************************************************************

        subroutine write_node_2d(it,nvar,cv)

        implicit none

        integer it
        integer nvar
        real cv(nvar)

        integer i

        write(89,*) it,(cv(i),i=1,nvar)

        end

c***************************************************************

        subroutine write_profile_c(it,i,ki,ke,lmax,ivar,h,z,c,hlv,hl)

        implicit none

        integer it,i,ki,ke
        integer lmax
        integer ivar
        real z,h
        real c(lmax)
        real hlv(lmax)
        real hl(lmax)

        logical bcenter
        integer l
        integer nlvaux,nsigma
        real hsigma
        real uv

        bcenter = .true.        !depth at center of layer ?

        call get_sigma_info(nlvaux,nsigma,hsigma)

        call get_layer_thickness(lmax,nsigma,hsigma,z,h,hlv,hl)
        call get_bottom_of_layer(bcenter,lmax,z,hl,hl)  !orig hl is overwritten

        write(2,*) it,i,ke,ki,lmax,ivar
        do l=1,lmax
          write(2,*) hl(l),c(l)
        end do

        end

c***************************************************************

	subroutine ilhe2k(nkn,nel,nen3v,ilhv,ilhkv)

	implicit none

	integer nkn,nel
	integer nen3v(3,nel)
	integer ilhv(nel)
	integer ilhkv(nkn)

	integer ie,ii,k,l

	ilhkv = 0

	do ie=1,nel
	  l = ilhv(ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    ilhkv(k) = max(l,ilhkv(k))
	  end do
	end do

	end

c***************************************************************

        subroutine open_shy_file(file,status,nunit)

c open SHY file
c
c nunit is 0 if no other file exists

	use clo

        implicit none

        character*(*) status
        integer nunit

        character*80 file
        integer ifileo

        nunit = 0
        if( file == ' ' ) return

        nunit = ifileo(0,file,'unform',status)

        if( nunit .le. 0 ) then
          write(6,*) 'file: ',trim(file)
          stop 'error stop open_next_shy_file: opening file'
        end if

        end

c***************************************************************

