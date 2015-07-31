
!====================================================
	module elabutil
!====================================================

	implicit none

	logical, save, private :: binitialized = .false.

	logical, save :: bout
	logical, save :: baverbas
	logical, save :: baver
	logical, save :: bsum
	logical, save :: bmin
	logical, save :: bmax
	logical, save :: bsumvar
	logical, save :: bsplit
	logical, save :: b2d

	logical, save :: bmem
	logical, save :: bask
	logical, save :: bverb
	logical, save :: bwrite
	logical, save :: bquiet
	logical, save :: bdate

	integer, save :: nodesp
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

	integer, save :: istep
	integer, save :: mode
	integer, save :: modeb

	integer, save :: date,time,datetime(2)

        character*80, save :: infile
        character*80, save :: stmin,stmax
        character*80, save :: nodefile

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
	else
	  write(6,*) 'type : ',trim(type)
	  stop 'error stop elabutil_set_options: unknown type'
	end if

        call clo_add_info('returns info on or elaborates a nos file')

        call clo_add_sep('what to do (only one of these may be given)')

        call clo_add_option('out',.false.,'writes new nos file')
        call clo_add_option('averbas',.false.,'average over basin')
        call clo_add_option('aver',.false.,'average over records')
        call clo_add_option('sum',.false.,'sum over records')
        call clo_add_option('min',.false.,'minimum of records')
        call clo_add_option('max',.false.,'maximum of records')
        call clo_add_option('sumvar',.false.,'sum over variables')
        call clo_add_option('split',.false.,'split file for variables')
	call clo_add_option('2d',.false.,'average vertically to 2d field')

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

	end subroutine elabutil_set_options

!************************************************************

	subroutine elabutil_get_options(type)

	use clo

	character*(*) type

	if( binitialized ) return

        call clo_get_option('out',bout)
        call clo_get_option('averbas',baverbas)
        call clo_get_option('aver',baver)
        call clo_get_option('sum',bsum)
        call clo_get_option('min',bmin)
        call clo_get_option('max',bmax)
        call clo_get_option('sumvar',bsumvar)
        call clo_get_option('split',bsplit)
        call clo_get_option('2d',b2d)

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
	  else
	    write(6,*) 'type : ',trim(type)
	    stop 'error stop elabutil_get_options: unknown type'
	  end if
        end if

        bnode = nodesp > 0
        bnodes = nodefile .ne. ' '

        boutput = bout .or. bsplit .or. b2d
        !btrans is added later

        bneedbasin = b2d .or. baverbas .or. bnode .or. bnodes

        modeb = 2
        if( bneedbasin ) modeb = 3

	end subroutine elabutil_get_options

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

	function elabutil_check_time(it,itnew,itold)

	logical elabutil_check_time
	integer it,itnew,itold

	logical btimew
	double precision dtime,atime
	double precision dtimenew,atimenew
	double precision dtimeold,atimeold

        btimew = .true.

        dtime = it
        call fem_file_convert_time(datetime,dtime,atime)

        if( btmin ) btimew = btimew .and. atime >= atmin
        if( btmax ) btimew = btimew .and. atime <= atmax

	elabutil_check_time = btimew

	if( .not. binclusive ) return

	!write(6,*) 'inclusive..........',btimew,it,itnew

        if( btmin ) then
          dtimenew = itnew
          call fem_file_convert_time(datetime,dtimenew,atimenew)
	  !write(6,*) 'checking min: ',atime,atmin,atimenew
	  btimew = btimew .or. (atime < atmin .and. atmin < atimenew)
	end if
        if( btmax ) then
          dtimeold = itold
          call fem_file_convert_time(datetime,dtimeold,atimeold)
	  btimew = btimew .or. (atimeold < atmax .and. atmax < atime)
	end if

	elabutil_check_time = btimew

	end function elabutil_check_time

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

!====================================================
	end module elabutil
!====================================================

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

        include 'param.h'

        integer nkn
        integer nlvddi
        integer ilhkv(1)

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

        include 'param.h'

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
        real c(1)
        real hlv(1)
        real hl(1)

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

