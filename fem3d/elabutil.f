!
! utility routines for shyelab: elabutil
!
! revision log :
!
! 15.07.2015	ggu	written from scratch
! 22.09.2015	ggu	new routine open_shy_file()
! 10.10.2015	ggu	code added to handle FLX routines
! 22.02.2016	ggu	handle catmode
! 15.04.2016	ggu	handle gis files with substitution of colon
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

	logical, save :: bmem		= .false.
	logical, save :: bask		= .false.
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

	!logical, save :: btmin
	!logical, save :: btmax
	logical, save :: binclusive
	!double precision, save :: atmin
	!double precision, save :: atmax

	logical, save :: bthreshold
	double precision, save :: threshold

	integer, save :: nodesp
	integer, save :: nnodes = 0
	integer, save, allocatable :: nodes(:)
	integer, save, allocatable :: nodese(:)

	real, save :: fact			= 1

	integer, save :: istep
	integer, save :: mode
	integer, save :: modeb

	!integer, save :: date = 0
	!integer, save :: time = 0
	!integer, save :: datetime(2) = 0

	integer, save :: catmode = 0

        character*80, save :: infile		= ' '
        character*80, save :: stmin		= ' '
        character*80, save :: stmax		= ' '
        character*80, save :: nodefile		= ' '
        character*10, save :: outformat		= ' '

!====================================================
	contains
!====================================================

	subroutine elabutil_init(type,what)

	use clo

	character*(*) type
	character*(*), optional :: what

	character*80 program

	program = 'shyelab'
	if( present(what) ) program = what

	call elabutil_set_options(type,program)
	call clo_parse_options
	call elabutil_get_options(type,program)

	binitialized = .true.

	end subroutine elabutil_init

!************************************************************

	subroutine elabutil_set_options(type,program)

	use clo

	character*(*) type
	character*(*) program

	if( binitialized ) return

	if( type == 'SHY' ) then
          call clo_init(program,'shy-file','3.0')
	else if( type == 'NOS' ) then
          call clo_init(program,'nos-file','3.0')
	else if( type == 'OUS' ) then
          call clo_init(program,'ous-file','3.0')
	else if( type == 'EXT' ) then
          call clo_init(program,'ext-file','3.0')
	else if( type == 'FLX' ) then
          call clo_init(program,'flx-file','3.0')
	else
	  write(6,*) 'type : ',trim(type)
	  stop 'error stop elabutil_set_options: unknown type'
	end if

        call clo_add_info('returns info on or elaborates a shy file')

        call clo_add_sep('what to do (only one of these may be given)')

        call clo_add_option('out',.false.,'writes new shy file')
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
	!call clo_add_option('mem',.false.,'if no file given use memory')
        !call clo_add_option('ask',.false.,'ask for simulation')
        call clo_add_option('verb',.false.,'be more verbose')
        call clo_add_option('write',.false.,'write min/max of values')
        call clo_add_option('quiet',.false.,'do not be verbose')

        call clo_add_sep('additional options')

        call clo_add_option('node n',0,'extract vars of node number n')
        call clo_add_option('nodes file',' '
     +				,'extract vars at nodes given in file')
	call clo_add_option('freq n',0.,'frequency for aver/sum/min/max')
        call clo_add_option('tmin time',' '
     +                  ,'only process starting from time')
        call clo_add_option('tmax time',' '
     +                  ,'only process up to time')
        call clo_add_option('inclusive',.false.
     +			,'output includes whole time period given')

        call clo_add_option('outformat form','native','output format')

        call clo_add_option('catmode cmode',0.,'concatenation mode')

	end subroutine elabutil_set_options

!************************************************************

	subroutine elabutil_get_options(type,program)

	use clo

	character*(*) type
	character*(*) program

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

        !call clo_get_option('mem',bmem)
        !call clo_get_option('ask',bask)
        call clo_get_option('verb',bverb)
        call clo_get_option('write',bwrite)
        call clo_get_option('quiet',bquiet)

        call clo_get_option('freq',ifreq)
        call clo_get_option('tmin',stmin)
        call clo_get_option('tmax',stmax)
        call clo_get_option('inclusive',binclusive)

        call clo_get_option('outformat',outformat)

        call clo_get_option('catmode',catmode)

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

        boutput = bout .or. b2d
	boutput = boutput .or. outformat /= 'native'
        !btrans is added later
	if( bsumvar ) boutput = .false.

        bneedbasin = b2d .or. baverbas .or. bnode .or. bnodes
	bneedbasin = bneedbasin .or. outformat == 'gis'
	bneedbasin = bneedbasin .or. ( type == 'OUS' .and. bsplit )

        modeb = 2
        if( bneedbasin ) modeb = 3

	bthreshold = ( threshold /= flag )

	end subroutine elabutil_get_options

!************************************************************

	subroutine elabutil_check_options

	integer ic

	ic = count( (/b2d,bsplit,bsumvar,btrans/) )

	if( ic > 1 ) then
	  write(6,*) 'Only one of the following options can be given:'
	  write(6,*) '-2d -split -sumvar'
	  write(6,*) '-aver -sum -min -max -std -rms'
	  write(6,*) '-threshold -averdir'
	  stop 'error stop elabutil_check_options: incompatible options'
	end if

	end subroutine elabutil_check_options

!************************************************************
!************************************************************
!************************************************************


!************************************************************
!************************************************************
!************************************************************

	subroutine elabutil_set_averaging(nvar)

	integer nvar

	integer ic

        mode = 0
        btrans = .false.

	ic = count( (/baver,bsum,bmin,bmax,bstd,brms
     +				,bthreshold,baverdir/) )

	if( ic > 1 ) then
	  write(6,*) 'Only one of the following options can be given:'
	  write(6,*) '-aver -sum -min -max -std -rms'
	  write(6,*) '-threshold -averdir'
	  stop 'error stop elabutil_set_averaging: incompatible options'
	end if

        if( baver ) mode = 1
        if( bsum )  mode = 2
        if( bmin )  mode = 3
        if( bmax )  mode = 4
        if( bstd )  mode = 5
        if( brms )  mode = 6
        if( bthreshold )  mode = 7
        if( baverdir ) mode = 8

        if( mode > 0 ) then             !prepare for averaging
          btrans = .true.
          if( bsumvar ) then            !sum over variables
	    if( ifreq /= 0 ) then
	      write(6,*) 'For option -sumvar cannot use value for -freq'
	      write(6,*) 'freq = ',ifreq
	      stop 'error stop elabutil_set_averaging: freq'
	    end if
            istep = 1
          else if( ifreq .ge. 0 ) then  !normal averaging
            istep = 1
          else                          !accumulate every -ifreq record
            istep = -ifreq
            ifreq = 0
          end if
	end if

        if( bsumvar ) then            !sum over variables
	  if( ifreq /= 0 ) then
	    write(6,*) 'For option -sumvar cannot use value for -freq'
	    write(6,*) 'freq = ',ifreq
	    stop 'error stop elabutil_set_averaging: freq'
	  end if
          istep = 1
	end if

	end subroutine elabutil_set_averaging

c***************************************************************

	subroutine handle_nodes

	integer i

          if( bnodes ) then
            nnodes = 0
            call get_node_list(nodefile,nnodes,nodes)
            allocate(nodes(nnodes))
            allocate(nodese(nnodes))
            call get_node_list(nodefile,nnodes,nodes)
          else if( bnode ) then
            nnodes = 1
            allocate(nodes(nnodes))
            allocate(nodese(nnodes))
            nodes(1) = nodesp
          end if

	  if( nnodes <= 0 ) return
 
	  nodese = nodes
          write(6,*) 'nodes: ',nnodes,(nodes(i),i=1,nnodes)
          call convert_internal_nodes(nnodes,nodes)

          if( bnode ) bnodes = .true.

	end subroutine handle_nodes

c***************************************************************

	subroutine write_nodes(dtime,ivar,cv3)

	use levels

	double precision dtime
	integer ivar
	real cv3(nlvdi,*)

	integer j,node,it

	if( nnodes <= 0 ) return

        do j=1,nnodes
          node = nodes(j)
          it = nint(dtime)
          call write_node(j,node,cv3,it,ivar)
        end do

	end subroutine write_nodes

c***************************************************************

	subroutine write_nodes_vel(dtime,znv,uprv,vprv)

	use levels
	use mod_depth

	double precision dtime
	real znv(*)
	real uprv(nlvdi,*)
	real vprv(nlvdi,*)

	integer j,ki,ke,lmax,it,l,k
	real z,h
	real hl(nlvdi)
	real u(nlvdi),v(nlvdi)

	if( nnodes <= 0 ) return

        do j=1,nnodes
          ki = nodes(j)
          ke = nodese(j)
	  lmax = ilhkv(ki)
          it = nint(dtime)
          write(79,*) it,j,ke,ki,lmax
          write(79,*) znv(ki)
          write(79,*) (uprv(l,ki),l=1,lmax)
          write(79,*) (vprv(l,ki),l=1,lmax)

          z = 0.
          z = znv(ki)
          h = hkv(ki)
	  u = uprv(:,ki)
	  v = vprv(:,ki)
          call write_profile_uv(it,j,ki,ke,lmax,h,z
     +				,u,v,hlv,hl)
        end do

	write(80,*) it,(znv(nodes(k)),k=1,nnodes)

	end subroutine write_nodes_vel

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

        subroutine shy_write_aver(dtime,ivar,cmin,cmax,cmed,vtot)

c writes basin average to file

        implicit none

        double precision dtime
        integer ivar
        real cmin,cmax,cmed,vtot

	integer it
        real totmass

	it = nint(dtime)
        totmass = cmed * vtot

        !write(6,1234) it,ivar,cmin,cmed,cmax,totmass
        write(100+ivar,1235) it,cmin,cmed,cmax,totmass
        write(100,1236) it,vtot

        write(6,2234) dtime,ivar,cmin,cmed,cmax,totmass
        write(200+ivar,2235) dtime,cmin,cmed,cmax,totmass
        write(200,2236) dtime,vtot

	return
 1234   format(i10,i10,3f12.4,e14.6)
 1235   format(i10,3f12.4,e14.6)
 1236   format(i10,e14.6)
 2234   format(f15.2,i10,3f12.4,e14.6)
 2235   format(f15.2,3f12.4,e14.6)
 2236   format(f15.2,e14.6)
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
        write(100,1236) it,vtot

	return
 1234   format(i10,i10,3f12.4,e14.6)
 1235   format(i10,3f12.4,e14.6)
 1236   format(i10,e14.6)
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

	write(1,*) it,cv3(1,ki)

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

        subroutine write_2d_all_nodes(nnodes,nodes,cv2,it,ivar)

        implicit none

	integer nnodes
	integer nodes(nnodes)
        real cv2(*)
        integer it
        integer ivar

        integer iunit,i

	if( nnodes <= 0 ) return

	iunit = 200 + ivar

        write(iunit,1000) it,(cv2(nodes(i)),i=1,nnodes)
 1000	format(i11,30g14.6)

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

	subroutine write_node_vel

	end

c***************************************************************

        subroutine write_profile_c(it,j,ki,ke,lmax,ivar,h,z,c,hlv,hl)

        implicit none

        integer it,j,ki,ke
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

        write(2,*) it,j,ke,ki,lmax,ivar
        do l=1,lmax
          write(2,*) hl(l),c(l)
        end do

        end

c***************************************************************

        subroutine write_profile_uv(it,j,ki,ke,lmax,h,z,u,v,hlv,hl)

        implicit none

        integer it,j,ki,ke
        integer lmax
        real z,h
        real u(1)
        real v(1)
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

        write(82,*) it,j,ke,ki,lmax,z
        do l=1,lmax
          uv = sqrt( u(l)**2 + v(l)**2 )
          write(82,*) hl(l),u(l),v(l),uv
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

	subroutine ilhk2e(nkn,nel,nen3v,ilhkv,ilhv)

c create ilhv -> result is not exact and must be adjusted

	implicit none

	integer nkn,nel
	integer nen3v(3,nel)
	integer ilhkv(nkn)
	integer ilhv(nel)

	integer ie,ii,k,lmax

	ilhv = 0

	do ie=1,nel
	  lmax = 0
	  do ii=1,3
	    k = nen3v(ii,ie)
	    lmax = max(lmax,ilhkv(k))
	  end do
	  ilhv(ie) = lmax
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
          stop 'error stop open_shy_file: opening file'
        end if

        end

c***************************************************************

	function concat_cycle(it,itold,itstart,nrec)

	use elabutil

c decides if with concatenation we have to use record or not

	implicit none

	logical concat_cycle
	integer it,itold,itstart
	integer nrec

	concat_cycle = .false.

        !write(66,*) 'ggu: ',it,itold,itstart,nrec

        if( catmode < 0 .and. nrec /= 1 ) then
          if( it <= itold ) then
            write(6,*) 'skipping record: ',it
            it = itold
	    concat_cycle = .true.
          end if
        else if( catmode > 0 .and. itstart /= -1 ) then
          if( it >= itstart ) then
            write(6,*) 'skipping record: ',it
	    concat_cycle = .true.
          end if
        end if

	end

c***************************************************************
c***************************************************************
c***************************************************************

        subroutine gis_write_record(nb,it,ivar,nlvddi,ilhkv,cv)

c writes one record to file nb (3D)

        use basin

        implicit none

        integer nb,it,ivar,nlvddi
        integer ilhkv(nlvddi)
        real cv(nlvddi,*)

        integer k,l,lmax
	integer nout
        real x,y
	character*80 format,name
	character*20 line,dateline
	character*3 var

	integer ifileo

	call dtsgf(it,dateline)
	call gis_subst_colon(dateline,line)
	call i2s0(ivar,var)

	name = 'extract_'//var//'_'//line//'.gis'
        nout = ifileo(60,name,'form','new')
	!write(6,*) 'writing: ',trim(name)

        write(nout,*) it,nkn,ivar,dateline

	lmax = 1

        do k=1,nkn
          if( nlvddi > 1 ) lmax = ilhkv(k)
          x = xgv(k)
          y = ygv(k)

	  write(format,'(a,i5,a)') '(i10,2g14.6,i5,',lmax,'g14.6)'
          write(nout,format) k,x,y,lmax,(cv(l,k),l=1,lmax)
        end do

	close(nout)

        end

c***************************************************************

        subroutine gis_write_hydro(it,nlvddi,ilhkv,zv,uv,vv)

c writes one record to file (3D)

        use basin

        implicit none

        integer it,nlvddi
        integer ilhkv(nlvddi)
        real zv(nkn)
        real uv(nlvddi,nkn)
        real vv(nlvddi,nkn)

        integer k,l,lmax,nn,i
	integer nout
        real x,y
	character*80 format,name
	character*20 line,dateline
	character*3 var

	integer ifileo

	call dtsgf(it,dateline)
	call gis_subst_colon(dateline,line)

	name = 'extract_hydro_'//line//'.gis'
        nout = ifileo(60,name,'form','new')
	!write(6,*) 'writing: ',trim(name)

        write(nout,*) it,nkn,0,dateline

	lmax = 1

        do k=1,nkn
          if( nlvddi > 1 ) lmax = ilhkv(k)
          x = xgv(k)
          y = ygv(k)

	  nn = 1 + 2*lmax
	  write(format,'(a,i5,a)') '(i10,2g14.6,i5,',nn,'g14.6)'
          write(nout,format) k,x,y,lmax,zv(k)
     +			,(uv(l,k),vv(l,k),l=1,lmax)
        end do

	close(nout)

        end

c***************************************************************

	subroutine gis_subst_colon(line_old,line_new)

	implicit none

	character*(*) line_old,line_new

	integer n,i

	n = min(len(line_old),len(line_new))
	line_new = line_old

	do i=1,n
	  if( line_new(i:i) == ':' ) line_new(i:i) = '_'
	  if( line_new(i:i) == ' ' ) line_new(i:i) = '_'
	end do

	end

c***************************************************************

        subroutine gis_write_connect

c writes connectivity

        use basin

        implicit none

	integer ie,ii

	open(1,file='connectivity.gis',form='formatted',status='unknown')

	write(1,*) nel
	do ie=1,nel
	  write(1,*) ie,(nen3v(ii,ie),ii=1,3)
	end do

	close(1)

	end

c***************************************************************
