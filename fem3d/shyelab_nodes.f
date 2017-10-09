!
! shyelab_nodes.f: utility for extracting nodes
!
! revision log :
!
! 07.10.2017    ggu     restructured
!
c***************************************************************

	subroutine handle_nodes

	use elabutil

	implicit none

	integer i
	integer nodes_dummy(1)

          if( bnodes ) then
            nnodes = 0
            call get_node_list(nodefile,nnodes,nodes_dummy)
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
c***************************************************************
c***************************************************************

	subroutine write_nodes_scal(dtime,nvar,nndim,ivars,cv3all)

! writes scalar values for single nodes

	use shyfem_strings

	use levels
	use elabutil
	use mod_depth

	implicit none

	double precision dtime
	integer nvar,nndim
	integer ivars(nvar)
	real cv3all(nlvdi,nndim,0:nvar)

	logical b3d
	integer j,node,it,i,ivar,iv
	integer iu
	integer ki,ke,lmax,l
	real h,z,s0
	real hl(nlvdi)
	real scal(nlvdi)
	integer isubs(nvar)
	character*80 format,name
	character*10 numb,short
	character*10 shorts(nvar)
	integer, save :: icall = 0
	integer, save, allocatable :: ius(:,:,:)
	integer, save :: iuall2d = 0
	integer, save :: iuall3d = 0

	real cv3(nlvdi,nndim)	!to be deleted

	if( nnodes <= 0 ) return
	if( .not. bcompat ) return

	iu = 0
	b3d = nlv > 1

!-----------------------------------------------------------------
! open files
!-----------------------------------------------------------------

	if( icall == 0 ) then
	  
	  allocate(ius(nvar,nnodes,3))
	  ius = 0

	  do iv=1,nvar
	    call strings_get_short_name(ivars(iv),shorts(iv),isubs(iv))
	  end do

	  iu = 100
	  do j=1,nnodes
            write(numb,'(i5)') j
            numb = adjustl(numb)
	    do iv=1,nvar
	      short = shorts(iv)
	      if( isubs(iv) > 0 ) call adjust_with_isub(short,isubs(iv))
	      !write(6,*) 'short: ',short,isubs(iv)
	      call make_iunit_name(short,'','2d',j,iu)
	      ius(iv,j,2) = iu
	      if( .not. b3d ) cycle
	      call make_iunit_name(short,'','3d',j,iu)
	      ius(iv,j,3) = iu
	      call make_iunit_name(short,'_profile','3d',j,iu)
	      ius(iv,j,1) = iu
	    end do
	  end do
	  name = 'all_scal_nodes.2d.txt'
	  iuall2d = iu + 1
          !write(6,*) 'opening file : ',iu,trim(name)
          open(iuall2d,file=name,form='formatted',status='unknown')
	  if( b3d ) then
	    name = 'all_scal_nodes.3d.txt'
	    iuall3d = iu + 2
            !write(6,*) 'opening file : ',iu,trim(name)
            open(iuall3d,file=name,form='formatted',status='unknown')
	  end if
	end if
	icall = icall + 1

!-----------------------------------------------------------------
! write files
!-----------------------------------------------------------------

        it = nint(dtime)

        do j=1,nnodes
          ki = nodes(j)
          ke = nodese(j)
	  lmax = ilhkv(ki)
          h = hkv(ki)
          z = cv3all(1,ki,0)
	  format = ' '
	  if( b3d ) write(format,'(a,i3,a)') '(i12,',lmax,'f8.3)'

	  do iv=1,nvar
	    ivar = ivars(iv)
	    scal = cv3all(:,ki,iv)
	    call average_vertical_node(lmax,hlv,z,h,scal,s0)
	    iu = ius(iv,j,2)
	    write(iu,*) it,s0
            write(iuall2d,*) it,j,ke,ki,lmax,ivar
            write(iuall2d,*) s0
	    if( .not. b3d ) cycle
	    iu = ius(iv,j,3)
	    write(iu,format) it,scal(1:lmax)
	    iu = ius(iv,j,1)
            call write_profile_c(iu,it,j,ki,ke,lmax,ivar,h,z,scal,hlv)
            write(iuall3d,*) it,j,ke,ki,lmax,ivar
            write(iuall3d,*) scal(1:lmax)
	  end do
        end do

	end subroutine write_nodes_scal

c***************************************************************

	subroutine write_nodes_vel(dtime,znv,uprv,vprv)

	use levels
	use mod_depth
	use elabutil

	implicit none

	double precision dtime
	real znv(*)
	real uprv(nlvdi,*)
	real vprv(nlvdi,*)

	logical b3d
	integer i,j,ki,ke,lmax,it,l,k
	integer iu
	real z,h,u0,v0,s0,d0
	real hl(nlvdi)
	real u(nlvdi),v(nlvdi),s(nlvdi),d(nlvdi)
	integer, save :: icall = 0
	integer, save :: iuall = 0
	integer, save, allocatable :: ius(:,:,:)
	integer, parameter :: niu = 6
        character*5 :: what(niu) = (/'velx ','vely ','zeta '
     +                          ,'speed','dir  ','all  '/)
	character*5 :: numb
	character*10 short
	character*80 name
	character*80 format

	if( nnodes <= 0 ) return
	if( .not. bcompat ) return

	b3d = nlv > 1

!-----------------------------------------------------------------
! open files
!-----------------------------------------------------------------

	if( icall == 0 ) then
	  allocate(ius(niu,nnodes,3))
	  ius = 0
	  iu = 100
	  do j=1,nnodes
	    do i=1,niu
	      short = what(i)
	      call make_iunit_name(short,'','2d',j,iu)
	      ius(i,j,2) = iu
	      if( .not. b3d ) cycle
	      if( i == niu ) cycle	!do not write all.3d.*
	      call make_iunit_name(short,'','3d',j,iu)
	      ius(i,j,3) = iu
	    end do
	    if( b3d ) then
	      call make_iunit_name('vel','_profile','3d',j,iu)
	      ius(1,j,1) = iu
	    end if
	  end do
	  name = 'all_nodes.3d.txt'
	  iuall = iu + 1
          !write(6,*) 'opening file : ',iu,trim(name)
          open(iuall,file=name,form='formatted',status='unknown')
	end if
	icall = icall + 1

!-----------------------------------------------------------------
! write files
!-----------------------------------------------------------------

        it = nint(dtime)

        do j=1,nnodes
          ki = nodes(j)
          ke = nodese(j)
	  lmax = ilhkv(ki)
          write(iuall,*) it,j,ke,ki,lmax
          write(iuall,*) znv(ki)
          write(iuall,*) (uprv(l,ki),l=1,lmax)
          write(iuall,*) (vprv(l,ki),l=1,lmax)

          z = znv(ki)
          h = hkv(ki)
	  u = uprv(:,ki)
	  v = vprv(:,ki)

	  call average_vertical_node(lmax,hlv,z,h,u,u0)
	  call average_vertical_node(lmax,hlv,z,h,v,v0)
	  call c2p_ocean(u0,v0,s0,d0)

	  iu = ius(1,j,2)
	  write(iu,*) it,u0
	  iu = ius(2,j,2)
	  write(iu,*) it,v0
	  iu = ius(3,j,2)
	  write(iu,*) it,z
	  iu = ius(4,j,2)
	  write(iu,*) it,s0
	  iu = ius(5,j,2)
	  write(iu,*) it,d0
          iu = ius(6,j,2)
          write(iu,'(i12,5f12.4)') it,z,u0,v0,s0,d0

	  if( .not. b3d ) cycle

	  write(format,'(a,i3,a)') '(i12,',lmax,'f8.3)'
	  do l=1,lmax
	    call c2p_ocean(u(l),v(l),s(l),d(l))
	  end do

	  iu = ius(1,j,3)
	  write(iu,format) it,u(1:lmax)
	  iu = ius(2,j,3)
	  write(iu,format) it,v(1:lmax)
	  iu = ius(3,j,3)
	  write(iu,format) it,z
	  iu = ius(4,j,3)
	  write(iu,format) it,s(1:lmax)
	  iu = ius(5,j,3)
	  write(iu,format) it,d(1:lmax)

	  iu = ius(1,j,1)
          call write_profile_uv(iu,it,j,ki,ke,lmax,h,z,u,v,hlv)
        end do

!-----------------------------------------------------------------
! end of routine
!-----------------------------------------------------------------

	end subroutine write_nodes_vel

c***************************************************************
c***************************************************************
c***************************************************************

	subroutine write_nodes_scalar(dtime,nvar,nndim,strings,cv3all)

! writes FEM file out.fem - version for hydro

	use levels
	use mod_depth
	use elabutil
	use elabtime

	implicit none

	double precision dtime
	integer nvar,nndim
	character*(*) strings(nvar)
	real cv3all(nlvdi,nndim,0:nvar)

	integer j,iv,node
	integer iformat,lmax,np,nvers,ntype
	real regpar(7)
	real cv3(nlvdi,nnodes,nvar)
	integer il(nnodes)
	real hd(nnodes)
	character*80 file,string
	integer, save :: iunit = 0
	integer ifileo

	if( nnodes <= 0 ) return

	if( iunit == 0 ) then
          file = 'out.fem'
          iunit = ifileo(60,file,'form','unknown')
          if( iunit <= 0 ) goto 74
	end if

	do iv=1,nvar
          do j=1,nnodes
            node = nodes(j)
	    cv3(:,j,iv) = cv3all(:,node,iv)
	  end do
        end do
        do j=1,nnodes
          node = nodes(j)
	  il(j) = ilhkv(node)
	  hd(j) = hkv_max(node)
	end do

        nvers = 0
	iformat = 1
	ntype = 1
        np = nnodes
        lmax = nlv
        call fem_file_write_header(iformat,iunit,dtime
     +                          ,nvers,np,lmax
     +                          ,nvar,ntype
     +                          ,nlvdi,hlv,datetime_elab,regpar)

	do iv=1,nvar
	  string = strings(iv)
          call fem_file_write_data(iformat,iunit
     +                          ,nvers,np,lmax
     +                          ,string
     +                          ,il,hd
     +                          ,nlvdi,cv3(:,:,iv))

	end do

	return
   74	continue
        write(6,*) 'error opening file ',trim(file)
        stop 'error stop write_nodes_scal: opening file'
	end

c***************************************************************

	subroutine write_nodes_hydro(dtime,znv,uprv,vprv)

! writes FEM file out.fem - version for scalars

	use levels
	use mod_depth
	use elabutil
	use elabtime
	!use shyelab_out

	implicit none

	double precision dtime
	real znv(*)
	real uprv(nlvdi,*)
	real vprv(nlvdi,*)

	integer j,iv,node,isub
	integer iformat,lmax,np,nvers,ntype
	integer ivar,nvar
	real regpar(7)
	real z(nnodes)
	real u(nlvdi,nnodes)
	real v(nlvdi,nnodes)
	integer il(nnodes)
	real hd(nnodes)
	character*80 file,string,stringx,stringy
	integer, save :: iunit = 0
	integer ifileo

	if( nnodes <= 0 ) return

	if( iunit == 0 ) then
          file = 'out.fem'
          iunit = ifileo(60,file,'form','unknown')
          if( iunit <= 0 ) goto 74
	end if

	do iv=1,nvar
          do j=1,nnodes
            node = nodes(j)
	    z(j) = znv(node)
	    u(:,j) = uprv(:,node)
	    v(:,j) = vprv(:,node)
	  end do
        end do
        do j=1,nnodes
          node = nodes(j)
	  il(j) = ilhkv(node)
	  hd(j) = hkv_max(node)
	end do

        nvers = 0
	iformat = 1
	ntype = 1
        np = nnodes
        lmax = nlv
	nvar = 3
        call fem_file_write_header(iformat,iunit,dtime
     +                          ,nvers,np,lmax
     +                          ,nvar,ntype
     +                          ,nlvdi,hlv,datetime_elab,regpar)

	ivar = 1
	lmax = 1
	call ivar2string(ivar,string,isub)

        call fem_file_write_data(iformat,iunit
     +                          ,nvers,np,lmax
     +                          ,string
     +                          ,il,hd
     +                          ,lmax,znv)

	ivar = 2
	lmax = nlv
	call ivar2string(ivar,string,isub)
	stringx = trim(string) // ' x'
	stringy = trim(string) // ' y'

        call fem_file_write_data(iformat,iunit
     +                          ,nvers,np,lmax
     +                          ,stringx
     +                          ,il,hd
     +                          ,nlvdi,uprv)

        call fem_file_write_data(iformat,iunit
     +                          ,nvers,np,lmax
     +                          ,stringy
     +                          ,il,hd
     +                          ,nlvdi,vprv)

	return
   74	continue
        write(6,*) 'error opening file ',trim(file)
        stop 'error stop write_nodes_scal: opening file'
	end


c***************************************************************
c***************************************************************
c***************************************************************

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

        subroutine write_2d_all_nodes(nnodes,nodes,cv2,it,ivar)

! for nos files - obsolecent

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
c***************************************************************
c***************************************************************

        subroutine write_profile_c(iu,it,j,ki,ke,lmax,ivar,h,z,c,hlv)

	use shyfem_strings

        implicit none

	integer iu
        integer it,j,ki,ke
        integer lmax
        integer ivar
        real z,h
        real c(lmax)
        real hlv(lmax)

        logical bcenter
        integer l
        integer nlvaux,nsigma
        real hsigma
        real uv
        real hd(lmax)
        real hl(lmax)

        bcenter = .true.        !depth at center of layer
        call get_sigma_info(nlvaux,nsigma,hsigma)
        call get_layer_thickness(lmax,nsigma,hsigma,z,h,hlv,hd)
        call get_depth_of_layer(bcenter,lmax,z,hd,hl)

        write(iu,*) it,j,ke,ki,lmax,ivar
        do l=1,lmax
          write(iu,*) hl(l),c(l)
        end do

        end

c***************************************************************

        subroutine write_profile_uv(iu,it,j,ki,ke,lmax,h,z,u,v,hlv)

        implicit none

	integer iu
        integer it,j,ki,ke
        integer lmax
        real z,h
        real u(lmax)
        real v(lmax)
        real hlv(lmax)

        logical bcenter
        integer l
        integer nlvaux,nsigma
        real hsigma
        real uv
        real hd(lmax)
        real hl(lmax)

        bcenter = .true.        !depth at center of layer
        call get_sigma_info(nlvaux,nsigma,hsigma)
        call get_layer_thickness(lmax,nsigma,hsigma,z,h,hlv,hd)
        call get_depth_of_layer(bcenter,lmax,z,hd,hl)

        write(iu,*) it,j,ke,ki,lmax,z
        do l=1,lmax
          uv = sqrt( u(l)**2 + v(l)**2 )
          write(iu,*) hl(l),u(l),v(l),uv
        end do

        end

c***************************************************************

	subroutine make_iunit_name(short,modi,dim,j,iu)

	implicit none

	character*(*) short,modi,dim
	integer j
	integer iu

	logical bopen
	character*5 numb
	character*80 dimen
	character*80 name

        write(numb,'(i5)') j
        numb = adjustl(numb)

	dimen = '.' // trim(dim) // '.'
	name = trim(short)//trim(modi)//trim(dimen)//numb

	inquire(file=name,opened=bopen)
	if( bopen ) goto 99
	iu = iu + 1
	inquire(unit=iu,opened=bopen)
	if( bopen ) goto 98

        !write(6,*) 'opening file : ',iu,trim(name)
        open(iu,file=name,form='formatted',status='unknown')

	return
   99	continue
	write(6,*) 'file already open: ',trim(name)
	stop 'error stop make_iunit_name: internal error (1)'
   98	continue
	write(6,*) 'unit already open: ',iu
	stop 'error stop make_iunit_name: internal error (2)'
	end

c***************************************************************

	subroutine adjust_with_isub(short,isub)

	implicit none

	character*(*) short
	integer isub

	integer i
	character*4 sub

	write(sub,'(i4)') isub

	do i=1,4
	  if( sub(i:i) == ' ' ) sub(i:i) = '0'
	end do
	sub(1:1) = '_'

	short = trim(short) // sub

	end

c***************************************************************

