
!--------------------------------------------------------------------------
!
!    Copyright (C) 2016-2020  Georg Umgiesser
!    Copyright (C) 2020  Leslie Aveytua
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

! shyelab_nodes.f: utility for extracting nodes
!
! revision log :
!
! 27.06.2016	ggu	changed VERS_7_5_16
! 09.09.2016	ggu	changed VERS_7_5_17
! 09.05.2017	ggu	changed VERS_7_5_26
! 07.10.2017	ggu	restructured
! 04.11.2017	ggu	changed VERS_7_5_34
! 22.02.2018	ggu	changed VERS_7_5_42
! 31.08.2018	ggu	changed VERS_7_5_49
! 15.10.2018	ggu	added option -coord (bcoord,scoord)
! 16.02.2019	ggu	changed VERS_7_5_60
! 18.03.2020	ggu&laa	better writing of scalar values
! 13.06.2020	ggu	bug fix in write_nodes_hydro_fem()
! 25.06.2021	ggu	better legend for node profiles
!
!***************************************************************

	subroutine initialize_nodes

	use elabutil
	use basin
	use mod_depth

	implicit none

	integer j,ki,ke
	integer nodes_dummy(1)

          if( bnodes ) then
            nnodes = 0
            call get_node_file(nodefile,nnodes,nodes_dummy)
            allocate(nodesi(nnodes))
            allocate(nodese(nnodes))
            call get_node_file(nodefile,nnodes,nodese)
          else if( bnode ) then
            nnodes = 0
            call get_node_list(nodelist,nnodes,nodes_dummy)
            allocate(nodesi(nnodes))
            allocate(nodese(nnodes))
            call get_node_list(nodelist,nnodes,nodese)
          else if( bcoord ) then
            nnodes = 1
            allocate(nodesi(nnodes))
            allocate(nodese(nnodes))
	    call get_closest_node(scoord,nodesi(1))
	    nodese = nodesi
            call convert_to_external_nodes(nnodes,nodese)
          end if

	  if( nnodes <= 0 ) return

          if( bnode ) bnodes = .true.
          if( bcoord ) bnodes = .true.
 
	  nodesi = nodese
          call convert_to_internal_nodes(nnodes,nodesi)

	  if( .not. bquiet ) then
	    write(6,*) 'Nodes to be extracted: ',nnodes
            !write(6,*) nodese
            write(6,*) '          j        node'
     +			//'            x                y'
     +			//'            depth'
	    do j=1,nnodes
	      ki = nodesi(j)
	      ke = nodese(j)
	      write(6,*) j,ke,xgv(ki),ygv(ki),hkv(ki)
	    end do
	  end if

	end subroutine initialize_nodes

!***************************************************************
!***************************************************************
!***************************************************************

	subroutine write_nodes(atime,ftype,nndim,nvar,ivars,cv3all)

! manages writing of nodes

	use basin
	use levels

	implicit none

	double precision atime
	integer ftype
	integer nndim
	integer nvar
	integer ivars(nvar)
	real cv3all(nlvdi,nndim,0:nvar)

	logical bhydro,bscalar
	real znv(nkn)
	real uprv(nlvdi,nkn)
	real vprv(nlvdi,nkn)

	bhydro = ( ftype == 1 )
	bscalar = ( ftype == 2 )

        if( bhydro ) then    !hydro output
          call prepare_hydro(.true.,nndim,cv3all,znv,uprv,vprv)
          call write_nodes_hydro_ts(atime,znv,uprv,vprv)
          call write_nodes_hydro_fem(atime,znv,uprv,vprv)
        else if( bscalar ) then
          call write_nodes_scalar_ts(atime,nndim,nvar,ivars,cv3all)
          call write_nodes_scalar_fem(atime,nndim,nvar,ivars,cv3all)
        end if

	end subroutine write_nodes

!***************************************************************

	subroutine write_nodes_final(ftype,nvar,ivars)

! final write to show what files have been written

	use elabutil
	use shyfem_strings

	implicit none

	integer ftype
	integer nvar
	integer ivars(nvar)

	logical bhydro,bscalar
	integer i,iv,ivar
	character*40 format,range
	character*20 filename
	character*15 short
	character*40 full

	integer, parameter :: niu = 6
        character*5, save :: what(niu) = (/'velx ','vely ','zeta '
     +                          ,'speed','dir  ','all  '/)
        character*26, save :: descrp(niu) = (/
     +           'velocity in x-direction   '
     +          ,'velocity in y-direction   '
     +          ,'water level               '
     +          ,'current speed             '
     +          ,'current direction         '
     +          ,'all hydrodynamic variables'
     +					/)

	bhydro = ( ftype == 1 )
	bscalar = ( ftype == 2 )

        write(6,*) 'output has been written to the following files: '
	filename = 'out.fem'
        write(6,*) '  ',filename,'all variables in fem format'

	if( bhydro ) then
          write(6,*) '  what.dim.point'
          write(6,*) 'what is one of the following:'
	  call write_special_vars(niu,what,descrp)      !write hydro variables
	  call write_special_vars(1,'vel_p'
     +				,'profile for velocities')
          write(6,*) 'dim is 2d or 3d'
          write(6,*) '  2d for depth averaged variables'
          write(6,*) '  3d for output at each layer'
	  call compute_range(nnodes,range)
	  write(6,'(a,a)') ' point is consecutive numbering: '
     +				,trim(range)
	  write(6,*) 'the 5 columns in vel_p.3d.* are:'
	  write(6,*) '  depth,velx,vely,speed,dir'
	end if

	if( bscalar ) then
	  write(6,*) '  what.dim.txt        surface values of all nodes'
          write(6,*) '  what.dim.point      values for single nodes'
          write(6,*) 'what is one of the following:'
	  call write_vars(nvar,ivars)
	  call write_extra_vars(nvar,ivars,'_p',' (profile)')
	  call write_extra_vars(nvar,ivars,'_s',' (surface)')
          write(6,*) 'dim is 2d or 3d'
          write(6,*) '  2d for depth averaged variables'
          write(6,*) '  3d for output at each layer'
	  call compute_range(nnodes,range)
	  write(6,'(a,a)') ' point is consecutive numbering: '
     +				,trim(range)
	  write(6,*) 'the 2 columns in *_p.3d.* are:'
	  write(6,*) '  depth,value'
	end if

	end subroutine write_nodes_final

!***************************************************************
!***************************************************************
!***************************************************************

        subroutine convert_to_external_nodes(n,nodes)

        use basin

        implicit none

        integer n
        integer nodes(n)

        integer ne,ni,i
        integer ipext

        if( n <= 0 ) return

	do i=1,n
	  ni = nodes(i)
          if( ni <= 0 ) goto 99
          ne = ipext(ni)
          if( ne <= 0 ) goto 98
	  nodes(i) = ne
        end do

	return
   98	continue
        write(6,*) 'cannot find internal node: ',ni
        stop 'error stop convert_to_external_nodes: no such node'
   99	continue
        write(6,*) 'cannot convert internal node: ',ni
        stop 'error stop convert_to_external_nodes: no such node'
        end

!***************************************************************

        subroutine convert_to_internal_nodes(n,nodes)

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
        stop 'error stop convert_to_internal_nodes: no such node'
   99	continue
        write(6,*) 'cannot convert node: ',ne
        stop 'error stop convert_to_internal_nodes: no such node'
        end

!***************************************************************

	subroutine get_node_list(list,n,nodes)

! for n == 0 only checks how many nodes to read
! for n > 0 reads nodes into nodes() (error if n is too small)

	implicit none

	character*(*) list
	integer n
	integer nodes(n)

	integer nscan,ndim
	double precision d(n)
	integer iscand

	ndim = n
	nscan = iscand(list,d,0)

	if( nscan < 0 ) then
	  goto 98
	else if( ndim <= 0 ) then		!only count
	  n = nscan
	else if( nscan > ndim ) then
	  goto 96
	else if( nscan == 0 ) then
	  goto 97
	else
	  nscan = iscand(list,d,ndim)
	  n = nscan
	  nodes(1:n) = nint(d(1:n))
	end if
	  
	return
   96	continue
	write(6,*) 'nscan,ndim :',nscan,ndim
	write(6,*) 'list: ',trim(list)
	stop 'error stop get_node_list: dimension error'
   97	continue
	write(6,*) 'no data in list: ',trim(list)
	stop 'error stop get_node_list: no data'
   98	continue
	write(6,*) 'nscan: ',nscan
	write(6,*) 'error in read of list: ',trim(list)
	stop 'error stop get_node_list: read error'
	end

!***************************************************************

	subroutine get_node_file(file,n,nodes)

! for n == 0 only checks how many nodes to read
! for n > 0 reads nodes into nodes() (error if n is too small)

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
	stop 'error stop get_node_file: dimension error'
   97	continue
	write(6,*) 'no data in file ',trim(file)
	stop 'error stop get_node_file: no data'
   98	continue
	write(6,*) 'read error in record ',n
	write(6,*) 'in file ',trim(file)
	stop 'error stop get_node_file: read error'
   99	continue
	write(6,*) 'file: ',trim(file)
	stop 'error stop get_node_file: cannot open file'
	end

!***************************************************************

	subroutine get_closest_node(scoord,node)

        use basin

        implicit none

	character*(*) scoord
	integer node

	integer icoord
        integer k,kclose
	real x,y,xp,yp,xydist,dist

        integer iscanf
	real f(3)

        icoord = iscanf(scoord,f,3)

        if( icoord > 0 .and. icoord /= 2 ) then
          write(6,*) 'parse error for -coord: need two numbers'
          write(6,*) '-coord:  ',trim(scoord)
          stop 'error stop get_closest_node: need 2 values'
        end if

        xp = f(1)
        yp = f(2)
	kclose = 0
        xydist = 1.e+30

	do k=1,nkn
	  x = xgv(k)
	  y = ygv(k)
          dist = (x-xp)**2 + (y-yp)**2
          if( dist < xydist ) then
	    kclose = k
            xydist = dist
          end if
        end do

	if( kclose == 0 ) then
	  write(6,*) 'kclose = ',kclose
	  stop 'error stop get_closest_node: internal error (1)'
	end if

	x = xgv(kclose)
	y = ygv(kclose)
        write(6,*) 'requested coords:      ',xp,yp
        write(6,*) 'extracted coords:      ',x,y

	node = kclose

	end

!***************************************************************
!***************************************************************
!***************************************************************

	subroutine write_nodes_scalar_ts(atime,nndim,nvar,ivars,cv3all)

! writes scalar values for single nodes

	use shyfem_strings

	use levels
	use elabutil
	use mod_depth

	implicit none

	double precision atime
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
	real scals(nnodes)
	integer isubs(nvar)
	character*80 format,name,line
	character*10 numb,short
	character*20 fname,fulln
	character*20 dline
	character*20 date
	character*20 filename(nvar),fullname(nvar)
	integer, save :: icall = 0
	integer, save :: iuall2d = 0
	integer, save :: iuall3d = 0
	integer, save, allocatable :: ius(:,:,:)

	real cv3(nlvdi,nndim)	!to be deleted

	if( nnodes <= 0 ) return
	if( .not. bcompat ) return

	iu = 0
	b3d = nlv > 1
	date = '#               date'

!-----------------------------------------------------------------
! open files
!-----------------------------------------------------------------

	if( icall == 0 ) then
	  
	  allocate(ius(nvar,nnodes,5))
	  ius = 0

	  do iv=1,nvar
	    call ivar2filename(ivars(iv),filename(iv))
	    call strings_get_full_name(ivars(iv),fullname(iv))
	  end do

	  do j=1,nnodes
            write(numb,'(i5)') j
            numb = adjustl(numb)
	    do iv=1,nvar
	      fulln = fullname(iv)
	      fname = filename(iv)
	      line = date//'  '//trim(fulln)
	      call make_iunit_name(fname,'','2d',j,iu)		!verically aver
	      write(iu,'(a)') trim(line)//' (vertically averaged)'
	      ius(iv,j,2) = iu
	      if( .not. b3d ) cycle
	      call make_iunit_name(fname,'','3d',j,iu)		!3d
	      write(iu,'(a)') trim(line)//' (all layers)'
	      ius(iv,j,3) = iu
	      call make_iunit_name(fname,'_p','3d',j,iu)	!profile
	      ius(iv,j,1) = iu
	      call make_iunit_name(fname,'_s','2d',j,iu)	!surface
	      write(iu,'(a)') trim(line)//' (surface layer)'
	      ius(iv,j,4) = iu
	    end do
	  end do

	  name = 'all_scal_nodes.2d.txt'
	  call get_new_unit(iuall2d)
          open(iuall2d,file=name,form='formatted',status='unknown')

	  if( b3d ) then
	    name = 'all_scal_nodes.3d.txt'
	    call get_new_unit(iuall3d)
            open(iuall3d,file=name,form='formatted',status='unknown')
	  end if

	  do iv=1,nvar			!all nodes in one file
	    fname = filename(iv)
	    call make_iunit_name(fname,'_s','2d',0,iu)
	    write(iu,'(a)') trim(line)//' (surface layer - all nodes)'
	    ius(iv,1,5) = iu
	  end do

	end if

	icall = icall + 1

!-----------------------------------------------------------------
! write files
!-----------------------------------------------------------------

        call dts_format_abs_time(atime,dline)

        do j=1,nnodes
          ki = nodesi(j)
          ke = nodese(j)
	  lmax = ilhkv(ki)
          h = hkv(ki)
          z = cv3all(1,ki,0)
	  format = ' '
	  write(format,'(a,i3,a)') '(a20,',lmax,'f8.3)'

	  do iv=1,nvar
	    ivar = ivars(iv)
	    scal = cv3all(:,ki,iv)
	    call average_vertical_node(lmax,hlv,z,h,scal,s0)
	    iu = ius(iv,j,2)
	    write(iu,format) dline,s0
            write(iuall2d,'(a20,5i10)') dline,j,ke,ki,lmax,ivar
            write(iuall2d,*) s0
	    if( .not. b3d ) cycle
	    iu = ius(iv,j,3)
	    write(iu,format) dline,scal(1:lmax)
	    iu = ius(iv,j,1)
	    call write_profile_c(iu,dline,j,ki,ke,lmax,ivar,h,z,scal,hlv)
	    iu = ius(iv,j,4)
	    write(iu,format) dline,scal(1)
            write(iuall3d,'(a20,5i10)') dline,j,ke,ki,lmax,ivar
            write(iuall3d,*) scal(1:lmax)
	  end do
        end do

	write(format,'(a,i4,a)') '(a20,',nnodes,'f8.3)'
	do iv=1,nvar
	  ivar = ivars(iv)
	  scals = cv3all(1,nodesi(1:nnodes),iv)
	  iu = ius(iv,1,5)
	  write(iu,format) dline,scals
	end do

	end subroutine write_nodes_scalar_ts

!***************************************************************

	subroutine write_nodes_hydro_ts(atime,znv,uprv,vprv)

! writes hydro values for single nodes

	use levels
	use mod_depth
	use elabutil

	implicit none

	double precision atime
	real znv(*)
	real uprv(nlvdi,*)
	real vprv(nlvdi,*)

	logical b3d,debug
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
	character*20 dline
	character*80 name
	character*80 format

	if( nnodes <= 0 ) return
	if( .not. bcompat ) return

	debug = .false.
	b3d = nlv > 1

!-----------------------------------------------------------------
! open files
!-----------------------------------------------------------------

	if( icall == 0 ) then
	  allocate(ius(niu,nnodes,3))
	  ius = 0
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
	      call make_iunit_name('vel','_p','3d',j,iu)
	      ius(1,j,1) = iu
	    end if
	  end do
	  if( debug ) then
	    name = 'all_nodes.3d.txt'
	    call get_new_unit(iuall)
            open(iuall,file=name,form='formatted',status='unknown')
	  end if
	end if
	icall = icall + 1

!-----------------------------------------------------------------
! write files
!-----------------------------------------------------------------

        call dts_format_abs_time(atime,dline)

        do j=1,nnodes
          ki = nodesi(j)
          ke = nodese(j)
	  lmax = ilhkv(ki)

	  if( iuall > 0 ) then
            write(iuall,'(a20,5i10)') dline,j,ke,ki,lmax
            write(iuall,*) znv(ki)
            write(iuall,*) (uprv(l,ki),l=1,lmax)
            write(iuall,*) (vprv(l,ki),l=1,lmax)
	  end if

          z = znv(ki)
          h = hkv(ki)
	  u = uprv(:,ki)
	  v = vprv(:,ki)

	  call average_vertical_node(lmax,hlv,z,h,u,u0)
	  call average_vertical_node(lmax,hlv,z,h,v,v0)
	  call c2p_ocean(u0,v0,s0,d0)

	  iu = ius(1,j,2)
	  write(iu,*) dline,u0
	  iu = ius(2,j,2)
	  write(iu,*) dline,v0
	  iu = ius(3,j,2)
	  write(iu,*) dline,z
	  iu = ius(4,j,2)
	  write(iu,*) dline,s0
	  iu = ius(5,j,2)
	  write(iu,*) dline,d0
          iu = ius(6,j,2)
          write(iu,'(a20,5f12.4)') dline,z,u0,v0,s0,d0

	  if( .not. b3d ) cycle

	  write(format,'(a,i3,a)') '(a20,',lmax,'f8.3)'
	  do l=1,lmax
	    call c2p_ocean(u(l),v(l),s(l),d(l))
	  end do

	  iu = ius(1,j,3)
	  write(iu,format) dline,u(1:lmax)
	  iu = ius(2,j,3)
	  write(iu,format) dline,v(1:lmax)
	  iu = ius(3,j,3)
	  write(iu,format) dline,z
	  iu = ius(4,j,3)
	  write(iu,format) dline,s(1:lmax)
	  iu = ius(5,j,3)
	  write(iu,format) dline,d(1:lmax)

	  iu = ius(1,j,1)
          call write_profile_uv(iu,dline,j,ki,ke,lmax,h,z,u,v,hlv)
        end do

!-----------------------------------------------------------------
! end of routine
!-----------------------------------------------------------------

	end subroutine write_nodes_hydro_ts

!***************************************************************
!***************************************************************
!***************************************************************

	subroutine write_nodes_scalar_fem(atime,nndim,nvar,ivars,cv3all)

! writes FEM file out.fem - version for scalars

	use levels
	use mod_depth
	use elabutil
	use elabtime
	use shyfem_strings

	implicit none

	double precision atime
	integer nvar,nndim
	integer ivars(nvar)
	real cv3all(nlvdi,nndim,0:nvar)

	integer j,iv,node
	integer iformat,lmax,np,nvers,ntype
	integer date,time,datetime(2)
	double precision dtime
	real regpar(7)
	real cv3(nlvdi,nnodes,nvar)
	integer il(nnodes)
	real hd(nnodes)
	character*80 file,string
	character*80 strings(nvar)
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
            node = nodesi(j)
	    cv3(:,j,iv) = cv3all(:,node,iv)
	  end do
	  call ivar2femstring(ivars(iv),strings(iv))
        end do
        do j=1,nnodes
          node = nodesi(j)
	  il(j) = ilhkv(node)
	  hd(j) = hkv_max(node)
	end do

        nvers = 0
	iformat = 1
	ntype = 1
        np = nnodes
        lmax = nlv

	dtime = 0.
        call dts_from_abs_time(date,time,atime)
	datetime = (/date,time/)

        call fem_file_write_header(iformat,iunit,dtime
     +                          ,nvers,np,lmax
     +                          ,nvar,ntype
     +                          ,nlvdi,hlv,datetime,regpar)

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

!***************************************************************

	subroutine write_nodes_hydro_fem(atime,znv,uprv,vprv)

! writes FEM file out.fem - version for hydro (velocities)

	use levels
	use mod_depth
	use elabutil
	use elabtime
	!use shyelab_out

	implicit none

	double precision atime
	real znv(*)
	real uprv(nlvdi,*)
	real vprv(nlvdi,*)

	integer j,iv,node,isub
	integer iformat,lmax,np,nvers,ntype
	integer ivar,nvar
	integer date,time,datetime(2)
	double precision dtime
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

        do j=1,nnodes
          node = nodesi(j)
	  z(j) = znv(node)
	  u(:,j) = uprv(:,node)
	  v(:,j) = vprv(:,node)
	end do

        do j=1,nnodes
          node = nodesi(j)
	  il(j) = ilhkv(node)
	  hd(j) = hkv_max(node)
	end do

        nvers = 0
	iformat = 1
	ntype = 1
        np = nnodes
        lmax = nlv
	nvar = 3

	dtime = 0.
        call dts_from_abs_time(date,time,atime)
	datetime = (/date,time/)

        call fem_file_write_header(iformat,iunit,dtime
     +                          ,nvers,np,lmax
     +                          ,nvar,ntype
     +                          ,nlvdi,hlv,datetime,regpar)

	ivar = 1
	lmax = 1
	call ivar2femstring(ivar,string)

        call fem_file_write_data(iformat,iunit
     +                          ,nvers,np,lmax
     +                          ,string
     +                          ,il,hd
     +                          ,lmax,z)

	ivar = 2
	lmax = nlv
	call ivar2femstring(ivar,string)
	stringx = trim(string) // ' x'
	stringy = trim(string) // ' y'

        call fem_file_write_data(iformat,iunit
     +                          ,nvers,np,lmax
     +                          ,stringx
     +                          ,il,hd
     +                          ,nlvdi,u)

        call fem_file_write_data(iformat,iunit
     +                          ,nvers,np,lmax
     +                          ,stringy
     +                          ,il,hd
     +                          ,nlvdi,v)

	return
   74	continue
        write(6,*) 'error opening file ',trim(file)
        stop 'error stop write_nodes_scal: opening file'
	end

!***************************************************************
!***************************************************************
!***************************************************************

        subroutine write_profile_c(iu,dline,j,ki,ke,lmax,ivar,h,z,c,hlv)

	use shyfem_strings

        implicit none

	integer iu
	character*20 dline
        integer j,ki,ke
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
	character*80 header1,header2

	header1 = '#      date_and_time     point'
     +			// '  ext-node  int-node    layers'
     +			// '      ivar'

        bcenter = .true.        !depth at center of layer
        call get_sigma_info(nlvaux,nsigma,hsigma)
        call get_layer_thickness(lmax,nsigma,hsigma,z,h,hlv,hd)
        call get_depth_of_layer(bcenter,lmax,z,hd,hl)

        !write(iu,'(a20,5i10)') dline,j,ke,ki,lmax,ivar
        !write(iu,'(a)') '#               date      lmax'
        write(iu,'(a)') header1
        write(iu,'(a20,5i10)') dline,j,ke,ki,lmax,ivar
        !write(iu,'(a20,5i10)') dline,lmax
        write(iu,'(a)') '#       depth       value'
        do l=1,lmax
          write(iu,*) hl(l),c(l)
        end do

        end

!***************************************************************

        subroutine write_profile_uv(iu,dline,j,ki,ke,lmax,h,z,u,v,hlv)

        implicit none

	integer iu
	character*20 dline
        integer j,ki,ke
        integer lmax
        real z,h
        real u(lmax)
        real v(lmax)
        real hlv(lmax)

        logical bcenter
        integer l
        integer nlvaux,nsigma
        real hsigma
        real uv,dir
        real hd(lmax)
        real hl(lmax)
	character*80 header1,header2

	header1 = '#      date_and_time     point'
     +			// '  ext-node  int-node    layers'
     +			// '        zeta'
	header2 = '#        depth         x-vel         y-vel'
     +			// '         speed     direction'

        bcenter = .true.        !depth at center of layer
        call get_sigma_info(nlvaux,nsigma,hsigma)
        call get_layer_thickness(lmax,nsigma,hsigma,z,h,hlv,hd)
        call get_depth_of_layer(bcenter,lmax,z,hd,hl)

        write(iu,'(a)') header1
        write(iu,'(a20,4i10,f12.3)') dline,j,ke,ki,lmax,z
        write(iu,'(a)') header2
        do l=1,lmax
	  call c2p_ocean(u(l),v(l),uv,dir)
          write(iu,'(5f14.3)') hl(l),u(l),v(l),uv,dir
        end do

        end

!***************************************************************

