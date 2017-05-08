
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

	subroutine write_nodes(dtime,ivar,cv3)

	use levels
	use elabutil
	use mod_depth

	implicit none

	double precision dtime
	integer ivar
	real cv3(nlvdi,*)

	integer j,node,it
	integer ki,ke,lmax,l
	real h,z
	real hl(nlvdi)
	character*80 format

	if( nnodes <= 0 ) return
	if( .not. bcompat ) return

        it = nint(dtime)

        do j=1,nnodes
          ki = nodes(j)
          ke = nodese(j)
	  lmax = ilhkv(ki)

	  write(1,*) it,cv3(1,ki)

          write(4,*) it,j,ke,ki,lmax,ivar
          write(4,*) (cv3(l,ki),l=1,lmax)

          write(3,*) it,j,ke,ki,lmax,ivar
          write(format,'(a,i4,a)') '(',lmax,'(f12.4))'
          write(3,format) (cv3(l,ki),l=1,lmax)

          z = 0.
          h = hkv(ki)
          call write_profile_c(it,j,ki,ke,lmax,ivar,h,z
     +				,cv3(1,ki),hlv,hl)
        end do

	end subroutine write_nodes

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

	integer j,ki,ke,lmax,it,l,k
	real z,h
	real hl(nlvdi)
	real u(nlvdi),v(nlvdi)

	if( nnodes <= 0 ) return
	if( .not. bcompat ) return

        it = nint(dtime)

        do j=1,nnodes
          ki = nodes(j)
          ke = nodese(j)
	  lmax = ilhkv(ki)
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

c***************************************************************
c***************************************************************
c***************************************************************

	subroutine write_nodes_scalar(dtime,nvar,nndim,strings,cv3all)

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

	!write(6,*) 'ggguuu: ',z,h,hl
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

