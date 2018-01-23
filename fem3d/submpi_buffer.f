
!******************************************************************
!******************************************************************
!******************************************************************

        subroutine count_buffer(n0,nlvddi,n,nc,il,nodes,nb)

        integer n0,nlvddi,n,nc
        integer il(n)
        integer nodes(nc)
        integer nb

        integer i,k,l,lmax

        if( nlvddi == 1 ) then
          nb = nc * (2-n0)
        else
          nb = 0
          do i=1,nc
            k = nodes(i)
            lmax = il(k)
            nb = nb + lmax - n0 + 1
          end do
        end if

        end subroutine count_buffer

!******************************************************************

        subroutine to_buffer_i(n0,nlvddi,n,nc,il,nodes,val,nb,buffer)

        integer n0,nlvddi,n,nc
        integer il(n)
        integer nodes(nc)
        integer val(n0:nlvddi,n)
        integer nb
        integer buffer(nb)

        integer i,k,l,lmax

        if( nlvddi == 1 .and. n0 == 1 ) then
          do i=1,nc
            k = nodes(i)
            buffer(i) = val(1,k)
          end do
          nb = nc
        else
          nb = 0
          do i=1,nc
            k = nodes(i)
            lmax = il(k)
            do l=n0,lmax
              nb = nb + 1
              buffer(nb) = val(l,k)
            end do
          end do
        end if

        end subroutine to_buffer_i

!******************************************************************

        subroutine from_buffer_i(n0,nlvddi,n,nc,il,nodes,val,nb,buffer)

        integer n0,nlvddi,n,nc
        integer il(n)
        integer nodes(nc)
        integer val(n0:nlvddi,n)
        integer nb
        integer buffer(nb)

        integer i,k,l,lmax

        if( nlvddi == 1 .and. n0 == 1 ) then
          do i=1,nc
            k = nodes(i)
            val(1,k) = buffer(i)
          end do
          nb = nc
        else
          nb = 0
          do i=1,nc
            k = nodes(i)
            lmax = il(k)
            do l=n0,lmax
              nb = nb + 1
              val(l,k) = buffer(nb)
            end do
          end do
        end if

        end subroutine from_buffer_i

!******************************************************************

        subroutine to_buffer_r(n0,nlvddi,n,nc,il,nodes,val,nb,buffer)

        integer n0,nlvddi,n,nc
        integer il(n)
        integer nodes(nc)
        real val(n0:nlvddi,n)
        integer nb
        real buffer(nb)

        integer i,k,l,lmax

        if( nlvddi == 1 .and. n0 == 1 ) then
          do i=1,nc
            k = nodes(i)
            buffer(i) = val(1,k)
          end do
          nb = nc
        else
          nb = 0
          do i=1,nc
            k = nodes(i)
            lmax = il(k)
            do l=n0,lmax
              nb = nb + 1
              buffer(nb) = val(l,k)
            end do
          end do
        end if

        end subroutine to_buffer_r

!******************************************************************

        subroutine from_buffer_r(n0,nlvddi,n,nc,il,nodes,val,nb,buffer)

        integer n0,nlvddi,n,nc
        integer il(n)
        integer nodes(nc)
        real val(n0:nlvddi,n)
        integer nb
        real buffer(nb)

        integer i,k,l,lmax

        if( nlvddi == 1 .and. n0 == 1 ) then
          do i=1,nc
            k = nodes(i)
            val(1,k) = buffer(i)
          end do
          nb = nc
        else
          nb = 0
          do i=1,nc
            k = nodes(i)
            lmax = il(k)
            do l=n0,lmax
              nb = nb + 1
              val(l,k) = buffer(nb)
            end do
          end do
        end if

        end subroutine from_buffer_r

!******************************************************************
!******************************************************************
!******************************************************************

