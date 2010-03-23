
        program scalar

        implicit none

        integer nlvdim,ntdim,npdim
        parameter (nlvdim=20)
        parameter (ntdim=200)
        parameter (npdim=10)   !number of points

        integer nlv,nt
        real data(nlvdim,ntdim)

        call qopen

        call read_matrix(nlvdim,ntdim,nlv,nt,data)
        call print_matrix(nlvdim,ntdim,nlv,nt,data)
        call plot_matrix(nlvdim,ntdim,nlv,nt,data)

        call qclose

        end

c******************************************************************

        subroutine plot_matrix(nlvdim,ntdim,nlv,nt,data)

        implicit none

        integer nlvdim,ntdim
        integer nlv,nt
        real data(nlvdim,ntdim)

        integer it,l
        real xmin,ymin,xmax,ymax
        real colmin,colmax,valmin,valmax
        real aux,value,col,fact
        real dx,dy,x0,y0

        colmin = 0.9
        colmax = 0.1
        valmin = 15.
        valmin = 20.
        valmax = 24.

        aux = (colmax-colmin) / (valmax-valmin)
        fact = 1.
        if( colmin .gt. colmax ) fact = -1

        call qstart

        call qgetvp(xmin,ymin,xmax,ymax)

        call qtxts(6)
        call qfont('Courier')

        xmin = xmin + 1.
        ymin = ymin + 1.
        xmax = xmax - 1.
        ymax = ymax - 15.

        call qsetvp(xmin,ymin+1.5,xmax,ymax)
        call qworld(xmin,ymin,xmax,ymax)

        call set_color_table( 0 )         !gray scale
        call qsetc(0.)          !black

        dx = (xmax-xmin) / nt
        dy = (ymax-ymin) / nlv

        call pbox(xmin,ymin,xmax,ymax)

        call set_color_table( 1 )         !hue scale

        do it=1,nt
          do l=1,nlv
            x0 = xmin + it*dx
            y0 = ymax - l*dy

            value = data(l,it)
            col = colmin + (value-valmin) * aux
            col = fact*max(fact*col,fact*colmin)
            col = fact*min(fact*col,fact*colmax)
            write(6,*) l,it,value,col

            call qsetc(col)
            call qrfill(x0-dx,y0,x0,y0+dy)
          end do
        end do

        call qend

        end

c******************************************************************

        subroutine read_matrix(nlvdim,ntdim,nlv,nt,data)

        implicit none

        integer nlvdim,ntdim
        integer nlv,nt
        real data(nlvdim,ntdim)

        integer inode,iwhat
        integer itime,i,kext,kint,lmax,ivar
        integer it,l
        real dummy

        inode = 4
        iwhat = 12

        it = 0
    1   continue
          read(5,*,end=2) itime,i,kext,kint,lmax,ivar

          if( i .eq. inode .and. ivar .eq. iwhat ) then
            it = it + 1
            if( l .gt. nlvdim ) stop 'error stop read_matrix: nlvdim'
            if( it .gt. ntdim ) stop 'error stop read_matrix: ntdim'
            read(5,*) (data(l,it),l=1,lmax)
          else
            read(5,*) (dummy,l=1,lmax)
          end if

          goto 1
    2   continue

        nlv = lmax
        nt = it

        end

c******************************************************************

        subroutine print_matrix(nlvdim,ntdim,nlv,nt,data)

        implicit none

        integer nlvdim,ntdim
        integer nlv,nt
        real data(nlvdim,ntdim)

        integer it,l

        write(6,*) nlv,nt
        do it=1,nt
          write(6,*) (data(l,it),l=1,nlv)
        end do

        end

c******************************************************************


