c
c $Id: subproj.f,v 1.5 2010-03-11 15:36:38 georg Exp $
c
c subroutines for handling projection
c
c revision log :
c
c 17.11.2011    ggu     written from scratch
c 10.01.2012    ggu     bug fix: c_param was real
c
c****************************************************************            

	subroutine handle_projection

c handles projection - converts x/y to lat/lon

	implicit none

        include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	real xgv(nkndim)
	common /xgv/xgv
	real ygv(nkndim)
	common /ygv/ygv
	real xgeov(nkndim)
	common /xgeov/xgeov
	real ygeov(nkndim)
	common /ygeov/ygeov

	integer mode,iproj,i
	double precision c_param(5)

        real getpar

	integer icall
	save icall
        data icall / 0 /
        
        if( icall .ne. 0 ) return

	icall = 1
       
c---------------------------------------------------------------
c initialization
c---------------------------------------------------------------

	mode = 1	!from cartesian to lat/lon

        iproj = nint(getpar('iproj'))

	if( iproj .eq. 0 ) then
	  !nothing
	else if( iproj .eq. 1 ) then
	  c_param(1) = getpar('c_fuse')
	  c_param(2) = getpar('c_x0')
	  c_param(3) = getpar('c_y0')
	else if( iproj .eq. 2 ) then
	  c_param(1) = getpar('c_zone')
	  c_param(2) = getpar('c_x0')
	  c_param(3) = getpar('c_y0')
	else if( iproj .eq. 3 ) then
	  c_param(1) = getpar('c_phi')
	  c_param(2) = getpar('c_lon0')
	  c_param(3) = getpar('c_lat0')
	else if( iproj .eq. 4 ) then
	  c_param(1) = getpar('c_lamb')
	  c_param(2) = getpar('c_x0')
	  c_param(3) = getpar('c_y0')
	  c_param(4) = getpar('c_skal')
	else
	  write(6,*) 'iproj = ',iproj
	  stop 'error stop handle_projection: value for iproj not allowed'
	end if

	call init_coords(iproj,c_param)
	call convert_coords(mode,nkn,xgv,ygv,xgeov,ygeov)

	write(6,*) 'start of handle_projection'
	write(6,*) 'mode  = ',mode
	write(6,*) 'iproj = ',iproj

	write(6,*) (xgv(i),i=1,5)
	write(6,*) (ygv(i),i=1,5)
	write(6,*) (xgeov(i),i=1,5)
	write(6,*) (ygeov(i),i=1,5)

	write(6,*) 'end of handle_projection'

	end

c****************************************************************            

