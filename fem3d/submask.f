c
c $Id: subreg.f,v 1.17 2009-09-14 08:20:58 georg Exp $
c
c routines for interpolation onto regular grid
c
c contents :
c
c subroutine setgeo(x0,y0,dx,dy,flag)
c		sets grid values for interpolation
c subroutine getgeo(x0,y0,dx,dy,flag)
c		gets grid values for interpolation
c subroutine getgeoflag(flag)
c		gets flag value for interpolation
c
c subroutine av2am(av,am,ip,jp)
c		interpolation of av onto a regular net
c subroutine av2amk(bwater,av,am,ip,jp)
c		interpolation of av onto a regular net
c function intri(x,y,xp,yp)
c		point in triangle or not
c function intrid(x,y,xp,yp)
c		point in triangle or not (double precision version)
c subroutine am2av(am,av,ip,jp)
c		interpolation of am onto finite element mesh
c function am2val(am,ip,jp,xx,yy)
c		interpolation of am onto finite element mesh
c subroutine ave2am(av,am,ip,jp)
c		interpolation of av (elementwise) onto a regular net
c
c subroutine mkmask(bwater,zv,href,hzoff)
c		makes mask in element
c
c subroutine mimareg(am,ip,jp,amin,amax)
c		computes min/max of regular matrix (without flag values)
c subroutine a2char(am,ac,ip,jp)
c		creates 1 char representation of matrix
c subroutine prchar(ac,ip,jp)
c		prints 1 char representation of matrix
c
c subroutine femintp(ie,z,xp,yp,zp)
c               interpolation in element (with ev)
c subroutine elemintp(x,y,z,xp,yp,zp)
c               interpolation in element (no ev)
c
c subroutine find_elem_from_old(ieold,xp,yp,ielem)
c		finds element for point (xp,yp) starting from ieold
c subroutine find_element(xp,yp,ielem)
c		finds element for point (xp,yp)
c function in_element(ie,xp,yp)
c		checks if point (xp,yp) is in element ie
c subroutine get_xy_elem(ie,x,y)
c		returns x,y of vertices of element ie
c
c revision log :
c
c 18.11.1998	ggu	routine commented
c 18.11.1998	ggu	routine setgeo introduced
c 19.11.1998	ggu	routines a2char, prchar added
c 19.10.1999	ggu	routine mkmask added from subutl
c 25.11.2004	ggu	new routines femintp and elemintp for interpolation
c 14.03.2005	ggu	new routines for interpolation in element
c 11.03.2009	ggu	new helper routine getgeoflag()
c 12.06.2009	ggu	passing to double precision, intrid, bug bug_f_64bit
c 26.01.2011	ggu&mb	handling extrapolation in am2av()
c 27.01.2011	ggu&ccf	bug fix in find_elem_from_old() BUG_27.01.2011
c 31.03.2011	ggu	new routine elemmask()
c 24.11.2011	ggu	new routine find_close_elem()
c 20.06.2012	ggu	new routine get_scal_elem()
c 07.10.2012	ggu	new routine av2fm()
c 10.10.2012	ggu	new routine fm2am2d() and fm2am3d()
c 26.10.2012	ggu	bug fix: do not access not existing storage
c 30.05.2014	ggu	in av2amk() do not interpolate for flag values
c 07.07.2014	ggu	new routine intp_reg()
c 25.09.2015	ggu	new routines intp_reg_nodes(), intp_reg_elems()
c 05.05.2016	ggu	file restructured (module)
c 14.05.2016	ggu	allow for extension of grid -> bregextend
c 23.06.2016	ggu	allow for eps in computing box
c 23.09.2016	ggu	allow for eps in computing box and reg intp
c 23.04.2017	ggu	new routine intp_reg_single_nodes()
c
c notes :
c
c pxareg,pyareg         coordinates of lower left point of matrix
c pxdreg,pydreg         grid size of matrix
c pzlreg                value of z for land points
c
c******************************************************
c******************************************************
c******************************************************

	subroutine set_dry_mask(bwater,zv,zev,href,hzoff)

c makes mask for dry and wet areas - zenv must be available
c
c bwater is elementwise mask:	true = water point

	use basin

	implicit none

c arguments
	logical bwater(nel)
	real zv(nkn)
	real zev(3,nel)
	real href,hzoff
c local
	integer itot,itot1
	integer ie,ii

        do ie=1,nel

          itot=0
          do ii=1,3
            if( hm3v(ii,ie)+zev(ii,ie)-href .gt. hzoff ) then
		itot=itot+1    !wet
	    end if
          end do

          itot1=0
          do ii=1,3
            if(zv(nen3v(ii,ie)).eq.zev(ii,ie)) itot1=itot1+1
          end do

          !if(itot.ne.3.or.itot1.ne.3)  bwater(ie) = .false.
          if(itot.ne.3.and.itot1.ne.3)  bwater(ie) = .false.

        end do

	end

c******************************************************

	subroutine set_level_mask(bwater,ilhv,level)

c makes mask for water points (level)
c
c bwater is elementwise mask:	true = water point

	use basin, only : nkn,nel,ngr,mbw

	implicit none

c arguments
	logical bwater(nel)
	integer ilhv(nel)
	integer level
c local
	integer ie,nedry

	nedry = 0

	!write(6,*) 'set_level_mask: ',level
	!write(6,*) (ilhv(ie),ie=1,nel,nel/10)

	do ie=1,nel
	  if( level .gt. ilhv(ie) ) then
	    bwater(ie) = .false.
	    nedry = nedry + 1
	  end if
	end do

	write(6,*) 'level exist  (dry/wet/total): ',nedry,nel-nedry,nel

	end

c******************************************************

	subroutine make_dry_node_mask(bwater,bkwater)

c makes node mask from element mask
c
c bwater is elementwise mask:	true = water point

	use basin

	implicit none

c arguments
	logical bwater(nel)
	logical bkwater(nkn)

	integer ie,ii,k
	integer nndry,nedry

	nndry = 0
	nedry = 0

	bkwater = .false.

	do ie=1,nel
	  if( bwater(ie) ) then
	    do ii=1,3
	      k = nen3v(ii,ie)
	      bkwater(k) = .true.
	    end do
	  else
	  end if
	end do

	end

c******************************************************

	subroutine make_dry_elem_mask(bwater,bkwater)

c makes elem mask from node mask
c
c bwater is elementwise mask:	true = water point

	use basin

	implicit none

c arguments
	logical bwater(nel)
	logical bkwater(nkn)
c local
	integer ie,ii,k
	integer nedry

        bwater = .true.

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    if( .not. bkwater(k) ) then
	      bwater(ie) = .false.
	    end if
	  end do
	end do

	end

c******************************************************

        subroutine info_dry_mask(bwater,bkwater)

        use basin

        implicit none

	logical bwater(nel)
	logical bkwater(nkn)
        
        integer nedry,nndry,newet,nnwet

        newet = count(bwater)
        nnwet = count(bkwater)
        nedry = nel - newet
        nndry = nkn - nnwet

	write(6,*) 'dry elements (dry/wet/total): ',nedry,newet,nel
	write(6,*) 'dry nodes    (dry/wet/total): ',nndry,nnwet,nkn

        end

c******************************************************
c******************************************************
c******************************************************

