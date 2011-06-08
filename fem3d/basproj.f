c
c $Id: nos_elab.f,v 1.2 2009-04-07 10:43:57 georg Exp $
c
c revision log :
c
c 24.01.2011    ggu     written from scratch
c
c****************************************************************

	program basproj

c projection of basin

	implicit none

	include 'param.h'
	include 'basin.h'

	character*50 gfile,nfile
	character*80 title

	real hkv(nkndim)
	real hev(neldim)

	logical berror
	integer ike
	integer mode,iproj,isphe
	double precision c_param(10)

c---------------------------------------------------------------
c open grd file
c---------------------------------------------------------------

        write(6,*)
        write(6,*) 'I need the name of the grid file '
        write(6,*)
        write(6,*) 'Enter file name: '
        read(5,'(a)') gfile
        if( gfile .eq. ' ' ) stop
        write(6,*) 'grid is read from file : ', gfile
        write(6,*)

	call read_grd(gfile,hkv,hev,ike)

	call check_spheric_ev
	call get_coords_ev(isphe)

	mode = 1		!+1: cart to geo  -1: geo to cart
	if( isphe .eq. 1 ) mode = -1
	write(6,*) 'isphe,mode: ',isphe,mode
	if( mode .eq. 1 ) then
	  write(6,*) 'converting from cartesian to geographical'
	else
	  write(6,*) 'converting from geographical to cartesian'
	end if

c---------------------------------------------------------------
c parameters for projection
c---------------------------------------------------------------


c mediterranean

c	iproj = 3	     	     !equidistant cylindrical
c        c_param(1) = 15.             !longitude of origin (lon0)
c        c_param(2) = 38.             !latitude of origin (lat0)
c        c_param(3) = 38.             !central latitude (phi)

c black sea

	iproj = 3		     !equidistant cylindrical
        c_param(1) = 34.             !longitude of origin (lon0)
        c_param(2) = 43.5            !latitude of origin (lat0)
        c_param(3) = 43.5            !central latitude (phi)

c Klaipeda

c	iproj = 4		     !UTM Lithuania
c        c_param(1) = 24.             !longitude of origin (lon0)
c        c_param(2) = -500000.        !false easting
c        c_param(3) = 0.              !false northing
c        c_param(4) = 0.9998          !scale factor

c---------------------------------------------------------------
c do projection
c---------------------------------------------------------------

	call init_coords(iproj,c_param)
	call convert_coords(mode,nkn,xgv,ygv,xgv,ygv)	!overwrite coords

c---------------------------------------------------------------
c write new grd file
c---------------------------------------------------------------

        nfile = 'proj.grd'
        open(1,file=nfile,status='unknown',form='formatted')
        call wrgrd(1,hkv,hev,ike)
        close(1)
        write(6,*) 'file has been written to ',nfile

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

	stop
	end

c***************************************************************


c***************************************************************

