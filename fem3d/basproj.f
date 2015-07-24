c
c $Id: nos_elab.f,v 1.2 2009-04-07 10:43:57 georg Exp $
c
c revision log :
c
c 24.01.2011    ggu     written from scratch
c 18.11.2011    ggu     adapted to new function call
c 16.03.2012    ggu     writes also transformed lines
c
c****************************************************************

	program basproj

c projection of basin

	use basin

	implicit none

	include 'param.h'

	character*50 gfile,nfile
	character*80 title

	real hkv(nkndim)
	real hev(neldim)

	logical berror
	integer nk,ne,nl,nne,nnl
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

        call grd_read(gfile)

        call grd_get_params(nk,ne,nl,nne,nnl)
        write(6,*) 'grid info: ',nk,ne,nl

        if( nk == 0 .or. ne == 0 ) then
          write(6,*) 'nk,ne: ',nk,ne
          stop 'error stop vp: no nodes or elements in basin'
        end if

        call grd_to_basin

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

c	1	Gauss-Boaga
c	2	UTM
c	3	equidistant cylindrical
c	4	UTM non standard

c---------------------------------------------------------------

c Mediterranean

c	iproj = 3	     	     !equidistant cylindrical
c        c_param(1) = 38.             !central latitude (phi)
c        c_param(2) = 15.             !longitude of origin (lon0)
c        c_param(3) = 38.             !latitude of origin (lat0)

c Nador

	iproj = 3	     	     !equidistant cylindrical
        c_param(1) = 35.             !central latitude (phi)
        c_param(2) = -3.             !longitude of origin (lon0)
        c_param(3) = 35.             !latitude of origin (lat0)

c Black Sea

c	iproj = 3		     !equidistant cylindrical
c        c_param(1) = 43.5            !central latitude (phi)
c        c_param(2) = 34.             !longitude of origin (lon0)
c        c_param(3) = 43.5            !latitude of origin (lat0)

c Klaipeda

c	iproj = 4		     !UTM Lithuania
c        c_param(1) = 24.             !longitude of origin (lon0)
c        c_param(2) = -500000.        !false easting
c        c_param(3) = 0.              !false northing
c        c_param(2) = -220000.        !false easting
c        c_param(3) = +6100000.              !false northing
c        c_param(4) = 0.9998          !scale factor

c ??

c	iproj = 2		     !UTM
c        c_param(1) = 34.             !zone
c        c_param(2) = -500000.        !false easting
c        c_param(3) = 0.              !false northing

c Laguna di Venezia

c	iproj = 1		     !Gauss-Boaga
c        c_param(1) = 2.              !fuse
c        c_param(2) = 2280000.        !shift in x
c        c_param(3) = 5000000.        !shift in y

c Laguna di Marano-Grado

c	iproj = 1		     !Gauss-Boaga
c        c_param(1) = 2.              !fuse
c        c_param(2) = 0.              !shift in x
c        c_param(3) = 0.              !shift in y

c Turkey lake for Ali

c	iproj = 2		     !UTM
c        c_param(1) = 36.             !zone
c        c_param(2) = -500000.        !false easting
c        c_param(3) = 0.              !false northing
c        c_param(2) = -0.13E+06       !false easting
c        c_param(3) = 0.418E+07       !false northing

c---------------------------------------------------------------
c do projection
c---------------------------------------------------------------

	call init_coords(iproj,c_param)
	call convert_coords(mode,nkn,xgv,ygv,xgv,ygv)	!overwrite coords

c---------------------------------------------------------------
c write new grd file
c---------------------------------------------------------------

        nfile = 'proj.grd'
        call basin_to_grd

	call grd_write(nfile)

        write(6,*) 'file has been written to ',nfile

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

	stop
	end

c***************************************************************

