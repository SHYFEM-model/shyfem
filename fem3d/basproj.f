!
! $Id: nos_elab.f,v 1.2 2009-04-07 10:43:57 georg Exp $
!
! revision log :
!
! 24.01.2011    ggu     written from scratch
! 18.11.2011    ggu     adapted to new function call
! 16.03.2012    ggu     writes also transformed lines
!
!****************************************************************

        program basproj

! projection of basin

        use basin
        use projection

        implicit none

        character(50) :: gfile,nfile
        character(80) :: title

        logical :: berror
        integer :: nk,ne,nl,nne,nnl
        integer :: mode,iproj,isphe
        double precision, dimension(9) :: c_param

!---------------------------------------------------------------
! open grd file
!---------------------------------------------------------------

        write(6,*)
        write(6,*) 'I need the name of the grid file '
        write(6,*)
        write(6,*) 'Enter file name: '
        read(5,'(a)') gfile
        if( gfile == ' ' ) stop
        write(6,*) 'grid is read from file : ', gfile
        write(6,*)

        call grd_read(gfile)

        call grd_get_params(nk,ne,nl,nne,nnl)
        write(6,*) 'grid info: ',nk,ne,nl

        if( nk == 0 .OR. ne == 0 ) then
            write(6,*) 'nk,ne: ',nk,ne
            stop 'error stop vp: no nodes or elements in basin'
        end if

        call grd_to_basin

        call check_spheric_ev
        call get_coords_ev(isphe)

        mode = 1		!+1: cart to geo  -1: geo to cart
        if( isphe == 1 ) mode = -1
        write(6,*) 'isphe,mode: ',isphe,mode
        if( mode == 1 ) then
            write(6,*) 'converting from cartesian to geographical'
        else
            write(6,*) 'converting from geographical to cartesian'
        end if

	c_param = 0.

!---------------------------------------------------------------
! parameters for projection
!---------------------------------------------------------------

!	1	Gauss-Boaga
!	2	UTM
!	3	equidistant cylindrical
!	4	UTM non standard

!---------------------------------------------------------------

! Mediterranean

!	iproj = 3	     	     !equidistant cylindrical
!        c_param(1) = 38.             !central latitude (phi)
!        c_param(2) = 15.             !longitude of origin (lon0)
!        c_param(3) = 38.             !latitude of origin (lat0)

! Nador

!	iproj = 3	     	     !equidistant cylindrical
!        c_param(1) = 35.             !central latitude (phi)
!        c_param(2) = -3.             !longitude of origin (lon0)
!        c_param(3) = 35.             !latitude of origin (lat0)

! Black Sea

!	iproj = 3		     !equidistant cylindrical
!        c_param(1) = 43.5            !central latitude (phi)
!        c_param(2) = 34.             !longitude of origin (lon0)
!        c_param(3) = 43.5            !latitude of origin (lat0)

! Klaipeda

!	iproj = 4		     !UTM Lithuania
!        c_param(1) = 24.             !longitude of origin (lon0)
!        c_param(2) = -500000.        !false easting
!        c_param(3) = 0.              !false northing
!        c_param(2) = -220000.        !false easting
!        c_param(3) = +6100000.       !false northing
!        c_param(4) = 0.9998          !scale factor

! ??

        iproj = 2		     !UTM
        c_param(1) = 33.             !zone
        c_param(2) = -500000.        !false easting
        c_param(3) = 0.              !false northing

! Laguna di Venezia

!	iproj = 1		     !Gauss-Boaga
!        c_param(1) = 2.              !fuse
!        c_param(2) = 2280000.        !shift in x
!        c_param(3) = 5000000.        !shift in y

! Laguna di Marano-Grado

!	iproj = 1		     !Gauss-Boaga
!        c_param(1) = 2.              !fuse
!        c_param(2) = 0.              !shift in x
!        c_param(3) = 0.              !shift in y

! Turkey lake for Ali

!	iproj = 2		     !UTM
!        c_param(1) = 36.             !zone
!        c_param(2) = -500000.        !false easting
!        c_param(3) = 0.              !false northing
!        c_param(2) = -0.13E+06       !false easting
!        c_param(3) = 0.418E+07       !false northing

!---------------------------------------------------------------
! do projection
!---------------------------------------------------------------

        call init_coords(iproj,c_param)
        call convert_coords(mode,nkn,xgv,ygv,xgv,ygv)	!overwrite coords

!---------------------------------------------------------------
! write new grd file
!---------------------------------------------------------------

        nfile = 'proj.grd'
        call basin_to_grd

        call grd_write(nfile)

        write(6,*) 'file has been written to ',nfile

!---------------------------------------------------------------
! end of routine
!---------------------------------------------------------------

        stop
        end program

!***************************************************************
