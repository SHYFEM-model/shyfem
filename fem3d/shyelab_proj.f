!
! projection handling routines
!
! contents :
!
! subroutine shy_proj
!
! revision log :
!
! 20.07.2018    dbf     projection routines for shyelab
!
!*****************************************************************

        subroutine shy_proj

        use clo
        use basin
        use shympi
        use projection
        use coordinates

        implicit none

        character*80 sproj
        integer iscand,jj,iproj1
        integer                         :: mode,i
        double precision		:: c_param(9)
        double precision		:: d(9)
        real                            :: xgeov1(nkn),ygeov1(nkn) 

	logical is_geographical
        
        !mode: +1: cart to geo  -1: geo to cart                
        !iproj1 1=GB 2=UTM 3=CPP 4=UTM non standard
        
        call clo_get_option('proj',sproj)

        jj=iscand(sproj,d,6)		!mode,iproj,c_param

        mode = nint(d(1))
        iproj1 = nint(d(2))
        c_param(1)=d(3)
        c_param(2)=d(4)
        c_param(3)=d(5)
        if( jj>=6 ) c_param(4)=d(6)

        !write(6,*)mode,iproj1,c_param

        call init_coords(iproj1,c_param)

        if( mode == 1) then
                call convert_coords(1,nkn,xgv,ygv,xgeov1,ygeov1)
        else
		if( .not. is_geographical(nkn,xgv,ygv) ) goto 99
                call convert_coords(-1,nkn,xgv,ygv,xgeov1,ygeov1)
        endif

        xgv = xgeov1
        ygv = ygeov1

	return
   99	continue
        write(6,*)'coordinates are not lat,lon'
        write(6,*)'  xmin,xmax: ',minval(xgv),maxval(xgv)
        write(6,*)'  ymin,ymax: ',minval(ygv),maxval(ygv)
        stop' error stop shy_proj: not lat/lon'
        end subroutine shy_proj
       
!*****************************************************************

