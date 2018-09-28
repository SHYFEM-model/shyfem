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
        double precision, dimension(9)  :: c_param
        double precision d(5)
        real                            :: xgeov1(nkn),ygeov1(nkn) 
        
        !mode: +1: cart to geo  -1: geo to cart                
        !iproj1 1=GB 2=UTM 3=CPP 4=UTM non standard
        
        call clo_get_option('proj',sproj)

        jj=iscand(sproj,d,6)!mode,iproj,c_param

        mode = d(1)
        iproj1 = nint(d(2))
        c_param(1)=d(3)
        c_param(2)=d(4)
        c_param(3)=d(5)
        if( jj>=6 ) c_param(4)=d(6)

        write(6,*)mode,iproj1,c_param


        call init_coords(iproj1,c_param)
        if( mode == 1) then
                call convert_coords(mode,nkn,xgv,ygv,xgeov1,ygeov1)
        else
                if( xgv(1) >= 10000 )then
                        write(6,*)'ERROR, coordinates are not lat,lon'
                        !stop
                else
                call convert_coords(-1,nkn,xgv,ygv,xgeov1,ygeov1)
                endif
        endif
        xgv = xgeov1
        ygv = ygeov1

        end subroutine shy_proj
       