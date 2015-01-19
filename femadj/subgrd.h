
        real xscale,yscale,zscale
        common /vscale/ xscale,yscale,zscale
        save /vscale/

        integer nin,iline,ianz
        common /grdcom_i/ nin,iline,ianz
        real f(80)
        common /grdcom_r/ f
        character*132 line
        common /grdcom_c/ line
        save /grdcom_i/, /grdcom_r/, /grdcom_c/

