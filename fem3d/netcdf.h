
        integer rec_varid
        integer dimids_2d(2)
        integer dimids_3d(3)
        integer coord_varid(6)

        common /nc_common/ rec_varid,dimids_2d,dimids_3d,coord_varid

        save /nc_common/

