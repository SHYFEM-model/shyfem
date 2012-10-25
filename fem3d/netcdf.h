
c has to be still commented ... FIXME

        integer rec_varid
        integer dimids_2d(5)
        integer dimids_3d(5)
        integer coord_varid(9)

        common /nc_common/ rec_varid,dimids_2d,dimids_3d,coord_varid

        save /nc_common/

