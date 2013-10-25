c
c netcdf utility routines - header
c

        integer dimids_2d(5)	!dimensions for 2D case
        integer dimids_3d(5)	!dimensions for 3D case

        integer rec_varid	!id for time
        integer coord_varid(9)	!ids for coordinates

        common /nc_common/ rec_varid,dimids_2d,dimids_3d,coord_varid
        save /nc_common/

