
c include file for connectivity simulations

        integer nconnect_dim
        parameter ( nconnect_dim = 50 )

        real lagr_connect_pps
        !parameter ( lagr_connect_pps = 1./200. )
        parameter ( lagr_connect_pps = 0. )

        real r_connect_radius
        parameter ( r_connect_radius = 2000. )
        !parameter ( r_connect_radius = 0. )

	integer lagr_connect_itmonth
        parameter ( lagr_connect_itmonth = 30.5*86400 )

	integer np_station
	common /int_conn/ np_station

        integer i_connect_elems(neldim)
        common /i_connect_elems/i_connect_elems

        integer i_connect_total(neldim)
        common /i_connect_total/i_connect_total

        real t_connect_total(neldim)
        common /t_connect_total/t_connect_total

        integer i_connect_released(nconnect_dim)
        common /i_connect_released/i_connect_released

        real a_connect_area(nconnect_dim)
        common /a_connect_area/a_connect_area

        integer i_connect(nconnect_dim,nconnect_dim)
        common /i_connect/i_connect

        real t_connect(nconnect_dim,nconnect_dim)
        common /t_connect/t_connect

        integer if_connect(nconnect_dim,nconnect_dim)
        common /if_connect/if_connect

        real tf_connect(nconnect_dim,nconnect_dim)
        common /tf_connect/tf_connect

	save /int_conn/
        save /i_connect_elems/

	save /i_connect/,/t_connect/
	save /if_connect/,/tf_connect/

	save /t_connect_total/,/i_connect_total/

        save /i_connect_released/,/a_connect_area/

