	
	

        integer nbdy         !numero totale di elementi
        common /nbdy/nbdy

        integer est(nbdydim) ! numero elemento iniziale contenente particella iesima
        common /est/est

        real xst(nbdydim) ! coordinata iniziale x particella iesima
        common /xst/xst

        real yst(nbdydim) ! coordinata iniziale y particella iesima
        common /yst/yst

	real zst(nbdydim) ! coordinata iniziale z particella iesima
	common /zst/zst

		
        integer ie_body(nbdydim) ! numero elemento contenente particella iesima
        common /ie_body/ie_body

        real x_body(nbdydim) ! coordinata x particella iesima
        common /x_body/x_body

        real y_body(nbdydim) ! coordinata y particella iesima
        common /y_body/y_body

	real z_body(nbdydim) ! coordinata z particella iesima
	common /z_body/z_body

	real tin(nbdydim) ! tempo di rilascio
	common /tin/tin

	real v_lag(neldim),u_lag(neldim) ! velocità lagrangiane per il termine advectiv
	common /v_lag/v_lag,/u_lag/u_lag

	
	real tdecay      ! tempo di dimezzamento
	common /tdecay/tdecay

        real fx(3,neldim)
        common /fx/fx

        real vel_ie(3,neldim) ! velocita lati di ie, la
        common /vel_ie/vel_ie
                              ! la numerazione indica il lato
                              ! opposto al vertice

        integer lt_body(nbdydim) ! e il lato interno 1,2 o 3 da cui parte il body
        common /lt_body/lt_body
                                 ! e a cui corrsponde la vel_ie(1-2-3,ie)
                                 ! variabile in output

        real dvert(3,neldim)
        common /dvert/dvert

	real fall  		! velocitá verticale
	common /fall/fall

	real alt(neldim)	! altezza media su elemento
	common /alt/alt

