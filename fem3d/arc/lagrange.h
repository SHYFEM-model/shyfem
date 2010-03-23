	
	

        integer nbdy         !numero totale di elementi
        common /nbdy/nbdy

        integer est(nbdydim) ! numero elemento iniziale contenente particella iesima
        common /est/est

        real*8 xst(nbdydim) ! coordinata iniziale x particella iesima
        common /xst/xst

        real*8 yst(nbdydim) ! coordinata iniziale y particella iesima
        common /yst/yst

        integer ie_body(nbdydim) ! numero elemento contenente particella iesima
        common /ie_body/ie_body

        real*8 x_body(nbdydim) ! coordinata x particella iesima
        common /x_body/x_body

        real*8 y_body(nbdydim) ! coordinata y particella iesima
        common /y_body/y_body

        real*8 fx(3,neldim)
        common /fx/fx

        real*8 vel_ie(3,neldim) ! velocita lati di ie, la
        common /vel_ie/vel_ie
                              ! la numerazione indica il lato
                              ! opposto al vertice

        integer lt_body(nbdydim) ! e il lato interno 1,2 o 3 da cui parte il body
        common /lt_body/lt_body
                                 ! e a cui corrsponde la vel_ie(1-2-3,ie)
                                 ! variabile in output

        real*8 dvert(3,neldim)
        common /dvert/dvert

