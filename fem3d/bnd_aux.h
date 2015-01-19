
        real ruv(nkndim), rvv(nkndim)       !momentum input (2D)
        common /ruv/ruv, /rvv/rvv

        real crad(neldim)                       !$$GWI (radiation)
        common /crad/crad

	save /ruv/,/rvv/,/crad/

