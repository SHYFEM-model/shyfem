
        integer nsidim
        parameter( nsidim = 20 )

        integer nbsig
        common /nbsig/nbsig
        save /nbsig/

        real siglev(0:nsidim)
        common /siglev/siglev
        save /siglev/

        real sigdep(nkndim)
        common /sigdep/sigdep
        save /sigdep/

