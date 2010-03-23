
        integer nsidim
        parameter( nsidim = 20 )

        integer nsigma
        common /nsigma/nsigma
        save /nsigma/

        real siglev(0:nsidim)
        common /siglev/siglev
        save /siglev/

        real sigdep(nkndim)
        common /sigdep/sigdep
        save /sigdep/

