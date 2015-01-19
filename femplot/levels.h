
        integer ilhkv(nkndim)
        common /ilhkv/ilhkv

        integer ilhv(neldim)
        common /ilhv/ilhv

        real hlv(nlvdim), hldv(nlvdim)
        common /hlv/hlv, /hldv/hldv

        integer ilmv(neldim)
        common /ilmv/ilmv

        integer ilmkv(nkndim)
        common /ilmkv/ilmkv

	save /ilhkv/,/ilhv/,/hlv/,/hldv/,/ilmv/,/ilmkv/

