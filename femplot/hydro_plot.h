
        real uv(nkndim)
        common /uv/ uv
        real vv(nkndim)
        common /vv/ vv

        real uvnv(neldim)
        common /uvnv/ uvnv
        real vvnv(neldim)
        common /vvnv/ vvnv

        real usnv(neldim)
        common /usnv/ usnv
        real vsnv(neldim)
        common /vsnv/ vsnv
        real wsnv(nkndim)
        common /wsnv/ wsnv

	save /uv/,/vv/,/uvnv/,/vvnv/,/usnv/,/vsnv/,/wsnv/

