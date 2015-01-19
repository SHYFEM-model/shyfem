
        real zov(nkndim), znv(nkndim)
        common /zov/zov, /znv/znv
        real zeov(3,neldim), zenv(3,neldim)       !$$ZEONV
        common /zeov/zeov, /zenv/zenv

	save /zov/,/znv/,/zeov/,/zenv/

        real utlov(nlvdim,neldim)
        common /utlov/utlov
        real utlnv(nlvdim,neldim)
        common /utlnv/utlnv
        real vtlov(nlvdim,neldim)
        common /vtlov/vtlov
        real vtlnv(nlvdim,neldim)
        common /vtlnv/vtlnv

	save /utlov/,/utlnv/,/vtlov/,/vtlnv/

