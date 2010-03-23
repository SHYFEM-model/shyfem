        call qopen
        call qstart
        call qline(1.,1.,10.,10.)
        call qend
        call qstart
        call qtext(2.,2.,'Hello World!')
        call qend
        call qclose
        end

