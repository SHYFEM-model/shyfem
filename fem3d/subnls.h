
        character*80 sname      !name of section
        integer snum            !number of section
        logical sread           !section has been read

        common /secsec/ snum,sread,sname
        save /secsec/

        integer unit
        common /nrdcom/ unit
	save /nrdcom/

	integer nlsdim
	parameter (nlsdim=1000)

	double precision dnlscom(nlsdim)
	common /dnlscom/dnlscom
	save /dnlscom/

