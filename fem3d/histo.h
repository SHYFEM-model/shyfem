
        integer ndim_histo,nid_histo
        parameter(ndim_histo=100,nid_histo=14758327)

        integer ncbin,ncid
        integer icount(ndim_histo+1)
        real abin(ndim_histo)

        common /ihisto/ncid,ncbin,icount
        common /rhisto/abin

        save /ihisto/,/rhisto/

