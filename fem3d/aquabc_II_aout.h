cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Parameters and arrays for writing state variables to ASCII files for defined nodes 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



C      Output of variables for water column

	      integer NBIOTSMX             !Maximal number of time series(nodes) for WC 
		  parameter(NBIOTSMX=100)    
	      integer NBIOTS               !Actual number of time series(nodes) for WC       
		  integer BIOTSFUN(NBIOTSMX)   !file units for WC state variables ASCII output	
		  integer BIOTSNOD(NBIOTSMX)   !nodes for WC state variables ASCII output

C      Output  variables for sediments
		   integer NBIOTSMX_sed               !Maximal number of time series(nodes) for BS
		   parameter(NBIOTSMX_sed=100)    
	       integer NBIOTS_sed                 !Actual number of time series(nodes) for BS	
		   integer BIOTSFUN_sed(NBIOTSMX_sed) !file units for BS state variables ASCII output 	
		   integer BIOTSNOD_sed(NBIOTSMX_sed) !nodes for BS state variables ASCII output



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Parameters and arrays for writing diagnostics(components of derivatives) to ASCII files for defined nodes 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C WATER COLUMN
		integer   NDIAGVAR           !Maximum number of different intermediate variables in output.  
		parameter (NDIAGVAR=26)      !Number of auxilary variables for output	
		integer   NDGTSMX            !Total number of nodes for intermediates of each state variable 
		parameter (NDGTSMX=100)
		integer NDGTS(nstate)        !Actual number of nodes for intermediates of each state variable
		integer DGTSFUN(NDGTSMX, nstate) !Output file units for EUTRO diagnostics time series(nodes) plots
		integer DGTSNOD(NDGTSMX, nstate) !Nodes for intermediate variables in output
		real    dgar(nstate,NDIAGVAR)    !dgar (Diagnostics data) for each sediments layer


C BOTTOM SEDIMENTS
		integer   NDIAGVAR_sed          !Maximum number of different intermediate variables in output. Chek biotser_write formats
		parameter (NDIAGVAR_sed=25)      !Number of auxilary variables for output	
		integer   NDGTSMX_sed           !Total number of nodes for intermediates of each state variable 
		parameter (NDGTSMX_sed=100)
		integer NDGTS_sed(nsstate)      !Actual number of nodes for intermediates of each state variable
		integer DGTSFUN_sed(NDGTSMX, nsstate) !Output file units for EUTRO diagnostics time series(nodes) plots
		integer DGTSNOD_sed(NDGTSMX, nsstate) !Nodes for intermediate variables in output
	    real    dgar_sed(NOSLAY,nsstate,NDIAGVAR_sed)    ! Diagnostics data, for temporary storage


