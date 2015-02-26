
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!             Main parameters of Water quality model          
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	integer nstate              !Chek biotser_write formats variables to be printed in one row               
	parameter(nstate =  24)     !Number of Water column state variables
	integer noutput             ! Number of output of WC variables(includes state and some auxilaries)
	parameter (noutput = 28)
	integer nconst
	parameter (nconst = 229)    ! Number of WC constants
	integer nsstate             !Chek biotser_write formats variables to be printed in one row 
	parameter(nsstate = 15)     !Number of Bottom Sediments state variables
	integer nsoutput            !Number of BS variables for output (includes state and some auxilaries). Just for ASCII in stations 
	parameter (nsoutput = 17) 
!	parameter (nsoutput = 12) 
	integer nsconst
	parameter (nsconst = 42)
	integer NOSLAY            
	parameter (NOSLAY =  6)     !Number of layers in bottom sediments
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    
