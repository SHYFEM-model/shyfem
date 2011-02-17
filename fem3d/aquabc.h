
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C             Main parameters of Water quality model          
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	integer nstate              !Chek biotser_write formats variables to be printed in one row               
	parameter(nstate =  25)     !Number of Water column state variables
	integer noutput             ! Number of output of WC variables(includes state and some auxilaries)
	parameter (noutput = 25)
	integer nsstate             !Chek biotser_write formats variables to be printed in one row 
	parameter(nsstate = 12)     !Number of Bottom Sediments state variables
	integer nsoutput            !Number of BS variables for output (includes state and some auxilaries). Just for ASCII in stations 
	parameter (nsoutput = 14) 
c	parameter (nsoutput = 12) 	
	integer NOSLAY            
	parameter (NOSLAY =  3)     !Number of layers in bottom sediments
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    
