! shypart.f, v 1.10 2016/02/16 CMCC 

!----------------------------------------------------------------

	program shypart_main

        use partitioning
        use constants_par

!-----------------------------------------------------------
! read STR file
!-----------------------------------------------------------

	call cstinit

	call cstfile				!read STR and basin

        call partMetis
	
        end


