
#---------------------------------------------------

default:

fem:

list:
	pwd; ls

#---------------------------------------------------

clean:
	recursivemake $@ skel

cleanall: clean cleantest cleanfem
	recursivemake $@ skel

cleantest:
	-rm -rf test_*

cleanfem:
	-rm -rf fem3d_*

#---------------------------------------------------

