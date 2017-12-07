
FEMDIR    = ..
FEMBIN    = $(FEMDIR)/fembin

FEMDIRS   = mpi

#-----------------------------------------------------------------------

default: mpi

fem: mpi

mpi:
	mpif90 mpi_test.f
	/usr/bin/mpirun -np 4 a.out

regress:

list:
	pwd; ls

depend:

clean:
	-rm -f a.out

cleanall: clean

save:

zip:

#-----------------------------------------------------------------------

