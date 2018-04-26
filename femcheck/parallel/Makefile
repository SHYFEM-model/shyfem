
#--------------------------------------------------

default: run

run: mpi omp mix

#--------------------------------------------------

mpi:
	./hw_mpi.sh

omp:
	./hw_omp.sh

mix:
	./hw_mix.sh

#--------------------------------------------------

clean:
	-rm -f a.out

cleanzip:
	-rm -f parallel.zip

zip: clean cleanzip
	zip parallel *

#--------------------------------------------------

