
#------------------------------------------------------

DIR = unf

#------------------------------------------------------

UNFOBJS = unformatted.o

EXES = unfnos unf_diff unf_info unf_joel \
	unf_64to32 unf_32to64 unf_info32 unf_info64

#------------------------------------------------------

default: all

all: $(EXES)

unfnos: unfnos.o $(UNFOBJS)
	cc -o $@ $@.o $(UNFOBJS)

unf_diff: unf_diff.o $(UNFOBJS)
	cc -o $@ $@.o $(UNFOBJS)

unf_info: unf_info.o $(UNFOBJS)
	cc -o $@ $@.o $(UNFOBJS)

unf_info32: unf_info32.o $(UNFOBJS)
	cc -o $@ $@.o $(UNFOBJS)

unf_info64: unf_info64.o $(UNFOBJS)
	cc -o $@ $@.o $(UNFOBJS)

unf_joel: unf_joel.o $(UNFOBJS)
	cc -o $@ $@.o $(UNFOBJS)

unf_64to32: unf_64to32.o $(UNFOBJS)
	cc -o $@ $@.o $(UNFOBJS)

unf_32to64: unf_32to64.o $(UNFOBJS)
	cc -o $@ $@.o $(UNFOBJS)

fortran_test:
	gfortran test.f
	./a.out

#------------------------------------------------------

save: zip
zip: cleanall
	zip $(DIR).zip *

clean:
	-rm -f $(EXES) a.out
	-rm -f *.o
	-rm -f *.nos *.con
	-rm -f ggg hhh
	-rm -f *.dat

cleanall: clean
	-rm -f $(DIR).zip

#------------------------------------------------------

