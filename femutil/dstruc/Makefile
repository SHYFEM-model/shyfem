
DIR = dstruc

#------------------------------------------------------

EXES = stack queue hash
F77 = gfortran

#------------------------------------------------------

default: all

all: $(EXES)

stack: stack.o
	$(F77) -o $@ $@.f

queue: queue.o
	$(F77) -o $@ $@.f

hash: hash.o
	$(F77) -o $@ $@.f

#------------------------------------------------------

clean:
	-rm -f *.o
	-rm -f *.mod
	-rm -f a.out
	-rm -f $(EXES)
	-rm -f ggg hhh

cleanall: clean
	-rm -f $(DIR).zip

save: zip
zip: cleanall
	zip $(DIR) *

copy:
	cp $(DIR).zip $(HOME)/C

#------------------------------------------------------

.f.o:
	$(F77) -c $(FFLAGS) $<

#------------------------------------------------------

