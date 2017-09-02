
DIR = dstruc

#------------------------------------------------------

EXES = stack queue
F77 = gfortran

#------------------------------------------------------

default: all

all: $(EXES)

stack: stack.o
	$(F77) -o $@ $@.f

queue: queue.o
	$(F77) -o $@ $@.f

#------------------------------------------------------

clean:
	-rm -f *.o
	-rm -f *.mod
	-rm -f a.out
	-rm -f $(EXES)

cleanall: clean
	-rm -f $(DIR).zip

save: zip
zip: cleanall
	zip $(DIR) *

#------------------------------------------------------

.f.o:
	$(F77) -c $(FFLAGS) $<

#------------------------------------------------------

