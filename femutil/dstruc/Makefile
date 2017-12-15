
DIR = dstruc

TOTAL = dsall.f

#------------------------------------------------------

EXES = stack queue hash tree sort spline
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

tree: tree.o
	$(F77) -o $@ $@.f

sort: sort.o
	$(F77) -o $@ $@.f

spline: spline.o
	$(F77) -o $@ $@.f

heap: heap.o
	$(F77) -o $@ $@.f

#------------------------------------------------------

test: all
	@for prog in $(EXES); do ./$$prog; done

compact:
	@rm -f $(TOTAL)
	@for prog in $(EXES); do ./compact.pl $$prog.f >> $(TOTAL) ; done

#------------------------------------------------------

clean:
	-rm -f *.o
	-rm -f *.mod
	-rm -f *.tmp
	-rm -f *.bak
	-rm -f a.out
	-rm -f $(EXES)
	-rm -f ggg hhh
	-rm -f out*.ps
	-rm -f fort.*

cleanzip: clean
	-rm -f $(DIR).zip

cleanall: clean cleanzip
	-rm -f $(TOTAL)

save: zip
zip: cleanzip
	zip $(DIR) *

copy:
	cp $(DIR).zip $(HOME)/C

#------------------------------------------------------

.f.o:
	$(F77) -c $(FFLAGS) $<

#------------------------------------------------------

