
DIR = profile

OBJSF = subfemfile.f subfil.f subdts.f

default: profile

all: profile

profile:
	gfortran -o $@ $@.f $(OBJSF)

clean:
	-rm -f a.out profile
	-rm -f *.mod
	-rm -f $(DIR).zip

cleanall: clean
	-rm -f *.fem

zip: clean
	zip $(DIR).zip *

