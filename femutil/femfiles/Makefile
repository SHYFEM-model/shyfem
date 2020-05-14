
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------


DIR = profile

OBJSF = subfemfile.f subfil.f subdts.f

default: profile

all: profile

profile: profile.o
	gfortran -o $@ $@.f $(OBJSF)

clean:
	-rm -f a.out profile
	-rm -f *.o *.mod
	-rm -f $(DIR).zip

cleanall: clean
	-rm -f *.fem

zip: clean
	zip $(DIR).zip *

