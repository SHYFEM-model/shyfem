
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------


DIR = femcheck

default: nothing

fem: nothing

nothing:
	@true

check:
	./check_software.sh

setup:
	./setup.sh

install:
	@echo "do not use this command..."
	#./install.sh

links:
	./link_compiler.sh

list:
	pwd; ls

depend: nothing

save:
	zip $(DIR) *

clean:
	-rm -f *.o *.out *.ext
	-rm -f ggg hhh
	-rm -f test.*
	-rm -f *.bak *.tmp
	-rm -f *~
	cd parallel; make clean
	cd copyright; make clean
	cd netcdf; make clean

cleanall: clean

