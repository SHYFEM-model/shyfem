
DIR = femcheck

default:

fem:

check:
	./check_setup.sh

setup:
	./setup.sh

install:
	./install.sh

links:
	./link_compiler.sh

list:
	pwd; ls

save:
	zip $(DIR) *

clean:
	-rm -f *.o *.out *.ext
	-rm -f ggg hhh
	-rm -f test.*
	-rm -f *.bak *.tmp
	-rm -f *~

cleanall: clean
	-rm -f CHECKLOG COMPILER_LOG

