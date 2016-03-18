
# makefile for fem/bin

DIR = fembin

DEFDIR  = $(HOME)
FEMDIR  = ..
DISTDIR  = $(FEMDIR)/femdist

REMOTEDIR = /net/lagoon/usr/users/georg/fem/$(DIR)
LOCALDIR  = $(FEMDIR)/$(DIR)
UPFILES   = *

#----------------------------------------------------------

fem:

#----------------------------------------------------------

clean:
	-rm -f *.zip
	-rm -f *.o a.out
	-rm -f tmp.gnu *.bak *.ps *.tmp
	-rm -f ggg hhh
	-rm -f fort.*
	-rm -f *~
	-rm -f errout.dat

cleanall: clean

cleantotal: cleanall
	-rm -f f77 gcc

save:	zip

zip:	clean
	zip $(DIR) *

list:
	pwd; ls

depend:

femsave: save
	mv $(DIR).zip $(DISTDIR)

#----------------------------------------------------------

update:
	@cd $(REMOTEDIR); diffs -d -n -s -u $(LOCALDIR) $(UPFILES)

diffs-remote:
	@diffs -d -n -s $(REMOTEDIR) $(UPFILES)

diffs:
	@cd $(REMOTEDIR); diffs -d -n -s $(LOCALDIR) $(UPFILES)

#----------------------------------------------------------

