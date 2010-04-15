
default:

all: anim
anim:
	./makeanim.sh

clean:
	-rm -f plot.*.ps
	-rm -f plot.*.eps
	-rm -f plot.*.eps.bak
	-rm -f tmp.gif

cleanall: clean
	-rm -f plot.ps
	-rm -f plot.*.gif

cleantotal: cleanall
	-rm -f anim.gif

save: zip
zip: clean
	-zip anim Makefile makeanim.sh rename.pl

