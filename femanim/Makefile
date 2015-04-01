
default:
	@echo "targets: gif avi clean cleanmore cleanall cleantotal zip"
	@./makeanim.sh -h

fem:
all:

anim: gif
gif:
	./makeanim.sh -gif

avi:
	./makeanim.sh -avi

clean:
	-rm -f plot.*.ps
	-rm -f plot.*.eps
	-rm -f plot.*.eps.bak
	-rm -f tmp.gif
	-rm -f .directory

cleanmore: clean
	-rm -f plot.*.gif
	-rm -f plot.*.jpg
	-rm -f plot.*.png

cleanall: cleanmore
	-rm -f plot.ps

cleantotal: cleanall
	-rm -f anim.gif anim.avi

save: zip
zip: clean
	-zip anim Makefile makeanim.sh rename.pl

