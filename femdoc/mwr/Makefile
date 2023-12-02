
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------


DIR = mwr

default:

discret:
	latex $@.tex
	dvips $@.dvi
	ps2pdf $@.ps
	cp $@.pdf ~/Downloads
	gv $@.ps

2d: discret2d
discret2d:
	latex $@.tex
	dvips $@.dvi
	ps2pdf $@.ps
	cp $@.pdf ~/Downloads
	gv $@.ps

3d: discret3d
discret3d:
	latex $@.tex
	dvips $@.dvi
	ps2pdf $@.ps
	cp $@.pdf ~/Downloads
	gv $@.ps

clean:
	-rm -f *.dvi *.ps *.aux *.log *.pdf

zip:
	zip $(DIR) *
