
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------


DIR = mwr

default:
	latex discret.tex
	dvips discret.dvi
	ps2pdf discret.ps
	cp discret.pdf ~/Downloads
	gv discret.ps

3d:
	latex discret3d.tex
	dvips discret3d.dvi
	ps2pdf discret3d.ps
	cp discret3d.pdf ~/Downloads
	gv discret3d.ps

clean:
	-rm -f *.dvi *.ps *.aux *.log *.pdf

zip:
	zip $(DIR) *
