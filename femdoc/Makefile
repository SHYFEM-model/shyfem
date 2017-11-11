
DOC = fem3d
DOC = shyfem

default: nothing

nothing:
	@true

.PHONY: clean ps final
.IGNORE: clean

#--------------------------------

latex:
	latex $(DOC)

bibtex:
	bibtex $(DOC)

tex:
	./femdoc.sh

check:
	tex --version

pdf:
	pdflatex $(DOC)

ps:
	dvips -o $(DOC).ps $(DOC)

final:
	cp shyfem.pdf final
	#cp shyfem.ps final

#--------------------------------

fem:    nothing
doc:	all
docs:	all
all:	dvi post final nothing

dvi:	check tex latex bibtex
	latex $(DOC)
	latex $(DOC)

post:	ps pdf

list:
	pwd; ls

depend: nothing

#--------------------------------

spell:
	ispelling4 *.tex

#--------------------------------

cleanall:	clean cleanps

clean:
	-rm -f *.aux *.dvi *.log *.toc *.bbl *.blg
	-rm -f S_*.tex P_*.tex
	-rm -f $(DOC)*.html
	-rm -f femdoc.zip
	-rm -f ggg hhh
	-rm -f *.tmp *.bak
	-rm -f *~

cleanps:
	-rm -f *.ps *.pdf

#--------------------------------

zip:	cleanall
	zip -r femdoc *

save:	zip

#--------------------------------

