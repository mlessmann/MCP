IMAGES := $(wildcard *.jpg)

Protokoll.pdf: Protokoll.tex Protokoll.aux $(IMAGES)
	pdflatex $<

Protokoll.aux: Protokoll.tex
	pdflatex $<

.PHONY: clean
clean:
	rm -f Protokoll.aux
	rm -f Protokoll.log
