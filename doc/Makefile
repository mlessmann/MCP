IMAGES := $(wildcard *.jpg)
IMAGES += $(wildcard *.eps)

.PHONY: all
all: Protokoll.pdf Praesentation.pdf

%.pdf: %.tex %.aux $(IMAGES)
	pdflatex $<

%.aux: %.tex
	pdflatex $<

.PHONY: clean
clean:
	rm -f *.aux
	rm -f *.log
	rm -f *-eps-converted-to.pdf
	rm -f *.nav *.out *.snm *.toc
