%.pdf: %.tex
	pdflatex $<

%.png: %.pdf
	convert -density 300 $< -quality 100 -resize x320 -depth 4 $@

all: binary_indexed_tree.png

clear:
	rm -f *.png *.pdf *.log *.aux *.fls *.fdb_latexmk

copy:
	cp *.png ../../../e-maxx-eng/img/
