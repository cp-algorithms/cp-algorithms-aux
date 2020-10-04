%.pdf: %.tex
	pdflatex $<

%.png: %.pdf
	convert -density 300 $< -quality 100 -resize 350x -depth 4 $@

all: second_best_mst_1.png second_best_mst_2.png

clear:
	rm -f *.png *.pdf *.log *.aux *.fls *.fdb_latexmk

copy:
	cp *.png ../../../e-maxx-eng/img/
