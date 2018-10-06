%.pdf: %.tex
	pdflatex $<

%.png: %.pdf
	convert -density 300 $< -quality 100 -resize 450x -depth 4 $@

all: point_location_goal.png point_location_events.png

clear:
	rm -f *.png *.pdf *.log *.aux *.fls *.fdb_latexmk

copy:
	cp *.png ../../../e-maxx-eng/img/
