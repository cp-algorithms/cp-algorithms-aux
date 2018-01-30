#!/bin/bash
pdflatex 2SAT_SCC.tex
convert -density 300 2SAT_SCC.pdf -quality 100 -resize 732 -depth 4 2SAT_SCC.png
