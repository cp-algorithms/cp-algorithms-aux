#!/bin/bash
pdflatex 2SAT.tex
convert -density 300 2SAT.pdf -quality 100 -resize 732 -depth 4 2SAT.png
