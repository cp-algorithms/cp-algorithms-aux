#!/bin/bash
pdflatex CartesianTreeLCA.tex
convert -density 300 CartesianTreeLCA.pdf -quality 100 -resize 550 -depth 4 CartesianTreeLCA.png
