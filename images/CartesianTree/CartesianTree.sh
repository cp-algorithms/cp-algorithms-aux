#!/bin/bash
pdflatex CartesianTree.tex
convert -density 300 CartesianTree.pdf -quality 100 -resize 550 -depth 4 CartesianTree.png
