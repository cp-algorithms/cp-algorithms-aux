#!/bin/bash
pdflatex LCA_Euler.tex
convert -density 300 LCA_Euler.pdf -quality 100 -resize 350 -depth 4 LCA_Euler.png
