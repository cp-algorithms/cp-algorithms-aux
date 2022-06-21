#!/bin/bash
for i in 1 2 3; do
  pdflatex topological_$i.tex
  convert -density 300 topological_$i.pdf -quality 100 -resize 320 -depth 4 topological_$i.png
done
