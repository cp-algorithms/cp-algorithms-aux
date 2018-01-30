#!/bin/bash
pdflatex sum-segment-tree-query.tex
convert -density 300 sum-segment-tree-query.pdf -quality 100 -resize x380 -depth 4 sum-segment-tree-query.png
