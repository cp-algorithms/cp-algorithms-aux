#!/bin/bash
pdflatex sum-segment-tree-update.tex
convert -density 300 sum-segment-tree-update.pdf -quality 100 -resize x380 -depth 4 sum-segment-tree-update.png
