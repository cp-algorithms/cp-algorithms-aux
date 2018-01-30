#!/bin/bash
pdflatex sum-segment-tree.tex
convert -density 300 sum-segment-tree.pdf -quality 100 -resize x380 -depth 4 sum-segment-tree.png
