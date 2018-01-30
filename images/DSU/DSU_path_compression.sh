#!/bin/bash
pdflatex DSU_path_compression.tex
convert -density 300 DSU_path_compression.pdf -quality 100 -resize 732 -depth 4 DSU_path_compression.png
