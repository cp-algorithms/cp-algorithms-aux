#!/bin/bash
for f in $(ls *.tex | cut -f 1 -d '.');
do
    pdflatex $f.tex
    convert -density 300 $f.pdf -quality 100 -resize 500 -depth 4 $f.png
done
