#!/bin/bash
for i in SA*.tex;
do
    pdflatex $i
    convert -density 300 ${i%.tex}.pdf -quality 100 -resize 50% -depth 4 ${i%.tex}.png
done
convert -density 300 SA_suffix_links.pdf -quality 100 -resize 42% -depth 4 SA_suffix_links.png
