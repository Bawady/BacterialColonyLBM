#!/bin/sh

gnuplot< plot.gnu
pdflatex All
pdflatex Nutrients

rm *.eps
rm *.tex
rm *.log
rm *-inc-eps-converted-to.pdf
rm *.aux

pdfcrop All.pdf
pdfcrop Nutrients.pdf
rm All.pdf
rm Nutrients.pdf
mv All-crop.pdf All.pdf
mv Nutrients-crop.pdf Nutrients.pdf