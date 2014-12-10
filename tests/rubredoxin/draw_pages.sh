#!/bin/sh
for i in page??.gnuplot; do gnuplot $i ; echo $i ; done
gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=plots.pdf -dBATCH page??.pdf

# Clean up
rm page??.pdf
