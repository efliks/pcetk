#!/bin/sh

# Render pages
for i in page??.gnuplot; do
    gnuplot $i
    echo $i
done

# Unite pages
pdfunite page??.pdf plots.pdf

# Clean up
rm page??.pdf
