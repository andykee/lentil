#!/bin/sh

# note: this requires poppler to be installed (brew install poppler)
xelatex logo-dark.tex
pdftocairo -svg logo-dark.pdf logo-dark.svg

xelatex logo-light.tex
pdftocairo -svg logo-light.pdf logo-light.svg

rm *.pdf
rm *.log
rm *.aux