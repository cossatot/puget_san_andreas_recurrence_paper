#!/bin/bash

OUT_TYPE=$1

OUTFILE_NAME="styron_et_al_recurrence_paper."

if [ "$OUT_TYPE" == "doc" ]; then
    pandoc recurrence_manuscript.md --from-markdown \
    -o $OUTFILE_NAME$OUT_TYPE \
    --filter pandoc-citeproc --bibliography=pug_saf_recurrence_manuscript.bib \
    --csl=/Users/itchy/Zotero/styles/american-geophysical-union.csl

elif [ "$OUT_TYPE" == "pdf" ]; then
    pandoc recurrence_manuscript.md \
    -o $OUTFILE_NAME$OUT_TYPE \
    --latex-engine=xelatex \
    --template=default.latex \
    --filter pandoc-citeproc --bibliography=pug_saf_recurrence_manuscript.bib \
    --csl=/Users/itchy/Zotero/styles/american-geophysical-union.csl

elif [ "$OUT_TYPE" == "tex" ]; then
    pandoc recurrence_manuscript.md \
    -o $OUTFILE_NAME$OUT_TYPE \
    --latex-engine=pdflatex \
    --filter pandoc-citeproc --bibliography=pug_saf_recurrence_manuscript.bib \
    --csl=/Users/itchy/Zotero/styles/american-geophysical-union.csl
#--from-markdown \
fi
