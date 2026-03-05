#!/bin/sh

${PYTHON:-python3.12} -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
ln -s ../PyGly/pygly
ln -s ../PyGly/scripts/glyres.py
./glyres.py GlyTouCanNoCache allcrossrefs | fgrep kegg_glycan | awk '{print $3,$1}' > kegg.txt
./glyres.py PubChemDownload allchebigtc | awk '{print $2,$1}'> gtc2chebi.pubchem.txt


