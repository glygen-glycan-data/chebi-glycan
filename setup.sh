#!/bin/sh

${PYTHON:-python3.12} -m venv .venv
.venv/bin/python -m pip -r requirements.txt
ln -s ../PyGly/pygly
ln -s ../PyGly/scripts/glyres.py
./glyres.py GlyTouCanNoCache allcrossrefs | fgrep kegg > kegg.txt
./glyres.py PubChemDownload allchebigtc > gtc2chebi.pubchem.txt


