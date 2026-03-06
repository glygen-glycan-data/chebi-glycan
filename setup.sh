#!/bin/sh
set -x
if [ ! -d .venv ]; then
    ${PYTHON:-python3.12} -m venv .venv
    .venv/bin/python -m pip install -r requirements.txt
fi
if [ ! -f chebi.obo.gz ]; then
    wget -O chebi.obo.gz 'https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.obo.gz'
fi
if [ ! -f GNOme.obo ]; then
    wget -O GNOme.obo 'https://purl.obolibrary.org/obo/gno.obo'
fi
if [ ! -d ${PYGLY:-../PyGly} ]; then
    echo "Please checkout the glygen-glycan-data/PyGly repository from GitHub in `dirname ${PWD}` or set PYGLY appropriately" 1>&2
    exit 1
fi
if [ ! -d pygly ]; then
    ln -s ${PYGLY:-../PyGly}/pygly
fi
if [ ! -f kegg.txt ]; then
    ${PYGLY:-../PyGly}/scripts/glyres.py GlyTouCanNoCache allcrossrefs | fgrep kegg_glycan | awk '{print $3,$1}' > kegg.txt
fi
if [ ! -f gtc2chebi.pubchem.txt ]; then
    ${PYGLY:-../PyGly}/scripts/glyres.py PubChemDownload allchebigtc > gtc2chebi.pubchem.txt
fi


