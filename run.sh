#!/bin/sh
./extract_xrefs.py --chebi chebi.obo.gz --gnome GNOme.obo --pubchem gtc2chebi.pubchem.txt --kegg kegg.txt > chebi-togtc.tsv
./addgnome.py --gnome GNOme.obo --chebi chebi.obo.gz chebi-togtc.tsv > chebi-togtc-gnome.tsv
./expand_parents.py chebi-togtc-gnome.tsv > chebi-togtc-gnome-parents-first.tsv
