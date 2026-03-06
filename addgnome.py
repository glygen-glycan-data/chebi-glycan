#!.venv/bin/python -u

import csv, sys, os
from collections import Counter, defaultdict
from pronto import Ontology
import argparse

parser = argparse.ArgumentParser(description="Add GNOme subusmption relationships")
parser.add_argument('--chebi', required=True, 
                    help = 'ChEBI ontology file (OBO format).')
parser.add_argument('--gnome', required=True,
                    help = 'GNOme ontology file (OBO format).')
parser.add_argument('input',
                    help = 'ChEBI ids with GlyTouCan accessions (TSV format).')
args = parser.parse_args()

from progress import *

with elapsed("Read GNOme ontology"):
    gnome = Ontology(args.gnome)

with elapsed("Read ChEBI ontology"):
    chebi = Ontology(args.chebi,encoding='utf8')

prog = elapsed("Read input")
reader = csv.DictReader(open(args.input),dialect='excel-tab')
headers = None
allchebi = set()
chebi2gtc = defaultdict(set)
for row in reader:
    if not headers:
        headers = reader.fieldnames
    chebiid = row["id"]
    allchebi.add(chebiid)
    gtcacc = row["gtcacc"]
    if not gtcacc:
        continue
    try:
        gtcterm = gnome.get_term("GNO:"+gtcacc)
    except KeyError:
        continue
    chebi2gtc[chebiid].add(gtcterm)
prog.done()

chebirel = defaultdict(set)
for c1 in progress("Compute relationships",chebi2gtc):
    for gtc1term in chebi2gtc[c1]:
        gtc1termsuper = gtc1term.superclasses(with_self=False).to_set()
        for c2 in chebi2gtc:
            for gtc2term in chebi2gtc[c2]:
                if gtc2term in gtc1termsuper:
                    chebirel[c1].add(c2)

with elapsed("Find shortcuts"):
    toremove = set()
    for c1 in chebirel:
        for c2 in chebirel[c1]:
            if c2 not in chebirel:
                continue
            for c3 in chebirel[c2]:
                if c3 in chebirel[c1]:
                    toremove.add((c1,c3))

with elapsed("Remove %d shortcuts"%(len(toremove),)):
    for f,t in toremove:
        chebirel[f].remove(t)

chebirel1 = defaultdict(set)
for c1 in progress("Compute relationships",allchebi):
    c1super = set(t.id for t in chebi.get_term(c1).superclasses(with_self=False))
    for c2 in allchebi & c1super:
        chebirel1[c1].add(c2)

with elapsed("Find shortcuts in ChEBI subsumption"):
    toremove = set()
    for c1 in chebirel1:
        for c2 in chebirel1[c1]:
            if c2 not in chebirel1:
                continue
            for c3 in chebirel1[c2]:
                if c3 in chebirel1[c1]:
                    toremove.add((c1,c3))

with elapsed("Remove %d shortcuts in ChEBI subsumption"%(len(toremove))):
    for f,t in toremove:
        chebirel1[f].remove(t)
        
i = headers.index("level")
headers.insert(i,"chebi_isa")
headers.insert(i,"subsumed_by")

prog = elapsed("Rewrite input...")
print("\t".join(headers))
reader = csv.DictReader(open(args.input),dialect='excel-tab')
for row in reader:
    chebiid = row["id"]
    row['subsumed_by'] = ";".join(sorted(chebirel[chebiid]))
    row['chebi_isa'] = ";".join(sorted(chebirel1[chebiid]))
    print("\t".join(row.get(h,"") for h in headers))
prog.done()





