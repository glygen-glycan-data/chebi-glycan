#!.venv/bin/python -u

import gzip, sys, os, re
from collections import defaultdict
import argparse
from pronto import Ontology, LiteralPropertyValue, ResourcePropertyValue

from progress import *

from pygly.GlycanResource import GlyTouCan
from pygly.GlycanMultiParser import GlycanMultiParser, GlycanParseError

parser = argparse.ArgumentParser(description="Analyze ChEBI Ontology Glycans")

parser.add_argument('--chebi', required=True, 
                    help = 'ChEBI ontology file (OBO format).')
parser.add_argument('--gnome', required=True,
                    help = 'GNOme ontology file (OBO format).')
parser.add_argument('--pubchem', required=True,
                    help = 'Mapping from ChEBI to GlyTouCan from PubChem.')
parser.add_argument('--kegg', required=True,
                    help = 'Mapping from KEGG to GlyTouCan from GlyTouCan.')
args = parser.parse_args()

prog = elapsed("Read PubChem-based ChEBI to GTC mapping")
pubchem_based = defaultdict(set)
for l in open(args.pubchem,'rt'):
    sl = l.split()
    pubchem_based[sl[0]].add(sl[1])
prog.done()

prog = elapsed("Read KEGG to GTC mapping")
kegg_togtc = defaultdict(set)
for l in open(args.kegg,'rt'):
    sl = l.split()
    kegg_togtc[sl[0]].add(sl[1])
prog.done()

with elapsed("Read GNOme ontology"):
    gnome = Ontology(args.gnome)

levelmap = {"GNO:00000013": "basecomposition",
            "GNO:00000014": "composition",
            "GNO:00000015": "topology",
            "GNO:00000016": "saccharide"}

headers = """
id term stars gtcacc monocount basecomposition composition level score fully_determined name glytoucan glygen
pubchem kegg glylookup_wurcs glylookup_iupac kegg.glycan wurcs iupac
""".split()

gtcaccheaders = """
name glytoucan glygen pubchem kegg glylookup_wurcs glylookup_iupac
""".split()

from glyomicsclient import GlyLookupClient
glylookup = GlyLookupClient()
gtc = GlyTouCan()
multiparser = GlycanMultiParser()

with elapsed("Read ChEBI ontology"):
    chebi = Ontology(args.chebi,encoding='utf8')

chebi_class = dict(
    carbohydrate_and_derivitives = chebi.get_term("CHEBI:78616"), #carbohydrates and carbohydrate derivatives
    carbohydrate = chebi.get_term("CHEBI:16646"), #carbohydrate
    polysaccharide = chebi.get_term("CHEBI:18154"), #polysaccharide
    oligosaccharide = chebi.get_term("CHEBI:50699"), #oligosaccharide
    glycan = chebi.get_term("CHEBI:167559"), #glycan
    partially_defined_glycan = chebi.get_term("CHEBI:146306"), #partially-defined glycan
    topology_glycan = chebi.get_term("CHEBI:167503"), #topology glycan
    composition_glycan = chebi.get_term("CHEBI:167502"), #composition glycan
    basecomposition_glycan = chebi.get_term("CHEBI:167481"), #basecomposition glycan
)

def distance(frm,to):
    try:
        frmsuper = frm.superclasses()
    except AttributeError:
        return -1
    if to not in frmsuper:
        return -1
    step = 0
    while True:
        if to in frm.superclasses(distance=step):
            return step
        step += 1
    return -1

headers.extend(chebi_class.keys())
headers.append("parents")
print("\t".join(headers))

for term,entry in progress("ChEBI terms",chebi.items()):
    # print(term,file=sys.stderr)

    d = dict(id=term,term=(entry.name if entry.name else ""))
    for ss in entry.subsets:
        if ss.endswith(":STAR"):
            d['stars'] = ss.split(':',1)[0]
    
    if term in pubchem_based:
        d['pubchem'] = ",".join(sorted(pubchem_based[term]))

    if entry.name:
        split_name = entry.name.split()
        if split_name[0] == "GlyTouCan":
            d['name'] = split_name[1]

    for xref in entry.xrefs:
        split_xref = xref.id.split(":",1)
        if split_xref[0] in ("glygen","glytoucan","kegg.glycan"):
            if split_xref[0] not in d:
                d[split_xref[0]] = set()
            d[split_xref[0]].add(split_xref[1])

    if d.get('kegg.glycan'):
        for kg in d['kegg.glycan']:
            if kg in kegg_togtc:
                if "kegg" not in d:
                    d["kegg"] = set()
                d['kegg'].update(kegg_togtc[kg])
    for k in ("glygen","glytoucan","kegg.glycan","kegg"):
        if d.get(k):
            d[k] = ",".join(sorted(d[k]))
        elif k in d:
            del d[k]

    for ann in entry.annotations:
        if ann.property == "chemrof:wurcs_representation":
            d['wurcs'] = ann.literal
            wurcs_acc = glylookup.get_accession_for_sequence(d['wurcs'])
            if wurcs_acc:
                d['glylookup_wurcs'] = wurcs_acc

    iupac = None
    if entry.name:
        skip = False
        try:
            multiparser.toGlycan(entry.name)
        except GlycanParseError:
            skip = True
        if not skip:
            d['iupac'] = entry.name
            iupac = glylookup.get_accession_for_sequence(entry.name)
            if iupac:
                d['glylookup_iupac'] = iupac

    if not iupac:
        for syn in entry.synonyms:
            # print(syn.scope,syn.type.id if syn.type else None,syn.description)
            if syn.scope == "EXACT" and syn.type and syn.type.id == "IUPAC:NAME":
                try:
                    multiparser.toGlycan(syn.description)
                except GlycanParseError:
                    continue
                d['iupac'] = syn.description
                iupac = glylookup.get_accession_for_sequence(syn.description)
                if iupac:
                    d['glylookup_iupac'] = iupac
                    break

    nglycocols = len(set(d.keys()).difference(set(["id","term","stars"])))
    undercarb = (distance(entry,chebi_class['carbohydrate_and_derivitives']) >= 0)
    
    if nglycocols == 0 and not undercarb:
        continue

    for cls in chebi_class:
        d[cls] = str(distance(entry,chebi_class[cls]))
    d["parents"] = ";".join(sorted([ "%s - %s"%(t.id,t.name) for t in entry.superclasses(distance=1,with_self=False)]))

    gtcaccs = set(filter(None,[d.get(h) for h in gtcaccheaders ]))
    gtcaccs1 = set(filter(None,[d.get(h) for h in filter(lambda h: 'iupac' not in h, gtcaccheaders) ]))
    gtcacc = None
    if len(gtcaccs) == 1:
        gtcacc = list(gtcaccs)[0]
    elif len(gtcaccs1) == 1:
        gtcacc = list(gtcaccs1)[0]
    if gtcacc and ',' not in gtcacc:
        d['gtcacc'] = gtcacc
    
    if d.get('gtcacc'):
        try:
            gnoentry = gnome.get_term("GNO:"+d['gtcacc'])
        except KeyError:
            gnoentry = None
        if gnoentry:
            for ann in gnoentry.annotations:
                if ann.property == "GNO:00000102": #score
                    d['score'] = ann.literal
                elif ann.property == "GNO:00000021": #level
                    d['level'] = levelmap[ann.resource]

        gly = gtc.getGlycan(d['gtcacc'])
        if gly:
            if gly.fully_determined():
                d['fully_determined'] = "Y"
            else:
                d['fully_determined'] = "N"
            if not gly.repeated():
                comp = gly.iupac_composition(aggregate_basecomposition=False)
                d['monocount'] = str(comp['Count'])
                d['composition'] = gly.composition_string()
                d['basecomposition'] = gly.basecomposition_string()
    
    print("\t".join(map(lambda h: d.get(h,""),headers)))
    sys.stdout.flush()
