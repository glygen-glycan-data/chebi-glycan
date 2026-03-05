#!.venv/bin/python

import csv, sys
from collections import Counter

if len(sys.argv) > 2:
    expr = sys.argv[2]
else:
    expr = 'True'

reader = csv.DictReader(open(sys.argv[1]),dialect='excel-tab')
headers = None
for row in reader:
    if not headers:
        headers = ["parent_id","parent_term"] + reader.fieldnames[:-1]
        print("\t".join(headers))
    if not eval(expr,row):
        continue
    parents = row['parents'].split(';')
    row1 = dict(row.items())
    for p in parents:
        pid,pterm = p.split(' - ',1)
        row1['parent_id'] = pid
        row1['parent_term'] = pterm
        print("\t".join(row1.get(h,"") for h in headers))
