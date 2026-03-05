#!.venv/bin/python

import csv, sys
from collections import Counter

if len(sys.argv) > 2:
    expr = sys.argv[2]
else:
    expr = 'True'

freq = Counter()
count = 0
for row in csv.DictReader(open(sys.argv[1]),dialect='excel-tab'):
    if not eval(expr,row):
        continue
    parents = row['parents'].split(';')
    chebiid = row['id']
    for p in parents:
        freq[p] += 1
    count += 1
total = sum(freq.values())

print("id\tterm\t# as parent\t% as parent")
for k,v in sorted(freq.items(),key=lambda t: -t[1]):
    chebiid,term = k.split(' - ',1)
    print(f"{chebiid}\t{term}\t{v}\t{round(100*v/count,2)}")
