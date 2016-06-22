#!/usr/local/bin/python3.5

import json

#with open('../names', 'r') as f:
#    names = f.read().splitlines()

with open('../hostimgs.json', 'r') as f:
    hostimgs = json.loads(f.read())

count = 0
newhostimgs = []
for ei, entry in enumerate(hostimgs):
    if entry[1] == 'SDSS':
    #if entry[0] not in names:
        count = count + 1
        newhostimgs.append(entry)

print(count)

jsonstring = json.dumps(newhostimgs, indent='\t', separators=(',',':'))
with open('../hostimgs.json', 'w') as f:
    f.write(jsonstring)
