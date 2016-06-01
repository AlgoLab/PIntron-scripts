#!/usr/bin/env python3

import sys
import os
import json
import operator
from operator import itemgetter

#lettura del file.json
jsonFile = open("full.json", "r")
data = json.load(jsonFile)
jsonFile.close()

def sort_introns(intron):
    return intron['chromosome start'] * 1000000 + intron['length']
    # Works only if no intron is longer than 1000000

introns = list(data['introns'].values())
sorted_introns = sorted(introns, key=sort_introns)

data['introns'] = dict(zip(range(1, len(introns)), sorted_introns))
# print(sorted_introns)


print(json.dumps(data, sort_keys=True, indent=4))
