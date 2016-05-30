import sys
import os
import json
import operator
from operator import itemgetter

#lettura del file.json
jsonFile = open("full.json", "r")
data = json.load(jsonFile)
jsonFile.close()

#modifica dei dati
def extract_start(json):
    try:
        return int(json['introns'][]['absolute_start'])
    except KeyError:
        return 0

data.sort(key=extract_start, reverse=False)

#apertura del file in scrittura
jsonFile = open("full.json", "w+")
jsonFile.write(json.dumps(data))
jsonFile.close()
