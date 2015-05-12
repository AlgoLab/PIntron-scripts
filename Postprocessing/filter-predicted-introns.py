#!/usr/bin/env python2.7

import logging
import sys
import argparse
import os
import json
import gzip

def main():
    logging.basicConfig(level=logging.INFO)
    logging.info('Program started.')

    #creazione del parser di input da linea di comando
    parser = argparse.ArgumentParser(prog = "predicted-introns-filtered",
                                      description = "Filtering intron from predicted-introns",
                                      formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-p', '--predicted-introns', help = "Predicted-introns file in JSON format.",
                         required = False, dest = 'pfile')
    parser.add_argument('-o', '--output-file', help = "Name of the output file in JSON format.",
                         required = False, dest = 'wfile')
    parser.add_argument('-n', '--filter-value', help = ' Filter value',
                         required = False, dest = 'nValue', type = int, default = 20)
    parser.add_argument('-m', '--minimum-canonic-reads', help = ' Minimum canonic reads and',
                         required = False, dest = 'minimum', type = int)
    args = parser.parse_args()

    # controllo che il predicted-introns e il file di output sia specificato dall'utente
    # sia che il controllo fallisca, sia che il controllo vada a buon fine scrive sul file di LOG
    # In entrambi i casi scrive l'indirizzo assoluto del file sul file di LOG'
    if not args.pfile :
        logging.error('No predicted-introns file given.\nAborting...')
        sys.exit(1)

    if not (args.pfile.endswith(".json") or args.pfile.endswith(".gz")):
        logging.error('Input file name must be a JSON (or GZ).\nAborting...')
        sys.exit(1)

    logging.info('Predicted-introns file: ' + args.pfile)

    if not args.wfile:
        logging.error('No output file specified.\nAborting...')
        sys.exit(1)

    if not (args.wfile.endswith(".json") or args.wfile.endswith(".gz")):
        logging.error('Output file name must be a JSON (or GZ).\nAborting...')
        sys.exit(1)

    logging.info('Output file: ' + args.wfile)

    # controllo su valori di taglio inseriti dall'utente (se non specificati sono N=20 e M=N/2=10)'
    if args.minimum > args.nValue :
        logging.error('Minimum value is greater than N the filter value .\nAborting...')
        sys.exit(1)

    logging.info('Filter threshold: ' + str(args.nValue))

    # inizializzo le variabili con i valori presi dal parser e apro i file input/output
    N = args.nValue
    M = 0
    if not args.minimum:
	M = N / 2
    else: 
    	M = args.minimum
    logging.info('Minimum reads: ' + str(M))

    jdata = {}
    if args.pfile.endswith(".gz"):
       with gzip.open(args.pfile, 'rb') as jfile:
           jdata = json.loads(jfile.read().decode("ascii"))
    else:
        with open(args.pfile, 'r') as jfile:
            jdata = json.load(jfile)
    # print jdata

    if args.wfile.endswith(".gz"):
        out = gzip.open(args.wfile, 'wb')
    else:
        out = open(args.wfile, 'w')

    num_kept_introns = 0
    num_tot_introns = 0
    introns = jdata['introns']
    new_introns = []
    for i in introns:
        num_tot_introns += 1
        num_reads = int(i['number_of_supporting_transcripts'])
        intron_type = int(i['type'])
        intron_pattern = i['pattern']
        if (num_reads >= N) or \
           ( \
             (num_reads >= M) and \
             (intron_type == 0 or intron_type == 1) and \
             (intron_pattern == "GTAG" or intron_pattern == "GCAG" or intron_pattern == "ATAC") \
         ):
            # Controllo dei parametri filtrati
            logging.debug(str(num_reads) + " " + str(intron_type) + " " + intron_pattern)
            num_kept_introns += 1
            new_introns.append(i)

    logging.info("Num. input introns: " + str(num_tot_introns))
    logging.info("Num. kept introns: " + str(num_kept_introns))
    jdata['introns'] = new_introns
    json.dump(jdata, out)
    
    logging.info('Program completed.')

if __name__ == "__main__":
    main()
