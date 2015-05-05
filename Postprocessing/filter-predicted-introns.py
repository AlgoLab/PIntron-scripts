#!/usr/bin/env python2.7

import logging
import sys, math
import argparse
import time
import os

def main():
    logging.basicConfig( level=logging.INFO )
    logging.info( 'Program started.' )

    #creazione del parser di input da linea di comando
    parser = argparse.ArgumentParser( prog = "predicted-introns-filtered",
                                      description = "Filtering intron from predicted-introns",
                                      formatter_class = argparse.ArgumentDefaultsHelpFormatter )
    parser.add_argument( '-p', '--predicted-introns', help = " Predicted-introns file",
                         required = False, dest = 'pfile' )
    parser.add_argument( '-o', '--where-save-output', help = ' Name of file where save the output',
                         required = False, dest = 'wfile')
    parser.add_argument( '-n', '--filter-value', help = ' Filter value',
                         required = False, dest = 'nValue', type = int, default = 20 )
    parser.add_argument( '-m', '--minimum-canonic-reads', help = ' Minimum canonic reads and',
                         required = False, dest = 'minimum', type = int )
    args = parser.parse_args(  )

    # controllo che il predicted-introns e il file di output sia specificato dall'utente
    # sia che il controllo fallisca, sia che il controllo vada a buon fine scrive sul file di LOG
    # In entrambi i casi scrive l'indirizzo assoluto del file sul file di LOG'
    if not args.pfile :
        logging.error( 'No predicted-introns file given.\nAborting...' )
        sys.exit( 1 )

    abspfile = os.path.abspath( args.pfile )
    logging.info( 'Predicted-introns file: ' + abspfile )

    if not args.wfile:
        logging.error( 'ERROR: No output file specified.\nAborting...' )
        sys.exit( 1 )

    abswfile = os.path.abspath( args.wfile )
    logging.info('Output file: ' + abswfile )

    # controllo su valori di taglio inseriti dall'utente (se non specificati sono N=20 e M=N/2=10)'
    if args.minimum > args.nValue :
        logging.error( 'ERROR: Minimum value is greater than N the filter value .\nAborting...' )
        sys.exit( 1 )

    logging.info( 'Filter threshold: ' + str( args.nValue ) )

    # inizializzo le variabili con i valori presi dal parser e apro i file input/output
    N = args.nValue
    M = 0
    if not args.minimum:
	M = N / 2
    else: 
    	M = args.minimum
    logging.info( 'Minimum reads: ' + str( M ) )
    
    with open( args.pfile, 'r' ) as f, open( args.wfile, 'w' ) as out:
        # legge la linea e la trasforma in un lista di caratteri indicizzati [0,1...N] con l'istruzione .split'
   
        # CHECK corretto inserimento valori di taglio
        logging.debug( "N: " + str( N ) + " M: " + str( M ) )
        for linea in f.readlines( ):
            lista = linea.split()

            # fa il filtraggio dei dati
            nReads = int( lista[5] )
            intron_type = int( lista[13] )
            if ( nReads >= N ) or \
               ( \
                 ( nReads >= M ) and \
                 ( intron_type == 0 or intron_type == 1 ) and \
                 ( lista[14] == "GTAG" or lista[14] == "GCAG" or lista[14] == "ATAC" ) \
             ):
                # CHECK controllo dei parametri filtrati,
                logging.debug( lista[5] + " " + lista[13] + " " + lista[14] )
                out.write( linea )
    logging.info( 'Program completed.' )

if __name__ == "__main__":
    main()
