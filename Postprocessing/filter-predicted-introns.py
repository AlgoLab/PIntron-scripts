# -*- coding: utf-8 -*-
import sys, math
import argparse
import time
import os

def main():

#creo/apro il file di LOG, che verra aperto/creato nella cartella dove si sta lavorando

    log=open('LOG.txt','a')
    logTime = time.localtime(time.time())
    FormatlogTime = time.strftime("%d/%m/%Y %H:%M:%S", logTime)
    log.write("**************************\n")
    log.write(FormatlogTime +'\n')

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
        print >> sys.stderr, "ERROR: No predicted-introns file given.\nAborting..."
        log.write('predicted-introns: ERROR: No predicted-introns file given.\n')
        sys.exit( 1 )

    abspfile=os.path.abspath(args.pfile)
    log.write('predicted-introns: '+abspfile+ '\n')

    if not args.wfile:
        print >> sys.stderr, "ERROR: No output file specified.\nAborting..."
        log.write('output-file: ERROR: No predicted-introns file given.\n')
        sys.exit( 1 )

    abswfile=os.path.abspath(args.wfile)
    log.write('output file: ' +abswfile+ '\n')

# controllo su valori di taglio inseriti dall'utente (se non specificati sono N=20 e M=N/2=10)'

    if args.minimum > args.nValue :
        print >> sys.stderr, "ERROR: Minimum value is greater than N the filter value .\nAborting..."
        log.write('filtering threshold: ERROR: Minimum value is greater than N the filter value.')
        sys.exit( 1 )

    log.write('filtering threshold: %d \n' %args.nValue)

# inizializzo le variabili con i valori presi dal parser e apro i file input/output

    N = args.nValue
    if not args.minimum:
	M = N/2
    else: 
    	M = args.minimum
    log.write('minimum reads: %d \n' %M)
    f=open( args.pfile,'r')
    out=open( args.wfile,'w')

# legge la linea e la trasforma in un lista di caratteri indicizzati [0,1...N] con l'istruzione .split'
   
    # CHECK corretto inserimento valori di taglio
    # print (N, M)

    while 1:
        linea = f.readline()
        if len(linea) == 0:
            break
        lista = linea.split()

# fa il filtraggio dei dati

        nReads = int(lista[5])
        intron_type = int(lista[13])
        if (nReads>=N) or \
		( \
		(nReads>=M) and \
		(intron_type==0 or intron_type==1) and \
		(lista[14]=="GTAG" or lista[14]=="GCAG" or lista[14]=="ATAC") \
		):
            # CHECK controllo dei parametri filtrati, 
            # print "{0}, {1}, {2}".format(lista[5], lista[13], lista[14])
            out.write(linea)	

if __name__ == "__main__":
    main()
