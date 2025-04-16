

import os, sys, argparse, shutil, gzip
from Bio.SeqIO.FastaIO import SimpleFastaParser

def main(argv):
    
    parser = argparse.ArgumentParser()

    parser.add_argument("contigs", help="contig file")
    
    #parser.add_argument("csv", help="output unitig coverage file (.csv)")
    
    args = parser.parse_args()

    contigFilename = args.contigs

    if(".gz" in contigFilename):
        fileHandle = gzip.open(contigFilename, "rt")
    else:
        fileHandle = open(contigFilename)

    nbBases = 0

    for header, seq in SimpleFastaParser(fileHandle):
        nbBases += len(seq)

    print(nbBases)

if __name__ == "__main__":
    main(sys.argv[1:])  