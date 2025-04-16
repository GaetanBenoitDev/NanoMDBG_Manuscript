

import os, sys, argparse, shutil, gzip
from Bio.SeqIO.FastaIO import SimpleFastaParser

def main(argv):
    
    parser = argparse.ArgumentParser()

    parser.add_argument("inputFilename", help="contig file")
    parser.add_argument("assembler", help="nanoMDBG, metaMDBG, hifiasm, metaflye")
    parser.add_argument("minLength", help="min contig length")
    parser.add_argument("--circ", help="count circular contigs only", action="store_true")
    parser.add_argument("--checkm2", help="path to checkm2 quality report")
    parser.add_argument("--plasmid", help="path to genomad plasmid annotation file")
    parser.add_argument("--virus", help="path to genomad virus annotation file")
    parser.add_argument("--out", help="output fasta filename")
    
    #parser.add_argument("csv", help="output unitig coverage file (.csv)")
    
    args = parser.parse_args()

    inputFilename = args.inputFilename
    assembler = args.assembler
    minContigLength = int(args.minLength)
    
    collectCircularContigs = False
    if args.circ:
        collectCircularContigs = True
    
    isLongContig, isContigCircular = loadCircularContigs(inputFilename, assembler, minContigLength, collectCircularContigs)

    outputFile = None
    if args.out:
        outputFile = open(args.out, "w")
    
    isNearCompleteContig = None
    isVirusContig = None
    isPlasmidContig = None

    if args.checkm2:
        isNearCompleteContig = loadNearCompleteContigs(args.checkm2)
        #print(contigQuality)
    elif args.virus:
        isVirusContig = loadVirusContigs(args.virus)
    elif args.plasmid:
        isPlasmidContig = loadPlasmidContigs(args.plasmid)

    n = 0

    if isNearCompleteContig is not None:

        if isContigCircular is not None:
            for contigName in isNearCompleteContig:
                if contigName in isContigCircular:
                    n += 1
        else:
            n = len(isNearCompleteContig)

    elif isVirusContig is not None:

        if isContigCircular is not None:
            for contigName in isVirusContig:
                if contigName in isContigCircular:
                    n += 1
        else:
            n = len(isVirusContig)

    elif isPlasmidContig is not None:

        if isContigCircular is not None:
            for contigName in isPlasmidContig:
                if contigName in isContigCircular:
                    n += 1
        else:
            n = len(isPlasmidContig)

    else:

        if isContigCircular:
            for contigName in isLongContig:
                if contigName in isContigCircular:
                    n += 1
        else:
            n = len(isLongContig)

    print(n)


def loadCircularContigs(inputFilename, assembler, minContigLength, collectCircularContigs):

    isCircularContig = set()
    isLongContig = set()

    fileHandle = None

    if(".gz" in inputFilename):
        fileHandle = gzip.open(inputFilename, "rt")
    else:
        fileHandle = open(inputFilename)
        
    if "metaMDBG" in assembler or "nanoMDBG" in assembler:

        for header, seq in SimpleFastaParser(fileHandle):

            contigName, lengthStr, coverageStr, circularStr = header.split(" ")
            
            if circularStr.split("=")[1] == "yes":
                isCircularContig.add(contigName)

            if len(seq) >= minContigLength:
                isLongContig.add(contigName)

    elif "hifiasm" in assembler:
        
        for header, seq in SimpleFastaParser(fileHandle):

            contigName = header.split(" ")[0]

            if header.endswith("c"):
                isCircularContig.add(contigName)

            if len(seq) >= minContigLength:
                isLongContig.add(contigName)

    elif "metaflye" in assembler:
        fileHandle.readline() #skip header

        for line in fileHandle:
            line = line.rstrip()
            if len(line) == 0: continue
            fields = line.split("\t")

            contigName = fields[0]

            isCircular = fields[3] == "Y"
            if isCircular:
                isCircularContig.add(contigName)

            contigLength = int(fields[1])
            if contigLength >= minContigLength:
                isLongContig.add(contigName)

    if collectCircularContigs:
        return isLongContig, isCircularContig
    else:
        return isLongContig, None

"""
def countHifiasm(inputFilename, minContigLength, isNearCompleteContig, isVirusContig, isPlasmidContig, outputFile):

    nbLinearContigs = 0
    nbCircularContigs = 0
    nbNearCompleteContigs = 0

    fileHandle = None

    if(".gz" in inputFilename):
        fileHandle = gzip.open(inputFilename, "rt")
    else:
        fileHandle = open(inputFilename)

    for header, seq in SimpleFastaParser(fileHandle):

        contigName = header.split(" ")[0]

        if isNearCompleteContig is not None:
            if contigName in isNearCompleteContig:



        if len(seq) < minContigLength: continue

        nbLinearContigs += 1
        
        if not header.endswith("c"): continue

        #if header.endswith("rc"): continue
        #if not header.endswith("c"): continue
        nbCircularContigs += 1
        if outputFile:
            outputFile.write(">" + header + "\n")
            outputFile.write(seq + "\n")
              
        #print(len(seq), header, file=sys.stderr)


    return nbCircularContigs, nbLinearContigs

#seq_name	length	cov.	circ.	repeat	mult.	alt_group	graph_path
#contig_5192	6127306	23	N	N	1	*	5192
def countMetaflye(inputFilename, minContigLength, isNearCompleteContig, isVirusContig, isPlasmidContig, outputFile):

    nbLinearContigs = 0
    nbCircularContigs = 0
    
    fileHandle = open(inputFilename)
    fileHandle.readline() #skip header

    for line in fileHandle:
        line = line.rstrip()
        fields = line.split("\t")

        contigLength = int(fields[1])
        if contigLength < minContigLength: continue

        contigName = fields[0]

        if not contigQuality is None:
            headerShorten = contigName.split(" ")[0]
            completeness, contamination = contigQuality[headerShorten]
            if completeness < 90: continue
            if contamination > 5: continue

        nbLinearContigs += 1

        isCircular = fields[3] == "Y"
        if isCircular:
            nbCircularContigs += 1
            #print(completeness, contamination)

    return nbCircularContigs, nbLinearContigs

"""

def loadNearCompleteContigs(checkm2ReportFilename):

    isNearCompleteContig = set()
    #contigQuality = {}

    file = open(checkm2ReportFilename)
    file.readline()

    for line in file:
        line = line.rstrip()
        if len(line) == 0: continue

        fields = line.split("\t")

        contigName = fields[0]
        completeness = float(fields[1])
        contamination = float(fields[2])

        #contigQuality[contigName] = (completeness, contamination)
        if completeness >= 90 and contamination <= 5:
            isNearCompleteContig.add(contigName)

    #print("Nb complete contigs: ", len(completeContigs))

    return isNearCompleteContig

def loadVirusContigs(genomadAnnotationFilename):

    isVirusContig = set()

    #filename = os.path.split(contigFilename)[1].split(".")[0]
    #resultFilename = genomadOutputDir + "/" + filename + "_summary/" + filename + "_virus_summary.tsv"
    
    resultFile = open(genomadAnnotationFilename)
    resultFile.readline()

    for line in resultFile:
        line = line.rstrip()
        if len(line) == 0: continue

        fields = line.split("\t")
        contigName = fields[0]

        isVirusContig.add(contigName)

    return isVirusContig
    
def loadPlasmidContigs(genomadAnnotationFilename):

    isPlasmidContig = set()

    #filename = os.path.split(contigFilename)[1].split(".")[0]
    #resultFilename = genomadOutputDir + "/" + filename + "_summary/" + filename + "_plasmid_summary.tsv"
    
    resultFile = open(genomadAnnotationFilename)
    resultFile.readline()

    for line in resultFile:
        line = line.rstrip()
        if len(line) == 0: continue

        fields = line.split("\t")
        contigName = fields[0]

        isPlasmidContig.add(contigName)

    return isPlasmidContig

if __name__ == "__main__":
    main(sys.argv[1:])  
