



import os, sys, argparse, shutil, glob, gzip
from Bio.SeqIO.FastaIO import SimpleFastaParser



def main(argv):

    parser = argparse.ArgumentParser()

    parser.add_argument("outputDir", help="")
    parser.add_argument("assemblyFilename", help="contig file")
    parser.add_argument("minContigLength", help="min contig length")
    parser.add_argument("nbCores", help="")
    
    args = parser.parse_args()

    outputDir = args.outputDir + "/checkm/"
    nbCores = int(args.nbCores)

    if not os.path.exists(outputDir):
        os.makedirs(outputDir)

    binDir = outputDir + "/bins/"
    if not os.path.exists(binDir):
        os.makedirs(binDir)
    
    createLongContigBins(binDir, args.assemblyFilename, int(args.minContigLength))

    exeDirName = os.path.dirname(os.path.realpath(__file__))

    command = "python3 " + exeDirName + "/checkm3.py " + binDir + " " + outputDir + " " + str(nbCores)
    #print(command)
    ret = os.system(command)
    if(ret != 0): sys.exit(1)


def createLongContigBins(binDir, contigFilename, minContigLength):


    fileHandle = None
    if(".gz" in contigFilename):
        fileHandle = gzip.open(contigFilename, "rt")
    else:
        fileHandle = open(contigFilename)

    for header, seq in SimpleFastaParser(fileHandle):
        if len(seq) < minContigLength: continue

        if " " in header: header = header.split(" ")[0]


        binFileCirc = open(binDir + "/" + header + ".fa", "w")

        binFileCirc.write(">" + header + "\n")
        binFileCirc.write(seq + "\n")

        binFileCirc.close()




if __name__ == "__main__":
    main(sys.argv[1:])  
