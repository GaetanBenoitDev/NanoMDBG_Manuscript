
import os, sys, argparse, shutil, glob, gzip
from Bio.SeqIO.FastaIO import SimpleFastaParser



def main(argv):


    parser = argparse.ArgumentParser()

    parser.add_argument("outputDir", help="")
    parser.add_argument("contigFilename", help="")
    parser.add_argument("contigInfoFilename", help="")
    parser.add_argument("assembler", help="")
    parser.add_argument("genomadSummaryFilename", help="")
    parser.add_argument("nbCores", help="")

    
    args = parser.parse_args()
    
    outputDir = args.outputDir
    contigFilename = args.contigFilename
    genomadSummaryFilename = args.genomadSummaryFilename
    nbCores = args.nbCores

    isCircularContig = loadCircularContigs(args.contigInfoFilename, args.assembler)

    outputFilename = outputDir + "/results.txt"
    outputDir += "/__results/"

    if os.path.exists(outputDir):
        #os.makedirs(outputDir)
        shutil.rmtree(outputDir)

    os.makedirs(outputDir)




    isVirusContig = getVirusContigs(genomadSummaryFilename)

    fileHandle = None

    if(".gz" in contigFilename):
        fileHandle = gzip.open(contigFilename, "rt")
    else:
        fileHandle = open(contigFilename)

    virusContigFilename = outputDir + "/virusContigs.fasta"
    virusContigFile  = open(virusContigFilename, "w")

    for header, seq in SimpleFastaParser(fileHandle):

        header = header.split(" ")[0]

        if header in isVirusContig and header in isCircularContig:
            virusContigFile.write(">" + header + "\n")
            virusContigFile.write(seq + "\n")

    virusContigFile.close()

    print("Nb virus predictions: ", len(isVirusContig))

    command = "conda run -n checkv checkv end_to_end " + virusContigFilename + " " + outputDir + " -t " + str(args.nbCores) + " -d ~/appa/data/checkv/checkv-db-v1.5/" 
    print(command)
    os.system(command)


    nbVirusPerQuality = {}

    checkvFilename = outputDir + "/quality_summary.tsv"
    checkvFile = open(checkvFilename)
    checkvFile.readline()


    qualities = ["Complete", "High-quality", "Medium-quality", "Low-quality", "Not-determined"]
    for quality in qualities:
        #nbVirusPerQuality[quality + "_linear"] = 0
        #nbVirusPerQuality[quality + "_circular"] = 0
        nbVirusPerQuality[quality] = 0
        
    for line in checkvFile:
        line = line.rstrip()
        if len(line) == 0: continue

        fields = line.split("\t")
        contigName = fields[0]
        checkvQuality = fields[7]

        if not contigName in isCircularContig: continue
        #if contigName in isCircularContigs:
        #    checkvQuality += "_circular"
        #else:
        #    checkvQuality += "_linear"
        #if not checkvQuality in nbVirusPerQuality:
        #    nbVirusPerQuality[checkvQuality] = 0
            
        nbVirusPerQuality[checkvQuality] += 1

    resultFile = open(outputFilename, "w")

    #print("Linear:")
    #resultFile.write("Linear:\n")

    for quality in qualities:
        #print("\t" + quality + ":", nbVirusPerQuality[quality + "_linear"])
        #resultFile.write("\t" + quality + ": " + str(nbVirusPerQuality[quality+ "_linear"]) + "\n")
        resultFile.write("Checkv_" + quality + ": " + str(nbVirusPerQuality[quality]) + "\n")

    #print("Circular:")
    #checkvResultFile.write("Circular:\n")

    #for quality in qualities:
    #    print("\t" + quality + ":", nbVirusPerQuality[quality + "_circular"])
    #    checkvResultFile.write("\t" + quality + ": " + str(nbVirusPerQuality[quality+ "_circular"]) + "\n")

    #print("")

    #print("> High Quality: ", nbVirusPerQuality["Complete"] + nbVirusPerQuality["High-quality"])
    #resultFile.write("High Quality linear: " + str(nbVirusPerQuality["Complete_linear"] + nbVirusPerQuality["High-quality_linear"]) + "\n")
    #resultFile.write("High Quality circular: " + str(nbVirusPerQuality["Complete_circular"] + nbVirusPerQuality["High-quality_circular"]) + "\n")
    resultFile.write("High Quality: " + str(nbVirusPerQuality["Complete"] + nbVirusPerQuality["High-quality"]) + "\n")
    resultFile.close()

def getVirusContigs(genomadAnnotationFilename):

    isVirusContig = set()
    
    resultFile = open(genomadAnnotationFilename)
    resultFile.readline()

    for line in resultFile:
        line = line.rstrip()
        if len(line) == 0: continue

        fields = line.split("\t")
        contigName = fields[0]

        isVirusContig.add(contigName)

    return isVirusContig

def loadCircularContigs(inputFilename, assembler):

    isCircularContig = set()

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

    elif "hifiasm" in assembler:
        
        for header, seq in SimpleFastaParser(fileHandle):

            contigName = header.split(" ")[0]

            if header.endswith("c"):
                isCircularContig.add(contigName)

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

    return isCircularContig
    
if __name__ == "__main__":
    main(sys.argv[1:])  
