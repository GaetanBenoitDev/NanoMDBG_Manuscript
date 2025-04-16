
import os, sys, argparse, time, glob, gzip


exeDirName = os.path.dirname(os.path.realpath(__file__))
from Bio.SeqIO.FastaIO import SimpleFastaParser
#os.chdir(exeDirName)

def main(argv):
    
    parser = argparse.ArgumentParser()

    parser.add_argument("outputDir", help="output dir")
    parser.add_argument("contigs", help="contigs fasta filename")
    parser.add_argument("experimentDir", help="")
    parser.add_argument("minimap2Preset", help="map-ont | map-hifi")
    parser.add_argument("nbCores", help="")
    parser.add_argument("reads", help="")

    args = parser.parse_args()

    outputDir = args.outputDir
    experimentDir = args.experimentDir
    minimap2Preset = args.minimap2Preset
    nbCores = args.nbCores
    outputFilename = outputDir + "/results.txt"

    if os.path.exists(outputDir):
        print("output dir exists")
        #exit(1)
    else:
        os.makedirs(outputDir)

    if(experimentDir[-1] == "/"):
        datasetSize = int(experimentDir.split("/")[-3].replace("G", ""))*1000000000
    else:
        datasetSize = int(experimentDir.split("/")[-2].replace("G", ""))*1000000000
    
    alignFilename = outputDir + "/align.paf"
    command = "minimap2 -c -x " + minimap2Preset + " -t " + str(nbCores) + " " + args.contigs + " " + args.reads  + " > " + alignFilename 
    print(command)
    os.system(command)

    totalBpsMapped = 0
    readsBestHit = {}
    contigNbBpsMapped = {}

    for line in open(alignFilename):
        line = line.rstrip()
        fields = line.split("\t")

        ctgName = fields[5]
        nbMatches = float(fields[9])
        alignLength = float(fields[10])
        readLength = float(fields[1])

        if nbMatches / alignLength < 0.95: continue
        if alignLength < 1000: continue


        if alignLength / readLength < 0.8: continue
        #if alignLength > 2500: continue

        readName = fields[0]
        readLength = int(fields[1])


        if not readName in readsBestHit:
            readsBestHit[readName] = (ctgName, alignLength, readLength)
        else:
            if alignLength > readsBestHit[readName][1]:
                readsBestHit[readName] = (ctgName, alignLength, readLength)

        
    for contigName, alignLength, readLength in readsBestHit.values():
        if not contigName in contigNbBpsMapped:
            contigNbBpsMapped[contigName] = 0

        contigNbBpsMapped[contigName] += readLength
        totalBpsMapped += alignLength


    quality_to_nbBases = {}
    quality_to_nbBases["single-complete"] = 0
    quality_to_nbBases["complete"] = 0
    quality_to_nbBases["high"] = 0
    quality_to_nbBases["med"] = 0

    binningCheckmDir = experimentDir + "/binning/checkm/__checkm/bins_/"
    binIndex_to_contigNames, contigName_to_binIndex = loadContigBin(experimentDir + "/binning/bins/contig_bins.tsv")
    binIndex_to_nbContigs = loadBinNbContigs(experimentDir + "/binning/bins/bins_info.tsv")

    for quality in ["complete", "high", "med"]:

        for binFilename in glob.glob(binningCheckmDir + "/" + quality + "/*.fa"):

            binIndex = int(os.path.split(binFilename)[1].replace("SemiBin_", "").replace(".fa", ""))

            for contigName in binIndex_to_contigNames[binIndex]:

                if not contigName in contigNbBpsMapped: continue

                if quality == "complete" and binIndex_to_nbContigs[binIndex] == 1:
                    quality_to_nbBases["single-complete"] += contigNbBpsMapped[contigName]
                else:
                    quality_to_nbBases[quality] += contigNbBpsMapped[contigName]




    outputFile = open(outputFilename, "w")
    
    for quality, nbBases in quality_to_nbBases.items():
        fractionCovered = float(nbBases) / datasetSize
        outputFile.write(quality + ": " + str(fractionCovered) + "\n")

    outputFile.write("Assembly: " + str(float(totalBpsMapped) / datasetSize) + "\n")

    outputFile.close()

    os.remove(alignFilename)

def loadBinNbContigs(binInfoFilename):

    fileHandle = open(binInfoFilename)
    fileHandle.readline()

    binIndex_to_nbContigs = {}

    for line in fileHandle:
        line = line.rstrip()
        if len(line) == 0: continue

        fields = line.split("\t")

        binFilename = fields[0]
        nbContigs = int(fields[2])

        binIndex = int(os.path.split(binFilename)[1].replace("SemiBin_", "").replace(".fa", ""))

        binIndex_to_nbContigs[binIndex] = nbContigs

    fileHandle.close()

    return binIndex_to_nbContigs

def loadContigBin(contigBinFilename):

    binIndex_to_contigNames = {}
    contigName_to_binIndex = {}

    fileHandle = open(contigBinFilename)
    fileHandle.readline()

    for line in fileHandle:
        line = line.rstrip()
        if len(line) == 0: continue

        fields = line.split("\t")

        contigName = fields[0]
        binIndex = int(fields[1])

        if not binIndex in binIndex_to_contigNames:
            binIndex_to_contigNames[binIndex] = []

        binIndex_to_contigNames[binIndex].append(contigName)

        contigName_to_binIndex[contigName] = binIndex

    fileHandle.close()

    return binIndex_to_contigNames, contigName_to_binIndex

"""
    binIndex_to_contigNames, contigName_to_binIndex = loadContigBin(contigBinFilename)
    contigName_to_coverage = loadContigCoverage(coverageFilename)
    contigName_to_contigLength = loadContigLengths(experimentDir + "/binning/bins/output_bins/")
    binIndex_to_quality = loadSingleContigQuality(singleContigCheckmDir, contigName_to_binIndex)

    if(experimentDir[-1] == "/"):
        datasetSize = int(experimentDir.split("/")[-3].replace("G", ""))*1000000000
        #print(datasetSize)
    else:
        datasetSize = int(experimentDir.split("/")[-2].replace("G", ""))*1000000000
        #print(datasetSize)

    #print(datasetSize)
    #print(binIndex_to_quality)

    binIndex_to_quality = loadBinQuality(binningCheckmDir, binIndex_to_quality)
    #binIndex_to_quality

    #print(binIndex_to_quality)

    quality_to_nbBases = {}
    quality_to_nbBases["single-complete"] = 0
    quality_to_nbBases["complete"] = 0
    quality_to_nbBases["high"] = 0
    quality_to_nbBases["med"] = 0

    for binIndex, contigNames in binIndex_to_contigNames.items():

        if not binIndex in binIndex_to_quality: continue
        
        quality = binIndex_to_quality[binIndex]
        #if not quality in quality_to_nbBases:
        #    quality_to_nbBases[quality] = 0

        for contigName in contigNames:
            coverage = contigName_to_coverage[contigName]
            contigLength = contigName_to_contigLength[contigName]

            nbBases = contigLength * coverage
            quality_to_nbBases[quality] += nbBases


    outputFile = open(outputFilename, "w")
    

    for quality, nbBases in quality_to_nbBases.items():
        fractionCovered = float(nbBases) / datasetSize
        outputFile.write(quality + ": " + str(fractionCovered) + "\n")

    outputFile.close()



def loadContigCoverage(coverageFilename):

    contigName_to_coverage = {}

    fileHandle = open(coverageFilename)
    fileHandle.readline()

    for line in fileHandle:
        line = line.rstrip()
        if len(line) == 0: continue

        fields = line.split(",")

        contigName = fields[0]
        coverage = float(fields[1])

        contigName_to_coverage[contigName] = coverage

    fileHandle.close()

    return contigName_to_coverage

def loadContigLengths(outputBinDir):

    contigName_to_contigLength = {}

    for binFilename in glob.glob(outputBinDir + "/*.fa"):

        f = open(binFilename)

        for header, seq in SimpleFastaParser(f):
            contigName = header.split(" ")[0]
            contigName_to_contigLength[contigName] = len(seq)

        f.close()

    return contigName_to_contigLength

def loadSingleContigQuality(checkmDir, contigName_to_binIndex):

    binIndex_to_quality = {}

    for quality in ["complete"]:

        for binFilename in glob.glob(checkmDir + "/" + quality + "/*.fa"):

            contigName = os.path.split(binFilename)[1].replace(".fa", "")
            binIndex = contigName_to_binIndex[contigName]
            binIndex_to_quality[binIndex] = "single-" + quality

    return binIndex_to_quality


"""

if __name__ == "__main__":
    main(sys.argv[1:])  

















