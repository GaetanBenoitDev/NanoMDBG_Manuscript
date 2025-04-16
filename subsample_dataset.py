



import os, sys, argparse, shutil, glob, gzip
from Bio import SeqIO

def main(argv):
    
    parser = argparse.ArgumentParser()

    parser.add_argument("readFilename", help="")
    parser.add_argument("outputDir", help="")
    parser.add_argument("nbCores", help="")
    
    args = parser.parse_args()

    readFilename = args.readFilename
    outputDir = args.outputDir
    nbCores = int(args.nbCores)

    sampleName = os.path.split(readFilename)[1].replace(".fastq.gz", "")

    #outputDir = outputDir + "/" + sampleName + "_subsample/"

    if not os.path.exists(outputDir):
        os.makedirs(outputDir)

    #sampleSizeBps = 50000000000
    #sampleSizeBps = 10000000000
    sampleSizeBps = 100000000000
    
    maxSubsampleSize = 500000000000

    while True:
        
        currentSampleSize = 0
        continueSampling = False


        sizeStr = str(int(maxSubsampleSize / 1000000000)) + "G"
        sampleFilename = outputDir + "/" + sampleName + "_" + str(sizeStr) + ".fastq"
        outputFile = open(sampleFilename, "w")

        print("Creating Subsample: ", maxSubsampleSize, sampleFilename)
        
        inputFile = gzip.open(readFilename, "rt")

        for record in SeqIO.parse(inputFile, "fastq"):

            currentSampleSize += len(record.seq)
            SeqIO.write(record, outputFile, "fastq")

            #print("\t" + str(currentSampleSize))

            if currentSampleSize > maxSubsampleSize:
                continueSampling = True
                maxSubsampleSize += sampleSizeBps
                break
        
        inputFile.close()
        outputFile.close()

        print("Compressing: " + sampleFilename)
        gzipCommand = "pigz -f -p " + str(nbCores) + " " + sampleFilename
        os.system(gzipCommand)

        if not continueSampling: break
        #if maxSubsampleSize == 50000000000: break

if __name__ == "__main__":
    main(sys.argv[1:])  
