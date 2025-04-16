
import os, sys, argparse, time, glob, gzip


exeDirName = os.path.dirname(os.path.realpath(__file__))
from Bio.SeqIO.FastaIO import SimpleFastaParser
#os.chdir(exeDirName)

def main(argv):
    
    parser = argparse.ArgumentParser()

    parser.add_argument("outputDir", help="output dir")
    parser.add_argument("contigs", help="contigs fasta filename")
    parser.add_argument("minimap2Preset", help="map-ont | map-hifi")
    parser.add_argument("nbCores", help="")
    parser.add_argument("reads", nargs='+', help="Read filename(s) (separated by space)")
    
    args = parser.parse_args()

    outputDir = args.outputDir
    #checkm2ReportFilename = args.checkm2ReportFilename
    contigFilename = args.contigs
    minimap2Preset = args.minimap2Preset
    nbCores = int(args.nbCores)
    binDir = outputDir + "/bins/"
    genomadOutputDir = outputDir + "/genomad/"

    readFilenames = []
    for path in args.reads:
        readFilenames.append(path)

    if os.path.exists(outputDir):
        print("output dir exists")
        #exit(1)
    else:
        os.makedirs(outputDir)


    allBamFilename = []

    for readFilename in readFilenames:

        filename = os.path.basename(readFilename)
        outputFilename = outputDir + "/" + filename + "_mapping.bam"

        #command = "conda run -n minimap2 python3 " + exeDirName + "/mapReadsJob.py " + readFilename + " " + contigFilenameFragmented + " " + outputFilename + " " + str(nbCores)
        #print(command)
        command = "minimap2 --split-prefix " + outputFilename+".split.tmp" + " --secondary=no --sam-hit-only -t " + str(nbCores) + " -a -x " + minimap2Preset + " " + contigFilename + " " + readFilename + " | samtools sort -@ " + str(nbCores) + " -o " + outputFilename
        print(command)

        allBamFilename.append(outputFilename)


        os.system(command)


    command = "conda run -n semibin2 SemiBin2 single_easy_bin --random-seed 42 --sequencing-type=long_read -p " + str(nbCores) + " -i " + contigFilename + " -b " + outputDir + "/*.bam -o " + binDir + " --self-supervised --compression=none --tmpdir " + binDir + "/tmp/"
    print(command)
    os.system(command)
    
    binDir = binDir + "/output_bins/"

    command = "python3 " + exeDirName + "/checkm3.py " + binDir + " " + outputDir + "/checkm" + " " + str(nbCores)
    print(command)
    os.system(command)


    #command = "genomad end-to-end --splits 8 --cleanup " + contigFilename + " " + genomadOutputDir + " ~/appa/data/genomad/genomad_db/ --threads " + str(nbCores)
    #print(command)
    #os.system(command)


    
    for filename in allBamFilename:
        if os.path.exists(filename): os.remove(filename)


if __name__ == "__main__":
    main(sys.argv[1:])  
