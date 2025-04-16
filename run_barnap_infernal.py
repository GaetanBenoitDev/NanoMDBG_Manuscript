

import os, sys, argparse, shutil, glob, gzip
from Bio.SeqIO.FastaIO import SimpleFastaParser

def main(argv):
    
    parser = argparse.ArgumentParser()

    parser.add_argument("outputDir", help="")
    parser.add_argument("binDir", help="")
    parser.add_argument("nbCores", help="")
    
    args = parser.parse_args()
    

    binDir = args.binDir
    outputDir = args.outputDir
    nbCores = args.nbCores


    outputFilename = outputDir + "/results.txt"
    outputFilenameMAG = outputDir + "/results_mags.tsv"
    outputDir += "/__results/"

    if os.path.exists(outputDir):
        #os.makedirs(outputDir)
        shutil.rmtree(outputDir)

    os.makedirs(outputDir)


    completeMAGs = 0

    resultFileMAG = open(outputFilenameMAG, "w")
    resultFileMAG.write("BinName\t# tRNA\t# 5S\t# 16S\t# 23S\n")

    for binFilename in glob.glob(args.binDir + "/*.fa"):

        print(binFilename)

        barrnapFilename = outputDir + "/results_barrnap.txt"
        infernalFilename = outputDir + "/results_infernal.txt"

        if os.path.exists(barrnapFilename): os.remove(barrnapFilename)
        if os.path.exists(infernalFilename): os.remove(infernalFilename)


        binName = os.path.split(binFilename)[1]
        binName = os.path.splitext(binName)[0]

        command = "conda run -n infernal cmscan --cpu " + args.nbCores + " --cut_ga --rfam --nohmmonly --fmt 2 --tblout " + infernalFilename + " --clanin ~/appa/data/rfam/Rfam.clanin ~/appa/data/rfam/Rfam.cm " + binFilename
        print(command)
        os.system(command)

        command = "conda run -n barrnap barrnap --threads 64 --evalue 0.01 " + binFilename + " > " + barrnapFilename
        print(command)
        os.system(command)


        tRNA = 0
        for line in open(infernalFilename):
            line = line.rstrip()
            if line.startswith("#"): continue

            fields = line.split()

            eValue = float(fields[17])

            #trunc = fields[12]

            if eValue > 0.01: continue

            type = fields[1]
            if(type != "tRNA"): continue

            contigName = fields[3]

            tRNA += 1

        rRNA_5S = 0
        rRNA_16S = 0
        rRNA_23S = 0

        for line in open(barrnapFilename):
            line = line.rstrip()
            if line == "": continue
            if line.startswith("#"): continue

            fields = line.split("\t")

            contigName = fields[0]

            eValue = float(fields[5])

            if eValue > 0.01: continue

            desc = fields[8]
            if "5S" in desc:
                rRNA_5S += 1
            elif "16S" in desc:
                rRNA_16S += 1
            elif "23S" in desc:
                rRNA_23S += 1

        resultFileMAG.write(binName + "\t" + str(tRNA) + "\t" + str(rRNA_5S) + "\t" + str(rRNA_16S) + "\t" + str(rRNA_23S) + "\n")

        if tRNA >= 18 and rRNA_5S > 0 and rRNA_16S > 0 and rRNA_23S > 0:
            completeMAGs += 1

        print(tRNA, rRNA_5S, rRNA_16S, rRNA_23S)

    resultFileMAG.close()

    resultFile = open(outputFilename, "w")
    resultFile.write("BarnapInfernalMags: " + str(completeMAGs) + "\n")
    resultFile.close()

if __name__ == "__main__":
    main(sys.argv[1:])  