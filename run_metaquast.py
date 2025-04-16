




import os, sys, argparse, glob





def main(argv):
    
    parser = argparse.ArgumentParser()

    parser.add_argument("contigFilename", help="")
    parser.add_argument("refFilename", help="")
    parser.add_argument("outDir", help="")
    parser.add_argument("nbCores", help="")

    
    args = parser.parse_args()

    contigFilename = args.contigFilename
    refFilename = args.refFilename
    outputDir = args.outDir

    referenceStr = ""

    for line in open(refFilename):
        line = line.rstrip()
        if len(line) == 0: continue

        if not os.path.exists(line): continue

        if(len(referenceStr) == 0):
            referenceStr += line
        else:
            referenceStr += "," + line

        
        
    #for filename in glob.glob(referenceDir + "/*.fa") + glob.glob(referenceDir + "/*.fna") + glob.glob(referenceDir + "/*.fasta"):
    #    if(len(referenceStr) == 0):
    #        referenceStr += filename
    #    else:
    #        referenceStr += "," + filename

    print(referenceStr)

    #--min-contig 50000 --no-html
    command = "conda run -n quast2 metaquast --unique-mapping --reuse-combined-alignments  --threads " + args.nbCores + " --no-check " + contigFilename + " -o " + outputDir + " -r " + referenceStr
    print(command)
    os.system(command)

if __name__ == "__main__":
    main(sys.argv[1:])  
