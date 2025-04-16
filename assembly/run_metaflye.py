



import os, sys, argparse, shutil, glob, gzip

def main(argv):
    
    parser = argparse.ArgumentParser()

    parser.add_argument("readFilename", help="")
    parser.add_argument("dataType", help="Hifi | Nanopore")
    parser.add_argument("outputDir", help="")
    parser.add_argument("nbCores", help="")
    
    args = parser.parse_args()

    outputDir = args.outputDir
    dataType = args.dataType

    if not os.path.exists(outputDir): os.makedirs(outputDir)

    inputArg = ""
    if dataType == "Hifi":
        inputArg = "--pacbio-hifi"
    elif dataType == "Nanopore":
        inputArg = "--nano-hq"

    command = "flye " + inputArg + " " + args.readFilename + " --out-dir " + outputDir + " --threads " + args.nbCores + " --meta"
    runCommand(command, "flye", outputDir)

    shutil.rmtree(outputDir + "/00-assembly")
    shutil.rmtree(outputDir + "/10-consensus")
    shutil.rmtree(outputDir + "/20-repeat")
    shutil.rmtree(outputDir + "/30-contigger")
    shutil.rmtree(outputDir + "/40-polishing")

    command = "pigz " + outputDir + "/assembly.fasta"
    os.system(command)

    command = "pigz " + outputDir + "/assembly_graph.gfa"
    os.system(command)

def runCommand(command, condaEnvName, experimentOutputDir):

    print(command)
    
    command = "{ /usr/bin/time -v " + command + "; } 2> " + experimentOutputDir + "/perf.txt"
    
    jobFilename = experimentOutputDir + "/job.sh"
    jobFile = open(jobFilename, "w")
    jobFile.write(command)
    jobFile.close()

    command = "conda run -n " + condaEnvName + " sh " + jobFilename

    ret = os.system(command)
    if ret != 0:
        print("Failed")
        sys.exit(1)




if __name__ == "__main__":
    main(sys.argv[1:])  
