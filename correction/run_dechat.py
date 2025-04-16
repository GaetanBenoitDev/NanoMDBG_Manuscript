



import os, sys, argparse, shutil, glob, gzip

def main(argv):
    
    parser = argparse.ArgumentParser()

    parser.add_argument("readFilename", help="")
    parser.add_argument("outputDir", help="")
    parser.add_argument("nbCores", help="")
    
    args = parser.parse_args()

    outputDir = args.outputDir

    if not os.path.exists(outputDir): os.makedirs(outputDir)

    readFilename = args.readFilename
    filename = os.path.split(readFilename)[1]
    filename = filename.replace(".fastq", "")
    filename = filename.replace(".gz", "")
    outputFilename = outputDir + "/" + filename + "_"

    command = "dechat -t " + args.nbCores + " -i " + args.readFilename + " -o " + outputFilename
    
    runCommand(command, "dechat", outputDir)


    #command = "pigz " + outputFilename
    #os.system(command)

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
