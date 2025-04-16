
import os, sys, argparse, shutil, glob, gzip

exeDirName = os.path.dirname(os.path.realpath(__file__))

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
        inputArg = "--in-hifi"
    elif dataType == "Nanopore":
        inputArg = "--in-ont"

    command = "hifiasm_meta -o " + outputDir + "/asm" + " -t " + args.nbCores + " " + args.readFilename
    runCommand(command, "hifiasm-meta", outputDir)

    if os.path.exists(outputDir + "/asm.ec.bin"): os.remove(outputDir + "/asm.ec.bin")
    if os.path.exists(outputDir + "/asm.ec.mt.bin"): os.remove(outputDir + "/asm.ec.mt.bin")
    if os.path.exists(outputDir + "/asm.ovecinfo.bin"): os.remove(outputDir + "/asm.ovecinfo.bin")
    if os.path.exists(outputDir + "/asm.ovlp.reverse.bin"): os.remove(outputDir + "/asm.ovlp.reverse.bin")
    if os.path.exists(outputDir + "/asm.ovlp.source.bin"): os.remove(outputDir + "/asm.ovlp.source.bin")

    command = "pigz " + outputDir + "/asm.p_utg.gfa"
    os.system(command)

    command = "pigz " + outputDir + "/asm.r_utg.gfa"
    os.system(command)

    command = "python3 " + exeDirName + "/gfaToFasta.py " + outputDir + "/asm.p_ctg.gfa" + " " + outputDir + "/asm.p_ctg.fasta"
    os.system(command)

    command = "pigz " + outputDir + "/asm.p_ctg.gfa"
    os.system(command)

    command = "pigz " + outputDir + "/asm.p_ctg.fasta"
    os.system(command)

    command = "pigz " + outputDir + "/asm.rescue.all.fa"
    os.system(command)

    command = "pigz " + outputDir + "/asm.rescue.fa"
    os.system(command)

    command = "pigz " + outputDir + "/asm.a_ctg.gfa"
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
