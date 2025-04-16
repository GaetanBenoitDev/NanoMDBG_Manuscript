

import os, sys, argparse, shutil, subprocess


#os.chdir(exeDirName)

def main(argv):

    parser = argparse.ArgumentParser()

    parser.add_argument("outputDir", help="")
    parser.add_argument("contigFilename", help="contig file")
    parser.add_argument("circFile", help="assembly file containing circular information")
    parser.add_argument("assembler", help="mdbg | hifiasm | metaflye")
    parser.add_argument("nbCores", help="")
    
    args = parser.parse_args()

    outputDir = args.outputDir
    contigFilename = args.contigFilename
    circFile = args.circFile
    assembler = args.assembler
    nbCores = int(args.nbCores)
    minContigLength = 1000000

    outputFilename = outputDir + "/results.txt"
    outputDir += "/__results/"

    if os.path.exists(outputDir):
        #os.makedirs(outputDir)
        shutil.rmtree(outputDir)

    os.makedirs(outputDir)

    exeDirName = os.path.dirname(os.path.realpath(__file__))

    command = "python3 " + exeDirName + "/computeAssemblySize3.py " + contigFilename
    ret = runCommand(command)
    assemblySize = int(ret)
    print("Assembly size:", assemblySize)
    
    command = "seqtk comp " + contigFilename + " | cut -f 2 | sort -rn | awk '{ sum += $0; print $0, sum }' | tac | awk 'NR==1 { halftot=$2/2 } lastsize>halftot && $2<halftot { print $1 } { lastsize=$2 }'"
    #command = "python3 " + exeDirName + "/computeAssemblyN50.py " + contigFilename 
    ret = runCommand(command)
    n50 = int(ret)
    print("N50:", n50)

    command = "python3 " + exeDirName + "/countCircularContigs3.py " + circFile + " " + assembler + " " + str(minContigLength) 
    ret = runCommand(command)
    nbLongContigs = int(ret)
    print("Long contigs: ", nbLongContigs)
    
    command = "python3 " + exeDirName + "/countCircularContigs3.py " + circFile + " " + assembler + " " + str(minContigLength)  + " --circ"
    ret = runCommand(command)
    nbCircularContigs = int(ret)
    print("Long circular contigs: ", nbCircularContigs)

    command = "python3 " + exeDirName + "/computeMAG_singleContigs3.py " + outputDir + " " + contigFilename + " " + str(minContigLength) + " " + str(nbCores) 
    ret = runCommand(command)

    nbLongContigs_complete = 0
    nbCircularContigs_complete = 0
    checkm2_filename = outputDir + "/checkm/__checkm/quality_report.tsv" 

    #print(outputDir + "/checkm/__checkm/binScore.csv")

    if os.path.exists(checkm2_filename):
        command = "python3 " + exeDirName + "/countCircularContigs3.py " + circFile + " " + assembler + " " + str(minContigLength) + " --checkm2 " + checkm2_filename
        ret = runCommand(command)
        nbLongContigs_complete = int(ret)
        print("Long near-complete contigs: ", nbLongContigs_complete)
        
        command = "python3 " + exeDirName + "/countCircularContigs3.py " + circFile + " " + assembler + " " + str(minContigLength)  + " --circ " + " --checkm2 " + checkm2_filename
        ret = runCommand(command)
        nbCircularContigs_complete = int(ret)
        print("Long near-complete circular contigs: ", nbCircularContigs_complete)


    #if not os.path.exists(outputDir + "/checkm/__checkm/binScore.csv"):
    #    singleContigs_results = "0" #no contigs > minLength
    #else:

    #    f = open(outputDir + "/checkm/__checkm/binScore.csv")
    #    singleContigs_results = f.read()
    #    f.close()
    #    print("Checkm:\n", singleContigs_results + "\n")
    
    outputFile = open(outputFilename, "w")
    
    outputFile.write("Assembly size: " + str(assemblySize) + "\n")
    outputFile.write("N50: " + str(n50) + "\n")
    outputFile.write("Long contigs: " +  str(nbLongContigs) + "\n")
    outputFile.write("Circular contigs: " +  str(nbCircularContigs) + "\n")
    outputFile.write("Long near-complete contigs: " +  str(nbLongContigs_complete) + "\n")
    outputFile.write("Circular near-complete contigs: " +  str(nbCircularContigs_complete) + "\n")
    outputFile.close()


    print("\n")
    for line in open(outputFilename):
        line = line.rstrip()
        print(line)

def runCommand(command):
    print(command)
    sp = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    if sp.returncode is not None:
        print("Error")
        sys.exit(1)
    subprocess_return = sp.stdout.read()
    return subprocess_return

if __name__ == "__main__":
    main(sys.argv[1:])  
