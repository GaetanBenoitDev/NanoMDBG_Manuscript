
import os, sys, argparse, subprocess, shutil

exeDirName = os.path.dirname(os.path.realpath(__file__))

def main(argv):
    
    parser = argparse.ArgumentParser()

    parser.add_argument("outputDir", help="output dir")
    parser.add_argument("contigFilename", help="contigs fasta filename")
    parser.add_argument("contigInfoFilename", help="contigs info")
    parser.add_argument("assembler", help="mdbg | hifiasm | metaflye")
    parser.add_argument("nbCores", help="")
    
    args = parser.parse_args()

    outputDir = args.outputDir
    #checkm2ReportFilename = args.checkm2ReportFilename
    contigFilename = args.contigFilename
    contigInfoFilename = args.contigInfoFilename
    assembler = args.assembler
    minContigLength = 0
    nbCores = int(args.nbCores)
    #genomadOutputDir = outputDir + "/genomad/"

    outputFilename = outputDir + "/results.txt"
    outputDir += "/__results/"

    if os.path.exists(outputDir):
        #os.makedirs(outputDir)
        shutil.rmtree(outputDir)

    os.makedirs(outputDir)


    if not os.path.exists(outputDir): os.makedirs(outputDir)

    command = "conda run -n genomad genomad end-to-end --splits 8 --conservative --cleanup " + contigFilename + " " + outputDir + " ~/appa/data/genomad/genomad_db/ --threads " + str(nbCores)
    print(command)
    os.system(command)

    filename = os.path.split(contigFilename)[1]
    filename = filename.replace(".fasta", "")
    filename = filename.replace(".gz", "")

    virusResultFilename = outputDir + "/" + filename + "_summary/" + filename + "_virus_summary.tsv"
    plasmidResultFilename = outputDir + "/" + filename + "_summary/" + filename + "_plasmid_summary.tsv"

    command = "python3 " + exeDirName + "/countCircularContigs3.py " + contigInfoFilename + " " + assembler + " " + str(minContigLength) + " --virus " + virusResultFilename
    ret = runCommand(command)
    nbVirus = int(ret)
    print("Nb virus: ", nbVirus)

    command = "python3 " + exeDirName + "/countCircularContigs3.py " + contigInfoFilename + " " + assembler + " " + str(minContigLength) + " --virus " + virusResultFilename + " --circ"
    ret = runCommand(command)
    nbCircularVirus = int(ret)
    print("Nb circular virus: ", nbCircularVirus)

    command = "python3 " + exeDirName + "/countCircularContigs3.py " + contigInfoFilename + " " + assembler + " " + str(minContigLength) + " --plasmid " + plasmidResultFilename
    ret = runCommand(command)
    nbPlasmids = int(ret)
    print("Nb plasmids: ", nbPlasmids)

    command = "python3 " + exeDirName + "/countCircularContigs3.py " + contigInfoFilename + " " + assembler + " " + str(minContigLength) + " --plasmid " + plasmidResultFilename + " --circ"
    ret = runCommand(command)
    nbCircularPlasmids = int(ret)
    print("Nb circular plasmids: ", nbCircularPlasmids)

    outputFile = open(outputFilename, "w")
    
    outputFile.write("Virus: " + str(nbVirus) + "\n")
    outputFile.write("Circular virus: " + str(nbCircularVirus) + "\n")
    outputFile.write("Plasmids: " + str(nbPlasmids) + "\n")
    outputFile.write("Circular plasmids: " + str(nbCircularPlasmids) + "\n")
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
