



import os, sys, argparse, shutil, glob

def main(argv):
    
    parser = argparse.ArgumentParser()

    parser.add_argument("binDir", help="")
    parser.add_argument("outputDir", help="")
    parser.add_argument("nbCores", help="")
    
    args = parser.parse_args()

    binDir = args.binDir
    outputDir = args.outputDir + "/__checkm/"
    nbCores = args.nbCores


    resultFilename = outputDir + "/quality_report.tsv"
    outputBinDir = outputDir + "/bins_/"
    completeDir = outputBinDir + "/complete"
    highDir = outputBinDir + "/high"
    medDir = outputBinDir + "/med"

    if not os.path.exists(outputDir):
        os.makedirs(outputDir)

    nbBins = 0
    for binFilename in glob.glob(binDir + "/*.fa"):
        nbBins += 1
    
    if nbBins > 0 and not os.path.exists(resultFilename):
        
        #shutil.rmtree(outputDir)

        outputFilename = outputDir + "/checkm_results.txt"

        command = "conda run -n checkm2 checkm2 predict --force --threads " + str(nbCores) + " -x fa -i " + binDir + " -o " + outputDir
        #command = "conda run -n checkm2 checkm lineage_wf -r -x fa " + binDir + " " + outputDir + " -t " + str(nbCores) + " --pplacer_threads " + str(nbCores) + " >  " + outputFilename
        print(command)
        os.system(command)

        #command = "conda run -n checkm checkm2 qa -o 2 --tab_table -f " + resultFilename + " " + outputDir + "/lineage.ms" + " " + outputDir  
        #print(command)
        #os.system(command)

        if os.path.exists(outputDir + "/diamond_output/"): shutil.rmtree(outputDir + "/diamond_output/")
        if os.path.exists(outputDir + "/protein_files/"): shutil.rmtree(outputDir + "/protein_files/")
        
    #else:
    #    print("Output dir exists!")

    if os.path.exists(completeDir): shutil.rmtree(completeDir)
    os.makedirs(completeDir)

    if os.path.exists(highDir): shutil.rmtree(highDir)
    os.makedirs(highDir)

    if os.path.exists(medDir): shutil.rmtree(medDir)
    os.makedirs(medDir)


    nbMags_nearComplete = 0
    nbMags_highQuality = 0
    nbMags_mediumQuality = 0
    nbMags_contaminated = 0

    if os.path.exists(resultFilename):
        f = open(resultFilename)
        f.readline() #skip header

        for line in f:
            line = line.rstrip()
            if line[0] == "-": continue

            fields = line.split("\t")

            binName = fields[0]
            binFilename = binName + ".fa"
            #binFilename = outputDir + "/dereplicated_genomes/" + binName
            #print(binFilename)
            #score = float(fields[1])
            completeness = float(fields[1])
            contamination = float(fields[2])
            #strain = float(fields[7])
            #strainContamination = contamination * strain/100
            #contamination -= strainContamination

            #print(line)
            #print(fields)
            #print(completeness, contamination, strain)
            
            #score = completeness - 5*contamination
            
            #command = ""

            if contamination > 10:
                print("Contaminated: ", contamination)
                #command = "ln -s " + binFilename + " " + contaminatedDir + "/" + binName
                nbMags_contaminated += 1
            
            if completeness >= 90 and contamination <= 5:
                #command = "ln -s " + binFilename + " " + completeDir + "/" + binName
                if os.path.exists(binDir + "/" + binFilename): shutil.copy2(binDir + "/" + binFilename, completeDir + "/" + binFilename)
                nbMags_nearComplete += 1
            elif completeness >= 70 and contamination <= 10:
                #command = "ln -s " + binFilename + " " + highDir + "/" + binName
                if os.path.exists(binDir + "/" + binFilename):  shutil.copy2(binDir + "/" + binFilename, highDir + "/" + binFilename)
                nbMags_highQuality += 1
            elif completeness >= 50 and contamination <= 10:
                #command = "ln -s " + binFilename + " " + medDir + "/" + binName
                if os.path.exists(binDir + "/" + binFilename):  shutil.copy2(binDir + "/" + binFilename, medDir + "/" + binFilename)
                nbMags_mediumQuality += 1

            #if command != "": os.system(command)

        f.close()

    #for circularStr in qualityScores.keys():
    #    print(circularStr, ": ", qualityScores[circularStr]["high"], qualityScores[circularStr]["med"], qualityScores[circularStr]["low"], "    ", qualityScores[circularStr]["contaminated"])

    f = open(outputDir + "/binScore.csv", "w")
    f.write("Near-complete: " + str(nbMags_nearComplete) + "\n")
    f.write("High-quality: " + str(nbMags_highQuality) + "\n")
    f.write("Medium-quality: " + str(nbMags_mediumQuality) + "\n")
    f.write("Contaminated: " + str(nbMags_contaminated) + "\n")
    #for circularStr in qualityScores.keys():
    #    f.write(circularStr + ": " + str(qualityScores[circularStr]["high"]) + " " + str(qualityScores[circularStr]["med"]) + " " + str(qualityScores[circularStr]["low"])+ "    " + str(qualityScores[circularStr]["contaminated"]) + "\n")
    f.close()


if __name__ == "__main__":
    main(sys.argv[1:])  
