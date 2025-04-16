

#python3 ~/zeus/scripts/nanopore/evaluation/computeReferenceCompleteness_mergeResults.py ~/appa/run/paper/experiments/Zymo_R10/Nanopore/50G/ ~/zeus/tmp/refCompResults.tsv

from mimetypes import init
import os, sys, argparse, shutil, glob, gzip

referenceName_to_abundance = {
"Akkermansia_muciniphila": "1.4%",
"Bacteroides_fragilis": "13.1%",
"Bifidobacterium_adolescentis": "1.3%",
"Candida_albican": "1.6%",
"Clostridioides_difficile": "1.8%",
"Clostridium_perfringens":  "0.0%",
"Enterococcus_faecalis": "0.0%",
"Escherichia_coli_B1109": "8.4%",
"Escherichia_coli_B3008": "8.2%",
"Escherichia_coli_B766": "7.8%",
"Escherichia_coli_JM109": "8.3%",
"Escherichia_coli_b2207": "8.3%",
"Faecalibacterium_prausnitzii": "14.4%",
"Fusobacterium_nucleatum": "3.8%",
"Lactobacillus_fermentum": "0.9%",
"Methanobrevibacter_smithii": "0.0%",
"Prevotella_corporis": "5.4%",
"Roseburia_hominis": "3.9%",
"Veillonella_rogosae": "11.0%",
#Candida albican	1.6%
#Saccharomyces cerevisiae	0.2%
#Salmonella enterica	0.0%
}






def main(argv):

    parser = argparse.ArgumentParser()

    parser.add_argument("inputDir", help="")
    parser.add_argument("outputFilename", help="")
    
    args = parser.parse_args()

    #tmpDir = args.tmpDir
    #if not os.path.exists(tmpDir):
    #    os.makedirs(tmpDir)


    results_nanoMDBG = loadResults(args.inputDir + "/nanoMDBG/refComp/results.txt")
    results_metaMDBG = loadResults(args.inputDir + "/metaMDBG/refComp/results.txt")
    results_metaflye = loadResults(args.inputDir + "/metaflye/refComp/results.txt")
    results_hifiasm = loadResults(args.inputDir + "/hifiasm/refComp/results.txt")

    #print(results_nanoMDBG)
    #print(results_metaMDBG)
    #print(results_metaflye)

    outputFile = open(args.outputFilename, "w")

    header1 = "\t\tNanoMDBG\t\t\t\tMetaMDBG\t\t\t\tMetaflye\t\t\t\tHifiasm-meta\t\t\n"
    header2 = "Species\tAbundance\tStatus\tCompleteness\tANI\tContig lenthgs\tStatus\tCompleteness\tANI\tContig lenthgs\tStatus\tCompleteness\tANI\tContig lenthgs\tStatus\tCompleteness\tANI\tContig lenthgs\n"
    outputFile.write(header1)
    outputFile.write(header2)

    for i in range(0, len(results_nanoMDBG)):
        referenceName = results_nanoMDBG[i]["referenceName"]

        print(referenceName)

        line = ""
        line += referenceName + "\t"
        line += referenceName_to_abundance[referenceName] + "\t"

        line += getStatusText(results_nanoMDBG[i]) + "\t"
        line += getCompletenessText(results_nanoMDBG[i]) + "\t"
        line += getAniText(results_nanoMDBG[i]) + "\t"
        line += getContigLengthsText(results_nanoMDBG[i]) + "\t"

        line += getStatusText(results_metaMDBG[i]) + "\t"
        line += getCompletenessText(results_metaMDBG[i]) + "\t"
        line += getAniText(results_metaMDBG[i]) + "\t"
        line += getContigLengthsText(results_metaMDBG[i]) + "\t"

        line += getStatusText(results_metaflye[i]) + "\t"
        line += getCompletenessText(results_metaflye[i]) + "\t"
        line += getAniText(results_metaflye[i]) + "\t"
        line += getContigLengthsText(results_metaflye[i]) + "\t"

        line += getStatusText(results_hifiasm[i]) + "\t"
        line += getCompletenessText(results_hifiasm[i]) + "\t"
        line += getAniText(results_hifiasm[i]) + "\t"
        line += getContigLengthsText(results_hifiasm[i]) + "\t"

        outputFile.write(line + "\n")

    outputFile.close()

def getStatusText(result):

    if result["status"] == "circular":
        return "circular"
    
    nbContigs = int(result["nbContigs"])

    if nbContigs == 1:
        return "1 contig"
    else:
        return str(nbContigs) + " contigs"

def getCompletenessText(result):

    completeness = float(result["completeness"])
    completeness *= 100

    return str(truncate(completeness, 2)) + "%" 

def getAniText(result):

    ani = float(result["ani"])
    ani *= 100

    return str(truncate(ani, 2)) + "%" 

def getContigLengthsText(result):
    
    contigLengths = result["contigLengths"]
    contigLengths.sort(reverse=True)

    print(contigLengths)

    size = min(10, len(contigLengths))
    s = ""

    for i in range(0, size):
        if i == size-1:
            s += str(contigLengths[i])
        else:
            s += str(contigLengths[i]) + ", "

    if len(contigLengths) > size:
        s += "..."
    return s


def truncate(f, n):
    s = '{}'.format(f)
    if 'e' in s or 'E' in s:
        return '{0:.{1}f}'.format(f, n)
    i, p, d = s.partition('.')
    return '.'.join([i, (d+'0'*n)[:n]])

def loadResults(resultFilename):

    results = []

    for line in open(resultFilename):
        line = line.rstrip()
        if len(line) == 0:
            results.append(result)
            #print(result)
            continue

        #print(line)

        if line[0] == "\t":
            fields = line.split(" ")
            contigLength = int(fields[1].replace("(", "").replace(")", ""))
            #print(fields)
            result["contigLengths"].append(contigLength)
        else:
            fields = line.split(" ")
            #print(fields)

            referenceName = fields[0]
            nbContigs = 0
            if len(fields[3]) != 0: nbContigs = int(fields[3])
            status = fields[4]
            completeness = float(fields[5])
            ani = float(fields[6])

            result = {
                "referenceName": referenceName, 
                "nbContigs": nbContigs, 
                "status": status, 
                "completeness": completeness,
                "ani": ani,
                "contigLengths": [],
            }


    return results

if __name__ == "__main__":
    main(sys.argv[1:])  
