


#python3 ~/zeus/scripts/nanopore/evaluation/computeReferenceCompleteness.py ~/appa/data/nanopore/mock/Zymo_D6331/references_filenames.txt ~/appa/run/paper/experiments/Zymo_R10/Nanopore/50G/nanoMDBG/asm/contigs.fasta.gz ~/appa/run/paper/experiments/Zymo_R10/Nanopore/50G/nanoMDBG/asm/contigs.fasta.gz metaMDBG ~/appa/tmp/refComp_2/ 0.99 32

from mimetypes import init
import os, sys, argparse, shutil, glob, gzip
from Bio.SeqIO.FastaIO import SimpleFastaParser

#reprise: circular seulemen si > minCOmpleteness

class Align:

    def __init__(self, contigName, referenceName, similarity):
        self.referenceName = referenceName
        #self.referenceLength = referenceLength
        self.contigName = contigName
        self.totalAlignLength = 0
        self.subAligns = []
        self.similarity = similarity

    def addSubAlign(self, start, end, nbMatches):

        #if self.contigName == "ctg945c":
        #    print("Add sub align: ", self.contigName, end-start, start, end)
        if self.isContainedOverlap(start, end):
            #print("failed")
            return

        identity = nbMatches / (end - start)
        self.totalAlignLength += ((end - start) * identity)
        self.subAligns.append((start, end))

    def __repr__(self):

        s = self.contigName
        s += " -> " + self.referenceName
        s += " (" + str(self.totalAlignLength) + ")"

        return s

    def isContainedOverlap(self, start, end):

        for existingStart, existingEnd in self.subAligns:
            if start >= existingStart and end <= existingEnd: return True

        return False
        """
        hits = [False] * self.referenceLength

        for start, end in self.subAligns:
            for i in range(start, end):
                hits[i] = True

        nbContainedHits = 0
        alignLength = 0

        for i in range(0, len(hits)):
            if hits[i]:
                nbContainedHits += 1
            alignLength += 1

        return float(nbContainedHits) / float(alignLength) > 0.5
        """

class Reference:

    def __init__(self, referenceName, fastaFilename):
        self.fastaFilename = fastaFilename
        self.referenceName = referenceName
        self.fragmentHits = {}

    def addFragment(self, fragmentName, fragmentLength):
        self.fragmentHits[fragmentName] = [False] * fragmentLength

    def applyAlign(self, alignObject):

        print(alignObject)
        isNewAlign = False

        hits = self.fragmentHits[alignObject.referenceName]

        for start, end in alignObject.subAligns:
            if(self.isContainedOverlap(hits, start, end)): continue
            for i in range(start, end):
                hits[i] = True
            isNewAlign = True

        return isNewAlign

    def isContainedOverlap(self, hits, start, end):

        nbContainedHits = 0

        for i in range(start, end):
            if hits[i]:
                nbContainedHits += 1

        return float(nbContainedHits) / float(end-start) > 0.9

    def computeCompleteness(self):

        referenceLength = 0
        nbHits = 0

        for fragmentName, hits in self.fragmentHits.items():
            for i in range(0, len(hits)):
                if hits[i]: nbHits += 1
            referenceLength += len(hits)

        return float(nbHits) / float(referenceLength)

def main(argv):

    parser = argparse.ArgumentParser()

    #/mnt/gpfs/gaetan/run/experiments/rust-mdbg/AD_origin/binning/bin*.fa
    parser.add_argument("references", help="Reference input file (one fasta filename per line)")
    parser.add_argument("contigs", help="Contig filename (.fasta)")
    parser.add_argument("assemblyFilename", help="Assembly file that contaisn circular information")
    parser.add_argument("assembler", help="mdbg, hifiasm, metaflye")
    #parser.add_argument("minContigLength", help="min contig length for circular contigs")
    parser.add_argument("tmpDir", help="")
    parser.add_argument("minCompleteness", help="Min completeness fraction for the reference to be considered complete (0-1)")
    parser.add_argument("nbCores", help="")
    
    args = parser.parse_args()

    tmpDir = args.tmpDir
    if not os.path.exists(tmpDir):
        os.makedirs(tmpDir)

    resultFile = open(tmpDir + "/results.txt", "w")
    bandageFile = open(tmpDir + "/contigAssignments.csv", "w")
    bandageFile.write("Name,Color\n")

    nbCores = int(args.nbCores)
    minCompleteness = float(args.minCompleteness)
    #combinedReferenceFilename = args.references + "_combined_.fasta"

    assembler = args.assembler

    isLongContig, circularContigs = loadCircularContigs(args.assemblyFilename, assembler, 1000000, True)
        

    #circularContigs = loadCircularContigs(args.assembler, args.assemblyFilename, 1000000)
    references, validReferenceNames, fragmentName_to_reference, alignFilenames = getReferenceGenomes(args.references, tmpDir, args.contigs, nbCores)
    #print(references)




    #exit(1)
    alignments = []
    contigLengths = {}

    for alignFilename in alignFilenames:
        for line in open(alignFilename):
            line = line.rstrip()

            fields = line.split("\t")

            #if fields[5] != refName: continue
            contigName = fields[0]
            referenceName = fields[5]
            contigLength = int(fields[1])
            alignStart = int(fields[7])
            alignEnd = int(fields[8])
            nbMatches = int(fields[9])
            alignLength = int(fields[10])

            #print(referenceName in validReferenceNames, len(validReferenceNames))
            if not referenceName in validReferenceNames: continue

            if not contigName in contigLengths:
                contigLengths[contigName] = contigLength

            similarity = float(fields[12].replace("gi:f:", ""))

            alignments.append((contigName, referenceName, alignStart, alignEnd, alignEnd-alignStart, nbMatches, similarity))

            #if not contigName in aligns:
            #    aligns[contigName] = [] #Align(contigName, len(referenceHits))

            #for al in aligns[contigName]:
            #    if al.referenceName == referenceName:
            #        al.addSubAlign(alignStart, alignEnd)

            #aligns.append([alignLength, alignStart, alignEnd, contigName])

    alignments.sort(key=lambda x: x[4], reverse=True)

    #for i in range(0, 100):
    #    print(alignments[i])

    contigAlignments = {}

    for alignment in alignments:

        contigName, referenceName, alignStart, alignEnd, alignLength, nbMatches, similarity = alignment

        key = contigName + "-" + referenceName
        if not key in contigAlignments:
            contigAlignments[key] = Align(contigName, referenceName, similarity)

        contigAlignments[key].addSubAlign(alignStart, alignEnd, nbMatches)
        """
        alignObject = None

        for al in contigAlignments[contigName]:
            if al.referenceName == referenceName:
                alignObject = al
                break

        if alignObject is None:
            alignObject = Align(contigName, referenceName)
            contigAlignments[contigName].append(alignObject)
            
        alignObject.addSubAlign(alignStart, alignEnd)
        """

    
    
    #print("-----------------------")
    contigAlignments_sorted = []
    for key, alignObject in contigAlignments.items():
        contigAlignments_sorted.append(alignObject)
        #print(contigAlignment)
        #contigAlignment.sort(key=lambda x: x.totalAlignLength, reverse=True)
        #contigAlignments_bestHits.append(contigAlignment[0])

    
    contigAlignments_sorted.sort(key=lambda x: x.totalAlignLength, reverse=True)

    #for i in range(0, 100):
    #    print(contigAlignments_sorted[i])
    #sys.exit(1)
    #print("-----------------------")

    referenceContigs = {}
    completeReferences = {}
    mappedContigs = {}

    for alignObject in contigAlignments_sorted:
        #print(alignObject)
        contigName = alignObject.contigName
        similarity = alignObject.similarity

        referenceName = fragmentName_to_reference[alignObject.referenceName]
        if (contigName in mappedContigs) and (mappedContigs[contigName] != referenceName): continue
        #if referenceName in completeReferences:

        referenceObject = references[referenceName]
        isNewAlign = referenceObject.applyAlign(alignObject)

        if isNewAlign:
            completeness = referenceObject.computeCompleteness()
            print(referenceName, completeness)
            mappedContigs[contigName] = referenceName

            if not referenceName in referenceContigs:
                referenceContigs[referenceName] = [0, []]

            referenceContigs[referenceName][0] = completeness
            referenceContigs[referenceName][1].append(contigName) #(contigName, similarity)

            if not referenceName in completeReferences and completeness >= minCompleteness: #(completeness >= minCompleteness or contigName in circularContigs):
                completeReferences[referenceName] = [completeness, referenceContigs[referenceName][1].copy()]

            
        #print(alignObject.contigName + " -> " + fragmentName_to_reference[alignObject.referenceName])
        #print(contigAlignment)

    for referenceName, stats in referenceContigs.items():
        if not referenceName in completeReferences:
            completeReferences[referenceName] = stats


    sortedReferenceName = list(completeReferences.keys())
    sortedReferenceName.sort()

    for referenceName in sortedReferenceName:
        stats = completeReferences[referenceName]
        mappedContigs = stats[1]
        referenceObject = references[referenceName]
        
        #aniScore = 0
        #for contigName, similarity in mappedContigs:
        #    aniScore += similarity

        #if len(mappedContigs) != 0:
        #    aniScore /= len(mappedContigs)
        aniScore = computeANI(tmpDir, referenceObject.fastaFilename, mappedContigs, args.contigs)
        #stats[0] = round(stats[0]*100, 2)

        if len(mappedContigs) == 1 and mappedContigs[0] in circularContigs:
            print(referenceName, "(" + str(len(referenceObject.fragmentHits)), "contigs)" + ": ", "circular", stats[0], aniScore)
            resultFile.write(referenceName + " (" + str(len(referenceObject.fragmentHits)) + " contigs)" + ": " + " circular " + str(stats[0]) + " " + str(aniScore) + "\n")
        else:
            print(referenceName, "(" + str(len(referenceObject.fragmentHits)), "contigs)" + ": ", len(mappedContigs), "contigs", stats[0], aniScore)
            resultFile.write(referenceName + " (" + str(len(referenceObject.fragmentHits)) + " contigs)" + ": " + str(len(mappedContigs)) + " contigs " + str(stats[0]) + " " + str(aniScore) + "\n")

        for contigName in mappedContigs:
            resultFile.write("\t" + contigName + " (" + str(contigLengths[contigName])  + ")\n")
            bandageFile.write(contigName + "," + referenceName + "\n")
        resultFile.write("\n")
            
    for referenceName in references.keys():
        if not referenceName in completeReferences:
            print(referenceName + ": " + "lost")
            resultFile.write(referenceName + ": " + "lost\n")

    resultFile.close()
    bandageFile.close()

def getReferenceGenomes(inputFilename, tmpDir, contigFilename, nbCores):

    
    contigFilenameBgzip = tmpDir + "/__tmp_contigs_bgzip.fasta.gz"
    command = "zcat " + contigFilename + " | bgzip --threads " + str(nbCores) + " -c > " + contigFilenameBgzip
    os.system(command)

    command = "samtools faidx " + contigFilenameBgzip
    os.system(command)
    

    #combinedReferenceFile = open(combinedReferenceFilename, "w")
    references = {}
    fragmentName_to_reference = {}
    validReferenceNames = set()
    alignFilenames = []
    #referenceName = ""
    #referenceLength = 0

    for line in open(inputFilename):
        filename = line.rstrip()
        print(filename)
        referenceName = os.path.splitext(os.path.basename(filename))[0]
        print(referenceName)
        
        if referenceName == "Salmonella_enterica": continue

        alignFilename = tmpDir + "/" + referenceName + ".paf"
        alignFilenames.append(alignFilename)

        
        
        refFilenameBgzip = tmpDir + "/__tmp_ref.fasta.gz"
        command = "cat " + filename + " | bgzip --threads " + str(nbCores) + " -c > " + refFilenameBgzip
        os.system(command)

        command = "samtools faidx " + refFilenameBgzip
        os.system(command)

        #if not "Escherichia_coli" in referenceName: continue
        #if referenceName != "Escherichia_coli_B766": continue
        #if referenceName != "Escherichia_coli_B3008": continue
        #if referenceName != "Porphyromonas_gingivalis_ATCC_33277": continue


        #command = "minimap2 -t " + str(nbCores) + " -cxasm20 " + filename + " " + contigFilename + " > " + alignFilename
        #command = "wfmash -t " + str(nbCores) + " " + filename + " " + contigFilename + " > " + alignFilename

        command = "wfmash " + refFilenameBgzip + " " + contigFilenameBgzip + " -t " + str(nbCores) + " > " + alignFilename

        print(command)
        ret = os.system(command)
        if ret != 0: sys.exit(ret)
        

        if not referenceName in references:
            references[referenceName] = Reference(referenceName, filename)

        for header, seq in SimpleFastaParser(open(filename)):
            if len(seq) < 500000: continue #Remove plasmids
            fragmentName = header.split(" ")[0]
            print(fragmentName, len(seq))
            references[referenceName].addFragment(fragmentName, len(seq))

            #combinedReferenceFile.write(">" + header + "\n")
            #combinedReferenceFile.write(seq + "\n")

            validReferenceNames.add(fragmentName)
            fragmentName_to_reference[fragmentName] = referenceName
            #print ("Reference: ", header, len(seq))
            #if len(seq) > referenceLength:
            #    referenceLength = len(seq)
            #    referenceName = header.split(" ")[0]

    #return referenceName, referenceLength
    #combinedReferenceFile.close()

    return references, validReferenceNames, fragmentName_to_reference, alignFilenames






def loadCircularContigs(inputFilename, assembler, minContigLength, collectCircularContigs):

    isCircularContig = set()
    isLongContig = set()

    fileHandle = None

    if(".gz" in inputFilename):
        fileHandle = gzip.open(inputFilename, "rt")
    else:
        fileHandle = open(inputFilename)
        
    if "metaMDBG" in assembler or "nanoMDBG" in assembler:

        for header, seq in SimpleFastaParser(fileHandle):

            contigName, lengthStr, coverageStr, circularStr = header.split(" ")
            
            if circularStr.split("=")[1] == "yes":
                isCircularContig.add(contigName)

            if len(seq) >= minContigLength:
                isLongContig.add(contigName)

    elif "hifiasm" in assembler:
        
        for header, seq in SimpleFastaParser(fileHandle):

            contigName = header.split(" ")[0]

            if header.endswith("c"):
                isCircularContig.add(contigName)

            if len(seq) >= minContigLength:
                isLongContig.add(contigName)

    elif "metaflye" in assembler:
        fileHandle.readline() #skip header

        for line in fileHandle:
            line = line.rstrip()
            if len(line) == 0: continue
            fields = line.split("\t")

            contigName = fields[0]

            isCircular = fields[3] == "Y"
            if isCircular:
                isCircularContig.add(contigName)

            contigLength = int(fields[1])
            if contigLength >= minContigLength:
                isLongContig.add(contigName)

    if collectCircularContigs:
        return isLongContig, isCircularContig
    else:
        return isLongContig, None

def computeANI(tmpDir, referenceFilename, mappedContigs, contigFilename):

    aniInputDir = tmpDir + "/__aniInputDir__"
    aniResultDir = tmpDir + "/__aniResultDir__"

    if os.path.exists(aniInputDir):
        shutil.rmtree(aniInputDir)
    os.makedirs(aniInputDir)

    if os.path.exists(aniResultDir):
        shutil.rmtree(aniResultDir)
    #os.makedirs(aniResultDir)

    shutil.copy2(referenceFilename, aniInputDir)

    contigFile = open(aniInputDir + "/mappedContigs.fasta", "w")

    if(".gz" in contigFilename):
        fileHandle = gzip.open(contigFilename, "rt")
    else:
        fileHandle = open(contigFilename)

    for header, seq in SimpleFastaParser(fileHandle):
        
        contigName = header.split(" ")[0]

        if not contigName in mappedContigs: continue
        
        contigFile.write(">" + header + "\n")
        contigFile.write(seq + "\n")

    contigFile.close()
    
    command = "conda run -n pyani average_nucleotide_identity.py -i " + aniInputDir + " -o " + aniResultDir
    os.system(command)

    resultFile = open(aniResultDir + "/ANIm_percentage_identity.tab")
    aniScore = float(resultFile.readlines()[1].split("\t")[2])
    
    return aniScore

if __name__ == "__main__":
    main(sys.argv[1:])  
