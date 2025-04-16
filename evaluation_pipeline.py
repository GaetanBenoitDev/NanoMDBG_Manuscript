


"""
- perf time a convertire en format readable

python3 ~/zeus/scripts/nanopore/evaluation/evaluation_pipeline.py metaflye ~/zeus/scripts/nanopore/evaluation/evaluation_pipeline_input_test.tsv ~/appa/tmp/test_eval 32
"""

import os, sys, argparse, shutil, glob, gzip
#from Bio import SeqIO

exeDirName = os.path.dirname(os.path.realpath(__file__))

def main(argv):


	parser = argparse.ArgumentParser()

	parser.add_argument("datasetName", help="")
	parser.add_argument("SequencingTechnology", help="")
	parser.add_argument("datasetSize", help="")
	parser.add_argument("assemblerName", help="")
	parser.add_argument("readFilename", help="")
	parser.add_argument("outputDir", help="")
	parser.add_argument("nbCores", help="")

	args = parser.parse_args()

	inputInfo = {}
	inputInfo["DatasetName"] = args.datasetName
	inputInfo["DatasetSize"] = args.datasetSize
	inputInfo["SequencingTechnology"] = args.SequencingTechnology
	inputInfo["AssemblerName"] = args.assemblerName
	inputInfo["ReadFilename"] = args.readFilename
	inputInfo["OutputDir"] = args.outputDir
	inputInfo["NbCores"] = args.nbCores


	if not os.path.exists(inputInfo["OutputDir"]): os.makedirs(inputInfo["OutputDir"])

	evaluateDataset(inputInfo)

def evaluateDataset(inputInfo):

	experimentOutputDir = inputInfo["OutputDir"] + "/" + inputInfo["DatasetName"] + "/" + inputInfo["SequencingTechnology"] + "/" + inputInfo["DatasetSize"] + "/" + inputInfo["AssemblerName"] + "/"
	if not os.path.exists(experimentOutputDir): os.makedirs(experimentOutputDir)

	print("Evaluating: ", inputInfo)
	print("Experiment dir: ", experimentOutputDir)

	runAssembly(inputInfo, experimentOutputDir)
	runSingleContigs(inputInfo, experimentOutputDir)
	runBinning(inputInfo, experimentOutputDir)
	run_fractionMappedReads(inputInfo, experimentOutputDir)


	if "Zymo_mock" in inputInfo["ReadFilename"]: #Mock evaluation #SRR17913199
		run_mockCompletenessEvaluation(inputInfo, experimentOutputDir)
		return


	run_virusPlasmids(inputInfo, experimentOutputDir)
	run_checkv(inputInfo, experimentOutputDir)
	run_barnap(inputInfo, experimentOutputDir)

def runAssembly(inputInfo, experimentOutputDir):

	assemblyDir = experimentOutputDir + "/asm/"
	if not os.path.exists(assemblyDir): os.makedirs(assemblyDir)

	if "metaflye" in inputInfo["AssemblerName"]:	
		command = "python3 " + exeDirName + "/assembly/run_metaflye.py " + inputInfo["ReadFilename"] + " " + inputInfo["SequencingTechnology"] + " " + assemblyDir + " " + str(inputInfo["NbCores"])
		runCommand(command, assemblyDir)
	elif "metaMDBG" in inputInfo["AssemblerName"]:
		command = "python3 " + exeDirName + "/assembly/run_metaMDBG.py " + inputInfo["ReadFilename"] + " " + inputInfo["SequencingTechnology"] + " " + assemblyDir + " " + str(inputInfo["NbCores"])
		runCommand(command, assemblyDir)
	elif "nanoMDBG" in inputInfo["AssemblerName"]:
		command = "python3 " + exeDirName + "/assembly/run_nanoMDBG.py " + inputInfo["ReadFilename"] + " " + inputInfo["SequencingTechnology"] + " " + assemblyDir + " " + str(inputInfo["NbCores"])
		runCommand(command, assemblyDir)
	elif "hifiasm" in inputInfo["AssemblerName"]:
		command = "python3 " + exeDirName + "/assembly/run_hifiasm.py " + inputInfo["ReadFilename"] + " " + inputInfo["SequencingTechnology"] + " " + assemblyDir + " " + str(inputInfo["NbCores"])
		runCommand(command, assemblyDir)



def runSingleContigs(inputInfo, experimentOutputDir):

	outputDir = experimentOutputDir + "/singleContigs/"
	if not os.path.exists(outputDir): os.makedirs(outputDir)

	contigFilename = getContigFilename(inputInfo["AssemblerName"], experimentOutputDir + "/asm/")
	assemblyInfoFilename = getContigInfoFilename(inputInfo["AssemblerName"], experimentOutputDir + "/asm/")

	command = "python3 " + exeDirName + "/run_singleContigs3.py " + outputDir + " " + contigFilename + " " + assemblyInfoFilename + " " + inputInfo["AssemblerName"] + " " + str(inputInfo["NbCores"])
	runCommand(command, outputDir)

def runBinning(inputInfo, experimentOutputDir):

	outputDir = experimentOutputDir + "/binning/"
	if not os.path.exists(outputDir): os.makedirs(outputDir)

	contigFilename = getContigFilename(inputInfo["AssemblerName"], experimentOutputDir + "/asm/")
	assemblyInfoFilename = getContigInfoFilename(inputInfo["AssemblerName"], experimentOutputDir + "/asm/")

	minimap2Preset = ""
	if inputInfo["SequencingTechnology"] == "Nanopore":
		minimap2Preset = "map-ont"
	elif inputInfo["SequencingTechnology"] == "Hifi":
		minimap2Preset = "map-hifi"

	command = "python3 " + exeDirName + "/computeMAG_binning3.py " + outputDir + " " + contigFilename + " " + minimap2Preset + " " + str(inputInfo["NbCores"]) + " " + inputInfo["ReadFilename"]
	runCommand(command, outputDir)



def run_virusPlasmids(inputInfo, experimentOutputDir):

	outputDir = experimentOutputDir + "/virusPlasmids/"
	if not os.path.exists(outputDir): os.makedirs(outputDir)

	contigFilename = getContigFilename(inputInfo["AssemblerName"], experimentOutputDir + "/asm/")
	assemblyInfoFilename = getContigInfoFilename(inputInfo["AssemblerName"], experimentOutputDir + "/asm/")

	command = "python3 ~/zeus/scripts/nanopore/evaluation/compute_virusPlasmids.py " + outputDir + " " + contigFilename + " " + assemblyInfoFilename + " " + inputInfo["AssemblerName"] + " " + str(inputInfo["NbCores"])
	runCommand(command, outputDir)
    
def run_checkv(inputInfo, experimentOutputDir):

	outputDir = experimentOutputDir + "/checkv/"
	if not os.path.exists(outputDir): os.makedirs(outputDir)

	contigFilename = getContigFilename(inputInfo["AssemblerName"], experimentOutputDir + "/asm/")
	assemblyInfoFilename = getContigInfoFilename(inputInfo["AssemblerName"], experimentOutputDir + "/asm/")

	filename = os.path.split(contigFilename)[1]
	filename = filename.replace(".fasta", "")
	filename = filename.replace(".gz", "")

	genomadSummaryFilename = experimentOutputDir + "/virusPlasmids/__results/" + filename + "_summary/" + filename + "_virus_summary.tsv"

	command = "python3 ~/zeus/scripts/nanopore/evaluation/run_checkv.py " + outputDir + " " + contigFilename + " " + assemblyInfoFilename + " " + inputInfo["AssemblerName"] + " " + genomadSummaryFilename + " " + str(inputInfo["NbCores"])
	runCommand(command, outputDir)

def run_barnap(inputInfo, experimentOutputDir):

	outputDir = experimentOutputDir + "/barnap/"
	if not os.path.exists(outputDir): os.makedirs(outputDir)

	binDir = experimentOutputDir + "/binning/checkm/__checkm/bins_/complete/"

	command = "python3 ~/zeus/scripts/nanopore/evaluation/run_barnap_infernal.py " + outputDir + " " + binDir + " " + str(inputInfo["NbCores"])
	runCommand(command, outputDir)

def run_fractionMappedReads(inputInfo, experimentOutputDir):

	outputDir = experimentOutputDir + "/fractionMappedReads/"
	if not os.path.exists(outputDir): os.makedirs(outputDir)

	contigFilename = getContigFilename(inputInfo["AssemblerName"], experimentOutputDir + "/asm/")
	assemblyInfoFilename = getContigInfoFilename(inputInfo["AssemblerName"], experimentOutputDir + "/asm/")

	minimap2Preset = ""
	if inputInfo["SequencingTechnology"] == "Nanopore":
		minimap2Preset = "map-ont"
	elif inputInfo["SequencingTechnology"] == "Hifi":
		minimap2Preset = "map-hifi"

	command = "python3 ~/zeus/scripts/nanopore/evaluation/computeFractionMappedReads.py " + outputDir + " " + contigFilename + " " + experimentOutputDir + " " + minimap2Preset + " " + str(inputInfo["NbCores"]) + " " + inputInfo["ReadFilename"]
	runCommand(command, outputDir)

def run_mockCompletenessEvaluation(inputInfo, experimentOutputDir):

	outputDir = experimentOutputDir + "/refComp/"
	if not os.path.exists(outputDir): os.makedirs(outputDir)

	contigFilename = getContigFilename(inputInfo["AssemblerName"], experimentOutputDir + "/asm/")
	assemblyInfoFilename = getContigInfoFilename(inputInfo["AssemblerName"], experimentOutputDir + "/asm/")

	referenceFilename = "~/appa/data/nanopore/mock/Zymo_D6331/references_filenames.txt"
	command = "python3 " + exeDirName + "/computeReferenceCompleteness.py " + referenceFilename + " " + contigFilename + " " + assemblyInfoFilename + " " + inputInfo["AssemblerName"] + " " + outputDir + " 0.99 " + str(inputInfo["NbCores"])
	runCommand(command, outputDir)


def getContigFilename(assemblerName, assemblyDir):
	if "metaflye" in assemblerName:
		return assemblyDir + "/assembly.fasta.gz"
	elif "metaMDBG" in assemblerName:
		return assemblyDir + "/contigs.fasta.gz"
	elif "nanoMDBG" in assemblerName:
		return assemblyDir + "/contigs.fasta.gz"
	elif "hifiasm" in assemblerName:
		return assemblyDir + "/asm.p_ctg.fasta.gz"
     
def getContigInfoFilename(assemblerName, assemblyDir):
	if "metaflye" in assemblerName:
		return assemblyDir + "/assembly_info.txt"
	elif "metaMDBG" in assemblerName:
		return assemblyDir + "/contigs.fasta.gz"
	elif "nanoMDBG" in assemblerName:
		return assemblyDir + "/contigs.fasta.gz"
	elif "hifiasm" in assemblerName:
		return assemblyDir + "/asm.p_ctg.fasta.gz"
     
def runCommand(command, experimentOutputDir):

	checkpointFilename = experimentOutputDir + "/_isDone"
	if os.path.exists(checkpointFilename): return

	print(command)

	command = "conda run -n magEvaluation " + command

	ret = os.system(command)
	if ret != 0:
		print("Failed")
		sys.exit(1)

	file = open(checkpointFilename, "w")
	file.close()




if __name__ == "__main__":
    main(sys.argv[1:])  
