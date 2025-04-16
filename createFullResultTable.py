
"""
python3 ~/zeus/scripts/nanopore/evaluation/createFullResultTable.py  ~/appa/run/paper/experiments/ ~/appa/run/paper/experiments_v2/ ~/zeus/scripts/nanopore/evaluation/results/
"""

import os, sys, argparse, shutil, subprocess


#os.chdir(exeDirName)

def main(argv):

	parser = argparse.ArgumentParser()

	parser.add_argument("experimentDir", help="")
	parser.add_argument("experimentDir_v2", help="")
	parser.add_argument("outputDir", help="")

	args = parser.parse_args()

	experimentDir = args.experimentDir
	experimentDir_v2 = args.experimentDir_v2
	outputDir = args.outputDir

	createFullResultTable(experimentDir, experimentDir_v2, outputDir)

def createFullResultTable(experimentDir, experimentDir_v2, outputDir):

	assemblerOrder = ["nanoMDBG", "metaMDBG", "metaflye", "hifiasm", "nanoMDBG-herro", "nanoMDBG-dechat", "metaMDBG-herro", "metaMDBG-dechat", "metaflye-herro", "metaflye-dechat"]
	sequencingTechnologyOrder = ["Nanopore", "Hifi"]
	datasetOrder = ["HumanGut", "ZymoFecalReference", "Soil"]
	#datasetOrder = ["Soil"]
	datasetSizeOrder = ["10G", "20G", "30G", "40G", "50G", "100G", "150G", "200G", "250G", "300G", "350G", "400G", "450G", "500G", "550G", "600G"]

	resultFileAll = open(outputDir + "/all.tsv", "w")
	resultFile = open(outputDir + "/mags.tsv", "w")
	correctionToolFile = open(outputDir + "/correctionTools.tsv", "w")
	perfFile = open(outputDir + "/perf.tsv", "w")
	fractionMappedReadsFile = open(outputDir + "/fractionMappedReads.tsv", "w")
	virusPlasmidsFile = open(outputDir + "/virusPlasmids.tsv", "w")

	header = ""
	header += "Dataset\t"
	header += "Sequencing Platform\t"
	header += "Dataset Size\t"
	header += "Assembler\t"
	header += "Assembly Size\t"
	header += "Assembly N50\t"
	header += ">1Mb contigs\t"
	header += ">1Mb near-complete contigs\t"
	header += ">1Mb circular contigs\t"
	header += ">1Mb near-complete circular contigs\t"
	header += "Near-complete MAGs\t"
	header += "High-quality MAGs\t"
	header += "Medium-quality MAGs\t"
	header += "Near-Complete tRNA rRNA MAGs"
	resultFile.write(header + "\n")

	header = ""
	header += "Dataset\t"
	header += "Sequencing Platform\t"
	header += "Dataset Size\t"
	header += "Assembler\t"
	header += "Assembly Size\t"
	header += "Assembly N50\t"
	header += ">1Mb contigs\t"
	header += ">1Mb near-complete contigs\t"
	header += ">1Mb circular contigs\t"
	header += ">1Mb near-complete circular contigs\t"
	header += "Near-complete MAGs\t"
	header += "High-quality MAGs\t"
	header += "Medium-quality MAGs\t"
	header += "Near-Complete tRNA rRNA MAGs\t"
	header += "Virus\t"
	header += "Circular Virus\t"
	header += "High-quality Circular Virus\t"
	header += "Plasmids\t"
	header += "Circular Plasmids"
	resultFileAll.write(header + "\n")

	header = ""
	header += "Dataset\t"
	header += "Sequencing Platform\t"
	header += "Dataset Size\t"
	header += "Assembler\t"
	header += "Correction\t"
	header += ">1Mb near-complete contigs\t"
	header += "Near-complete MAGs\t"
	header += "High-quality MAGs\t"
	header += "Medium-quality MAGs"
	correctionToolFile.write(header + "\n")

	header = ""
	header += "Dataset\t"
	header += "Sequencing Platform\t"
	header += "Dataset Size\t"
	header += "Assembler\t"
	header += "Runtime\t"
	header += "Peak Memory (GB)"
	perfFile.write(header + "\n")

	header = ""
	header += "Dataset\t"
	header += "Sequencing Platform\t"
	header += "Dataset Size\t"
	header += "Assembler\t"
	header += "Near-complete contigs\t"
	header += "Near-complete MAGs\t"
	header += "High-quality MAGs\t"
	header += "Medium-quality MAGs\t"
	fractionMappedReadsFile.write(header + "\n")

	header = ""
	header += "Dataset\t"
	header += "Sequencing Platform\t"
	header += "Dataset Size\t"
	header += "Assembler\t"
	header += "Virus\t"
	header += "Circular Virus\t"
	header += "High-quality Circular Virus\t"
	header += "Plasmids\t"
	header += "Circular Plasmids"
	virusPlasmidsFile.write(header + "\n")

	for datasetName in datasetOrder:
		for sequencingTechnology in sequencingTechnologyOrder:
			for datasetSize in datasetSizeOrder:
				for assemblerName in assemblerOrder:

					results = collectResults(experimentDir, experimentDir_v2, datasetName, sequencingTechnology, datasetSize, assemblerName, resultFile)
					if results is None: continue
					writeResults(datasetName, sequencingTechnology, datasetSize, assemblerName, results, resultFile)
					writeResultsAll(datasetName, sequencingTechnology, datasetSize, assemblerName, results, resultFileAll)
					writeCorrectionTool(datasetName, sequencingTechnology, datasetSize, assemblerName, results, correctionToolFile)
					writePerf(datasetName, sequencingTechnology, datasetSize, assemblerName, results, perfFile)
					writeFractionMappedReads(datasetName, sequencingTechnology, datasetSize, assemblerName, results, fractionMappedReadsFile)
					writeVirusPlasmids(datasetName, sequencingTechnology, datasetSize, assemblerName, results, virusPlasmidsFile)


	resultFileAll.close()
	resultFile.close()
	correctionToolFile.close()
	perfFile.close()
	fractionMappedReadsFile.close()
	virusPlasmidsFile.close()

def collectResults(experimentDir, experimentDir_v2, datasetName, sequencingTechnology, datasetSize, assemblerName, resultFile):

	results = {}

	if assemblerName == "nanoMDBG" or (assemblerName == "metaMDBG" and sequencingTechnology == "Hifi"):
		resultDir = experimentDir_v2 + "/" + datasetName + "/" + sequencingTechnology + "/" + datasetSize + "/" + assemblerName + "/"
	else:
		resultDir = experimentDir + "/" + datasetName + "/" + sequencingTechnology + "/" + datasetSize + "/" + assemblerName + "/"

	if not os.path.exists(resultDir):
		print("Experiment not started: ", resultDir)
		return None
	
	performancesResults = getPerformances(resultDir + "/asm/")
	if performancesResults is None: return None

	results["PeakMemory"] = performancesResults["PeakMemory"]
	results["CPU"] = performancesResults["CPU"]
	results["WallclockTime"] = performancesResults["WallclockTime"]
	results["CpuTime"] = performancesResults["CpuTime"]

	singleContigsResults = getResults(resultDir + "/singleContigs/", resultDir + "/singleContigs/results.txt")
	if singleContigsResults is None: return None

	results["Assembly size"] = singleContigsResults["Assemblysize"]
	results["Assembly N50"] = singleContigsResults["N50"]
	results["Long contigs"] = singleContigsResults["Longcontigs"]
	results["Long circular Contigs"] = singleContigsResults["Circularcontigs"]
	results["Long near-complete contigs"] = singleContigsResults["Longnear-completecontigs"]
	results["Long circular near-complete contigs"] = singleContigsResults["Circularnear-completecontigs"]

	binningResults = getResults(resultDir + "/binning/", resultDir + "/binning/checkm/__checkm/binScore.csv")
	if binningResults is None: return None

	results["Near-complete MAGs"] = binningResults["Near-complete"]
	results["High-quality MAGs"] = binningResults["High-quality"]
	results["Medium-quality MAGs"] = binningResults["Medium-quality"]
	results["Contaminated MAGs"] = binningResults["Contaminated"]

	virusPlasmidsResults = getResults(resultDir + "/virusPlasmids/", resultDir + "/virusPlasmids/results.txt")
	if virusPlasmidsResults is not None:

		results["Virus"] = virusPlasmidsResults["Virus"]
		results["Circular virus"] = virusPlasmidsResults["Circularvirus"]
		results["Plasmids"] = virusPlasmidsResults["Plasmids"]
		results["Circular plasmids"] = virusPlasmidsResults["Circularplasmids"]

	checkvResults = getResults(resultDir + "/checkv/", resultDir + "/checkv/results.txt")
	if checkvResults is not None:

		results["Virus High-quality"] = checkvResults["HighQuality"]

	barnapResults = getResults(resultDir + "/barnap/", resultDir + "/barnap/results.txt")
	if barnapResults is not None:

		results["Near-Complete tRNA rRNA MAGs"] = barnapResults["BarnapInfernalMags"]

	fractionMappedReadsResults = getResults(resultDir + "/fractionMappedReads/", resultDir + "/fractionMappedReads/results.txt")
	if fractionMappedReadsResults is not None:

		#print(resultDir)
		results["Coverage Near-complete contigs"] = str(round(float(fractionMappedReadsResults["single-complete"]), 2))
		results["Coverage Near-complete MAGs"] = str(round(float(fractionMappedReadsResults["complete"]), 2))
		results["Coverage High-quality MAGs"] = str(round(float(fractionMappedReadsResults["high"]), 2))
		results["Coverage Medium-quality MAGs"] = str(round(float(fractionMappedReadsResults["med"]), 2))

	return results

def getResults(resultDir, resultFilename):

	if not isExperimentFinished(resultDir):
		print("Experiment is not finished: ", resultDir)
		return None

	results = {}

	for line in open(resultFilename):

		line = line.rstrip()
		if len(line) == 0: continue
		line = line.replace(" ", "")
		fieldName, value = line.split(":")
		
		results[fieldName] = value

	return results

def getPerformances(assemblyDir):

	if not isExperimentFinished(assemblyDir):
		print("Experiment is not finished: ", assemblyDir)
		return None
	
	perfFilename = assemblyDir + "/perf.txt"
	performances = {}

	peakMemoryStr = "Maximumresidentsetsize(kbytes):"
	cpuStr = "PercentofCPUthisjobgot:"
	wallclockTimeStr = "Elapsed(wallclock)time(h:mm:ssorm:ss):"
	cpuTimeStr = "Usertime(seconds):"

	for line in open(perfFilename):
		line = line.rstrip()
		line = line.replace("\t", "")
		line = line.replace(" ", "")
		if len(line) == 0: continue

		if peakMemoryStr in line:
			performances["PeakMemory"] = str(float(line.replace(peakMemoryStr, ""))/1000000)
		elif cpuStr in line:
			performances["CPU"] = line.replace(cpuStr, "")
		elif wallclockTimeStr in line:
			performances["WallclockTime"] = line.replace(wallclockTimeStr, "")
		elif cpuTimeStr in line:
			performances["CpuTime"] = line.replace(cpuTimeStr, "")

	return performances

def isExperimentFinished(dir):
	return os.path.exists(dir + "/_isDone")

def writeResults(datasetName, sequencingTechnology, datasetSize, assemblerName, results, resultFile):

	if "dechat" in assemblerName: return
	if "herro" in assemblerName: return
	if datasetName == "HumanGut" and datasetSize != "50G": return
	if datasetName == "ZymoFecalReference" and datasetSize != "200G": return
	if datasetName == "Soil" and (datasetSize != "250G" and datasetSize != "400G"): return

	line = ""
	line += datasetName + "\t"
	line += sequencingTechnology + "\t"
	line += datasetSize.replace("G", "") + "\t"
	line += assemblerName + "\t"
	line += results["Assembly size"] + "\t"
	line += results["Assembly N50"] + "\t"
	line += results["Long contigs"] + "\t"
	line += results["Long near-complete contigs"] + "\t"
	line += results["Long circular Contigs"] + "\t"
	line += results["Long circular near-complete contigs"] + "\t"
	line += results["Near-complete MAGs"] + "\t"
	line += results["High-quality MAGs"] + "\t"
	line += results["Medium-quality MAGs"] + "\t"

	if "Near-Complete tRNA rRNA MAGs" in results:
		line += results["Near-Complete tRNA rRNA MAGs"] 
	else:
		line += "-"

	resultFile.write(line + "\n")
		
def writeResultsAll(datasetName, sequencingTechnology, datasetSize, assemblerName, results, resultFile):

	if "dechat" in assemblerName: return
	if "herro" in assemblerName: return

	line = ""
	line += datasetName + "\t"
	line += sequencingTechnology + "\t"
	line += datasetSize.replace("G", "") + "\t"
	line += assemblerName + "\t"
	line += results["Assembly size"] + "\t"
	line += results["Assembly N50"] + "\t"
	line += results["Long contigs"] + "\t"
	line += results["Long near-complete contigs"] + "\t"
	line += results["Long circular Contigs"] + "\t"
	line += results["Long circular near-complete contigs"] + "\t"
	line += results["Near-complete MAGs"] + "\t"
	line += results["High-quality MAGs"] + "\t"
	line += results["Medium-quality MAGs"] + "\t"

	if "Near-Complete tRNA rRNA MAGs" in results:
		line += results["Near-Complete tRNA rRNA MAGs"] + "\t"
	else:
		line += "-" + "\t"

	line += results["Virus"] + "\t"
	line += results["Circular virus"] + "\t"

	if "Virus High-quality" in results:
		line += results["Virus High-quality"] + "\t"
	else:
		line += "0\t"
		
	line += results["Plasmids"] + "\t"
	line += results["Circular plasmids"]

	resultFile.write(line + "\n")

def writeCorrectionTool(datasetName, sequencingTechnology, datasetSize, assemblerName, results, resultFile):

	if datasetName != "HumanGut": return
	if sequencingTechnology != "Nanopore": return
	if datasetSize != "50G": return


	if "dechat" in assemblerName:
		correctionMethod = "Dechat corrected"
	elif "herro" in assemblerName:
		correctionMethod = "Herro corrected"
	else:
		correctionMethod = "No correction"

	assemblerName = assemblerName.replace("-herro", "").replace("-dechat", "")

	line = ""
	line += datasetName + "\t"
	line += sequencingTechnology + "\t"
	line += datasetSize.replace("G", "") + "\t"
	line += assemblerName + "\t"
	line += correctionMethod + "\t"
	line += results["Long near-complete contigs"] + "\t"
	line += results["Near-complete MAGs"] + "\t"
	line += results["High-quality MAGs"] + "\t"
	line += results["Medium-quality MAGs"] 

	resultFile.write(line + "\n")

def writePerf(datasetName, sequencingTechnology, datasetSize, assemblerName, results, resultFile):

	if "dechat" in assemblerName: return
	if "herro" in assemblerName: return
	if datasetName == "HumanGut" and datasetSize != "50G": return
	if datasetName == "ZymoFecalReference" and datasetSize != "200G": return
	if datasetName == "Soil" and (datasetSize != "250G" and datasetSize != "400G"): return

	line = ""
	line += datasetName + "\t"
	line += sequencingTechnology + "\t"
	line += datasetSize.replace("G", "") + "\t"
	line += assemblerName + "\t"
	line += results["WallclockTime"].split(":")[0] + "\t"
	line += results["PeakMemory"] 

	resultFile.write(line + "\n")

def writeFractionMappedReads(datasetName, sequencingTechnology, datasetSize, assemblerName, results, resultFile):

	if "dechat" in assemblerName: return
	if "herro" in assemblerName: return
	if not "Coverage Near-complete contigs" in results: return
	if not "Nanopore" in sequencingTechnology: return
	if datasetName == "HumanGut" and datasetSize != "50G": return
	if datasetName == "ZymoFecalReference" and datasetSize != "200G": return
	if datasetName == "Soil" and datasetSize != "400G": return

	line = ""
	line += datasetName + "\t"
	line += sequencingTechnology + "\t"
	line += datasetSize.replace("G", "") + "\t"
	line += assemblerName + "\t"
	line += results["Coverage Near-complete contigs"] + "\t"
	line += results["Coverage Near-complete MAGs"] + "\t"
	line += results["Coverage High-quality MAGs"] + "\t"
	line += results["Coverage Medium-quality MAGs"] 

	resultFile.write(line + "\n")

def writeVirusPlasmids(datasetName, sequencingTechnology, datasetSize, assemblerName, results, resultFile):


	if "dechat" in assemblerName: return
	if "herro" in assemblerName: return
	if datasetName == "HumanGut" and datasetSize != "50G": return
	if datasetName == "ZymoFecalReference" and datasetSize != "200G": return
	if datasetName == "Soil" and (datasetSize != "250G" and datasetSize != "400G"): return
	
	"""
	if "dechat" in assemblerName: return
	if "herro" in assemblerName: return
	if not "Virus" in results: return
	if assemblerName == "hifiasm": return
	if datasetName == "HumanGut" and datasetSize != "50G": return
	if datasetName == "ZymoFecalReference" and datasetSize != "200G": return
	if "Nanopore" in sequencingTechnology:
		if datasetName == "Soil" and (datasetSize != "400G" and datasetSize != "250G"): return
	if "Hifi" in sequencingTechnology:
		if datasetName == "Soil" and datasetSize != "250G": return
	"""


	line = ""
	line += datasetName + "\t"
	line += sequencingTechnology + "\t"
	line += datasetSize.replace("G", "") + "\t"
	line += assemblerName + "\t"
	line += results["Virus"] + "\t"
	line += results["Circular virus"] + "\t"

	if "Virus High-quality" in results:
		line += results["Virus High-quality"] + "\t"
	else:
		line += "0\t"
		
	line += results["Plasmids"] + "\t"
	line += results["Circular plasmids"]

	
	resultFile.write(line + "\n")

if __name__ == "__main__":
	main(sys.argv[1:])  
