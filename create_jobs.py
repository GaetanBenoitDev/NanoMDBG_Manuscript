

"""
python3 ~/zeus/scripts/nanopore/evaluation/create_jobs.py ~/appa/run/paper/jobs/ ~/appa/run/paper/experiments 32

"""

import os, sys, argparse, shutil, glob, gzip

exeDirName = os.path.dirname(os.path.realpath(__file__))

def main(argv):
		
	parser = argparse.ArgumentParser()

	parser.add_argument("outputDirJobs", help="")
	parser.add_argument("outputDirResults", help="")
	parser.add_argument("nbCores", help="")

	args = parser.parse_args()

	outputDirJobs = args.outputDirJobs
	outputDirResults = args.outputDirResults
	nbCores = args.nbCores

		
	nanoporeDatasets = []
	nanoporeDatasets.append(("Zymo_R10", "Nanopore", "50G", "~/appa/data/paper/nanopore/zymo/Zymo_mock.fastq.gz"))
	nanoporeDatasets.append(("Zymo", "Nanopore", "50G", "~/appa/data/nanopore/mock/Zymo_D6331/SRR17913199_Q20/SRR17913199.1.fastq.gz"))
	nanoporeDatasets.append(("Zymo_R9", "Nanopore", "50G", "~/appa/data/nanopore/mock/Zymo_D6331/SRR17913200.1.fastq.gz"))
	nanoporeDatasets.append(("AD_R9", "Nanopore", "50G", "~/appa/data/nanopore/R9/AD/ERR7014876.fastq.gz"))
	nanoporeDatasets.append(("HumanGut", "Nanopore", "10G", "~/appa/data/paper/nanopore/humanGut/HuFilterd_GA_Input_pass_subsample_10G/HuFilterd_GA_Input_pass_10G.fastq.gz"))
	nanoporeDatasets.append(("HumanGut", "Nanopore", "20G", "~/appa/data/paper/nanopore/humanGut/HuFilterd_GA_Input_pass_subsample_10G/HuFilterd_GA_Input_pass_20G.fastq.gz"))
	nanoporeDatasets.append(("HumanGut", "Nanopore", "30G", "~/appa/data/paper/nanopore/humanGut/HuFilterd_GA_Input_pass_subsample_10G/HuFilterd_GA_Input_pass_30G.fastq.gz"))
	nanoporeDatasets.append(("HumanGut", "Nanopore", "40G", "~/appa/data/paper/nanopore/humanGut/HuFilterd_GA_Input_pass_subsample_10G/HuFilterd_GA_Input_pass_40G.fastq.gz"))
	nanoporeDatasets.append(("HumanGut", "Nanopore", "50G", "~/appa/data/paper/nanopore/humanGut/HuFilterd_GA_Input_pass_subsample/HuFilterd_GA_Input_pass_50G.fastq.gz"))
	nanoporeDatasets.append(("HumanGut", "Nanopore", "100G", "~/appa/data/paper/nanopore/humanGut/HuFilterd_GA_Input_pass_subsample/HuFilterd_GA_Input_pass_100G.fastq.gz"))
	nanoporeDatasets.append(("HumanGut", "Nanopore", "150G", "~/appa/data/paper/nanopore/humanGut/HuFilterd_GA_Input_pass_subsample/HuFilterd_GA_Input_pass_150G.fastq.gz"))
	nanoporeDatasets.append(("HumanGut", "Nanopore", "200G", "~/appa/data/paper/nanopore/humanGut/HuFilterd_GA_Input_pass_subsample/HuFilterd_GA_Input_pass_200G.fastq.gz"))
	nanoporeDatasets.append(("ZymoFecalReference", "Nanopore", "50G", "~/appa/data/nanopore/zymo_fecal_reference/zymoFecalReference_subsample/zymoFecalReference_50G.fastq.gz"))
	nanoporeDatasets.append(("ZymoFecalReference", "Nanopore", "100G", "~/appa/data/nanopore/zymo_fecal_reference/zymoFecalReference_subsample/zymoFecalReference_100G.fastq.gz"))
	nanoporeDatasets.append(("ZymoFecalReference", "Nanopore", "150G", "~/appa/data/nanopore/zymo_fecal_reference/zymoFecalReference_subsample/zymoFecalReference_150G.fastq.gz"))
	nanoporeDatasets.append(("ZymoFecalReference", "Nanopore", "200G", "~/appa/data/nanopore/zymo_fecal_reference/zymoFecalReference_subsample/zymoFecalReference_200G.fastq.gz"))
	nanoporeDatasets.append(("Soil", "Nanopore", "50G", "~/appa/data/paper/nanopore/soil/CEHSoil_ONT_subsample/CEHSoil_ONT_50G.fastq.gz"))
	nanoporeDatasets.append(("Soil", "Nanopore", "100G", "~/appa/data/paper/nanopore/soil/CEHSoil_ONT_subsample/CEHSoil_ONT_100G.fastq.gz"))
	nanoporeDatasets.append(("Soil", "Nanopore", "150G", "~/appa/data/paper/nanopore/soil/CEHSoil_ONT_subsample/CEHSoil_ONT_150G.fastq.gz"))
	nanoporeDatasets.append(("Soil", "Nanopore", "200G", "~/appa/data/paper/nanopore/soil/CEHSoil_ONT_subsample/CEHSoil_ONT_200G.fastq.gz"))
	nanoporeDatasets.append(("Soil", "Nanopore", "250G", "~/appa/data/paper/nanopore/soil/CEHSoil_ONT_subsample/CEHSoil_ONT_250G.fastq.gz"))
	nanoporeDatasets.append(("Soil", "Nanopore", "300G", "~/appa/data/paper/nanopore/soil/CEHSoil_ONT_subsample/CEHSoil_ONT_300G.fastq.gz"))
	nanoporeDatasets.append(("Soil", "Nanopore", "400G", "~/appa/data/paper/nanopore/soil/CEHSoil_ONT_subsample/CEHSoil_ONT_400G.fastq.gz"))
	nanoporeDatasets.append(("Soil", "Nanopore", "500G", "~/appa/data/paper/nanopore/soil/full/CEHSoil_ONT_filtered_subsample/CEHSoil_ONT_filtered_500G.fastq.gz"))
	nanoporeDatasets.append(("Soil", "Nanopore", "600G", "~/appa/data/paper/nanopore/soil/full/CEHSoil_ONT_filtered_subsample/CEHSoil_ONT_filtered_600G.fastq.gz"))


	hifiDatasets = []
	hifiDatasets.append(("Zymo_Hifi", "Hifi", "50G", "~/appa/data/nanopore/mock/Zymo_D6331/hifi/SRR13128014.std.fastq.gz"))
	hifiDatasets.append(("HumanGut", "Hifi", "10G", "~/appa/data/paper/hifi/humanGut/hu_HiFi_filt_trim_subsample_10G/hu_HiFi_filt_trim_10G.fastq.gz"))
	hifiDatasets.append(("HumanGut", "Hifi", "20G", "~/appa/data/paper/hifi/humanGut/hu_HiFi_filt_trim_subsample_10G/hu_HiFi_filt_trim_20G.fastq.gz"))
	hifiDatasets.append(("HumanGut", "Hifi", "30G", "~/appa/data/paper/hifi/humanGut/hu_HiFi_filt_trim_subsample_10G/hu_HiFi_filt_trim_30G.fastq.gz"))
	hifiDatasets.append(("HumanGut", "Hifi", "40G", "~/appa/data/paper/hifi/humanGut/hu_HiFi_filt_trim_subsample_10G/hu_HiFi_filt_trim_40G.fastq.gz"))
	hifiDatasets.append(("HumanGut", "Hifi", "50G", "~/appa/data/paper/hifi/humanGut/hu_HiFi_filt_trim_subsample/hu_HiFi_filt_trim_50G.fastq.gz"))
	hifiDatasets.append(("ZymoFecalReference", "Hifi", "50G", "~/appa/data/hifi/zymo_human/zymoFecalReferenceHifi_subsample/zymoFecalReferenceHifi_50G.fastq.gz"))
	hifiDatasets.append(("ZymoFecalReference", "Hifi", "100G", "~/appa/data/hifi/zymo_human/zymoFecalReferenceHifi_subsample/zymoFecalReferenceHifi_100G.fastq.gz"))
	hifiDatasets.append(("ZymoFecalReference", "Hifi", "150G", "~/appa/data/hifi/zymo_human/zymoFecalReferenceHifi_subsample/zymoFecalReferenceHifi_150G.fastq.gz"))
	hifiDatasets.append(("ZymoFecalReference", "Hifi", "200G", "~/appa/data/hifi/zymo_human/zymoFecalReferenceHifi_subsample/zymoFecalReferenceHifi_200G.fastq.gz"))
	hifiDatasets.append(("Soil", "Hifi", "50G", "~/appa/data/paper/hifi/soil/CEH_hifi_sample_trimmed_subsample/CEH_hifi_sample_trimmed_50G.fastq.gz"))
	hifiDatasets.append(("Soil", "Hifi", "100G", "~/appa/data/paper/hifi/soil/CEH_hifi_sample_trimmed_subsample/CEH_hifi_sample_trimmed_100G.fastq.gz"))
	hifiDatasets.append(("Soil", "Hifi", "150G", "~/appa/data/paper/hifi/soil/CEH_hifi_sample_trimmed_subsample/CEH_hifi_sample_trimmed_150G.fastq.gz"))
	hifiDatasets.append(("Soil", "Hifi", "200G", "~/appa/data/paper/hifi/soil/CEH_hifi_sample_trimmed_subsample/CEH_hifi_sample_trimmed_200G.fastq.gz"))
	hifiDatasets.append(("Soil", "Hifi", "250G", "~/appa/data/paper/hifi/soil/CEH_hifi_sample_trimmed_subsample/CEH_hifi_sample_trimmed_250G.fastq.gz"))

	dechatDatasets = []
	dechatDatasets.append(("HumanGut", "Nanopore", "50G", "~/appa/run/paper/experiments/dechat/HumanGut/Nanopore/50G/HuFilterd_GA_Input_pass_50G_.ec.fa"))
	
	herroDatasets = []
	herroDatasets.append(("HumanGut", "Nanopore", "50G", "~/appa/run/paper/experiments/herro/HumanGut/Nanopore/50G/HuFilterd_GA_Input_pass_50G_doradoCorrected.fasta.gz"))

	for assemblerName in ["nanoMDBG", "metaMDBG", "metaflye", "hifiasm"]:
		for datasetInfo in nanoporeDatasets:
			createJob(outputDirJobs, outputDirResults, assemblerName, datasetInfo, nbCores)

	for assemblerName in ["metaMDBG", "metaflye", "hifiasm"]:
		for datasetInfo in hifiDatasets:
			createJob(outputDirJobs, outputDirResults, assemblerName, datasetInfo, nbCores)

	for assemblerName in ["nanoMDBG-dechat", "metaMDBG-dechat", "metaflye-dechat", "hifiasm-dechat"]:
		for datasetInfo in dechatDatasets:
			createJob(outputDirJobs, outputDirResults, assemblerName, datasetInfo, nbCores)

	for assemblerName in ["nanoMDBG-herro", "metaMDBG-herro", "metaflye-herro", "hifiasm-herro"]:
		for datasetInfo in herroDatasets:
			createJob(outputDirJobs, outputDirResults, assemblerName, datasetInfo, nbCores)

	for datasetInfo in nanoporeDatasets:
		createJob_herro(outputDirJobs + "/herro/", outputDirResults + "/herro/", datasetInfo, nbCores)

	for datasetInfo in nanoporeDatasets:
		createJob_dechat(outputDirJobs + "/dechat/", outputDirResults + "/dechat/", datasetInfo, nbCores)

def createJob(outputDir, outputDirResults, assemblerName, datasetInfo, nbCores):

	jobDir = outputDir + "/" + datasetInfo[0] + "/" + datasetInfo[1] + "/" + datasetInfo[2] + "/" + assemblerName + "/"
	if not os.path.exists(jobDir): os.makedirs(jobDir)

	readFilename = datasetInfo[3]

	jobName = "job_" + assemblerName + "_" + datasetInfo[0] + "_" + datasetInfo[1] + "_" + datasetInfo[2] + ".sh"

	print(jobDir, jobName)
	jobFile = open(jobDir + "/" + jobName, "w")

	command = "python3 " + exeDirName + "/evaluation_pipeline.py " + datasetInfo[0] + " " + datasetInfo[1] + " " + datasetInfo[2] + " " + assemblerName + " " + readFilename + " " + outputDirResults + " " + nbCores

	jobFile.write("#!/bin/bash\n")
	jobFile.write("\n")
	jobFile.write(command + "\n")
	jobFile.write("\n")

	jobFile.close()

def createJob_herro(outputDir, outputDirResults, datasetInfo, nbCores):

	jobDir = outputDir + "/" + datasetInfo[0] + "/" + datasetInfo[1] + "/" + datasetInfo[2] + "/"
	if not os.path.exists(jobDir): os.makedirs(jobDir)

	readFilename = datasetInfo[3]

	jobName = "job_herro_" + datasetInfo[0] + "_" + datasetInfo[1] + "_" + datasetInfo[2] + ".sh"

	print(jobDir, jobName)
	jobFile = open(jobDir + "/" + jobName, "w")

	outputDir = outputDirResults + "/" + datasetInfo[0] + "/" + datasetInfo[1] + "/" + datasetInfo[2] + "/" 
	command = "python3 " + exeDirName + "/correction/run_herro.py " + readFilename + " " + outputDir + " " + nbCores

	jobFile.write("#!/bin/bash\n")
	jobFile.write("\n")
	jobFile.write(command + "\n")
	jobFile.write("\n")

	jobFile.close()

def createJob_dechat(outputDir, outputDirResults, datasetInfo, nbCores):

	jobDir = outputDir + "/" + datasetInfo[0] + "/" + datasetInfo[1] + "/" + datasetInfo[2] + "/"
	if not os.path.exists(jobDir): os.makedirs(jobDir)

	readFilename = datasetInfo[3]

	jobName = "job_dechat_" + datasetInfo[0] + "_" + datasetInfo[1] + "_" + datasetInfo[2] + ".sh"

	print(jobDir, jobName)
	jobFile = open(jobDir + "/" + jobName, "w")

	outputDir = outputDirResults + "/" + datasetInfo[0] + "/" + datasetInfo[1] + "/" + datasetInfo[2] + "/" 
	command = "python3 " + exeDirName + "/correction/run_dechat.py " + readFilename + " " + outputDir + " " + nbCores

	jobFile.write("#!/bin/bash\n")
	jobFile.write("\n")
	jobFile.write(command + "\n")
	jobFile.write("\n")

	jobFile.close()

if __name__ == "__main__":
	main(sys.argv[1:])  
