# NanoMDBG_Manuscript

A collection of scripts and workflows to assess HiFi metagenomics assembly quality.

# Dependencies
- biopython
- semiBin2
- checkm2
- minimap2 2.24+
- wfmash
- pyani

# Assess circular contigs
Assess contig quality in an assembly, and their quality with checkM.
```
MetaMDBG:
python3 ./run_singleContigs3.py outputDir contigs.fasta.gz contigs.fasta.gz metaMDBG nbCores

Metaflye:
python3 ./run_singleContigs3.py outputDir contigs.fasta.gz assembly_info.txt metaflye nbCores
```

After successful execution, the file outputDir/results.txt contains the results with the following format:
```
Assembly size: 5353314150
N50: 19080
Long contigs: 46
Circular contigs: 1
Long near-complete contigs: 0
Circular near-complete contigs: 0
```

# Assess non-circular MAGs (binning)
Reconstruct non-circular MAGs using semiBin2, and assess their quality with checkM.
```
python3 ./computeMAG_binning3.py outputDir contigs.fasta.gz map-ont nbCores reads_1.fastq.gz reads_2.fastq.gz...
```

After successful execution, the file "outputDir/checkm/\_\_checkm/binScore.csv" contains the results with the following format:
Near-complete: 5
High-quality: 20
Medium-quality: 47
Contaminated: 72

# Assess assembly completeness and fragmentation with reference sequences
Map contigs to references and compute the number of contigs required to cover at least 99% of the references.

```
python3 ./computeReferenceCompleteness.py referenceFile contigs.fasta.gz contigs.fasta.gz metaMDBG tmpDir 0.99 nbCores
```
"referenceFile" contains the list of reference filenames, one filename per line.
The result file tmpDir/results.txt provides the following information for each reference:
```
ReferenceName (number of contig of reference): [Assembly status: circular| number of contigs] completeness ANI_with_reference
Staphylococcus_aureus_ATCC_BAA_1556 (1 contigs):  circular 0.9991514199561156 0.999987844786959
```
