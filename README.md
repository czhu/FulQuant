# FulQuant
FulQuant method to identify and quantify transcript isoforms based on Nanopore long-read data. It integrates accurate full-length read identification, transcript quantification and visualization from ONT cDNA sequences. This new method allows *de novo* transcript annotation and employs stringent and complex criteria for filtering reads, alignments, and transcripts, to rule out artifacts from reverse transcription and sequencing errors, and groups reads into transcripts based on their splice sites, generating a set of highly confident transcripts with abundance estimates.


# Software dependency
To run all the scripts you will need:

1. R with dependent packages
2. samtools
3. bedtools
4. minimap2
5. java
6. python with dependent packages
7. custom R [nanopore](https://github.com/czhu/R_nanopore/) package

# Genome dependency
in genome folder
1. indexed genome for minimap2 (setup in sw/run_minimap2.sh)
2. genome file for chromosome lengths (example file included for hg38)
3. polyA file (example included for hg38)
4. reference transcriptome for classification (example included for hg38)

# How to run the pipeline
Run the scripts in the scripts folder step by step. You need to adapt project relevant trivial parameters such as `infolder`, `outfolder`, `ncpu` in each folder.

You fastq data should have both adapters (untrimmed) because the pipeline needs this information to identify full-length reads. Also please make sure you give the correct sequencing adapters in `01_identify_ispcr.R`.

Depending on the depth of the data, please change the filtering parameters in filter_clustering.R script. Currently, we only consider transcripts with at least 3 read count.

We are in the process of wrapping all scripts into R package. Some important functions are included in the [custom nanopore package](https://github.com/czhu/R_nanopore/).

Please do not hesitate to contact me in case of questions: (czhu5 [at] stanford.edu).