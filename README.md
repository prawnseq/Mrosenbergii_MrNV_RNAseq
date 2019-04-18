# Transcriptomic analysis (RNAseq) of Macrobrachium rosenbergii post-larvae in response to Macrobrachium rosenbergii nodavirus (MrNV) infection.

## Publication information

## Data availibility

## About
This pipeline was used for the transcriptome assembly and differential expression analysis of M.rosenbergii in response to MrNV infection. The pipeline was written using **[Snakemake](https://snakemake.readthedocs.io/en/stable/)** tool. The pipeline contained three major sections including raw data pre-processing, transcriptome assembly, and post-processing of the transcriptome (see below). 

![alt text](https://github.com/prawnseq/Mrosenbergii_MrNV_RNAseq/blob/master/AnalysisPipeline.png "analysis pipeline")

## Processing steps
Raw data pre-processing:
1. Quality assessment (**[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)** v 0.11.5 and **[MultiQC](https://multiqc.info)** v 1.8)
2. Quality trimming (**[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)** v 0.36)
3. Merge

Transcriptome assembly:
1. De novo transcriptome assembly (**[Trinity](https://github.com/trinityrnaseq/trinityrnaseq)** v 2.8.0)
2. Clustering redundant transcript (**[CD-HIT](http://weizhongli-lab.org/cd-hit/)**)

Post-processing:
1. Transcriptome quality assessment
