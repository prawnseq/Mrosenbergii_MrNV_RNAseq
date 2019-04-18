# Transcriptomic analysis (RNAseq) of Macrobrachium rosenbergii post-larvae in response to Macrobrachium rosenbergii nodavirus (MrNV) infection.

## Publication information

## Data availibility

## About
This pipeline was used for the transcriptome assembly and differential expression analysis of M.rosenbergii in response to MrNV infection.  **[Snakemake](https://snakemake.readthedocs.io/en/stable/)** tool was used to create reproducible and scalable automated pipeline. The pipeline contained three major sections including raw data pre-processing, transcriptome assembly, and post-processing of the transcriptome (see below). 

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
1. Transcriptome quality assessment:
- Fragment mapping rates  (**[Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)**)
- Orthologs completeness against arthropoda_odb9 database (**[BUSCO](https://busco.ezlab.org)**)
- **[ExN50](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Transcriptome-Contig-Nx-and-ExN50-stats)** statistics
2. Transctiprome annotation:
- Homology search against **[UniProt](https://data.broadinstitute.org/Trinity/Trinotate_v3_RESOURCES/uniprot_sprot.pep.gz)** and **[non-redundant arthropods](https://ftp.ncbi.nlm.nih.gov/blast/db/v5/nr_v5.*.tar.gz)** database (**[Blastx](https://www.ncbi.nlm.nih.gov/BLAST/)** v 2.8.0)
- Obtain functional annotation from 
