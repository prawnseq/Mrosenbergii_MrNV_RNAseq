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
- Fragment mapping rates  (**[Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)** v 2.3.0)
- Orthologs completeness against arthropoda_odb9 database (**[BUSCO](https://busco.ezlab.org)** v 3)
- **[ExN50](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Transcriptome-Contig-Nx-and-ExN50-stats)** statistics
2. Transctiprome annotation:
- Homology search against **[UniProt](https://data.broadinstitute.org/Trinity/Trinotate_v3_RESOURCES/uniprot_sprot.pep.gz)** database (**[BlastX](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)** v 2.6.0)
- Homology search against **[non-redundant arthropods](https://ftp.ncbi.nlm.nih.gov/blast/db/v5/nr_v5.*.tar.gz)** database (**[BlastX](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)** v 2.8.0)
- Obtain functional annotation using BlastX UniProt results from **[EggNOG](http://eggnogdb.embl.de/#/app/home)** (Evolutionary Genealogy of Genes: Non-supervised Orthologous Groups), **[KEGG](https://www.kegg.jp)** (Kyoto Encyclopedia of Genes and Genomes), and **[GO](http://geneontology.org)** (Gene Ontology) database (**[Trinotate](https://github.com/Trinotate/Trinotate.github.io/wiki)** v 3.0.2)
3. Differential expression analysis:
- Abundance estimation (**[RSEM](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-323)**)
- Differential expression analysis (**[EdgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)**)
