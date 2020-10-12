# Transcriptomic analysis (RNAseq) of *Macrobrachium rosenbergii* post-larvae in response to *Macrobrachium rosenbergii* nodavirus (*Mr*NV) infection.

## Publication information
Phongthana Pasookhush, Charles Hindmarch, Siwaporn Longyant, Paisarn Sithigorngul, William G. Bendena, Parin Chaivisuthangkura. **Transcriptomic analysis of *Macrobrachium rosenbergii* (giant fresh water prawn) postlarvae in response to *M. rosenbergii* nodavirus (*Mr*NV) infection: *de novo* assembly and functional annotation**. *BMC Genomics* volume 20, Article number: 762 (2019), DOI: 10.1186/s12864-019-6102-6 **[Article URL](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-019-6102-6)**

## Data availibility
Raw data has been uploaded to the National Centre for Biotechnology Information Sequence Read Archive (SRA) under the  BioProject number: PRJNA550272.
## About
This pipeline was used for the transcriptome assembly and differential expression analysis of *M.rosenbergii* in response to *Mr*NV infection.  **[Snakemake](https://snakemake.readthedocs.io/en/stable/)** tool was used to create reproducible and scalable automated pipeline. The pipeline contained three major sections including raw data pre-processing, transcriptome assembly, and post-processing of the transcriptome (See below). 

![alt text](https://github.com/prawnseq/Mrosenbergii_MrNV_RNAseq/blob/master/AnalysisPipeline.png "analysis pipeline")

## Processing steps
Raw data pre-processing:
1. Perform quality assessment using **[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)** v 0.11.5 and **[MultiQC](https://multiqc.info)** v 1.8 to examine **[Phred](https://en.wikipedia.org/wiki/Phred_quality_score)** quality score, GC content, adaptor contamination, size distribution, and N base ratio
2. Perform quality trimming using **[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)** v 0.36 to trim low quality base, remove N base, and discard short read
3. Concatenate all the left and right reads data into universal left and right reads data using merge command

Transcriptome assembly:
1. De novo transcriptome assembly using **[Trinity](https://github.com/trinityrnaseq/trinityrnaseq)** v 2.8.0
2. Remove the redundancy using **[CD-HIT](http://weizhongli-lab.org/cd-hit/)** (95% similarity threshold was used)

Post-processing:
1. Transcriptome quality assessment:
- Calculate fragment mapping rates by mapping reads back to the transcripts using **[Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)** v 2.3.0 
- Examine orthologs completeness against 1,066 complete universal single copy orthologous gene from **[arthropoda_odb9 database](https://busco.ezlab.org/datasets/arthropoda_odb9.tar.gz)** using **[BUSCO](https://busco.ezlab.org)** v 3
- Calculate ExN50 statistics (Top x% most expressed transcripts that have at least N50 length) (See **[details](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Transcriptome-Contig-Nx-and-ExN50-stats)**)
2. Transcriptome annotation:
- Homology search against **[UniProt](https://data.broadinstitute.org/Trinity/Trinotate_v3_RESOURCES/uniprot_sprot.pep.gz)** database using **[BlastX](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)** v 2.6.0
- Homology search against **[non-redundant arthropods](https://ftp.ncbi.nlm.nih.gov/blast/db/v5/nr_v5.*.tar.gz)** database using **[BlastX](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)** v 2.8.0
- Obtain functional annotation using BlastX UniProt results from **[EggNOG](http://eggnogdb.embl.de/#/app/home)** (Evolutionary Genealogy of Genes: Non-supervised Orthologous Groups), **[KEGG](https://www.kegg.jp)** (Kyoto Encyclopedia of Genes and Genomes), and **[GO](http://geneontology.org)** (Gene Ontology) database using **[Trinotate](https://github.com/Trinotate/Trinotate.github.io/wiki)** v 3.0.2 
3. Differential expression analysis:
- Quantify transcripts using alignment based abundance estimation method using **[RSEM](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-323)**
- Perform differential expression analysis using **[EdgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)** (The trimmed mean of M-values normalization method (TMM) )

## Contact informations
**Phongthana Pasookhush, Ph.D.**
- Researcher, Division of Bioinfiormatics and Data Management for Research, Faculty of Medicine, Siriraj Hospital, Mahidol University
- Former Ph.D. student in Biotechnology, Department of Biology, Srinakharinwirot University
- [ResearchGate](https://www.researchgate.net/profile/Phongthana_Pasookhush2)| phongthana.pas@mahidol.edu

**Charles Hindmarch, Ph.D.**
- Adjunct Assistant Professor, Department of Medicine, Queen's University

**Siwaporn Longyant, Ph.D.**
- Associate professor, Department of Biology, Srinakharinwirot University

**Paisarn Sithigorngul, Ph.D.**
- Professor, Department of Biology, Srinakharinwirot University

**William Bendena, Ph.D.**
- Professor of Molecular Physiology and Behaviour, Department of Biology, Queen's University

**Parin Chaivisuthangkura, Ph.D.**
- Associate professor, Department of Biology, Srinakharinwirot University
