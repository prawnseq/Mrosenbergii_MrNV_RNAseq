#### Command ####
# snakemake -np
# snakemake --dag | dot -Tsvg > dag.svg
# snakemake --report report.html

#### Load config file ####
configfile: "config.yaml"

#### Set up from config file ####
# Software #
fastqc = config["software"]["fastqc"]
multiqc = config["software"]["multiqc"]
trimmomatic = config["software"]["trimmomatic"]
trinity = config["software"]["trinity"]
cdhit = config["software"]["cdhit"]
bowtie2build = config["software"]["bowtie2build"]
bowtie2 = config["software"]["bowtie2"]
busco = config["software"]["busco"]
rsem =  config["software"]["rsem"]
trinitystat = config["software"]["trinitystat"]
exn50 = config["software"]["exn50"]
exn50plot = config["software"]["exn50plot"]
matrix = config["software"]["matrix"]
deg = config["software"]["deg"]
pca = config["software"]["pca"]
heatmap = config["software"]["heatmap"]
blastx = config["software"]["blastx"]
blastxnr = config["software"]["blastxnr"]
blastp = config["software"]["blastp"]
hmmscan = config["software"]["hmmscan"]
signalp = config["software"]["signalp"]
tmhmm = config["software"]["tmhmm"]
transdecoder_longorfs = config["software"]["transdecoder_longorfs"]
transdecoder_predict = config["software"]["transdecoder_predict"]
genetransmap = config["software"]["genetransmap"]
trinotate = config["software"]["trinotate"]
trinotategoextract = config["software"]["trinotategoextract"]

#### Set up output dir ####
fastqc_raw_out="output/fastqc_raw/"
multiqc_raw_out="output/multiqc_raw/"
trimmomatic_out="output/trim/"
fastqc_trim_out="output/fastqc_trim/"
multiqc_trim_out="output/multiqc_trim/"
merge_out="output/merge/"
trinity_out="output/trinity/"
cdhit_out="output/cdhit/"
bowtie2_out="output/bowtie2/"
rsem_out="output/rsem/"
exn50stat_out="output/exn50stat/"
deg_out="output/deg_edgeR/"
blastx_trinotate_out="output/blastx_trinotate/"
blastx_nr_out="output/blastx_nr/"
blastp_trinotate_out="output/blastp_trinotate/"

log_out="output/log/"

#### all ####
rule target:
    input:
        multiqc_raw_out + "multiqc_report.html",
        multiqc_trim_out + "multiqc_report.html",
        directory("run_BUSCO/"),
        bowtie2_out + "align_stats.txt",
        exn50stat_out + "ExN50.stats.plot.pdf",
        exn50stat_out + "TrinityStats100.txt",
        exn50stat_out + "TrinityStats95.txt",
        "Trinotate.xls",
        directory(deg_out),
        blastx_nr_out + "nr.blastx.outfmt6"

#### FASTQC on raw reads ####
rule fastqc_raw:
    input:
        lambda wildcards: config["samplesname"][wildcards.sample]
    output:
        fastqc_raw_out+ "{sample}/{sample}.f_fastqc.zip",
        fastqc_raw_out+ "{sample}/{sample}.r_fastqc.zip"
    log:
        log_out + "fastqc_raw.{sample}.log"
    threads: config["threads"]["fastqc"]
    shell:
        "{fastqc} {input} -o {fastqc_raw_out}{wildcards.sample}/ -t {threads}"

#### MULTIQC on fastqc_raw ####
rule multiqc_raw:
    input:
        expand("output/fastqc_raw/{sample}/{sample}.f_fastqc.zip", sample=config["samplesname"]),
        expand("output/fastqc_raw/{sample}/{sample}.r_fastqc.zip", sample=config["samplesname"])
    output:
        multiqc_raw_out + "multiqc_report.html"
    log:
        log_out + "multiqc_raw.log"
    shell:
        "{multiqc} --force {input} -o {multiqc_raw_out} "

#### Trim adaptor on raw reads using Trimmomatic ####
rule trimmomatic:
    input:
        lambda wildcards: config["samplesname"][wildcards.sample]
    output:
        trimmomatic_out + "{sample}/{sample}.forward.paired.fastq.gz",
        trimmomatic_out + "{sample}/{sample}.reverse.paired.fastq.gz"
    log:
        log_out + "trimmomatic.{sample}.log"
    threads: config["threads"]["trimmomatic"]
    params: config["trimmomatic_params"]
    shell:
        "{trimmomatic} PE -threads {threads} {input} \
        {trimmomatic_out}{wildcards.sample}/{wildcards.sample}.forward.paired.fastq.gz \
        {trimmomatic_out}{wildcards.sample}/{wildcards.sample}.forward.unpaired.fastq.gz \
        {trimmomatic_out}{wildcards.sample}/{wildcards.sample}.reverse.paired.fastq.gz \
        {trimmomatic_out}{wildcards.sample}/{wildcards.sample}.reverse.unpaired.fastq.gz \
        {params}"
        ## Quality trimmed reads in paired (P)
        ## Quality trimmed paired reads which one was discarded and one was kept (U)

#### FASTQC on trim reads ####
rule fastqc_trim:
    input:
        trimmomatic_out + "{sample}/{sample}.forward.paired.fastq.gz",
        trimmomatic_out + "{sample}/{sample}.reverse.paired.fastq.gz"
    output:
        fastqc_trim_out + "{sample}.forward.paired_fastqc.zip",
        fastqc_trim_out + "{sample}.reverse.paired_fastqc.zip",
    log:
        log_out + "fastqc_trim.{sample}.log"
    threads: config["threads"]["fastqc"]
    shell:
        "{fastqc} {trimmomatic_out}{wildcards.sample}/{wildcards.sample}.forward.paired.fastq.gz \
        {trimmomatic_out}{wildcards.sample}/{wildcards.sample}.reverse.paired.fastq.gz -o {fastqc_trim_out} -t {threads}"

#### MULTIQC on fastqc_trim ####
rule multiqc_trim:
    input:
        expand("output/fastqc_trim/{sample}.forward.paired_fastqc.zip", sample=config["samplesname"]),
        expand("output/fastqc_trim/{sample}.reverse.paired_fastqc.zip", sample=config["samplesname"])
    output:
        multiqc_trim_out + "multiqc_report.html"
    log:
        log_out + "multiqc_trim.log"
    shell:
        "{multiqc} --force {input} -o {multiqc_trim_out} "

#### Merge trim left reads for Trinity ####
rule merge_left:
    input:
        left=expand("output/trim/{sample}/{sample}.forward.paired.fastq.gz", sample=config["samplesname"])
    output:
        merge_out + "merged.forward.paired.fastq.gz"
    log:
        log_out + "merge_left.log"
    shell:
        "cat {input} > {output}"

#### Merge trim right reads for Trinity ####
rule merge_right:
    input:
        right=expand("output/trim/{sample}/{sample}.reverse.paired.fastq.gz", sample=config["samplesname"])
    output:
        merge_out + "merged.reverse.paired.fastq.gz"
    log:
        log_out + "merge_right.log"
    shell:
        "cat {input} > {output}"

#### de novo assembly using Trinity ####
rule trinity:
    input:
        left= merge_out + "merged.forward.paired.fastq.gz",
        right= merge_out + "merged.reverse.paired.fastq.gz"
    output:
        trinity_out + "Trinity.fasta"
    threads: config["threads"]["trinity"]
    params:
        max_memory = config["trinity_params"]["max_memory"],
        SS_lib_type = config["trinity_params"]["SS_lib_type"]
    log:
        log_out + "trinity.log"
    shell:
        "{trinity} --seqType fq --CPU {threads} \
        --max_memory {params.max_memory} --SS_lib_type {params.SS_lib_type} \
        --left {input.left} \
        --right {input.right} \
        --output {trinity_out}"

#### remove redundancy using CD-HIT ####
rule cdhit:
    input:
        transcriptome100= trinity_out + "Trinity.fasta",
    output:
        cdhit_out + "Trinity95.fasta"
    threads: config["threads"]["cdhit"]
    params:
        c = config["cdhit_params"]["c"]
    log:
        log_out + "cdhit.log"
    shell:
        "{cdhit} -i {input} -o {output} -c {params.c} -p 1 -d 0 -b 3 -T {threads}"

#### transcriptome quality assessment using Bowtie2 ####
rule bowtie2:
    input:
        transcriptome95= cdhit_out + "Trinity95.fasta",
    output:
        bowtie2_out + "align_stats.txt"
    params:
        left= merge_out + "merged.forward.paired.fastq.gz",
        right= merge_out + "merged.reverse.paired.fastq.gz"
    log:
        log_out + "bowtie2.log"
    shell:
        """
        {bowtie2build} {input} {input}

        {bowtie2} -p 10 -q --no-unal -k 20 \
        -x {input} -1 {params.left} -2 {params.right}  \
        2>{bowtie2_out}/align_stats.txt| samtools view -@10 -Sb -o bowtie2.bam

        cat 2>&1 {bowtie2_out}/align_stats.txt
        """

#### transcriptome quality assessment using BUSCO ####
    #### transcripts completeness using Arthropoda odb9 database from busco website ####
rule busco:
    input:
        transcriptome95= cdhit_out + "Trinity95.fasta",
    output:
        directory("run_BUSCO/")
    ## run_XXXXX, XXXXX has to match config["busco_params"]["out_name"] ##
    threads: config["threads"]["busco"]
    params:
        mode = config["busco_params"]["mode"],
        lineage = config["busco_params"]["lineage"],
        out_name = config["busco_params"]["out_name"]
    log:
        log_out + "busco.log"
    shell:
        "{busco} --in {input} --out {params.out_name} --force\
        --lineage_path {params.lineage} --mode {params.mode} \
        --cpu {threads}"

#### Abundance estimation using RSEM ####
rule rsem:
    input:
        lambda wildcards: config["samplesname"][wildcards.sample],
        transcriptome95= cdhit_out + "Trinity95.fasta",
    output:
        "output/rsem/{sample}/RSEM.isoforms.results"
    log:
        log_out + "rsem.{sample}.log"
    shell:
        "{rsem} --seqType fq --est_method RSEM --aln_method bowtie2 \
        --trinity_mode  --prep_reference\
        --transcripts {input.transcriptome95} \
        --left {trimmomatic_out}{wildcards.sample}/{wildcards.sample}.forward.paired.fastq.gz \
        --right {trimmomatic_out}{wildcards.sample}/{wildcards.sample}.reverse.paired.fastq.gz \
        --output_dir {rsem_out}{wildcards.sample}"

rule matrix:
    input:
        expand("output/rsem/{sample}/RSEM.isoforms.results", sample=config["samplesname"])
    output:
        "RSEM.isoform.counts.matrix",
        "RSEM.isoform.TMM.EXPR.matrix"
    log:
        log_out + "matrix.log"
    shell:
        "{matrix} --est_method RSEM  --gene_trans_map none --name_sample_by_basedir {input}"

#### transcriptome quality assessment using ExN50 and TrinityStats ####
rule exn50stat:
    input:
        transcriptome100= trinity_out + "Trinity.fasta",
        transcriptome95= cdhit_out + "Trinity95.fasta",
        matrix= "RSEM.isoform.TMM.EXPR.matrix"
    output:
        exn50stat_out + "ExN50.stats.plot.pdf",
        exn50stat_out + "TrinityStats100.txt",
        exn50stat_out + "TrinityStats95.txt"
    log:
        log_out + "exn50stat.log"
    shell:
        """
        {trinitystat} {input.transcriptome100} > {exn50stat_out}TrinityStats100.txt

        {trinitystat} {input.transcriptome95} > {exn50stat_out}TrinityStats95.txt

        {exn50} {input.matrix} {input.transcriptome95} | tee {exn50stat_out}ExN50.stats

        {exn50plot} {exn50stat_out}ExN50.stats
        """

#### Perform differential expression analysis using edgeR ####
#### Extract DEG and generate PCA plot ####
rule deg:
    input:
        counts="RSEM.isoform.counts.matrix",
        tmm="RSEM.isoform.TMM.EXPR.matrix"
    output:
        directory(deg_out)
    params:
        samples = config["deg_params"]["samples"],
        pc = config["deg_params"]["pc"],
        p = config["deg_params"]["p"],
        c = config["deg_params"]["c"],
        maxDEgene = config["deg_params"]["maxDEgene"]
    log:
        log_out + "deg.log"
    shell:
        """
        {deg} --matrix {input.counts} --method edgeR --samples_file {params.samples}\
        --output {output}

        {pca} --matrix {input.counts} --samples {params.samples} \
        --min_rowSums 10 --log2 --CPM --center_rows --prin_comp {params.pc}

        cd {output}

        {heatmap} --matrix {input.tmm} --samples {params.samples} \
        -P {params.p} -C {params.c} --max_DE_genes_per_comparison {params.maxDEgene}

        cd ~
        """
#### Extract potential protein-coding regions in transcripts ####
#### using Transdecoder longorfs and predict function ####
rule transdecoder:
    input:
        transcriptome95= cdhit_out + "Trinity95.fasta",
    output:
        "Trinity95.fasta.transdecoder.pep"
    threads: config["threads"]["transdecoder"]
    params:
        min_length = config["transdecoder_params"]["min_length"]
    log:
        log_out + "transdecoder.log"
    shell:
        """
        {transdecoder_longorfs} -t {input} -m {params.min_length} &> {log}
        {transdecoder_predict} -t {input} --cpu {threads} &>> {log}
        {genetransmap} {input} > {trinity_out}Trinity.fasta.gene_trans_map
        """
#### Transcript annotation using BLASTx 2.8 against Nr Arthropod database ####
rule blastx_nr:
    input:
        transcriptome95= cdhit_out + "Trinity95.fasta",
    output:
        blastx_nr_out + "nr.blastx.outfmt6"
    threads: config["threads"]["blastx"]
    params:
        database = config["blastx_nr_params"]["database"],
        evalue = config["blastx_nr_params"]["evalue"],
        max_target_seqs = config["blastx_nr_params"]["max_target_seqs"],
        outfmt = config["blastx_nr_params"]["outfmt"],
        taxidlist = config["blastx_nr_params"]["taxidlist"]
    log:
        log_out + "blastx_nr.log"
    shell:
#### MUST DOWNLOAD NR DATABASE V5 FIRST !!! ####
        "{blastxnr} -query {input} -db {params.database} \
        -out {blastx_nr_out}nr.blastx.outfmt6 \
        -evalue {params.evalue} -num_threads {threads} \
        -taxidlist {params.taxidlist}\
        -max_target_seqs {params.max_target_seqs} -outfmt '{params.outfmt}'"

#### Trinotate transcripts annotation using BLASTx against Uniprot_pep ####
rule blastx_trinotate:
    input:
        transcriptome95= cdhit_out + "Trinity95.fasta",
    output:
        blastx_trinotate_out + "swissprot.blastx.outfmt6"
    threads: config["threads"]["blastx"]
    params:
        database = config["blastx_trinotate_params"]["database"],
        evalue = config["blastx_trinotate_params"]["evalue"],
        max_target_seqs = config["blastx_trinotate_params"]["max_target_seqs"],
        outfmt = config["blastx_trinotate_params"]["outfmt"]
    log:
        log_out + "blastx_trinotate.log"
    shell:
#### MUST PERFORM makeblastdb FIRST !!! ####
## using this command "makeblastdb -in PATHtoDB -dbtype prot" ##
#### wget 'https://data.broadinstitute.org/Trinity/Trinotate_v3_RESOURCES/uniprot_sprot.pep.gz' ####
        "{blastx} -query {input} -db {params.database} \
        -out {blastx_trinotate_out}swissprot.blastx.outfmt6 \
        -evalue {params.evalue} -num_threads {threads} \
        -max_target_seqs {params.max_target_seqs} -outfmt {params.outfmt}"

#### Trinotate transcripts annotation using BLASTp against Uniprot_pep ####
rule blastp_trinotate:
    input:
        "Trinity95.fasta.transdecoder.pep",
    output:
        blastp_trinotate_out + "swissprot.blastp.outfmt6"
    threads: config["threads"]["blastp"]
    params:
        database = config["blastx_trinotate_params"]["database"],
        evalue = config["blastx_trinotate_params"]["evalue"],
        max_target_seqs = config["blastx_trinotate_params"]["max_target_seqs"],
        outfmt = config["blastx_trinotate_params"]["outfmt"]
        ### USE SAME PARAMS AS BLASTx ###
    log:
        log_out + "blastp_trinotate.log"
    shell:
        "{blastp} -query {input} -db {params.database} \
        -out {blastp_trinotate_out}swissprot.blastp.outfmt6 \
        -evalue {params.evalue} -num_threads {threads} \
        -max_target_seqs {params.max_target_seqs} -outfmt {params.outfmt}"

#### Protein family alignment using hmmscan against Pfam database ####
rule hmmscan:
    input:
        "Trinity95.fasta.transdecoder.pep",
    output:
        "TrinotatePFAM.out"
    threads: config["threads"]["hmmscan"]
    params:
        database = config["hmmscan_params"]["database"],
    shell:
#### MUST PERFORM "hmmpress Pfam-A.hmm" FIRST !! ####
#### wget 'https://data.broadinstitute.org/Trinity/Trinotate_v3_RESOURCES/Pfam-A.hmm.gz'
        "{hmmscan} --cpu {threads} --domtblout {output} \
        {params.database} {input}"

#### Signal sequence prediction using SipnalP ####
rule signalp:
    input:
        "Trinity95.fasta.transdecoder.pep",
    output:
        "signalp.out"
    shell:
        "{signalp} -f short -n signalp.out {input}"

#### Transmembrane sequence prediction using tmhmm ####
rule tmhmm:
    input:
        "Trinity95.fasta.transdecoder.pep",
    output:
        "tmhmm.out"
    shell:
        "{tmhmm} --short < {input} > tmhmm.out"

#### Obtain the Trinotate.xls annotation file ####
rule trinotate:
    input:
        blastx=blastx_trinotate_out + "swissprot.blastx.outfmt6",
        blastp=blastp_trinotate_out + "swissprot.blastp.outfmt6",
        hmmscan="TrinotatePFAM.out",
        signalp="signalp.out",
        tmhmm="tmhmm.out"
    output:
        "Trinotate.xls"
    params:
        transcripts=cdhit_out + "Trinity95.fasta",
        pep="Trinity95.fasta.transdecoder.pep",
        genetransmap=trinity_out + "Trinity.fasta.gene_trans_map",
        sqlite = config["trinotate_params"]["sqlite"]
    shell:
#### MUST HAVE Trinotate_V3.sqlite extracted in home_dir FIRST !!!! ####
#### wget 'https://data.broadinstitute.org/Trinity/Trinotate_v3_RESOURCES/Trinotate_v3.sqlite.gz'
        """
        {trinotate} {params.sqlite} init \
        --gene_trans_map {params.genetransmap} \
        --transcript_fasta {params.transcripts} \
        --transdecoder_pep {params.pep}

        {trinotate} {params.sqlite} LOAD_swissprot_blastx {input.blastx}

        {trinotate} {params.sqlite} LOAD_swissprot_blastp {input.blastp}

        {trinotate} {params.sqlite} LOAD_pfam {input.hmmscan}

        {trinotate} {params.sqlite} LOAD_signalp {input.signalp}

        {trinotate} {params.sqlite} LOAD_tmhmm {input.tmhmm}

        {trinotate} {params.sqlite} report > Trinotate.xls

        {trinotategoextract} --Trinotate_xls Trinotate.xls -G > go_annotations.txt
        """
