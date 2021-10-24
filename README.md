# Whole-Exome sequencing data analysis

#### Note
This code/pipeline was written by me in 2018-2019 while working in the [RayaGen startup](http://rayagn.com/). The main field where it could be employed is automating whole-exome sequencing data analysis from raw (fastq) to variant call format (vcf). This vcf could be further pass to a variant annotation tool of your choice. This pipeline was written by ``` snakemake``` workflow management system.  This is how original website describe ``` snakemake```  
>The Snakemake workflow management system is a tool to create reproducible and scalable data analyses. Workflows are described via a human readable, Python based language. They can be seamlessly scaled to server, cluster, grid and cloud environments, without the need to modify the workflow definition. Finally, Snakemake workflows can entail a description of required software, which will be automatically deployed to any execution environment.
For full documentation and installation please refere to the original [website](https://snakemake.readthedocs.io/en/stable/).

In this workflow, several address and path was defined in the ```snakefile```:
```shell
SAMPLES, = glob_wildcards("raw/{sample}_1.fastq.gz") # your file should be in a folder in this way: "/wes/raw"
READS=["1", "2"]
BASE_DIR = "/path/to/wes"
ADAPTORS = BASE_DIR + "/software/trimmomatic-0.38/adapters/TruSeq3-PE.fa"
REF_GENOME = "/path/to/wes/ref/human_g1k_v37_decoy.fasta"
BWA_INDEX = BASE_DIR + "/ref"
KNOWEN_SITE_SNP = BASE_DIR + "/ref/dbsnp_138.b37.vcf"
KNOWEN_SITE_INDEL = BASE_DIR + "/ref/Mills_and_1000G_gold_standard.indels.b37.vcf"
GATK3 = "/path/to/wes/software/GenomeAnalysisTK.jar"
PICARD = "/path/to/wes/software/picard.jar"



rule all:
  input:  expand("VCF/{sample}_filtered_final.vcf",sample=SAMPLES),
          #expand("VCF/{sample}_filtered_INDEL.vcf",sample=SAMPLES),
          #expand("VCF/{sample}_filtered_SNP.vcf",sample=SAMPLES),
          #expand("VCF/{sample}_raw_INDEL.vcf",sample=SAMPLES),
          #expand("VCF/{sample}_raw_SNP.vcf",sample=SAMPLES),
          #expand("VCF/{sample}_raw.vcf",sample=SAMPLES),
          #expand("BQSR/recal_reads_{sample}.bam",sample=SAMPLES),
          #"BQSR/recal_data.table",
          #expand("alignment/marked_duplicates_{sample}.bam",sample=SAMPLES),
          #expand("alignment/sorted_{sample}.bam",sample=SAMPLES),
         # expand("alignment/{sample}.bam",sample=SAMPLES),
          #expand("interleaved_fastq/{sample}_interleaved.fastq.gz", sample=SAMPLES),
          expand("multi_qc/{sample}.multiqc.html", sample=SAMPLES),
          #expand("seqs_trimmed/{sample}_trimmed_paired_{read}.fastq.gz", sample=SAMPLES, read=READS),
          expand("qc_trimmed/{sample}_trimmed_paired_{read}_fastqc.html", sample=SAMPLES, read=READS),
          expand("qc_trimmed/{sample}_trimmed_paired_{read}_fastqc.zip", sample=SAMPLES, read=READS),
          #expand("seqs_trimmed/{sample}_trimmed_paired_R1.fastq.gz", sample=SAMPLES),
          #expand("seqs_trimmed/{sample}_trimmed_unpaired_R1.fastq.gz", sample=SAMPLES),
          #expand("seqs_trimmed/{sample}_trimmed_paired_R2.fastq.gz", sample=SAMPLES),
          #expand("seqs_trimmed/{sample}_trimmed_unpaired_R2.fastq.gz", sample=SAMPLES),
          expand("qc_raw/{sample}_{read}_fastqc.html", sample=SAMPLES, read=READS),
          expand("qc_raw/{sample}_{read}_fastqc.zip", sample=SAMPLES, read=READS),

rule fastqc_raw:
    input:  expand("raw/{sample}_{read}.fastq.gz", sample=SAMPLES, read=READS)
    output: expand("qc_raw/{sample}_{read}_fastqc.html", sample=SAMPLES, read=READS),
            expand("qc_raw/{sample}_{read}_fastqc.zip", sample=SAMPLES, read=READS)
    log:    expand("logs/fastqc/{sample}_{read}_raw.fastqc.err", sample=SAMPLES, read=READS)
    priority:1
    message: "Executing Quality check (FASTQC)  on the following files {input}."
    shell: """
      mkdir -p logs/fastqc
      mkdir -p qc_raw
      fastqc --outdir qc_raw  --thread 50 --nogroup {input} 2>>{log}
        """
rule trimming:
     input:
        fwd = "raw/{sample}_R1.fastq.gz",
        rev = "raw/{sample}_R2.fastq.gz",
     output:
        fwd_paired = "seqs_trimmed/{sample}_trimmed_paired_R1.fastq.gz",
        fwd_unpaired = "seqs_trimmed/{sample}_trimmed_unpaired_R1.fastq.gz",
        rev_paired = "seqs_trimmed/{sample}_trimmed_paired_R2.fastq.gz",
        rev_unpaired = "seqs_trimmed/{sample}_trimmed_unpaired_R2.fastq.gz"
     params:
       TRIMMOMATIC_JAR = "~/RayaGene/wes/software/trimmomatic-0.38/trimmomatic-0.38.jar",
       PAIR = "PE SE".split(),
       threads=50,
       TRIM_OPTS = "-phred33 ILLUMINACLIP:" + ADAPTORS + ":2:30:10 LEADING:3 TRAILING:3  SLIDINGWINDOW:4:20 MINLEN:36"
     log:  "logs/trimmomatic/{sample}.paired.trimomatic.log"
     priority:2
     message: "Executing Trimming (trimmomatic) on the following files {input.fwd} and {input.rev}."
     shell: """
        mkdir -p logs/trimmomatic
        java -jar {params.TRIMMOMATIC_JAR} \
           {params.PAIR[0]} -threads {params.threads} {input.fwd} {input.rev} \
           {output.fwd_paired} {output.fwd_unpaired} \
           {output.rev_paired} {output.rev_unpaired} \
           {params.TRIM_OPTS}  2>{log}
        """
rule fastqc_trimmed:
    input:  expand("seqs_trimmed/{sample}_trimmed_paired_{read}.fastq.gz", sample=SAMPLES, read=READS)
    output: expand("qc_trimmed/{sample}_trimmed_paired_{read}_fastqc.html", sample=SAMPLES, read=READS),
            expand("qc_trimmed/{sample}_trimmed_paired_{read}_fastqc.zip", sample=SAMPLES, read=READS)
    log:    expand("logs/fastqc/{sample}_{read}_trimmed.fastqc.err", sample=SAMPLES, read=READS)
    priority:3
    message: "Executing Quality Checks (FASTQC)  on the following files {input}."
    shell: """
       fastqc --outdir qc_trimmed --thread 80 --nogroup {input} \
             >> qc_trimmed/fastqc.log 2>>{log}
          """
rule multiqc:
    input:  raw_fwd=    expand("qc_raw/{sample}_R1_fastqc.zip", sample= SAMPLES ),
            raw_rev=    expand("qc_raw/{sample}_R2_fastqc.zip", sample= SAMPLES ),
            trimed_fwd= expand("qc_trimmed/{sample}_trimmed_paired_R1_fastqc.zip", sample= SAMPLES ),
            trimed_rev= expand("qc_trimmed/{sample}_trimmed_paired_R2_fastqc.zip", sample= SAMPLES ),
    output: expand("multi_qc/{sample}.multiqc.html", sample=SAMPLES)
    log:    expand("logs/multi_qc/{sample}.log", sample=SAMPLES)
    priority:4
    shell:"""
        mkdir -p multi_qc
        multiqc {input.raw_fwd} {input.raw_rev} {input.trimed_fwd} {input.trimed_rev} > {output} 2>{log}
        mv multiqc_report.html ./multi_qc
        mv multiqc_data multi_qc
        """
rule seqtk:
    input:  fwd= expand("seqs_trimmed/{sample}_trimmed_paired_R1.fastq.gz", sample= SAMPLES ),
            rev= expand("seqs_trimmed/{sample}_trimmed_paired_R2.fastq.gz", sample= SAMPLES ),
    output: expand("interleaved_fastq/{sample}_interleaved.fastq.gz", sample=SAMPLES),
    log:    expand("logs/seqtk/{sample}.log", sample=SAMPLES),
    priority:5
    shell:"""
        mkdir -p interleaved_fastq
       seqtk mergepe {input.fwd} {input.rev} > {output} 2>{log}
       """
rule bwa_map:
    input:
        ref=REF_GENOME,
        read=expand("interleaved_fastq/{sample}_interleaved.fastq.gz", sample=SAMPLES),
    output:    bam = expand("alignment/{sample}.bam",sample=SAMPLES),
               #sorted = expand("alignment/sorted_{sample}.bam",sample=SAMPLES),
               #bai = expand("alignment/sorted_{sample}.bam.bai",sample=SAMPLES),
    log:expand("logs/alignment/{sample}.log", sample=SAMPLES),
    params:
    rg = expand("'@RG\\tID:C6C0TANXX_2\\tSM:ZW177\\tLB:ZW177lib\\tPL:ILLUMINA'")
    priority:6
    shell:
       r"""
           mkdir -p alignment
           mkdir -p logs/alignment
           bwa mem -M -t 50 -p {input.ref} {input.read} -R {params.rg}\
            |samtools view -Sb - > {output.bam} 2>log
            """
rule bam_sort_index:
            input:expand("alignment/{sample}.bam",sample=SAMPLES),
            output:sorted = expand("alignment/sorted_{sample}.bam",sample=SAMPLES),
                   bai = expand("alignment/sorted_{sample}.bam.bai",sample=SAMPLES),
            shell:
               r"""
               samtools sort {input} > {output.sorted}
               samtools index {output.sorted} > {output.bai}
               """
rule mark_dup:
    input: expand("alignment/sorted_{sample}.bam",sample=SAMPLES),
    output:    bam = expand("alignment/marked_duplicates_{sample}.bam",sample=SAMPLES),
               metrics ="alignment/marked_dup_metrics.txt",
               bai =expand("alignment/marked_duplicates_{sample}.bam.bai",sample=SAMPLES),
    log:expand("logs/mark_dup/{sample}.log", sample=SAMPLES),
    params:
       PICARD_JAR = PICARD,
    threads: 8
    priority:7
    shell:
       r"""
            mkdir -p logs/mark_dup
            java -jar -Xmx64G {params.PICARD_JAR} MarkDuplicates I={input} O={output.bam} M={output.metrics} TMP_DIR=/tmp
            samtools index {output.bam} > {output.bai}
            """
rule BQSR_1:
    input: bam = expand("alignment/marked_duplicates_{sample}.bam",sample=SAMPLES),
           ref= REF_GENOME,
           known_site_snp= KNOWEN_SITE_SNP,
           known_site_indel = KNOWEN_SITE_INDEL,
    output:"BQSR/recal_data.table",
    log:expand("logs/BQSR/{sample}.log", sample=SAMPLES),
    params:
       GATK3_JAR = GATK3,
    priority:8
    shell:
       r"""
            mkdir -p logs/BQSR
            java -jar -Xmx64G {params.GATK3_JAR} -nct 32 -T BaseRecalibrator \
            -R {input.ref} -I {input.bam} \
            -knownSites {input.known_site_snp} \
            -knownSites {input.known_site_indel} \
            -o {output}
            """
rule BQSR_2:
    input: bam = expand("alignment/marked_duplicates_{sample}.bam",sample=SAMPLES),
           ref= REF_GENOME,
           recal_table = "BQSR/recal_data.table",
    output:expand("BQSR/recal_reads_{sample}.bam",sample=SAMPLES),
    log:expand("logs/BQSR/{sample}_2.log", sample=SAMPLES),
    params:
       GATK3_JAR = GATK3,
    priority:9
    shell:
       r"""
            java -jar -Xmx64G {params.GATK3_JAR} -nct 64 -T PrintReads \
            -R {input.ref} \
            -I {input.bam} \
            -BQSR {input.recal_table} \
            -o {output} 2>log
            """
rule Variant_calling:
    input: bam = expand("BQSR/recal_reads_{sample}.bam",sample=SAMPLES),
           ref= REF_GENOME,
           known_site_snp= KNOWEN_SITE_SNP,
           recal_table = "BQSR/recal_data.table",
           known_site_indel = KNOWEN_SITE_INDEL
    output:expand("VCF/{sample}_raw.vcf",sample=SAMPLES),
    log:expand("logs/VCF/{sample}_raw_vcf.log", sample=SAMPLES),
    params:
       GATK3_JAR = GATK3,
    priority:10
    shell:
       r"""
            java -jar -Xmx80G {params.GATK3_JAR} -nct 64 -T HaplotypeCaller  \
            -R {input.ref} \
            -I {input.bam} \
            -stand_call_conf 30 \
            --genotyping_mode DISCOVERY \
            -o {output} 2>log
            """
rule Select_SNP:
    input: vcf = expand("VCF/{sample}_raw.vcf",sample=SAMPLES),
           ref= REF_GENOME,
    output:expand("VCF/{sample}_raw_SNP.vcf",sample=SAMPLES),
    log:expand("logs/VCF/{sample}_raw_SNP_vcf.log", sample=SAMPLES),
    params:
       GATK3_JAR = GATK3,
    threads: 8
    priority:11
    shell:
       r"""
            java -jar -Xmx80G {params.GATK3_JAR} -T SelectVariants  \
            -R {input.ref} \
            -V {input.vcf} \
            -selectType SNP \
            -o {output} 2>log
            """
rule Select_INDEL:
    input: vcf = expand("VCF/{sample}_raw.vcf",sample=SAMPLES),
           ref= REF_GENOME,
    output:expand("VCF/{sample}_raw_INDEL.vcf",sample=SAMPLES),
    log:expand("logs/VCF/{sample}_raw_INDEL_vcf.log", sample=SAMPLES),
    params:
       GATK3_JAR = GATK3,
    threads: 8
    priority:12
    shell:
       r"""
            java -jar -Xmx80G {params.GATK3_JAR} -T SelectVariants  \
            -R {input.ref} \
            -V {input.vcf} \
            -selectType INDEL \
            -o {output} 2>log
            """
rule SNP_FILTER:
    input: vcf = expand("VCF/{sample}_raw_SNP.vcf",sample=SAMPLES),
           ref= REF_GENOME,
    output:expand("VCF/{sample}_filtered_SNP.vcf",sample=SAMPLES),
    log:expand("logs/VCF/{sample}_filtered_SNP_vcf.log", sample=SAMPLES),
    params:
       GATK3_JAR = GATK3,
    threads: 8
    priority:13
    shell:
       r"""
            java -jar -Xmx80G {params.GATK3_JAR} -T VariantFiltration  \
            -R {input.ref} \
            -V {input.vcf} \
            --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
            --filterName "Low_Qual" \
            -o {output} 2>log
            """
rule INDEL_FILTER:
    input: vcf = expand("VCF/{sample}_raw_INDEL.vcf",sample=SAMPLES),
           ref= REF_GENOME,
    output:expand("VCF/{sample}_filtered_INDEL.vcf",sample=SAMPLES),
    log:expand("logs/VCF/{sample}_filtered_INDEL_vcf.log", sample=SAMPLES),
    params:
       GATK3_JAR = GATK3,
    threads: 8
    priority:13
    shell:
       r"""
            java -jar -Xmx80G {params.GATK3_JAR} -T VariantFiltration  \
            -R {input.ref} \
            -V {input.vcf} \
            --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
            --filterName "Low_Qual" \
            -o {output} 2>log
            """
rule combine_variants:
    input:snp =   expand("VCF/{sample}_filtered_SNP.vcf",sample=SAMPLES),
          indel = expand("VCF/{sample}_filtered_INDEL.vcf",sample=SAMPLES),
           ref= REF_GENOME,
    output:expand("VCF/{sample}_filtered_final.vcf",sample=SAMPLES),
    log:expand("logs/VCF/{sample}_filtered_final_vcf.log", sample=SAMPLES),
    params:
       GATK3_JAR = GATK3,
    priority:13
    shell:
       r"""
            java -jar -Xmx80G {params.GATK3_JAR} -T CombineVariants  \
            -R {input.ref} \
            -V {input.snp} \
            -V {input.indel} \
            --genotypemergeoption UNSORTED \
            -o {output} 2>log
            
```
There is a configuration ```yaml``` file in this directory that could be of help in making same virtual environment in ```python```to what I had configured, ```wes```. By ```conda``` you can make same environment to ```wes```. 
 ### Pipeline steps
 This pipeline consisted of several steps:
 * QC by [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
 * read trimming by [trimmomatics](https://github.com/timflutre/trimmomatic)
 * summarizing qc parameters by [multiqc](https://multiqc.info/)
 * read mapping by [bwa](http://bio-bwa.sourceforge.net/)
 * marking duplicated by [picard](https://broadinstitute.github.io/picard/)
 * base quality recalibration by [gatk](https://gatk.broadinstitute.org/hc/en-us)
 * variant calling by [gatk](https://gatk.broadinstitute.org/hc/en-us)
 

 *fnal note:*
 In the repositiry there is a ```Snakefile``` file. The content of this file is copied in ```Snakefile.sh``` for inspection purpose, however for pipeline running ```Snakefile``` should be downloaded and located in your ```wes``` directory. 
 
