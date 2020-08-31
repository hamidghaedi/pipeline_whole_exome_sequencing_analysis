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
 In the repositiry there is a ```Snakefile``` file. The content of this file is copied in ```Snakefile.sh``` for inspection purpose, hwoever for pipeline running ```Snakefile``` should be downloaded and located in your ```wes``` directory. 
 
