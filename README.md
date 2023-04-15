# Whole-genome sequencing data analysis for calling genomic variants
Variant calling pipeline for identifying gremline and somatic mutations in patient genomes. 

## Preparations 
### Download human reference genome
cd hg_wgs_variant_calling

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz inputs/fasta/ & unzip inputs/hg38.fa.gz

### Download vcf files for human known variants
1. Visit https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false

2. Download files into inputs/vcf/ folder:

  1000G_omni2.5.hg38.vcf.gz
  1000G_omni2.5.hg38.vcf.gz.tbi

  1000G_phase1.snps.high_confidence.hg38.vcf.gz

1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi

Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz

Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi

Homo_sapiens_assembly38.dbsnp138.vcf.gz

Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi

Homo_sapiens_assembly38.known_indels.vcf.gz

Homo_sapiens_assembly38.known_indels.vcf.gz.tbi

Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi

hapmap_3.3.hg38.vcf.gz

hapmap_3.3.hg38.vcf.gz.tbi


### install conda environment
conda env create -f envs/environment.yaml -p envs/smake_wgs

conda activate envs/smake_wgs

## Run pipeline
The pipeline is writen as a snakemake workflow, and can be run within HPC cluster (default is using sbatch command for submission of jobs). 
To test the pipeline, do a dry run: 

cd hg_wgs_variant_calling

snakemake -n


To run in cluster with nohup:

nohup snakemake --profile profile &

## Custom script for calculating edit distance between variant site and gRNA
The custom script "off_target.pl" takes and input file containing a list of variants (eg. called from the pipeline above). The input file should contain 5 columns: 1. chromosome id; 2. variant start position; 3. variant end position; 4. reference allele; 5. alternative allele. 

Example input file: 

Chr	Start	End	Ref	Alt

chr1	148542190	148542190	G	T

chr9	41098700	41098700	T	C


To get the edit distance of each variant to the gRNA sequence, Run:

./off_target_search/off_target.pl -in off_target_search/example.input.txt -out off_target_search/example.output.txt -genome inputs/fasta/hg38.fa -seq ACTCACGGTGGATCCCGCTG
