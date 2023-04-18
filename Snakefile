import os
import pandas as pd

configfile: "inputs/input.yaml"
units = pd.read_table(config["samplesheet"], dtype=str).set_index(["sample"], drop=False)
glen = pd.read_table(config["genome"] + ".bed", dtype=str, header=0, names=['chr','len']).set_index(["chr"], drop=False)
SAMPLES = units.index.get_level_values('sample').unique().tolist()
CHRS = glen.index.get_level_values('chr').unique().tolist()
TMPDIR = config["tmp_dir"]
RESULTDIR = config["result_dir"]

def get_read1(wildcards):
    if len(units.loc[wildcards.sample,"read1fq"].split(","))==1:
        return units.loc[wildcards.sample,"read1fq"]
    else:
        return units.loc[wildcards.sample,"read1fq"].split(",")[int(wildcards.subfile)]

def get_read2(wildcards):
    if len(units.loc[wildcards.sample,"read2fq"].split(","))==1:
        return units.loc[wildcards.sample,"read2fq"]
    else:
        return units.loc[wildcards.sample,"read2fq"].split(",")[int(wildcards.subfile)]

def get_bam(wildcards):
    if len(units.loc[wildcards.sample,"read1fq"].split(","))==1:
        return TMPDIR + "/" + wildcards.sample + ".part0.sort.bam"
    else:
        return expand(TMPDIR + "/" + wildcards.sample + ".part{subfile}.sort.bam",subfile = list(range(0,len(units.loc[wildcards.sample,"read1fq"].split(",")))))

def get_gvcfs_input(wildcards):
    return expand(TMPDIR + "/" + "{sample}.HC.g.vcf", sample = wildcards.samples.split("."))

def get_bams_input(wildcards):
    return expand(RESULTDIR + "/bam/" + "{sample}.sort.dedup.recal.bam", sample = wildcards.samples.split("."))

def get_bai_input(wildcards):
    return expand(RESULTDIR + "/bam/" + "{sample}.sort.dedup.recal.bam.bai", sample = wildcards.samples.split("."))

def get_chr_vcfgzs(wildcards):
    return expand(TMPDIR + "/" + wildcards.sample + ".SC.{chr}.fil.vcf.gz", chr = CHRS)

def get_chr_len(wildcards):
    return glen.loc[wildcards.chr,"len"]

def rchop(s, suffix):
    if suffix and s.endswith(suffix):
        return s[:-len(suffix)]
    return s

onstart:
    print("##### Creating profile pipeline #####\n") 
    print("\t Creating jobs output subfolders...\n")
    shell("mkdir -p jobs/all")
    shell("mkdir -p jobs/indexgenome")
    shell("mkdir -p jobs/cutadapt")
    shell("mkdir -p jobs/bwamem")
    shell("mkdir -p jobs/sortbam")
    shell("mkdir -p jobs/markduplicate")
    shell("mkdir -p jobs/indexbam_result")
    shell("mkdir -p jobs/indexbam_temp")
    shell("mkdir -p jobs/indexfasta")
    shell("mkdir -p jobs/BQSR")
    shell("mkdir -p jobs/HaplotypeCaller")
    shell("mkdir -p jobs/JointGenotyping")
    shell("mkdir -p jobs/VQSR")
    shell("mkdir -p jobs/Mutect2")
    shell("mkdir -p jobs/Lofreq")
    shell("mkdir -p jobs/Scalpel")
    shell("mkdir -p jobs/Bgzip")
    shell("mkdir -p jobs/MergeVcfs")

rule all:
    input:
        expand(RESULTDIR + "/vcf/{targetfile}", targetfile=config["targetfiles"].split(","))

rule indexgenome:
    input:
        config["genome"]
    output:
        config["genome"] + ".bwt"
    threads: 1
    shell:
        "bwa index {input}"

rule cutadapt:
    input:
        read1 = get_read1,
        read2 = get_read2
    output:
        read1 = temp(TMPDIR + "/" + "{sample}.part{subfile,[0-9]}.r1.fastq"),
        read2 = temp(TMPDIR + "/" + "{sample}.part{subfile,[0-9]}.r2.fastq")
    threads: 8
    shell:
        "cutadapt -g AGATGTGTATAAGAGACAG -G AGATGTGTATAAGAGACAG "
        "-a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT "
        "-g TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -G TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG "
        "-a CTGTCTCTTATACACATCTGACGCTGCCGACGA -A CTGTCTCTTATACACATCTGACGCTGCCGACGA "
        "-g GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -G GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG "
        "-a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -A CTGTCTCTTATACACATCTCCGAGCCCACGAGAC "
        "-j {threads} -o {output.read1} -p {output.read2} -q 20 --trim-n -m 20 --max-n=0 {input.read1} {input.read2}"

rule bwamem:
    input:
        read1 = TMPDIR + "/" + "{sample}.part{subfile}.r1.fastq",
        read2 = TMPDIR + "/" + "{sample}.part{subfile}.r2.fastq",
        genome = config["genome"],
        index = config["genome"] + ".bwt"
    output:
        temp(TMPDIR + "/" + "{sample}.part{subfile,[0-9]}.bam")
    threads: 8
    shell:
        "bwa mem -t {threads} {input.genome} -R \"@RG\\tID:{wildcards.sample}.{wildcards.subfile}\\tSM:{wildcards.sample}\\tPL:illumina\\tLB:{wildcards.sample}\\tPU:{wildcards.sample}.{wildcards.subfile}\" {input.read1} {input.read2} | samtools view -bS > {output}"

rule sortbam:
    input:
        TMPDIR + "/" + "{sample}.part{subfile}.bam",
    output:
        temp(TMPDIR + "/" + "{sample}.part{subfile,[0-9]}.sort.bam")
    threads: 1
    params:
        TMPDIR
    shell:
        "gatk SortSam -I {input} -O {output} -SO coordinate --VALIDATION_STRINGENCY SILENT --TMPDIR {params}"

rule markduplicate:
    input:
        get_bam
    output:
        temp(TMPDIR + "/" + "{sample}.sort.dedup.bam")
    threads: 1
    params:
        input = lambda wildcards: " ".join(expand(" -I {inputs} ",inputs = get_bam(wildcards))),
        tmpdir = TMPDIR
    shell:
        "gatk MarkDuplicates {params.input} -O {output} -M {params.tmpdir}/{wildcards.sample}.metrics --VALIDATION_STRINGENCY SILENT --TMP_DIR {params.tmpdir}"

rule indexbam_result:
    input:
        RESULTDIR + "/bam/" + "{sample}.bam"
    output:
        RESULTDIR + "/bam/" + "{sample}.bam.bai"
    threads: 1
    shell:
        "samtools index {input}"

rule indexbam_temp:
    input:
        TMPDIR + "/" + "{sample}.bam"
    output:
        temp(TMPDIR + "/" + "{sample}.bam.bai")
    threads: 1
    shell:
        "samtools index {input}"

rule indexfasta:
    input:
        config["genome"]
    output:
        config["genome"] + ".fai",
        rchop(rchop(rchop(config["genome"],".gz"),".fasta"),".fa") + ".dict"
    threads: 1
    shell:
        "samtools faidx {input} && "
        "gatk CreateSequenceDictionary -R {input}"

rule BQSR:
    input:
        bam = TMPDIR + "/" + "{sample}.sort.dedup.bam",
        bai = TMPDIR + "/" + "{sample}.sort.dedup.bam.bai",
        genome = config["genome"],
        fai = config["genome"] + ".fai",
        dict = rchop(rchop(rchop(config["genome"],".gz"),".fasta"),".fa") +  ".dict"
    output:
        table = temp(TMPDIR + "/" + "{sample}.recal.table"),
        bam = RESULTDIR + "/bam/" + "{sample}.sort.dedup.recal.bam",
        bai = RESULTDIR + "/bam/" + "{sample}.sort.dedup.recal.bai"
    params:
        " ".join(expand(" --known-sites {knownvcfs} ",knownvcfs = config['knownvcfsBQSR'].split(";")))
    threads: 1
    shell:
        "gatk BaseRecalibrator -R {input.genome} -I {input.bam} {params} -O {output.table} --use-original-qualities && "
        "gatk ApplyBQSR --create-output-bam-md5 --add-output-sam-program-record -R {input.genome} -I {input.bam} -O {output.bam} -bqsr {output.table} --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-qual 30"

rule HaplotypeCaller:
    input:
        bam = RESULTDIR + "/bam/" + "{sample}.sort.dedup.recal.bam",
        bai = RESULTDIR + "/bam/" + "{sample}.sort.dedup.recal.bam.bai",
        genome = config["genome"]
    output:
        temp(TMPDIR + "/" + "{sample}.HC.g.vcf")
    threads: 1
    shell:
        "gatk HaplotypeCaller -R {input.genome} -I {input.bam} -O {output} -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 -ERC GVCF"

rule JointGenotyping:
    input:
        gvcf = get_gvcfs_input,
        genome = config["genome"]
    output:
        gvcf = temp(TMPDIR + "/" + "{samples}.HC.joint.g.vcf"),
        gvcfidx = temp(TMPDIR + "/" + "{samples}.HC.joint.g.vcf.idx"),
        vcf = temp(RESULTDIR + "/vcf/" + "{samples}.HC.raw.vcf"),
        vcfidx = temp(RESULTDIR + "/vcf/" + "{samples}.HC.raw.vcf.idx")
    threads: 1
    params:
        lambda wildcards: " ".join(expand("--variant " + TMPDIR + "/" + "{sample}.HC.g.vcf", sample=wildcards.samples.split(".")))
    shell:
        "gatk CombineGVCFs -R {input.genome} {params} -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation -O {output.gvcf} && "
        "gatk GenotypeGVCFs -R {input.genome} -V {output.gvcf} -O {output.vcf} -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation"

rule VQSR:
    input:
        vcf = "{sample}.HC.raw.vcf",
        vcfidx = "{sample}.HC.raw.vcf.idx",
        genome = config["genome"],
        dict = rchop(rchop(rchop(config["genome"],".gz"),".fasta"),".fa") + ".dict"
    output:
        tmpvcf = temp("{sample}.HC.indelrecal.vcf"),
        tmpvcfidx = temp("{sample}.HC.indelrecal.vcf.idx"),
        vcf = "{sample}.HC.fil.vcf",
        idx = temp("{sample}.HC.fil.vcf.idx"),
        tmp = temp("{sample}.tmp.sites_only.vcf"),
        tmpidx = temp("{sample}.tmp.sites_only.vcf.idx"),
        indelrecal = temp("{sample}.indel.recal"),
        indelrecalidx = temp("{sample}.indel.recal.idx"),
        indeltranche = temp("{sample}.indel.tranche"),
        snprecal = temp("{sample}.snp.recal"),
        snprecalidx = temp("{sample}.snp.recal.idx"),
        snptranche = temp("{sample}.snp.tranche")
#        RESULTDIR + "/vcf/" + "{samples}.HC.fil.vcf"
    threads: 1
    params:
        indel = " ".join(config["paramsVQSRindel"].split(";")),
        snp = " ".join(config["paramsVQSRsnp"].split(";"))
    shell:
        "gatk MakeSitesOnlyVcf -I {input.vcf} -O {output.tmp} && "
        "gatk VariantRecalibrator -V {output.tmp} -O {output.indelrecal} --tranches-file {output.indeltranche} --trust-all-polymorphic "
        "-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 "
        "-an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP "
        "-AS -mode INDEL --max-gaussians 4 {params.indel} && "
        "gatk VariantRecalibrator -V {output.tmp} -O {output.snprecal} --tranches-file {output.snptranche} --trust-all-polymorphic "
        "-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 "
        "-an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP "
        "-AS -mode SNP --max-gaussians 6 {params.snp} && "
        "gatk ApplyVQSR -V {input.vcf} -O {output.tmpvcf} --recal-file {output.indelrecal} --tranches-file {output.indeltranche} --truth-sensitivity-filter-level 99.0 "
        "--create-output-variant-index true -mode INDEL"
        "gatk ApplyVQSR -V {output.tmpvcf} -O {output.vcf} --recal-file {output.snprecal} --tranches-file {output.snptranche} --truth-sensitivity-filter-level 99.7 "
        "--create-output-variant-index true -mode SNP -AS"

rule Mutect2:
    input:
        bam = get_bams_input, 
        bai =  get_bai_input, 
        genome = config["genome"]
    output:
        raw = temp(TMPDIR + "/{samples}.MT.raw.vcf"),
        fil = RESULTDIR + "/vcf/" + "{samples}.MT.fil.vcf",
        rawidx = temp(TMPDIR + "/vcf/" + "{samples}.MT.raw.vcf.idx"),
        filidx = temp(RESULTDIR + "/vcf/" + "{samples}.MT.fil.vcf.idx"),
        filstat = temp(RESULTDIR + "/vcf/" + "{samples}.MT.fil.vcf.filteringStats.tsv"),
        rawstat = temp(TMPDIR + "/vcf/" + "{samples}.MT.raw.vcf.stats"),
    params:
        lambda wildcards: " ".join(expand("-I" + RESULTDIR + "/bam/" + "{sample}.sort.dedup.recal.bam", sample=wildcards.samples.split(".")))
    threads: 8
    shell:
        "gatk Mutect2 -R {input.genome} {params} -O {output.raw} && "
        "gatk FilterMutectCalls -V {output.raw} -R {input.genome} -O {output.fil}"

rule Lofreq:
    input:
        bam = RESULTDIR + "/bam/" + "{sample}.sort.dedup.recal.bam",
        bai = RESULTDIR + "/bam/" + "{sample}.sort.dedup.recal.bam.bai",
        genome = config["genome"]
    output:
        RESULTDIR + "/vcf/" + "{sample}.LF.fil.vcf"
    threads: 8
    params:
        " ".join(expand(" -S {knownvcfs} ",knownvcfs = config['knownvcfsBQSR'].split(";")))
    shell:
        "lofreq call-parallel --pp-threads {threads} -f {input.genome} -o {output} -s {params} {input.bam}"

rule Scalpel:
    input:
        bam = RESULTDIR + "/bam/" + "{sample}.sort.dedup.recal.bam",
        bai = RESULTDIR + "/bam/" + "{sample}.sort.dedup.recal.bam.bai",
        genome = config["genome"]
    output:
        vcf = temp(TMPDIR + "/" + "{sample}.SC.{chr}.fil.vcf"),
        dir = temp(directory(TMPDIR + "/" + "{sample}.{chr}"))
    params:
        chrlen = get_chr_len,
        window = 1000
    threads: 8
    shell:
        "scalpel-discovery --single --bam {input.bam} --ref {input.genome} --bed {wildcards.chr}:1-{params.chrlen} --window {params.window} --numprocs {threads} --dir {output.dir} && "
        "scalpel-export --single --db {output.dir}/variants.db --bed {wildcards.chr}:1-{params.chrlen} --ref {input.genome} > {output.vcf}"


rule Bgzip:
    input:
        "{sample}.vcf"
    output:
        gz = temp("{sample}.vcf.gz"),
        tbi = temp("{sample}.vcf.gz.tbi")
    shell:
        "bgzip {input} && tabix {output.gz}"

rule MergeVcfs:
    input:
        get_chr_vcfgzs
    output:
        RESULTDIR + "/vcf/" + "{sample}.SC.fil.vcf"
    shell:
        "vcf-concat {input} > {output}"
