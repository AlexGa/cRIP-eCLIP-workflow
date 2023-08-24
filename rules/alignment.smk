ruleorder: star_index_gtf > star_index_fasta
ruleorder: extract_subgenome > mark_duplicates

def get_fq(wildcards):
    u = units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
    path = "results/trimmed/round2/" + wildcards.sample + "/" + wildcards.unit + "/"
    return([path + u["fq1"], path + u["fq2"]])

def get_mem_star(wildcards):
    host_virus = config['ref']['host'] + "_" + config['ref']['virus']
    return config["params"]["star-map-" + host_virus + "-ram"] / 1e6

def is_in_sample_comp(wildcards):
    return wildcards.unit in sample_comp['unit2sample']

rule star_index_gtf:
    input:
        fasta = lambda wildcards: config["ref"]["fasta"][wildcards.organism],
        gtf = lambda wildcards: config["ref"]["gtf"][wildcards.organism]
    output:
        files = expand("references/STAR_INDEX/{{organism}}/{suf}", 
                suf = ["chrNameLength.txt", "Log.out", "chrStart.txt"]),
    log:
        "logs/STAR/{organism}/star_index.log"
    threads: 20
    params:
        ram = lambda wildcards: config["params"]["star-map-" + wildcards.organism + "-ram"],
        star_usr = lambda wildcards: config["params"]["star-index-" + wildcards.organism],
        dir = directory("references/STAR_INDEX/{organism}")
    wrapper:
        "file:cRIP-eCLIP-workflow/wrapper/star_index"

rule star_index_fasta:
    input:
        fasta = lambda wildcards: config["ref"]["fasta"][wildcards.organism]
    output:
        files = expand("references/STAR_INDEX/{{organism}}/{suf}", 
                suf = [ "chrNameLength.txt", "Log.out", "chrStart.txt"]),
    log:
        "logs/STAR/{organism}/{organism}_index.log"
    threads: 20
    params:
        ram = lambda wildcards: config["params"]["star-map-" + wildcards.organism + "-ram"],
        star_usr = lambda wildcards: config["params"]["star-index-" + wildcards.organism],
        dir = directory("references/STAR_INDEX/{organism}")
    wrapper:
        "file:cRIP-eCLIP-workflow/wrapper/star_index"

rule star_pe_multi:
    input:
        fq1 = "results/trimmed/round2/{sample}-{unit}/{sample}_{unit}_R1.fq.gz",
        fq2 = "results/trimmed/round2/{sample}-{unit}/{sample}_{unit}_R2.fq.gz",
        idx = expand("references/STAR_INDEX/" + config['ref']['host'] + "_" + config['ref']['virus'] + "/{suf}", suf = ["chrNameLength.txt", "Log.out", "chrStart.txt"])
    output:
        # see STAR manual for additional output files
        bam = "results/star/" + config['ref']['host'] + "_" + config['ref']['virus'] + "/{sample}-{unit}/{sample}-{unit}_aligned.bam",
        log = "results/star/" + config['ref']['host'] + "_" + config['ref']['virus'] + "/{sample}-{unit}/Log.out",
        log_final = "results/star/" + config['ref']['host'] + "_" + config['ref']['virus'] + "/{sample}-{unit}/Log.final.out",
        unmapped = expand("results/star/" + config['ref']['host'] + "_" + config['ref']['virus'] + "/{{sample}}-{{unit}}/Unmapped.out.mate{read}", read = ["1", "2"])
    log:
        "logs/STAR_align/" + config['ref']['host'] + "_" + config['ref']['virus'] + "/{sample}-{unit}.log"
    resources:
        mem_mb = get_mem_star
    threads: 16
    params:
        # path to STAR reference genome index
        idx="references/STAR_INDEX/" + config['ref']['host'] + "_" + config['ref']['virus'] + "/",
        # optional parameters
        extra=lambda wildcards: config["params"]["star-map-" + config['ref']['host'] + "_" + config['ref']['virus']],
        tmpdir = "results/star/" + config['ref']['host'] + "_" + config['ref']['virus'] + "/{sample}-{unit}"
    wrapper:
        "file:cRIP-eCLIP-workflow/wrapper/star_align"

rule reanno:
    input: "results/star/{organism}/{sample}-{unit}/{sample}-{unit}_aligned.bam"
    output: "results/star/{organism}/{sample}-{unit}/{sample}-{unit}_aligned_reanno.bam"
    conda: "../envs/clip.yaml"   
    resources:
        mem_mb=100000
    shell:
        "java -jar -Xmx{resources.mem_mb}m cRIP-eCLIP-workflow/java/ReannoBarCode.jar -i {input[0]} -o {output[0]}"

rule mark_duplicates:
    input:
        bams="results/star/" + config['ref']['host'] + "_" + config['ref']['virus'] + "/{sample}-{unit}/{sample}-{unit}_aligned_reanno.bam"
    output:
        bam="results/dedup/" + config['ref']['host'] + "_" + config['ref']['virus'] + "/{sample}-{unit}/{sample}-{unit}_aligned.bam",
        metrics="results/dedup/" + config['ref']['host'] + "_" + config['ref']['virus'] + "/{sample}-{unit}/{sample}-{unit}.metrics.txt"
    log:
        "logs/dedup/" + config['ref']['host'] + "_" + config['ref']['virus'] + "/{sample}-{unit}_marked_duplicates.log",
    params:
        extra="--REMOVE_DUPLICATES true --BARCODE_TAG XB"
    resources:
        mem_mb=100000
    wrapper:
        "file:cRIP-eCLIP-workflow/wrapper/picard_dedup"

rule bam_index:
    input: "results/{sample}.bam"
    output: "results/{sample}.bai"
    log:
        "logs/index/{sample}.log",
    wrapper:
        "file:cRIP-eCLIP-workflow/wrapper/samtools_index"

rule split_by_strand:
    input: "results/dedup/{organism}/{sample}-{unit}/{sample}-{unit}_aligned.bam"
    output: fwd = "results/dedup/{organism}/{sample}-{unit}/{sample}-{unit}_fwd.bam", 
            rev = "results/dedup/{organism}/{sample}-{unit}/{sample}-{unit}_rev.bam",
            paired = "results/dedup/{organism}/{sample}-{unit}/{sample}-{unit}_paired.bam"
    log:
        "log/samtools/split_by_strand/{organism}/{sample}-{unit}.log"
    params:
        extra = "",  # optional params string
        tmpdir = temp("results/dedup/{organism}/{sample}-{unit}/tmp")
    threads: 1
    wrapper:
        "file:cRIP-eCLIP-workflow/wrapper/samtools_split_strand"

rule coverage_rpkm_bw:
    input:
        bam="results/dedup/{organism}/{sample}-{unit}/{sample}-{unit}_aligned.bam",
        bai="results/dedup/{organism}/{sample}-{unit}/{sample}-{unit}_aligned.bai"
    output:
        bigwig="results/bigwig/{organism}/{sample}-{unit}/{sample}-{unit}-{strand}.bw",
    log:
        "logs/bigwig/{organism}/{sample}-{unit}-{strand}.log",
    params:
        extra="--binSize 1 --normalizeUsing RPKM",
        strand = "{strand}"
    wrapper:
         "file:cRIP-eCLIP-workflow/wrapper/deeptools_coverage"

rule coverage_counts_bw:
    input:
        bam="results/dedup/{organism}/{sample}-{unit}/{sample}-{unit}_aligned.bam",
        bai="results/dedup/{organism}/{sample}-{unit}/{sample}-{unit}_aligned.bai"
    output:
        bigwig="results/bigwig/{organism}/{sample}-{unit}/{sample}-{unit}.bw",
    log:
        "logs/bigwig/{organism}/{sample}-{unit}.log",
    params:
        extra = "--binSize 1",
        strand = "" 
    wrapper:
         "file:cRIP-eCLIP-workflow/wrapper/deeptools_coverage"

rule coverage_counts_strand:
    input:
        bam="results/dedup/{organism}/{sample}-{unit}/{sample}-{unit}_{strand}.bam",
        bai="results/dedup/{organism}/{sample}-{unit}/{sample}-{unit}_{strand}.bai"
    output:
        bigwig="results/bigwig/{organism}/{sample}-{unit}/{sample}-{unit}_{strand}.bw",
    log:
        "logs/bigwig/{organism}/{sample}-{unit}_{strand}.log",
    params:
        extra = "--binSize 1",
        strand = "" 
    wrapper:
         "file:cRIP-eCLIP-workflow/wrapper/deeptools_coverage"

rule coverage_counts_bw_reverse:
    input:
        bam="results/dedup/{organism}/{sample}-{unit}/{sample}-{unit}_aligned_{strand}.bam",
        bai="results/dedup/{organism}/{sample}-{unit}/{sample}-{unit}_aligned_{strand}.bai"
    output:
        bigwig="results/bigwig/{organism}/reverse/{sample}-{unit}/{sample}-{unit}_strand.bw",
    log:
        "logs/bigwig/{organism}/reverse/{sample}-{unit}.log",
    params:
        extra = "--binSize 1 --strand {strand}",
        strand = "" 
    wrapper:
         "file:cRIP-eCLIP-workflow/wrapper/deeptools_coverage"

rule coverage_bw:
    input:
        bam="results/R2/{organism}/{sample}-{unit}/{sample}-{unit}_{strand}.bam",
        bai="results/R2/{organism}/{sample}-{unit}/{sample}-{unit}_{strand}.bai"
    output:
        bigwig="results/bigwig/{organism}/{sample}-{unit}/{sample}-{unit}_{strand}.bw",
    log:
        "logs/bigwig/{organism}/{sample}-{unit}_{strand}.log",
    params:
        extra="--binSize 1"
    wrapper:
         "file:cRIP-eCLIP-workflow/wrapper/deeptools_coverage"

rule fastq_sort:
    input:
        "{star_dir}/{sample}-{unit}/Unmapped.out.mate{read}"
    output:
        temp("{star_dir}/{sample}-{unit}/{sample}_{unit}_R{read}.fq")
    conda: "../envs/fastq_tools.yaml"
    shell: "fastq-sort --id {input} > {output}"


rule pack_reads:
    input:
        "{star_dir}/{sample}-{unit}/{sample}_{unit}_R{read}.fq"
    output:
        "{star_dir}/{sample}-{unit}/{sample}_{unit}_R{read}.fq.gz"
    threads: 1
    conda: "../envs/fastq_tools.yaml"
    shell: "gzip {input}"

rule extract_subgenome:
    input: "results/dedup/" + config['ref']['host'] + "_{virus}/{sample}-{unit}/{sample}-{unit}_aligned.bam"
    output: bam_sorted = temp("results/dedup/" + config['ref']['host'] + "_{virus}/{sample}-{unit}/{sample}-{unit}_aligned_sorted.bam"),
            bam = "results/dedup/{virus}/{sample}-{unit}/{sample}-{unit}_aligned.bam"
    log:
        "logs/extract_chromosomes_{virus}/{sample}/{unit}.log",
    params:
        region= lambda wildcards: config['subgenome'][wildcards.virus] #"MN908947.3",  # optional params string
    threads: 1
    wrapper:
        "file:cRIP-eCLIP-workflow/wrapper/samtools_view"

rule ratio_viral_human:
    input: all = "results/dedup/" + config['ref']['host'] + "_{suffix}/{sample}-{unit}/{sample}-{unit}_{suffix}.bam",
           viral = "results/dedup/{virus}/{sample}-{unit}/{sample}-{unit}_{suffix}.bam",
    output: txt = "results/host_virus_ratio_{virus}/{sample}-{unit}.txt"
    log: "logs/ratio/ratio_host_virus_{virus}/{sample}-{unit}.log"
    threads: 4
    params:
        read_type = "paired" # single for single end read
    wrapper:
        "file:cRIP-eCLIP-workflow/wrapper/ratio_virus_human"

rule extract_R2:
    input: bam = "results/dedup/{organism}/{sample}-{unit}/{sample}-{unit}_{suffix}.bam"
    output: bam = "results/R2/{organism}/{sample}-{unit}/{sample}-{unit}_{suffix}_R2.bam",
            bai = "results/R2/{organism}/{sample}-{unit}/{sample}-{unit}_{suffix}_R2.bam.bai"
    threads: 8
    log: "logs/extract_R2/{organism}/{sample}-{unit}_{suffix}.log"
    wrapper:
        "file:cRIP-eCLIP-workflow/wrapper/samtools_extract_R2"



