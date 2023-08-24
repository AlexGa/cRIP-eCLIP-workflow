def get_fastq_files(wildcards):
    fq_files = units.loc[wildcards.sample].loc[wildcards.unit]
    return([fq_files['fq1'], fq_files['fq2']])

rule clipTrim:
  input: get_fastq_files
  output: r1="results/clipTrim/{sample}-{unit}/{sample}_{unit}_R1.fq.gz",
          r2="results/clipTrim/{sample}-{unit}/{sample}_{unit}_R2.fq.gz"
  params: outDir = "results/clipTrim/{sample}-{unit}/{sample}_{unit}_R"
  conda: "../envs/clip.yaml"
  shell:
    "java -jar cRIP-eCLIP-workflow/java/DeBarCode.jar -fq1 {input[0]} -fq2 {input[1]} -p {params.outDir}"


rule cutadapt_pe_round1:
    input:
        r1="results/clipTrim/{sample}-{unit}/{sample}_{unit}_R1.fq.gz",
        r2="results/clipTrim/{sample}-{unit}/{sample}_{unit}_R2.fq.gz"
    output:
        fastq1="results/trimmed/round1/{sample}-{unit}/{sample}_{unit}_R1.fq.gz",
        fastq2="results/trimmed/round1/{sample}-{unit}/{sample}_{unit}_R2.fq.gz",
        qc="results/trimmed/round1/{sample}-{unit}/{sample}_{unit}.paired.qc.txt"
    log:
        "logs/cutadapt/round1/{sample}-{unit}/{sample}_{unit}.log",
    params:
        others=config["params"]["cutadapt-pe-round1"],
    threads: 8
    wrapper:
        "file:cRIP-eCLIP-workflow/wrapper/cutadapt_v0.30.0/pe"

rule cutadapt_pe_round2:
    input:
        r1="results/trimmed/round1/{sample}-{unit}/{sample}_{unit}_R1.fq.gz",
        r2="results/trimmed/round1/{sample}-{unit}/{sample}_{unit}_R2.fq.gz"
    output:
        fastq1="results/trimmed/round2/{sample}-{unit}/{sample}_{unit}_R1.fq.gz",
        fastq2="results/trimmed/round2/{sample}-{unit}/{sample}_{unit}_R2.fq.gz",
        qc="results/trimmed/round2/{sample}-{unit}/{sample}_{unit}.paired.qc.txt",
    log:
        "logs/cutadapt/round2/{sample}-{unit}/{sample}_{unit}.log",
    params:
        others=config["params"]["cutadapt-pe-round2"],
    threads: 8
    wrapper:
        "file:cRIP-eCLIP-workflow/wrapper/cutadapt_v0.30.0/pe"