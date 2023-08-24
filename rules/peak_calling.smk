def get_peak_file(wildcards):
    path = "results/peak_intersect/"+wildcards.organism+"/"
    return(path + sample_comp['unit2sample'][wildcards.unit] + "-IP-SM_"+wildcards.strand+"_peaks.bed")

def query_bed_files_for_DBA(wildcards):
    path = "results/confirmed_peaks/DBA/"+wildcards.organism+"/"+ sample_comp['bed-intersect']['query'][wildcards.sample] +"-"+wildcards.ip+"-"+wildcards.sm+"/"
    return(path + sample_comp['bed-intersect']['query'][wildcards.sample] + "_"+wildcards.ip+"_"+wildcards.sm+"_"+wildcards.strand+"_signif.bed")

def target_bed_files_for_DBA(wildcards):
    return(expand("results/confirmed_peaks/DBA/"+wildcards.organism+"/{unit_prefix}-"+wildcards.ip+"-"+wildcards.sm+"/{unit_prefix}_"+wildcards.ip+"_"+wildcards.sm+"_"+wildcards.strand+"_signif.bed", unit_prefix = sample_comp['bed-intersect']['target'][wildcards.sample]))

def query_bed_files(wildcards):
    path = "results/macs2/"+wildcards.organism+"/"+ sample_comp['bed-intersect']['query'][wildcards.sample] +"-"+wildcards.ip+"-"+wildcards.sm+"/"
    return(path + wildcards.strand+"_peaks.narrowPeak")

def target_bed_files(wildcards):
    return(expand("results/macs2/"+wildcards.organism+"/{unit_prefix}-"+wildcards.ip+"-"+wildcards.sm+"/" +wildcards.strand+"_peaks.narrowPeak", unit_prefix = sample_comp['bed-intersect']['target'][wildcards.sample]))

rule callpeak:
    input:
        treatment="results/dedup/{organism}/{ip}-{unit}/{ip}-{unit}_{strand}.bam",   
        control="results/dedup/{organism}/{sm}-{unit}/{sm}-{unit}_{strand}.bam"
    output:
        multiext("results/macs2/{organism}/{unit}-{ip}-{sm}/",
                 "{strand}_peaks.xls", 
                 "{strand}_peaks.narrowPeak",
                 "{strand}_summits.bed"
                 )
    log:
        "logs/macs2_stranded/{organism}/{unit}/{unit}-{ip}-vs-{sm}/{strand}/callpeak.log"
    params:
        lambda wildcards: "-f BAM " + config["params"]["macs-" + wildcards.organism]
    wrapper:
        "file:cRIP-eCLIP-workflow/wrapper/macs_peak"
        
rule peakCoverage:
    input: bam = "results/dedup/{organism}/{sample}-{unit}/{sample}-{unit}_{strand}.bam",
           bed = "results/macs2/{organism}/{unit}-IP-SM/{strand}_peaks.narrowPeak"
    output: "results/peaks_coverage/single/{organism}/{unit}-IP-SM/{sample}-{unit}_{strand}.bedgraph"
    conda: "../envs/clip.yaml"
    shell: "java -jar -Xmx20G cRIP-eCLIP-workflow/java/PeakCoverage2Bedgraph.jar -i {input.bam} -b {input.bed} -o {output[0]}"

rule peakCoverage_DBA:
    input: bam = "results/dedup/{organism}/{sample}-{unit}/{sample}-{unit}_{strand}.bam",
           bed = get_peak_file
    output: "results/peaks_coverage/DBA/{organism}/{unit}-IP-SM/{sample}-{unit}_{strand}.bedgraph"
    conda: "../envs/clip.yaml"
    shell: "java -jar -Xmx20G cRIP-eCLIP-workflow/java/PeakCoverage2Bedgraph.jar -i {input.bam} -b {input.bed} -o {output[0]}"

rule confirm_peaks_single:
    input:
        ip = "results/dedup/{organism}/{ip}-{unit}/{ip}-{unit}_{strand}.bam",   
        sm = "results/dedup/{organism}/{sm}-{unit}/{sm}-{unit}_{strand}.bam",
        ip_bed = "results/peaks_coverage/single/{organism}/{unit}-{ip}-{sm}/{ip}-{unit}_{strand}.bedgraph",   
        sm_bed = "results/peaks_coverage/single/{organism}/{unit}-{ip}-{sm}/{sm}-{unit}_{strand}.bedgraph",
        peaks = "results/macs2/{organism}/{unit}-{ip}-{sm}/{strand}_peaks.narrowPeak"
    output:
        raw = "results/confirmed_peaks/single/{organism}/{unit}-{ip}-{sm}/{unit}_{ip}_{sm}_{strand}_raw.tsv",
        filtered = "results/confirmed_peaks/single/{organism}/{unit}-{ip}-{sm}/{unit}_{ip}_{sm}_{strand}_signif.tsv",
        filtered_fc1 = "results/confirmed_peaks/single/{organism}/{unit}-{ip}-{sm}/{unit}_{ip}_{sm}_{strand}_signif_fc1.tsv",
        filtered_fc2 = "results/confirmed_peaks/single/{organism}/{unit}-{ip}-{sm}/{unit}_{ip}_{sm}_{strand}_signif_fc2.tsv",
        narrow = "results/confirmed_peaks/single/{organism}/{unit}-{ip}-{sm}/{unit}_{ip}_{sm}_{strand}_signif.bed",
        narrow_fc1 = "results/confirmed_peaks/single/{organism}/{unit}-{ip}-{sm}/{unit}_{ip}_{sm}_{strand}_signif_fc1.bed",
        narrow_fc2 = "results/confirmed_peaks/single/{organism}/{unit}-{ip}-{sm}/{unit}_{ip}_{sm}_{strand}_signif_fc2.bed"
    params: count_threshold = 20,
            use_thresholds = True
    wrapper:
        "file:cRIP-eCLIP-workflow/wrapper/confirm_peaks"

rule confirm_peaks_DBA:
    input:
        ip = "results/dedup/{organism}/{ip}-{unit}/{ip}-{unit}_{strand}.bam",   
        sm = "results/dedup/{organism}/{sm}-{unit}/{sm}-{unit}_{strand}.bam",
        ip_bed = "results/peaks_coverage/DBA/{organism}/{unit}-{ip}-{sm}/{ip}-{unit}_{strand}.bedgraph",   
        sm_bed = "results/peaks_coverage/DBA/{organism}/{unit}-{ip}-{sm}/{sm}-{unit}_{strand}.bedgraph",
        peaks = get_peak_file
    output:
        raw = "results/confirmed_peaks/DBA/{organism}/{unit}-{ip}-{sm}/{unit}_{ip}_{sm}_{strand}_raw.tsv",
        filtered = "results/confirmed_peaks/DBA/{organism}/{unit}-{ip}-{sm}/{unit}_{ip}_{sm}_{strand}_signif.tsv",
        filtered_fc1 = "results/confirmed_peaks/DBA/{organism}/{unit}-{ip}-{sm}/{unit}_{ip}_{sm}_{strand}_signif_fc1.tsv",
        filtered_fc2 = "results/confirmed_peaks/DBA/{organism}/{unit}-{ip}-{sm}/{unit}_{ip}_{sm}_{strand}_signif_fc2.tsv",
        narrow = "results/confirmed_peaks/DBA/{organism}/{unit}-{ip}-{sm}/{unit}_{ip}_{sm}_{strand}_signif.bed",
        narrow_fc1 = "results/confirmed_peaks/DBA/{organism}/{unit}-{ip}-{sm}/{unit}_{ip}_{sm}_{strand}_signif_fc1.bed",
        narrow_fc2 = "results/confirmed_peaks/DBA/{organism}/{unit}-{ip}-{sm}/{unit}_{ip}_{sm}_{strand}_signif_fc2.bed"
    params: count_threshold = 20,
            use_thresholds = True
    wrapper:
        "file:cRIP-eCLIP-workflow/wrapper/confirm_peaks"

rule peak_intersect:
    input: query = query_bed_files,
           target = target_bed_files
    output: "results/peak_intersect/{organism}/{sample}-{ip}-{sm}_{strand}_peaks.bed"
    threads: 1
    log: "logs/peak_intersect/{organism}/{sample}-{ip}-{sm}_{strand}_peaks.log"
    wrapper: "file:cRIP-eCLIP-workflow/wrapper/peak_intersect"

rule peak_merge_DBA:
    input: query = query_bed_files_for_DBA,
           target = target_bed_files_for_DBA
    output: "results/peak_merge_DBA/{organism}/{sample}-{ip}-{sm}_{strand}_peaks.bed"
    threads: 1
    log: "logs/peak_merge_DBA/{organism}/{sample}-{ip}-{sm}_{strand}_peaks.log"
    script:
        "../R/merge_narrow_files.R"

rule DBA:
    input: peaks = "results/peak_merge_DBA/{organism}/{unit}-IP-SM_{strand}_peaks.bed",
           ctrl = lambda wildcards: "results/dedup/{organism}/IP-" + sample_comp['DBA'][wildcards.unit]['ctrl'] + "/IP-" + sample_comp['DBA'][wildcards.unit]['ctrl'] + "_{strand}.bam",
           ko = lambda wildcards: "results/dedup/{organism}/IP-" + sample_comp['DBA'][wildcards.unit]['ko'] + "/IP-" + sample_comp['DBA'][wildcards.unit]['ko'] + "_{strand}.bam"
    output: raw = "results/DBA/{organism}/{unit}/{unit}_{strand}_raw.tsv",
            filtered = "results/DBA/{organism}/{unit}/{unit}_{strand}_signif.tsv",
            bed = "results/DBA/{organism}/{unit}/{unit}_{strand}_signif.bed"
    threads: 1
    log: "logs/DBA/{organism}/{unit}_{strand}_peaks.log"
    wrapper: "file:cRIP-eCLIP-workflow/wrapper/DBA"

rule DBA_strand:
    input: peaks = "results/confirmed_peaks/single/{organism}/{unit}-IP-SM/{unit}_IP_SM_{strand}_signif.bed",
           fwd = "results/dedup/{organism}/IP-{unit}/IP-{unit}_fwd.bam",
           rev = "results/dedup/{organism}/IP-{unit}/IP-{unit}_rev.bam"
    output: raw = "results/DBA_strand/{organism}/{unit}/{unit}_{strand}_peaks_raw.tsv",
            filtered = "results/DBA_strand/{organism}/{unit}/{unit}_{strand}_peaks_signif.tsv",
    threads: 1
    params: strand = "{strand}",
            count_threshold = 10
    log: "logs/DBA_strand/{organism}/{unit}_{strand}_peaks.log"
    wrapper: "file:cRIP-eCLIP-workflow/wrapper/DBA_strand"

rule bw_tracks_norm_peaks:
    input: ip_bam = expand("results/dedup/{{organism}}/IP-{{unit}}/IP-{{unit}}_{strand}.bam", strand = ['fwd', 'rev']),
           sm_bam = expand("results/dedup/{{organism}}/SM-{{unit}}/SM-{{unit}}_{strand}.bam", strand = ['fwd', 'rev']),
           ip_bw = expand("results/bigwig/{{organism}}/IP-{{unit}}/IP-{{unit}}_{strand}.bw", strand = ['fwd', 'rev']),
           sm_bw = expand("results/bigwig/{{organism}}/SM-{{unit}}/SM-{{unit}}_{strand}.bw", strand = ['fwd', 'rev']),
    output: raw = expand("results/normed_tracks_peaks/{{organism}}/Raw/{{unit}}/{{unit}}_RAW_IP_{strand}.bw", strand = ['fwd', 'rev']),
            rpm = expand("results/normed_tracks_peaks/{{organism}}/RPM/{{unit}}/{{unit}}_RPM_IP_{strand}.bw", strand = ['fwd', 'rev']),
            subtract = expand("results/normed_tracks_peaks/{{organism}}/Subtract/{{unit}}/{{unit}}_SUBTRACT_{strand}.bw", strand = ['fwd', 'rev']),
            kl = expand("results/normed_tracks_peaks/{{organism}}/KL/{{unit}}/{{unit}}_KL_{strand}.bw", strand = ['fwd', 'rev']),
            lfc = expand("results/normed_tracks_peaks/{{organism}}/LFC/{{unit}}/{{unit}}_LFC_{strand}.bw", strand = ['fwd', 'rev'])
    threads: 1
    params: raw = "results/normed_tracks_peaks/{organism}/Raw/{unit}/{unit}_RAW",
            rpm = "results/normed_tracks_peaks/{organism}/RPM/{unit}/{unit}_RPM",
            subtract = "results/normed_tracks_peaks/{organism}/Subtract/{unit}/{unit}_SUBTRACT",
            kl = "results/normed_tracks_peaks/{organism}/KL/{unit}/{unit}_KL",
            lfc = "results/normed_tracks_peaks/{organism}/LFC/{unit}/{unit}_LFC"
    log: "logs/normed_tracks_peaks/{organism}/{unit}_xl_sites.log"
    wrapper: "file:cRIP-eCLIP-workflow/wrapper/norm_tracks"