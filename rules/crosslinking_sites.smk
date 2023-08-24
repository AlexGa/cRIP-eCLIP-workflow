rule pileup_xlsites:
    input: "results/dedup/{organism}/{sample}-{unit}/{sample}-{unit}_{strand}.bam"
    output: "results/XL_sites_bedgraph/raw/{organism}/{sample}-{unit}/{sample}-{unit}_{strand}.bedgraph"
    conda: "../envs/clip.yaml"
    shell: "java -jar -Xmx60G cRIP-eCLIP-workflow/java/Pileup_XLsites.jar -i {input[0]} -o {output[0]}"

rule pileup_xlsites_single:
    input: bam = "results/dedup/{organism}/{sample}-{unit}/{sample}-{unit}_{strand}.bam",
           bed = "results/confirmed_peaks/single/{organism}/{unit}-IP-SM/{unit}_IP_SM_{strand}_signif.bed",
    output: "results/XL_sites_bedgraph/single/{organism}/{sample}-{unit}/{sample}-{unit}_{strand}.bedgraph"
    conda: "../envs/clip.yaml"
    shell: "java -jar -Xmx60G cRIP-eCLIP-workflow/java/Pileup_XLsites.jar -i {input.bam} -b {input.bed} -o {output[0]}"

rule pileup_xlsites_DBA:
    input: bam = "results/dedup/{organism}/{sample}-{unit}/{sample}-{unit}_{strand}.bam",
           bed = lambda wildcards: "results/DBA/{organism}/" + sample_comp['unit2sample'][wildcards.unit] + "/" + sample_comp['unit2sample'][wildcards.unit] + "_{strand}_signif.bed"
    output: "results/XL_sites_bedgraph/DBA/{organism}/{sample}-{unit}/{sample}-{unit}_{strand}.bedgraph"
    conda: "../envs/clip.yaml"
    shell: "java -jar -Xmx60G cRIP-eCLIP-workflow/java/Pileup_XLsites.jar -i {input.bam} -b {input.bed} -o {output[0]}"

rule XL_sites:
    input: ip = "results/dedup/{organism}/IP-{unit}/IP-{unit}_{strand}.bam",   
           sm = "results/dedup/{organism}/SM-{unit}/SM-{unit}_{strand}.bam",
           ip_bed = "results/XL_sites_bedgraph/raw/{organism}/IP-{unit}/IP-{unit}_{strand}.bedgraph",   
           sm_bed = "results/XL_sites_bedgraph/raw/{organism}/SM-{unit}/SM-{unit}_{strand}.bedgraph"
    output: tsv_raw = "results/XL_sites/all/{organism}/{unit}-IP-SM/{unit}_IP_SM_{strand}_raw.tsv",
            tsv_filtered = "results/XL_sites/all/{organism}/{unit}-IP-SM/{unit}_IP_SM_{strand}_signif.tsv",
            bed = "results/XL_sites/all/{organism}/{unit}-IP-SM/{unit}_IP_SM_{strand}_signif.bed",
            bw_ip = "results/XL_sites/all/{organism}/{unit}-IP-SM/IP-{unit}_{strand}_raw.bw",
            bw_sm = "results/XL_sites/all/{organism}/{unit}-IP-SM/SM-{unit}_{strand}_raw.bw",
            bw_filtered = "results/XL_sites/all/{organism}/{unit}-IP-SM/IP-{unit}_{strand}_signif.bw",
    params: orientation = "{strand}",
            count_threshold = 10
    threads: 1
    log: "logs/XL_sites/all/{organism}/{unit}_IP_SM_{strand}_peaks.log"
    wrapper: "file:cRIP-eCLIP-workflow/wrapper/XL_sites"

rule XL_sites_single:
    input: ip = "results/dedup/{organism}/IP-{unit}/IP-{unit}_{strand}.bam",   
           sm = "results/dedup/{organism}/SM-{unit}/SM-{unit}_{strand}.bam",
           ip_bed = "results/XL_sites_bedgraph/single/{organism}/IP-{unit}/IP-{unit}_{strand}.bedgraph",   
           sm_bed = "results/XL_sites_bedgraph/single/{organism}/SM-{unit}/SM-{unit}_{strand}.bedgraph",
           bed = "results/confirmed_peaks/single/{organism}/{unit}-IP-SM/{unit}_IP_SM_{strand}_signif.bed"
    output: tsv_raw = "results/XL_sites/single/{organism}/{unit}-IP-SM/{unit}_IP_SM_{strand}_raw.tsv",
            tsv_filtered = "results/XL_sites/single/{organism}/{unit}-IP-SM/{unit}_IP_SM_{strand}_signif.tsv",
            bed = "results/XL_sites/single/{organism}/{unit}-IP-SM/{unit}_IP_SM_{strand}_signif.bed",
            bw_ip = "results/XL_sites/single/{organism}/{unit}-IP-SM/IP-{unit}_{strand}_raw.bw",
            bw_sm = "results/XL_sites/single/{organism}/{unit}-IP-SM/SM-{unit}_{strand}_raw.bw",
            bw_filtered = "results/XL_sites/single/{organism}/{unit}-IP-SM/IP-{unit}_{strand}_signif.bw",
    params: orientation = "{strand}",
            count_threshold = 10
    threads: 1
    log: "logs/XL_sites/single/{organism}/{unit}_IP_SM_{strand}_peaks.log"
    wrapper: "file:cRIP-eCLIP-workflow/wrapper/XL_sites"

rule XL_sites_DBA:
    input: ip = "results/dedup/{organism}/IP-{unit}/IP-{unit}_{strand}.bam",   
           sm = "results/dedup/{organism}/SM-{unit}/SM-{unit}_{strand}.bam",
           ip_bed = "results/XL_sites_bedgraph/DBA/{organism}/IP-{unit}/IP-{unit}_{strand}.bedgraph",   
           sm_bed = "results/XL_sites_bedgraph/DBA/{organism}/SM-{unit}/SM-{unit}_{strand}.bedgraph",
           bed = lambda wildcards: "results/DBA/{organism}/" + sample_comp['unit2sample'][wildcards.unit] + "/" + sample_comp['unit2sample'][wildcards.unit] + "_{strand}_signif.bed"
    output: tsv_raw = "results/XL_sites/DBA/{organism}/{unit}-IP-SM/{unit}_IP_SM_{strand}_raw.tsv",
            tsv_filtered = "results/XL_sites/DBA/{organism}/{unit}-IP-SM/{unit}_IP_SM_{strand}_signif.tsv",
            bed = "results/XL_sites/DBA/{organism}/{unit}-IP-SM/{unit}_IP_SM_{strand}_signif.bed",
            bw_ip = "results/XL_sites/DBA/{organism}/{unit}-IP-SM/IP-{unit}_{strand}_raw.bw",
            bw_sm = "results/XL_sites/DBA/{organism}/{unit}-IP-SM/SM-{unit}_{strand}_raw.bw",
            bw_filtered = "results/XL_sites/DBA/{organism}/{unit}-IP-SM/IP-{unit}_{strand}_signif.bw",
    params: orientation = "{strand}",
            count_threshold = 10
    threads: 1
    log: "logs/XL_sites/DBA/{organism}/{unit}_IP_SM_{strand}_peaks.log"
    wrapper: "file:cRIP-eCLIP-workflow/wrapper/XL_sites"

rule DBA_XL_sites:
    input: ctrl = lambda wildcards: "results/dedup/{organism}/IP-" + sample_comp['DBA'][wildcards.unit]['ctrl'] + "/IP-" + sample_comp['DBA'][wildcards.unit]['ctrl'] + "_{strand}.bam",
           ko = lambda wildcards: "results/dedup/{organism}/IP-" + sample_comp['DBA'][wildcards.unit]['ko'] + "/IP-" + sample_comp['DBA'][wildcards.unit]['ko'] + "_{strand}.bam",
           ctrl_bed = lambda wildcards: "results/XL_sites_bedgraph/DBA/{organism}/IP-" + sample_comp['DBA'][wildcards.unit]['ctrl'] + "/IP-" + sample_comp['DBA'][wildcards.unit]['ctrl'] + "_{strand}.bedgraph",   
           ko_bed = lambda wildcards: "results/XL_sites_bedgraph/DBA/{organism}/IP-" + sample_comp['DBA'][wildcards.unit]['ko'] + "/IP-" + sample_comp['DBA'][wildcards.unit]['ko'] + "_{strand}.bedgraph"
    output: tsv_raw = "results/DBA_XL_sites/{organism}/{unit}/{unit}_{strand}_raw.tsv",
            tsv_filtered = "results/DBA_XL_sites/{organism}/{unit}/{unit}_{strand}_signif.tsv",
            bed = "results/DBA_XL_sites/{organism}/{unit}/{unit}_{strand}_signif.bed",
            bw_filtered_ctrl = "results/DBA_XL_sites/{organism}/{unit}/CTRL-{unit}_{strand}_signif.bw",
            bw_filtered_ko = "results/DBA_XL_sites/{organism}/{unit}/KO-{unit}_{strand}_signif.bw"
    params: orientation = "{strand}",
            count_threshold = 10
    threads: 1
    log: "logs/XL_sites_from_pub_peaks/{organism}/{unit}_IP_SM_{strand}_peaks.log"
    wrapper: "file:cRIP-eCLIP-workflow/wrapper/DEG_XL_sites"

rule bw_tracks_norm_XL_sites:
    input: ip_bam = expand("results/dedup/{{organism}}/IP-{{unit}}/IP-{{unit}}_{strand}.bam", strand = ['fwd', 'rev']),
           sm_bam = expand("results/dedup/{{organism}}/SM-{{unit}}/SM-{{unit}}_{strand}.bam", strand = ['fwd', 'rev']),
           ip_bw = expand("results/XL_sites/{{suffix}}/{{organism}}/{{unit}}-IP-SM/IP-{{unit}}_{strand}_raw.bw", strand = ['fwd', 'rev']),
           sm_bw = expand("results/XL_sites/{{suffix}}/{{organism}}/{{unit}}-IP-SM/SM-{{unit}}_{strand}_raw.bw", strand = ['fwd', 'rev']),
    output: raw = expand("results/normed_tracks_XL_sites/{{suffix}}/{{organism}}/Raw/{{unit}}/{{unit}}_RAW_IP_{strand}.bw", strand = ['fwd', 'rev']),
            rpm = expand("results/normed_tracks_XL_sites/{{suffix}}/{{organism}}/RPM/{{unit}}/{{unit}}_RPM_IP_{strand}.bw", strand = ['fwd', 'rev']),
            subtract = expand("results/normed_tracks_XL_sites/{{suffix}}/{{organism}}/Subtract/{{unit}}/{{unit}}_SUBTRACT_{strand}.bw", strand = ['fwd', 'rev']),
            kl = expand("results/normed_tracks_XL_sites/{{suffix}}/{{organism}}/KL/{{unit}}/{{unit}}_KL_{strand}.bw", strand = ['fwd', 'rev']),
            lfc = expand("results/normed_tracks_XL_sites/{{suffix}}/{{organism}}/LFC/{{unit}}/{{unit}}_LFC_{strand}.bw", strand = ['fwd', 'rev'])
    threads: 1
    params: raw = "results/normed_tracks_XL_sites/{suffix}/{organism}/Raw/{unit}/{unit}_RAW",
            rpm = "results/normed_tracks_XL_sites/{suffix}/{organism}/RPM/{unit}/{unit}_RPM",
            subtract = "results/normed_tracks_XL_sites/{suffix}/{organism}/Subtract/{unit}/{unit}_SUBTRACT",
            kl = "results/normed_tracks_XL_sites/{suffix}/{organism}/KL/{unit}/{unit}_KL",
            lfc = "results/normed_tracks_XL_sites/{suffix}/{organism}/LFC/{unit}/{unit}_LFC",
            onlyXL = True
    log: "logs/bw_tracks_norm_XL_sites_{suffix}/{organism}/{unit}_xl_sites.log"
    wrapper: "file:cRIP-eCLIP-workflow/wrapper/norm_tracks"